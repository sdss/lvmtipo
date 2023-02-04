# -*- coding: utf-8 -*-
#
# @Author: Florian Briegel (briegel@mpia.de)
# @Date: 2021-08-18
# @Filename: scraper.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import asyncio

from datetime import datetime as dt
import aio_pika as apika

from clu.client import AMQPClient, AMQPReply
from cluplus.proxy import flatten, Client


class ScraperDataStore(object):
    def __init__(self, config={}):
        self.actor_key_maps = config
        self.data = {}

    def copy(self):
        o = type(self).__new__(self.__class__)
        o.actor_key_maps = self.actor_key_maps.copy()
        o.data = self.data.copy()
        return o

    def __repr__(self):
        return self.data.__repr__()

    def keys(self):
         return self.data.keys()

    def __getitem__(self, key):
        return self.get(key)

    def __setitem__(self, key, value):
        return self.set(key, value)

    def actors(self):
        return list(self.actor_key_maps.keys()) if self.actor_key_maps else []

    def set(self, key, val, timestamp=dt.utcnow()):
        self.data[key] = (val, timestamp)

    def get(self, key, default=None):
        return self.data.get(key, (default, None))[0]

    def update(self, data:dict, timestamp=dt.utcnow()):
        self.data.update({k:(v, timestamp) for k, v in data.items()})

    def update_with_actor_key_maps(self, actor, data:dict, timestamp=dt.utcnow()):
        akm = self.actor_key_maps.get(actor, None)
        self.data.update({akm[k]:(v, timestamp) for k, v in data.items() if k in akm.keys()})

    def items(self):
        return self.data.items()


class Scraper(Client):

    def __init__(
        self,
        scraper:dict,
        **kwargs
    ):

        super().__init__(**kwargs)
        self.scraper_store = ScraperDataStore(scraper)

    async def handle_reply(self, message: apika.IncomingMessage) -> AMQPReply:
        """Handles a reply received from the exchange.
        """
        reply = AMQPReply(message, log=self.log)
        if reply.sender in self.scraper_store.actors() and reply.headers.get("message_code", None) in ":i":
            timestamp = apika.message.decode_timestamp(message.timestamp) if message.timestamp else datetime.utcnow()
            self.scraper_store.update_with_actor_key_maps(reply.sender, flatten(reply.body), timestamp)

        return reply

