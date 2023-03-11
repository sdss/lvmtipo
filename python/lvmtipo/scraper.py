# -*- coding: utf-8 -*-
#
# @Author: Florian Briegel (briegel@mpia.de)
# @Date: 2021-08-18
# @Filename: scraper.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import asyncio
import json
from typing import Callable, Optional

from datetime import datetime as dt
import aio_pika as apika

from clu.client import AMQPClient, AMQPReply, CommandStatus
from cluplus.proxy import Client, ProxyDict, unpack, flatten


class ScraperDataStore(object):
    """
    Scrapers data store initialized with a dictofdicts or yaml

    config =
           {'lvm.sci.pwi':
                {'ra_j2000_hours': 'sci:ra_j2000_h',
                 'dec_j2000_degs': 'sci:dec_j2000_d',
                 'altitude_degs': 'sci:altitude_d',
                 'azimuth_degs': 'sci:azimuth_d'},
            'lvm.sci.foc':
                {'Position': 'sci:foc_dt'},
            'lvm.sci.km':
                {'Position': 'sci:km_d', 'SkyPA': 'sci:sky_d'},
            'lvm.tel':
                {'temperature': 'bentemp',
                 'humidity': 'benhum',
                 'pressure': 'benpress'},
            'lvm.sci.agcam':
                {'east.temperature': 'sci:east.temp',
                 'east.filename': 'sci:east.file',
                 'west.temperature': 'sci:west.temp',
                 'west.filename': 'sci:west.file',
                 'center.temperature': 'sci:center.temp',
                 'center.filename': 'sci:center.file'}}

      # or
        import yaml

        config = "
        lvm.foc:
            Position: foc_dt

        lvm.km:
            Position: km_d
            SkyPA: sky_d

        lvm.tel:
            temperature: bentemp
            humidity: benhum
            pressure: benpress
        "
        yaml.safe_load(config_string)

        # or not storing data, only for callback usage
        config = "
        lvm.sci.ag: null
        "

        # or with a lambda funktion
        tel="sci"

        config = lambda tel: f"
        lvm.{tel}.foc:
            Position: foc_dt

        lvm.{tel}.km:
            Position: km_d
            SkyPA: sky_d

        lvm.{tel}.tel:
            temperature: bentemp
            humidity: benhum
            pressure: benpress

        lvm.{tel}.ag0: null
        "
        yaml.safe_load(config(tel))

    """

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
        if isinstance(akm, dict):
            self.data.update({akm[k]:(v, timestamp) for k, v in data.items() if k in akm.keys()})

    def items(self):
        return self.data.items()


class Scraper(Client):
    """ Scraper listening for actor messages """

    def __init__(
        self,
        scraper: dict,
        callback: Optional[Callable[[ProxyDict], None]] = None,
        **kwargs
    ):

        super().__init__(**kwargs)


        self.scraper_store = ScraperDataStore(scraper)
        self.callback = callback

    async def handle_reply(self, message: apika.IncomingMessage) -> AMQPReply:
        """Handles a reply received from the exchange.
        """
        reply = await super().handle_reply(message)

#        reply = AMQPReply(message, log=self.log)
        if reply.sender in self.scraper_store.actors() and reply.headers.get("message_code", None) in ":i":
            timestamp = apika.message.decode_timestamp(message.timestamp) if message.timestamp else datetime.utcnow()
            self.scraper_store.update_with_actor_key_maps(reply.sender, flatten(reply.body), timestamp)

            if self.callback:
                msg = ProxyDict(json.loads(reply.message.body))
                msg.sender = reply.sender
                msg.command_status = CommandStatus.code_to_status(reply.message_code)
                msg.timestamp = timestamp
                self.callback(msg)

        return reply

