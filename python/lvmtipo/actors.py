# -*- coding: utf-8 -*-
#
# @Author: Florian Briegel (briegel@mpia.de)
# @Date: 2021-08-18
# @Filename: actors.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


import asyncio
from typing import Union

from cluplus.proxy import Proxy, invoke, unpack, flatten
from clu import BaseClient

class lvm:
    nps = Proxy("lvmnps")
    ieb = Proxy("lvmieb")
    ecp = Proxy("lvmecp")
    scp = Proxy("lvmscp")
#    scraper = Proxy("lvm.scraper")
    tel = Proxy("lvm.tel")


    class sci:
        foc = Proxy("lvm.sci.foc")
        km = Proxy("lvm.sci.km")
        pwi = Proxy("lvm.sci.pwi")
        agc = Proxy("lvm.sci.agcam")
        ag = Proxy("lvm.sci.ag")

        async def start(amqpc:BaseClient = None):
            await asyncio.gather(
                lvm.sci.foc.start(amqpc),
                lvm.sci.km.start(amqpc),
                lvm.sci.pwi.start(amqpc),
                lvm.sci.agc.start(amqpc),
                lvm.sci.ag.start(amqpc),
                return_exceptions=True
            )
            return lvm.sci

    class skye:
        foc = Proxy("lvm.skye.foc")
        km = Proxy("lvm.skye.km")
        pwi = Proxy("lvm.skye.pwi")
        agc = Proxy("lvm.skye.agcam")
        ag = Proxy("lvm.skye.ag")
        
        async def start(amqpc:BaseClient = None):
            await asyncio.gather(
                lvm.skye.foc.start(amqpc),
                lvm.skye.km.start(amqpc),
                lvm.skye.pwi.start(amqpc),
                lvm.skye.agc.start(amqpc),
                lvm.skye.ag.start(amqpc),
                return_exceptions=True
            )
            return lvm.skye

    class skyw:
        foc = Proxy("lvm.skyw.foc")
        km = Proxy("lvm.skyw.km")
        pwi = Proxy("lvm.skyw.pwi")
        agc = Proxy("lvm.skyw.agcam")
        ag = Proxy("lvm.skyw.ag")
        async def start(amqpc:BaseClient = None):
            await asyncio.gather(
                lvm.skyw.foc.start(amqpc),
                lvm.skyw.km.start(amqpc),
                lvm.skyw.pwi.start(amqpc),
                lvm.skyw.agc.start(amqpc),
                lvm.skyw.ag.start(amqpc),
                return_exceptions=True
            )
            return lvm.skyw

    class spec:
        foc = Proxy("lvm.spec.foc")
        fibsel = Proxy("lvm.spec.fibsel")
        pwi = Proxy("lvm.spec.pwi")
        agc = Proxy("lvm.spec.agcam")
        ag = Proxy("lvm.spec.ag")
        async def start(amqpc:BaseClient = None):
            await asyncio.gather(
                lvm.spec.foc.start(amqpc),
                lvm.spec.fibsel.start(amqpc),
                lvm.spec.pwi.start(amqpc),
                lvm.spec.agc.start(amqpc),
                lvm.spec.ag.start(amqpc),
                return_exceptions=True
            )
            return lvm.spec

    def from_string(subsys: str):
        if subsys == 'sci': 
            return lvm.sci
        elif subsys == 'skye': 
            return lvm.skye
        elif subsys == 'skyw': 
            return lvm.skyw
        elif subsys == 'spec': 
            return lvm.spec
        else: return None

    TelSubSystem = Union[sci, skye, skyw, spec]

