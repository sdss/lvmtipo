# -*- coding: utf-8 -*-
#
# @Author: Florian Briegel (briegel@mpia.de)
# @Date: 2021-08-18
# @Filename: lvm/actors.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

from cluplus.proxy import Proxy


class lvm:
    class sci:
        foc = Proxy( "lvm.sci.foc")
        km = Proxy( "lvm.sci.km")
        pwi = Proxy( "lvm.sci.pwi")
        agc = Proxy( "lvm.sci.agcam")
        ag = Proxy( "lvm.sci.ag")

        async def start():
            rc = await asyncio.gather(
                lvm.sci.foc.start(),
                lvm.sci.km.start(),
                lvm.sci.pwi.start(),
                lvm.sci.agc.start(),
                lvm.sci.ag.start(),
                return_exceptions=True
            )
            return lvm.sci


    class skye:
        foc = Proxy( "lvm.skye.foc")
        km = Proxy( "lvm.skye.km")
        pwi = Proxy( "lvm.skye.pwi")
        agc = Proxy( "lvm.skye.agcam")
        ag = Proxy( "lvm.skye.ag")
        
        async def start():
            await invoke(
                lvm.skye.foc.start(),
                lvm.skye.km.start(),
                lvm.skye.pwi.start(),
                lvm.skye.agc.start(),
                lvm.skye.ag.start(),
                return_exceptions=True
            )
            return lvm.skye


    class skyw:
        foc = Proxy( "lvm.skyw.foc")
        km = Proxy( "lvm.skyw.km")
        pwi = Proxy( "lvm.skyw.pwi")
        agc = Proxy( "lvm.skyw.agcam")
        ag = Proxy( "lvm.skyw.ag")
        async def start():
            await lvm_amqpc.start()
            await invoke(
                lvm.skyw.foc.start(),
                lvm.skyw.km.start(),
                lvm.skyw.pwi.start(),
                lvm.skyw.agc.start(),
                lvm.skyw.ag.start(),
                return_exceptions=True
            )
            return lvm.skyw


    class spec:
        foc = Proxy( "lvm.spec.foc")
        fibsel = Proxy( "lvm.spec.fibsel")
        pwi = Proxy( "lvm.spec.pwi")
        agc = Proxy( "lvm.spec.agcam")
        ag = Proxy( "lvm.spec.ag")
        async def start():
            await lvm_amqpc.start()
            await invoke(
                lvm.spec.foc.start(),
                lvm.spec.fibsel.start(),
                lvm.spec.pwi.start(),
                lvm.spec.agc.start(),
                lvm.spec.ag.start(),
                return_exceptions=True
            )
            return lvm.spec

    def from_string(subsys: str):
        if subsys == 'sci': 
            return lvm.sci.start()
        elif subsys == 'skye': 
            return lvm.skye.start()
        elif subsys == 'skyw': 
            return lvm.skyw.start()
        elif subsys == 'spec': 
            return lvm.spec.start()
        else: return None

