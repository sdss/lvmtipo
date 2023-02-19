import asyncio
import argparse

import yaml
from cluplus.proxy import ProxyDict, flatten, unpack
from lvmtipo.scraper import Scraper
from clu.client import CommandStatus

from PyQt5.QtCore import QObject, QThread, pyqtSignal, QEventLoop
from asyncqt import QEventLoop


def default(tel):
  return f"""
lvm.{tel}.pwi:
    ra_j2000_hours: {tel}:ra_j2000_h
    dec_j2000_degs: {tel}:dec_j2000_d
    altitude_degs: {tel}:altitude_d
    azimuth_degs: {tel}:azimuth_d

lvm.{tel}.foc:
    Position: {tel}:foc_dt

lvm.{tel}.km:
    Position: {tel}:km_d
    SkyPA: {tel}:sky_d

lvm.tel:
    temperature: bentemp
    humidity: benhum
    pressure: benpress

lvm.{tel}.agcam:
    east.temperature: {tel}:east.temp
    east.filename: {tel}:east.file
    west.temperature: {tel}:west.temp
    west.filename: {tel}:west.file
    center.temperature: {tel}:center.temp
    center.filename: {tel}:center.file
"""

tel = ['sci', 'skye', 'skyw', 'spec']


class QtLvmScraper(QObject):
    data_rcvd_signal = pyqtSignal(ProxyDict)

    def handle(self, data):
        # gets executed on scraper_event
        self.data_rcvd_signal.emit(data)

    def start(self, app, qtslot):
        loop = QEventLoop(app)
        asyncio.set_event_loop(loop)

        parser = argparse.ArgumentParser()
        parser.add_argument("-v", '--verbose', action='store_true', help="print some notes to stdout")
        parser.add_argument("-H", '--rmqhost', type=str, default="localhost", help="rabbitmq server")
        parser.add_argument("-t", '--telsubsys', type=str, default="sci", help="Telescope subsystem: sci, skye, skyw or spec")

        args = parser.parse_args()

        config = yaml.safe_load(default(args.telsubsys))

        scraper = Scraper(config, callback=self.handle, host=args.rmqhost)

        self.data_rcvd_signal.connect(qtslot)

        with loop:
            loop.run_until_complete(scraper.start())
            loop.run_forever()






