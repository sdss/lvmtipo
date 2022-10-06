from keyreader import KeyReader

async def keep_alive():
  keyreader = KeyReader(echo=False, block=False)
  while not keyreader.getch():
     await asyncio.sleep(0.42)

