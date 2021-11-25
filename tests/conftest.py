# encoding: utf-8
#
# conftest.py

"""
Here you can add fixtures that will be used for all the tests in this
directory. You can also add conftest.py files in underlying subdirectories.
Those conftest.py will only be applies to the tests in that subdirectory and
underlying directories. See https://docs.pytest.org/en/2.7.3/plugins.html for
more information.
"""

import asyncio
import uuid

import pytest

from proto.actor.actor import ProtoActor

@pytest.fixture(scope="session")
def event_loop():
    return asyncio.get_event_loop()


@pytest.fixture(scope="session")
async def proto_test_actor(event_loop):

    actor = ProtoActor(name=f"proto_test-{uuid.uuid4().hex[:8]}")
    await actor.start()
   
    yield actor

    await actor.stop()

