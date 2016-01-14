# Copyright 2013 Netherlands eScience Center
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import anydbm
from intbitset import intbitset


def fastdump(bitsets, out):
    for bid, bitset in bitsets.iteritems():
        out[bid] = bitset.fastdump()
    return out


def dump_bitsets(bitsets, out_fn):
    db = anydbm.open(out_fn, 'c')
    fastdump(bitsets, db)
    db.close()


def fastload(db, bitsets):
    for bid, serial_bitset in db.iteritems():
        bs = intbitset()
        bs.fastload(serial_bitset)
        bitsets[bid] = bs
    return bitsets


def load_bitsets(fn):
    bitsets = {}
    db = anydbm.open(fn)
    fastload(db, bitsets)
    return bitsets


def dump_keys(fn):
    return anydbm.open(fn).keys()
