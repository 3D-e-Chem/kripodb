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

import UserDict
import anydbm
from intbitset import intbitset


class IntbitsetDict(UserDict.DictMixin):
    """Dict of bitsets.

    Similar to Shelve in shelve package,
    but with intbitset fastload/fastdump instead of pickle
    """
    number_of_bits = None

    def __init__(self, dict, number_of_bits):
        self.dict = dict
        self.number_of_bits = number_of_bits

    def keys(self):
        return self.dict.keys()

    def __iter__(self):
        for k in self.dict.keys():
            yield k

    def __len__(self):
        return len(self.dict)

    def has_key(self, key):
        return key in self.dict

    def __contains__(self, key):
        return key in self.dict

    def get(self, key, default=None):
        if key in self.dict:
            return self[key]
        return default

    def __getitem__(self, key):
        serialized_value = self.dict[key]
        value = intbitset()
        value.fastload(serialized_value)
        return value

    def __setitem__(self, key, value):
        self.dict[key] = value.fastdump()

    def __delitem__(self, key):
        del self.dict[key]

    def close(self):
        try:
            self.dict.close()
        except AttributeError:
            pass

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()


class IntbitsetDictDbm(IntbitsetDict):
    """Bitset dictionary using dbm as storage instead of memory

    """
    def __init__(self, filename, number_of_bits, flag='c'):
        IntbitsetDict.__init__(self,
                               anydbm.open(filename, flag),
                               number_of_bits)
