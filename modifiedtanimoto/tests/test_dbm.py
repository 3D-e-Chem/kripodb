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

from intbitset import intbitset
from mock import Mock
from modifiedtanimoto.dbm import IntbitsetDict


class TestIntbitsetDict:
    bitsetA = None
    samples = None
    bitsets = None
    number_of_bits = 100

    def setup(self):
        self.bitsetA = intbitset([2, 5, 8, 23])
        self.samples = {
            'a': self.bitsetA.fastdump()
        }
        self.number_of_bits = 100

        self.bitsets = IntbitsetDict(self.samples, self.number_of_bits)

    def test_keys(self):
        result = self.bitsets.keys()

        expected = ['a']
        assert result == expected

    def test_len(self):
        assert len(self.bitsets) == 1

    def test_contains(self):
        assert 'a' in self.bitsets

    def test_contains_keyerror(self):
        assert 'b' not in self.bitsets

    def test_number_of_bits(self):
        assert self.bitsets.number_of_bits == self.number_of_bits

    def test_get(self):
        result = self.bitsets.get('a')
        assert result == self.bitsetA

    def test_get_default(self):
        result = self.bitsets.get('b', self.bitsetA)
        assert result == self.bitsetA

    def test_getitem(self):
        result = self.bitsets['a']
        assert result == self.bitsetA

    def test_set(self):
        self.bitsets['b'] = self.bitsetA

        assert self.bitsets['b'] == self.bitsetA

    def test_del(self):
        del(self.bitsets['a'])

        assert len(self.bitsets) == 0

    def test_close_dbm(self):
        d = Mock()
        z = IntbitsetDict(d, self.number_of_bits)

        z.close()

        assert d.close.called

    def test_close_dict(self):
        z = IntbitsetDict({}, self.number_of_bits)

        try:
            z.close()
        except AttributeError:
            assert False
        assert True

    def test_iter(self):
        result = [k for k in self.bitsets]
        expected = ['a']
        assert result == expected

    def test_with(self):
        d = Mock()

        with IntbitsetDict(d, self.number_of_bits):
            pass

        assert d.close.called
