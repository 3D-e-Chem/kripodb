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

import StringIO
from intbitset import intbitset
import modifiedtanimoto.pairs as pairs


class Testpairs(object):
    bitsets = None
    pairs = None
    number_of_bits = None

    def setup(self):
        self.bitsets = {
            'a': intbitset([1, 2, 3]),
            'b': intbitset([1, 2, 4, 5, 8]),
            'c': intbitset([1, 2, 4, 8])
        }
        self.pairs = [
            ('a', 'c', 0.7523),
            ('b', 'c', 0.8342)]
        self.number_of_bits = 8

    def test_dump_pairs_tsv(self):
        out = StringIO.StringIO()

        pairs.dump_pairs_tsv(self.pairs, out)
        result = out.getvalue()

        expected = "a\tc\t0.7523\nb\tc\t0.8342\n"
        assert result == expected

    def test_dump_pairs_tsv_compact(self):
        label2id = {
            'a': 1,
            'b': 2,
            'c': 3
        }
        precision = 100
        out = StringIO.StringIO()

        pairs.dump_pairs_tsv_compact(self.pairs, label2id, precision, out)
        result = out.getvalue()

        expected = "1\t3\t75\n2\t3\t83\n"
        assert result == expected
