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
import modifiedtanimoto.matrix as matrix


def test_dump_matrix_tsv():
    bitsets = {
        'a': intbitset([1, 2, 3]),
        'b': intbitset([1, 2, 4, 5, 8]),
        'c': intbitset([1, 2, 4, 8])
    }
    fp_size = 8
    out = StringIO.StringIO()

    matrix.dump_matrix_tsv(bitsets, bitsets, fp_size, 0.4, 0.5, 0.01, out)
    result = out.getvalue()
    print result

    expected = "a\tc\t0.0766666666667\n"

    assert result == expected


def test_dump_matrix_tsv_numbered():
    bitsets = {
        'a': intbitset([1, 2, 3]),
        'b': intbitset([1, 2, 4, 5, 8]),
        'c': intbitset([1, 2, 4, 8])
    }
    fp_size = 8
    label2id = {
        'a': 1,
        'b': 2,
        'c': 3
    }
    out = StringIO.StringIO()

    matrix.dump_matrix_tsv_numbered(bitsets, bitsets, fp_size, 0.4, 0.5, 0.01, label2id, out)
    result = out.getvalue()
    print result

    expected = "1\t3\t0.0766666666667\n"
    assert result == expected
