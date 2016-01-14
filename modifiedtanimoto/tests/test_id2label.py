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

from StringIO import StringIO
from modifiedtanimoto.id2label import read, swap_label2id


def test_read():
    infile = StringIO('1\ta\n2\tb\n3\tc\n')

    result = read(infile)

    expected = {
        1: 'a',
        2: 'b',
        3: 'c'
    }
    assert result == expected


def test_swap_label2id():
    id2label = {
        1: 'a',
        2: 'b',
        3: 'c'
    }

    result = swap_label2id(id2label)

    expected = {
        'a': 1,
        'b': 2,
        'c': 3
    }
    assert result == expected
