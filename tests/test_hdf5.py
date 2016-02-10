# Copyright 2016 Netherlands eScience Center
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

from nose.tools import eq_

from kripodb.hdf5 import DistanceMatrix


class TestPairsTable(object):
    matrix = None

    def setUp(self):
        self.matrix = DistanceMatrix('data/distances.h5')

    def tearDown(self):
        self.matrix.close()

    def test_find_1(self):
        frag_id = self.matrix.labels().by_label('2n6i_4FU_frag1')

        result = self.matrix.pairs().find(frag_id, 0.98)

        expected = {357742L: 1.0, 357686L: 1.0}
        eq_(result, expected)
