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


class TestDistanceMatrix(object):
    matrix = None

    def setUp(self):
        self.matrix = DistanceMatrix('data/distances.h5')

    def tearDown(self):
        self.matrix.close()

    def test_find_1(self):
        result = self.matrix.find('2n6i_4FU_frag1', 0.98)

        expected = [('2n6i_4FU_frag2', 1.0), ('2n6i_4FU_frag6', 1.0)]
        eq_(list(result), expected)

    def test_iter_first2(self):
        myiter = iter(self.matrix)

        result = [myiter.next(), myiter.next()]

        expected = [('2mlm_2W7_frag1', '2mlm_2W7_frag2', 0.5877164873731594), ('2mlm_2W7_frag2', '3wvm_STE_frag1', 0.4633096818493935)]
        eq_(result, expected)

    def test_iter_last(self):
        myiter = iter(self.matrix)

        result = None
        for row in myiter:
            result = row

        expected = ('3wyl_3KB_frag20', '3wyl_3KB_frag21', 0.999496452277409)
        eq_(result, expected)

