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
from pyramid import testing
from webtest import TestApp
from kripodb import webservice
from kripodb.hdf5 import DistanceMatrix
from kripodb.version import __version__


class TestWebservice(object):
    def setUp(self):
        self.matrix = DistanceMatrix('data/distances.h5')
        self.request = testing.DummyRequest()
        self.request.registry.settings = {'matrix': self.matrix}

    def tearDown(self):
        self.matrix.close()

    def test_similar(self):
        self.request.swagger_data = {
            'fragment_id': '3j7u_NDP_frag24',
            'cutoff': 0.85,
        }

        result = webservice.similar(self.request)

        expected = [
            {'query_frag_id': '3j7u_NDP_frag24', 'hit_frag_id': '3j7u_NDP_frag23', 'score': 0.8991},
        ]
        eq_(result, expected)

    def test_version(self):
        result = webservice.version(self.request)

        expected = {'version': __version__}
        eq_(result, expected)

    def test_wsgi_app(self):
        app = webservice.wsgi_app(self.request.registry.settings)
        testapp = TestApp(app)

        result = testapp.get('/kripo/version')

        expected = {'version': __version__}
        eq_(result.json_body, expected)