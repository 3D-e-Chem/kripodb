# Copyright 2016 Netherlands eScience Center
#
# Licensed under the Apache License, Version 2.0 (the 'License");
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
import requests


class WebserviceClient(object):
    def __init__(self, base_url):
        self.base_url = base_url

    def similar_fragments(self, fragment_id, cutoff, limit=1000):
        url = self.base_url + '/fragments/{fragment_id}/similar'.format(fragment_id=fragment_id)
        params = {'cutoff': cutoff, 'limit': limit}
        response = requests.get(url, params)
        response.raise_for_status()
        return response.json()