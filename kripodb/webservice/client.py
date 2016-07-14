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
"""Module for Client for kripo web service"""
from __future__ import absolute_import
import requests


class WebserviceClient(object):
    """Client for kripo web service

    Example:
        >>> client = WebserviceClient('http://localhost:8084/kripo')
        >>> client.similar_fragments('3j7u_NDP_frag24', 0.85)
        [{'query_frag_id': '3j7u_NDP_frag24', 'hit_frag_id': '3j7u_NDP_frag23', 'score': 0.8991}]

    Args:
        base_url (str): Base url of web service. e.g. http://localhost:8084/kripo
    """
    def __init__(self, base_url):
        self.base_url = base_url

    def similar_fragments(self, fragment_id, cutoff, limit=1000):
        """Find similar fragments to query.

        Args:
            fragment_id (str): Query fragment identifier
            cutoff (float): Cutoff, similarity scores below cutoff are discarded.
            limit (int): Maximum number of hits. Default is None for no limit.

        Returns:
            List(Dict()): Query fragment identifier, hit fragment identifier and similarity score
        """
        url = self.base_url + '/fragments/{fragment_id}/similar'.format(fragment_id=fragment_id)
        params = {'cutoff': cutoff, 'limit': limit}
        response = requests.get(url, params)
        response.raise_for_status()
        return response.json()
