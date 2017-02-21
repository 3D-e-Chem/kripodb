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
from rdkit.Chem.AllChem import MolFromMolBlock


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
            list[dict]: Query fragment identifier, hit fragment identifier and similarity score
        """
        url = self.base_url + '/fragments/{fragment_id}/similar'.format(fragment_id=fragment_id)
        params = {'cutoff': cutoff, 'limit': limit}
        response = requests.get(url, params)
        response.raise_for_status()
        return response.json()

    def fragments_by_pdb_codes(self, pdb_codes, chunk_size=450):
        """Retrieve fragments by their PDB code

        Args:
            pdb_codes (List[str]): List of PDB codes
            chunk_size (int): Number of PDB codes to retrieve in a single http request

        Returns:
            list[dict]: List of fragment information

        Raises:
            requests.HTTPError: When one of the PDB codes could not be found.
        """
        return self._fetch_chunked_fragments('pdb_codes', pdb_codes, chunk_size)

    def fragments_by_id(self, fragment_ids, chunk_size=100):
        """Retrieve fragments by their identifier

        Args:
            fragment_ids (List[str]): List of fragment identifiers
            chunk_size (int): Number of PDB codes to retrieve in a single http request

        Returns:
            list[dict]: List of fragment information

        Raises:
            requests.HTTPError: When one of the identifiers could not be found.
        """
        return self._fetch_chunked_fragments('fragment_ids', fragment_ids, chunk_size)

    def _fetch_chunked_fragments(self, idtype, ids, chunk_size):
        fragments = []
        for start in range(0, len(ids), chunk_size):
            stop = chunk_size + start
            fragments += self._fetch_fragments(idtype, ids[start:stop])
        return fragments

    def _fetch_fragments(self, idtype, ids):
        url = self.base_url + '/fragments?{idtype}={ids}'.format(idtype=idtype, ids=','.join(ids))
        response = requests.get(url)
        response.raise_for_status()
        fragments = response.json()
        # Convert molblock string to RDKit Mol object
        for fragment in fragments:
            if fragment['mol'] is not None:
                fragment['mol'] = MolFromMolBlock(fragment['mol'])
        return fragments
