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
from requests import HTTPError


class Incomplete(Exception):
    def __init__(self, message, absent_identifiers):
        super(Incomplete, self).__init__(message)
        self.absent_identifiers = absent_identifiers


class IncompleteFragments(Incomplete):

    def __init__(self, absent_identifiers, fragments):
        """List of fragments and list of identifiers for which no information could be found

        Args:
            absent_identifiers (List[str]): List of identifiers for which no information could be found
            fragments (List[dict]): List of fragment information that could be retrieved
        """
        message = 'Some identifiers could not be found'
        super(IncompleteFragments, self).__init__(message, absent_identifiers)
        self.fragments = fragments


class IncompletePharmacophores(Incomplete):
    def __init__(self, absent_identifiers, pharmacophores):
        """List of fragments and list of identifiers for which no information could be found

        Args:
            absent_identifiers (List[str]): List of identifiers for which no information could be found
            pharmacophores (List[dict]): List of pharmacophores that could be retrieved
        """
        message = 'Some identifiers could not be found'
        super(IncompletePharmacophores, self).__init__(message, absent_identifiers)
        self.pharmacophores = pharmacophores


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

        Raises:
            request.HTTPError: When fragment_id could not be found
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
            chunk_size (int): Number of fragment to retrieve in a single http request

        Returns:
            list[dict]: List of fragment information

        Raises:
            IncompleteFragments: When one or more of the identifiers could not be found.
        """
        return self._fetch_chunked_fragments('fragment_ids', fragment_ids, chunk_size)

    def _fetch_chunked_fragments(self, idtype, ids, chunk_size):
        fragments = []
        absent_identifiers = []
        for start in range(0, len(ids), chunk_size):
            stop = chunk_size + start
            (chunk_fragments, chunk_absent_identifiers) = self._fetch_fragments(idtype, ids[start:stop])
            fragments += chunk_fragments
            absent_identifiers += chunk_absent_identifiers
        if chunk_absent_identifiers:
            raise IncompleteFragments(absent_identifiers, fragments)
        return fragments

    def _fetch_fragments(self, idtype, ids):
        url = self.base_url + '/fragments?{idtype}={ids}'.format(idtype=idtype, ids=','.join(ids))
        absent_identifiers = []
        try:
            response = requests.get(url)
            response.raise_for_status()
            fragments = response.json()
        except HTTPError as e:
            if e.response.status_code == 404:
                body = e.response.json()
                fragments = body['fragments']
                absent_identifiers = body['absent_identifiers']
            else:
                raise e
        # Convert molblock string to RDKit Mol object
        for fragment in fragments:
            if fragment['mol'] is not None:
                fragment['mol'] = MolFromMolBlock(fragment['mol'])
        return fragments, absent_identifiers

    def pharmacophores(self, fragment_ids):
        absent_identifiers = []
        pharmacophores = []
        for fragment_id in fragment_ids:
            url = self.base_url + '/fragments/{0}.phar'.format(fragment_id)
            try:
                response = requests.get(url)
                response.raise_for_status()
                pharmacophore = response.text
                pharmacophores.append(pharmacophore)
            except HTTPError as e:
                if e.response.status_code == 404:
                    pharmacophores.append(None)
                    absent_identifiers.append(fragment_id)
                else:
                    raise e
        if absent_identifiers:
            raise IncompletePharmacophores(absent_identifiers, pharmacophores)
        return pharmacophores
