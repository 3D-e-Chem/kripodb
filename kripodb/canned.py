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
"""Module with functions which use pandas DataFrame as input and output.

For using Kripo data files inside KNIME (http://www.knime.org)
"""

from __future__ import absolute_import

import numpy as np
import pandas as pd
from requests import HTTPError

from .db import FragmentsDb
from .pairs import similar, open_similarity_matrix
from .pharmacophores import PharmacophoresDb, as_phar
from .webservice.client import WebserviceClient, IncompleteFragments, IncompletePharmacophores


class IncompleteHits(Exception):
    def __init__(self, absent_identifiers, hits):
        """List of hits and list of identifiers for which no information could be found

        Args:
            absent_identifiers (List[str]): List of identifiers for which no information could be found
            hits (pandas.DataFrame): Data frame with query_fragment_id, hit_frag_id and score columns
        """
        message = 'Some query fragment identifiers could not be found'
        super(IncompleteHits, self).__init__(message)
        self.absent_identifiers = absent_identifiers
        self.hits = hits


def similarities(queries, similarity_matrix_filename_or_url, cutoff, limit=1000):
    """Find similar fragments to queries based on similarity matrix.

    Args:
        queries (List[str]): Query fragment identifiers
        similarity_matrix_filename_or_url (str): Filename of similarity matrix file or base url of kripodb webservice
        cutoff (float): Cutoff, similarity scores below cutoff are discarded.
        limit (int): Maximum number of hits for each query.
            Default is 1000. Use is None for no limit.

    Examples:
        Fragments similar to '3j7u_NDP_frag24' fragment.

        >>> import pandas as pd
        >>> from kripodb.canned import similarities
        >>> queries = pd.Series(['3j7u_NDP_frag24'])
        >>> hits = similarities(queries, 'data/similaritys.h5', 0.55)
        >>> len(hits)
        11

        Retrieved from web service instead of local similarity matrix file.
        Make sure the web service is running,
        for example by `kripodb serve data/similarities.h5 data/fragments.sqlite data/pharmacophores.h5`.

        >>> hits = similarities(queries, 'http://localhost:8084/kripo', 0.55)
        >>> len(hits)
        11

    Returns:
        pandas.DataFrame: Data frame with query_fragment_id, hit_frag_id and score columns

    Raises:
        IncompleteHits: When one or more of the identifiers could not be found.
    """
    hits = []
    absent_identifiers = []
    if similarity_matrix_filename_or_url.startswith('http'):
        client = WebserviceClient(similarity_matrix_filename_or_url)
        for query in queries:
            try:
                qhits = client.similar_fragments(query, cutoff, limit)
                hits.extend(qhits)
            except HTTPError as e:
                if e.response.status_code == 404:
                    absent_identifiers.append(query)
    else:
        similarity_matrix = open_similarity_matrix(similarity_matrix_filename_or_url)
        for query in queries:
            try:
                for query_id, hit_id, score in similar(query, similarity_matrix, cutoff, limit):
                    hit = {'query_frag_id': query_id,
                           'hit_frag_id': hit_id,
                           'score': score,
                           }
                    hits.append(hit)
            except KeyError:
                absent_identifiers.append(query)

        similarity_matrix.close()

    if absent_identifiers:
        if len(hits) > 0:
            df = pd.DataFrame(hits, columns=['query_frag_id', 'hit_frag_id', 'score'])
        else:
            # empty hits array will give dataframe without columns
            df = pd.DataFrame({'query_frag_id': pd.Series(dtype=str),
                               'hit_frag_id': pd.Series(dtype=str),
                               'score': pd.Series(dtype=np.double)
                               }, columns=['query_frag_id', 'hit_frag_id', 'score'])
        raise IncompleteHits(absent_identifiers, df)

    return pd.DataFrame(hits, columns=['query_frag_id', 'hit_frag_id', 'score'])


def fragments_by_pdb_codes(pdb_codes, fragments_db_filename_or_url, prefix=''):
    """Retrieve fragments based on PDB codes.

    See http://www.rcsb.org/pdb/ for PDB structures.

    Args:
        pdb_codes (List[str]): List of PDB codes
        fragments_db_filename_or_url (str): Filename of fragments db or base url of kripodb webservice
        prefix (str): Prefix for output columns

    Examples:
        Fetch fragments of '2n2k' PDB code

        >>> from kripodb.canned import fragments_by_pdb_codes
        >>> pdb_codes = pd.Series(['2n2k'])
        >>> fragments = fragments_by_pdb_codes(pdb_codes, 'data/fragments.sqlite')
        >>> len(fragments)
        3

        Retrieved from web service instead of local fragments db file.
        Make sure the web service is running,
        for example by `kripodb serve data/similarities.h5 data/fragments.sqlite data/pharmacophores.h5`.

        >>> fragments = fragments_by_pdb_codes(pdb_codes, 'http://localhost:8084/kripo')
        >>> len(fragments)
        3

    Returns:
        pandas.DataFrame: Data frame with fragment information

    Raises:
        IncompleteFragments: When one or more of the identifiers could not be found.
    """
    if fragments_db_filename_or_url.startswith('http'):
        client = WebserviceClient(fragments_db_filename_or_url)
        try:
            fragments = client.fragments_by_pdb_codes(pdb_codes)
        except IncompleteFragments as e:
            df = pd.DataFrame(e.fragments)
            df.rename(columns=lambda x: prefix + x, inplace=True)
            raise IncompleteFragments(e.absent_identifiers, df)
    else:
        fragmentsdb = FragmentsDb(fragments_db_filename_or_url)
        fragments = []
        absent_identifiers = []
        for pdb_code in pdb_codes:
            try:
                for fragment in fragmentsdb.by_pdb_code(pdb_code):
                    fragments.append(fragment)
            except LookupError as e:
                absent_identifiers.append(pdb_code)
        if absent_identifiers:
            df = pd.DataFrame(fragments)
            df.rename(columns=lambda x: prefix + x, inplace=True)
            raise IncompleteFragments(absent_identifiers, df)

    df = pd.DataFrame(fragments)
    df.rename(columns=lambda x: prefix + x, inplace=True)
    return df


def fragments_by_id(fragment_ids, fragments_db_filename_or_url, prefix=''):
    """Retrieve fragments based on fragment identifier.

    Args:
        fragment_ids (List[str]): List of fragment identifiers
        fragments_db_filename_or_url (str): Filename of fragments db or base url of kripodb webservice
        prefix (str): Prefix for output columns

    Examples:
        Fetch fragments of '2n2k_MTN_frag1' fragment identifier

        >>> from kripodb.canned import fragments_by_id
        >>> fragment_ids = pd.Series(['2n2k_MTN_frag1'])
        >>> fragments = fragments_by_id(fragment_ids, 'data/fragments.sqlite')
        >>> len(fragments)
        1

        Retrieved from web service instead of local fragments db file.
        Make sure the web service is running,
        for example by `kripodb serve data/similarities.h5 data/fragments.sqlite data/pharmacophores.h5`.

        >>> fragments = fragments_by_id(fragment_ids,, 'http://localhost:8084/kripo')
        >>> len(fragments)
        1

    Returns:
        pandas.DataFrame: Data frame with fragment information

    Raises:
        IncompleteFragments: When one or more of the identifiers could not be found.
    """
    if fragments_db_filename_or_url.startswith('http'):
        client = WebserviceClient(fragments_db_filename_or_url)
        try:
            fragments = client.fragments_by_id(fragment_ids)
        except IncompleteFragments as e:
            df = pd.DataFrame(e.fragments)
            df.rename(columns=lambda x: prefix + x, inplace=True)
            raise IncompleteFragments(e.absent_identifiers, df)
    else:
        fragmentsdb = FragmentsDb(fragments_db_filename_or_url)
        fragments = []
        absent_identifiers = []
        for frag_id in fragment_ids:
            try:
                fragments.append(fragmentsdb[frag_id])
            except KeyError:
                absent_identifiers.append(frag_id)
        if absent_identifiers:
            df = pd.DataFrame(fragments)
            df.rename(columns=lambda x: prefix + x, inplace=True)
            raise IncompleteFragments(absent_identifiers, df)

    df = pd.DataFrame(fragments)
    df.rename(columns=lambda x: prefix + x, inplace=True)
    return df


def pharmacophores_by_id(fragment_ids, pharmacophores_db_filename_or_url):
    """Fetch pharmacophore points by fragment identifiers

    Args:
        fragment_ids (pd.Series): List of fragment identifiers
        pharmacophores_db_filename_or_url: Filename of pharmacophores db or base url of kripodb webservice

    Returns:
        pandas.Series: Pandas series with pharmacophores as string in phar format.
                       Fragment without pharmacophore will return None

    Examples:
        Fragments similar to '3j7u_NDP_frag24' fragment.

        >>> from kripodb.canned import pharmacophores_by_id
        >>> fragment_ids = pd.Series(['2n2k_MTN_frag1'], ['Row0'])
        >>> pharmacophores = pharmacophores_by_id(fragment_ids, 'data/pharmacophores.h5')
        >>> len(pharmacophores)
        1

        Retrieved from web service instead of local pharmacophores db file.
        Make sure the web service is running,
        for example by `kripodb serve data/similarities.h5 data/fragments.sqlite data/pharmacophores.h5`.

        >>> pharmacophores = pharmacophores_by_id(fragment_ids, 'http://localhost:8084/kripo')
        >>> len(pharmacophores)
        1
    """
    pphors = pd.Series([], dtype=str)
    if pharmacophores_db_filename_or_url.startswith('http'):
        client = WebserviceClient(pharmacophores_db_filename_or_url)
        try:
            pphorsarray = client.pharmacophores(fragment_ids)
            pphors = pd.Series(pphorsarray, fragment_ids.index, dtype=str)
        except IncompletePharmacophores as e:
            pphors = pd.Series(e.pharmacophores, fragment_ids.index, dtype=str)
            raise IncompletePharmacophores(e.absent_identifiers, pphors)
    else:
        with PharmacophoresDb(pharmacophores_db_filename_or_url) as pharmacophoresdb:
            absent_identifiers = []
            for row_id, frag_id in fragment_ids.iteritems():
                try:
                    phar = as_phar(frag_id, pharmacophoresdb[frag_id])
                    pphors[row_id] = phar
                except KeyError:
                    pphors[row_id] = None
                    absent_identifiers.append(frag_id)
            if absent_identifiers:
                raise IncompletePharmacophores(absent_identifiers, pphors)
    return pphors
