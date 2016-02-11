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

For using Kripo data files inside Knime (http://www.knime.org)
"""

import pandas as pd

from .hdf5 import DistanceMatrix
from .pairs import similar
from .db import FragmentsDb


def similarities(queries, distance_matrix_filename, cutoff, limit=1000):
    """Find similar fragments to queries based on distance matrix.

    Args:
        queries (List[str]): Query fragment identifiers
        distance_matrix_filename (str): Filename of distance matrix file
        cutoff (float): Cutoff, distance scores below cutoff are discarded.
        limit (int): Maximum number of hits for each query.
            Default is 1000. Use is None for no limit.

    Examples:
        Fragments similar to '3j7u_NDP_frag24' fragment.

        >>> from kripodb.canned import similarities
        >>> queries = pd.Series(['3j7u_NDP_frag24'])
        >>> hits = similarities(queries, 'data/distances.h5', 0.55)
        >>> len(hits)
        5

    Returns:
        pandas.DataFrame: Data frame with query_fragment_id, hit_frag_id and score columns
    """
    distance_matrix = DistanceMatrix(distance_matrix_filename)
    hits = []
    for query in queries:
        for query_id, hit_id, score in similar(query, distance_matrix, cutoff, limit):
            hit = {'query_frag_id': query_id,
                   'hit_frag_id': hit_id,
                   'score': score,
                   }
            hits.append(hit)

    distance_matrix.close()

    return pd.DataFrame(hits)


def fragments_by_pdb_codes(pdb_codes, fragments_db_filename):
    """Retrieve fragments based on PDB codes.

    See http://www.rcsb.org/pdb/ for PDB structures.

    Args:
        pdb_codes (List[str]): List of PDB codes
        fragments_db_filename (str): Filename of fragments db

    Examples:
        Fetch fragments of '2n2k' PDB code

        >>> from kripodb.canned import fragments_by_pdb_codes
        >>> pdb_codes = pd.Series(['2n2k'])
        >>> fragments = fragments_by_pdb_codes(pdb_codes, 'data/fragments.sqlite')
        >>> len(fragments)
        3

    Returns:
        pandas.DataFrame: Data frame with fragment information
    """
    fragmentsdb = FragmentsDb(fragments_db_filename)
    fragments = []
    for pdb_code in pdb_codes:
        for fragment in fragmentsdb.by_pdb_code(pdb_code):
            fragments.append(fragment)
    return pd.DataFrame(fragments)


def fragments_by_id(fragment_ids, fragments_db_filename):
    """Retrieve fragments based on fragment identifier.

    Args:
        fragment_ids (List[str]): List of fragment identifiers
        fragments_db_filename (str): Filename of fragments db

    Examples:
        Fetch fragments of '2n2k' PDB code

        >>> from kripodb.canned import fragments_by_id
        >>> fragment_ids = pd.Series(['2n2k_MTN_frag1'])
        >>> fragments = fragments_by_id(fragment_ids, 'data/fragments.sqlite')
        >>> len(fragments)
        1

    Returns:
        pandas.DataFrame: Data frame with fragment information
    """
    fragmentsdb = FragmentsDb(fragments_db_filename)
    fragments = [fragmentsdb[frag_id] for frag_id in fragment_ids]
    return pd.DataFrame(fragments)
