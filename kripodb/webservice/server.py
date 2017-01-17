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
"""Kripo datafiles wrapped in a webservice"""
from __future__ import absolute_import

from kripodb.db import FragmentsDb
from pkg_resources import resource_filename
import logging

from rdkit.Chem.AllChem import Mol
from rdkit.Chem.AllChem import MolToMolBlock
from rdkit.Chem.Draw import rdMolDraw2D
from six.moves.urllib_parse import urlparse

import connexion
from flask import current_app, abort
from flask.json import JSONEncoder

from ..version import __version__
from ..pairs import open_similarity_matrix

LOGGER = logging.getLogger(__name__)


class KripodbJSONEncoder(JSONEncoder):
    """JSON encoder for KripoDB object types

    Copied from http://flask.pocoo.org/snippets/119/"""
    def default(self, obj):
        try:
            if isinstance(obj, Mol):
                return MolToMolBlock(obj).encode()
            iterable = iter(obj)
        except TypeError:
            pass
        else:
            return list(iterable)
        return JSONEncoder.default(self, obj)


def get_similar_fragments(fragment_id, cutoff, limit):
    """Find similar fragments to query.

    Args:
        fragment_id (str): Query fragment identifier
        cutoff (float): Cutoff, similarity scores below cutoff are discarded.
        limit (int): Maximum number of hits. Default is None for no limit.

    Returns:
        List(Dict()): Query fragment identifier, hit fragment identifier and similarity score

    Raises:
        werkzeug.exceptions.NotFound: When the fragments_id could not be found

    """
    similarity_matrix = current_app.config['matrix']
    query_id = fragment_id
    hits = []
    try:
        raw_hits = similarity_matrix.find(query_id, cutoff, limit)
        # add query column
        for hit_id, score in raw_hits:
            hits.append({'query_frag_id': query_id, 'hit_frag_id': hit_id, 'score': score})
    except LookupError as e:
        abort(404, 'Fragment with identifier \'{0}\' not found'.format(fragment_id))
    return hits


def get_fragments(fragment_ids=None, pdb_codes=None):
    """Retrieve fragments based on their identifier or PDB code.

    Args:
        fragment_ids (List[str]): List of fragment identifiers
        pdb_codes (List[str]): List of PDB codes

    Returns:
        List(Dict()): List of fragment information

    Raises:
        werkzeug.exceptions.NotFound: When one of the fragments_ids or pdb_code could not be found
    """
    fragments_db_filename = current_app.config['db_fn']
    with FragmentsDb(fragments_db_filename) as fragmentsdb:
        fragments = []
        if fragment_ids:
            for frag_id in fragment_ids:
                try:
                    fragments.append(fragmentsdb[frag_id])
                except LookupError:
                    abort(404, 'Fragment with identifier \'{0}\' not found'.format(frag_id))
        if pdb_codes:
            for pdb_code in pdb_codes:
                try:
                    for fragment in fragmentsdb.by_pdb_code(pdb_code):
                        fragments.append(fragment)
                except LookupError:
                    abort(404, 'Fragments with PDB code \'{0}\' not found'.format(pdb_code))
        # TODO if fragment_ids and pdb_codes are both None then return paged list of all fragments
        return fragments


def mol2svg(mol, width, height):
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg


def get_fragment_svg(fragment_id, width, height):
    fragments_db_filename = current_app.config['db_fn']
    with FragmentsDb(fragments_db_filename) as fragmentsdb:
        try:
            fragment = fragmentsdb[fragment_id]
            LOGGER.warning([fragment_id, width, height])
            mol = fragment['mol']
            return mol2svg(mol, width, height)
        except LookupError:
            abort(404, 'Fragment with identifier \'{0}\' not found'.format(fragment_id))


def get_version():
    """
    Returns:
        Dict: Version of web service
    """
    # TODO check if matrix is usable
    return {'version': __version__}


def wsgi_app(sim_matrix, frags_db_fn, external_url='http://localhost:8084/kripo'):
    """Create wsgi app

    Args:
        sim_matrix (SimilarityMatrix): Similarity matrix to use in webservice
        frags_db_fn (FragmentsDb): Fragment database filename
        external_url (str): URL which should be used in Swagger spec

    Returns:
        connexion.App
    """
    app = connexion.App(__name__)
    url = urlparse(external_url)
    swagger_file = resource_filename(__name__, 'swagger.json')
    app.add_api(swagger_file, base_path=url.path, arguments={'hostport': url.netloc, 'scheme': url.scheme})
    app.app.json_encoder = KripodbJSONEncoder
    app.app.config['matrix'] = sim_matrix
    app.app.config['db_fn'] = frags_db_fn
    return app


def serve_app(matrix, db, internal_port=8084, external_url='http://localhost:8084/kripo'):
    """Serve webservice forever

    Args:
        matrix: Filename of similarity matrix hdf5 file
        db: Filename of fragments database file
        internal_port: TCP port on which to listen
        external_url (str): URL which should be used in Swagger spec
    """
    sim_matrix = open_similarity_matrix(matrix)
    app = wsgi_app(sim_matrix, db, external_url)
    LOGGER.setLevel(logging.INFO)
    LOGGER.addHandler(logging.StreamHandler())
    LOGGER.info(' * Swagger spec at {}/swagger.json'.format(external_url))
    LOGGER.info(' * Swagger ui at {}/ui'.format(external_url))
    try:
        app.run(port=internal_port)
    finally:
        sim_matrix.close()
