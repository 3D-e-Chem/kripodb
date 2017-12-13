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

import logging

import connexion
import flask
from flask import current_app
from flask.json import JSONEncoder
from pkg_resources import resource_filename
from rdkit.Chem.AllChem import Mol
from rdkit.Chem.AllChem import MolToMolBlock
from rdkit.Chem.Draw import rdMolDraw2D
from six.moves.urllib_parse import urlparse

from kripodb.pharmacophores import as_phar, PharmacophoresDb
from ..db import FragmentsDb
from ..pairs import open_similarity_matrix
from ..version import __version__

LOGGER = logging.getLogger(__name__)


class KripodbJSONEncoder(JSONEncoder):
    """JSON encoder for KripoDB object types

    Copied from http://flask.pocoo.org/snippets/119/"""
    def default(self, obj):
        try:
            if isinstance(obj, Mol):
                return MolToMolBlock(obj)
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
        list[dict]: List of dict with query fragment identifier, hit fragment identifier and similarity score

    Raises:
        werkzeug.exceptions.NotFound: When the fragments_id could not be found

    """
    similarity_matrix = current_app.config['similarities']
    query_id = fragment_id
    hits = []
    try:
        raw_hits = similarity_matrix.find(query_id, cutoff, limit)
        # add query column
        for hit_id, score in raw_hits:
            hits.append({'query_frag_id': query_id, 'hit_frag_id': hit_id, 'score': score})
    except LookupError:
        return fragment_not_found(fragment_id)
    return hits


def fragment_not_found(fragment_id):
    title = 'Not Found'
    description = 'Fragment with identifier \'{0}\' not found'.format(fragment_id)
    ext = {'identifier': fragment_id}
    return connexion.problem(404, title, description, ext=ext)


def get_fragments(fragment_ids=None, pdb_codes=None):
    """Retrieve fragments based on their identifier or PDB code.

    Args:
        fragment_ids (List[str]): List of fragment identifiers
        pdb_codes (List[str]): List of PDB codes

    Returns:
        list[dict]: List of fragment information

    Raises:
        werkzeug.exceptions.NotFound: When one of the fragments_ids or pdb_code could not be found
    """
    fragments_db_filename = current_app.config['fragments']
    with FragmentsDb(fragments_db_filename) as fragmentsdb:
        fragments = []
        missing_ids = []
        if fragment_ids:
            for frag_id in fragment_ids:
                try:
                    fragments.append(fragmentsdb[frag_id])
                except LookupError:
                    missing_ids.append(frag_id)

        if pdb_codes:
            for pdb_code in pdb_codes:
                try:
                    for fragment in fragmentsdb.by_pdb_code(pdb_code):
                        fragments.append(fragment)
                except LookupError:
                    missing_ids.append(pdb_code)
        # TODO if fragment_ids and pdb_codes are both None then return paged list of all fragments
        if missing_ids:
            title = 'Not found'
            label = 'identifiers'
            if pdb_codes:
                label = 'PDB codes'
            description = 'Fragments with {1} \'{0}\' not found'.format(','.join(missing_ids), label)
            # connexion.problem is using json.dumps instead of flask custom json encoder, so performing convert myself
            # TODO remove mol2string conversion when https://github.com/zalando/connexion/issues/266 is fixed
            for fragment in fragments:
                if fragment['mol']:
                    fragment['mol'] = MolToMolBlock(fragment['mol'])
            ext = {'absent_identifiers': missing_ids, 'fragments': fragments}
            return connexion.problem(404, title, description, ext=ext)
        return fragments


def mol2svg(mol, width, height):
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg


def get_fragment_svg(fragment_id, width, height):
    """2D drawing of fragment in SVG format

    Args:
        fragment_id (str): Fragment identifier
        width (int): Width of SVG in pixels
        height (int): Height of SVG in pixels

    Returns:
        flask.Response|connexion.lifecycle.ConnexionResponse: SVG document|problem
    """
    fragments_db_filename = current_app.config['fragments']
    with FragmentsDb(fragments_db_filename) as fragmentsdb:
        try:
            fragment = fragmentsdb[fragment_id]
            if not fragment['mol']:
                title = 'Not Found'
                description = 'Fragment with identifier \'{0}\' has no molblock'.format(fragment_id)
                ext = {'identifier': fragment_id}
                return connexion.problem(404, title, description, ext=ext)
            mol = fragment['mol']
            svg = mol2svg(mol, width, height)
            return flask.Response(svg, mimetype='image/svg+xml')
        except LookupError:
            return fragment_not_found(fragment_id)


def get_fragment_phar(fragment_id):
    """Pharmacophore in phar format of fragment

    Args:
        fragment_id (str): Fragment identifier

    Returns:
        flask.Response|connexion.lifecycle.ConnexionResponse: Pharmacophore|problem

    """
    pharmacophores_db = current_app.config['pharmacophores']
    try:
        points = pharmacophores_db[fragment_id]
        phar = as_phar(fragment_id, points)
        return flask.Response(phar, mimetype='text/plain')
    except LookupError:
        return fragment_not_found(fragment_id)


def get_version():
    """
    Returns:
        dict[version]: Version of web service
    """
    # TODO check if matrix is usable
    return {'version': __version__}


def wsgi_app(similarities, fragments, pharmacophores, external_url='http://localhost:8084/kripo'):
    """Create wsgi app

    Args:
        similarities (SimilarityMatrix): Similarity matrix to use in webservice
        fragments (FragmentsDb): Fragment database filename
        pharmacophores: Filename of pharmacophores hdf5 file
        external_url (str): URL which should be used in Swagger spec

    Returns:
        connexion.App
    """
    app = connexion.App(__name__)
    url = urlparse(external_url)
    swagger_file = resource_filename(__name__, 'swagger.yaml')
    app.app.json_encoder = KripodbJSONEncoder
    app.app.config['similarities'] = similarities
    app.app.config['fragments'] = fragments
    app.app.config['pharmacophores'] = pharmacophores
    arguments = {'hostport': url.netloc, 'scheme': url.scheme, 'version': __version__}
    # Keep validate_responses turned off, because of conflict with connexion.problem
    # see https://github.com/zalando/connexion/issues/266
    app.add_api(swagger_file, base_path=url.path, arguments=arguments)
    return app


def serve_app(similarities, fragments, pharmacophores, internal_port=8084, external_url='http://localhost:8084/kripo'):
    """Serve webservice forever

    Args:
        similarities: Filename of similarity matrix hdf5 file
        fragments: Filename of fragments database file
        pharmacophores: Filename of pharmacophores hdf5 file
        internal_port: TCP port on which to listen
        external_url (str): URL which should be used in Swagger spec
    """
    sim_matrix = open_similarity_matrix(similarities)
    pharmacophores_db = PharmacophoresDb(pharmacophores)
    app = wsgi_app(sim_matrix, fragments, pharmacophores_db, external_url)
    LOGGER.setLevel(logging.INFO)
    LOGGER.addHandler(logging.StreamHandler())
    LOGGER.info(' * Swagger spec at {}/swagger.json'.format(external_url))
    LOGGER.info(' * Swagger ui at {}/ui'.format(external_url))
    try:
        app.run(port=internal_port)
    finally:
        sim_matrix.close()
