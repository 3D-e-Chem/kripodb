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
from pkg_resources import resource_filename
import logging
from six.moves.urllib_parse import urlparse


import connexion
from flask import current_app, abort

from ..version import __version__
from ..pairs import open_distance_matrix

LOGGER = logging.getLogger(__name__)


def get_similar_fragments(fragment_id, cutoff, limit):
    """Find similar fragments to query.

    Args:
        fragment_id (str): Query fragment identifier
        cutoff (float): Cutoff, distance scores below cutoff are discarded.
        limit (int): Maximum number of hits. Default is None for no limit.

    Returns:
        List(Dict()): Query fragment identifier, hit fragment identifier and distance score
    """
    distance_matrix = current_app.config['matrix']
    query_id = fragment_id
    hits = []
    try:
        raw_hits = distance_matrix.find(query_id, cutoff, limit)
        # add query column
        for hit_id, score in raw_hits:
            hits.append({'query_frag_id': query_id, 'hit_frag_id': hit_id, 'score': score})
    except KeyError:
        abort(404)
    return hits


def get_version():
    """
    Returns:
        Dict: Version of web service
    """
    # TODO check if matrix is usable
    return {'version': __version__}


def wsgi_app(dist_matrix, external_url='http://localhost:8084/kripo'):
    """Create wsgi app

    Args:
        dist_matrix (DistanceMatrix): Distance matrix to use in webservice
        external_url (str): URL which should be used in Swagger spec

    Returns:
        connexion.App
    """
    app = connexion.App(__name__)
    url = urlparse(external_url)
    swagger_file = resource_filename(__name__, 'swagger.json')
    app.add_api(swagger_file, base_path=url.path, arguments={'hostport': url.netloc, 'scheme': url.scheme})
    app.app.config['matrix'] = dist_matrix
    return app


def serve_app(matrix, internal_port=8084, external_url='http://localhost:8084/kripo'):
    """Serve webservice forever

    Args:
        matrix: Filename of distance matrix hdf5 file
        internal_port: TCP port on which to listen
        external_url (str): URL which should be used in Swagger spec
    """
    dist_matrix = open_distance_matrix(matrix)
    app = wsgi_app(dist_matrix, external_url)
    LOGGER.setLevel(logging.INFO)
    LOGGER.addHandler(logging.StreamHandler())
    LOGGER.info(' * Swagger spec at {}/swagger.json'.format(external_url))
    LOGGER.info(' * Swagger ui at {}/ui'.format(external_url))
    try:
        app.run(port=internal_port)
    finally:
        dist_matrix.close()
