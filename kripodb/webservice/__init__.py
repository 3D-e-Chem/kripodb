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
from pyramid.config import Configurator
from pyramid.view import view_config
from waitress import serve
from kripodb.version import __version__
from kripodb.hdf5 import DistanceMatrix


@view_config(route_name='fragments.similar', renderer='json')
def similar(request):
    distance_matrix = request.registry.settings['matrix']
    query_id = request.swagger_data['fragment_id']
    cutoff = request.swagger_data['cutoff']
    limit = None
    raw_hits = distance_matrix.find(query_id, cutoff, limit)
    # add query column
    hits = []
    for hit_id, score in raw_hits:
        hits.append({'query_frag_id': query_id, 'hit_frag_id': hit_id, 'score': score})
    return hits


@view_config(route_name='version', renderer='json')
def version(request):
    # TODO check if matrix is usable
    return {'version': __version__}


def wsgi_app(settings):
    # must be same as basePath in api_docs/swagger.json
    # and same as prefix in reverse proxy
    url_prefix = '/kripo'
    settings['pyramid_swagger.exclude_paths'] = [r'^/kripo/static', r'^/kripo/swagger.json']
    config = Configurator(settings=settings)
    config.include('pyramid_swagger', route_prefix=url_prefix)
    config.add_route('fragments.similar', url_prefix + '/fragments/{fragment_id}/similar')
    config.add_route('version', url_prefix + '/version')
    # Unpack dist directory from swagger ui (https://github.com/swagger-api/swagger-ui) into kripodb/static.
    config.add_static_view(url_prefix + '/static', path='kripodb:static')
    config.scan()
    return config.make_wsgi_app()


def serve_app(matrix, port, host):
    dist_matrix = DistanceMatrix(matrix)
    app = wsgi_app({'matrix': dist_matrix})
    print('Api docs on http://{host}:{port}/static/?url=http://{host}:{port}/kripo/swagger.json'.format(host=host, port=port))
    try:
        serve(app, host=host, port=port)
    finally:
        dist_matrix.close()
