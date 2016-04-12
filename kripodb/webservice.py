from kripodb.hdf5 import DistanceMatrix
from pyramid.config import Configurator
from waitress import serve


def similar(request):
    distance_matrix = request.registry.settings['matrix']
    query_id = request.swagger_data['fragment_id']
    cutoff = request.swagger_data['cutoff']
    limit = None
    raw_hits = distance_matrix.find(query_id, cutoff, limit)
    # add query column
    hits = []
    for hit_id, score in raw_hits:
        hits.append({'query_id': query_id, 'hit_id': hit_id, 'score': score})
    return hits


def wsgi_app(distancesfn, prefix):
    settings = {'matrix': DistanceMatrix(distancesfn)}
    config = Configurator(settings=settings)
    config.include('pyramid_swagger')
    config.add_route('fragments.similar', prefix + '/fragments/{fragment_id}/similar')
    config.add_view(similar, route_name='fragments.similar', renderer='json')
    return config.make_wsgi_app()


def serve_app(matrix, prefix, port, host):
    app = wsgi_app(matrix, prefix)
    serve(app, host=host, port=port)
