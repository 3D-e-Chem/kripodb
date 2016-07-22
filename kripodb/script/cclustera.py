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
from __future__ import absolute_import

from ..db import FragmentsDb


def cclustera_sphere(inputfile, outputfile, onlyfrag1):
    import json, math
    frags_db = FragmentsDb(inputfile)
    nodes = {}

    # distribute fragments evenly on sphere using Fibonacci sphere algorithm
    # from http://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
    samples = len(frags_db)

    sql = 'SELECT frag_id, pdb_code, het_code FROM fragments'
    if onlyfrag1:
        sql += ' WHERE frag_id LIKE "%_frag1"'
        frags_db.cursor.execute('SELECT count(*) FROM fragments WHERE frag_id LIKE "%_frag1"')
        samples = frags_db.cursor.fetchone()[0]

    rnd = 1.
    offset = 2. / samples
    increment = math.pi * (3. - math.sqrt(5.));

    frag_ids = frags_db.cursor.execute(sql)
    for i, frag in enumerate(frag_ids):
        y = ((i * offset) - 1) + (offset / 2);
        r = math.sqrt(1 - pow(y, 2))

        phi = ((i + rnd) % samples) * increment

        x = math.cos(phi) * r
        z = math.sin(phi) * r

        node_info = {
            'Path': [],
            'Coordinates': [x, y, z],
            'Categories': [frag[1], frag[2]],
            'Properties': []
        }
        nodes[frag[0]] = node_info

    json.dump(nodes, outputfile)