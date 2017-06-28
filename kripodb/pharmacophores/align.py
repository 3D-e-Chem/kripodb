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
from __future__ import absolute_import

import copy
from math import sqrt

import numpy as np
from rdkit.Numerics import Alignment
from scipy.spatial.distance import pdist, squareform

"""

Code was copied from Kripo yasara plugin written by Dave Wood, August 8th 2008
"""


def distances(points):
    array = np.array([p[1:] for p in points])
    return squareform(pdist(array)).tolist()


def transform_point(pt, data):
    x = data[0][0] * pt[0] + data[0][1] * pt[1] + data[0][2] * pt[2] + data[0][3]
    y = data[1][0] * pt[0] + data[1][1] * pt[1] + data[1][2] * pt[2] + data[1][3]
    z = data[2][0] * pt[0] + data[2][1] * pt[1] + data[2][2] * pt[2] + data[2][3]
    return x, y, z


class NoOverlapFound(ValueError):
    pass


class Aligner(object):
    def __init__(self, reference, probe):
        """Generates transformation matrix to align points of a reference pharmacophore to points of a probe pharmacophore

        Args:
            reference( list): Reference pharmacophore points. List of points where each point is (key,x,y,z)
            probe (list): Probe pharamcophore points. List of points where each point is (key,x,y,z)
        """
        self.reference = reference
        self.probe = probe

        # CONTAINERS FOR ALL POSSIBLE FEATURE PAIRS BETWEEN THE PHARMACOPHORES
        self.pairs = {}
        self.nodes = {}

        # CLIQUE DETECTION VARIABLES
        self.max_clique = 0
        self.clique_results = []
        self.complete = 0
        self.clique_count = 0

    def calculate_possible_pairs(self, cutoff=1.0):
        reference_distances = distances(self.reference)
        probe_distances = distances(self.probe)

        # IDENTIFY EACH PAIR IN PHARMACOPHORE A
        for a in range(len(self.reference) - 1):
            for b in range(a + 1, len(self.probe)):
                # RECORD THE FEATURE TYPES AND THE DISTANCE BETWEEN THE FEATURES
                a_type_a = self.reference[a][0]
                a_type_b = self.probe[b][0]
                a_distance = reference_distances[a][b]

                # CHECK WHETHER FEATURE PAIR IS SYMMETRICAL
                symetrical = a_type_a == a_type_b

                # IDENTIFY FIRST FEATURE FOR POSSIBLE PAIR IN PHARMACOPHORE B
                for c in range(len(self.probe) - 1):
                    # CHECK WHETHER THE FIRST FEATURE MATCHES EITHER OR BOTH
                    # OF THE PPHORE A FEATURE PAIR. CONTINUE IF NOT
                    b_type_a = self.probe[c][0]
                    if b_type_a == a_type_a:
                        match_a = 1
                    elif b_type_a == a_type_b:
                        match_a = 2
                    else:
                        continue

                    # IDENITFY SECOND FEATURE FOR PPHORE B FEATURE PAIR
                    for d in range(c + 1, len(self.probe)):
                        # CHECK WHETHER FEATURES MATCH THE FEATURE PAIR IN PHARMACOPHORE A
                        # IF NOT THEN CONTINUE
                        b_type_b = self.probe[d][0]
                        if not ((match_a == 1 and a_type_b == b_type_b) or (match_a == 2 and a_type_a == b_type_b)):
                            continue

                        # CALCULATE DISTANCE BETWEEN THE PPHORE B FEATURE PAIR
                        b_distance = probe_distances[c][d]

                        # CHECK WHETHER DISTANCE MATCHES
                        if abs(a_distance - b_distance) < cutoff:

                            # IF IT DOES, RECORD THE PAIR OF PAIRS AS BEING COMPATIBLE

                            # MAKE SURE YOU PAIR THE FEATURES THE RIGHT WAY ROUND
                            if match_a == 1:

                                # PPHORE A PAIR HAS NOT BEEN SEEN BEFORE MAKE A NEW ENTRY
                                if (a, b) not in self.pairs:
                                    self.pairs[a, b] = [(c, d)]
                                    self.pairs[b, a] = [(d, c)]

                                # OTHERWISE ADD PPHORE B PAIR TO THE COMPATIBLE PAIRS FOR PPHORE A PAIR
                                else:
                                    self.pairs[a, b].append((c, d))
                                    self.pairs[b, a].append((d, c))

                                # RECORD THE FEATURE PAIRS
                                self.nodes[a, c] = 0
                                self.nodes[b, d] = 0
                                if symetrical:
                                    self.nodes[a, d] = 0
                                    self.nodes[b, c] = 0

                            # FEATURES PAIRED THE OTHER WAY ROUND
                            elif match_a == 2:

                                # PPHORE A PAIR HAS NOT BEEN SEEN BEFORE MAKE A NEW ENTRY
                                if (a, b) not in self.pairs:
                                    self.pairs[a, b] = [(d, c)]
                                    self.pairs[b, a] = [(c, d)]

                                # PPHORE A PAIR HAS BEEN SEEN. ADD PPHORE B PAIR TO PPHORE A PAIR LIST
                                else:
                                    self.pairs[a, b].append((d, c))
                                    self.pairs[b, a].append((c, d))

                                # RECORD THE FEATURE PAIRS
                                self.nodes[a, d] = 0
                                self.nodes[b, c] = 0

    def _detect_cliques(self, r, p, x, break_num_cliques=3000, complete=False):
        """BronKerbosch clique detection algorithm
        Based on the pseudo code in 'Reporting maximal cliques: new
        insights into an old problem' Cazaals and Chinmay 2007"""

        # IF SATISFACTORY CLIQUE HAS BEEN FOUND THEN BREAK OUT OF THE LOOP
        if self.complete:
            return

        # SET COMPLETION IF BREAK POINT REACHED
        elif break_num_cliques != 0 and self.clique_count > break_num_cliques:
            self.complete = 1
            return

        # COUNT THE NUMBER OF CLIQUES RECORDED SO FAR
        elif len(p) == 0 and len(x) == 0:
            self.clique_count += 1

            # RECORD CLIQUE IF IT IS THE LARGEST FOUND SO FAR
            if len(r) > self.max_clique:
                self.clique_results.append(copy.copy(r))
                self.max_clique = len(r)

                # IF CLIQUE IS LARGER THAN HALF THE SIZE OF THE SMALLEST PHARMACOPHORE SET COMPLETE
                if (not complete) and len(r) >= (min(len(self.reference), len(self.probe)) / 2):
                    self.complete = 1

        # DO THE BRON KERBOSCH THING
        else:

            # FOR EACH POSSIBLE PAIR OF FEATURES
            i = 0
            while i < len(p):
                u_i = p.pop(i)

                # MAKE COPY OF MATCHED PAIRS
                r_new = copy.copy(r)
                r_new.append(u_i)

                # MAKE COPY OF UNTESTED PAIRS
                p_new = copy.copy(p)

                # CHECK WHETHER UNTESTED FEATURE PAIRS ARE COMPATIBLE WITH CURRENT CLIQUE
                j = 0
                while j < len(p_new):
                    p_i = p_new[j]
                    if ((u_i[0], p_i[0]) not in self.pairs) or not self.pairs[u_i[0], p_i[0]].count(
                            (u_i[1], p_i[1])):
                        p_new.pop(j)
                        j -= 1
                    j += 1

                # CHECK WHETHER NEW CLIQUE IS PART OF PREVIOUSLY TESTED CLIQUE
                x_new = copy.copy(x)
                j = 0
                while j < len(x_new):
                    x_i = x_new[j]
                    if ((u_i[0], x_i[0]) not in self.pairs) or not self.pairs[u_i[0], x_i[0]].count(
                            (u_i[1], x_i[1])):
                        x_new.pop(j)
                        j -= 1
                    j += 1

                # PROGRESS A LEVEL IN THE SEARCH TREE
                self._detect_cliques(r_new, p_new, x_new, break_num_cliques, complete)

                # IF THE MAX CLIQUE HAS BEEN FOUND THEN BACKTRACK THROUGH THE SEARCH TREE
                if self.complete:
                    break

                # ADD FEATURE PAIR TO CURRENT CLIQUE
                x.append((u_i[0], u_i[1]))

    def cliques(self, break_num_cliques=3000):
        # BRON KERBOSCH CLIQUE DETECTION BETWEEN PHARMACOPHORES
        r = []
        p = list(self.nodes.keys())
        x = []

        self.clique_results = []
        self._detect_cliques(r, p, x, break_num_cliques)
        return self.clique_results

    def cliqued_pharmacophores(self):
        # PREPARE THE CLIQUE PHARMACOPHORES
        clique_reference = []
        clique_probe = []
        if not self.clique_results:
            raise NoOverlapFound()
        for pair in self.clique_results[-1]:
            clique_reference.append(self.reference[int(pair[0])])
            clique_probe.append(self.probe[int(pair[1])])
        return clique_reference, clique_probe

    def transformation(self, cutoff=1.0, break_num_cliques=3000):
        """Calculate the transformation matrix to align probe pharmacophore to reference pharmacophore.

        Args:
            cutoff (float): Tolerance threshold for considering two distances to be equivalent
            break_num_cliques (int): Break when {break_num_cliques} cliques found, default 3000.
                Set to zero (0) for a complete search.

        Returns:
            tuple: First item is the RSMD value for the alignment and
                    second item is the 4x4 transform matrix

        """
        self.calculate_possible_pairs(cutoff)
        self.cliques(break_num_cliques)

        (reference, probe) = self.cliqued_pharmacophores()
        reference_points = [r[1:] for r in reference]
        probe_points = [r[1:] for r in probe]
        ssd, matrix = Alignment.GetAlignmentTransform(reference_points, probe_points)
        rmsd = sqrt(ssd / len(reference_points))
        return rmsd, matrix


def align(pharmacophore, matrix):
    """Applies alignment transformation matrix on pharmacophore

    Args:
        pharmacophore (list): List of points where each point is (key,x,y,z)
        matrix (list): 4x4 transform matrix

    Returns:
        list: List of aligned points where each point is (key,x,y,z)
    """
    aligned_probe = []
    for d in pharmacophore:
        p = transform_point(d[1:], matrix)
        aligned_probe.append((
            d[0],
            p[0],
            p[1],
            p[2],
        ))
    return aligned_probe
