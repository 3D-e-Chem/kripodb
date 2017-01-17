# coding=utf-8
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
"""Module to calculate modified tanimoto similarity"""

from __future__ import absolute_import
from math import fsum
import six


def calc_mean_onbit_density(bitsets, number_of_bits):
    """Calculate the mean density of bits that are on in bitsets collection.

    Args:
        bitsets (list[intbitset.intbitset]): List of fingerprints
        number_of_bits: Number of bits for all fingerprints

    Returns:
        float: Mean on bit density

    """
    all_nr_onbits = [len(v) for v in bitsets]
    mean_onbit = fsum(all_nr_onbits) / float(len(all_nr_onbits))
    density = mean_onbit / number_of_bits
    return float(density)


def corrections(mean_onbit_density):
    """Calculate corrections

    See :func:`similarity` for explanation of corrections.

    Args:
        mean_onbit_density (float): Mean on bit density

    Returns:
        float: S\ :sub:`T` correction, S\ :sub:`T0` correction
    """
    p0 = mean_onbit_density
    corr_st = (2 - p0) / 3
    corr_sto = (1 + p0) / 3
    return corr_st, corr_sto


def similarity(bitset1, bitset2, number_of_bits, corr_st, corr_sto):
    """Calculate modified Tanimoto similarity between two fingerprints

    Given two fingerprints of length n with a and b bits set in each fingerprint,
    respectively, and c bits set in both fingerprint,
    selected from a data set of fingerprint with a mean bit density of ρ\ :sub:`0`,
    the modified Tanimoto similarity S\ :sub:`MT` is calculated as

    .. math::

        S_{MT} = (\\frac{2 - ρ_0}{3}) S_T + (\\frac{1 + ρ_0}{3}) S_{T0}

    where ST is the standard Tanimoto coefficient

    .. math::

        S_T = \\frac{c}{a + b - c}

    and Sr0 is the inverted Tanimoto coefficient

    .. math::

        S_{T0} = \\frac{n - a - b - c}{n -c}

    Args:
        bitset1 (intbitset.intbitset): First fingerprint
        bitset2 (intbitset.intbitset): Second fingerprint
        number_of_bits (int): Number of bits for all fingerprints
        corr_st (float): St correction
        corr_sto (float): Sto correction

    Returns:
        float: modified Tanimoto similarity
    """
    a = len(bitset1)
    b = len(bitset2)
    c = len(bitset1 & bitset2)
    n = number_of_bits
    st = float(c) / (a + b - c)
    st0 = (n - a - b - + c) / float(n - c)
    smt = corr_st * st + corr_sto * st0
    return smt


def similarities(bitsets1, bitsets2, number_of_bits, corr_st, corr_sto, cutoff, ignore_upper_triangle=False):
    """Calculate modified tanimoto similarity between two collections of fingerprints

    Excludes similarity of the same fingerprint.

    Args:
        bitsets1 (Dict{str, intbitset.intbitset}): First dict of fingerprints
            with fingerprint label as key and intbitset as value
        bitsets2 (Dict{str, intbitset.intbitset}): Second dict of fingerprints
            with fingerprint label as key and intbitset as value
        number_of_bits (int): Number of bits for all fingerprints
        corr_st (float): St correction
        corr_sto (float): Sto correction
        cutoff (float): Cutoff, similarity scores below cutoff are discarded.
        ignore_upper_triangle (Optional[bool]): When true returns similarity where label1 > label2,
            when false returns all similarities

    Yields:
        (fingerprint label 1, fingerprint label2, similarity score)

    """
    for (label1, bs1) in six.iteritems(bitsets1):
        for (label2, bs2) in six.iteritems(bitsets2):
            if label1 == label2:
                # always skip self
                continue
            if ignore_upper_triangle and label1 > label2:
                continue

            score = similarity(bs1, bs2, number_of_bits, corr_st, corr_sto)

            if score >= cutoff:
                yield label1, label2, score
