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
from six import StringIO, BytesIO

from mock import patch
import pytest

from kripodb.pdb import PdbReport


@pytest.fixture
def mock_fetch_response():
    mresponse = BytesIO()
    mresponse.write(b'structureId,chainId,structureTitle,compound,ecNo,uniprotAcc,uniprotRecommendedName\n')
    mresponse.write(
        b'"104L","B","HOW AMINO-ACID INSERTIONS ARE ALLOWED IN AN ALPHA-HELIX OF T4 LYSOZYME","T4 LYSOZYME","3.2.1.17","P00720","Endolysin"\n')
    mresponse.write(b'"12E8","H","2E8 FAB FRAGMENT","IGG1-KAPPA 2E8 FAB (HEAVY CHAIN)","","",""\n')
    mresponse.seek(0)
    return mresponse


class TestPdbReport(object):

    def test_url_default(self):
        pdb_report = PdbReport()

        url = pdb_report.url

        expected = 'http://www.rcsb.org/pdb/rest/customReport?' \
                   'pdbids=*&' \
                   'customReportColumns=structureTitle,compound,ecNo,uniprotAcc,uniprotRecommendedName&' \
                   'format=csv&service=wsfile'
        assert url == expected

    def test_url_custom(self):
        pdbids = ['1kvm', '2mbs']
        fields = ['resolution']
        pdb_report = PdbReport(pdbids, fields)

        url = pdb_report.url

        expected = 'http://www.rcsb.org/pdb/rest/customReport?' \
                   'pdbids=1kvm,2mbs&' \
                   'customReportColumns=resolution&' \
                   'format=csv&service=wsfile'
        assert url == expected

    @patch('kripodb.pdb.urlopen')
    def test_fetch(self, mocked_urlopen, mock_fetch_response):
        mocked_urlopen.return_value = mock_fetch_response

        pdb_report = PdbReport(['104L', '12E8'])

        pdbs = list(pdb_report.fetch())

        expected = [{
            'chainId': 'B',
            'structureId': '104L',
            'structureTitle': 'HOW AMINO-ACID INSERTIONS ARE ALLOWED IN AN ALPHA-HELIX OF T4 LYSOZYME',
            'ecNo': '3.2.1.17',
            'uniprotAcc': 'P00720',
            'compound': 'T4 LYSOZYME',
            'uniprotRecommendedName': 'Endolysin'
        }, {
            'chainId': 'H',
            'structureId': '12E8',
            'structureTitle': '2E8 FAB FRAGMENT',
            'ecNo': None,
            'uniprotAcc': None,
            'compound': 'IGG1-KAPPA 2E8 FAB (HEAVY CHAIN)',
            'uniprotRecommendedName': None
        }]
        assert pdbs == expected




