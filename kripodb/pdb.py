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
import codecs
import csv
from six.moves import zip
from six.moves.urllib.request import urlopen


def parse_csv_file(thefile):
    """Parse csv file, yielding rows as dictionary.

    The csv file should have an header.

    Args:
        thefile (file): File like object

    Yields:
         dict: Dictionary with column header name as key and cell as value

    """
    reader = csv.reader(codecs.iterdecode(thefile, 'ISO-8859-1'))
    # read header
    colnames = next(reader)
    # data rows
    for row in reader:
        pdb = {}
        for k, v in zip(colnames, row):
            if v is '':
                v = None
            pdb[k] = v

        yield pdb


class PdbReport(object):
    """Client for the Custom Report Web Services of the RCSB PDB website

    See http://www.rcsb.org/pdb/software/wsreport.do for more information.

    Args:
        pdbids (List[str]): List of pdb identifiers to fetch. Default is ['*'] which fetches all.
        fields: (List[str]: List of fields to fetch.
            Default is ['structureTitle', 'compound', 'ecNo', 'uniprotAcc', 'uniprotRecommendedName']
            See http://www.rcsb.org/pdb/results/reportField.do for possible fields.

    Attributes:
        url (str): Url of report, based on pdbids and fields.

    """
    url_tpl = 'http://www.rcsb.org/pdb/rest/customReport?pdbids={pdbids}' \
              '&customReportColumns={fields}&format=csv&service=wsfile'

    def __init__(self, pdbids=None, fields=None):
        self.pdbids = ['*']
        self.fields = ['structureTitle', 'compound', 'ecNo', 'uniprotAcc', 'uniprotRecommendedName']
        if pdbids is not None:
            self.pdbids = pdbids
        if fields is not None:
            self.fields = fields

    @property
    def url(self):
        return self.url_tpl.format(pdbids=','.join(self.pdbids), fields=','.join(self.fields))

    def fetch(self):
        """Fetch report from PDB website

        Yields:
            dict: Dictionary with keys same as ['structureId', 'chainID'] + self.fields

        """
        response = urlopen(self.url)
        return parse_csv_file(response)
