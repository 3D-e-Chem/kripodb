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

import csv


def read_id2label(filename):
    with open(filename) as f:
        return read(f)


def read(infile):
    id2label = {}
    reader = csv.reader(infile, delimiter="\t")
    for row in reader:
        id2label[int(row[0])] = row[1]
    return id2label


def swap_label2id(id2label):
    return {v: k for k, v in id2label.iteritems()}
