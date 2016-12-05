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

import argparse
import sys

from .fingerprints import make_fingerprints_parser
from .fragments import make_fragments_parser
from .similarities import make_similarities_parser
from .dive import make_dive_parsers
from kripodb.version import __version__


def make_parser():
    """Creates a parser with sub commands

    Returns:
        argparse.ArgumentParser: parser with sub commands
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version', version=__version__)
    subparsers = parser.add_subparsers(dest='subcommand')

    make_fingerprints_parser(subparsers)

    make_fragments_parser(subparsers)

    make_similarities_parser(subparsers)

    make_dive_parsers(subparsers)

    return parser


def main(argv=sys.argv[1:]):
    """Main script function.

    Calls run method of selected sub commandos.

    Args:
        argv (list[str]): List of command line arguments

    """
    parser = make_parser()
    args = parser.parse_args(argv)
    fargs = vars(args)
    if 'func' in fargs:
        func = args.func
        del(fargs['subcommand'])
        del(fargs['func'])
        func(**fargs)
    else:
        if 'subcommand' in args:
            parser.parse_args([args.subcommand, '--help'])
        else:
            parser.print_help()
