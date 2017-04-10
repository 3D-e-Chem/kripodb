#!/usr/bin/env bash

cp ../current/fragments.sqlite .
kripodb fragments shelve fragments.shelve fragments.sqlite
kripodb fragments sdf fragments.sd fragments.sqlite
kripodb fragments pdb fragments.sqlite