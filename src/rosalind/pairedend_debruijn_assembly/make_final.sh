#!/bin/bash

SOURCES="paired_kmer.py graph.py path.py string_builder.py cycle_finder.py cycle_generator.py path_finder.py main.py"

TARGET=pairedend_debruijn_assembly.py

[[ -f $TARGET ]] && rm $TARGET

for f in $SOURCES
do
  (cat "$f"; echo '') >> $TARGET
done
