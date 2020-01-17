#!/bin/bash

for nTrees in 1 10 100 1000 10000
  do
  for treeDepth in 1 2 3 4 5
  do
    #echo $nTrees $treeDepth
    python scripts/TrainingPandoraBdt.py $treeDepth $nTrees
  done
done
