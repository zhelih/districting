#!/bin/bash

echo "Generating config for Grid $1, model $2"

CFG_NAME="grid_$1.config"

echo "Writing down $CFG_NAME..."

echo "dimacs grid_$1_$1.dimacs" > $CFG_NAME # deletes
echo "distance grid_$1_$1_distances.csv" >> $CFG_NAME
echo "population grid_$1_$1.population" >> $CFG_NAME
echo "L $1" >> $CFG_NAME
echo "U $1" >> $CFG_NAME
echo "k $1" >> $CFG_NAME
echo "model $2" >> $CFG_NAME
echo "output grid_output.csv" >> $CFG_NAME
echo "Done"
