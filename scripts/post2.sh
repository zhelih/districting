#!/bin/bash

abbr=$1

rootdir="XXXXX/src/data/districting_census2020/"

if [ -z "$abbr" ]
then
  echo "Specify abbreviation as a first argument"
  exit 1
fi

echo "Postprocess for ($abbr)"
echo -n "Creating folders..."

mkdir -p $rootdir$abbr
mkdir -p $rootdir$abbr/counties
mkdir -p $rootdir$abbr/counties/graph
mkdir -p $rootdir$abbr/counties/maps
mkdir -p $rootdir$abbr/census_tracts
mkdir -p $rootdir$abbr/census_tracts/graph
mkdir -p $rootdir$abbr/census_tracts/maps

echo "done"

echo -n "Moving files..."
mv $rootdir$abbr/counties/maps/*.csv $rootdir$abbr/counties/graph/
mv $rootdir$abbr/census_tracts/maps/*.csv $rootdir$abbr/census_tracts/graph/
echo "done"

## Compute_all.native
echo -n "Running counties compute_all.native..."
./compute_all.native "$rootdir$abbr/counties/graph/$abbr"
echo -n "Running census_tracts compute_all.native..."
./compute_all.native "$rootdir$abbr/census_tracts/graph/$abbr"

## Distances
echo -n "Sorting and truncating the counties distance file..."
./main.py "$rootdir$abbr/counties/graph/$abbr.hash" "$rootdir$abbr/counties/graph/${abbr}_d.csv" "$rootdir$abbr/counties/graph/${abbr}_distances.csv"
echo "done"

echo -n "Sorting and truncating the census tracts distance file..."
./main.py "$rootdir$abbr/census_tracts/graph/$abbr.hash" "$rootdir$abbr/census_tracts/graph/${abbr}_d.csv" "$rootdir$abbr/census_tracts/graph/${abbr}_distances.csv"
echo "done"

## Sorting pop and hash
echo -n "Sorting pop and hash files..."
sort -n "$rootdir$abbr/counties/graph/$abbr.hash" > "$rootdir$abbr/counties/graph/$abbr.hash_sorted"
mv "$rootdir$abbr/counties/graph/$abbr.hash_sorted" "$rootdir$abbr/counties/graph/$abbr.hash"
sort -n "$rootdir$abbr/counties/graph/$abbr.population" > "$rootdir$abbr/counties/graph/$abbr.pop_sorted"
mv "$rootdir$abbr/counties/graph/$abbr.pop_sorted" "$rootdir$abbr/counties/graph/$abbr.population"
./swap.py "$rootdir$abbr/counties/graph/$abbr.population"

sort -n "$rootdir$abbr/census_tracts/graph/$abbr.hash" > "$rootdir$abbr/census_tracts/graph/$abbr.hash_sorted"
mv "$rootdir$abbr/census_tracts/graph/$abbr.hash_sorted" "$rootdir$abbr/census_tracts/graph/$abbr.hash"
sort -n "$rootdir$abbr/census_tracts/graph/$abbr.population" > "$rootdir$abbr/census_tracts/graph/$abbr.pop_sorted"
mv "$rootdir$abbr/census_tracts/graph/$abbr.pop_sorted" "$rootdir$abbr/census_tracts/graph/$abbr.population"
./swap.py "$rootdir$abbr/census_tracts/graph/$abbr.population"
echo "done"

echo -n "Moving features file..."
mv "$rootdir$abbr/counties/graph/${abbr}_features.csv" "$rootdir$abbr/counties/${abbr}_features.csv"
mv "$rootdir$abbr/census_tracts/graph/${abbr}_features.csv" "$rootdir$abbr/census_tracts/${abbr}_features.csv"
echo "done"

echo "All done!"
