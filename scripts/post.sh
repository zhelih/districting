#!/bin/bash

folder=$1
abbr=$2

rootdir="/home/zhelih/html/files/public/gerry/"

if [ -z "$folder" ]
then
  echo "Specify folder as a first argument"
  exit 1
fi

if [ -z "$abbr" ]
then
  echo "Specify abbreviation as a second argument"
  exit 1
fi

echo "Postprocess for $folder ($abbr) and uploading to mango"
echo -n "Creating folders..."

ssh mango "mkdir $rootdir$abbr"
ssh mango "mkdir $rootdir$abbr/counties"
ssh mango "mkdir $rootdir$abbr/counties/graph"
ssh mango "mkdir $rootdir$abbr/counties/maps"
ssh mango "mkdir $rootdir$abbr/census_tracts"
ssh mango "mkdir $rootdir$abbr/census_tracts/graph"
ssh mango "mkdir $rootdir$abbr/census_tracts/maps"

echo "done"

## Comput_all.native
echo -n "Running counties compute_all.native..."
./compute_all.native "$folder/counties/$abbr"
echo -n "Running census_tracts compute_all.native..."
./compute_all.native "$folder/$abbr"

## Distances
echo -n "Sorting and truncating the counties distance file..."
./main.py "$folder/counties/$abbr.hash" "$folder/counties/${abbr}_d.csv" "$folder/counties/${abbr}_distances.csv"
echo "done"

echo -n "Sorting and truncating the census tracts distance file..."
./main.py "$folder/$abbr.hash" "$folder/${abbr}_d.csv" "$folder/${abbr}_distances.csv"
echo "done"

## Sorting pop and hash
echo -n "Sorting pop and hash files..."
sort -n "$folder/counties/$abbr.hash" > "$folder/counties/$abbr.hash_sorted"
mv "$folder/counties/$abbr.hash_sorted" "$folder/counties/$abbr.hash"
sort -n "$folder/counties/$abbr.population" > "$folder/counties/$abbr.pop_sorted"
mv "$folder/counties/$abbr.pop_sorted" "$folder/counties/$abbr.population"
./swap.py "$folder/counties/$abbr.population"

sort -n "$folder/$abbr.hash" > "$folder/$abbr.hash_sorted"
mv "$folder/$abbr.hash_sorted" "$folder/$abbr.hash"
sort -n "$folder/$abbr.population" > "$folder/$abbr.pop_sorted"
mv "$folder/$abbr.pop_sorted" "$folder/$abbr.population"
./swap.py "$folder/$abbr.population"

echo "done"

## Upload maps
# dbf prj qpj shp shx
# _tracts _centers _centers_utm
echo -n "Uploading maps..."
scp "$folder/counties/${abbr}_counties.dbf" mango:$rootdir/$abbr/counties/maps/
scp "$folder/counties/${abbr}_counties.prj" mango:$rootdir/$abbr/counties/maps/
scp "$folder/counties/${abbr}_counties.qpj" mango:$rootdir/$abbr/counties/maps/
scp "$folder/counties/${abbr}_counties.shp" mango:$rootdir/$abbr/counties/maps/
scp "$folder/counties/${abbr}_counties.shx" mango:$rootdir/$abbr/counties/maps/

scp "$folder/counties/${abbr}_centers.dbf" mango:$rootdir/$abbr/counties/maps/
scp "$folder/counties/${abbr}_centers.prj" mango:$rootdir/$abbr/counties/maps/
scp "$folder/counties/${abbr}_centers.qpj" mango:$rootdir/$abbr/counties/maps/
scp "$folder/counties/${abbr}_centers.shp" mango:$rootdir/$abbr/counties/maps/
scp "$folder/counties/${abbr}_centers.shx" mango:$rootdir/$abbr/counties/maps/

scp "$folder/counties/${abbr}_centers_proj.dbf" mango:$rootdir/$abbr/counties/maps/
scp "$folder/counties/${abbr}_centers_proj.prj" mango:$rootdir/$abbr/counties/maps/
scp "$folder/counties/${abbr}_centers_proj.qpj" mango:$rootdir/$abbr/counties/maps/
scp "$folder/counties/${abbr}_centers_proj.shp" mango:$rootdir/$abbr/counties/maps/
scp "$folder/counties/${abbr}_centers_proj.shx" mango:$rootdir/$abbr/counties/maps/

# for tracts
scp "$folder/${abbr}_tracts.dbf" mango:$rootdir/$abbr/census_tracts/maps/
scp "$folder/${abbr}_tracts.prj" mango:$rootdir/$abbr/census_tracts/maps/
scp "$folder/${abbr}_tracts.qpj" mango:$rootdir/$abbr/census_tracts/maps/
scp "$folder/${abbr}_tracts.shp" mango:$rootdir/$abbr/census_tracts/maps/
scp "$folder/${abbr}_tracts.shx" mango:$rootdir/$abbr/census_tracts/maps/
scp "$folder/${abbr}_tracts.cpg" mango:$rootdir/$abbr/census_tracts/maps/

scp "$folder/${abbr}_centers.dbf" mango:$rootdir/$abbr/census_tracts/maps/
scp "$folder/${abbr}_centers.prj" mango:$rootdir/$abbr/census_tracts/maps/
scp "$folder/${abbr}_centers.qpj" mango:$rootdir/$abbr/census_tracts/maps/
scp "$folder/${abbr}_centers.shp" mango:$rootdir/$abbr/census_tracts/maps/
scp "$folder/${abbr}_centers.shx" mango:$rootdir/$abbr/census_tracts/maps/

scp "$folder/${abbr}_centers_proj.dbf" mango:$rootdir/$abbr/census_tracts/maps/
scp "$folder/${abbr}_centers_proj.prj" mango:$rootdir/$abbr/census_tracts/maps/
scp "$folder/${abbr}_centers_proj.qpj" mango:$rootdir/$abbr/census_tracts/maps/
scp "$folder/${abbr}_centers_proj.shp" mango:$rootdir/$abbr/census_tracts/maps/
scp "$folder/${abbr}_centers_proj.shx" mango:$rootdir/$abbr/census_tracts/maps/

echo "done"

## Upload graph data
echo -n "Uploading graph data..."
scp "$folder/counties/${abbr}.dimacs" mango:$rootdir/$abbr/counties/graph/
scp "$folder/counties/${abbr}.hash" mango:$rootdir/$abbr/counties/graph/
scp "$folder/counties/${abbr}.population" mango:$rootdir/$abbr/counties/graph/
scp "$folder/counties/${abbr}_distances.csv" mango:$rootdir/$abbr/counties/graph/

scp "$folder/${abbr}.dimacs" mango:$rootdir/$abbr/census_tracts/graph/
scp "$folder/${abbr}.hash" mango:$rootdir/$abbr/census_tracts/graph/
scp "$folder/${abbr}.population" mango:$rootdir/$abbr/census_tracts/graph/
scp "$folder/${abbr}_distances.csv" mango:$rootdir/$abbr/census_tracts/graph/

echo "done"

# upload features
scp "$folder/counties/${abbr}_features.csv" mango:$rootdir/$abbr/counties/
scp "$folder/${abbr}_features.csv" mango:$rootdir/$abbr/census_tracts/

echo "All done!"
