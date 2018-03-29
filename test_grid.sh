#!/bin/bash

for i in `seq 2 15`;
do
  res=$(./gerry_assign grid $i)
  obj=$(echo "$res" | grep Obj | tail -n 1 | awk '{ print $2}')
  tim=$(echo "$res" | grep Explored | awk '{ print $8 }')
  laz=$(echo "$res" | grep Lazy | awk '{ print $3 }')
  echo "$i, $obj, $tim, $laz"
#  sleep 0.5 # w8 for shared resources
done

wait
