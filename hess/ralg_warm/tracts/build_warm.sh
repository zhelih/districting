#!/bin/bash

set -x

states="
AK
AL
AR
AZ
CA
CO
CT
DE
FL
GA
HI
IA
ID
IL
IN
KS
KY
LA
MA
MD
ME
MI
MN
MO
MS
MT
NC
ND
NE
NH
NJ
NM
NV
NY
OH
OK
OR
PA
RI
SC
SD
TN
TX
UT
VA
VT
WA
WI
WV
WY
"

for state in $states; do
  echo "Running $state"
  /home/lykhovyd/progs/districting/hess/districting myconfig.txt $state
done
