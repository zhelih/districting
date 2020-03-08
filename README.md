# Imposing contiguity constraints in political districting models
C++ code and experimental results to accompany the paper "Imposing contiguity constraints in political districting models" by Hamidreza Validi, Austin Buchanan, and Eugene Lykhovyd. [(pdf)](http://www.optimization-online.org/DB_HTML/2020/01/7582.html) [(slides)](https://github.com/zhelih/districting/blob/master/Districting_slides.pdf) :earth_americas:

#### Build:
```
cd src/
make
```
Environment variable `GUROBI_HOME` must point to current GUROBI installation. The current `Makefile` is written for version 8.1, corresponding changes might be needed to correctly link with another version.
#### Binaries

- `ralg_hot_start` computes good starting point for the r-algorithm, e.g., computes Lagrangian Dual bound. This is important step to fix as many variables as possible.

- `districting` main binary: computes Lagrangian Dual, heuristic, fixes variables and finds the districting partition.

- `translate` converts results of districting to GEO mapping

- `sol_to_png.py` converts GEO mapping to .png using QGIS


**TODO more details about binaries and format**
#### Config format
```
# ####### config file for "districting"
# the database from online source, see the Link in the end of this README
database /path/to/db
# level must say tracts or counties
level change_to_tracts_or_counties
# Two-letter state identifier. Optional, can be passed using cmd arguments.
state AA
# If need to run for specific graph files and not from the database. Alternative configuration.
# Remove if db is used.
dimacs /path/to/dimacs
distance /path/to/dist
population /path/to/pop
# L,U,k - interger parameters for the mode. Use auto if using the db. Can be any number.
L 10
U auto
k auto
# see available models running ./districting
model hess
# Optional hot start for r-algorithm. Can be passed with cmd arguments.
ralg_hot_start /path/to/file
# Resulting CSV file. Appends comma-separated computational results
output /path/to/output.csv
```
#### Online Database!
Districting database can be found [here](https://lykhovyd.com/files/public/districting).
