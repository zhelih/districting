# districting
C++ code and experimental results to accompany the paper "Imposing contiguity constraints in political districting models" by Hamidreza Validi, Austin Buchanan, and Eugene Lykhovyd. :earth_americas:

#### Build:
```
cd src/
make
```
#### Binaries

- `ralg_hot_start` computes good starting point for the r-algorithm, which compute Lagrangian Dual bound. This is important step to fix as much variables as possible.

- `districting` main binary: compute Lagrangian Dual, heuristic, fixes variables and find distrcting partition.

- `translate` converts result of districting to GEO mapping

- `sol_to_png.py` converts GEO mapping to .png using QGIS


**TODO more details about binaries and format**
**TODO config format description**

Districting database can be found [here](https://lykhovyd.com/files/public/districting).
