This is the decription of the data files (example for texas state)

All census tracts are indexed by GEOID (later should be translated into vertex labeling)


All data is contained in the folder texas/

**Shapfile of centroids for census tracts**

texas_centers.prj  texas_centers.shp 
texas_centers.qpj  texas_centers.shx
texas_centers.dbf

**Shapefile of census tracts as polygons**

texas_c.dbf  texas_c.prj  texas_c.shp
texas_c.qpj  texas_c.shx

**Distance matrix between centroids as CSV table (huge)**

texas_d.csv
texas_d.zip (zipped version of .csv)

**Data table listing population and neigbors**

texas_features.xls 
texas_features.csv

**QGIS project files**

texas.qgs
texas.qgs~

** Graph generated file **

texas.dimacs
texas.hash
texas.population


Algorithm to craete these files:
1) open QGIS
2) load Profile-County GDB
3) filter layer: right-click in layers-panel, then filter, expression: "GEOID10" LIKE '45%'
4) save as shapefile: right-click in layers-panel, save-as
5) centroids: menu vector -> geometry tools -> polygon centroids
6) distance matrix: menu vector -> analysis -> distance matrix : choose full
7) generate features: plugins -> python console -> load script -> run (don't forget to select correct active layer)
8) run a programming tool to get .dimacs, .hash and .population
9) get a glass of beer to enjoy the evening
