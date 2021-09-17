#NB: run from QGis console

# display internal processing algorithms
for alg in QgsApplication.processingRegistry().algorithms():
    print(alg.id(), "->", alg.displayName())

# display algorithm args and how to run it
processing.algorithmHelp("native:centroids")
processing.algorithmHelp("qgis:distancematrix")
