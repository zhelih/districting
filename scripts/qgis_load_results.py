#  Copyright 2018 Eugene Lykhovyd
#  Select the shapefile layer
#  Run the script
#  Go to layer properties, select categorized fill, then press classify
from qgis.utils import iface
from PyQt4.QtCore import QVariant

import random

NEW_DISTRICT_FIELD = "district"
DISTRICT_FILE = "/home/lykhovyd/progs/gerry/lp2/districting/hess/translation.out"

# firstly read the file
dict = {} # dictionary
with open(DISTRICT_FILE) as f:
    for line in f:
        splitLine = line.split()
        dict[splitLine[0]] = int(splitLine[1])

layer = iface.activeLayer()
layer.startEditing()
layer.dataProvider().addAttributes(
        [QgsField(NEW_DISTRICT_FIELD, QVariant.Int)])
layer.updateFields()

# loop through all GEOIDs
for f in layer.getFeatures():
    geoid = f["GEOID10"]
    if geoid in dict:
        # add district number as attribute
        f[NEW_DISTRICT_FIELD] = dict[geoid]
    else:
        print "Missing district for geoid %s" % geoid
    layer.updateFeature(f)
    print "Done for %s" % geoid
layer.commitChanges()
print "Done"
