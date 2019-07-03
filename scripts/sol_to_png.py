#!/bin/python3

from qgis.utils import iface
from PyQt5.QtCore import QVariant, QSize
from PyQt5.QtGui import QColor

from random import randrange

from qgis.core import *
from qgis.gui import *

DATABASE_PATH = "/home/lykhovyd/progs/gerry/lp2/districting/hess/data/districting/"
DISTRICT_FILE = "/home/lykhovyd/progs/gerry/lp2/districting/hess/translation.out"
QGIS_PATH = "/usr"

# state = two letter abbreviation
# level = [ counties | tracts ]
state="OK"
level="tracts"


########################################################################################

NEW_DISTRICT_FIELD = "district"

# setup Qgs environment without ``running'' QGIS
QgsApplication.setPrefixPath(QGIS_PATH, True)
app = QgsApplication([], True)
app.initQgis()

print("QGIS environment loaded")


folder=""
if level == "tracts":
    folder = "census_tracts"
elif level == "counties":
    folder = level
else:
    print("Unknown level %s, must be counties or tracts!" % level)
    exit(1)

SHP_FILE = DATABASE_PATH + state + "/" + folder + "/maps/" + state + "_" + level + ".shp"

layer = QgsVectorLayer(SHP_FILE, "district_map", "ogr")
#layer = iface.addVectorLayer(SHP_FILE, "Global", "ogr")

if not layer.isValid():
    print("Layer failed to load!")
    exit(1)

# firstly read the file
dict = {} # dictionary
with open(DISTRICT_FILE) as f:
    for line in f:
        splitLine = line.split()
        dict[splitLine[0]] = int(splitLine[1])

layer.startEditing()
if -1 == layer.fields().indexFromName(NEW_DISTRICT_FIELD):
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
        print("Missing district for geoid %s" % geoid)
    layer.updateFeature(f)
    print("Done for %s" % geoid)
layer.commitChanges()
print("Done adding district feature")

# setup categorized renderer
fni = layer.fields().indexFromName(NEW_DISTRICT_FIELD)
unique_values = layer.uniqueValues(fni)
categories = []
palette = [ '255, 0, 0', '0, 255, 0', '0, 0, 255', '255, 255, 0', '0, 255, 255', '255, 0, 255', '128, 128, 128', '128, 0, 0', '0, 128, 0', '128, 128, 0', '0, 128, 128', '0, 0, 128' ]
for i, unique_value in enumerate(unique_values):
    # initialize the default symbol for this geometry type
    symbol = QgsSymbol.defaultSymbol(layer.geometryType())

    # configure a symbol layer
    layer_style = {}
    # layer_style['color'] = '%d, %d, %d' % (randrange(0, 256), randrange(0, 256), randrange(0, 256))
    layer_style['color'] = palette[ i % len(palette) ]
    layer_style['outline'] = '#000000'
    symbol_layer = QgsSimpleFillSymbolLayer.create(layer_style)

    # replace default symbol layer with the configured one
    if symbol_layer is not None:
        symbol.changeSymbolLayer(0, symbol_layer)

    # create renderer object
    category = QgsRendererCategory(unique_value, symbol, str(unique_value))
    # entry for the list of category items
    categories.append(category)
layer.setRenderer( QgsCategorizedSymbolRenderer ( NEW_DISTRICT_FIELD, categories ) )
layer.triggerRepaint()

print("Categorizing done")

# render image
QgsProject.instance().addMapLayer(layer)
options = QgsMapSettings()
options.setLayers([layer]) #"district_map"])
options.setBackgroundColor(QColor(255, 255, 255))
options.setOutputSize(QSize(500, 500))
options.setExtent(layer.extent())

render = QgsMapRendererParallelJob(options)
def finished():
    img = render.renderedImage()
    img.save(state + "_" + level + ".png", "png")
    print("Saving image done")
render.finished.connect(finished)
render.start()
render.waitForFinished()

app.exitQgis()
