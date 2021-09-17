#!/usr/bin/python3

import sys

if len(sys.argv) < 4:
    print("Usage : %s <shapefile> <blockfile> <output.png>" % sys.argv[0])
    exit(0)

from qgis.utils import iface
from PyQt5.QtCore import QVariant, QSize
from PyQt5.QtGui import QColor

from random import randrange

from qgis.core import *
from qgis.gui import *


SHP_FILE = sys.argv[1]
DISTRICT_FILE = sys.argv[2]
OUTPUT_FILE = sys.argv[3]
QGIS_PATH = "/usr"

# state = two letter abbreviation
# level = [ counties | tracts ]
#state="OK"
#level="tracts"


########################################################################################

NEW_DISTRICT_FIELD = "district"

# setup Qgs environment without ``running'' QGIS
QgsApplication.setPrefixPath(QGIS_PATH, True)
app = QgsApplication([], True)
app.initQgis()

print("QGIS environment loaded")


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
    geoid = f["GEOID20"]
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
# black -> gray, gradient
outline = [ '#000000', '#696969', '#808080', '#A9A9A9', '#C0C0C0' ]

palette_prism = [ '#5F4690', '#1D6996', '#38A6A5', '#0F8554', '#73AF48', '#EDAD08', '#E17C05', '#CC503E', '#94346E', '#6F4070', '#994E95', '#666666' ]
palette_pastel = [ '#66C5CC', '#F6CF71', '#F89C74', '#DCB0F2', '#87C55F', '#9EB9F3', '#FE88B1', '#C9DB74', '#8BE0A4', '#B497E7', '#D3B484', '#B3B3B3',
'#855C75', '#D9AF6B', '#AF6458', '#736F4C', '#526A83', '#625377', '#68855C', '#9C9C5E', '#A06177', '#8C785D', '#467378', '#7C7C7C' ]
palette_random = QgsLimitedRandomColorRamp.randomColors(len(unique_values))

outline_color = outline[ 0 ]
print("Setting outline color to %s" % outline_color)

for i, unique_value in enumerate(unique_values):
    # initialize the default symbol for this geometry type
    symbol = QgsSymbol.defaultSymbol(layer.geometryType())

    # configure a symbol layer
    layer_style = {}
    layer_style['color'] = palette_pastel[ i % len(palette_pastel) ] # palette_random[ i ].name()
    layer_style['strokeColor'] = outline_color
    layer_style['strokeWidth'] = "1.0"
    symbol_layer = QgsSimpleFillSymbolLayer.create(layer_style)
    symbol_layer.setStrokeColor(QColor(layer_style['strokeColor'])) # must be done due to QGIS bug

    # replace default symbol layer with the configured one
    if symbol_layer is not None:
        symbol.changeSymbolLayer(0, symbol_layer)

    # create renderer object
    category = QgsRendererCategory(unique_value, symbol, str(unique_value))
    # entry for the list of category items
    categories.append(category)
layer.setRenderer( QgsCategorizedSymbolRenderer ( NEW_DISTRICT_FIELD, categories ) )
layer.triggerRepaint(False)

print("Categorizing done")

# render image
QgsProject.instance().addMapLayer(layer)
options = QgsMapSettings()
options.setLayers([layer]) #"district_map"])
options.setBackgroundColor(QColor(255, 255, 255))
options.setOutputDpi(50.0)
options.setOutputSize(QSize(1024, 1024))
options.setExtent(layer.extent())

render = QgsMapRendererParallelJob(options)
def finished():
    img = render.renderedImage()
    img.save(OUTPUT_FILE, "png")
    print("Saving image done")
render.finished.connect(finished)
render.start()
render.waitForFinished()

app.exitQgis()
