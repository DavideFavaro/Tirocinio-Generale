
from __future__ import print_function
from builtins import str
from builtins import range
from qgis.PyQt import QtCore, QtGui

from PyQt5.QtCore import QSettings, QTranslator, QCoreApplication, Qt, QObject, pyqtSignal, pyqtRemoveInputHook 
from PyQt5.QtGui import QIcon, QStandardItem,QStandardItemModel
from PyQt5.QtWidgets import QAction, QDialog, QFormLayout, QMenu, QComboBox, QTableWidgetItem, QHBoxLayout, QLineEdit, QPushButton, \
                            QWidget, QSpinBox, QTableWidgetItem, QMessageBox, QFileDialog

import sys

# Initialize Qt resources from file resources.py
#import resources


# Import the code for the dialog

#from open_risk_dialog import OpenRiskDialog
import os.path
try:
  import sqlite3
except:
  # fix_print_with_import
  print("librerie per la connessione al database sqlite non trovate")

import qgis
from qgis.core import *
from qgis.gui import *
from qgis.utils import iface


import numpy as np
import math



from osgeo import gdal,ogr,osr

import pdb

from configuration_dialog import ConfigurationDialog 

sys.path.append( os.path.dirname(__file__)+"/../library" )


class Dialog(ConfigurationDialog):
    
    def __init__(self, iface):
        QDialog.__init__(self, iface.mainWindow())
        self.setupUi(self)

        

        # servizi=['https://idt2.regione.veneto.it/geoportal/csw','https://idt2-geoserver.regione.veneto.it/geoserver/ows','https://idt2.regione.veneto.it/gwc/service/wmts'];

        # for i in servizi:
        #     item = QStandardItem(i)
        #     item.setCheckable(True)
        #     model.appendRow(item)


