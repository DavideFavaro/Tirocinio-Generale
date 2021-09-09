# -*- coding: utf-8 -*-
"""
/***************************************************************************
 EnviFate
                                 A QGIS plugin
 EnviFate: Open source tool for environmental risk analysis
                              -------------------
        begin                : 2016-07-15
        git sha              : $Format:%H$
        copyright            : (C) 2016 by Francesco Geri
        email                : francescogeri@tim.it
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
from __future__ import print_function
from __future__ import absolute_import
from builtins import object
from PyQt5.QtCore import QSettings, QTranslator, QCoreApplication, Qt, QObject, pyqtSignal, pyqtRemoveInputHook
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QAction, QDialog, QFormLayout, QMenu, QTreeWidget, QTreeWidgetItem
import sys

# Initialize Qt resources from file resources.py
from . import resources

from .envifate_dockwidget import envifateDockWidget



grass7path = r'C:\OSGeo4W\apps\grass\grass-7.2.svn'
grass7bin_win = r'C:\OSGeo4W\bin\grass72svn.bat'
# Linux
#grass7bin_lin = 'grass72'
grass7bin_lin = 'grass'
# MacOSX
grass7bin_mac = '/Applications/GRASS/GRASS-7.1.app/'

# Import the code for the dialog


import os.path
try:
  import sqlite3
except:
  # fix_print_with_import
  print("librerie per la connessione al database sqlite non trovate")

#from qgis.core import QgsMapLayerRegistry

#import pdb
sys.path.append( os.path.dirname(__file__)+"/library" )
sys.path.append(os.path.dirname(__file__) + '/tools')
import subprocess

import pdb


# if sys.platform.startswith('linux'):
#     # we assume that the GRASS GIS start script is available and in the PATH
#     # query GRASS 7 itself for its GISBASE
#     grass7bin = grass7bin_lin
# elif sys.platform.startswith('win'):
#     grass7bin = grass7bin_win
# else:
#     grass7bin = grass7bin_mac

# startcmd = grass7bin + ' --config path'


# p = subprocess.Popen(startcmd, shell=True,
#                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# out, err = p.communicate()
# if p.returncode != 0:
#     status_grass=0


######per evitare l'errore con il porting di grass ########

status_grass=0

import functions, leaching,do_daf,do_plume,do_river,do_lake,do_noise,do_wnoise, do_light, do_sediment, do_runoff, do_thermic

if status_grass==1:
    import do_transport

#from daf_dialog import DafDialog

class EnviFate(object):
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        self.canvas=self.iface.mapCanvas()
        #self.registry = QgsMapLayerRegistry.instance()
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'EnviFate{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)
        self.dockwidget = None

        # Create the dialog (after translation) and keep reference
        #self.dlg = OpenRiskDialog()

        # Declare instance attributes
        # self.actions = []
        # self.menu = self.tr(u'&EnviFate')

        # self.toolbar = self.iface.addToolBar(u'EnviFate')
        # self.toolbar.setObjectName(u'EnviFate')

    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('EnviFate', message)




    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""


        #icon_path = ':/plugins/OpenRisk/icon1.png'
        self.envifate_menu = QMenu(QCoreApplication.translate("EnviFate", "&EnviFate"))
        self.envifate_menu.setIcon(QIcon(":/plugins/envifate/icons/EnviFate_mini.png"))

        self.acquifer_item = QAction(QIcon(":/plugins/envifate/icons/acquifer_icon.png"),
                                        QCoreApplication.translate("EnviFate", "Dispersione in falda"), self.iface.mainWindow())

        if status_grass==1:
            self.transport_item = QAction(QIcon(":/plugins/envifate/icons/acquifer_icon.png"),QCoreApplication.translate("EnviFate", "Trasporto soluti in falda"), self.iface.mainWindow())


        self.lake_item = QAction(QIcon(":/plugins/envifate/icons/lake_icon.png"),
                                        QCoreApplication.translate("EnviFate", "Dispersione in laghi e bacini"), self.iface.mainWindow())

        self.river_item = QAction(QIcon(":/plugins/envifate/icons/river_icon.png"),
                                        QCoreApplication.translate("EnviFate", "Dispersione fluviale"), self.iface.mainWindow())

        self.atm_item = QAction(QIcon(":/plugins/envifate/icons/atm_icon.png"),
                                        QCoreApplication.translate("EnviFate", "Dispersione atmosferica"), self.iface.mainWindow())

        self.noise_item = QAction(QIcon(":/plugins/envifate/icons/noise_icon.png"),
                                        QCoreApplication.translate("EnviFate", "Analisi rumore"), self.iface.mainWindow())

        self.wnoise_item = QAction(QIcon(":/plugins/envifate/icons/wnoise_icon.png"),
                                        QCoreApplication.translate("EnviFate", "Analisi rumore in acqua"), self.iface.mainWindow())


        self.light_item = QAction(QIcon(":/plugins/envifate/icons/light_icon.png"),
                                        QCoreApplication.translate("EnviFate", "Analisi inquinamento luminoso"), self.iface.mainWindow())


        self.sediment_item = QAction(QIcon(":/plugins/envifate/icons/sediment_icon.png"),
                                        QCoreApplication.translate("EnviFate", "Analisi sedimentazione marina"), self.iface.mainWindow())


        self.runoff_item = QAction(QIcon(":/plugins/envifate/icons/runoff_icon.png"),
                                        QCoreApplication.translate("EnviFate", "Analisi ruscellamento"), self.iface.mainWindow())

        self.thermic_item = QAction(QIcon(":/plugins/envifate/icons/thermic_icon.png"),
                                        QCoreApplication.translate("EnviFate", "Analisi inquinamento termico"), self.iface.mainWindow())

        self.acquifer_item.triggered.connect(self.run_leaching)
        if status_grass==1:
            self.transport_item.triggered.connect(self.run_transport)
        self.lake_item.triggered.connect(self.run_lake)
        self.river_item.triggered.connect(self.run_river)
        self.atm_item.triggered.connect(self.run_plume)
        self.noise_item.triggered.connect(self.run_noise)
        self.wnoise_item.triggered.connect(self.run_wnoise)
        self.light_item.triggered.connect(self.run_light)
        self.sediment_item.triggered.connect(self.run_sediment)
        self.runoff_item.triggered.connect(self.run_runoff)
        self.thermic_item.triggered.connect(self.run_thermic)

        #self.envifate_menu.addActions([self.CreateReceiverPoints_item, self.CalculateNoiseLevels_item, self.AssignLevelsToBuildings_item,self.ApplyNoiseSymbology_item, self.Informations_item])

        self.envifate_menu.addActions([self.acquifer_item,self.lake_item,self.river_item,self.atm_item,self.noise_item,self.wnoise_item,self.light_item,self.sediment_item,self.runoff_item,self.thermic_item])
        if status_grass==1:
            self.envifate_menu.addActions([self.transport_item])
        self.menu = self.iface.pluginMenu()
        self.menu.addMenu( self.envifate_menu )
        if self.dockwidget == None:
            # Create the dockwidget (after translation) and keep reference
            self.dockwidget = envifateDockWidget()

        # connect to provide cleanup on closing of dockwidget
        self.dockwidget.closingPlugin.connect(self.onClosePlugin)

        # show the dockwidget
        # TODO: fix to allow choice of dock location
        self.iface.addDockWidget(Qt.LeftDockWidgetArea, self.dockwidget)
        # items = [QTreeWidgetItem("item: {}".format(i)) for i in range(10)]
        #
        # self.treeWidget.insertTopLevelItems(0, items)
        self.dockwidget.treeWidget.itemDoubleClicked.connect(self.treeMedia_doubleClicked)
        self.dockwidget.show()



    def onClosePlugin(self):
        """Cleanup necessary items here when plugin dockwidget is closed"""

        #print "** CLOSING prova"

        # disconnects
        self.dockwidget.closingPlugin.disconnect(self.onClosePlugin)

        # remove this statement if dockwidget is to remain
        # for reuse if plugin is reopened
        # Commented next statement since it causes QGIS crashe
        # when closing the docked window:
        # self.dockwidget = None

        self.pluginIsActive = False

    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        # for action in self.actions:
        #     self.iface.removePluginMenu(
        #         self.tr(u'&EnviFate'),
        #         action)
        #     self.iface.removeToolBarIcon(action)

        # del self.toolbar
        self.iface.removePluginMenu("&EnviFate", self.acquifer_item)
        if status_grass==1:
            self.iface.removePluginMenu("&EnviFate", self.transport_item)
        self.iface.removePluginMenu("&EnviFate", self.lake_item)
        self.iface.removePluginMenu("&EnviFate", self.river_item)
        self.iface.removePluginMenu("&EnviFate", self.atm_item)
        self.iface.removePluginMenu("&EnviFate", self.noise_item)
        self.iface.removePluginMenu("&EnviFate", self.wnoise_item)
        self.iface.removePluginMenu("&EnviFate", self.light_item)
        self.iface.removePluginMenu("&EnviFate", self.sediment_item)
        self.iface.removePluginMenu("&EnviFate", self.runoff_item)
        self.iface.removePluginMenu("&EnviFate", self.thermic_item)







    def treeMedia_doubleClicked(self,index):
        getSelected = self.dockwidget.treeWidget.selectedItems()
        if getSelected:
            baseNode = getSelected[0]
            getChildNode = baseNode.text(0)
            if baseNode.text(0)=='Dispersione in falda':
                self.run_leaching()
            elif baseNode.text(0)=='Dispersione in laghi e bacini':
                self.run_lake()
            elif baseNode.text(0)=='Analisi rumore':
                self.run_noise()
            elif baseNode.text(0)=='Analisi rumore in acqua':
                self.run_wnoise()
            elif baseNode.text(0)=='Dispersione fluviale':
                self.run_river()
            elif baseNode.text(0)=='Dispersione atmosferica':
                self.run_plume()
            elif baseNode.text(0)=='Analisi inquinamento luminoso':
                self.run_light()
            elif baseNode.text(0)=='Analisi sedimentazione marina':
                self.run_sediment()
            elif baseNode.text(0)=='Analisi ruscellamento':
                self.run_runoff()
            elif baseNode.text(0)=='Analisi inquinamento termico':
                self.run_thermic()
        # item = self.dockwidget.treeView.selectedIndexes()[0]
        # print(item.model().itemFromIndex(index).text())


    def run_leaching(self):
        """Run method that performs all the real work"""
        # show the dialog
        d = do_daf.Dialog(self.iface)
        d.show()
        d.exec_()

    def run_transport(self):
        """Run method that performs all the real work"""
        # show the dialog
        d = do_transport.Dialog(self.iface)
        d.show()
        d.exec_()

    def run_river(self):
        """Run method that performs all the real work"""
        # show the dialog
        d = do_river.Dialog(self.iface)
        d.show()
        d.exec_()


    def run_lake(self):
        """Run method that performs all the real work"""
        # show the dialog
        d = do_lake.Dialog(self.iface)
        d.show()
        d.exec_()

    def run_noise(self):
        """Run method that performs all the real work"""
        # show the dialog
        d = do_noise.Dialog(self.iface)
        d.show()
        d.exec_()

    def run_wnoise(self):
        """Run method that performs all the real work"""
        # show the dialog
        d = do_wnoise.Dialog(self.iface)
        d.show()
        d.exec_()

    def run_plume(self):
        """Run method that performs all the real work"""
        # show the dialog
        p = do_plume.Dialog(self.iface)
        p.show()
        p.exec_()


    def run_light(self):
        """Run method that performs all the real work"""
        # show the dialog
        p = do_light.Dialog(self.iface)
        p.show()
        p.exec_()


    def run_sediment(self):
        """Run method that performs all the real work"""
        # show the dialog
        p = do_sediment.Dialog(self.iface)
        p.show()
        p.exec_()


    def run_runoff(self):
        """Run method that performs all the real work"""
        # show the dialog
        p = do_runoff.Dialog(self.iface)
        p.show()
        p.exec_()

    def run_thermic(self):
        """Run method that performs all the real work"""
        # show the dialog
        p = do_thermic.Dialog(self.iface)
        p.show()
        p.exec_()

        # self.dlg_dialog1=DafDialog()

        # self.dlg_dialog1.show()

        #self.dlg.show()
        #self.popolacombo()
        # Run the dialog event loop
        #result = self.dlg.exec_()

        # See if OK was pressed
        # if result:
        #     pass
