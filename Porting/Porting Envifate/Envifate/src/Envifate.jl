module Envifate


# -*- coding: utf-8 -*-
"""
/***************************************************************************
 EnviFate
                                 A QGIS plugin
 EnviFate: Open source tool for environmental risk analysis
                              -------------------
        begin                : 2016-07-15
        git sha              : $Format:%H
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


""" imports
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
"""

using Parameters


@with_kw mutable struct EnviFate
    """QGIS Plugin Implementation."""

    iface
    canvas
    # registry
    plugin_dir
    pluginIsActive
    translator = nothing
    dockwidget = nothing
    # dlg
    # actions
    # menu
    # toolbar
    envifate_menu
    acquifier_item
    lake_item
    river_item
    atm_item
    noise_item
    wnoise_item
    light_item
    sediment_ite
    runoff_item
    thermic_item
    run_leaching
    
    """Constructor.

    :param iface: An interface instance that will be passed to this class
        which provides the hook by which you can manipulate the QGIS
        application at run time.
    :type iface: QgsInterface
    """

    # Save reference to the QGIS interface, initialize plugin directory
    function EnviFate(iface)
        # initialize locale
        locale = QSettings().value("locale/userLocale")[0:2]
        locale_path = joinpath( ef.plugin_dir, "i18n", "EnviFate$locale.qm" )

        if ispath(locale_path)
            translator = QTranslator()
            load( translator, locale_path)
    
            if qVersion() > "4.3.3"
                QCoreApplication.installTranslator(translator)
            end
        end

        # Create the dialog (after translation) and keep reference
        # dlg = OpenRiskDialog()

        # Declare instance attributes
        # actions = []
        # menu = tr( u"&EnviFate" )

        # toolbar = addToolBar( iface, u"EnviFate" )
        # setObjectName( toolbar, u"EnviFate" )

        return new( iface, iface.mapCanvas(), #=QgsMapLayerRegistry.instance(),=# pwd() )
    end   
end



# noinspection PyMethodMayBeStatic
function tr( message )
    """Get the translation for a string using Qt translation API.

    We implement this ourselves since we do not inherit QObject.

    :param message: String for translation.
    :type message: str, QString

    :returns: Translated version of message.
    :rtype: QString
    """
    # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
    return QCoreApplication.translate("EnviFate", message)
end



function initGui( ef::EnviFate )
    """Create the menu entries and toolbar icons inside the QGIS GUI."""


    #icon_path = ':/plugins/OpenRisk/icon1.png'
    ef.envifate_menu = QMenu( QCoreApplication.translate("EnviFate", "&EnviFate") )
    setIcon( ef.envifate_menu, QIcon(":/plugins/envifate/icons/EnviFate_mini.png") )

    ef.acquifer_item = QAction( QIcon(":/plugins/envifate/icons/acquifer_icon.png"),
                                QCoreApplication.translate("EnviFate", "Dispersione in falda"),
                                mainWindow( ef.iface) )
    if status_grass == 1
        ef.transport_item = QAction( QIcon(":/plugins/envifate/icons/acquifer_icon.png"),
                                     QCoreApplication.translate("EnviFate", "Trasporto soluti in falda"),
                                     mainWindow( ef.iface) )
    end
    ef.lake_item = QAction( QIcon(":/plugins/envifate/icons/lake_icon.png"),
                              QCoreApplication.translate("EnviFate", "Dispersione in laghi e bacini"),
                              mainWindow( ef.iface) )
    elf.river_item = QAction( QIcon(":/plugins/envifate/icons/river_icon.png"),
                               QCoreApplication.translate("EnviFate", "Dispersione fluviale"),
                               mainWindow( ef.iface) )
    ef.atm_item = QAction( QIcon(":/plugins/envifate/icons/atm_icon.png"),
                             QCoreApplication.translate("EnviFate", "Dispersione atmosferica"),
                             mainWindow( ef.iface) )
    ef.noise_item = QAction( QIcon(":/plugins/envifate/icons/noise_icon.png"),
                               QCoreApplication.translate("EnviFate", "Analisi rumore"),
                               mainWindow( ef.iface) )
    ef.wnoise_item = QAction( QIcon(":/plugins/envifate/icons/wnoise_icon.png"),
                                QCoreApplication.translate("EnviFate", "Analisi rumore in acqua"),
                                mainWindow( ef.iface) )
    ef.light_item = QAction( QIcon(":/plugins/envifate/icons/light_icon.png"),
                               QCoreApplication.translate("EnviFate", "Analisi inquinamento luminoso"),
                               mainWindow( ef.iface) )
    ef.sediment_item = QAction( QIcon(":/plugins/envifate/icons/sediment_icon.png"),
                                  QCoreApplication.translate("EnviFate", "Analisi sedimentazione marina"),
                                  mainWindow( ef.iface) )
    ef.runoff_item = QAction( QIcon(":/plugins/envifate/icons/runoff_icon.png"),
                                QCoreApplication.translate("EnviFate", "Analisi ruscellamento"),
                                mainWindow( ef.iface) )
    ef.thermic_item = QAction( QIcon(":/plugins/envifate/icons/thermic_icon.png"),
                                 QCoreApplication.translate("EnviFate", "Analisi inquinamento termico"),
                                 mainWindow( ef.iface) )

    connect.(
        [ef.acquifer_item.triggered,ef.lake_item.triggered,ef.river_item.triggered,ef.atm_item.triggered,ef.noise_item.triggered,
         ef.wnoise_item.triggered,ef.light_item.triggered,ef.sediment_item.triggered,ef.runoff_item.triggered,ef.thermic_item.triggered],
        [ef.run_leaching,ef.run_lake,ef.run_river,ef.run_plume,ef.run_noise,ef.run_wnoise,ef.run_light,ef.run_sediment,ef.run_runoff,
         ef.run_thermic]
    )
    if status_grass == 1
        connect( ef.transport_item.triggered, ef.run_transport )
    end

    #ef.envifate_menu.addActions([ef.CreateReceiverPoints_item, ef.CalculateNoiseLevels_item, ef.AssignLevelsToBuildings_item,ef.ApplyNoiseSymbology_item, ef.Informations_item])

    addActions(ef.envifate_menu, [ef.acquifer_item,ef.lake_item,ef.river_item,ef.atm_item,ef.noise_item,ef.wnoise_item,ef.light_item,ef.sediment_item,ef.runoff_item,ef.thermic_item])
    if status_grass == 1
        addActions( ef.envifate_menu, [ef.transport_item] )
    ef.menu = pluginMenu( ef.iface )
    addMenu( ef.menu, ef.envifate_menu )
    if isnothing(ef.dockwidget)
        # Create the dockwidget (after translation) and keep reference
        ef.dockwidget = envifateDockWidget()

    # connect to provide cleanup on closing of dockwidget
    connect( ef.dockwidget.closingPlugin, ef.onClosePlugin )

    # show the dockwidget
    # TODO: fix to allow choice of dock location
    addDockWidget( ef.iface, Qt.LeftDockWidgetArea, ef.dockwidget )
    # items = [QTreeWidgetItem("item: {}".format(i)) for i in range(10)]
    #
    # ef.treeWidget.insertTopLevelItems(0, items)
    connect( ef.dockwidget.treeWidget.itemDoubleClicked, ef.treeMedia_doubleClicked )
    show( ef.dockwidget )
end


function onClosePlugin(ef::EnviFate)
    """Cleanup necessary items here when plugin dockwidget is closed"""

    #print "** CLOSING prova"

    # disconnects
    disconnect( ef.dockwidget.closingPlugin, ef.onClosePlugin )

    # remove this statement if dockwidget is to remain
    # for reuse if plugin is reopened
    # Commented next statement since it causes QGIS crashe
    # when closing the docked window:
    # ef.dockwidget = None

    ef.pluginIsActive = false
end

function unload( ef::EnviFate )
    """Removes the plugin menu item and icon from QGIS GUI."""
    # for action in ef.actions:
    #     ef.iface.removePluginMenu(
    #         ef.tr(u'&EnviFate'),
    #         action)
    #     ef.iface.removeToolBarIcon(action)

    # del ef.toolbar


    if status_grass == 1
        removePluginMenu( ef.iface, "&EnviFate", ef.transport_item)
    removePluginMenu.(Ref(ef.iface),Ref("&EnviFate"),[ef.acquifer_item,ef.lake_item,ef.river_item,ef.atm_item,ef.noise_item,ef.wnoise_item,
                      ef.light_item,ef.sediment_item,ef.runoff_item,ef.thermic_item] )
end


function treeMedia_doubleClicked(ef::EnviFate,index)
    getSelected = selectedItems( ef.dockwidget.treeWidget )
    if getSelected
        baseNode = getSelected[0]
        getChildNode = text( baseNode, 0 )
        if baseNode.text(0) == "Dispersione in falda"
            run_leaching(ef)
        elseif baseNode.text(0)=="Dispersione in laghi e bacini"
            run_lake(ef)
        elseif baseNode.text(0)=="Analisi rumore"
            run_noise(ef)
        elseif baseNode.text(0)=="Analisi rumore in acqua"
            run_wnoise(ef)
        elseif baseNode.text(0)=="Dispersione fluviale"
            run_river(ef)
        elseif baseNode.text(0)=="Dispersione atmosferica"
            run_plume(ef)
        elseif baseNode.text(0)=="Analisi inquinamento luminoso"
            run_light(ef)
        elseif baseNode.text(0)=="Analisi sedimentazione marina"
            run_sediment(ef)
        elseif baseNode.text(0)=="Analisi ruscellamento"
            run_runoff(ef)
        elseif baseNode.text(0)=="Analisi inquinamento termico"
            run_thermic(ef)
    # item = selectedIndexes(ef.dockwidget.treeView)[0]
    # print(text(itemFromIndex(model(item),index)))


function run_leaching(ef::EnviFate)
    """Run method that performs all the real work"""
    # show the dialog
    d = Dialog( do_daf, ef.iface )
    show(d)
    exec_(d)
end

function run_transport(ef::EnviFate)
    """Run method that performs all the real work"""
    # show the dialog
    d = Dialog(do_transport, ef.iface)
    show(d)
    exec_(d)
end

function run_river(ef::EnviFate)
    """Run method that performs all the real work"""
    # show the dialog
    d = Dialog(do_river, ef.iface)
    show(d)
    exec_(d)
end

function run_lake(ef::EnviFate)
    """Run method that performs all the real work"""
    # show the dialog
    d = Dialog(do_lake, ef.iface)
    show(d)
    exec_(d)
end

function run_noise(ef::EnviFate)
    """Run method that performs all the real work"""
    # show the dialog
    d = Dialog(do_noise, ef.iface)
    show(d)
    exec_(d)
end

function run_wnoise(ef::EnviFate)
    """Run method that performs all the real work"""
    # show the dialog
    d = Dialog(do_wnoise, ef.iface)
    show(d)
    exec_(d)
end

function run_plume(ef::EnviFate)
    """Run method that performs all the real work"""
    # show the dialog
    p = Dialog(do_plume, ef.iface)
    show(p)
    exec_(p)
end


function run_light(ef::EnviFate)
    """Run method that performs all the real work"""
    # show the dialog
    p = Dialog(do_light, ef.iface)
    show(p)
    exec_(p)
end


function run_sediment(ef::EnviFate)
    """Run method that performs all the real work"""
    # show the dialog
    p = Dialog(do_sediment, ef.iface)
    show(p)
    exec_(p)
end


function run_runoff(ef::EnviFate)
    """Run method that performs all the real work"""
    # show the dialog
    p = Dialog(do_runoff, ef.iface)
    show(p)
    exec_(p)
end

function run_thermic(ef::EnviFate)
    """Run method that performs all the real work"""
    # show the dialog
    p = Dialog(do_thermic, ef.iface)
    show(p)
    exec_(p)
end

    # ef.dlg_dialog1 = DafDialog()

    # show(ef.dlg_dialog1)

    #show(ef.dlg)
    #popolacombo(ef)
    # Run the dialog event loop
    #result = exec_(ef.dlg)

    # See if OK was pressed
    # if result:
    #     pass

end # module