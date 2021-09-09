# -*- coding: utf-8 -*-
"""
/***************************************************************************
 EnviFate
                                 A QGIS plugin
 EnviFate: Open source tool for environmental risk analysis
                             -------------------
        begin                : 2016-07-15
        copyright            : (C) 2016 by Francesco Geri
        email                : fgeri@icloud.com
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load EnviFate class from file EnviFate.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .envifate import EnviFate
    return EnviFate(iface)
