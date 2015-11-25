#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  sans titre.py
#  
#  Copyright 2015 pdefromont <pdefromont@advil.obspm.fr>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

from DEUSProfile import *

plotter = DEUSProfile()
plotter.setDataPath("/home/pdefromont/advil/data/home/pdefromont/code/DEUS-EX-Machina/ProfileTracer/data/output/")
plotter.set_power_spectrum_path("/home/pdefromont/advil/data_bingo/Babel/data/")
plotter.load("new_work_2592_2048_output60")
# plotter.globalSmoothSpectrum(10.0/sqrt(10.0), "exp")
plotter.PlotMeanProfile(30.)
plotter.PlotMeanProfile(40.)
plotter.PlotMeanProfile(15.)
plotter.PlotMeanProfile(20.)
