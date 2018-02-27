#!/usr/bin/env python
# Using the method in Frenkel + Smit "Understanding Molecular Simulation"
#   the Flyvbjerg & Petersen method (J. of Chem. Phys., 1989)

# Standard Library
import sys, os
import math
import argparse

# 3rd-Party Library
import numpy as np
import matplotlib.pyplot as plt 

import utils.parserIO as io
import utils.blocking
from   utils.welford import Welford

#---------------------------------#
#   Get Command line options      #
#---------------------------------#
parser = argparse.ArgumentParser(description="Perform block averaging  \
                                using method of Flyvbjerg + Petersen (JCP, 1989) \
                                found in Appendix D of Frenkel + Smit \
                                'Understanding Molecular Simulation' \
                                of a specified property. ")
parser.add_argument( "filetype", default="LAMMPS_Log", type=str,
                      help="File type / format of FILENAME to use for \
                            reading in data. Default is LAMMPS file format")
parser.add_argument( "filename", type=str,
                      help="Filename that contains data." )
parser.add_argument( "property", type=str,
                     help="String of property in FILENAME that we would \
                           like to perform block averaging over.")
parser.add_argument( "--plot", action="store_true",
                     help="Turn on plotting using PyPlot" )

args = parser.parse_args()

#-----------------------------------------------#
#   Error checking of command line options      #
#-----------------------------------------------#
if not args.filename:
    raise Exception("No filename specified.")

if not args.property:
    raise Exception("No property string specified.")


#-----------------------------#
#   Read in Log File          #
#-----------------------------#
f = io.selectFileType( args.filetype )

thermo_string = args.property
data_dict = f.readData( args.filename, thermo_string )

#-------------------------------#
#    Loop over data             #
#-------------------------------#
for i in data_dict:

    print "Data file chunk: ", i+1
    data = data_dict[i]

    #-------------------------------#
    #    Initialization - Blocking  #
    #-------------------------------#
    xlist = []
    ylist = []
    ylisthi = []
    ylistlo = []

    stats = Welford()
    stats(data)
    totalMean = stats.mean
    totalStd  = stats.std
    totalVar  = stats.std * stats.std

    #--------------------------------#
    #   Perform blocking eval on PE  #
    #--------------------------------#
  
    bData = blocking.selectBlockMethod( "Flyvbjerg+Petersen" )
    output = bData.scanBlocking( data )

    for M in output:
        xlist.append( M )
        ylist.append(   output[M][5] )
        ylisthi.append( output[M][6] )
        ylistlo.append( output[M][7] )
    

    #-----------------------------#
    #   Plot Std. Dev. and its    #
    #   error bars to determine   #
    #   plateau of blocking       #
    #-----------------------------#
    if args.plot:
        x = np.asarray( xlist )
        y = np.asarray( ylist )
        yhi = np.asarray( ylisthi )
        ylo = np.asarray( ylistlo )

        ytop = yhi - y
        ybot = y - ylo

        plt.errorbar( x, y, yerr=(ytop,ybot), fmt='-o' )
    plt.show()

