import numpy             as np
import matplotlib.pyplot as plt
import sys,os, math
from welford import Welford
#----------------------------------------------------------------------#
def selectBlockMethod( method ):
    if      method == "Flyvbjerg+Petersen" or method == "FP" \
         or method == "F+P":
        return blockFlyvbjergPetersen( method )

    elif    method == "SetBlocks" or method == "setBlocks" \
         or method == "Standard"  or method == "standard":
        return blockSetBlocks( method )

    else:
        raise Exception('Blocking method "' + method + '" not supported (...yet)')

#----------------------------------------------------------------------#
#   Parent Blocking Class - inherited by all others and calls others   #
#----------------------------------------------------------------------#
class BlockingMethod( object ):
    def __init__( self, blockingMethod ):
        self.type = blockingMethod

    def getBlockedData( self, nblocks, data ):
        return self.block( nblocks, data ) 

#----------------------------------------------------------------------#
#   Using the method in Frenkel + Smit Textbook                        #
#   "Understanding Molecular Simulation"                               #
#   the Flyvbjerg & Petersen method (J. of Chem. Phys., 1989)          #
#----------------------------------------------------------------------#
class blockFlyvbjergPetersen( BlockingMethod ):
    @property
    def method(self): return "Flyvbjerg & Petersen Method:" \
                             "   J. of Chem. Phys., Vol. 91 (1), 461-466"

    def blockData( self, data ):
        dataBlocked = []
        x = iter( data )
        for i in x:
            iPlus1 = next( x, None )
            if iPlus1 == None: break
            dataBlocked.append( (i + iPlus1) / 2.0 )
        return dataBlocked

    def block( self, nblocks, data ):
        blocks = len(data)
        while blocks > nblocks:
            data = self.blockData( data )
            blocks = len(data)
        return data

    def blockNtimes( self, ntimes, data ):
        n = 0
        while n != ntimes:
            data = self.blockData( data )
            n +=1
        return data

    def scanBlocking( self, data ):
        self.blockDict = {}

        print "M #Blcks  MEAN  VAR  VAR+  VAR-  STD   STD+    STD-"
        print "---------------------------------------------------"
        M = 0
        while len(data) >= 2:
            M += 1
            data = self.blockNtimes( 1, data )
            stats = Welford()
            stats(data)

            #mean   = sum(data) / float( len(data) )
            #meansq = mean * mean           
            #var = 0.00
            #for x in data:
            #    var += x*x - meansq
            #var = var / (len(data) - 1.0)
            #print "CHECK -  welford:", stats.std * stats.std, \
            #      "        standard:", var, \
            #      "            diff:", stats.std*stats.std - var
        
            if len(data) < 2: break

            L = float( len(data) )
            Lminus1 = L - 1.0
            
            var      = stats.std*stats.std
            var      = var / Lminus1
            varStd   = math.sqrt( 2.0 / Lminus1 )
            varPlus  = var * ( 1.0 + varStd )
            varMinus = var * ( 1.0 - varStd )

            std      = math.sqrt( var )
            stdStd   = 1.0 / math.sqrt( 2.0 * Lminus1 )
            stdPlus  = std * ( 1.0 + stdStd ) 
            stdMinus = std * ( 1.0 - stdStd ) 
 
            out = [len(data), stats.mean,\
                   var, varPlus, varMinus,\
                   std, stdPlus, stdMinus]
            self.blockDict[M] = out
           

            print "{0: <2} {1: <4} {2:.3f} " \
                  "{3:.3f} {4:.3f} {5:.3f} " \
                  "{6:.3f} {7:.3f} {8:.3f} ".format(  \
                  M, len(data), stats.mean, \
                  var, varPlus, varMinus, \
                  std, stdPlus, stdMinus )
    


        #-----------------------------#
        #   Plot Std. Dev. and its    #
        #   error bars to determine   #
        #   plateau of blocking       #
        #-----------------------------#
        def printScanBlocking( self ):
            xlist = []
            xlist = []
            ylist = []
            ylisthi = []
            ylistlo = []

            for M in self.blockDict:
                xlist.append( M )
                ylist.append(   output[M][5] )
                ylisthi.append( output[M][6] )
                ylistlo.append( output[M][7] )


            x = np.asarray( xlist )
            y = np.asarray( ylist )
            yhi = np.asarray( ylisthi )
            ylo = np.asarray( ylistlo )

            ytop = yhi - y
            ybot = y - ylo

            plt.errorbar( x, y, yerr=(ytop,ybot), fmt='-o' )
            print "PRINTING PLOT"
            plt.show()
            return

#----------------------------------------------------------------------#
#   Using the standard method by picking # of blocks                   #
#----------------------------------------------------------------------#
class blockSetBlocks( BlockingMethod ):
    @property
    def method(self): return "Standard Method:" \
                             "   Explicitly set # of blocks"

    def block( self, nblocks, data ):
        blockLength = len(data) / nblocks
        i=0
        blockAvgs = []
        while i < nblocks:  
            blockData = data[i*blockLength:(i+1)*blockLength]
            stats = Welford()
            stats(blockData)
            blockAvgs.append(stats.mean)
            i += 1
        return blockAvgs

