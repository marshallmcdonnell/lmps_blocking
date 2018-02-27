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
 
            out = { "length"    : len(data), 
                    "mean"      : stats.mean,
                    "var"       : var,
                    "var_plus"  : varPlus,
                    "var_minus" : varMinus,
                    "std"       : std, 
                    "std_plus"  : stdPlus, 
                    "std_minus" : stdMinus
            }
            self.blockDict[M] = out
            
           

   
    def printScanBlocking(self) :
        print "M #Blcks  MEAN  VAR  VAR+  VAR-  STD   STD+    STD-"
        print "---------------------------------------------------"
        for k, v in self.blockDict.items():
            print "{M: <2} {length: <4} {mean:.3f} " \
                  "{var:.3f} {var_plus:.3f} {var_minus:.3f} " \
                  "{std:.3f} {std_plus:.3f} {std_minus:.3f} ".format(M=k, **v )

    #-----------------------------#
    #   Plot Std. Dev. and its    #
    #   error bars to determine   #
    #   plateau of blocking       #
    #-----------------------------#
    def plotScanBlocking(self):
        xlist = []
        xlist = []
        ylist = []
        ylisthi = []
        ylistlo = []

        for k, v in self.blockDict.items():
            xlist.append(k)
            ylist.append(v["std"])
            ylisthi.append(v["std_plus"])
            ylistlo.append(v["std_minus"])


        x = np.asarray(xlist)
        y = np.asarray(ylist)
        yhi = np.asarray(ylisthi)
        ylo = np.asarray(ylistlo)

        ytop = yhi - y
        ybot = y - ylo

        plt.errorbar( x, y, yerr=(ytop,ybot), fmt='-o' )
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

