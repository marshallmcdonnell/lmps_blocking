
import re


def setTypes( data ):
    if isinstance(data, list):
        data = [data]

    data_dict=dict()
    for x in data:
        for key in x[:-1]:
            if x[-1] == "float":
                data_dict[key] = float
            if x[-1] == "str":
                data_dict[key] = str
            if x[-1] == "int":
                data_dict[key] = int
    return data_dict

class FileClass( object ):
    def __init__( self, fileProps ):
        self.type  = fileProps[0]
        self.name  = fileProps[1]
        self.xdataType = fileProps[3]
        self.ydataType = fileProps[5]
        self.xdata = [ fileProps[2], self.xdataType ]
        self.ydata = [ fileProps[4], self.ydataType ]


#---------------------------------------------------------------------------------#
def selectFileType( fileType ):

    lmpsPairCoeffNames= ['LAMMPS_pairCoeff', 'lammps_pairCoeff', 
                         'LAMMPSpairCoeff', 'lammpspairCoeff', 
                         'lmpspairCoeff', 'lmps_pairCoeff'
                         'LAMMPS_pairCoeff', 'lammps_pairCoeff', 
                         'LAMMPSpairCoeff', 'lammpspairCoeff', 
                         'lmp_pairCoeff', 'lmppairCoeff' ]
    lmpsLogNames  = ['LAMMPS_Log', 'lammps_Log', 'LAMMPSLog', 'lammpsLog', \
                      'log.lammps', 'log' ]
    lmpsDumpNames = ['LAMMPS_Dump', 'lammps_Dump', 'LAMMPSDump', 'lammpsDump', \
                     'dump' ]
    lmpsDumpBox   = ['LAMMPS_DumpBox', 'lammps_DumpBox', 'LAMMPSDumpBox', \
                     'lammpsDumpBox',  'dumpBox' ]
    lmpsDumpAtoms   = ['LAMMPS_DumpAtoms', 'lammps_DumpAtoms', \
                       'LAMMPSDumpAtoms',  'lammpsDumpAtoms',  'dumpAtoms' ]
    lmpsComputeNames  = ['LAMMPS_compute', 'lammps_compute',
                         'LAMMPSCompute', 'lammpsCompute',
                         'lammps_compute', 'lmps_compute'   ]
    xyzFileNames  = [ 'xyz', 'XYZ', 'xyzFile', 'XYZFile' ]
    customChemPot = ['Widom', 'ChemPot', 'MyWidom', 'MyChemPot', \
                     'WidomChemPot', 'MyWidomChemPot' ]
    simpleNames   = [ 'columns', 'simple', 'Simple', 'Columns' ]
 
    if      fileType in lmpsPairCoeffNames: return LAMMPSPairCoeff(fileType)

    elif    fileType in lmpsLogNames:     return LAMMPSLog(fileType)

    elif    fileType in lmpsDumpNames:    return LAMMPSDump(fileType)

    elif    fileType in lmpsDumpBox:      return  LAMMPSDumpBoxBounds(FileType)

    elif    fileType in lmpsDumpAtoms:      return  LAMMPSDumpAtomsBounds(FileType)

    elif    fileType in lmpsComputeNames: return  LAMMPSCompute(FileType)

    elif    fileType in xyzFileNames:     return XYZ(fileType)

    elif    fileType in customChemPot:    return ChemicalPot(fileType)

    elif    fileType in simpleNames:      return SimpleColumns(fileType) 

    else:
        raise Exception('File format "' + fileType + '" is not supported.')

#---------------------------------------------------------------------------------#
#   Parent File Type Class - Inherited by all others and calls other classes      #
#---------------------------------------------------------------------------------#
class FileType(object):
    def __init__(self, fileType):
        self.fileType = fileType

    def readData(self, fileName=None, xcol=None, ycol=None, \
                       section=None,  key=None):
        if not fileName:
            raise Exception("Must pass in filename for readData(...).")

        '''
        if not xcol and not ycol:
            raise Exception("Must pass in either the X column string " \
                            "or Y column string(s) " \
                            "at the very least to readData(...)." )
        '''

        if xcol and not isinstance(xcol, dict):
            raise Exception("In readData(...), xcol must be " \
                            "dictionaries w/ column string as key and "\
                            "data type as value.")

        if ycol and not isinstance(ycol,dict):
            raise Exception("In readData(...), ycol must be " \
                            "dictionaries w/ column string as key and "\
                            "data type as value.")

        if xcol and not all( isinstance(type(xcol[k]),type) for k in xcol):
            raise Exception("In readData(...), xcol is a " \
                            "dictionary but you must have the "\
                            "column string as key and data type as value. "\
                            "Exception was raised due to value not being " \
                            "data type.")

        if ycol and not all( isinstance(type(ycol[k]),type) for k in ycol):
            raise Exception("In readData(...), ycol is a " \
                            "dictionary but you must have the "\
                            "column string as key and data type as value. "\
                            "Exception was raised due to value not being " \
                            "data type.")
        if (section and not key) or (not section and key):
            raise Exception("In readData(...), section and key must both be " \
                            "specified to get a list back for that section " \
                            "and key in the dictionary." )

        self.fileName  = fileName
        data = self.data( xcol, ycol )
        if section and key:
            data = getListFor( data, section, key )
        return data
 
    # gets list from the dictionary of dictionaries for given section and key in
    # that section: 
    #   - The first dictionary contains sections / timesteps as keys
    #     and a dictionary of the data as the value
    #   - The second dictionary / value of the 1st dictionary has the 
    #     headers of the data as the key and a list of the data as the values
    # this function  
    def getListFor( self, data, section, key ):
        return list( data[section][key] )
       
    @property
    def file_obj(self): 
        if self.fileName:
            return open(self.fileName, 'r')
        else:
            raise Exception("Filename is not set.")

    def getData( self, header, rows, xcol, ycol ):
        xID = None
        yID = None
        # X column header
        if xcol: 
            if set(xcol).issubset(set(header)):
                xID = dict( (key, header.index(key)) for key in xcol)
            else:
                raise Exception( "ERROR: " + str(set(xcol)) + " NOT Found in "\
                                  + self.fileName + " LAMMPS log file "\
                                  "header section.")

        # Y column header(s)
        if ycol:
            if set(ycol).issubset(set(header)):
                yID = dict( (key, header.index(key) ) for key in ycol)
            else:
                raise Exception( "ERROR: One of the following column headers: " \
                                  + str(set(ycol)) + " was NOT Found in "\
                                  + self.fileName + " LAMMPS log file "\
                                  "header section.")

        return self.convertData( rows, xcol=xcol, xID=xID, \
                                       ycol=ycol, yID=yID )
          
    def convertData( self, rows, xcol=None, xID=None, ycol=None, yID=None ):
        # Get data, where data{key=column header, value=data list}
        # and column header = [x column, ycolumns=(y column 1, y column 2, ...) ]
        data = dict()
        for row in rows:
            if row:
                values = row.split()
                if xcol:
                    for key in xcol:
                        convertType = xcol[key]
                        index       = xID[key]
                        if key in data:
                            value = float( values[index] )
                            data[key].append( convertType(value) )
                        else:
                            value = float( values[index] )
                            data[key] = [convertType(value)]
                if ycol:
                    for key in ycol:
                        convertType = ycol[key]
                        index       = yID[key]

                        #--- string case ---#

                        if convertType == str:
                            if key in data:
                                data[key].append( convertType(values[index]) )
                            else:
                                data[key] = [convertType(values[index])]

                        #--- number case ---#

                        else:
                            if key in data:
                                value = float( values[index] )
                                data[key].append( convertType(value) )
                            else:
                                value = float( values[index] )
                                data[key] = [convertType(value)]
        return data


#---------------------------------------------------------------------------------#
#   LAMMPS Log File Class                                                         #
#---------------------------------------------------------------------------------#
class LAMMPSPairCoeff(FileType):
    
    @property
    def format(self): 
        descript = ("LAMMPS-style Pair Coeff File - "
                    "Large-scale Atomic / Molecular " 
                    "Massively Parallel Simulator.")
        return descript

    def data(self, xcol, ycol):
        eps=1
        sig=2
        data = dict()
        with self.file_obj as f:
            logfile = f.read()
            pattern = re.compile(r'pair_coeff\s+(.*?)\s+^(?!pair_coeff$)', re.DOTALL)
            for i, lines in enumerate( re.findall(pattern, logfile) ):
                data[i] = self.getDataCols( lines[0], eps, sig )

        return data

    def getDataCols( self, logData, xcol, ycol ):
        logData = logData.split("\n")
        rows    = logData[0:]
        if xcol:
            xID     = dict( (key, int(key)-1) for key in xcol)
        if ycol:
            yID     = dict( (key, int(key)-1) for key in ycol)
        
        if xcol and ycol:
            data = self.convertData( rows, xcol=xcol, xID=xID, \
                                           ycol=ycol, yID=yID )
        elif xcol and not ycol:
            data = self.convertData( rows, xcol=xcol, xID=xID ) 
        
        elif not xcol and ycol:
            data = self.convertData( rows, ycol=ycol, yID=yID ) 

        return data

#---------------------------------------------------------------------------------#
#   LAMMPS Log File Class                                                         #
#---------------------------------------------------------------------------------#
class LAMMPSLog(FileType):
    
    @property
    def format(self): 
        descript = ("LAMMPS Log File - Large-scale Atomic / Molecular " 
                    "Massively Parallel Simulator.")
        return descript

    def data(self, xcol=None, ycol=None):
        data = dict()
        with self.file_obj as f:
            logfile = f.read()

            pattern = re.compile(r'Mbytes\s+(.*?)\s+($|Loop|WARNING)', re.DOTALL)
            for i, lines in enumerate( re.findall(pattern, logfile) ):
                data[i] = self.getDataLAMMPS( lines[0], xcol, ycol )

        return data
 
    def getDataLAMMPS( self, logData, xcol, ycol):
        logData = logData.split("\n")
        header  = logData[0].split()
        rows    = logData[1:]
        data    = self.getData( header, rows, xcol, ycol )
        return data

#----------------------------------------------------------------------------#
#   LAMMPS Dump File Class                                                   #
#----------------------------------------------------------------------------#
class LAMMPSDump(FileType):
    @property
    def format(self): 
        descript = ( 'LAMMPS Dump File - Large-scale Atomic / Molecular ' 
                    'Massively Parallel Simulator.' )
        return descript

    def data(self, xcol, ycol):
        data = dict()
        with self.file_obj as f:
            dumpfile = f.read()
            #pattern = re.compile(r'TIMESTEP(.*?)ITEM:.*?ITEM:\s+ENTRIES(.*?)($|ITEM)', re.DOTALL)
            pattern = re.compile(r'TIMESTEP(.*?)ITEM:.*?ITEM:\s+ATOMS(.*?)($|ITEM)', re.DOTALL)
            for lines in re.findall(pattern, dumpfile):
                data[int(lines[0])] = self.getDataLAMMPS( lines[1], xcol, ycol )
            
        return data
 
    def getDataLAMMPS( self, dumpData, xcol, ycol ):
        dumpData = dumpData.split("\n")
        header  = dumpData[0].split()
        rows    = dumpData[1:]
        data    = self.getData( header, rows, xcol, ycol )
        return data

#----------------------------------------------------------------------------#
#   LAMMPS Dump File Class - Specifically get box bounds                     #
#----------------------------------------------------------------------------#
class LAMMPSDumpBoxBounds(FileType):
    @property
    def format(self): 
        descript = ( 'LAMMPS Dump File - Large-scale Atomic / Molecular ' 
                    'Massively Parallel Simulator. Get Box Bounds sections' )
        return descript

    def data(self, xcol=None, ycol=None):
        data = dict()
        with self.file_obj as f:
            dumpfile = f.read()
            #pattern = re.compile(r'TIMESTEP(.*?)ITEM:.*?ITEM:\s+ENTRIES(.*?)($|ITEM)', re.DOTALL)
            pattern = re.compile(r'TIMESTEP(.*?)ITEM:.*?ITEM:\s+BOX\sBOUNDS(.*?)($|ITEM)', re.DOTALL)
            for lines in re.findall(pattern, dumpfile):
                data[int(lines[0])] = self.getDataLAMMPSBox( lines[1] )
            
        return data
 
    def getDataLAMMPSBox( self, dumpData):
        dumpData = dumpData.split("\n")
        header  = dumpData[0].split()
        rows    = dumpData[1:]
        data = dict()
        data['x'] = [ header[0], float(rows[0].split()[0]), float(rows[0].split()[1]) ]
        data['y'] = [ header[1], float(rows[1].split()[0]), float(rows[1].split()[1]) ]
        data['z'] = [ header[2], float(rows[2].split()[0]), float(rows[2].split()[1]) ]
        return data

#----------------------------------------------------------------------------#
#   LAMMPS Dump File Class - Specifically get # of atoms                     #
#----------------------------------------------------------------------------#
class LAMMPSDumpAtomsBounds(FileType):
    @property
    def format(self): 
        descript = ( 'LAMMPS Dump File - Large-scale Atomic / Molecular ' 
                    'Massively Parallel Simulator. Get NUMBER OF ATOMS sections' )
        return descript

    def data(self, xcol=None, ycol=None):
        data = dict()
        with self.file_obj as f:
            dumpfile = f.read()
            #pattern = re.compile(r'TIMESTEP(.*?)ITEM:.*?ITEM:\s+ENTRIES(.*?)($|ITEM)', re.DOTALL)
            pattern = re.compile(r'TIMESTEP(.*?)ITEM:\s+NUMBER\sOF\sATOMS(.*?)($|ITEM)', re.DOTALL)
            for lines in re.findall(pattern, dumpfile):
                data[int(lines[0])] = self.getDataLAMMPSAtoms( lines[1] )
            
        return data
 
    def getDataLAMMPSAtoms( self, dumpData):
        return dumpData.split()[0]

#----------------------------------------------------------------------------#
#   LAMMPS Compute Dumps                                                     #
#----------------------------------------------------------------------------#
class LAMMPSCompute(FileType):
    @property
    def format(self): 
        descript = ( 'LAMMPS Compute Dump File - Large-scale Atomic / Molecular ' 
                    'Massively Parallel Simulator. Get Box Bounds sections' )
        return descript

    def data(self, xcol=None, ycol=None):
        data = dict()
        with self.file_obj as f:
            # pre-process file to pull out header (3rd line) and body (5th to end)
            dumpfile = f.read().split('\n')
            lines = str(dumpfile[2].split("# ")[1] + '\n')
            cols=len(lines.split())
            for line in dumpfile[4:]:
                if len(line.split()) != cols: continue
                lines += str(line + '\n')
            data[0] = self.getDataLAMMPSCompute( lines, xcol, ycol )
        return data
 
    def getDataLAMMPSCompute( self, logData, xcol, ycol ):
        logData = logData.split("\n")
        header  = logData[0].split()
        rows    = logData[1:]
        data    = self.getData( header, rows, xcol, ycol )
        return data

#---------------------------------------------------------------------------------#
#   XYZ Trajectory File Class                                                     #
#---------------------------------------------------------------------------------#
class XYZ(FileType):
    
    @property
    def format(self): 
        descript = ("XYZ - Chemical File Format ")
        return descript

    def data(self, xcol, ycol):
        data = dict()
        with self.file_obj as f:
            logfile = f.read()
            match = re.search("(.*)\n", logfile )
            natoms = str(match.group(1))
            pattern =  r"" + re.escape(natoms) + "\D*(.*)"
            pattern = re.compile( pattern, re.DOTALL)
            match = re.findall(pattern, logfile )
            for i, lines in enumerate(match):
                data[i] = self.getDataXYZ( lines, xcol, ycol )
        return data
 
    def getDataXYZ( self, logData, xcol, ycol):
        logData = logData.split("\n")
        header  = ["element", "x", "y", "z"] 
        rows    = logData[0:] 
        data    = self.getData( header, rows, xcol, ycol )
        return data

#---------------------------------------------------------------------------------#
#   Chemical Potential Log File Class                                       #
#---------------------------------------------------------------------------------#
class ChemicalPot(FileType):
    
    @property
    def format(self): 
        descript = ('Chemical Potential Custom Format  - for output'
                    'of chemical potential for Widom Particle Insertion code'
                    'using LAMMPS as a shared library OR '
                    'for Boulougoris + Theodorou'
                    'Particle Deletion code using pe/atom from LAMMPS'
                    'accessible volume FORTRAN code.')
        return descript
    def data(self, xcol, ycol):
        data = dict()
        with self.file_obj as f:
            logfile = f.read()
            pattern = re.compile(r'------+\s+(.*)', re.DOTALL )
            for i, lines in enumerate( re.findall( pattern, logfile ) ):
                data[i] = self.getDataChemPot( lines, xcol, ycol )
        return data

    def getDataChemPot( self, logData, xcol, ycol ):
        logData = logData.split("\n")
        header  = logData[0].split()
        rows    = logData[1:]
        data    = self.getData( header, rows, xcol, ycol )
        return data

#---------------------------------------------------------------------------------#
#   Simple Column Data (No Header) File Class                                     #
#---------------------------------------------------------------------------------#
class SimpleColumns(FileType):
    
    @property
    def format(self): 
        descript = ('Simple Column Format  - for output'
                    'of columns with no headers')
        return descript

    def data(self, xcol, ycol):
        data = dict()
        with self.file_obj as f:
            lines = f.read()
            data[0] = self.getDataCols(lines, xcol, ycol )
        return data

    def getDataCols( self, logData, xcol, ycol ):
        logData = logData.split("\n")
        rows    = logData[0:]
        if xcol:
            xID     = dict( (key, int(key)-1) for key in xcol)
        if ycol:
            yID     = dict( (key, int(key)-1) for key in ycol)
        
        if xcol and ycol:
            data = self.convertData( rows, xcol=xcol, xID=xID, \
                                           ycol=ycol, yID=yID )
        elif xcol and not ycol:
            data = self.convertData( rows, xcol=xcol, xID=xID ) 
        
        elif not xcol and ycol:
            data = self.convertData( rows, ycol=ycol, yID=yID ) 

        return data
