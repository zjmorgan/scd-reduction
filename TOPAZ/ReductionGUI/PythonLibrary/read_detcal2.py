class ReadDetCal:
    """ Version 2 only reads the detector calibration parameters.
    It is does not write a header for an integration file."""
    
    def __init__(self):
        # Detector calibration parameters.
        self.L1 = 0.0
        self.t0_shift = 0.0
        self.nod = 0
        # initialize lists
        self.detNum = []
        self.nRows = []
        self.nCols = []
        self.width = []
        self.height = []
        self.depth = []
        self.detD = []
        self.centerX = []
        self.centerY = []
        self.centerZ = []
        self.baseX = []
        self.baseY = []
        self.baseZ = []
        self.upX = []
        self.upY = []
        self.upZ = []

    def read_detcal2(self, input):
        """
         Function to read the .DetCal file, write the header for
         an integrate file, and return the detector calibartion 
         parameters.
           A. J. Schultz,   February, 2012
           
        Version 2 only reads the detector calibration parameters.
        It is does not write a header for an integration file.
        May, 2012
        
        --------------------------------------------------------
        
         Comments from Fortran subroutine:
            This subroutine will read the first lines of a peaks or
            integrate file with the SNS format.
            A. J. Schultz	June, 2008
        
            The first variable of each record is ISENT.
        
            Linux version, January 2002, A.Schultz
        
            ISENT = 0 --> variable name list for describing the histogram
                  = 1 --> variable values of above list
                  = 2 --> variable name list for a reflection
                  = 3 --> variable values of above list
                  = 4 --> variable name list for parameters for detectors
                  = 5 --> variable values of parameters for detectors
                  = 6 --> variable name list: L1    T0_SHIFT
                  = 7 --> variable values: L1    T0_SHIFT
        """
            
        # input is the .DetCal file which is already open.    
        lineString = input.readline() # read first header line from integrate file.
        
        nod = 0     # number of detectors
        
        # begin reading the calibration file
        while True:
            lineString = input.readline()
            lineList = lineString.split()
            j = len(lineList)
            if j == 0: break                    # len = 0 if EOf
                    
            if lineList[0] == '#': continue
            formatFlag = int(lineList[0])       # test for line type
            
            if formatFlag == 7:
                self.L1 = float(lineList[1])         # L1 in centimeters
                self.t0_shift = float(lineList[2])   # t-zero offset in microseconds

            elif formatFlag == 5:
                self.nod += 1
                self.detNum.append( int(lineList[1]) )    # store parameters in lists
                self.nRows.append( int(lineList[2]) )
                self.nCols.append( int(lineList[3]) )
                self.width.append( float(lineList[4]) )
                self.height.append( float(lineList[5]) )
                self.depth.append( float(lineList[6]) )
                self.detD.append( float(lineList[7]) )
                self.centerX.append( float(lineList[8]) )
                self.centerY.append( float(lineList[9]) )
                self.centerZ.append( float(lineList[10]) )
                self.baseX.append( float(lineList[11]) )
                self.baseY.append( float(lineList[12]) )
                self.baseZ.append( float(lineList[13]) )
                self.upX.append( float(lineList[14]) )
                self.upY.append( float(lineList[15]) )
                self.upZ.append( float(lineList[16]) )
