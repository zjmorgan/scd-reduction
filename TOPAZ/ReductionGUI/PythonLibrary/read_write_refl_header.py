#--------------------------------------------------------
#               function readrefl_header
#--------------------------------------------------------
# Function to read the header of a peaks or integrate
# file and return the detector calibartion parameters.
#--------------------------------------------------------
# Jython version:
#   A. J. Schultz,   November, 2009
#--------------------------------------------------------
#
# Comments from Fortran subroutine:
# !!!	This subroutine will read the first lines of a peaks or
# !!!	integrate file with the SNS format.
# !!!	A. J. Schultz	June, 2008
#
# !	The first variable of each record is ISENT.
#
# !	Linux version, January 2002, A.Schultz
#
# !	ISENT = 0 --> variable name list for describing the histogram
# !	      = 1 --> variable values of above list
# !	      = 2 --> variable name list for a reflection
# !	      = 3 --> variable values of above list
# !	      = 4 --> variable name list for parameters for detectors
# !	      = 5 --> variable values of parameters for detectors
# !	      = 6 --> variable name list: L1    T0_SHIFT
# !	      = 7 --> variable values: L1    T0_SHIFT
#

# read_write_refl_header
# This version reads the header from one file and writes the header
# to another file.
# August, 2012

def read_write_refl_header(input, output):
    "Returns the detector calibration info from the header of a peaks or integrate file."
    
# input is the input peaks or integrate file which is already open.    
    lineString = input.readline()           # read first header line from integrate file
    output.write(lineString)
    print('Reflection file header: ' + lineString)
    
    nod = 0     # number of detectors
    
# define lists
    detNum = []
    nRows = []
    nCols = []
    width = []
    height = []
    depth = []
    detD = []
    centerX = []
    centerY = []
    centerZ = []
    baseX = []
    baseY = []
    baseZ = []
    upX = []
    upY = []
    upZ = []
    
# begin reading the peaks or integrate file
    while True:
        lineString = input.readline()
        output.write(lineString)
        lineList = lineString.split()
        j = len(lineList)
        if j == 0: break                    # len = 0 if EOf
        
        formatFlag = int(lineList[0])       # test for line type
        if formatFlag == 0: break           # finished header, start of peaks
        
        if formatFlag == 7:
            L1 = float(lineList[1])         # L1 in centimeters
            t0_shift = float(lineList[2])   # t-zero offset in microseconds

        elif formatFlag == 5:
            nod = nod + 1
            i = nod - 1                     # array index starts at 0
            detNum.append(int(lineList[1]))    # store parameters in arrays
            nRows.append(int(lineList[2]))
            nCols.append(int(lineList[3]))
            width.append(float(lineList[4]))
            height.append(float(lineList[5]))
            depth.append(float(lineList[6]))
            detD.append(float(lineList[7]))
            centerX.append(float(lineList[8]))
            centerY.append(float(lineList[9]))
            centerZ.append(float(lineList[10]))
            baseX.append(float(lineList[11]))
            baseY.append(float(lineList[12]))
            baseZ.append(float(lineList[13]))
            upX.append(float(lineList[14]))
            upY.append(float(lineList[15]))
            upZ.append(float(lineList[16]))

# finished
    # return L1, t0_shift, nod,       \ i = 0, 1, 2
        # detNum, nRows, nCols,       \ i = 3, 4, 5
        # width, height, depth, detD, \ i = 6, 7, 8, 9
        # centerX, centerY, centerZ,  \ i = 10, 11, 12
        # baseX, baseY, baseZ,        \ i = 13, 14, 15
        # upX, upY, upZ                 i = 16, 17, 18
        
    return L1, t0_shift, nod,       \
        detNum, nRows, nCols,       \
        width, height, depth, detD, \
        centerX, centerY, centerZ,  \
        baseX, baseY, baseZ,        \
        upX, upY, upZ
