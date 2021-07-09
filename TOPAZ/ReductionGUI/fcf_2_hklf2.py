#--------------------------------------------------------------------
#                             fcf_2_hkl.py
#--------------------------------------------------------------------

#   X.P. Wang
#   July, 2014
#
#   Script to rewrite a fcf file (LIST 4) to an hkl file (HKLF 2 Laue format)
#   w/ options to remove outlier reflections
#
#   First run SHELX with MERG 0, HKLF 2 and LIST 4.
#
#   Based on ref2hkl.f and fcf_2_hkl.py by  A. J. Schultz, August, 2011

import os
import sys
import numpy
import math

def inputCL(number, *types):
    for t in types:
        try:
            return t(number)
        except ValueError as TypeError:
            continue
    return number

def readFormatted_text(s, *args):
    position = 0
    for length in args:
        yield s[position:position + length]
        position += length

#$TOPAZ HKLF2 data with direction consines
#line ='   9   4   2   72.18   48.08  47 0.70021 0.04477  0.60802  0.09183 -0.74489  0.83642  0.61242 -0.67559 22406  76419 0.9053  39  0.79255  0.9076 176.00  20.00
#          4   4   4       8       8   4       8       8        9        9        9        9        9        9     6      7      7   4        9       8      7      7  
#       1234123412341234567812345678123412345678123456781234567891234567891234567891234567891234567891234567891234561234567123456712341234567891234567812345671234567
#line ='   3  -9  -4 4139.45  126.41   1 1.91195 0.13041 -0.29280  0.74848 -0.70844 -0.65601 -0.64217 -0.09708  8071      1 0.9545  17  1.88381  1.1811 101.00  26.00\n'

print('')
print('Script to rewrite a fcf file (LIST 4) to an hkl file (HKLF 2 Laue format).')
print('First run SHELX with MERG 0, HKLF 2 and LIST 4.')
print('')

#Read the name of the input hkl file
raw_inputFile = input('File name of input hkl file: ')

raw_hkl =[]
dir_cos_1 = [None] * 3
dir_cos_2 = [None] * 3

with open(raw_inputFile, "r") as raw_hkl_input:
    for line in raw_hkl_input:
        line.strip("\n")
        HKL = list(readFormatted_text(line, 4,4,4,8,8,4,8,8,9,9,9,9,9,9,6,7,7,4,9,8,7,7))
        h = int(HKL[0])
        k = int(HKL[1])
        l = int(HKL[2])
        if h==0 and k==0 and l==0: break
        fsq = float(HKL[3])
        sigfsq = float(HKL[4])
        hstnum = int(HKL[5])
        wl = float(HKL[6])
        tbar = float(HKL[7])
        dir_cos_1[0] = float(HKL[8])
        dir_cos_2[0] = float(HKL[9])
        dir_cos_1[1] = float(HKL[10])
        dir_cos_2[1] = float(HKL[11])
        dir_cos_1[2] = float(HKL[12])
        dir_cos_2[2] = float(HKL[13])
        curhst = int(HKL[14])
        seqnum = int(HKL[15])
        transmission = float(HKL[16])
        dn = int(HKL[17])
        twoth = float(HKL[18])
        dsp = float(HKL[19])
        col = float(HKL[20])
        row = float(HKL[21])

        raw_hkl.append([ h, k, l, fsq, sigfsq, hstnum, wl, tbar,
            dir_cos_1[0], dir_cos_2[0], dir_cos_1[1], dir_cos_2[1], 
            dir_cos_1[2], dir_cos_2[2], curhst, seqnum, transmission,
            dn, twoth, dsp, col, row])


# Read the name of the input fcf file
inputFile = input('File name of input fcf file: ')
hkl_input = open(inputFile, 'r')

# Read the name of the output hkl file
outputFile = input('File name of output hkl file: ')
outputFile.rstrip()
hkl_output = open(outputFile, 'w')

iHKL = []
hkl_all =[]
hkl_omit = []
hkl_used = []
raw_index = 0
foc  = 0
fco  = 0
fabs = 0
while True:
    lineString = hkl_input.readline()
    lineList = lineString.split()
    if len( lineList ) == 0: continue
    if lineList[0] == '_refln_observed_status': break

#'INPUT MAXIMUM FOSQ/FCSQ RATIO <0 FOR NO TEST>: '
r_foc = inputCL(eval(input('INPUT MAXIMUM FOSQ/FCSQ RATIO <0 FOR NO TEST>: ')), int, float)
#'INPUT MAXIMUM FCSQ/FOSQ RATIO <0 FOR NO TEST>: '
r_fco = inputCL(eval(input('INPUT MAXIMUM FCSQ/FOSQ RATIO <0 FOR NO TEST>: ')), int, float)
#'INPUT MAXIMUM ABS(FOSQ - FCSQ)/SIGFOSQ RATIO <0 FOR NO TEST>: '
r_fabs = inputCL(eval(input('INPUT MAXIMUM ABS(FOSQ - FCSQ)/SIGFOSQ RATIO <0 FOR NO TEST>: ')), int, float)
#'INPUT OVERALL SCALE FACTOR <1.0>: '

while True:
    lineString = hkl_input.readline()
    lineList = lineString.split()
    if len(lineList) == 0: break
    
    h = int( lineList[0] )
    k = int( lineList[1] )
    l = int( lineList[2] )
    fsq_calc = float("%0.2f" %float(lineList[3]))
    fsq_obs = float("%0.2f" %float(lineList[4]))
    sigfsq = float("%0.2f" %float(lineList[5]))
    iScale = 1

    if fsq_calc == 0: 
        foc = 1000.0 #System absence violations
        fco = 0.0
    elif fsq_obs == 0: 
        fco = 1000.0
    else:
        foc = fsq_obs/fsq_calc
        fco = fsq_calc/fsq_obs
    
    fabs = abs(fsq_obs-fsq_calc)/float(sigfsq)

    hkl_all.append([h,k,l,fsq_calc,fsq_obs,sigfsq,foc,fabs, raw_index])
    if (r_foc>0 and foc>r_foc) or  (r_foc>0 and fco>r_fco) or (r_fabs>0 and fabs>r_fabs):
        seqnum = raw_hkl[raw_index][15]
        hkl_omit.append([h,k,l,fsq_calc,fsq_obs,sigfsq,foc,fco,fabs, seqnum, raw_index+1])     
    else:
        iHKL = raw_hkl[raw_index]
        # output reflections to hkl file                
        hkl_output.write( 3*'%4d' % (iHKL[0],iHKL[1],iHKL[2])
            + 2*'%8.2f' % (iHKL[3],iHKL[4])
            + '%4d' % iHKL[5] 
            + 2*'%8.5f' % (iHKL[6], iHKL[7])
            + 6*'%9.5f' % (iHKL[8],iHKL[9],iHKL[10],iHKL[11],iHKL[12],iHKL[13])
            + '%6d' % (iHKL[14])
            + '%7d' % (iHKL[15])
            + '%7.4f%4d%9.5f%8.4f' % (iHKL[16],iHKL[17],iHKL[18],iHKL[19])
            + 2*'%7.2f' % (iHKL[20], iHKL[21])
            + '\n' )
    hkl_used.append(iHKL)                      
    raw_index += 1  
        #hkl_output.write('%4d%4d%4d%8.2f%8.2f%4d\n' \
        #    % (h, k, l, fsq_obs, sigfsq, iScale)

# last record all zeros for shelx
iz = 0
rz = 0.0
hkl_output.write( 3*'%4d' % (iz, iz, iz)
    + 2*'%8.2f' % (rz, rz)
    + '%4d' % iz 
    + 2*'%8.5f' % (rz, rz)
    + 6*'%9.5f' % (rz, rz, rz, rz, rz, rz)
    + '%6d' % (iz)
    + '%7d' % (iz)
    + '%7.4f%4d%9.5f%8.4f' % ( rz, iz, rz, rz )
    + '%7.2f%7.2f' % ( rz, rz )
    + '\n' )
    
hkl_input.close()

hkl_output.close()


#Save omitted reflections to file
flist = 'fcf_2hkl.lst'
flist_output = open(flist,'w')
flist_output.write('  H   K    L    F2calc     F2obs    SigF2  F2obs/F2calc  F2calc/F2obs |F2obs-F2calc|/SigF2 Seqnum. Line No.\n') 
for HKL in hkl_omit:
    flist_output.write( 3*'%4d' % (HKL[0],HKL[1],HKL[2])
                         + 2*'%10.2f' % (HKL[3],HKL[4])
                         + '%9.2f' %(HKL[5])
                         + 2*'%12.2f  ' %(HKL[6],HKL[7])
                         + '%16.2f     ' %(HKL[8])
                         + '%7d' %(HKL[9])
                         + '%9d' %(HKL[10])
                         + '\n')

flist_output.close()
                        
print('\n***   %g reflections were rejected out of a total of  %g. ***\n'%(len(hkl_omit),len(hkl_all)))

#
a =numpy.array(hkl_used)
print('number of peaks %g\n'%(len(a)))
# Output cif for TOPAZ experiment
cif = 'fcf4_hklf2.cif'
cif_output = open(cif, 'w')
cif_output.write('# TOPAZ CIF updates after outlier removal\n'
    + '# Cut-off criteria for outliers:\n'
    + '# MAXIMUM FOSQ/FCSQ RATIO = %d\n' % (r_foc)
    + '# MAXIMUM FCSQ/FOSQ RATIO = %d\n' % (r_fco)
    + '# MAXIMUM ABS(FOSQ - FCSQ)/SIGFOSQ RATIO = %d\n' % (r_fabs)
    + '# MINIUM D-SPACING %7.3f Angs'%(numpy.min(a[:,19]))
    + '\n'
    + 'data_topaz \n'
    + ' \n'
    + '_diffrn_ambient_temperature          ?\n'
    + '_diffrn_radiation_wavelength         "%4.2f-%5.2f"\n'%(numpy.min(a[:,6]),numpy.max(a[:,6]))
    + '_diffrn_radiation_wavelength_details  "Based on neutron time of flight"\n'
    + '_diffrn_radiation_type               neutrons\n'
    + '_diffrn_source                       "The ORNL Spallation Neutron Source"\n'
    + '_diffrn_measurement_device_type      TOPAZ\n'
    + '_diffrn_measurement_method           "time-of-flight Laue"\n'
    + '_diffrn_detector_area_resol_mean     ?\n'
    + '\n'
    + '_diffrn_reflns_theta_min         %7.3f\n'%(numpy.min(a[:,18])*180.0/math.pi/2.0)
    + '_diffrn_reflns_theta_max         %7.3f\n'%(numpy.max(a[:,18])*180.0/math.pi/2.0)
    + ' \n'
    + '_computing_data_collection        "SNS PyDas"\n'
    + '_computing_cell_refinement         Mantidplot\n'
    + '_computing_data_reduction          Mantidplot\n')

cif_output.close()

print('The End')

    
        
