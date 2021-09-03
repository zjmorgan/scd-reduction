
# File: topaz_reduction_mod_combinefiles.py
#
# Modified and combined python scripts from 
# lin_avs_coef2.py, ReduceSCD_Parallel.py and anvred3_topi.py
# 
# X. P. Wang, May 2018
#
# =========================================================================================
# Part One 
#   Calculate sample linear absorption coefficients and density lin_avs_coef2.py
#
# Part Two 
#  Peak Integration to generate  individual and combined .integrate and .mat files
#
# Part Three
#  Peak intensity normalization to generate .hkl data files
#
# =========================================================================================
#
import os
import sys
import threading
import time
import math
import numpy as np
from operator import itemgetter
from itertools import groupby

sys.path.insert(0,"/opt/mantidnightly/bin")
sys.path.insert(0,"/opt/mantidnightly/lib")

import ReduceDictionary

from mantid.simpleapi import *
from mantid.kernel import *
from mantid.geometry import PointGroupFactory, SpaceGroupFactory, UnitCell
#print '\n'.join(sys.path)

if os.path.exists('/SNS/TOPAZ/shared/PythonPrograms/Python3Library'):
    sys.path.append('/SNS/TOPAZ/shared/PythonPrograms/Python3Library')
    sys.path.append('/SNS/TOPAZ/shared/PythonPrograms/ISAW_PythonSources/Lib')

from readrefl_header import *
from readrefl_SNS import *
from readSpecCoef import *
from spectrumCalc import *
from spectrum2 import *
from absor_sphere import *


print("API Version")
print(apiVersion())

start_time = time.time()

# Get the config file name from the command line
#
if (len(sys.argv) < 2):
  print("You MUST give the config file name on the command line")
  exit(0)

config_file_name = sys.argv[1]

#
# Load the parameter names and values from the specified configuration file 
# into a dictionary and set all the required parameters from the dictionary.
#

params_dictionary = ReduceDictionary.LoadDictionary( config_file_name )

formulaString         = params_dictionary[ "formulaString" ]
zParameter            = float(params_dictionary[ "zParameter" ])
unitCellVolume        = float(params_dictionary[ "unitCellVolume" ])
nistdatafname         = params_dictionary["nistdatafname"]
calcRadius            = params_dictionary[ "calcRadius" ]
weight                = float(params_dictionary.get("weight",'0.0'))
sampleRadius          = float(params_dictionary.get("sampleRadius",'1.0'))

exp_name              = params_dictionary[ "exp_name" ]
output_directory      = params_dictionary[ "output_directory" ]
reduce_one_run_script = params_dictionary[ "reduce_one_run_script" ]
slurm_queue_name      = params_dictionary[ "slurm_queue_name" ] 
max_processes         = int(params_dictionary[ "max_processes" ])
min_d                 = params_dictionary[ "min_d" ]
max_d                 = params_dictionary[ "max_d" ]
tolerance             = params_dictionary[ "tolerance" ]
cell_type             = params_dictionary[ "cell_type" ] 
centering             = params_dictionary[ "centering" ]
run_nums              = params_dictionary[ "run_nums" ]

use_cylindrical_integration = params_dictionary[ "use_cylindrical_integration" ]
instrument_name       = params_dictionary[ "instrument_name" ]

read_UB               = params_dictionary[ "read_UB" ]
UB_filename           = params_dictionary[ "UB_filename" ]

iSpec                  = int(params_dictionary.get("iSpec",'0'))
specCoeffFile          = params_dictionary.get("specCoeffFile",str(''))
spectraFile            = params_dictionary["spectraFile"]
normToWavelength       = float(params_dictionary["normToWavelength"])
minIsigI               = float(params_dictionary["minIsigI"])
numBorderCh            = int(params_dictionary["numBorderCh"])
intiMin                = float(params_dictionary["intiMin"])
dMin                   = float(params_dictionary["dMin"])
iIQ                    = int(params_dictionary.get("iIQ",'1'))
scaleFactor            = float(params_dictionary["scaleFactor"])
wlMin                  = float(params_dictionary["wlMin"])
wlMax                  = float(params_dictionary["wlMax"])
abs_correc_type        = params_dictionary["abs_correc_type"]
pg_symbol              = params_dictionary["pg_symbol"]
z_score                = float(params_dictionary["z_score"])
starting_batch_number  = int(params_dictionary.get("starting_batch_number",'1'))

max_order              = params_dictionary.get('max_order', "0")
cross_terms            = params_dictionary.get("cross_terms",'False')
save_mod_info          = params_dictionary.get("save_mod_info",'False')
mod_vector1            = params_dictionary.get("mod_vector1", '0,0,0')
mod_vector2            = params_dictionary.get("mod_vector2", '0,0,0')
mod_vector3            = params_dictionary.get("mod_vector3", '0,0,0')
tolerance_satellite    = params_dictionary.get("tolerance_satellite", '0.10')
commonUB_for_all       = params_dictionary.get("commonUB_for_all",False)

reduce_one_run_script  = os.path.dirname(os.path.realpath(__file__)) +'/' + reduce_one_run_script

if not os.path.exists(output_directory):
    os.makedirs(output_directory)
print("\nWorking dir : %s" % os.getcwd())
os.chdir(output_directory)
print("\nOutput  dir : %s" % output_directory)

#*****************Part One *************************************************
#
#           lin_abs_coef v3.
#***************************************************************************

# Program to calculate linear absorption coefficients and density.
# Version 1 requires ISAW. Version 2 is a stand alone script.

# Version 2:
# A.J. Schultz, August 2015
# Version 3:
# X. P. Wang, May 2018
#

smu =0.0
amu = 0.0
radius = 0.0
ub_matrix_file  =''
anvred_integrate_fname =''

formulaList = formulaString.split()
numberOfIsotopes = len(formulaList)     # the number of elements or isotopes in the formula        

sumScatXs = 0.0
sumAbsXs = 0.0
sumAtWt = 0.0

logFileName = 'lin_abs_coef.log'
logFile = open( logFileName, 'w' )
logFile.write('Output from lin_abs_coef.py script:\n\n')

logFile.write('Chemical formula: ' + formulaString + '\n')
logFile.write('Number of formula units in the unit cell (Z): %6.3f\n' % zParameter)
logFile.write('Unit cell volume (A^3): %8.2f\n' % unitCellVolume)

logFile.write('\nCross sections in units of barns ( 1 barn = 1E-24 cm^2)\n')
logFile.write('Absorption cross section for 2200 m/s neutrons (wavelength = 1.8 A)\n')
logFile.write('For further information and references, see .../SNS/TOPAZ/IPTS-9918/shared/2016A/reduction/\n')

print('\nAtom      ScatXs      AbsXs')	# print headings
print('----      ------      -----')
logFile.write('\nAtom      ScatXs      AbsXs\n')
logFile.write(  '----      ------      -----\n')

# Except for hydrogen, cross-section values are from the NIST web site:
# http://www.ncnr.nist.gov/resources/n-lengths/list.html
# which are from:
# V. F. Sears, Neutron News, Vol. 3, No. 3, 1992, pp. 29-37.
# Hydrogen cross-sections are from:
# Howard, J. A. K.; Johnson, O.; Schultz, A. J.; Stringer, A. M.
#	J. Appl. Cryst. 1987, 20, 120-122.
        
#filename = XsecDirectory + 'NIST_cross-sections.dat'
nistdatafname = '/SNS/TOPAZ/shared/calibrations/NIST_cross-sections.dat'

# begin loop through each atom in the formula
for i in range(numberOfIsotopes):

    lenAtom = len(formulaList[i])   # length of symbol plus number in formula
    
    # begin test for number of characters in the isotope name or symbol
    for j in range(lenAtom):          
        lenSymbol = lenAtom - j - 1
        if formulaList[i][lenSymbol].isalpha(): break
    lenSymbol = lenSymbol + 1
    
    input = open(nistdatafname, 'r')         # this has the effect of rewinding the file
    lineString = input.readline()       # read the first comment line
    while lineString[0] == '#':         # search for the end of the comments block
        lineString = input.readline()

    # Begin to search the table for element/isotope match.
    
    lineList = lineString.split()       # this should be the H atom

    while formulaList[i][0:lenSymbol] != lineList[0]:
        lineString = input.readline()
        lineList = lineString.split()

    scatteringXs = float(lineList[1])   # the total scattering cross section
    absorptionXs = float(lineList[2])   # the true absorption cross section at 1.8 A
    atomicWeight = float(lineList[4])   # atomic weight
    number = float(formulaList[i][lenSymbol:])   # the number of this nuclei in the formula
    
    print('%-5s %10.5f %10.5f' % (lineList[0], scatteringXs, absorptionXs))
    logFile.write('%-5s %10.5f %10.5f\n' % (lineList[0], scatteringXs, absorptionXs))
    
    sumScatXs = sumScatXs + ( number * scatteringXs )
    sumAbsXs = sumAbsXs + ( number * absorptionXs )
    sumAtWt = sumAtWt + ( number * atomicWeight )
    input.close()
# end loop

print('\nMolecular weight for %s: %10.2f'%(formulaString,sumAtWt))

# Calculate the linear absorption coefficients in units of cm^-1
muScat = sumScatXs * zParameter / unitCellVolume
muAbs = sumAbsXs * zParameter / unitCellVolume

# Calculate the density of the crystal in g/cc
density = (sumAtWt / 0.6022) * zParameter / unitCellVolume

# Print the results.
print('\n')
print('The linear absorption coefficent for total scattering is %6.3f cm^-1' % muScat)
print('The linear absorption coefficent for true absorption at 1.8 A is %6.3f cm^-1' % muAbs)
print('The calculated density is %6.3f grams/cm^3' % density)
logFile.write('\n')
logFile.write('The linear absorption coefficent for total scattering is %6.3f cm^-1\n' % muScat)
logFile.write('The linear absorption coefficent for true absorption at 1.8 A is %6.3f cm^-1\n' % muAbs)
logFile.write('\nThe calculated density is %6.3f grams/cm^3\n' % density)

# if calcRadius:
if weight != 0.0:
    crystalVolume = weight / (density)   # sample volume in mm^3
    print('For a weight of %6.3f mg, the crystal volume is %6.3f mm^3' % (weight, crystalVolume))
    logFile.write('\nFor a weight of %6.3f mg, the crystal volume is %6.3f mm^3\n' % (weight, crystalVolume))
    crystalRadius = ( crystalVolume / ((4.0/3.0)*math.pi) )**(1.0/3.0)   # radius in mm
    print('The crystal radius is %6.3f mm, or %6.4f cm' % (crystalRadius, crystalRadius/10.))
    logFile.write('The crystal radius is %6.3f mm, or %6.4f cm\n' % (crystalRadius, crystalRadius/10.))
    # volCalc = (4.0/3.0) * math.pi * crystalRadius**3
    # print 'volCalc = %6.3f' % volCalc
else:
    print('The crystal radius from input is %6.3f mm, or %6.4f cm' % (sampleRadius, sampleRadius/10.))
    logFile.write('The crystal radius from input is %6.3f mm, or %6.4f cm\n' % (sampleRadius, sampleRadius/10.))

    weight = density *(4.0/3.0)*math.pi * (sampleRadius)**(3.0)
    print('For a radius of %6.3f mm, the crystal weight is %6.3f mg' % (sampleRadius, weight))
    logFile.write('\nFor a radius of %6.3f mm, the crystal weight is %6.3f mg' % (sampleRadius, weight))

logFile.close()
print('\nCompleted calculation of linear absorption coefficients\n')

#

print("\n*********************************************************************************")
# **************************Part Three ***********************************************************
# Modified Version 2.0 ReduceSCD_CombineFiles.py
#
# First combine all of the integrated files, by reading the separate files and
# appending them to a combined output file.
#
# X. P. Wang May 2018
#
print("\n**************************************************************************************")
print("************** Starting to Combine Results *******************************************")
print("**************************************************************************************\n")

niggli_name = output_directory + "/" + exp_name + "_Niggli"
niggli_integrate_file = niggli_name + ".integrate"
niggli_matrix_file = niggli_name + ".mat"

if not use_cylindrical_integration:
  first_time = True
  seqNum = 0
  fout = open(niggli_integrate_file,"w")
  for r_num in run_nums:
    one_run_file = output_directory + '/' + str(r_num) + '_Niggli.integrate'
    f = open(one_run_file,"r")
    lines = f.readlines()
    f.close()
    for line in lines:
      if line[0] =="0" or line[0] == "1" or line[0] == "2":
        fout.write(line)
      elif line[0] !="4" and line[0] != "5" and line[0] != "6" and line[0] != "7" and line[0] != "V":
        strCount = "%6d"%seqNum
        line = line[:2] + strCount + line[8:]
        fout.write(line)
        seqNum = seqNum + 1
      elif first_time:
        fout.write(line)
    first_time = False
  fout.close()

#
# Load the combined file and re-index all of the peaks together. 
# Save them back to the combined Niggli file (Or selcted UB file if in use...)
#
  peaks_ws = LoadIsawPeaks( Filename=niggli_integrate_file )

#
# Find a Niggli UB matrix that indexes the peaks in the combined peaks workspace
# Option to load UB instead of using FFT
# Index peaks using UB from UB of initial orientation run/or combined runs from first iteration of crystal orientation refinement

  if read_UB:
    LoadIsawUB(InputWorkspace=peaks_ws, Filename=UB_filename)
    #Refine UB and reindex peaks
    IndexPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance, ToleranceForSatellite=tolerance_satellite,CommonUBForAll=False, RoundHKLs=False)
    FindUBUsingIndexedPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance)
    #OptimizeCrystalPlacement(PeaksWorkspace=peaks_ws,ModifiedPeaksWorkspace=peaks_ws,
    #  FitInfoTable='CrystalPlacement_info',MaxIndexingError=tolerance)

  else:
      FindUBUsingFFT( PeaksWorkspace=peaks_ws, MinD=min_d, MaxD=max_d, Tolerance=tolerance,  Iterations=30)

  #Use q-vector from user input  
  if save_mod_info is True: 
    print('max_order :', str(max_order)) 
    IndexPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance, ToleranceForSatellite=tolerance_satellite, 
                  RoundHKLs=False, 
                  ModVector1=mod_vector1, ModVector2=mod_vector2, ModVector3=mod_vector3, 
                  MaxOrder = max_order,
                  CrossTerms=True, 
                  SaveModulationInfo=True)

# test to see if peaks have satellites
  if peaks_ws.sample().hasOrientedLattice():
      maxOrder = peaks_ws.sample().getOrientedLattice().getMaxOrder()
  else:
      maxOrder = int(max_order)

  if maxOrder > 0:
    IndexPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance, ToleranceForSatellite=tolerance_satellite,CommonUBForAll=False, RoundHKLs=False)

  #Use the common UB to index peaks in the combined peaks workspace
  IndexPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance, ToleranceForSatellite=tolerance_satellite,CommonUBForAll=commonUB_for_all, RoundHKLs=False)
  SaveIsawPeaks( InputWorkspace=peaks_ws, AppendFile=False, Filename=niggli_integrate_file, RenumberPeaks=True)
  SaveIsawUB( InputWorkspace=peaks_ws, Filename=niggli_matrix_file )

  anvred_integrate_fname = niggli_integrate_file
  ub_matrix_file = niggli_matrix_file

#
# If requested, also switch to the specified conventional cell and save the
# corresponding matrix and integrate file
#
if not use_cylindrical_integration:
  if (not cell_type is None) and (not centering is None) :
    conv_name = output_directory + "/" + exp_name + "_" + cell_type + "_" + centering
    conventional_integrate_file = conv_name + ".integrate"
    conventional_matrix_file = conv_name + ".mat"
    SelectCellOfType( PeaksWorkspace=peaks_ws, CellType=cell_type, Centering=centering,
                      Apply=True, Tolerance=tolerance, AllowPermutations=True )
    #Do IndexPeaks first
    IndexPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance, ToleranceForSatellite=tolerance_satellite, RoundHKLs=False)
    FindUBUsingIndexedPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance)

    # Least squares refinement of conventional cell after transfromation
    if (cell_type=='Rhombohedral' and centering=='R'):
      cell_type = 'Hexagonal'  # The R lattice is hexagonal after unit cell transformation
    #elif (cell_type=='Monoclinic'):
    #  cell_type = 'Monoclinic ( b unique )'
    OptimizeLatticeForCellType(PeaksWorkspace=peaks_ws, CellType=cell_type,Apply='1', Tolerance=tolerance)
#    OptimizeCrystalPlacement(PeaksWorkspace=peaks_ws,ModifiedPeaksWorkspace=peaks_ws,
#      FitInfoTable='CrystalPlacement_info',MaxIndexingError=tolerance)
    SaveIsawPeaks( InputWorkspace=peaks_ws, AppendFile=False, Filename=conventional_integrate_file,RenumberPeaks=True)
    SaveIsawUB( InputWorkspace=peaks_ws, Filename=conventional_matrix_file )
    #ANVRED
    anvred_integrate_fname = conventional_integrate_file
    ub_matrix_file = conventional_matrix_file

if use_cylindrical_integration: 
  if (not cell_type is None) or (not centering is None):
    print("WARNING: Cylindrical profiles are NOT transformed!!!")
  # Combine *.profiles files
  filename = output_directory + '/' + exp_name + '.profiles'
  output = open( filename, 'w' )

  # Read and write the first run profile file with header.
  r_num = run_nums[0]
  filename = output_directory + '/' + instrument_name + '_' + r_num + '.profiles'
  input = open( filename, 'r' )
  file_all_lines = input.read()
  output.write(file_all_lines)
  input.close()
  os.remove(filename)

  # Read and write the rest of the runs without the header.
  for r_num in run_nums[1:]:
      filename = output_directory + '/' + instrument_name + '_' + r_num + '.profiles'
      input = open(filename, 'r')
      for line in input:
          if line[0] == '0': break
      output.write(line)
      for line in input:
          output.write(line)
      input.close()
      os.remove(filename)

  # Remove *.integrate file(s) ONLY USED FOR CYLINDRICAL INTEGRATION!
  for file in os.listdir(output_directory):
    if file.endswith('.integrate'):
      os.remove(file)



print("\n**************************************************************************************")
print("****************************** DONE PROCESSING ALL RUNS ******************************")
print("**************************************************************************************\n")

print("\n**************************************************************************************")
print("****************************** ANVRED Reduction **************************************")
print("**************************************************************************************\n")

# Modified anvred3.py
# Ouptput topi for JANA2006
# Outlier removal based on crystal symmetry statistics
# Use integrated config file
#
# X. P. Wang, May, 2018
#--------------------------------------------------------------------
#                             anvred3.py
#
#
# Includes the option to perform a spherical absorption correction
# or calculate direction cosines for a polyhedral correction
# using the gaussian.f program.
#
# Includes a Tkinter gui interface.
#
#     A. J. Schultz, May 2014
#
#--------------------------------------------------------------------
#                             anvred2x.py
#
# Does not run from Isaw. User input read from anvred2x.inp file.
#   A. J. Schultz, February 2012
#--------------------------------------------------------------------
#
# Data reduction program:
#   Input is raw integrated intensities.
#   Output is relative Fsq's.
#
# Jython version:
#    A. J. Schultz, started December 2009
#
# anvred_py.py
#    Each spectrum is a separate file. The spectra files are created
#    by "TOPAZ_spectrum_multiple_banks.iss".
#
# anvred2.py:
#    This version reads one spectrum file containing spectra for
#    each detector. The spectra are created by "TOPAZ_spectrum.py".
#
# Modfications by Xiaoping Wang, April 2011
# Added Selection of neutron wavelengths limits wlMin, wlMax
# Omit zero intensity peaks in integrate file XP Wang 03/21/2011
# Changed to >=0 and used absolute value for minium I/sing(I) = 0  XP Wang 02/24/2011
# Added detector scale factors for vanadium/niobium spectrum XP Wang 09/24/2013
# Changed output format for seqnum from %6d to %7d, 05/24/2015
#
#
# 
# Comments from Fortran source:
# C**************************   ANVRED  ******************************
# C
# C ARGONNE NATIONAL LABORATORY VARIABLE WAVELENGTH DATA REDUCTION PROGRAM
# C
# C		Major contributions from:
# C			P. C. W. Leung
# C			A. J. Schultz
# C			R. G. Teller
# C			L. R. Falvello
# C
# C     The data output by this program are corrected  for  variations  in
# C  spectral distribution, variations in detector efficiency  across  the
# C  face  of  the  detector,  and the pertinent geometric factors such as
# C  (SIN(THETA))**2 and LAMBDA**4.
# C
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !	Linux version:	A. Schultz   January, 2003                    !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# !	Version 4:		A. Schultz		June, 2003
# !		1. Process data from two detectors
# !		2. Does not use an x-file.
# !		3. Gets MONCNT from integrate file. If zero, sets CMONX = 1.0.
# !		4. Corrected ALPHAP for calculation of SPECT1.
#
# !	Version 5:		A. Schultz		July, 2003
# !		This version outputs a expnam.hkl file which can be input
# !		into SHELX with HKL 2.
# !	Version 5a:
# !		Cleaned-up and removed a lot of unused code.
# !		Added a test for dmin.
# !
# !	Version 6:		L. Falvello		January, 2004
# !		Polyhedral absorption correction with two detectors.
# !
# !	Version 7:		A. Schultz		2007
# !		Use spectrum obtained from each SCD detector
# !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !	ANVRED_SNS:		A. Schultz		2008                                     !
# !		Process SNS data.                                                    !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !	ANVRED_SNS_v2:		A. Schultz		June, 2008
# !		New spherical absorption correction. Removed all
# !		of the old correction code.
# !	ANVRED_SNS-v2.1: read detector parameters from integrate file.
# !	ANVRED_SNS-v2.2: get filename for spectrum file.
# !       ANVRED_SNS-v2.3: everything included in one file.  8/17/2008
# !       anvredSNS_2.4: made compatible with gfortran including removal
# !                      of FREIN3 and READ133.         10/8/2008
# !	anvredSNS_2.5: the datacom_SNS.inc file is no longer used. Cleaned
# !			up the code using ftnchek.    10/13/08
# !	anvredSNS_2.6: assign a common scale factor for each
# !			crystal setting, or for each detector. 1/29/09
# !
# !	4/13/09	Number of possible spectra increased from 2 to 100.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 
def huq(h, k, l, UB):
    "Multiplies hkl times UB matrix to return q-vector"
    hh = [h, k, l]
    q = np.zeros((3))
    q = np.dot(hh, UB)
    return q

def rotate_matrix(UB, omega, chi, phi, SNS_or_IPNS):
    "Rotates UB matrix by setting angles"
    fmat = rotation_matrix(omega, chi, phi, SNS_or_IPNS)
    newmat = np.zeros((3,3))
    newmat = np.dot(UB, fmat)
    return newmat

def rotation_matrix(omega, chi, phi, SNS_or_IPNS):
    "Returns rotation matrix from setting angles"
    rad = 180. / math.pi

    ph = phi / rad
    cp = math.cos(ph)
    sp = math.sin(ph)
    R_phi = np.zeros((3,3))
    R_phi[0,0] = cp
    R_phi[0,1] = sp
    R_phi[1,0] = -sp
    R_phi[1,1] = cp
    R_phi[2,2] = 1.0

    ch = chi / rad        #changed -chi to chi, 8/23/07
    cc = math.cos(ch)
    sc = math.sin(ch)
    R_chi = np.zeros((3,3))
    R_chi[0,0] = 1.0
    R_chi[1,1] = cc
    R_chi[1,2] = sc
    R_chi[2,1] = -sc
    R_chi[2,2] = cc

    if SNS_or_IPNS == 'IPNS':
        om = -omega / rad   # for IPNS data set omega to -omega
    if SNS_or_IPNS == 'SNS':
        om = omega / rad      # for SNS data
    co = math.cos(om)
    so = math.sin(om)
    R_om = np.zeros((3,3))
    R_om[0,0] = co
    R_om[0,1] = so
    R_om[1,0] = -so
    R_om[1,1] = co
    R_om[2,2] = 1.0

    fmat = np.zeros((3,3))
    fmat = np.dot(R_phi, R_chi)
    fmat = np.dot(fmat, R_om)

    return fmat

def Rvec(twoth, az, L2):
    "Return R vectors for a peak in peaks file."
    
    R_IPNS = np.zeros((3))
    R_SNS = np.zeros((3))
    
    # IPNS axes with x as the beam direction and z is vertically upward
    R_IPNS[0] = math.cos(twoth) * L2
    R_IPNS[1] = math.cos(az) * math.sin(twoth) * L2
    R_IPNS[2] = math.sin(az) * math.sin(twoth) * L2
    
    # SNS axes with z as the beam direction and y is vertically upward
    R_SNS[0] = math.cos(az) * math.sin(twoth) * L2
    R_SNS[1] = math.sin(az) * math.sin(twoth) * L2
    R_SNS[2] = math.cos(twoth) * L2
    
    return R_IPNS, R_SNS
    

#------------------ Begin ------------------

smu = muScat 
amu = muAbs 
radius = sampleRadius / 10.0  #Crystal radius in cm

# Check that one of the two absorption correction types has been selected
if abs_correc_type != 'spherical' and abs_correc_type != 'polyhedral':
    print('')
    print('**********************************************************')
    print("Absorption correction type is not 'spherical' or 'polyhedral'.")
    print('Check your spelling.')
    print('**********************************************************')
    print('')

# Read UB matrix if polyhedral absorption
# if abs_correc_type == 'polyhedral':
if True:
    # Open matrix file
    UB_input = open(ub_matrix_file,'r')

    # Initialize UB_IPNS matrix
    UB_IPNS = np.zeros((3,3))
    print('\n Input from matrix file ' + ub_matrix_file + ':\n')

    # Read matrix file into UB_IPNS matrix
    for i in range(3):
        linestring = UB_input.readline()
        print(linestring.strip('\n'))
        linelist = linestring.split()
        for j in range(3):
            UB_IPNS[i,j] = float(linelist[j])           
    # Read next 2 lines containing lattice constants and esd's
    for i in range(2):
        linestring = UB_input.readline()
        print(linestring.strip('\n'))
    print('\n')
    # End of reading and printing matrix file
                
#
#detScale={13:1.0,14:1.0,16:1.0,17:1.0,18:1.0,19:1.0,\
#  20:1.0,21:1.0,22:1.0,23:1.0,26:1.0,27:1.0,28:1.0,29:1.0,\
#  33:1.0,34:1.0,36:1.0,37:1.0,38:1.0,39:1.0,\
#  46:1.0,47:1.0,48:1.0,49:1.0,
#  56:1,57:1,58:1,59:1}                               

#Scolecite 2013B, XP Wang Sept 23, 2013 
#detScale = {17:1.115862021,18:0.87451341,\
#      22:1.079102931,26:1.087379072,27:1.064563992,28:0.878683269, \
#      36:1.15493377,37:1.010047685,38:1.046416037,39:0.83264528, \
#      47:1.06806776,48:0.872542083,\
#      58:0.915242691}

#Scolecite 2015/08/31
#detScale={13:1.044824,14:1.06399,16:0.95552,17:1.09583,18:0.84606,19:1.00154,\
# 22:1.015696,23:1.05283,26:1.05720,27:1.04025,28:0.84767,29:1.03148,\
#  33:1.101041,34:0.98625,36:1.02813,37:1.10629,38:0.87417,39:0.86848,\
# 46:0.897062,47:1.23349,48:0.82479,49:1.02740,\
#  58:0.921672}
## Scolecite 09/02/2015 Use New V/Nb spectrum TOPAZ_12387_12389
#detScale={13:1.029577,14:1.03623,16:0.93879,17:1.07941,18:0.83849,19:1.02553,\
#  22:1.000527,23:1.03267,26:1.03788,27:1.02859,28:0.85372,29:1.10750,\
#  33:1.082904,34:0.96428,36:1.12968,37:0.97101,38:1.00890,39:0.86381,\
#  46:0.979966,47:1.02720,48:0.85466,49:1.10869,\
#  58:0.98899}
#Si test crystal 10/14/2015 Rev. Use New V/Nb spectrum TOPAZ_13383_13384
#detScale={13:1.068455,14:1.00000,16:1.04368,17:1.15773,18:0.83238,19:0.95370,\
#        20:1.00000,22:0.998234,23:1.10762,26:1.13599,27:1.07148,28:0.81808,29:0.93032,\
#        33:1.120385,34:1.03951,36:1.22393,37:1.02077,38:0.99489,39:0.75315,\
#        46:1.057146,47:1.02733,48:0.81607,49:0.82916}

#Si test crystal 04/11/2016
#detScale={13:1.110083,14:1.230529,16:1.05020,17:1.15631,18:0.84376,19:0.95072,\
#      20:0.692254,22:1.043647,23:1.11743,26:1.15772,27:1.10177,28:0.80728,29:0.85565,\
#      33:1.120911,34:1.07563,36:1.20957,37:1.03539,38:1.00286,39:0.71187,\
#      46:1.092741,47:1.05848,48:0.81366,49:0.76153}
#No Scale
#detScale={13:1.0,14:1.0,16:1.0,17:1.0,18:1.0,19:1.0,\
#      20:1.0,22:1.0,23:1.0,26:1.0,27:1.0,28:1.0,29:1.0,\
#      33:1.0,34:1.0,36:1.0,37:1.0,38:1.0,39:1.0,\
#      46:1.0,47:1.0,48:1.0,49:1.0}

#Scolecite 04/12/2016
#detScale={13:1.115765,14:1.229136,16:1.05960,17:1.15974,18:0.84514,19:0.94292,\
#     20:0.679666,22:1.046624,23:1.12315,26:1.16409,27:1.10448,28:0.80435,29:0.83439,\
#     33:1.121278,34:1.09044,36:1.21097,37:1.03847,38:1.00295,39:0.70865,\
#     46:1.102081,47:1.06180,48:0.81028,49:0.74404}

#No Scale
#detScale={13:1.0,14:1.0,16:1.0,17:1.0,18:1.0,19:1.0,\
#      20:1.0,21:1.0, 22:1.0,23:1.0,26:1.0,27:1.0,28:1.0,29:1.0,\
#      33:1.0,34:1.0,36:1.0,37:1.0,38:1.0,39:1.0,\
#      46:1.0,47:1.0,48:1.0,49:1.0}

#DetScale Scolecite	08-13-16															
#detScale={
#        13:1.090472,14:1.27426,16:1.00603,17:1.15163,18:0.84610,19:0.93433,\
#        20:0.705760,21:0.81827,22:1.05983,23:1.11589,26:1.12619,27:1.11085,28:0.80871,29:0.88591,\
#        33:1.155496,34:1.03284,36:1.23397,37:1.04722,38:1.03315,39:0.75703,\
#        46:1.077561,47:1.09341,48:0.82596,49:0.80910
#          }

#Det Scale Si Crystal 2017A 2-3 mm BN Aperture
#detScale={13:1.025627,14:1.246955,16:1.03047,17:1.14505,18:0.84294,19:1.01058,\
#        20:0.809640,21:0.816426,22:1.00095,23:1.07553,26:1.10310,27:1.09329,28:0.82710,29:1.00269,\
#         33:1.073669,34:0.98532, 36:1.15649,37:0.99906,38:1.03150,39:0.81993,\
#         43:1.145198,46:1.025365,47:1.01475,48:0.83253,49:0.88582}

# Scolecite 2017B 2-2 mm BN Aperture July 19.2017
#detScale={13:1.081525,14:1.233261,16:0.99724,17:1.14080,18:0.84160,19:0.93578,\
#         20:0.709277,21:0.816502,22:1.04281,23:1.11220,26:1.12005,27:1.11828,28:0.81039,29:0.86783,\
#          33:1.131051,34:1.02575, 36:1.21839,37:1.04092,38:1.01037,39:0.74954,\
#          43:1.199689,46:1.066310,47:1.08897,48:0.82660,49:0.81487}
                                                                    
#Det Scale,Scolecite 2018B 3-3 mm BN Aperture
#detScale={16:1.192051,17:0.929431,18:0.89680,19:0.87531,\
#         22:1.059484,26:0.935278,27:0.88258,28:0.92069,29:0.86587,\
#         33:1.084642,36:1.13744,37:1.10597,38:0.92006,39:0.89196,\
#         46:1.130662,47:1.00365,48:1.09350,49:1.07460}

# Det Scale, Si 2-3 mm BN Aperture Oct. 1, 2018
#detScale={16:0.985425,17:1.012897,18:0.96796,19:0.88418,\
#        20:0.770015,22:1.108947,26:1.04216,27:0.95799,28:0.93083,29:0.80043,\
#        33:1.183736,36:1.04361,37:1.16065,38:0.94089,39:0.79518,\
#        46:1.281827,47:1.06494,48:1.09645,49:0.97189}

# Scolecite 2-3 mm BN Aperture, Jan. 28, 2019
#detScale={13:1.194192,16:1.009037,17:1.02146,18:0.96004,19:0.87209,\
#          20:0.801821,22:1.179643,26:1.02576,27:0.94814,28:0.92824,29:0.82279,\
#          33:1.098516,36:1.06507,37:1.08936,38:0.97915,39:0.83023,\
#          46:1.277106,47:1.09529,48:0.97434,49:0.84137}

# Scolecite 2-3 mm BN Aperture, Feb. 13, 2019
#detScale={14:1.154117,16:0.982386,17:1.00935,18:0.96693,19:0.88913,\
#          20:0.820713,22:1.225087,26:1.00445,27:0.94385,28:0.94622,29:0.84920,\
#          33:1.087622,36:1.03464,37:1.07357,38:0.98321,39:0.84419,\
#          46:1.254725,47:1.09052,48:0.98651,49:0.86519}

# Si 2mm sphere / 2-3 mm BN Aperture July 1, 2019
#detScale={13:1.00838,14:1.116629,16:1.109081,17:0.97126,18:0.95387,19:0.92162,\
#    20:0.869095,22:1.111812,26:1.15751,27:0.91691,28:0.97282,29:0.81954,\
#    33:1.039446,36:1.18936,37:1.02026,38:0.96910,39:0.89641,\
#    46:1.180095,47:0.92802,48:0.99738,49:0.85607,\
#    56:1.13937 ,57:0.92770,58:0.96573,59:0.96253}

# Si 2mm sphere / 2-3 mm BN Aperture Sept 14, 2019
#detScale={13:1.00540,14:1.137220,16:1.112335,17:0.94446,18:0.90887,19:0.85586,\
#    20:0.950860,22:1.10691,26:1.14553,27:0.88001,28:0.92701,29:0.81353,\
#    33:1.017205,36:1.19528,37:0.99549,38:1.01062,39:0.89886,\
#    46:1.170565,47:0.90727,48:0.99052,49:0.87208,\
#    56:1.208060,57:0.97762,58:0.99840,59:0.97004}

# Si 2 mm sphere / 2-3 mm BN Aperture Nov 7, 2019
#detScale={13:1.06535,14:1.159507,16:1.091353,17:0.94523,18:0.90369,19:0.86009,\
#  20:0.826992,22:1.160543,26:1.14849,27:0.91142,28:0.95817,29:0.77649,\
#  33:1.069831,36:1.18313,37:1.01611,38:0.95483,39:0.83659,\
#  46:1.196279,47:0.94534,48:1.01646,49:0.83405,\
#  56:1.16650,57:0.93950,58:1.01668,59:1.01737}

#Jan 21 2020:Si / 2-3 mm BN Aperture
#detScale={13:1.05253,14:1.179764,16:1.155927,17:0.98914,18:0.95470,19:0.87329,\
#    20:0.765826,22:1.120419,26:1.19930,27:0.92780,28:0.95789,29:0.75228,\
#    33:1.087870,36:1.24715,37:1.02962,38:0.97085,39:0.85218,\
#    46:1.186007,47:0.91658,48:0.95868,49:0.80979,\
#    56:1.16396,57:0.92676,58:0.99084,59:0.93086}

#Jun/05/2020:Nd5Pb3O 130K 3-3 mm BN Aperture:
detScale={13:1.16688,14:1.286851,16:1.200336,17:0.96340,18:0.91786,19:0.77843,\
    20:0.702411,22:1.198750,26:1.25762,27:0.89396,28:0.81221,29:0.63821,\
    33:1.170489,36:1.36634,37:1.10622,38:0.92717,39:0.73494,\
    46:1.280152,47:0.92700,48:0.88473,49:0.63999,\
    56:1.27779,57:1.02779,58:0.94930,59:0.89119} 

#Feb/27/2021 Scolecite AG 3-3 mm BN Aperture
#detScale={13:1.124640,14:1.263114,16:1.33895,17:1.01740,18:0.91066,19:0.75849,\
#          20:0.687544,22:1.162599,26:1.23753,27:0.92761,28:0.83093,29:0.63115,\
#          33:1.149815,36:1.314750,37:1.09043,38:0.91098,39:0.69582,\
#          46:1.300390,47:0.959460,48:0.86798,49:0.60712,\
#          56:1.348560,57:1.041110,58:0.98660,59:0.83638}

#June 6/2021 Bixbyite AG 3-3 mm BN Aperture:
#detScale={13:1.05321,14:1.066755,16:0.973064,17:0.92406,18:0.93731,19:0.90281,\
#          20:0.927286,22:1.195631,26:1.09975,27:0.90359,28:0.95096,29:0.82056,\
#          33:1.145274,36:1.07773,37:1.01502,38:0.97703,39:0.88045,\
#          46:1.150778,47:0.94977,48:1.00642,49:0.84127,\
#          56:1.18183,57:0.95950,58:1.04685,59:1.01310}

# open the anvred.log file in the working directory
fileName = output_directory + '/anvred3.log'
logFile = open( fileName, 'w' )

# open the hkl file in the working directory
hklFileName = os.path.splitext(anvred_integrate_fname)[0] + '.hkl'
if abs_correc_type == 'polyhedral':
    hklFileName = os.path.splitext(anvred_integrate_fname)[0] + '.hkl_no_abs'
hklFile = open( hklFileName, 'w' )

# echo the input in the log file
logFile.write('\n********** anvred **********\n')
logFile.write('\nWorking directory: ' + output_directory)
logFile.write('\nExperiment name: ' + os.path.splitext(os.path.basename(anvred_integrate_fname))[0] + '\n')

logFile.write('\nTotal scattering linear absorption coefficient: %6.3f cm^-1' % smu )
logFile.write('\nTrue absorption linear absorption coefficient: %6.3f cm^-1' % amu )
logFile.write('\nRadius of spherical crystal: %6.3f cm\n' % radius )

# logFile.write('\nIncident spectrum and detector efficiency correction.')
# logFile.write('\n    iSpec = 1. Spectrum fitted to 11 coefficient GSAS Type 2 function')
# logFile.write('\n    iSpec = 0. Spectrum data read from a spectrum file.')
# logFile.write('\niSpec: %i\n' % iSpec)

if iSpec == 1:   # spectrum is fitted to equation with 12 coefficients
    logFile.write('\nFile with spectrum coefficients: ' + specCoeffFile + '\n' )
    
if iSpec == 0:   # spectrum is read as TOF vs. counts
    logFile.write('\nFile with spectra: ' + spectraFile + '\n' )

logFile.write('\nNormalize spectra to a wavelength of %4.2f' % normToWavelength)
logFile.write('\nThe minimum I/sig(I) ratio: %i' % minIsigI )
logFile.write('\nWidth of border: %i channels' % numBorderCh )
logFile.write('\nMinimum integrated intensity: %i' % intiMin )
logFile.write('\nMinimum d-spacing : %4.2f Angstroms\n' % dMin )

# logFile.write('\nScale factor identifier:' )
# logFile.write('\n     IQ = 1. Scale factor per crystal setting.' )
# logFile.write('\n     IQ = 2. Scale factor for each detector in each setting.')
# logFile.write('\n     IQ = 3. Scale factor for each detector for all settings.')
# logFile.write('\nIQ: %i\n' % iIQ )

logFile.write('\nMultiply FSQ and sig(FSQ) by: %f\n' % scaleFactor )

logFile.write('\nMinimum wavelength: %f\n' % wlMin )
logFile.write('Maximum wavelength: %f\n' % wlMax )

logFile.write( '\n***** absorption correction type: ' + abs_correc_type + '\n' )
if abs_correc_type == 'polyhedral':
    logFile.write( '\nUB matrix file: ' + ub_matrix_file + '\n' )

# C
# C  CHECK ON THE EXISTANCE OF THE integrate FILE
# C

# anvred_integrate_fname = output_directory + '/' + expName + '.integrate'
integFile = open(anvred_integrate_fname, 'r')

# !  Initial read of integrate file to get instrument and detectors calibration.
calibParam = readrefl_header( integFile )
L1 = float(calibParam[0])       # initial flight path length in cm
t0_shift = float(calibParam[1]) # t-zero offest in microseconds
nod = int(calibParam[2])    # number of detectors
print('********** nod = ', nod)

logFile.write('\nInitial flight path length: %10.4f cm' % L1 )
logFile.write('\nT-zero offset: %8.3f microseconds' % t0_shift )
logFile.write('\nNumber of detectors: %i' % nod )

# Initial values.
transmin = 1.0
transmax = 0.0
hom = 0.39559974    # Planck's constant divided by neutron mass

# Read spectrum coefficients if iSpec = 1
if iSpec == 1:
    specInput = open( specCoeffFile, 'r')
    # pj is a list of lists with dimensions (nod, 11)
    pj = readSpecCoef(specInput, logFile, nod)
    
# Read spectrum for each detector bank if iSpec = 0
if iSpec == 0:
    # spectra is an array of arrays containing the spectra in the
    # Spectrum_run1_run2.dat file.
    specInput = open( spectraFile, 'r' )
    
    for i in range(8):   # skip the first 8 lines
        lineString = specInput.readline()
        # print lineString
    
    # "spectra" is an array spectra[i][j] where i is the number
    # of the detector bank starting at zero, and j = 0 for
    # the array of times and j = 1 for the array of counts
    spectra = []
    
    lineString = specInput.readline()   # read "Bank 1" line
    
    for i in range( nod ):
        # set arrays to zero
        time = []
        counts = []
        
        print('Reading spectrum for ' + lineString[0:-1])
        while True:
            lineString = specInput.readline()
            lineList = lineString.split()
            if len(lineList) == 0: break
            if lineList[0] == 'Bank': break
            time.append( float( lineList[0] ) )
            counts.append( float( lineList[1] ) )
            
        spectra.append( [time, counts] )
    
    specInput.close()
    

# C-----------------------------------------------------------------------
# C  Calculate spectral correction at normToWavelength to normalize
# C  spectral correction factors later on.
spect1 = []     # spectrum value at normToWavelength for each detector
dist = []       # sample-to-detector distance
xtof = []       # = (L1+dist)/hom; TOF = wl * xtof

wavelength = normToWavelength
one = 1.0       # denominator in spectrum to calculate spect1

for id in range(nod):

    if iSpec == 1:  # The spectrum is calculated from coefficients
        
        spect = spectrumCalc(wavelength, calibParam, pj, id)
        spect1.append(spect)
        
    else:           # iSpec = 2           
        
        dist.append(calibParam[9][id])
        xtof.append((L1 + dist[id]) / hom)
        
        # spectra[id][0] are the times-of-flight
        # spectra[id][1] are the counts
        spectx = spectrum2( wavelength, xtof[id], \
            one, spectra[id][0], spectra[id][1] )
        spect = spectx[0]            # the spectral normalization parameter
        relSigSpect = spectx[1]      # the relative sigma of spect
        if spect == 0.0:
            print('*** Wavelength for normalizing to spectrum is out of range.')
        spect1.append(spect)
                          
# C-----------------------------------------------------------------------

# C
# C  SET THE CURRENT HISTOGRAM NUMBER TO 0 AND INITIALIZE THE MONITOR COUN
# C
curhst = 0
idet = 0
hstnum = int(starting_batch_number) - 1    # Set the starting batch number
cmon = 9.89E+5                             # Scale proton charge to 1 MW-hr
ncntr = 0                                  # Number of processed reflections

nrun = 0
dn = 0
chi = 0.0
phi = 0.0
omega = 0.0
moncnt = 1000000.
eof = 999
hkllists =[]   # List of reflections, XP Wang, May 3013
#Get peak index [start = 0] 
peak_index = -1
fsqmax = 0.0
# C
# C   SET UP LOOP TO PROCESS THE REFLECTION DATA
# C
while True:

    peak = readrefl_SNS( integFile, eof, nrun, dn, chi, phi, omega,\
        moncnt)
    peak_index = peak_index + 1
    eof = peak[22]
    if eof == 0: break
    
    nrun = peak[0]
    dn = peak[1]
    chi = float( peak[2] )
    phi = float( peak[3] )
    omega = float( peak[4] )
    moncnt = peak[5]
    seqnum = peak[6]
    h = peak[7]
    k = peak[8]
    l = peak[9]
    col = peak[10]
    row = peak[11]
    chan = peak[12]
    L2 = peak[13]
    twoth = peak[14]  # radians
    az = peak[15]  # azimuthal angle in radians
    wl = peak[16]
    dsp = peak[17]
    ipkobs = peak[18]
    inti = peak[19]
    sigi = abs(peak[20])
    reflag = peak[21]
    #
    if maxOrder > 0:
        m = peak[23]
        n = peak[24]
        p = peak[25]
    else:
        m = 0
        n = 0
        p = 0

    if dn !=99:
        if (nrun != curhst or dn != idet):
            if nrun != curhst:
                curhst = nrun
                if iIQ != 2: hstnum = hstnum + 1
                
                # Rotate UB matrix if polyhedral absorption correction
                # if abs_correc_type == 'polyhedral':
                if True:
                    SNS_or_IPNS = 'SNS'
                    # Using UB_IPNS with SNS rotation angles
                    newmat = rotate_matrix(UB_IPNS, omega, chi, phi, SNS_or_IPNS)
                    
                    # Calculate direction cosines for reversed incident beam vector.
                    # IPNS coordinates.
                    R_reverse_incident = [ -L2, 0.0, 0.0 ]
                    dir_cos_1 = [ 0, 0, 0 ]
                    # Begin loop through a-star, b-star and c-star
                    for i in range(3):
                        hkl = [ 0, 0, 0 ]
                        hkl[i] = 1
                        q_abc_star = huq( hkl[0], hkl[1], hkl[2], newmat )
                        length_q_abc_star = math.sqrt( np.dot( q_abc_star, q_abc_star) )                    
                        dir_cos_1[i] = ( np.dot( R_reverse_incident, q_abc_star ) 
                            / ( L2 * length_q_abc_star ) )
                    
            idet = dn  #IDET and DN is the arbitrary detector number.
                       #ID is a sequential number in the order they are listed.
         
            for id in range(nod):
                detNum = calibParam[3][id]
                if detNum == dn: break
                
            if iIQ == 2: hstnum = hstnum + 1
            
            mnsum = moncnt
            
            if mnsum == 0:
                cmonx = 1.0
            else:
                cmonx = cmon / mnsum
                if cmonx == 0: cmonx = 1.0
            
            logFile.write('\n\nHISTOGRAM NUMBER %5d' % nrun)      
            logFile.write('\nDETECTOR BANK NUMBER %2d     DETECTOR SEQUENTIAL NUMBER %2d'\
                % (dn, id))
            logFile.write('\nANGLES ARE CHI =%7.2f   PHI =%7.2f   OMEGA=%7.2f\n'\
                % ( chi, phi, omega ))                        
            logFile.write('TOTAL MONITOR COUNTS ELAPSED%10d   CMONX =%10.4f\n'\
                % ( mnsum, cmonx ))
            logFile.write('* DATA SCALED TO 100 MILLION MONITOR COUNTS *\n')
            logFile.write('CORREC = SCALEFACTOR * CMONX * SINSQT /' + \
                '( SPECT * (DET EFF) * WL4 * ABTRANS )\n')
            logFile.write('\n    H   K   L   M   N   P       FSQ     SIG     WL      INTI' + \
                '     SIG   SPECT  SINSQT       l1       l2       m1       m2       n1       n2\n')
        # end of set-up for new run or detector
       
        # Omit zero intensity peaks from integrate file XP Wang 03/21/2011
        # Changed to >=0 and absolute value  XP Wang 02/24/2011

        # Omit peaks not indexed XP Wang, March, 2013
        if (h==0 and k==0 and l==0):
            logFile.write(' %4d *** Peak not indexed for run %4d det %4d   \n' \
                % (seqnum,nrun,dn))        
            continue  
        if inti == 0.0 :
            logFile.write(' %4d%4d%4d%4d%4d%4d *** intI = 0.0 \n' \
                % (h, k, l, m, n, p))
            continue  

        if isnan(sigi) == True:
            logFile.write(' %4d%4d%4d%4d%4d%4d *** sigi = nan \n' \
                % (h, k, l, m, n, p))
            continue  

        if minIsigI >= 0 and inti < abs(minIsigI * sigi):
            logFile.write(' %4d%4d%4d%4d%4d%4d *** inti < (minIsigI * sigi) \n' \
                % (h, k, l, m, n, p))
            continue
            
        if inti < intiMin:
            logFile.write(' %4d%4d%4d%4d%4d%4d *** inti < intiMin \n' \
                % (h, k, l, m, n, p))
            continue

        # Set-up limits for neutron wavelentgh XP Wang 02/24/2011
        if wl < wlMin:
            logFile.write(' %4d%4d%4d%4d%4d%4d *** wl < wlMin \n' \
                % (h, k, l, m, n, p))
            continue

        if wl > wlMax:
            logFile.write(' %4d%4d%4d%4d%4d%4d *** wl > wlMax \n' \
                % (h, k, l, m, n, p))
            continue

        nRows = calibParam[4][id]
        nCols = calibParam[5][id]
        
        if col < numBorderCh:
            logFile.write(' %4d%4d%4d%4d%4d%4d *** col < numBorderCh \n' \
                % (h, k, l, m, n, p))
            continue
            
        if col > (nCols - numBorderCh):
            logFile.write(' %4d%4d%4d%4d%4d%4d *** col > (nCols - numBorderCh)\n' \
                % (h, k, l, m, n, p))
            continue
            
        if row < numBorderCh:
            logFile.write(' %4d%4d%4d%4d%4d%4d *** row < numBorderCh \n' \
                % (h, k, l, m, n, p))
            continue
            
        if row > (nRows - numBorderCh):
            logFile.write(' %4d%4d%4d%4d%4d%4d *** row > (nRows - numBorderCh)\n' \
                % (h, k, l, m, n, p))
            continue
                            
        if dsp < dMin:
            logFile.write(' %4d%4d%4d%4d%4d%4d *** dsp < dMin \n' \
                % (h, k, l, m, n, p))
            continue
        
        ncntr = ncntr + 1
        
        if iSpec == 1:
            spect = spectrumCalc(wl, calibParam, pj, id)
            spect = spect / spect1[id]
        
        if iSpec == 0:
            spectx = spectrum2( wl, xtof[id], \
              spect1[id], spectra[id][0], spectra[id][1] )
            spect = spectx[0]
            relSigSpect = spectx[1]
        if spect == 0.0:
            logFile.write(' %4d%4d%4d%4d%4d%4d *** spect == 0.0 \n' \
                % (h, k, l, m, n, p))
            continue
        
        # correct for the slant path throught the scintillator glass
        mu = (9.614 * wl) + 0.266    # mu for GS20 glass
        depth = calibParam[8][id]
        eff_center = 1.0 - exp(-mu * depth)  # efficiency at center of detector
        cosA = dist[id] / L2
        pathlength = depth / cosA
        eff_R = 1.0 - exp(-mu * pathlength)   # efficiency at point R
        sp_ratio = eff_center / eff_R  # slant path efficiency ratio
        
        sinsqt = ( wl / (2.0*dsp) )**2
        wl4 = wl**4    
        correc = sinsqt * sp_ratio / (wl4 * spect ) * detScale[detNum]
            
        # absorption correction
        # trans[0] is the transmission
        # trans[1] is tbar
        if abs_correc_type == 'spherical':
            trans = absor_sphere(smu, amu, radius, twoth, wl)
            transmission = trans[0]
            if trans[0] < transmin: transmin = trans[0]
            if trans[0] > transmax: transmax = trans[0]
            
            correc = correc / trans[0]
        
        fsq = inti * correc * cmonx  
        sigfsq = sigi * correc * cmonx

        # Include instrument background constant in sigma        
        sigfsq = sqrt( sigfsq**2 + (relSigSpect*fsq)**2 + 1.94)   
    
        fsq = fsq * scaleFactor
        sigfsq = sigfsq * scaleFactor

        # Calculates direction cosines for scattered beam vector 
        R_IPNS, R_SNS = Rvec(twoth, az, L2)
        # Begin loop through a-star, b-star and c-star
        dir_cos_2 = [ 0, 0, 0 ]   # dir_cos_2 is the scattered beam with a*, b*, c*
        for i in range(3):
            abc = [ 0, 0, 0 ]
            abc[i] = 1
            q_abc_star = huq( abc[0], abc[1], abc[2], newmat )
            len_q_abc_star = math.sqrt( np.dot( q_abc_star, q_abc_star) )
            dir_cos_2[i] = np.dot( R_IPNS, q_abc_star ) / ( L2 * len_q_abc_star )
        
        # Write out to hkl file for spherical absorption correction.    
        if abs_correc_type == 'spherical':
            # tbar is the Coppen's tbar
            tbar = trans[1]
            
            # output reflection to log file and to hkl file
            if maxOrder > 0:
                logFile.write( ' '
                    + 6*'%4d' % (h,k,l,m,n,p)
                    + '%10.2f%8.2f%7.3f%10.2f%8.2f' % (fsq, sigfsq, wl, inti, sigi)
                    + 4*'%8.4f' % (spect, sinsqt, trans[0], tbar)
                    + '\n' )
            
                hklFile.write( 6*'%4d' % (h,k,l,m,n,p)
                    + f'{round(fsq,5):>8g}{round(sigfsq,5):>8g}' 
                    + '%4d' % hstnum 
                    + 2*'%8.5f' % (wl, tbar)
                    + 6*'%9.5f' % (dir_cos_1[0], dir_cos_2[0], dir_cos_1[1], dir_cos_2[1], dir_cos_1[2], dir_cos_2[2])
                    + '%6d' % (curhst)
                    + '%7d' % (seqnum)
                    + '%7.4f%4d%9.5f%8.4f%7.2f%7.2f\n' % (transmission, dn, twoth, dsp, col, row) )            
            else:
                logFile.write( ' '
                    + 3*'%4d' % (h,k,l)
                    + '%10.2f%8.2f%7.3f%10.2f%8.2f' % (fsq, sigfsq, wl, inti, sigi)
                    + 4*'%8.4f' % (spect, sinsqt, trans[0], tbar)
                    + '\n' )
            
                hklFile.write( 3*'%4d' % (h,k,l)
                    + f'{round(fsq,5):>8g}{round(sigfsq,5):>8g}' 
                    + '%4d' % hstnum 
                    + 2*'%8.5f' % (wl, tbar)
                    + 6*'%9.5f' % (dir_cos_1[0], dir_cos_2[0], dir_cos_1[1], dir_cos_2[1], dir_cos_1[2], dir_cos_2[2])
                    + '%6d' % (curhst)
                    + '%7d' % (seqnum)
                    + '%7.4f%4d%9.5f%8.4f%7.2f%7.2f\n' % (transmission, dn, twoth, dsp, col, row) )            
        # Write out to hkl file for polyhedral absorption correction.    
        if abs_correc_type == 'polyhedral':
        
            # output reflection to log file and to hkl file
            tbar = 0.0
            transmission = 1.0

            if maxOrder > 0:
                logFile.write(' '
                    + f'{h:>4d}{k:>4d}{l:4d}{m:4d}{n:4d}{p:4d}'
                    + '%10.2f%8.2f%7.3f%10.2f%8.2f' % (fsq, sigfsq, wl, inti, sigi)
                    + f'{round(spect,5):>8g}{round(sinsqt,5):>8g}'
                    + 6*'%9.5f' % (dir_cos_1[0], dir_cos_2[0], dir_cos_1[1], dir_cos_2[1], dir_cos_1[2], dir_cos_2[2])
                    + '\n' )
            
                hklFile.write( f'{h:>4d}{k:>4d}{l:4d}{m:4d}{n:4d}{p:4d}'
                    + f'{round(fsq,5):>8g}{round(sigfsq,5):>8g}'
                    + '%4d' % hstnum 
                    + 2*'%8.5f' % (wl, tbar)
                    + 6*'%9.5f' % (dir_cos_1[0], dir_cos_2[0], dir_cos_1[1], dir_cos_2[1], dir_cos_1[2], dir_cos_2[2])
                    + '%6d' % (curhst)
                    + '%7d' % (seqnum)
                    + '%7.4f%4d%9.5f%8.4f%7.2f%7.2f\n' % (transmission, dn, twoth, dsp, col, row) )
            else:
                logFile.write(' '
                    + 3*'%4d' % (h, k, l)
                    + '%10.2f%8.2f%7.3f%10.2f%8.2f' % (fsq, sigfsq, wl, inti, sigi)
                    + f'{round(spect,5):>8g}{round(sinsqt,5):>8g}'
                    + 6*'%9.5f' % (dir_cos_1[0], dir_cos_2[0], dir_cos_1[1], dir_cos_2[1], dir_cos_1[2], dir_cos_2[2])
                    + '\n' )

                hklFile.write( 3*'%4d' % (h, k, l)
                    + f'{round(fsq,5):>8g}{round(sigfsq,5):>8g}'
                    + '%4d' % hstnum 
                    + 2*'%8.5f' % (wl, tbar)
                    + 6*'%9.5f' % (dir_cos_1[0], dir_cos_2[0], dir_cos_1[1], dir_cos_2[1], dir_cos_1[2], dir_cos_2[2])
                    + '%6d' % (curhst)
                    + '%7d' % (seqnum)
                    + '%7.4f%4d%9.5f%8.4f%7.2f%7.2f\n' % (transmission, dn, twoth, dsp, col, row) )
        #Add to hkl list XP Wang May 2013
        #Mantid uses Inelastic convention, k_i - k_f, as the default setting
        if config['Q.convention']=='Inelastic':
            #Change to diffraction convention, k_f - k_i
            q=peaks_ws.cell(peak_index,16)*(-1.)
        else:
            q=peaks_ws.cell(peak_index,16)
        #
        Q_sample =np.array([q.getX(),q.getY(),q.getZ()])
        if maxOrder > 0:
            hkllists.append([ h, k, l, fsq, sigfsq, hstnum, wl, tbar,
                dir_cos_1[0], dir_cos_2[0], dir_cos_1[1], dir_cos_2[1], 
                dir_cos_1[2], dir_cos_2[2], curhst, seqnum, transmission,
                dn, twoth, dsp, col, row, Q_sample[0],Q_sample[1],Q_sample[2], m, n, p])
        else:
            hkllists.append([ h, k, l, fsq, sigfsq, hstnum, wl, tbar,
                dir_cos_1[0], dir_cos_2[0], dir_cos_1[1], dir_cos_2[1], 
                dir_cos_1[2], dir_cos_2[2], curhst, seqnum, transmission,
                dn, twoth, dsp, col, row, Q_sample[0],Q_sample[1],Q_sample[2]])
        if fsq > fsqmax: fsqmax = fsq        
            
if abs_correc_type == 'spherical':
    print('\nMinimum and maximum transmission = %6.4f, %6.4f' % (transmin, transmax))

logFile.write('\n\n***** Minimum and maximum transmission = %6.4f, %6.4f' \
    % (transmin, transmax))
    
if fsqmax > (100000.00 -1):
    print()
    print('################################################################')
    print('     Maximum FOSQ is', fsqmax)
    print('     Re-run anvred with a scale factor of', 0.1*scaleFactor)
    print('################################################################')

# last record all zeros for shelx
iz = 0
rz = 0.0
if maxOrder > 0:
    hklFile.write( f'{iz:>4d}{iz:>4d}{iz:>4d}{iz:>4d}{iz:>4d}{iz:>4d}'
        + f'{rz:8.2f}{rz:8.2f}'
        + f'{iz:>4d}'
        + f'{rz:8.5f}{rz:8.5f}'
        + f'{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}'
        + f'{iz:>6d}'
        + f'{iz:>7d}'
        + '%7.4f%4d%9.5f%8.4f%7.2f%7.2f\n' % (rz, iz, rz, rz, rz, rz) )
else:
    hklFile.write( 3*'%4d' % (iz, iz, iz)
        + f'{rz:8.2f}{rz:8.2f}'
        + f'{iz:>4d}'
        + f'{rz:8.5f}{rz:8.5f}'
        + f'{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}'
        + f'{iz:>6d}'
        + f'{iz:>7d}'
        + '%7.4f%4d%9.5f%8.4f%7.2f%7.2f\n' % (rz, iz, rz, rz, rz, rz) )
# C-----------------------------------------------------------------------

hklFile.close()

#
a =np.array(hkllists)
#Absorption coefficient at 1.0 Angstrom in mm^-1
mu_1A = (smu + amu / 1.7982) *1.0 / 10.0
if abs_correc_type=='polyhedral':
   abs_type ='gaussian'
   abs_ref = str('  P. Coppens, L. Leiserowitz, D Rabinovich, Acta Cryst. \n'
    + '  (1965), 18, 1035-1038. \n'
    + '  P. R. Mallinson and K. W. Muir, J. Appl. Cryst. \n'
    + '  (1985). 18, 51-53.')
else:
  abs_type = 'sphere'
  abs_ref = str('C. W. Dwiggins, Jr., Acta Cryst. A31, 395 (1975).') 

# Output cif for TOPAZ experiment
cif = 'topaz_anvred3.cif'
cif_output = open(cif, 'w')
cif_output.write('# CIF for TOPAZ neutron time of flight Laue Experiment\n'
    + '# Data directory \n'
    + '# %s\n' % (output_directory)
    + '# Raw data file name  \n'
    + '# %s\n' % (os.path.splitext(os.path.basename(anvred_integrate_fname))[0] +'.integrate')
    + '\n'
    + 'data_topaz \n'
    + ' \n'
    + '_diffrn_ambient_temperature          ?\n'
    + '_diffrn_radiation_wavelength         "%4.2f-%5.2f"\n'%(np.min(a[:,6]),np.max(a[:,6]))
    + '_diffrn_radiation_wavelength_details  "Based on neutron time of flight"\n'
    + '_diffrn_radiation_type               neutrons\n'
    + '_diffrn_source                       "The ORNL Spallation Neutron Source"\n'
    + '_diffrn_measurement_device_type      TOPAZ\n'
    + '_diffrn_measurement_method           "time-of-flight Laue"\n'
    + '_diffrn_detector_area_resol_mean     ?\n'
    + '\n'
    + '_diffrn_reflns_theta_min         %7.3f\n'%(np.min(a[:,18])*180.0/pi/2.0)
    + '_diffrn_reflns_theta_max         %7.3f\n'%(np.max(a[:,18])*180.0/pi/2.0)
    + '_diffrn_reflns_min_d_Angs        %7.3f'%(np.min(a[:,19]))
    + ' \n'
    + '_exptl_absorpt_correction_type    %s\n' % (abs_type)
    + '_exptl_absorpt_coefficient_mu     %7.4f\n' % (mu_1A)
    + '_exptl_absorpt_correction_T_min   %7.4f\n' % (np.min(a[:,16]))
    + '_exptl_absorpt_correction_T_max   %7.4f\n' % (np.max(a[:,16]))
    + '_exptl_absorpt_process_details    \n'
    + '; \n'
    + '%s\n'%(abs_ref)
    + '; \n'
    + '\n'
    + '_exptl_absorpt_special_details    \n'
    + '; \n'
    + 'Neutron linear absorption coefficient is wavelength dependent.\n'
    + 'For each peak, the absorption coefficient mu is calculated as:\n'
    + '          mu = %6.4f + %6.4f/1.7982 * wavelength  [cm^-1^]\n' % (smu,amu)
    + 'where \n'
    + ' %6.4f cm^-1 is the absorption coefficient for total scattering;\n' %(smu) 
    + ' %6.4f cm^-1 is the absorption coefficient for true absorption at 1.7982 A.\n' % (amu)        
    + ' \n'
    + 'The value of %6.4f (mm^-1^) shown in _exptl_absorpt_coefficient_mu \n' % (mu_1A)
    + 'is the sample absorption coefficient for neutron wavelength at 1.0 \%A.\n'
    + ';\n'
    + '\n'
    + '_computing_data_collection        "SNS PyDas"\n'
    + '_computing_cell_refinement         Mantidplot\n'
    + '_computing_data_reduction          Mantidplot\n')

cif_output.close()
#
#Print theta angle and d-spacing range  
print('min d-spacing : %8.4f'%(np.min(a[:,19])))
print('max d-spacing : %8.4f'%(np.max(a[:,19])))
#
# Set scale ID equal to detector number.
# This code is from scale_by_detnum.py.
# Save reflections by DetNum  XP Wang May, 2013

if iIQ == 1 or iIQ == 3:
    hklFileName1 = hklFileName + '_dn'
    hkl_output = open(hklFileName1, 'w')
    
    #Sort and save the result per module, XP Wang, March 2013
    hkllists.sort(key=itemgetter(17,0,1,2))  # sort by detector number
    nDet = 0
    for iDet, iGroup in groupby(hkllists, itemgetter(17)):
        nDet = nDet + 1
        for iHKL in iGroup:
            iHKL[5] = nDet                
            # output reflection sorted by detector number to hkl file
            if maxOrder > 0:
                hkl_output.write(f'{iHKL[0]:>4d}{iHKL[1]:>4d}{iHKL[2]:4d}{iHKL[25]:4d}{iHKL[26]:4d}{iHKL[27]:4d}'
                    + f'{round(iHKL[3],5):>8g}{round(iHKL[4],5):>8g}'
                    + '%4d' % iHKL[5] 
                    + f'{iHKL[6]:8.5f}{iHKL[7]:8.5f}'
                    + f'{iHKL[8]:9.5f}{iHKL[9]:9.5f}{iHKL[10]:9.5f}{iHKL[11]:9.5f}{iHKL[12]:9.5f}{iHKL[13]:9.5f}'
                    + f'{iHKL[14]:>6d}'
                    + f'{iHKL[15]:>7d}'
                    + f'{iHKL[16]:7.4f}{iHKL[17]:>4d}{iHKL[18]:9.5f}{iHKL[19]:8.4f}'
                    + f'{iHKL[20]:7.2f}{iHKL[21]:7.2f}'
                    + '\n' )
            else:
                hkl_output.write( 3*'%4d' % (iHKL[0],iHKL[1],iHKL[2])
                    + f'{round(iHKL[3],5):>8g}{round(iHKL[4],5):>8g}'
                    + '%4d' % iHKL[5] 
                    + f'{iHKL[6]:8.5f}{iHKL[7]:8.5f}'
                    + f'{iHKL[8]:9.5f}{iHKL[9]:9.5f}{iHKL[10]:9.5f}{iHKL[11]:9.5f}{iHKL[12]:9.5f}{iHKL[13]:9.5f}'
                    + f'{iHKL[14]:>6d}'
                    + f'{iHKL[15]:>7d}'
                    + f'{iHKL[16]:7.4f}{iHKL[17]:>4d}{iHKL[18]:9.5f}{iHKL[19]:8.4f}'
                    + f'{iHKL[20]:7.2f}{iHKL[21]:7.2f}'
                    + '\n' )
                                
    # last record all zeros for shelx
    if maxOrder > 0:
        hkl_output.write( f'{iz:>4d}{iz:>4d}{iz:>4d}{iz:>4d}{iz:>4d}{iz:>4d}'
            + f'{rz: 8.2f}{rz: 8.2f}'
            + f'{iz:>4d}'
            + f'{rz:8.5f}{rz:8.5f}'
            + f'{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}'
            + f'{iz:>6d}'
            + f'{iz:>7d}'
            + f'{rz:7.4f}{iz:4d}{rz:9.5f}{rz:8.4f}'
            + f'{rz:7.2f}{rz:7.2f}'
            + '\n' )
    else:
        hkl_output.write( 3*'%4d' % (iz, iz, iz)
            + f'{rz: 8.2f}{rz: 8.2f}'
            + f'{iz:>4d}'
            + f'{rz:8.5f}{rz:8.5f}'
            + f'{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}'
            + f'{iz:>6d}'
            + f'{iz:>7d}'
            + f'{rz:7.4f}{iz:4d}{rz:9.5f}{rz:8.4f}'
            + f'{rz:7.2f}{rz:7.2f}'
            + '\n' )

hkl_output.close()

hkllists.sort(key=itemgetter(15))  # sort peaks by sequence number
peak_num_anvred=np.array(hkllists)[:,15].astype(np.int)

#
jana_fname = os.path.splitext(anvred_integrate_fname)[0] + '.topi'
print('Writing JANA topi file ...\n' + jana_fname)
topi_f=open(jana_fname,'w')
#
#Q_sample in IPNS coordinates system, x-along neutron beam; y-perpendicular to neutron beam; z-up
to_IPNS=np.array([[0,0,1],[1,0,0],[0,1,0]])
#
hkllists.sort(key=itemgetter(15))  # sort peaks by sequence number

#Output JANA topi file with batch number sorted by runs 
nBatch = 0
for iBatch, iGroup in groupby(hkllists, itemgetter(14)):                
    nBatch = nBatch + 1
    for runHKL in iGroup:
        runHKL[5] = nBatch               
        #Use IPNS coordinates
        Q_sample =np.array([runHKL[22],runHKL[23],runHKL[24]])
        Q_sample=np.dot(to_IPNS,Q_sample)
        #
        #Use normalized intensity and sigma after anvred correction
        #Use normalized intensity and sigma after anvred correction
        intensity=runHKL[3]
        error=runHKL[4]
        wl=runHKL[6]
        hstnum=runHKL[5]
        topi_f.write(2*'%9.2f' % (intensity,error)
                   + 3*'%13.8f' % (Q_sample[0],Q_sample[1],Q_sample[2])
                   + '%8.5f' % (wl) + '%4d' % (hstnum) + '\n')
topi_f.close()

# Peak statistics

#Remove peaks from anvred
#Replace peak intensity and sigma with normalized fsqr and sigma_fsqr from anvred
for i in range(peaks_ws.getNumberPeaks()):
    pk=peaks_ws.getPeak(i)
    pki=pk.getPeakNumber()
    anvred_index=np.where(np.array(peak_num_anvred)==pki)[0].astype(int)
    if len(anvred_index) == 0:
        pk.setIntensity(0)
        pk.setSigmaIntensity(0)
        continue
    else:
        anvred_index=int(anvred_index[0])
        pk.setIntensity(hkllists[anvred_index][3])
        pk.setSigmaIntensity(hkllists[anvred_index][4])
        #print(hkl, pki, index,pk.getIntensity(), pk.getSigmaIntensity())

peaks_ws = FilterPeaks(InputWorkspace = peaks_ws,  FilterVariable = 'Intensity',  FilterValue = 0,  Operator = '!=')

#
print('\nNumber of Peaks from anvred correction : {0}'.format(peaks_ws.getNumberPeaks()))

logFile.write('\n\nNumber of Peaks from anvred correction : {0}'.format(peaks_ws.getNumberPeaks()))

point_group = PointGroupFactory.createPointGroup(pg_symbol)

if centering == "P":
     lattice_centering='Primitive'
if centering == "C":
     lattice_centering='C-face centred'
if centering == "A":
     lattice_centering='A-face centred'
if centering == "B":
     lattice_centering='B-face centred'
if centering == "I":
     lattice_centering='Body centred'
if centering == "F":
     lattice_centering='All-face centred'
if centering == "R":
     lattice_centering='Rhombohedrally centred, obverse'

print('\nLattice Type: {}'.format(lattice_centering))

#
if str(point_group.getCrystalSystem()) =='Monoclinic' or str(point_group.getCrystalSystem()) =='Trigonal':
    pg_symbol =str(point_group.getHMSymbol())
else:
    pg_symbol =str(str(point_group.getHMSymbol()) +' (' +str(point_group.getCrystalSystem())+')')
print('Point Group symmetry: {}'.format(pg_symbol))
print('Z score: ', z_score)


# Run the SortHKL algorithm
sorted, statistics_table, eq_peaks_list = StatisticsOfPeaksWorkspace(InputWorkspace=peaks_ws, 
                                              PointGroup=pg_symbol, 
                                              LatticeCentering= lattice_centering, 
                                              OutputWorkspace='OutputPeaks', 
                                              SortBy='Overall', 
                                              SigmaCritical=z_score, 
                                              EquivalentIntensities="Median",
                                              WeightedZScore=True,
                                              EquivalentsWorkspace='peaks_eq')

statistics = statistics_table.row(0)

print('\nCrystal symmetry')
print('               Point Group: {}'.format(point_group.getHMSymbol()))
print('          Lattice Sysstem: {}'.format(point_group.getCrystalSystem()))
print('        Lattice Centering: {}'.format(lattice_centering))
print('\nPeak Statistics')
print('          Number of Peaks: {0}'.format(sorted.getNumberPeaks()))
print('             Multiplicity: {0}'.format(round(statistics['Multiplicity'],2)))
print('        Data Completeness: {0}%'.format(round(statistics['Data Completeness'],2)))
print('           Resolution Min: {0}'.format(round(statistics['Resolution Min'],2)))
print('           Resolution Max: {0}'.format(round(statistics['Resolution Max'],2)))
print('No. of Unique Reflections: {0}'.format(statistics['No. of Unique Reflections']))
print('         Mean ((I)/sd(I)): {0}'.format(round(statistics['Mean ((I)/sd(I))'],2)))
print('                   Rmerge: {0}%'.format(round(statistics['Rmerge'],2)))
print('                     Rpim: {0}%'.format(round(statistics['Rpim'],2)))
print('\n')

logFile.write('\nLattice Type: {}\n'.format(lattice_centering))
logFile.write('Point Group symmetry: {}'.format(pg_symbol))
logFile.write('Z score: {}'.format(z_score))
logFile.write('\nCrystal symmetry')
logFile.write('\n               Point Group: {}'.format(point_group.getHMSymbol()))
logFile.write('\n          Lattice Sysstem: {}'.format(point_group.getCrystalSystem()))
logFile.write('\n        Lattice Centering: {}'.format(lattice_centering))
logFile.write('\nPeak Statistics')
logFile.write('\n          Number of Peaks: {0}'.format(sorted.getNumberPeaks()))
logFile.write('\n             Multiplicity: {0}'.format(round(statistics['Multiplicity'],2)))
logFile.write('\n        Data Completeness: {0}%'.format(round(statistics['Data Completeness'],2)))
logFile.write('\n           Resolution Min: {0}'.format(round(statistics['Resolution Min'],2)))
logFile.write('\n           Resolution Max: {0}'.format(round(statistics['Resolution Max'],2)))
logFile.write('\nNo. of Unique Reflections: {0}'.format(statistics['No. of Unique Reflections']))
logFile.write('\n         Mean ((I)/sd(I)): {0}'.format(round(statistics['Mean ((I)/sd(I))'],2)))
logFile.write('\n                   Rmerge: {0}%'.format(round(statistics['Rmerge'],2)))
logFile.write('\n                     Rpim: {0}%'.format(round(statistics['Rpim'],2)))
logFile.write('\n')

rn=peaks_ws.column(0)
pix=peaks_ws.column(1)
intensity=peaks_ws.column(9)
sigma=peaks_ws.column(10)
peak_seqnum = []

for i in range(peaks_ws.getNumberPeaks()):
    pk =peaks_ws.getPeak(i)
    peak_seqnum.extend([pk.getPeakNumber()])

hkl_out=[]    
peak_num_anvred=np.array(hkllists)[:,15].astype(np.int)

print('Number of peaks after outlier removal: {0}'.format(mtd['OutputPeaks'].getNumberPeaks()))

logFile.write('\nNumber of peaks after outlier removal: {0}'.format(mtd['OutputPeaks'].getNumberPeaks()))
print('\nSaving result ...')

for i in range(mtd['OutputPeaks'].getNumberPeaks()):
    pk=mtd['OutputPeaks'].getPeak(i)
    hkl=pk.getHKL()
    #rni=pk.getRunNumber()
    #pixi=pk.getDetectorID()
    pki=pk.getPeakNumber()
    #index=int(np.where(np.logical_and(np.array(rn)==rni,np.array(detpixel)==pixi))[0])
    #get index number from peaks workspace
    index=int(np.where(np.array(peak_seqnum)==pki)[0])
    #match index number from peak list after anvred correction
    anvred_index=np.where(np.array(peak_num_anvred)==pki)[0].astype(int)
    #print(hkl,  pki, index,pk.getIntensity(), pk.getSigmaIntensity())
    pk.setIntensity(intensity[index])
    pk.setSigmaIntensity(sigma[index])
    #print(hkl, pki, index,pk.getIntensity(), pk.getSigmaIntensity())
    if len(anvred_index) == 0:     
       continue
    else:
        #print hkl_index
        #print hkllists[hkl_index][0:6],hkllists[hkl_index][15]
        anvred_index=int(anvred_index[0])
        hkl_out.extend([hkllists[anvred_index]])

#
# Save reflections after removing outliers by run numbers XP Wang Feb, 2018
print(hklFileName)
if iIQ == 1 or iIQ == 3:
    hklFileName2 = os.path.splitext(hklFileName)[0] + '_symm.hkl'
    hkl_output2 = open(hklFileName2, 'w')
    
    #Sort and save the result per module, XP Wang, March 2013
    hkl_out.sort(key=itemgetter(14,15))  # sort by run and seqnum number
    nBatch = int(starting_batch_number) - 1
    for iBatch, iGroup in groupby(hkl_out, itemgetter(14)):
        nBatch = nBatch + 1
        for iHKL in iGroup:
            iHKL[5] = nBatch               
            # output reflection sorted by run number to hkl file
            if maxOrder > 0:
                hkl_output2.write(f'{iHKL[0]:>4d}{iHKL[1]:>4d}{iHKL[2]:4d}{iHKL[25]:4d}{iHKL[26]:4d}{iHKL[27]:4d}'
                    + f'{round(iHKL[3],5):>8g}{round(iHKL[4],5):>8g}'
                    + '%4d' % iHKL[5] 
                    + f'{iHKL[6]:8.5f}{iHKL[7]:8.5f}'
                    + f'{iHKL[8]:9.5f}{iHKL[9]:9.5f}{iHKL[10]:9.5f}{iHKL[11]:9.5f}{iHKL[12]:9.5f}{iHKL[13]:9.5f}'
                    + f'{iHKL[14]:>6d}'
                    + f'{iHKL[15]:>7d}'
                    + f'{iHKL[16]:7.4f}{iHKL[17]:>4d}{iHKL[18]:9.5f}{iHKL[19]:8.4f}'
                    + f'{iHKL[20]:7.2f}{iHKL[21]:7.2f}'
                    + '\n' )
            else:
                hkl_output2.write( 3*'%4d' % (iHKL[0],iHKL[1],iHKL[2])
                    + f'{round(iHKL[3],5):>8g}{round(iHKL[4],5):>8g}'
                    + '%4d' % iHKL[5] 
                    + f'{iHKL[6]:8.5f}{iHKL[7]:8.5f}'
                    + f'{iHKL[8]:9.5f}{iHKL[9]:9.5f}{iHKL[10]:9.5f}{iHKL[11]:9.5f}{iHKL[12]:9.5f}{iHKL[13]:9.5f}'
                    + f'{iHKL[14]:>6d}'
                    + f'{iHKL[15]:>7d}'
                    + f'{iHKL[16]:7.4f}{iHKL[17]:>4d}{iHKL[18]:9.5f}{iHKL[19]:8.4f}'
                    + f'{iHKL[20]:7.2f}{iHKL[21]:7.2f}'
                    + '\n' )
                                
    # last record all zeros for shelx
    if maxOrder > 0:
        hkl_output2.write( f'{iz:>4d}{iz:>4d}{iz:>4d}{iz:>4d}{iz:>4d}{iz:>4d}'
            + f'{rz:8.2f}{rz:8.2f}'
            + f'{iz:>4d}'
            + f'{rz:8.5f}{rz:8.5f}'
            + f'{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}'
            + f'{iz:>6d}'
            + f'{iz:>7d}'
            + f'{rz:7.4f}{iz:4d}{rz:9.5f}{rz:8.4f}'
            + f'{rz:7.2f}{rz:7.2f}'
            + '\n' )
    else:
        hkl_output2.write( 3*'%4d' % (iz, iz, iz)
            + f'{rz: 8.2f}{rz: 8.2f}'
            + f'{iz:>4d}'
            + f'{rz:8.5f}{rz:8.5f}'
            + f'{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}'
            + f'{iz:>6d}'
            + f'{iz:>7d}'
            + f'{rz:7.4f}{iz:4d}{rz:9.5f}{rz:8.4f}'
            + f'{rz:7.2f}{rz:7.2f}'
            + '\n' )

hkl_output2.close()

print(hklFileName2)
#
# Save reflections after removing outliers by detector numbers XP Wang Feb, 2018
if iIQ == 1 or iIQ == 3:
    hklFileName2 = os.path.splitext(hklFileName)[0] + '_symm.hkl_dn'
    hkl_output2 = open(hklFileName2, 'w')
    
    #Sort and save the result per module, XP Wang, March 2013
    hkl_out.sort(key=itemgetter(17,0,1,2,3))  # sort by detector and hkl
    nBatch = 0
    for iBatch, iGroup in groupby(hkl_out, itemgetter(17)):
        nBatch = nBatch + 1
        for iHKL in iGroup:
            iHKL[5] = nBatch               
            # output reflection sorted by detector number to hkl file
            if maxOrder > 0:
                hkl_output2.write(f'{iHKL[0]:>4d}{iHKL[1]:>4d}{iHKL[2]:4d}{iHKL[25]:4d}{iHKL[26]:4d}{iHKL[27]:4d}'
                    + f'{round(iHKL[3],5):>8g}{round(iHKL[4],5):>8g}'
                    + '%4d' % iHKL[5] 
                    + f'{iHKL[6]:8.5f}{iHKL[7]:8.5f}'
                    + f'{iHKL[8]:9.5f}{iHKL[9]:9.5f}{iHKL[10]:9.5f}{iHKL[11]:9.5f}{iHKL[12]:9.5f}{iHKL[13]:9.5f}'
                    + f'{iHKL[14]:>6d}'
                    + f'{iHKL[15]:>7d}'
                    + f'{iHKL[16]:7.4f}{iHKL[17]:>4d}{iHKL[18]:9.5f}{iHKL[19]:8.4f}'
                    + f'{iHKL[20]:7.2f}{iHKL[21]:7.2f}'
                    + '\n' )
            else:
                hkl_output2.write( 3*'%4d' % (iHKL[0],iHKL[1],iHKL[2])
                    + f'{round(iHKL[3],5):>8g}{round(iHKL[4],5):>8g}'
                    + '%4d' % iHKL[5] 
                    + f'{iHKL[6]:8.5f}{iHKL[7]:8.5f}'
                    + f'{iHKL[8]:9.5f}{iHKL[9]:9.5f}{iHKL[10]:9.5f}{iHKL[11]:9.5f}{iHKL[12]:9.5f}{iHKL[13]:9.5f}'
                    + f'{iHKL[14]:>6d}'
                    + f'{iHKL[15]:>7d}'
                    + f'{iHKL[16]:7.4f}{iHKL[17]:>4d}{iHKL[18]:9.5f}{iHKL[19]:8.4f}'
                    + f'{iHKL[20]:7.2f}{iHKL[21]:7.2f}'
                    + '\n' )
                                
    # last record all zeros for shelx
    if maxOrder > 0:
        hkl_output2.write( f'{iz:>4d}{iz:>4d}{iz:>4d}{iz:>4d}{iz:>4d}{iz:>4d}'
            + f'{rz: 8.2f}{rz: 8.2f}'
            + f'{iz:>4d}'
            + f'{rz:8.5f}{rz:8.5f}'
            + f'{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}'
            + f'{iz:>6d}'
            + f'{iz:>7d}'
            + f'{rz:7.4f}{iz:4d}{rz:9.5f}{rz:8.4f}'
            + f'{rz:7.2f}{rz:7.2f}'
            + '\n' )
    else:
        hkl_output2.write( 3*'%4d' % (iz, iz, iz)
            + f'{rz: 8.2f}{rz: 8.2f}'
            + f'{iz:>4d}'
            + f'{rz:8.5f}{rz:8.5f}'
            + f'{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}{rz:9.5f}'
            + f'{iz:>6d}'
            + f'{iz:>7d}'
            + f'{rz:7.4f}{iz:4d}{rz:9.5f}{rz:8.4f}'
            + f'{rz:7.2f}{rz:7.2f}'
            + '\n' )

hkl_output2.close()

print("\n**************************************************************************************")
print("****************************** All DONE **********************************************")
print("**************************************************************************************\n")

print('Config file used for data reduction : ' + config_file_name) 

logFile.write('\nConfig file used for data reduction : ' + config_file_name) 

logFile.close()

