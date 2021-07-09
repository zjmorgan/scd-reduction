# File: ReduceOneSCD_Run.py
#
# Version 2.0, modified to work with Mantid's new python interface.
#
# This script will reduce one SCD run.  The configuration file name and
# the run to be processed must be specified as the first two command line 
# parameters.  This script is intended to be run in parallel using the 
# ReduceSCD_Parallel.py script, after this script and configuration file has 
# been tested to work properly for one run. This script will load, find peaks,
# index and integrate either found or predicted peaks for the specified run.  
# Either sphere integration or the Mantid PeakIntegration algorithms are 
# currently supported, but it may be updated to support other integration 
# methods.  Users should make a directory to hold the output of this script, 
# and must specify that output directory in the configuration file that 
# provides the parameters to this script.
#
# NOTE: All of the parameters that the user must specify are listed with 
# instructive comments in the sample configuration file: ReduceSCD.config.
#

#
# _v1: December 3rd 2013. Mads Joergensen
# This version now includes the posibility to use the 1D cylindrical integration method
# and the posibility to load a UB matrix which will be used for integration of the individual
# runs and to index the combined file (Code from Xiapoing).
#
# _v2: December 5 2013. A. J. Schultz
# This version includes the Boolean parameter use_monitor_counts to allow
# the use of either monitor counts (True) or proton charge (False) for
# scaling.
#
# _v3: September 14,2017. Xiaoping Wang
# This version includes the following new keys:

# Apply correction for Goniometer x and z offset 
#
# z_offset  
# x_offset  

# User defined peak threshold
# split_threshold

# User Qmax for peak search and integration
# Qmax

#
# read_UB can read UB for individual runs obtained from previous integration.
#        
#Perform instrument background subtraction
# subtract_bkg  True

# adaptive_Q_background  
# adaptive_Q_multiplier
#   
# Peaks 16 pixels from detector edge are excluded.
# 
import os
import sys
import shutil
import time
import ReduceDictionary
import sys,os,glob,subprocess
import struct
from collections import defaultdict

sys.path.insert(0,"/opt/mantid42/bin")

from mantid.simpleapi import FilterBadPulses,ConvertToMD,LoadIsawUB,PredictPeaks,CentroidPeaks,SaveIsawUB,SaveIsawPeaks,IntegratePeaksMD,Rebin,PeakIntegration,IntegrateEllipsoids,CopySample,CropWorkspace,MDNormSCD,BinMD,AddSampleLog,SetGoniometer
from mantid.api import mtd

def write_prenexus(ws, instrument, run, output_directory):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    writeOrAppend = "ab"
    if not os.path.exists(os.path.join(output_directory, instrument+"_"+run+"_neutron_event.dat")):
      writeOrAppend = "wb"
    newFile = open(os.path.join(output_directory, instrument+"_"+run+"_neutron_event.dat"), writeOrAppend)
    data_dict = defaultdict(list)
    for i in range(ws.getNumberHistograms()):
        tofs = ws.getSpectrum(i).getTofs()
        pid = ws.getDetector(i).getID()
        for j in range(len(tofs)):
            data_dict[tofs[j]].append(pid)
    sorted(data_dict)
    events = []
    for tof in data_dict:
        pids = data_dict.get(tof)
        for pid in pids:
            events.append(int(10*tof))
            events.append(pid)
    
    newFile.write(struct.pack('%si' % len(events), *events))
    newFile.close()

start_time = time.time()

#
# Get the last config file name and the run number to process from the command line
# config file is last one written in current directory
#

config_files=glob.glob(os.path.join('.','*.config'))
config_file_name = ""
if config_files!=[]:
    config_file_name=max(config_files,key=os.path.getmtime)


#
# Load the parameter names and values from the specified configuration file 
# into a dictionary and set all the required parameters from the dictionary.
#

params_dictionary = ReduceDictionary.LoadDictionary(config_file_name )

instrument_name           = params_dictionary[ "instrument_name" ]
output_directory          = params_dictionary[ "output_directory" ]
run = str(input.getRunNumber())

write_prenexus(input, instrument_name, run, output_directory)
