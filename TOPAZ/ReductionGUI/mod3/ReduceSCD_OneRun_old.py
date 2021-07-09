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
sys.path.insert(0,"/opt/mantidnightly/bin")
#sys.path.insert(0,"/opt/mantidnightly/release/bin")

from mantid.simpleapi import *
from mantid.api import *

print("API Version")
print((apiVersion()))

start_time = time.time()

#
# Get the config file name and the run number to process from the command line
#
if (len(sys.argv) < 3):
  print("You MUST give the config file name and run number on the command line")
  exit(0)

config_file_name = sys.argv[1]
run              = sys.argv[2]

#
# Load the parameter names and values from the specified configuration file 
# into a dictionary and set all the required parameters from the dictionary.
#
config_file_name = os.path.dirname(os.path.realpath(__file__)) +'/' + config_file_name 


params_dictionary = ReduceDictionary.LoadDictionary(config_file_name )

instrument_name           = params_dictionary[ "instrument_name" ]
calibration_file_1        = params_dictionary[ "calibration_file_1" ]
calibration_file_2        = params_dictionary[ "calibration_file_2" ]
data_directory            = params_dictionary[ "data_directory" ]
output_directory          = params_dictionary[ "output_directory" ]
min_tof                   = params_dictionary[ "min_tof" ] 
max_tof                   = params_dictionary[ "max_tof" ]
use_monitor_counts        = params_dictionary[ "use_monitor_counts" ]
min_monitor_tof           = params_dictionary[ "min_monitor_tof" ] 
max_monitor_tof           = params_dictionary[ "max_monitor_tof" ] 
monitor_index             = params_dictionary[ "monitor_index" ] 
cell_type                 = params_dictionary[ "cell_type" ] 
centering                 = params_dictionary[ "centering" ]
num_peaks_to_find         = params_dictionary[ "num_peaks_to_find" ]
min_d                     = params_dictionary[ "min_d" ]
max_d                     = params_dictionary[ "max_d" ]
tolerance                 = params_dictionary[ "tolerance" ]
integrate_predicted_peaks = params_dictionary[ "integrate_predicted_peaks" ]
min_pred_wl               = params_dictionary[ "min_pred_wl" ]
max_pred_wl               = params_dictionary[ "max_pred_wl" ]
min_pred_dspacing         = params_dictionary[ "min_pred_dspacing" ]
max_pred_dspacing         = params_dictionary[ "max_pred_dspacing" ]

use_sphere_integration    = params_dictionary[ "use_sphere_integration" ]
use_ellipse_integration   = params_dictionary[ "use_ellipse_integration" ]
use_fit_peaks_integration = params_dictionary[ "use_fit_peaks_integration" ]
use_cylindrical_integration = params_dictionary[ "use_cylindrical_integration" ]

peak_radius               = params_dictionary[ "peak_radius" ]
bkg_inner_radius          = params_dictionary[ "bkg_inner_radius" ]
bkg_outer_radius          = params_dictionary[ "bkg_outer_radius" ]
integrate_if_edge_peak    = params_dictionary[ "integrate_if_edge_peak" ]

rebin_step                = params_dictionary[ "rebin_step" ]
preserve_events           = params_dictionary[ "preserve_events" ] 
use_ikeda_carpenter       = params_dictionary[ "use_ikeda_carpenter" ]
n_bad_edge_pixels         = params_dictionary[ "n_bad_edge_pixels" ]

rebin_params = min_tof + "," + rebin_step + "," + max_tof

ellipse_region_radius     = params_dictionary[ "ellipse_region_radius" ]
ellipse_size_specified    = params_dictionary[ "ellipse_size_specified" ]

cylinder_radius           = params_dictionary[ "cylinder_radius" ]
cylinder_length           = params_dictionary[ "cylinder_length" ]

read_UB                   = params_dictionary[ "read_UB" ]
UB_filename               = params_dictionary[ "UB_filename" ]

x_offset                  = params_dictionary.get("x_offset",'0.0')
z_offset                  = params_dictionary.get("z_offset",'0.0')
adaptive_Q_background   = params_dictionary.get("adaptive_Q_background",'False')
adaptive_Q_multiplier   = params_dictionary.get("adaptive_Q_multiplier",'0.000')

integrate_in_HKL_space     = params_dictionary[ "integrate_in_HKL_space" ]

subtract_bkg               = params_dictionary[ "subtract_bkg" ]
no_sample_event_nxs_fname  = params_dictionary[ "no_sample_event_nxs_fname" ]
BN_aperture_size           = params_dictionary[ "BN_aperture_size" ]

split_threshold            = params_dictionary[ "split_threshold" ]

Qmax                       = params_dictionary.get('Qmax', "20")

maxQ ='%s,%s,%s'%(Qmax,Qmax,Qmax)
minQ ='-%s,-%s,-%s'%(Qmax,Qmax,Qmax)


# Get the fully qualified input run file name, either from a specified data 
# directory or from findnexus
#
short_filename = "%s_%s" % (instrument_name, str(run))
if data_directory is not None:
  full_name = data_directory + "/" + short_filename + ".nxs.h5"
  if not os.path.exists(full_name):
    full_name = data_directory + "/" + short_filename + "_event.nxs"
else:
  full_name = short_filename

print(("\nProcessing File: " + full_name + " ......\n"))

#
# Name the files to write for this run
#
run_niggli_matrix_file = output_directory + "/" + run + "_Niggli.mat"
run_niggli_integrate_file = output_directory + "/" + run + "_Niggli.integrate"

#
# Load the run data and find the total monitor counts
#
event_ws = Load( Filename=full_name, 
                           FilterByTofMin=min_tof, FilterByTofMax=max_tof)

#
omega = event_ws.getRun()['omega'].value[0]
omega_modified = float(omega) - float(z_offset)
AddSampleLog(Workspace=event_ws, LogName='omega', 
    LogText='%.4f'%omega_modified, LogType='Number Series')
AddSampleLog(Workspace=event_ws, LogName='omegaRequest', 
    LogText='%.4f'%omega_modified, LogType='Number Series')
#
chi = event_ws.getRun()['chi'].value[0]
chi_modified = float(chi) - float(x_offset)
AddSampleLog(Workspace=event_ws, LogName='chi', 
    LogText='%.4f'%chi_modified, LogType='Number Series')
#
SetGoniometer(Workspace='event_ws', Goniometers='Universal')
print(('TOPAZ_%s \nOmega =%9.4f, Omega_modified = %s \nChi = %9.4f, Chi_modified = %s\n'%(run, omega, omega_modified, chi, chi_modified)))

#
event_ws = FilterBadPulses(InputWorkspace=event_ws, LowerCutoff = 95)
#
# Load calibration file(s) if specified.  NOTE: The file name passed in to LoadIsawDetCal
# can not be None.  TOPAZ has one calibration file, but SNAP may have two.
#
if (calibration_file_1 is not None ) or (calibration_file_2 is not None):
  if (calibration_file_1 is None ):
    calibration_file_1 = ""
  if (calibration_file_2 is None ):
    calibration_file_2 = ""
  LoadIsawDetCal( event_ws, 
                  Filename=calibration_file_1, Filename2=calibration_file_2 )  

proton_charge = event_ws.getRun().getProtonCharge() * 1000.0  # get proton charge
print(("\n", run, " has integrated proton charge x 1000 of", proton_charge, "\n"))

try:
  monitor_ws = LoadNexusMonitors( Filename=full_name )
  integrated_monitor_ws = Integration( InputWorkspace=monitor_ws, 
                                       RangeLower=min_monitor_tof, RangeUpper=max_monitor_tof, 
                                       StartWorkspaceIndex=monitor_index, EndWorkspaceIndex=monitor_index )
  
  monitor_count = integrated_monitor_ws.dataY(0)[0]
except:
  monitor_count = 0
print(("\n", run, " has integrated monitor count", monitor_count, "\n"))

# Subtract no-sample background
if subtract_bkg is True:
  print('\nPerform instrument background subtraction using no sample data :')
  print(('\n%s'%no_sample_event_nxs_fname))
  #Load no-sample event nexus 
  if mtd.doesExist('bkg')==False: 
    bkg = LoadEventNexus( Filename=no_sample_event_nxs_fname,
                      FilterByTofMin=min_tof, FilterByTofMax=max_tof)
    bkg_proton_charge = bkg.getRun().getProtonCharge()* 1000.0 
    bkg = FilterBadPulses(InputWorkspace=bkg, LowerCutoff = 95)
    #
    LoadIsawDetCal( bkg, 
                Filename=calibration_file_1, Filename2=calibration_file_2 )  

    #
  #Subtract instrument background
  temp=bkg*(proton_charge/bkg_proton_charge) *(float(BN_aperture_size)/2.0)  
  # Ratio of background counts to 2 mm diameter BN aperture
  event_ws=Minus(LHSWorkspace=event_ws, RHSWorkspace=temp)


print(("Qmax = %s"%maxQ, "\n"))
#
# Make MD workspace using Lorentz correction, to find peaks 
#
MDEW = ConvertToMD( InputWorkspace=event_ws, QDimensions="Q3D",
                    dEAnalysisMode="Elastic", QConversionScales="Q in A^-1",
   	            LorentzCorrection='1', MinValues=minQ, MaxValues=maxQ,
                    SplitInto='2', SplitThreshold=split_threshold,MaxRecursionDepth='11' )
#
# Find the requested number of peaks.  Once the peaks are found, we no longer
# need the weighted MD event workspace, so delete it.
#
distance_threshold = 0.9 * 6.28 / float(max_d)
peaks0_ws = FindPeaksMD( MDEW, MaxPeaks=50,DensityThresholdFactor='5', 
                        PeakDistanceThreshold=distance_threshold )

#Remove peaks on detector edge
peaks_on_edge=[]
for i in range(peaks0_ws.getNumberPeaks()):
  pi=peaks0_ws.getPeak(i)
  if pi.getRow()<16 or pi.getRow()>240 or pi.getCol()<16 or pi.getCol()>240:
      peaks_on_edge.append(i)
DeleteTableRows(TableWorkspace=peaks0_ws,Rows=peaks_on_edge)
#
# Read or find UB for the run
# Read or find UB for the run
if read_UB:
  # Read orientation matrix from file
  ubpath=os.path.dirname(UB_filename)
  ubrunnum = run
  if os.path.exists(ubpath +'%s_Niggli.mat'%(run)):
    LoadIsawUB(InputWorkspace=peaks0_ws, Filename=ubpath +'%s_Niggli.mat'%(run))
    print(('Use UB: ', ubpath +'%s_Niggli.mat'%(run)))
    IndexPeaks( PeaksWorkspace=peaks0_ws, CommonUBForAll=True, Tolerance=tolerance)
    FindUBUsingIndexedPeaks(PeaksWorkspace=peaks0_ws, Tolerance=tolerance)

  else:
    LoadIsawUB( InputWorkspace=peaks0_ws, Filename=UB_filename)
    IndexPeaks( PeaksWorkspace=peaks0_ws, CommonUBForAll=False, Tolerance=tolerance)
    FindUBUsingIndexedPeaks(PeaksWorkspace=peaks0_ws, Tolerance=tolerance)
    IndexPeaks( PeaksWorkspace=peaks0_ws, CommonUBForAll=True, Tolerance=tolerance)

else:
  # Find a Niggli UB matrix that indexes the peaks in this run
  FindUBUsingFFT( PeaksWorkspace=peaks0_ws, MinD=min_d, MaxD=max_d, Tolerance=tolerance )
  IndexPeaks( PeaksWorkspace=peaks0_ws, CommonUBForAll=False, Tolerance=tolerance)
  
peaks_ws = FindPeaksMD( MDEW, MaxPeaks=num_peaks_to_find,DensityThresholdFactor='5', 
                        PeakDistanceThreshold=distance_threshold )
CopySample(InputWorkspace = peaks0_ws,OutputWorkspace = peaks_ws,CopyName = '0',CopyMaterial = '0',CopyEnvironment = '0',CopyShape = '0')
IndexPeaks( PeaksWorkspace=peaks_ws, Tolerance=tolerance)
IndexPeaks( PeaksWorkspace=peaks_ws, CommonUBForAll=True, Tolerance=tolerance)
#
# Get complete list of peaks to be integrated and load the UB matrix into
# the predicted peaks workspace, so that information can be used by the
# PeakIntegration algorithm.
#
if integrate_predicted_peaks:
  print("PREDICTING peaks to integrate....")
  peaks_ws = PredictPeaks( InputWorkspace=peaks_ws,
                WavelengthMin=min_pred_wl, WavelengthMax=max_pred_wl,
                MinDSpacing=min_pred_dspacing, MaxDSpacing=max_pred_dspacing, 
                ReflectionCondition='Primitive' )
  #Remove peaks on detector edge
  peaks_on_edge=[]
  for i in range(peaks_ws.getNumberPeaks()):
    pi=peaks_ws.getPeak(i)
    if pi.getRow()<16 or pi.getRow()>240 or pi.getCol()<16 or pi.getCol()>240:
        peaks_on_edge.append(i)
  DeleteTableRows(TableWorkspace=peaks_ws,Rows=peaks_on_edge)
  #
  #Find peak centroids from predicted peak position on detector face in event workspace
  peaks_ws = CentroidPeaks(InPeaksWorkspace=peaks_ws, InputWorkspace=event_ws, PeakRadius=4, EdgePixels=16)
  #Second option to find peak centroids from predicted peak position on detector face in MD workspace
  #CentroidPeaksMD(InputWorkspace=MDEW,, PeakRadius=tolerance, PeaksWorkspace=peaks_ws)

  #
  #Refine UB and reindex peaks
  FindUBUsingIndexedPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance)
  IndexPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance, RoundHKLs=False)
  FindUBUsingIndexedPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance)

else:
  print("Only integrating FOUND peaks ....")

#
#Delete MD workspace
AnalysisDataService.remove( MDEW.getName() )

#
# Save UB and peaks file, so if something goes wrong latter, we can at least 
# see these partial results
#
SaveIsawUB( InputWorkspace=peaks_ws,Filename=run_niggli_matrix_file )
SaveIsawPeaks( InputWorkspace=peaks_ws, AppendFile=False,
               Filename=run_niggli_integrate_file )

#
# Set the monitor counts for all the peaks that will be integrated
#
#
# Set the monitor counts for all the peaks that will be integrated
#
num_peaks = peaks_ws.getNumberPeaks()
for i in range(num_peaks):
  peak = peaks_ws.getPeak(i)
  if use_monitor_counts:
    peak.setMonitorCount( monitor_count )
  else:
    peak.setMonitorCount( proton_charge )
if use_monitor_counts:
  print('\n*** Beam monitor counts used for scaling.')
else:
  print('\n*** Proton charge x 1000 used for scaling.\n')
    
if use_sphere_integration:
#
# Integrate found or predicted peaks in Q space using spheres, and save 
# integrated intensities, with Niggli indexing.  First get an un-weighted 
# workspace to do raw integration (we don't need high resolution or 
# LorentzCorrection to do the raw sphere integration )
#
  MDEW = ConvertToMD( InputWorkspace=event_ws, QDimensions="Q3D",
                    dEAnalysisMode="Elastic", QConversionScales="Q in A^-1",
                    LorentzCorrection='0', MinValues=minQ, MaxValues=maxQ,
                    SplitInto='2', SplitThreshold=split_threshold,MaxRecursionDepth='11' )

  peaks_ws = IntegratePeaksMD( InputWorkspace=MDEW, PeakRadius=peak_radius,
                  CoordinatesToUse="Q (sample frame)",
	            BackgroundOuterRadius=bkg_outer_radius, 
                  BackgroundInnerRadius=bkg_inner_radius,
	            PeaksWorkspace=peaks_ws, 
                  IntegrateIfOnEdge=integrate_if_edge_peak )

elif use_fit_peaks_integration:
  event_ws = Rebin( InputWorkspace=event_ws,
                    Params=rebin_params, PreserveEvents=preserve_events )
  peaks_ws = PeakIntegration( InPeaksWorkspace=peaks_ws, InputWorkspace=event_ws, 
                              IkedaCarpenterTOF=use_ikeda_carpenter,
                              MatchingRunNo=True,
                              NBadEdgePixels=n_bad_edge_pixels )

elif use_ellipse_integration:
  peaks_ws= IntegrateEllipsoids( InputWorkspace=event_ws, 
                                 PeaksWorkspace = peaks_ws,
                                 IntegrateInHKL=integrate_in_HKL_space,
                                 RegionRadius = ellipse_region_radius,
                                 SpecifySize = ellipse_size_specified,
                                 PeakSize = peak_radius,
                                 BackgroundOuterSize = bkg_outer_radius,
                                 BackgroundInnerSize = bkg_inner_radius,
                                 AdaptiveQBackground=adaptive_Q_background,
                                 AdaptiveQMultiplier=adaptive_Q_multiplier)

elif use_cylindrical_integration:
  profiles_filename = output_directory + "/" + instrument_name + '_' + run + '.profiles'
  MDEW = ConvertToMD( InputWorkspace=event_ws, QDimensions="Q3D",
                    dEAnalysisMode="Elastic", QConversionScales="Q in A^-1",
                    LorentzCorrection='0', MinValues=minQ, MaxValues=maxQ,
                    SplitInto='2', SplitThreshold='100',MaxRecursionDepth='10' )

  peaks_ws = IntegratePeaksMD( InputWorkspace=MDEW, PeakRadius=cylinder_radius,
                  CoordinatesToUse="Q (sample frame)", 
                  Cylinder='1', CylinderLength = cylinder_length, 
                  PercentBackground = '20', ProfileFunction = 'NoFit',
                  ProfilesFile = profiles_filename,
	          PeaksWorkspace=peaks_ws, 
                  )
  if (not cell_type is None) and (not centering is None):
    print("WARNING: Cylindrical profiles are NOT transformed!!!")

#
# Save the final integrated peaks, using the Niggli reduced cell.  
# This is the only file needed, for the driving script to get a combined
# result.
#
SaveIsawPeaks( InputWorkspace=peaks_ws, AppendFile=False, 
               Filename=run_niggli_integrate_file )

# Print warning if user is trying to integrate using the cylindrical method and transorm the cell
if use_cylindrical_integration: 
  if (not cell_type is None) or (not centering is None):
    print("WARNING: Cylindrical profiles are NOT transformed!!!")
#
# If requested, also switch to the specified conventional cell and save the
# corresponding matrix and integrate file
#
if not use_ellipse_integration:
  if (not cell_type is None) and (not centering is None) :
    run_conventional_matrix_file = output_directory + "/" + run + "_" +    \
                                   cell_type + "_" + centering + ".mat"
    run_conventional_integrate_file = output_directory + "/" + run + "_" + \
                                      cell_type + "_" + centering + ".integrate"
    SelectCellOfType( PeaksWorkspace=peaks_ws, 
                      CellType=cell_type, Centering=centering, 
                      Apply=True, Tolerance=tolerance )
    SaveIsawPeaks( InputWorkspace=peaks_ws, AppendFile=False, 
                   Filename=run_conventional_integrate_file )
    SaveIsawUB( InputWorkspace=peaks_ws, Filename=run_conventional_matrix_file )

end_time = time.time()
print(('\nReduced run ' + str(run) + ' in ' + str(end_time - start_time) + ' sec'))
print(('using config file ' + config_file_name)) 

#
# Try to get this to terminate when run by ReduceSCD_Parallel.py, from NX session
#
sys.exit(0)

