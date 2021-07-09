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

#
# Add IntegratePeakProfileFitting, Jan 12, 2020, X. P. Wang
#

#
# Add option for cryogoniometer with 1 axis of rotation in omega, March 2, 2020, X. P. Wang
# 

import os
import sys
import shutil
import time
import ReduceDictionary
sys.path.insert(0,"/opt/mantid50/bin")
sys.path.insert(0,"/opt/mantid50/lib")

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
recenter                  = params_dictionary.get("recenter",'False')
use_sphere_integration    = params_dictionary[ "use_sphere_integration" ]
use_ellipse_integration   = params_dictionary[ "use_ellipse_integration" ]
use_fit_peaks_integration = params_dictionary[ "use_fit_peaks_integration" ]
use_cylindrical_integration = params_dictionary[ "use_cylindrical_integration" ]

peak_radius               = params_dictionary[ "peak_radius" ]
bkg_inner_radius          = params_dictionary[ "bkg_inner_radius" ]
bkg_outer_radius          = params_dictionary[ "bkg_outer_radius" ]
integrate_if_edge_peak    = params_dictionary[ "integrate_if_edge_peak" ]

#IntegratePeaksProfileFitting parameter
moderator_coefficients_file = params_dictionary.get("moderator_coefficients_file","/SNS/TOPAZ/shared/ProfileFitting/franz_coefficients_2017.dat")
minppl_frac                 = float(params_dictionary.get("minppl_frac",'0.9'))
maxppl_frac                 =float(params_dictionary.get("maxppl_frac",'1.1'))
frac_stop                   =float(params_dictionary.get("frac_stop",'0.06')) 
intensity_cutoff            =int(params_dictionary.get("intensity_cutoff",'50')) 
dq_max                      = float(params_dictionary.get("dq_max",'0.16')) 
#
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

Qmax                       = params_dictionary.get('Qmax', "16")

maxQ ='%s,%s,%s'%(Qmax,Qmax,Qmax)
minQ ='-%s,-%s,-%s'%(Qmax,Qmax,Qmax)


# Get the fully qualified input run file name, either from a specified data 
# directory or from findnexus
#
if os.path.isfile(str(run)):
    full_name = str(run)
else:
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
run_niggli_fit_integrate_file = output_directory + "/" + run + "_Niggli_fit.integrate"
run_niggli_fit_profile_file = output_directory + "/" + run + "_Niggli_profile.nxs"

#
# Load the run data and find the total monitor counts
#
event_ws = Load( Filename=full_name, 
                           FilterByTofMin=min_tof, FilterByTofMax=max_tof)
event_ws = FilterBadPulses(InputWorkspace=event_ws, LowerCutoff = 75)
#
#Check for goniometer settings
if event_ws.getRun()['gonio'].value[0] is None:
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

else:
  gonio_type =event_ws.getRun()['gonio'].value[0]
  if gonio_type==2:
    #MaskBTP(Workspace=event_ws,Bank="1-12,15,21,24-15,31-32,34-35,41-45,51-55")
    SetGoniometer(event_ws,Axis0="omega:Axis1,0,1,0,1")

#
event_ws = FilterBadPulses(InputWorkspace=event_ws, LowerCutoff = 75)
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
  #Subtract instrument background scaled down by 95%
  temp=0.95*bkg*(proton_charge/bkg_proton_charge) *(float(BN_aperture_size)/2.0)  
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
peaks_ws = FindPeaksMD( MDEW, MaxPeaks=500,DensityThresholdFactor=100, 
                        PeakDistanceThreshold=distance_threshold,
                        EdgePixels=n_bad_edge_pixels)

#
# Read or find UB for the run
if read_UB:
  # Read orientation matrix from file
  ubpath=os.path.dirname(UB_filename)
  print(('Use UB in %s'%(ubpath)))
  ubrunnum = run
  if os.path.exists(ubpath +'%s_Niggli.mat'%(run)):
    LoadIsawUB(InputWorkspace=peaks_ws, Filename=ubpath +'%s_Niggli.mat'%(run))
    print(('Use UB: ', ubpath +'%s_Niggli.mat'%(run)))
    IndexPeaks( PeaksWorkspace=peaks_ws, CommonUBForAll=True, Tolerance=tolerance)
    FindUBUsingIndexedPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance)

  else:
    LoadIsawUB( InputWorkspace=peaks_ws, Filename=UB_filename)
    IndexPeaks( PeaksWorkspace=peaks_ws, CommonUBForAll=False, Tolerance=tolerance)
    FindUBUsingIndexedPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance)
    IndexPeaks( PeaksWorkspace=peaks_ws, CommonUBForAll=True, Tolerance=tolerance)

else:
  # Find a Niggli UB matrix that indexes the peaks in this run
  FindUBUsingFFT( PeaksWorkspace=peaks_ws, MinD=min_d, MaxD=max_d, Tolerance=tolerance )
  IndexPeaks( PeaksWorkspace=peaks_ws, CommonUBForAll=False, Tolerance=tolerance)
  FindUBUsingIndexedPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance)
  IndexPeaks( PeaksWorkspace=peaks_ws, CommonUBForAll=True, Tolerance=tolerance)

#
# Save UB and peaks file from find peaks
#
SaveIsawUB( InputWorkspace=peaks_ws,Filename=run_niggli_matrix_file )
SaveIsawPeaks( InputWorkspace=peaks_ws, AppendFile=False,
               Filename=run_niggli_integrate_file )

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
                EdgePixels=n_bad_edge_pixels,
                ReflectionCondition='Primitive' )

  #
  #Find peak centroids from predicted peak position on detector face in event workspace
  #peaks_ws = CentroidPeaks(InPeaksWorkspace=peaks_ws, InputWorkspace=event_ws, PeakRadius=8, EdgePixels=n_bad_edge_pixels)
  #Second option to find peak centroids from predicted peak position on detector face in MD workspace
  if recenter:
    CentroidPeaksMD(InputWorkspace=MDEW, PeakRadius=bkg_inner_radius, PeaksWorkspace=peaks_ws, OutputWorkspace=peaks_ws)

  #
  #Refine UB and reindex peaks
  IndexPeaks(PeaksWorkspace=peaks_ws, CommonUBForAll=True, Tolerance=tolerance, RoundHKLs=False)
  FindUBUsingIndexedPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance)
  IndexPeaks(PeaksWorkspace=peaks_ws, CommonUBForAll=True, Tolerance=tolerance, RoundHKLs=False)

  #
  #Initial pass to get an estimate of peak intens using ellipse integration
  peaks_ws= IntegrateEllipsoids(InputWorkspace=event_ws, 
                             PeaksWorkspace = peaks_ws,
                             IntegrateInHKL=integrate_in_HKL_space,
                             RegionRadius = ellipse_region_radius,
                             SpecifySize = ellipse_size_specified,
                             PeakSize = peak_radius,
                             BackgroundOuterSize = bkg_outer_radius,
                             BackgroundInnerSize = bkg_inner_radius,
                             AdaptiveQBackground=adaptive_Q_background,
                             AdaptiveQMultiplier=adaptive_Q_multiplier,
                             UseOnePercentBackgroundCorrection = False)

else:
  print("Only integrating FOUND peaks ....")
  #Search for more peaks using split_threshold
  peaks_ws = FindPeaksMD( MDEW, MaxPeaks=num_peaks_to_find,
                        DensityThresholdFactor=split_threshold, 
                        PeakDistanceThreshold=distance_threshold,
                        EdgePixels=n_bad_edge_pixels )
  #
  LoadIsawUB( InputWorkspace=peaks_ws, Filename=run_niggli_matrix_file)
  IndexPeaks( PeaksWorkspace=peaks_ws, CommonUBForAll=True, Tolerance=tolerance)
  FindUBUsingIndexedPeaks(PeaksWorkspace=peaks_ws, Tolerance=tolerance)
  IndexPeaks( PeaksWorkspace=peaks_ws, CommonUBForAll=True, Tolerance=tolerance)
  peaks_ws= FilterPeaks(InputWorkspace=peaks_ws, FilterVariable='h^2+k^2+l^2', FilterValue=0, Operator='>')

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
#  event_ws = Rebin( InputWorkspace=event_ws,
#                    Params=rebin_params, PreserveEvents=preserve_events )
#  peaks_ws = PeakIntegration( InPeaksWorkspace=peaks_ws, InputWorkspace=event_ws, 
#                              IkedaCarpenterTOF=use_ikeda_carpenter,
#                              MatchingRunNo=True,
#                              NBadEdgePixels=n_bad_edge_pixels )
#
  #Prepare MD workspace for 3D profile fitting
  Qlab_md=ConvertToMD(InputWorkspace=event_ws, 
    QDimensions='Q3D', 
    dEAnalysisMode='Elastic', 
    Q3DFrames='Q_lab', 
    LorentzCorrection=False, 
    MinValues=minQ, 
    MaxValues=maxQ, 
    SplitInto='2', 
    SplitThreshold=split_threshold, MinRecursionDepth=7,MaxRecursionDepth='11')
  #
  mtd.remove('event_ws')

  #Perform 3D profile fitting
  IntegratePeaksProfileFitting(OutputPeaksWorkspace='peaks_out_%s'%(run), OutputParamsWorkspace='parameters_%s'%(run), 
    InputWorkspace=Qlab_md, 
    PeaksWorkspace=peaks_ws,
    UBFile=run_niggli_matrix_file,
    ModeratorCoefficientsFile= moderator_coefficients_file,
    MinpplFrac=minppl_frac, 
    MaxpplFrac=maxppl_frac, 
    FracStop=frac_stop, 
    IntensityCutoff=intensity_cutoff, 
    DQMax=dq_max,
    EdgeCutoff=n_bad_edge_pixels,
    StrongPeakParamsFile=None)

  #
  SaveIsawPeaks(InputWorkspace='peaks_out_%s'%(run), AppendFile=False, 
               Filename=output_directory + '/' + run + '_Niggli_fit.integrate')
  SaveNexus(InputWorkspace='parameters_%s'%(run), Filename=output_directory + '/' + run + '_Niggli_profile.nxs')

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
if (not use_fit_peaks_integration):
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

# ============
#
# If requested, also switch to the specified conventional cell and save the
# corresponding matrix and integrate file
#
if not use_cylindrical_integration:
  if (not cell_type is None) and (not centering is None) :
    conv_name = output_directory + "/" + run + "_" + cell_type + "_" + centering
    conventional_integrate_file = conv_name + ".integrate"
    conventional_matrix_file = conv_name + ".mat"
    SelectCellOfType( PeaksWorkspace=peaks_ws, CellType=cell_type, Centering=centering,
                      Apply=True, Tolerance=tolerance )
    # Least squares refinement of conventional cell after transfromation
    if (cell_type=='Rhombohedral' and centering=='R'):
      cell_type = 'Hexagonal'  # The R lattice is hexagonal after unit cell transformation
    #elif (cell_type=='Monoclinic'):
    #  cell_type = 'Monoclinic ( b unique )'
    OptimizeLatticeForCellType(PeaksWorkspace=peaks_ws, CellType=cell_type,Apply='1', Tolerance=tolerance)
    OptimizeCrystalPlacement(PeaksWorkspace=peaks_ws,ModifiedPeaksWorkspace=peaks_ws,
      FitInfoTable='CrystalPlacement_info',MaxIndexingError=tolerance)
    SaveIsawPeaks( InputWorkspace=peaks_ws, AppendFile=False, Filename=conventional_integrate_file )
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


#+++++++++++++

end_time = time.time()
print(('\nReduced run ' + str(run) + ' in ' + str(end_time - start_time) + ' sec'))
print(('using config file ' + config_file_name)) 

#
# Try to get this to terminate when run by topaz_reduction.py, from NX session
#
sys.exit(0)

