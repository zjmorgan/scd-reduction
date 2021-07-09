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

from mantid.simpleapi import ConvertToMD,LoadIsawUB,PredictPeaks,CentroidPeaks,SaveIsawUB,SaveIsawPeaks,IntegratePeaksMD,Rebin,PeakIntegration,IntegrateEllipsoids,CopySample,CropWorkspace,MDNormSCD,BinMD,AddSampleLog,SetGoniometer
from mantid.api import mtd

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
event_ws =CloneWorkspace('tmp')
run = str(event_ws.getRunNumber())

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

#
# Name the files to write for this run
#
run_niggli_matrix_file = output_directory + "/" + run + "_Niggli.mat"
run_niggli_integrate_file = output_directory + "/" + run + "_Niggli.integrate"

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
phi = event_ws.getRun()['phi'].value[0]
#
SetGoniometer(Workspace='event_ws', Goniometers='Universal')
print('TOPAZ_%s Events =%20d Phi =%9.4f Chi = %9.4f Omega = %9.4f\n'%(run, event_ws.getNumberEvents(), phi, chi_modified, omega_modified))

SaveNexus(InputWorkspace='event_ws', Filename=os.path.join(output_directory, instrument+"_"+run+".nxs.h5"))

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
#print "\n", run, " has integrated proton charge x 1000 of", proton_charge, "\n"

# Subtract no-sample background
if subtract_bkg is True:
  print('\nPerform instrument background subtraction using no sample data :')
  print('\n%s'%no_sample_event_nxs_fname)
  #Load no-sample event nexus 
  if mtd.doesExist('bkg')==False: 
    bkg = LoadEventNexus( Filename=no_sample_event_nxs_fname,
                      FilterByTofMin=min_tof, FilterByTofMax=max_tof)
    bkg_proton_charge = bkg.getRun().getProtonCharge()* 1000.0 
    #
    LoadIsawDetCal( bkg, 
                Filename=calibration_file_1, Filename2=calibration_file_2 )  

    #
  #Subtract instrument background
  temp=bkg*(proton_charge/bkg_proton_charge) *(float(BN_aperture_size)/2.0)  
  # Ratio of background counts to 2 mm diameter BN aperture
  event_ws=Minus(LHSWorkspace=event_ws, RHSWorkspace=temp)


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
peaks_ws = FindPeaksMD( MDEW, MaxPeaks=num_peaks_to_find,DensityThresholdFactor='5', 
                        PeakDistanceThreshold=distance_threshold )

#Remove peaks on detector edge
peaks_on_edge=[]
for i in range(peaks_ws.getNumberPeaks()):
  pi=peaks_ws.getPeak(i)
  if pi.getRow()<16 or pi.getRow()>240 or pi.getCol()<16 or pi.getCol()>240:
      peaks_on_edge.append(i)
DeleteTableRows(TableWorkspace=peaks_ws,Rows=peaks_on_edge)
#
# Read or find UB for the run
# Read or find UB for the run
try:
  if read_UB:
    # Read orientation matrix from file
    ubpath=os.path.dirname(UB_filename)
    ubrunnum = run
    if os.path.exists(ubpath +'%s_Niggli.mat'%(run)):
      LoadIsawUB(InputWorkspace=peaks_ws, Filename=ubpath +'%s_Niggli.mat'%(run))
      print('Use UB: ', ubpath +'%s_Niggli.mat'%(run))
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
  
  print(peaks_ws.sample().getOrientedLattice())
  indexed=IndexPeaks( PeaksWorkspace=peaks_ws, Tolerance=tolerance)
  print(("Number of Indexed Peaks: {:d}".format(indexed[0])))
  
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
  AnalysisDataService.remove( MDEW.name() )
  
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
  print('\nReduced run ' + str(run) + ' in ' + str(end_time - start_time) + ' sec')
  print('using config file ' + config_file_name) 
  
  
  
  import mantid.plots.helperfunctions as hf
  import numpy as np
  from configparser import ConfigParser
  import matplotlib
  import matplotlib.pyplot as plt
  from matplotlib import colors
  
  def get_colorscale_minimum(arr):
      x=arr[np.isfinite(arr)]
      x=x[x>0]
      xc=x[np.argsort(x)][len(x)*0.02] #ignore the bottom 2%
      return xc
  
  def dim2array(d,center=True):
      """
      Create a numpy array containing bin centers along the dimension d
      input: d - IMDDimension
      return: numpy array, from min+st/2 to max-st/2 with step st  
      """
      dmin=d.getMinimum()
      dmax=d.getMaximum()
      dstep=d.getX(1)-d.getX(0)
      if center:
          return np.arange(dmin+dstep/2,dmax,dstep)
      else:
          return np.linspace(dmin,dmax,d.getNBins()+1)
  
  def Plot2DMD(ax,ws,NumEvNorm=False,**kwargs):
      """
      Plot a 2D slice from an MDHistoWorkspace (assume all other dimensions are 1)
      input: ax - axis object
      input: ws - handle to the workspace
      input: NumEvNorm - must be set to true if data was converted to MD from a histo workspace (like NXSPE) and no MDNorm... algorithms were used
      input: kwargs - arguments that are passed to plot, such as plotting symbol, color, label, etc.
      """
      dims=ws.getNonIntegratedDimensions()
      if len(dims)!=2:
          raise ValueError("The workspace dimensionality is not 2")
      dimx=dims[0]
      x=dim2array(dimx,center=False)
      dimy=dims[1]
      y=dim2array(dimy,center=False)
      intensity=ws.getSignalArray()*1.
      if NumEvNorm:
          nev=ws.getNumEventsArray()
          intensity/=nev
      intensity=intensity.squeeze()
      intensity=np.ma.masked_where(np.isnan(intensity),intensity)
      XX,YY=np.meshgrid(x,y,indexing='ij')
      pcm=ax.pcolorfast(XX,YY,intensity,**kwargs)
      ax.set_xlabel(dimx.getName())
      ax.set_ylabel(dimy.getName())
      ax.set(aspect=1)
      for tick in ax.get_xticklabels():
          tick.set_rotation(90)
      return pcm
  
  
  #try to do HKL plots
  plot_config_files=glob.glob(os.path.join('.','*.ini'))
  plot_latest_config = ""
  if plot_config_files!=[]:
      plot_latest_config=max(plot_config_files,key=os.path.getmtime)
  plot_config_file=plot_latest_config
  try:
      with open(plot_config_file) as fh:
          config_data=fh.read().decode('utf-8')
  
      config = ConfigParser(allow_no_value=True)
      config.read_string(config_data)
  
      sa_filename = config['NORMALIZATION']['safile']
      flux_filename = config['NORMALIZATION']['fluxfile']
      cal_filename = config['REDUCTION']['calfile']
      ub_matrix_dir = config['REDUCTION']['ubdirectory']
      plot_params=[]
      for k in list(config.keys()):
          if 'PLOT' in k:
              plot_param={}
              for item in list(config[k].items()):
                  plot_param[item[0]] = item[1]
              plot_params.append(plot_param)
  except Exception as e:
      logger.warning("Could not read plotting configuration: "+str(e))
      sys.exit()
  
  #copy the UB
  CopySample(InputWorkspace = peaks_ws,OutputWorkspace = event_ws,CopyName = '0',CopyMaterial = '0',CopyEnvironment = '0',CopyShape = '0')
  event_ws=ConvertUnits(event_ws,Target='Momentum')
  have_van=True
  try:
      sa=Load(sa_filename)
      fl=Load(flux_filename)
  except:
      have_van=False
  XMin=1.
  XMax=10.
  if have_van:
      MaskDetectors(event_ws, MaskedWorkspace=sa)
      XMin = sa.getXDimension().getMinimum()
      XMax = sa.getXDimension().getMaximum()
  event_ws = CropWorkspace(InputWorkspace=event_ws,XMin=XMin,XMax=XMax)
  event_ws=Rebin(event_ws,Params="{},{},{}".format(XMin,XMax-XMin,XMax),PreserveEvents=True)
  mde = ConvertToMD(InputWorkspace=event_ws,
                    QDimensions='Q3D',
                    dEAnalysisMode='Elastic',
                    Q3DFrames='HKL',
                    QConversionScales='HKL')
  for plot_param in plot_params:
      AxisNames={'H':'[H,0,0]','K':'[0,K,0]','L':'[0,0,L]'}
      xaxis=mde.getDimension(mde.getDimensionIndexByName(AxisNames[plot_param['axis1']]))
      try:
          plot_xmin=float(plot_param['xmin'])
      except ValueError:
          plot_xmin=xaxis.getMinimum()
      try:
          plot_xmax=float(plot_param['xmax'])
      except ValueError:    
          plot_xmax=xaxis.getMaximum()
      try:
          plot_xbins=int(plot_param['xsteps'])
      except ValueError:
          plot_xbins=400
      AlignedDim0='{},{},{},{}'.format(AxisNames[plot_param['axis1']],plot_xmin,plot_xmax,plot_xbins)
      yaxis=mde.getDimension(mde.getDimensionIndexByName(AxisNames[plot_param['axis2']]))
      try:
          plot_ymin=float(plot_param['ymin'])
      except ValueError:
          plot_ymin=yaxis.getMinimum()
      try:
          plot_ymax=float(plot_param['ymax'])
      except ValueError:    
          plot_ymax=yaxis.getMaximum()
      try:
          plot_ybins=int(plot_param['ysteps'])
      except ValueError:
          plot_ybins=400
      AlignedDim1='{},{},{},{}'.format(AxisNames[plot_param['axis2']],plot_ymin,plot_ymax,plot_ybins)
      zaxis=mde.getDimension(mde.getDimensionIndexByName(AxisNames[plot_param['axis3']]))
      try:
          plot_zmin=float(plot_param['zmin'])
      except ValueError:
          plot_zmin=zaxis.getMinimum()
      try:
          plot_zmax=float(plot_param['zmax'])
      except ValueError:    
          plot_zmax=zaxis.getMaximum()
      plot_zbins=1
      AlignedDim2='{},{},{},{}'.format(AxisNames[plot_param['axis3']],plot_zmin,plot_zmax,plot_zbins)
      if have_van:
          data,norm = MDNormSCD(InputWorkspace=mde,
                                FluxWorkspace=fl,
                                SolidAngleWorkspace=sa,
                                SkipSafetyCheck=True,
                                AlignedDim0=AlignedDim0,
                                AlignedDim1=AlignedDim1,
                                AlignedDim2=AlignedDim2)
          mdh = data/norm 
      else:
          mdh = BinMD(InputWorkspace=mde,
                      AlignedDim0=AlignedDim0,
                      AlignedDim1=AlignedDim1,
                      AlignedDim2=AlignedDim2,
                      AxisAligned=True)
  
     
      intensity=mdh.getSignalArray()
      vmin = 1
      vmax = intensity.max()
      if vmax > 1:
        fig, ax = plt.subplots()
        ax.set_title('{}=[{},{}]'.format(AxisNames[plot_param['axis3']],plot_zmin,plot_zmax))
        logNorm = colors.LogNorm(vmin = vmin, vmax = vmax)
        cm = plt.cm.get_cmap('rainbow')
        pcm = Plot2DMD(ax,mdh,NumEvNorm=False,norm=logNorm)
        fig.colorbar(pcm,ax=ax)
        #plotSlice(mdws, xydim=[plot_param['axis1'],plot_param['axis2']], slicepoint=[0,0,0.5*(plot_param['zmin']+plot_param['zmax'])], colorscalelog=True)
        plt.show()
except Exception as e:
    logger.error("Could not find UB: "+str(e))
