Load(Filename='/SNS/TOPAZ/IPTS-22255/nexus/TOPAZ_35058.nxs.h5', OutputWorkspace='TOPAZ_35058.nxs_event', FilterByTofMin=500, FilterByTofMax=16666)
FilterBadPulses(InputWorkspace='TOPAZ_35058.nxs_event', OutputWorkspace='TOPAZ_35058.nxs_event', LowerCutoff=25)
LoadIsawDetCal(InputWorkspace='TOPAZ_35058.nxs_event', Filename='/SNS/TOPAZ/IPTS-22255/shared/2019Nov/calibration/TOPAZ_2019B_Nov_5.DetCal')
ConvertToMD(InputWorkspace='TOPAZ_35058.nxs_event', QDimensions='Q3D', dEAnalysisMode='Elastic', Q3DFrames='Q_sample', LorentzCorrection=True, OutputWorkspace='TOPAZ_35058.nxs_md', MinValues='-16,-16,-16', MaxValues='16,16,16', SplitInto='2', SplitThreshold=50, MaxRecursionDepth=13, MinRecursionDepth=7)
FindPeaksMD(InputWorkspace='TOPAZ_35058.nxs_md', PeakDistanceThreshold=0.12286956521739131, MaxPeaks=400, DensityThresholdFactor=100, OutputWorkspace='TOPAZ_35058.nxs_peaks')
#
LoadIsawUB(InputWorkspace='TOPAZ_35058.nxs_peaks', Filename='/SNS/TOPAZ/IPTS-22255/shared/2019Nov/n2_pea_0V/find_peaks/35058_Niggli.mat')
#
CopySample(InputWorkspace='TOPAZ_35058.nxs_peaks', OutputWorkspace='TOPAZ_35058.nxs_event', CopyName=False, CopyMaterial=False, CopyEnvironment=False, CopyShape=False)
IndexPeaks(PeaksWorkspace='TOPAZ_35058.nxs_peaks', Tolerance=0.16, RoundHKLs=False)
ShowPossibleCells(PeaksWorkspace='TOPAZ_35058.nxs_peaks', BestOnly=False, AllowPermutations=False)
IntegrateEllipsoids(InputWorkspace='TOPAZ_35058.nxs_event', PeaksWorkspace='TOPAZ_35058.nxs_peaks', RegionRadius=0.25, SpecifySize=True, PeakSize=0.089999999999999997, BackgroundInnerSize=0.11, BackgroundOuterSize=0.14000000000000001, OutputWorkspace='TOPAZ_35058.nxs_peaks')
OptimizeLatticeForCellType(PeaksWorkspace='TOPAZ_35058.nxs_peaks', CellType='Triclinic', Apply=True, Tolerance=0.16, EdgePixels=19, OutputDirectory='/SNS/TOPAZ/IPTS-22255/shared/2019Nov/CrystalPlan')
LoadIsawUB(InputWorkspace='TOPAZ_35058.nxs_peaks', Filename='/SNS/TOPAZ/IPTS-22255/shared/2019Nov/n2_pea_0V/find_peaks/35058_Niggli.mat')
CopySample(InputWorkspace='TOPAZ_35058.nxs_peaks', OutputWorkspace='TOPAZ_35058.nxs_event', CopyName=False, CopyMaterial=False, CopyEnvironment=False, CopyShape=False)
IndexPeaks(PeaksWorkspace='TOPAZ_35058.nxs_peaks', Tolerance=0.16, RoundHKLs=False)
IndexPeaks(PeaksWorkspace='TOPAZ_35058.nxs_peaks', Tolerance=0.25, RoundHKLs=False)
FindUBUsingIndexedPeaks(PeaksWorkspace='TOPAZ_35058.nxs_peaks', Tolerance=0.16)
IndexPeaks(PeaksWorkspace='TOPAZ_35058.nxs_peaks', Tolerance=0.16, RoundHKLs=False)
CopySample(InputWorkspace='TOPAZ_35058.nxs_peaks', OutputWorkspace='TOPAZ_35058.nxs_event', CopyName=False, CopyMaterial=False, CopyEnvironment=False, CopyShape=False)
IntegrateEllipsoids(InputWorkspace='TOPAZ_35058.nxs_event', PeaksWorkspace='TOPAZ_35058.nxs_peaks', RegionRadius=0.25, SpecifySize=True, PeakSize=0.089999999999999997, BackgroundInnerSize=0.11, BackgroundOuterSize=0.14000000000000001, OutputWorkspace='TOPAZ_35058.nxs_peaks')
proton_charge = mtd['TOPAZ_35058.nxs_event'].getRun().getProtonCharge() * 1000.0  # get proton charge
use_monitor_counts='False'
num_peaks = mtd['TOPAZ_35058.nxs_peaks'].getNumberPeaks()
for i in range(num_peaks):
  peak = mtd['TOPAZ_35058.nxs_peaks'].getPeak(i)
  peak.setMonitorCount( proton_charge )
SaveIsawPeaks(InputWorkspace='TOPAZ_35058.nxs_peaks',Filename='/SNS/TOPAZ/IPTS-22255/shared/2019Nov/n2_pea_0V/find_peaks/35058_Niggli.integrate')
SaveIsawUB(InputWorkspace='TOPAZ_35058.nxs_peaks',Filename='/SNS/TOPAZ/IPTS-22255/shared/2019Nov/n2_pea_0V/find_peaks/35058_Niggli.mat')