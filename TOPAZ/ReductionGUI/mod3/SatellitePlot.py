import os
import sys
sys.path.insert(0,"/opt/mantidnightly/bin")
sys.path.insert(0,"/opt/mantidnightly/lib")

#sys.path.insert(0,"/SNS/users/vel/mantid/release/bin")
from mantid.simpleapi import *
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import LogFormatter
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import proj3d
import numpy as np
import math
from mantid.kernel import ConfigService
from mantid.kernel import V3D
from mantid.geometry import OrientedLattice
from mantid import config
#ConfigService.setLogLevel(2)
config['Q.convention'] = 'Crystallography'
from scipy.spatial import cKDTree


class SatellitePlot():
    def __init__(self):
        self.tolerance = 0.5
        self.minIntensity = 10
        self.mod = ""

    def orthogonal_proj(self, zfront, zback):
        a = (zfront+zback)/(zfront-zback)
        b = -2*(zfront*zback)/(zfront-zback)
        return np.array([[1,0,0,0],
                            [0,1,0,0],
                            [0,0,a,b],
                            [0,0,-0.0001,zback]])


    def load_peaks(self, peaksFile, UBFile):
            #End Input from GUI
            self.peaks_ws = LoadIsawPeaks( Filename = peaksFile, OutputWorkspace='peaks')
            LoadIsawUB(self.peaks_ws, UBFile)

    def plot_Qpeaks(self):
            proj3d.persp_transformation = self.orthogonal_proj
            try:
            # Find size of screen
                import curses
                stdscr = curses.initscr()
                self.screen_x,self.screen_y = stdscr.getmaxyx()
                self.screen_x = min(self.screen_x, self.screen_y)
                self.screen_y = self.screen_x
            except:
                self.screen_x = 40
                self.screen_y = 40
            plt.rcParams.update({'font.size': 10})
            if len(self.mod) == 0:
                IndexPeaks( PeaksWorkspace = self.peaks_ws, Tolerance = self.tolerance, RoundHKLs = False, CommonUBForAll=False)
            elif len(self.mod) > 0:
                #IndexPeaksWithSatellites(PeaksWorkspace=self.peaks_ws, RoundHKLs=False, ModVector1=self.mod, MaxOrder=1)
                IndexPeaks(PeaksWorkspace=self.peaks_ws, 
                  RoundHKLs=False, 
                  ModVector1=self.mod, 
                  MaxOrder = 2, CommonUBForAll =False)

            npeaksTotal = self.peaks_ws.getNumberPeaks()
            self.c = []
            self.x = []
            self.y = []
            self.z = []
            self.s = []
            for i in range(npeaksTotal):
                    peak = self.peaks_ws.getPeak(i)
                    hklp = V3D(peak.getH()-math.floor(peak.getH()),peak.getK()-math.floor(peak.getK()),peak.getL()-math.floor(peak.getL()))
                    if peak.getSigmaIntensity() !=  0.0:
                        IsigI = peak.getIntensity()/peak.getSigmaIntensity()
                    else:
                        IsigI = 0.0
                    #if IsigI < 50:
                        #continue
                    intensity = max(self.minIntensity, peak.getIntensity())
                    for i in range(-1,1):
                        for j in range(-1,1):
                            for k in range(-1,1):
                                if (abs(hklp[0]+i) < 1.0 and abs(hklp[1]+j) < 1.0 and abs(hklp[2]+k) < 1.0):
                                    self.x.append(hklp[0]+i)
                                    self.y.append(hklp[1]+j)
                                    self.z.append(hklp[2]+k)
                                    self.c.append(intensity)
                                    self.s.append(5)
            nPts = 100
            x_grid = np.linspace(-1, 1, nPts+1)
            x_half = np.zeros(nPts)
            for j in range(nPts):
                x_half[j] = 0.5*(x_grid[j]+x_grid[j+1])
            sumIntensityX = np.zeros(nPts)
            for i in range(len(self.x)):
                for j in range(nPts):
                  if self.x[i] > x_grid[j] and self.x[i] < x_grid[j+1]:
                      sumIntensityX[j] = sumIntensityX[j] + self.c[i]
            sumIntensityY = np.zeros(nPts)
            for i in range(len(self.y)):
                for j in range(nPts):
                  if self.y[i] > x_grid[j] and self.y[i] < x_grid[j+1]:
                      sumIntensityY[j] = sumIntensityY[j] + self.c[i]
            sumIntensityZ = np.zeros(nPts)
            for i in range(len(self.z)):
                for j in range(nPts):
                  if self.z[i] > x_grid[j] and self.z[i] < x_grid[j+1]:
                      sumIntensityZ[j] = sumIntensityZ[j] + self.c[i]
            # Plot figure with subplots of different sizes
            fig = plt.figure("Modulated Structures",figsize = (self.screen_x, self.screen_y))
            # set up subplot grid
            gridspec.GridSpec(3,3)
            ax = plt.subplot2grid((3,3), (0,2))
            plt.plot(x_half,sumIntensityX)
            maxGrid = "max\n"
            for i in range(2,nPts-2):
                if abs(x_half[i]) > 0.55 or abs(x_half[i]) < 0.05:
                    continue
                if sumIntensityX[i] > sumIntensityX[i-1] and sumIntensityX[i] > sumIntensityX[i+1] and sumIntensityX[i-1] > sumIntensityX[i-2] and sumIntensityX[i+1] > sumIntensityX[i+2]:
                    maxGrid += "{:f}\n".format(x_half[i])
            ax.text(0.05, 0.60, maxGrid,
                fontsize = 10, transform = ax.transAxes)
            ax.set_yscale('log')
            ax.set_xlabel('H')
            ax.set_ylabel(r'$\Sigma$I')
            ax = plt.subplot2grid((3,3), (1,2))
            plt.plot(x_half,sumIntensityY)
            maxGrid = "max\n"
            for i in range(2,nPts-2):
                if abs(x_half[i]) > 0.5 or abs(x_half[i]) < 0.05:
                    continue
                if sumIntensityY[i] > sumIntensityY[i-1] and sumIntensityY[i] > sumIntensityY[i+1] and sumIntensityY[i-1] > sumIntensityY[i-2] and sumIntensityY[i+1] > sumIntensityY[i+2]:
                    maxGrid += "{:f}\n".format(x_half[i])
            ax.text(0.05, 0.60, maxGrid,
                fontsize = 10, transform = ax.transAxes)
            ax.set_yscale('log')
            ax.set_xlabel('K')
            ax.set_ylabel(r'$\Sigma$I')
            ax = plt.subplot2grid((3,3), (2,2))
            plt.plot(x_half,sumIntensityZ)
            maxGrid = "max\n"
            for i in range(2,nPts-2):
                if abs(x_half[i]) > 0.5 or abs(x_half[i]) < 0.05:
                    continue
                if sumIntensityZ[i] > sumIntensityZ[i-1] and sumIntensityZ[i] > sumIntensityZ[i+1] and sumIntensityZ[i-1] > sumIntensityZ[i-2] and sumIntensityZ[i+1] > sumIntensityZ[i+2]:
                    maxGrid += "{:f}\n".format(x_half[i])
            ax.text(0.05, 0.60, maxGrid,
                fontsize = 10, transform = ax.transAxes)
            ax.set_yscale('log')
            ax.set_xlabel('L')
            ax.set_ylabel(r'$\Sigma$I')
          
            vmin = min(self.c)
            vmax = max(self.c)
            logNorm = colors.LogNorm(vmin = vmin, vmax = vmax)
            self.axP = plt.subplot2grid((3,3), (0,0), colspan=2, rowspan=3, projection='3d')
            sp = self.axP.scatter(self.x, self.y, self.z, c = self.c, vmin = vmin, vmax = vmax, norm = logNorm, s = self.s, cmap='rainbow', picker = True, alpha=0.2)
               

            self.props = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 1.0)


            #cid = fig.canvas.mpl_connect('pick_event', self.onpick3)
            pid = fig.canvas.mpl_connect('key_press_event', self.onpress)

            try:
                lattice = self.peaks_ws.sample().getOrientedLattice()
            except:
                lattice = OrientedLattice(1,1,1)
            astar = V3D(0.5,0,0)
            bstar = V3D(0,0.5,0)
            cstar = V3D(0,0,0.5)
            hkls = np.array([np.array([x, y, z]) for x,y,z in zip(self.x,self.y,self.z)])
            labels, centroids = self.dbscan(hkls, eps=.0425, min_points=6)
            
            i = 0
            modTest = str(self.peaks_ws.getPeak(0).getRunNumber())
            lastRun = str(self.peaks_ws.getPeak(npeaksTotal-1).getRunNumber())
            if lastRun != modTest:
                modTest = modTest + "-" +lastRun
            modTest = modTest + " Type a, b, or c for axis alignment"
            for x in centroids:
                self.axP.plot([x[0]], [x[1]], [x[2]], 'wo', ms=40, markerfacecolor="None", markeredgecolor='red', markeredgewidth=3)
                if x[0]*x[0]+x[1]*x[1]+x[2]*x[2] > 0.01 and x[0]>=-0.5 and x[1]>=-0.5 and x[2]>=-0.5 and x[0] <=0.5 and x[1] <= 0.5 and x[2] <=0.5:
                    modTest += "\nModulation Vector = "+"{:.3f}".format(x[0]) + ", " + "{:.3f}".format(x[1]) +", " + "{:.3f}".format(x[2])
                    if len(self.mod) == 0:
                        self.mod = "{:.3f}".format(x[0]) + "," + "{:.3f}".format(x[1]) +"," + "{:.3f}".format(x[2])
                        #IndexPeaksWithSatellites(PeaksWorkspace=self.peaks_ws, RoundHKLs=False, ModVector1=self.mod, MaxOrder=1)
                        IndexPeaks(PeaksWorkspace=self.peaks_ws, 
                                   RoundHKLs=False, ModVector1=self.mod, MaxOrder = 1, CommonUBForAll =False)
            
            self.axP.plot([0,astar[0]],[0,astar[1]],[0,astar[2]], color='r')
            self.axP.text(astar[0],astar[1],astar[2],  '%s' % ('a*'), size=20, zorder=1,  color='r') 
            self.axP.plot([0,bstar[0]],[0,bstar[1]],[0,bstar[2]], color='b')
            self.axP.text(bstar[0],bstar[1],bstar[2],  '%s' % ('b*'), size=20, zorder=1,  color='b') 
            self.axP.plot([0,cstar[0]],[0,cstar[1]],[0,cstar[2]], color='g')
            self.axP.text(cstar[0],cstar[1],cstar[2],  '%s' % ('c*'), size=20, zorder=1,  color='g') 

            self.axP.text2D(0.05, 0.80, modTest +
                "\nLattice = " + " " + "{:.3f}".format(lattice.a()) + " " + "{:.3f}".format(lattice.b()) + " " + "{:.3f}".format(lattice.c()) + " " +
                "{:.3f}".format(lattice.alpha()) + " " + "{:.3f}".format(lattice.beta()) + " " + "{:.3f}".format(lattice.gamma()) +
                "\nError = " + " " + "{:.3E}".format(lattice.errora()) + " " + "{:.3E}".format(lattice.errorb()) + " " + "{:.3E}".format(lattice.errorc()) + " " +
                "{:.3E}".format(lattice.erroralpha()) + " " + "{:.3E}".format(lattice.errorbeta()) + " " + "{:.3E}".format(lattice.errorgamma()) ,
                fontsize = 10, transform = self.axP.transAxes)
            self.axP.set_xlabel('H')
            self.axP.set_ylabel('K')
            self.axP.set_zlabel('L')
    
    
            # Create cubic bounding box to simulate equal aspect ratio
            max_range = max([max(self.x)-min(self.x), max(self.y)-min(self.y), max(self.z)-min(self.z)])
            Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(max(self.x)+min(self.x))
            Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(max(self.y)+min(self.y))
            Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(max(self.z)+min(self.z))
            # Comment or uncomment following both lines to test the fake bounding box:
            for xb, yb, zb in zip(Xb, Yb, Zb):
               self.axP.plot([xb], [yb], [zb], 'w')
            plt.show()
            #fig.canvas.mpl_disconnect(cid)
            fig.canvas.mpl_disconnect(pid)

    def dbscan(self, points, eps=1.0, min_points=5):
        """DBSCAN Clustering algorithm
        
        @param points
        """
        tree = cKDTree(points)
        neighbours = tree.query_ball_point(points, eps)
        #neighbours = tree.query_ball_point(points, eps, n_jobs=-1)
        neighbours_count = np.array([len(x) for x in neighbours])
    
        labels = np.zeros(points.shape[0])
        cluster_index = 1
    
        indices = np.arange(points.shape[0])
        min_pts_mask = (neighbours_count > min_points)
        while True:
            possible_seeds = indices[min_pts_mask & (labels == 0)]
    
            if len(possible_seeds) == 0:
                break
    
            seed_index = np.random.choice(possible_seeds)
            search_queue = set([seed_index])
    
            # breadth first search of each neighbour of this cluster
            while len(search_queue) > 0:
                current_index = search_queue.pop()
    
                if labels[current_index] > 0:
                    continue
    
                if neighbours_count[current_index] < min_points:
                    # Ok, it's a leaf
                    labels[current_index] = cluster_index
                else:
                    # It's a branch
                    labels[current_index] = cluster_index
                    search_queue.update(neighbours[current_index])
    
            cluster_index += 1
    
        # Now we're done clustering find the centroids of the clusters
    
        def compute_centroid(cluster_index):
            subset = points[labels == cluster_index]
            return np.mean(subset, axis=0)
    
        n_clusters = cluster_index
        centroids = list(map(compute_centroid, list(range(1, n_clusters))))
    
        return labels, np.array(list(centroids))

    def onpress(self, event):
        if event.key not in ('a', 'b', 'c'):
            return
        if event.key == 'a':
            self.axP.view_init(elev=0, azim=0)
        elif event.key == 'b':
            self.axP.view_init(elev=0, azim=-90)
        elif event.key == 'c':
            self.axP.view_init(elev=90, azim=-90)
        plt.show()
                
    def onpick3(self, event):
        ind = event.ind
        zdirs = None
        xp = np.take(self.x, ind)
        yp = np.take(self.y, ind)
        zp = np.take(self.z, ind)
        label = 'h,k,l:'+'%.3f %.3f %.3f'%(xp[0],yp[0],zp[0])
        global txt
        try:
            txt.remove()
        except:
            print ("First peak picked")
        txt = self.axP.text(xp[0],yp[0],zp[0],label,zdirs, bbox = self.props)

#===============================================================================================

if __name__ == '__main__':  # if we're running file directly and not importing it
    print ("Wait for peaks in HKL")
    test = SatellitePlot()  # run the main function
    test.load_peaks('/SNS/TOPAZ/IPTS-22357/shared/2020/FeCuTe2_3x3x6_sat_find_peaks/FeCuTe2_3x3x6f_sat_Niggli.integrate', 
                    '/SNS/TOPAZ/IPTS-22357/shared/2020/FeCuTe2_3x3x6_sat_find_peaks/FeCuTe2_3x3x6f_sat_Niggli.mat')

#    test.load_peaks('/SNS/TOPAZ/IPTS-25887/shared/CrystalPlan/5K_36936_Monoclinic_I_sat.integrate',
#                     '/SNS/TOPAZ/IPTS-25887/shared/CrystalPlan/5K_36936_Monoclinic_I_sat.mat')
#    test.load_peaks( "/SNS/users/vel/shared/MOD/MOD_Triclinic_P.integrate","/SNS/users/vel/shared/MOD/MOD_Triclinic_P.mat")
#    test.load_peaks( "/SNS/users/vel/CORELLI_59429.peaks","/SNS/users/vel/CORELLI_59429.mat")
#    test.load_peaks( "/SNS/TOPAZ/IPTS-22342/shared/PNO_PSI_Mag1_all_bg/PNO_PSI_Mag1_all_bg_Orthorhombic_F.integrate", "/SNS/TOPAZ/IPTS-22342/shared/PNO_PSI_Mag1_all_bg/PNO_PSI_Mag1_all_bg_Orthorhombic_F.mat")
#    test.load_peaks('/SNS/TOPAZ/IPTS-10003/shared/integration/8071_P.integrate',"/SNS/TOPAZ/IPTS-10003/shared/integration/8071_P.mat")
#    test.load_peaks('/SNS/TOPAZ/IPTS-18749/shared/7147A_load_peaks/24281_Niggli.integrate',"/SNS/TOPAZ/IPTS-18749/shared/7147A_load_peaks/24281_Niggli.mat")
    test.plot_Qpeaks()
    print ("Plotting finished for current data")






