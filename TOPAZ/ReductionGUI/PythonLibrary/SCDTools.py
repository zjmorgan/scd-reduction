#"""*WIKI* 
#Library of python functions/methods to facilitate python access for data stored in Mantid Workspace
#*WIKI*"""
import os
import sys
import numpy as np
from numpy.linalg import norm, inv, lstsq
#import ReduceDictionary
from SCDTools import *


#=== Neutron Single Crystal Library ========================
#
# Library of python functions/methods to facilitate python access
# for data stored in Mantid Workspace
#
# 1. Methods to facilitate python access to data stored in Mantid Workspace
#        getUBworkspace(workspace)
#        getUworkspace(workspace)
#        getBworkspace(workspace)
#        getGstarworkspace(workspace)
#        getGworkspace(workspace)
#        getCellworkspace(workspace)
#        getPeakInfo(self,workspace, pkSquenceNo)
#
# 2.  Conversions in orientation matrix
#        readUB(filename)
#        toMantidUB(iswub)
#        toIsawUB(mantidub)
#        getGstarfromUB(UB)
#        getGfromUB(UB)
#        getCellfromG(G)
#        matrixToList(inputUB)    
#
# 3. Conversions between rotation matrix and angles 
#        matrix_from_XYZ_rotations(rotX,rotY,rotZ)
#        angles_from_rotation_matrix(mat)
#        matrix_from_rotation_about_axis(rotangle,rotvector)
#        rotation_angle_about_one_axis(mat)
#
# 4. Correction methods for pseudo rotation of orientation matrix
#        align_matrix(testU, refU, one_axis_rotation = 'True')
#        round_angle(angle)
#        sINT(abc)
#
#
#    Xiaoping Wang, December, 2012
#
class SCDTools():
#class SCDTools(PythonAlgorithm):
#    def category(self):
#        """
#        """
#        return "Crystal;PythonAlgorithms"
#
#    def name(self):
#        """
#        """
#        return "SCDTools"
#
    def __init__(self,*args):
        pass

    # Return UB matrix from an ISAW .mat file
    def readUB(self,filename):
        matfile = open(filename);
        UB = []
        item = 0
        for line in matfile.readlines():
            if item > 2: break
            y = [float(value) for value in line.split()]
            UB.append( y )
            item = item + 1
        matfile.close()
        return np.array(UB)

    # Transform ISAW UB to Mantid UB
    def toMantidUB(self,isawub):
        swapcolumn =np.array([[0.,1.,0.],[0.,0.,1.],[1.,0.,0.]])
        #mantidUB = np.dot(swapcolumn, isawUB.T)
        return np.dot(swapcolumn, isawub.T)
    
    # Transform Mantid UB to ISAW UB
    def toIsawUB(self,mantidub):
         swapcolumn =np.array([[0.,1.,0.],[0.,0.,1.],[1.,0.,0.]])
         #isawUB = np.dot(mantidUB.transpose(), swapcolumn) 
         return np.dot(mantidub.T, swapcolumn) 
     
    # Return G* from UB
    def getGstarfromUB(self,UB):
        Gstar = np.dot(UB.T,UB)
        return Gstar
    
    # Return G from UB
    def getGfromUB(self,UB):
        Gstar = np.dot(UB.T,UB)
        return inv(Gstar)
     
    # Return Cell parameters from G as a list of six float numbers, and unit cell volume
    def getCellfromG(self,G):
        cell_a = np.sqrt(G[0][0])
        cell_b = np.sqrt(G[1][1])
        cell_c = np.sqrt(G[2][2])
        alpha = np.degrees(np.arccos(G[1][2]/cell_b/cell_c))
        beta = np.degrees(np.arccos(G[0][2]/cell_a/cell_c))
        gamma =np.degrees(np.arccos(G[0][1]/cell_a/cell_b))
        V = sqrt(norm(G))
        return ([cell_a, cell_b, cell_c, alpha, beta, gamma], V)

    # Return UB in peak workspace    
    def getUBworkspace(self,workspace):
        sample = workspace.mutableSample()
        lattice = sample.getOrientedLattice()
        return lattice.getUB()
     
    # Get mantid workspace U matrix 
    def getUworkspace(self,workspace):
        sample = workspace.mutableSample()
        lattice = sample.getOrientedLattice()
        return lattice.getU()

    # Get mantid workspace B matrix 
    def getBworkspace(self,workspace):
        sample = workspace.mutableSample()
        lattice = sample.getOrientedLattice()
        return lattice.getB()
    
    # Return mantid workspace Gstar matrix
    def getGstarworkspace(self,workspace):
        sample = workspace.mutableSample()
        lattice = sample.getOrientedLattice()
        return lattice.getGstar()

    # Return mantid workspace G matrix 
    def getGworkspace(self,workspace):
        sample = workspace.mutableSample()
        lattice = sample.getOrientedLattice()
        return lattice.getG()
    
    # Return mantid workspace unit cell parameters
    def getCellworkspace(self,workspace):
        sample = workspace.mutableSample()
        lattice = sample.getOrientedLattice()
        return (lattice.a(), lattice.b(), lattice.c(), lattice.alpha(), lattice.beta(), lattice.gamma())

    # Convert the transMatrixation matrix to 9 comma separated number as a ingle-quoted string
    def matrixToList(self,inputUB):    
        outMatrix_str = inputUB.reshape(9)
        outMatrix_str = ','.join(str(i) for i in outMatrix_str)
        print(outMatrix_str) 
        return outMatrix_str

    # Rotation matrix R = Rz.Ry.Rx represents three rotations, first about the x-axis,
    # then the y-axis and finally the z-axis
    def matrix_from_XYZ_rotations(self,rotX,rotY,rotZ):
        res = np.identity(3, dtype = float)
        psi   = np.radians(rotX)
        theta = np.radians(rotY)
        phi   = np.radians(rotZ)
        cX = np.cos(psi)
        cY = np.cos(theta)
        cZ = np.cos(phi)
        sX = np.sin(psi)
        sY = np.sin(theta)
        sZ = np.sin(phi)
        res = [[cY*cZ,  sX*sY*cZ - cX*sZ,  cX*sY*cZ + sX*sZ],
               [cY*sZ,  sX*sY*sZ + cX*cZ,  cX*sY*sZ - sX*cZ], 
               [ -sY,   sX*cY,             cX*cY           ]]
        return np.array(res)

    # Return Euler angles from rotation first about the x-axis, psi, then the y-axis, theta, and
    # finally the z-axis, phi.
    def angles_from_rotation_matrix(self,mat):
        r00 = mat[0][0]; r01 = mat[0][1]; r02 = mat[0][2];
        r10 = mat[1][0]; r11 = mat[1][1]; r12 = mat[1][2];
        r20 = mat[2][0]; r21 = mat[2][1]; r22 = mat[2][2];
        if (r20 != 1 or r20 != -1) :
            theta = -1*np.arcsin(r20)
            theta1 = np.pi - theta
            psi = np.arctan2(r21/np.cos(theta),r22/np.cos(theta))
            psi1 = np.arctan2(r21/np.cos(theta1),r22/np.cos(theta1))
            phi = np.arctan2(r10/np.cos(theta),r00/np.cos(theta))
            phi1 = np.arctan2(r10/np.cos(theta1),r00/np.cos(theta1))
        else:
            # Gimbal lock: psi and phi are linked at theta equals +- 90 degrees.
            # Set phi to zero and compute psi
            phi = 0.0
            delta = np.arctan2(r01, r02)
            if (r20 == -1):
                theta = np.pi/2.0
                psi = phi + delta
            else:
                theta = -1*np.pi/2.0
                psi = -phi + delta
                # Output two sets of rotation angles
                # psi, theta, phi and psi1, theta1, phi1
        return ([np.degrees(psi), np.degrees(theta), np.degrees(phi), 
             np.degrees(psi1), np.degrees(theta1), np.degrees(phi1)])

    # Find the axis-angle representation of a rotation. 
    # Return a rotation matrix from a rotation about the Euler axis
    def matrix_from_rotation_about_axis(self,rotangle,rotvector):
        ca = np.cos(np.radians(rotangle))
        sa = np.sin(np.radians(rotangle))
        aa = np.zeros([3,3], dtype = float)
        a = np.zeros([3,3], dtype = float)
        r = np.zeros([3,3], dtype = float)
        I = np.identity(3)

        #Apply the Euler-Rodrigues tensor formula in cartesian coordinates
        rotvector = rotvector / norm(rotvector)
        aa[0][0] = rotvector[0]*rotvector[0]
        aa[1][1] = rotvector[1]*rotvector[1]
        aa[2][2] = rotvector[2]*rotvector[2]
        aa[0][1] = rotvector[0]*rotvector[1]; aa[1][0] = aa[0][1]
        aa[0][2] = rotvector[0]*rotvector[2]; aa[2][0] = aa[0][2]
        aa[1][2] = rotvector[1]*rotvector[2]; aa[2][1] = aa[1][2]
        a[0][1] = -rotvector[2]; a[1][0] = -a[0][1]
        a[0][2] =  rotvector[1]; a[2][0] = -a[0][2]
        a[1][2] = -rotvector[0]; a[2][1] = -a[1][2]
        rmat = ca*I + np.dot(I*(1-ca),aa) + sa*a
        #print "Euler - Rodrigues Matrix: \n", rmat.round(4)
        return rmat

    # Return the Euler angle and rotation axis from a rotation matrix
    def rotation_angle_about_one_axis(self,mat):
        epslon = np.float(1.0E-9)     # naught
        # Cosine of rotation  angle
        cRot = 0.5*(mat[0][0] + mat[1][1] + mat[2][2] - 1.0)
        print("cos(theta) = ", cRot)
        
        if abs(cRot - 1) < epslon:
           #CRot = 1
           #Set the rotation angle to zero
           theta = 0
           eRot = [1.0, 0.0, 0.0]
           ieRot = eRot
           ilength = norm(ieRot)
        
        elif abs(cRot + 1) < epslon:
           #CRot = -1
           theta = 180
           length = norm(np.array(mat)+ np.identity(len(mat)))
           eRot = np.array(mat + np.identity(len(mat)))/length
           ieRot = eRot
           ilength = norm(ieRot)
        else:
           theta = np.arccos(cRot)
           e1 = (mat[2][1] - mat[1][2]) / 2/np.sin(theta)
           e2 = (mat[0][2] - mat[2][0]) / 2/np.sin(theta)
           e3 = (mat[1][0] - mat[0][1]) / 2/np.sin(theta)
           length = np.sqrt(e1*e1 + e2*e2 + e3*e3)
           eRot = np.array([e1, e2, e3])/norm(np.array([e1, e2, e3]))
           ieRot = eRot/(np.max(np.absolute(eRot)))*2*(np.sqrt(6))
           s2 = 10.0*np.mod(ieRot,(np.sqrt(2)))
           s3 = 10.0*np.mod(ieRot,(np.sqrt(3)))
           if (np.absolute(s2).all) <= 1:
             ieRot = ieRot /np.sqrt(2)
             ieRot = ieRot/self.sINT(ieRot)
           elif (np.absolute(s3).all) <= 1:
             ieRot = ieRot /np.sqrt(3)
             ieRot = ieRot/self.sINT(ieRot)
           else:
             ieRot = eRot/(np.max(np.absolute(eRot)))*6
             ieRot = ieRot/self.sINT(ieRot)
           ilength = norm(ieRot)
        # 0.  theta is the Euler angle
        # 1.  erot is a unit vector for the rotation axis
        # 2.  ieRot is the axis vector
        # 3. ilength is the length of the axis vector
        return ([np.degrees(theta), eRot, ieRot, ilength])

    # Return the smallest whole number in an array
    def sINT(self,abc):
        small = 1
        abc = np.rint(np.absolute(abc))
        abc.sort()
        for item in abc:
          if (item != 0):
              small = item
              break
        return small


    #Round a rotation angle
    def round_angle(self,angle):
        angle = float(angle)
        outangle = angle
        print('Rotation angle from input = ', angle)
        if abs(angle -180.0) < 9.0:
          outangle = 180.0
        if abs(angle + 180) < 9.0:
          outangle = -180.0
        if abs(angle -150.0) < 2.5:
          outangle = 150.0
        if abs(angle + 150) < 2.5:
          outangle = -150.0
        if abs(angle - 144.73) < 2.63:
          outangle = 144.73333
        if abs(angle + 144.73) < 2.63:
          outangle = -144.73333
        if abs(angle - 135) < 3.0:
          outangle = 135
        if abs(angle + 135) < 3.0:
          outangle = -135
        if abs(angle - 125.27) < 2.63:
          outangle = 125.266667
        if abs(angle + 125.27) < 2.63:
          outangle = -125.266667
        if abs(angle - 120) < 2.63:
          outangle = 120.0
        if abs(angle + 120) < 2.63:
          outangle = -120.0
        if abs(angle - 109.47) < 3:
          outangle = 109.466667
        if abs(angle + 109.47) < 3:
          outangle = -109.466667
        if abs(angle - 90) < 5.0:
          outangle = 90.0
        if abs(angle + 90) < 5.0:
          outangle = -90.0
        if abs(angle - 75) < 2.2:
          outangle = 75.0
        if abs(angle + 75) < 2.2:
          outangle = -75.0
        if abs(angle - 70.53) < 2.2:
          outangle = 70.533333
        if abs(angle + 70.53) < 2.2:
          outangle = -70.533333
        if abs(angle - 60) < 2.63:
          outangle = 60.0
        if abs(angle + 60) < 2.63:
          outangle = -60.0
        if abs(angle - 54.73) < 2.63:
          outangle = 54.733333
        if abs(angle + 54.73) < 2.63:
          outangle = -54.733333
        if abs(angle - 45) < 5.0:
          outangle = 45
        if abs(angle + 45) < 5.0:
          outangle = -45
        if abs(angle - 35.27) < 2.63:
          outangle = 35.266667
        if abs(angle + 35.27) < 2.63:
          outangle = -35.266667
        if abs(angle - 30) < 2.63:
          outangle = 30.0
        if abs(angle + 30) < 2.63:
          outangle = -30.0
        if abs(angle) < 12.0:
          outangle = 0.0 
        print('Rotation angle changed to = ', outangle)
        return outangle

    def adjust_raxis(self, vnorm):
        print('Vector from imput: ', vnorm)
        vnorm = vnorm/norm(vnorm)
        signv =[x/np.absolute(x) for x in vnorm]
        r_vnorm =np.zeros(3)
        for i,vi in enumerate(vnorm):
          if np.absolute(vi) < 0.10:
            r_vnorm[i] = 0.0
          else:
            print('1/vi**2 for '+str(i)+': ', vi, 3.0/ (vi*vi),24.0/(vi*vi))
            rsqr = 24.0/(vi*vi)
            if np.absolute(rsqr - 28) < 2:
              rsqr = 28.0
            if np.absolute(rsqr - 36) < 2:
              rsqr = 36.0
            if np.absolute(rsqr - 42) < 2:
              rsqr = 42.0
            if np.absolute(rsqr - 56) < 2:
              rsqr = 56.0
            if np.absolute(rsqr - 63)<4 and (np.absolute(vnorm).max() - np.absolute(vi))>=0.08:
              rsqr = 63
            if np.absolute(rsqr - 64)<4 and (np.absolute(vnorm).max() - np.absolute(vi))<0.08:
              rsqr = 64
            if np.absolute(rsqr - 72) < 4:
              rsqr = 72.0
            if np.absolute(rsqr - 84) < 4:
              rsqr = 84.0
            if np.absolute(rsqr - 144) < 16:
              rsqr = 144.0
            if np.absolute(rsqr - 336) < 16:
              rsqr = 336
            if np.absolute(rsqr - 504) < 35:
              rsqr = 504
            r_vnorm[i] = signv[i] * np.sqrt(24.0/rsqr)
          print(r_vnorm[i])
          print('adjusted V'+str(i)+': ', r_vnorm[i])
        print('Returned vector: ', r_vnorm/norm(r_vnorm))
        return np.array(r_vnorm/norm(r_vnorm))


    # Align the orientation matrices for runs at different setting angles
    # Return the orientation matrix corrected for pseudo-rotation
    def align_matrix(self,testU, refU, one_axis_rotation = 'True'):
        #Find the rotation matrix between two orientation matrices
        rotmat = np.dot(testU.T,refU)
        print("Difference in orientation U from that of the reference frame \n", rotmat)
        print("trace =", np.trace(rotmat))
        if one_axis_rotation != 'True':
          #Method using rotations about X, Y and Z axes
          rotangle = np.round(self.angles_from_rotation_matrix(rotmat),1)
          print(rotangle)
          xRot = self.round_angle(float(rotangle[0]))
          yRot = self.round_angle(float(rotangle[1]))
          zRot = self.round_angle(float(rotangle[2]))
          #Rotate testU to match the orientation of the refU
          alignedU = np.dot(testU, self.matrix_from_XYZ_rotations(xRot, yRot, zRot) )
        else:
          #Method using rotation about one axis
          angle_axis = self.rotation_angle_about_one_axis(rotmat)
          angle = float((angle_axis[0]).round(4))
          print("Euler angle = ",  angle)
          print("Euler angle rounded to 360/n degrees: ", angle)
          angle= self.round_angle(angle)
          print("Unitary rotation vector \n", angle_axis[1])
          rotation_axis = self.adjust_raxis(angle_axis[1])
          print("Rotation axis \n", rotation_axis)
          alignedU = np.dot(testU, self.matrix_from_rotation_about_axis(angle, rotation_axis))
          print('alignedU \n', alignedU)
          #print "Rotation axis \n", angle_axis[2].round(4)
          #TODO Test the effect of rounding the axis vector to integer numbers
          #print "Indices of rotation axis rounded to integer numbers"
          #print  angle_axis[2].round(0)
          #Rotate testU to match the orientation of the refU
          #alignedU = np.dot(testU, self.matrix_from_rotation_about_axis(angle, angle_axis[2].round(0)))
        return  alignedU

    #Retrieve peak parameters from Mantid PeakWorkspace
    def getPeakInfo(self,workspace, pkSquenceNo):
        onePeak = workspace.getPeak(pkSquenceNo)
        hkl = onePeak.getHKL()
        h = hkl.X()
        k = hkl.Y()
        l = hkl.Z()
        pkInt = np.float(onePeak.getIntensity())
        pkSig = np.float(onePeak.getSigmaIntensity())
        wl = onePeak.getWavelength()
        twoth = onePeak.getScattering()
        L2 = onePeak.getL2()
        dsp =np.float(onePeak.getDSpacing())
        qsample = onePeak.getQSampleFrame()
        qlab = onePeak.getQLabFrame()
        onepkQSample =[qsample.X(),qsample.Y(),qsample.Z()]
        onepkQLab =[qlab.X(),qlab.Y(),qlab.Z()]
        detCol =onePeak.getCol()
        detRow =onePeak.getRow()
        tof = onePeak.getTOF()
        detPos =onePeak.getDetPos()
        detX = detPos.X()
        detY = detPos.Y()
        detZ = detPos.Z()
        detID =onePeak.getDetectorID()
        run = onePeak.getRunNumber()
        #Output parameters
        #0  [h,k,l]      
        #1  pkInt        
        #2  pkSig        
        #3  wl           
        #4  L2           
        #5  twoth        
        #6  dsp          
        #7  onepkQSample 
        #8  onepkQLab    
        #9  detCol       
        #10 detRow       
        #11 [detX, detY, detZ]   
        #12 run       
        #13 detID   
        return ([h,k,l], pkInt, pkSig, wl, L2, twoth, dsp,
                onepkQSample, onepkQLab, detCol, detRow,[detX, detY, detZ], run, detID)
   
#mtd.registerPyAlgorithm(SCDTools())
 
