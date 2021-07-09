#--------------------------------------#
# Cystal function library.             #
# A. J. Schultz     from October, 2006 #
#--------------------------------------#

"""
List of functions (June 5, 2012):

abc(UB)
    Input UB matrix.
    Return unit cell parameters.

calc_2th_wl_IPNS(q)
    Input q vector in IPNS coordinate system in units of 1/d.
    Return two-theta (deg) and wavelength.

calc_2th_wl_SNS(q)
    Input q vector in SNS coordinate system in units of 1/d.
    Return two-theta (deg) and wavelength.
    
center(h, k, l, center_type)
    Input hkl and the lattice centering type.
    Return True if hkl is allowed.
    Return False if hkl violates centering and is not allowed.

det_coord(q, wl, deta1, deta2, detd, det_rot_ang)
    *** This function is outdated.
    Input q vector, wl, and various detector parameters.
    Return detector coordinates xcm, ycm.
    
Gmatrix(a, b, c, alpha, beta, gamma)
    Input unit cell parameters in Angstroms and angles.
    Return real metric tensor G.
    
huq(h, k, l, UB)
    Input hkl and UB matrix.
    Return q vector in units of 1/d.

polar_paramters (q)
    Input q vector.
    Return polar coordinates for polar plot.

qvec(twoth, az, wl)
    Input the two-theta and azimuthal angles (rad) and the wavelength.
    Return the q vector.
    
rotate_matrix(UB, omega, chi, phi, SNS_or_IPNS)
    Input the UB matrix and the omega, chi and phi angles (deg).
    Return a rotated orientation matrix.
    
rotate_vector(q, omega, chi, phi, SNS_or_IPNS)
    Input q vector and 3 rotation angles (deg).
    Return a rotated q vector.
    
rotation_matrix(omega, chi, phi, SNS_or_IPNS):
    Input rotation angles.
    Return the rotation matrix.
    
UB_IPNS_2_SNS(UB_IPNS)
    Input the IPNS UB matrix.
    Return the SNS UB matrix.
"""

import math
import numpy
import numpy.linalg as linalg

#--------------------------------------------------------
#               function abc
#--------------------------------------------------------
# Function to obtain unit cell parameters from UB matrix.
# Same as the Fortran SUBROUTINE ABC.
#--------------------------------------------------------
def abc(UB):
    """Returns unit cell parameters calculated from UB matrix."""
    
    UBt = numpy.transpose(UB) #UBt the the transpose of UB

    Gi = numpy.dot(UB,UBt)    #Gi is the reciprocal metric tensor
                        #Gi is the product of UB and the transpose of UB

    G = linalg.inv(Gi)         #G is the real metric tensor

    # Obtain unit cell parameters
    rad = 180./math.pi

    a = numpy.zeros((7))       
    a[0] = math.sqrt(G[0,0])                     #a
    a[1] = math.sqrt(G[1,1])                     #b
    a[2] = math.sqrt(G[2,2])                     #c
    a[3] = math.acos(G[1,2] / (a[1] * a[2])) * rad     #alpha in degrees
    a[4] = math.acos(G[0,2] / (a[0] * a[2])) * rad     #beta in degrees
    a[5] = math.acos(G[0,1] / (a[0] * a[1])) * rad     #gamma in degrees
    a[6] = math.sqrt(linalg.det(G))                     #Volume

    return a
    

#---------------------------------------------
#           function calc_2th_wl_IPNS
#---------------------------------------------
# Calculate two-theta and wavelength from
# q vector in the IPNS coordinate system
# with x along the beam.
#---------------------------------------------
def calc_2th_wl_IPNS(q):
    "Returns two-theta angle and wavelength for q-vector"

    rad = 180. / math.pi

    dstar = math.sqrt(q[0]**2 + q[1]**2 + q[2]**2)
    d = 1.0 / dstar

    b = numpy.zeros((2))

    theta = 90.0 - math.acos(-q[0] / dstar) * rad    # theta in degrees
    b[0] = 2.0 * theta    # two-theta
	
    b[1] = 2.0 * d * math.sin((theta) * math.pi / 180.0) # wavelength
    b[1] = abs(b[1])

    return b
    
#---------------------------------------------
#           function calc_2th_wl_SNS
#---------------------------------------------
# Calculate two-theta and wavelength from
# q vector in the SNS coordinate system
# with z along the beam.
#---------------------------------------------  
def calc_2th_wl_SNS(q):
    "Returns two-theta angle and wavelength for q-vector"

    rad = 180. / math.pi

    dstar = math.sqrt(q[0]**2 + q[1]**2 + q[2]**2)
    d = 1.0 / dstar

    b = numpy.zeros((2))

    theta = 90.0 - (math.acos(-q[2] / dstar) * rad)    # theta
    b[0] = 2.0 * theta    # two-theta
	
    b[1] = 2.0 * d * math.sin((theta) * math.pi / 180.0) # wavelength
    b[1] = abs(b[1])

    return b


#--------------------------------------------------------
#               function center
#--------------------------------------------------------
#  Test for centering.
#  Return True if peak is allowed.
#  Return False if peaks is not allowed.
#--------------------------------------------------------
def center(h, k, l, center_type):
    
    if center_type == 'P':
        return True
    
    if center_type == 'A':
        sum = k + l
        if (sum % 2) == 0: return True
        return False
        
    if center_type == 'B':
        sum = h + l
        if (sum % 2) == 0: return True
        return False
        
    if center_type == 'C':
        sum = h + k
        if (sum % 2) == 0: return True
        return False
        
    if center_type == 'F':
        sum = h + k
        if (sum % 2) != 0: return False
        sum = h + l
        if (sum % 2) != 0: return False
        sum = k + l
        if (sum % 2) != 0: return False
        return True
        
    if center_type == 'I':
        sum = h + k + l
        if (sum % 2) != 0: return False
        return True
        
    if center_type == 'R':
        sum = -h + k + l
        if (sum % 3) != 0: return False
        return True
        
    print('Centering type not P, A, B, C, F, I or R.')

    
#---------------------------------------------
#           function det_coord
#---------------------------------------------
# Calculate detector positions (xcm,ycm)
# for the IPNS coordinate system.
#---------------------------------------------
def det_coord(q, wl, deta1, deta2, detd, det_rot_ang):
    "Returns detector coordinates xcm,ycm"

    # The detector will be assumed to be centered at zero angles,
    # with coordinates in q-space of (1,0,0). So the procedure is to rotate
    # the q-vector for the hkl peak by the two detector angles.

    # In ISAW, the coordinate system is x parallel
    # to the beam, with positive x pointing downstream.
    # Looking down the x-axis, +y is to the right and
    # +z points up. This means that the IPNS SCD detectors
    # are located at -120 and -75 degrees.

    #   beam stop <... crystal ...> source
    #            +x <----X
    #                    |
    #                    |
    #                    v +y

    # First translate the origin from the reciprocal
    # lattice origin to the center of the sphere of
    # reflection.
    xdp = q[0] + (1.0 / wl)

    # Angles and rotations as defined here:
    # http://mathworld.wolfram.com/SphericalCoordinates.html
    # and here:
    # http://mathworld.wolfram.com/RotationMatrix.html

    # Rotate to a deta1 to zero

    ang = radians(-deta1)               #Change sign of deta, 8/30/07
    if ang < 0.0: ang = ang + 2.0 * math.pi    #This ensures a counterclockwise rotation around x
    xt = xdp
    yt = q[1]
    xdp = xt * math.cos(ang) - yt * math.sin(ang)
    ydp = xt * math.sin(ang) + yt * math.cos(ang)

    ang = math.radians(-deta2)               #Again counterclockwise rotation around y
    xt = xdp
    zt = q[2]
    xdp = xt * math.cos(ang) - zt * math.sin(ang)
    zdp = xt * math.sin(ang) + zt * math.cos(ang)    

    # Calculate xcm and ycm

    xcm0 = -(ydp / xdp) * detd
    ycm0 = (zdp / xdp) * detd

    #  If detector is rotated (usually by 45 deg.),
    #  calculate detector coordinates.
    ang = math.radians(det_rot_ang)
    xcm = xcm0 * math.cos(ang) - ycm0 * math.sin(ang)
    ycm = xcm0 * math.sin(ang) + ycm0 * math.cos(ang)
    
    
    # If xdp is negative, then diffracted ray is in the opposite direction
    # of the detector.

    if xdp <= 0.0:
        xcm = ycm = 99.

    dc = numpy.zeros((2))
    dc = [xcm, ycm]

    return dc

#---------------------------------------------
#           function Gmatrix
#---------------------------------------------
# Calculate the metric tensor G
# from the lattice constants.
#---------------------------------------------

def Gmatrix(a, b, c, alpha, beta, gamma):
    """Create the G real space metric tensor."""
    G = numpy.zeros((3, 3))     # real metric tensor
    
    alpha = math.radians(alpha)
    beta = math.radians(beta)
    gamma = math.radians(gamma)

    G[0][0] = a * a
    G[0][1] = a * b * math.cos(gamma)
    if math.abs(G[0][1]) < 1.0e-10: G[0][1] = 0.0
    G[0][2] = a * c * math.cos(beta)
    if math.abs(G[0][2]) < 1.0e-10: G[0][2] = 0.0
    G[1][0] = G[0][1]
    G[1][1] = b * b
    G[1][2] = b * c * math.cos(alpha)
    if math.abs(G[1][2]) < 1.0e-10: G[1][2] = 0.0
    G[2][0] = G[0][2]
    G[2][1] = G[1][2]
    G[2][2] = c * c
    
    return G
 
 
#---------------------------------------------
#           function huq
#---------------------------------------------
#   Multiply hkl times UB matrix to obtain
#   qx,qy,qz.
#---------------------------------------------
def huq(h, k, l, UB):
    "Multiplies hkl times UB matrix to return q-vector"

    hh = [h, k, l]

    q = numpy.zeros((3))
    q = numpy.dot(hh, UB)

    return q

#---------------------------------------------
#       function polar_parameters
#---------------------------------------------
# Calculate coordinates for polar plots.
#---------------------------------------------
def polar_paramters (q):
    "Return coordinates for polar plot."


    pc = numpy.zeros((4))

    if q[0]!=0.0:
        pc[0] = math.degrees(math.atan(q[1] / q[0])) # polar plot theta angle
        if q[0]<0: pc[0] = pc[0] + 180.
    else:
        pc[0] = 90.0
        if q[1]<0: pc[0]=270.0
    
    #  Spherical angle phi, which will be the polar coordinate R
    pc[1] = math.degrees(acos(q[2]/sqrt(q[0]**2 + q[1]**2 + q[2]**2)))
    
    #  Cartesian coordinates xp,yp for plotting with some programs
    pc[2] = pc[1] * math.cos(math.radians(pc[0]))
    pc[3] = pc[1] * math.sin(math.radians(pc[0]))
    
    return pc

#---------------------------------------------
#       function qvec
#---------------------------------------------
# Calculate the q vector from two-theta, az 
# and wl for a peak.
#---------------------------------------------
def qvec(twoth, az, wl):
    "Return q vector for a peak in peaks file."
    
    q = numpy.zeros((3))
    
    # IPNS axes with x is the beam direction and z is vertically upward
    # q[0] = cos(az)*sin(twoth)/wl
    # q[1] = sin(az)*sin(twoth)/wl
    # q[2] = (cos(twoth) - 1.0)/wl
    
    # SNS axes with z is the beam direction and y is vertically upward
    q[0] = (math.cos(twoth) - 1.0) / wl
    q[1] = math.cos(az) * math.sin(twoth) / wl
    q[2] = math.sin(az) * math.sin(twoth) / wl
    
    return q
    
    
#--------------------------------------------------------
#               function rotate_matrix
#--------------------------------------------------------
# Same as Fortran SUBROUTINE NEWROT.
#C
#C   ROTATES A 3X3 MATRIX FOR WHICH ALL ANGLES ARE ZERO TO
#C   A 3X3 MATRIX FOR WHICH ALL ANGLES ARE NON ZERO. SEE W.
#C   C. HAMILTON, INT. TABLES. IV, PP 275-281
#--------------------------------------------------------
def rotate_matrix(UB, omega, chi, phi, SNS_or_IPNS):
    "Rotates UB matrix by setting angles"

    fmat = rotation_matrix(omega, chi, phi, SNS_or_IPNS)

    newmat = numpy.zeros((3,3))
    newmat = numpy.dot(UB, fmat)

    return newmat

#---------------------------------------------
#       function rotate_vector
#---------------------------------------------
# Rotate q vector by angles omega, chi and phi
#---------------------------------------------    
def rotate_vector(q, omega, chi, phi, SNS_or_IPNS):
    "Rotates q-vector by setting angles"

    fmat = rotation_matrix(omega, chi, phi, SNS_or_IPNS)
    
    newq = numpy.zeros((3))
    newq = numpy.dot(q, fmat)

    return newq

#---------------------------------------------
#       function rotation_matrix
#---------------------------------------------
# Calculate rotation matrix for omega, chi and phi angles
#---------------------------------------------
def rotation_matrix(omega, chi, phi, SNS_or_IPNS):
    "Returns rotation matrix from setting angles"

    rad = 180. / math.pi

    ph = phi / rad
    cp = math.cos(ph)
    sp = math.sin(ph)
    R_phi = numpy.zeros((3,3))
    R_phi[0,0] = cp
    R_phi[0,1] = sp
    R_phi[1,0] = -sp
    R_phi[1,1] = cp
    R_phi[2,2] = 1.0

    ch = chi / rad        #changed -chi to chi, 8/23/07
    cc = math.cos(ch)
    sc = math.sin(ch)
    R_chi = numpy.zeros((3,3))
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
    R_om = numpy.zeros((3,3))
    R_om[0,0] = co
    R_om[0,1] = so
    R_om[1,0] = -so
    R_om[1,1] = co
    R_om[2,2] = 1.0

    fmat = numpy.zeros((3,3))
    fmat = numpy.dot(R_phi, R_chi)
    fmat = numpy.dot(fmat, R_om)

    return fmat

#--------------------------------------------------------
#               function UB_IPNS_2_SNS
#--------------------------------------------------------
#  Transform the IPNS UB matrix usually store in the mat
#  file to the SNS matrix.
#  These are actually the UB tramspose matrices. The
#  transformation is:
#  IPNS:
#     col1  col2  col3
#  SNS:
#     col2  col3  col1
#--------------------------------------------------------
def UB_IPNS_2_SNS(UB_IPNS):
    "Transform the IPNS UB matrix to the SNS matrix."

    UB_SNS = numpy.zeros((3,3))
    for i in range(3):
        UB_SNS[i,0] = UB_IPNS[i,1]
        UB_SNS[i,1] = UB_IPNS[i,2]
        UB_SNS[i,2] = UB_IPNS[i,0]
    
    return UB_SNS

#################################################################
# test
if __name__ == '__main__':

    # nickel 5678_Niggle.mat
    UB = numpy.zeros((3,3))
    UB[0,] = [ 0.28403938,  0.01181680, -0.00120478 ]
    UB[1.] = [ 0.01218532, -0.28391245,  0.00737687 ]
    UB[2,] = [ 0.00095225, -0.00789868, -0.28362876 ]
    print('\nUB:')
    print(UB)
    
    omega = 144.0
    chi = 134.8
    phi = -0.02
    
    SNS_or_IPNS = 'SNS'
    
    newmat = rotate_matrix(UB, omega, chi, phi, SNS_or_IPNS)
    print('\nnewmat:')
    print(newmat)
    
    
    

        
