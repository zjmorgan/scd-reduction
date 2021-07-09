import math
def integrate_1d_peak( y, peakMin, peakMax, bkgMin, bkgMax):
    """
    Obtain the integrated intensity and the sigma of a peak
    in Q space.
    """
    peak = 0.0                # peak counts
    bkg = 0.0                 # background counts
    numPeakCh = peakMax - peakMin + 1            # number of peak channels
    numBkgCh =(bkgMax - bkgMin + 1) - numPeakCh  # number of background channels
            
    for i in range(bkgMin, peakMin):
        bkg = bkg + y[i]
    
    for i in range(peakMin, peakMax+1):
        peak = peak + y[i]
        
    for i in range(peakMax+1, bkgMax+1):
        bkg = bkg + y[i]

    ratio = float(numPeakCh)/float(numBkgCh)
    
    intI = peak - bkg*(ratio)
    sigI = math.sqrt(peak + bkg*ratio**2)
    
    return intI, sigI

# test
if __name__ == '__main__':

    import math
    
    input = open('profile_0_-3_-12.txt', 'r')
    y = []
    for i in range(10):
        lineString = input.readline()
        lineList = lineString.split()
        if len(lineList) == 0: break
        for j in range(10):
            y.append(int(lineList[j]))

    # print y
    ysum = sum(y)
    print(ysum)
    ysum = math.fsum(y)
    print(ysum)
    
    peakMin = 20
    peakMax = 79
    bkgMin = 0
    bkgMax = 99
    
    intI, sigI = integrate_1d_peak( y, peakMin, peakMax, bkgMin, bkgMax)
    
    print(intI, sigI)
    