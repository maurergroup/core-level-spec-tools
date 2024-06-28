import numpy as np

###### BROADENING PARAMETER ###################################

X_START = 285. # Start Value
X_STOP = 305. # End Value
BROAD = 0.7 # Broadening value for first section 
MIX = 0.3 # G/L mix ratio

###### INPUT PARAMETERS ########################################

ELEMENT = 'C'
ATOM_CONTRIBUTE = False

################################################################

def main():

# Read in the XPS peaks in generated with python script
    XPS_DATA = np.loadtxt(ELEMENT+'_XPS_peaks.txt')
    print(XPS_DATA)

# Apply the broadening
    x, y = dos_binning(XPS_DATA, BROADENING=BROAD, mix=MIX, START=X_START, STOP=X_STOP,
                COEFFS = None)

# Write out the spectrum to a text file
    spec_file = open(ELEMENT+'_XPS_spectrum.txt', 'w')
    for (xi, yi) in zip(x,y):
        SPEC_DAT = str(xi) + ' ' + str(yi) + '\n'
        spec_file.write(SPEC_DAT)
    spec_file.close()

    if ATOM_CONTRIBUTE == True:
        xs = []
        ys = []

        for z in range(len(XPS_DATA)):
            peak = []
            peak.append(XPS_DATA[z])
            x_tmp, y_tmp = dos_binning(peak, BROADENING=BROAD, mix=MIX, START=X_START, STOP=X_STOP,
                        COEFFS = None,)
            xs.append(x_tmp)
            ys.append(y_tmp)

            ATOM_file = open(ELEMENT+'_XPS_spectrum_'+ELEMENT+str(z)+'.txt','w')
            for (xsz, ysz) in zip(x_tmp, y_tmp):
                ATOM_DATA = str(xsz) + ' ' + str(ysz) + '\n'
                ATOM_file.write(ATOM_DATA)
            ATOM_file.close()
        else:
            quit()

##############################################################

def gaussian(X, X_MEAN, BROADENING):
    
    GAUSSIAN_VAL = np.sqrt((4*np.log(2))/(np.pi*(BROADENING**2)))* np.exp(-((4*np.log(2))/(BROADENING**2))*(X-X_MEAN)**2)
    return GAUSSIAN_VAL

def lorentzian(X, X_MEAN, BROADENING):

    LORENTZIAN_VAL = (1/(2*np.pi))* (BROADENING)/(((BROADENING/2)**2)+(X-X_MEAN)**2)
    return LORENTZIAN_VAL

def PseudoVoigt(X, X_MEAN, BROADENING, MIXING):
    """ 
    Combines gaussian and lorentzian schemes together
    """
    return (1-MIXING)*gaussian(X, X_MEAN, BROADENING)+MIXING*lorentzian(X, X_MEAN, BROADENING)

def dos_binning(EIGENVALUES, BROADENING=0.75, BIN_WIDTH=0.01, MIX=0.,
        COEFFS=None, START=0.0, STOP=10.0):
    """ 
    performs binning for a given set of EIGENVALUES and 
    optionally weight COEFFS.
    """
    if COEFFS is None:
        COEFFS = np.ones(len(EIGENVALUES))

    LOWEST_E = START
    HIGHEST_E = STOP
    NUM_BINS = int((HIGHEST_E-LOWEST_E)/BIN_WIDTH)
    X_AXIS = np.zeros([NUM_BINS])
    data = np.zeros([NUM_BINS])
    #setting up x-axis
    for i in range(NUM_BINS):
        X_AXIS[i] = LOWEST_E + i * BIN_WIDTH
    #get DOS

    SIGMA = BROADENING
    MIXING = MIX

    for i in range(NUM_BINS):
        PSEUDO_VOIGT_VEC = np.zeros((len(EIGENVALUES)))
        PSEUDO_VOIGT_VEC=PseudoVoigt(X_AXIS[i],EIGENVALUES,SIGMA,MIXING)*COEFFS
        data[i]= np.sum(PSEUDO_VOIGT_VEC)
    return X_AXIS, data

main()