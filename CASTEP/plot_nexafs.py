import numpy as np

###### BROADENING PARAMETERS ###################################

# Start and end values of spectra
X_START = 275.
X_STOP = 330.
# Broadeining parameters for the first and last ranges
BROAD_1 = 0.75
BROAD_2 = 2.0
# Energy of the lead peak
FIRST_PEAK = 290.0 
# Set the start and end point of the linearly increasing broadening
# change from the first and last ranges with respect to the leading
# peak
EWID_1 = FIRST_PEAK + 5.0
EWID_2 = FIRST_PEAK + 15.0
# Set the Gaussian/Lorentzian mixing ratio for the ranges
MIX_1 = 0.2
MIX_2 = 0.8

###### SYSTEM PARAMETERS ########################################

# Setting of the system being investigated
MOLECULE = 'azulene'
METAL = 'Ag'
ELEMENT = 'C'
# Index range of the atom directories created by autoscript.py 
NUM_START = 48
NUM_END = 57
# Type of NEXAFS spectrum to output
# 1 for total summed NEXAFS, 2 for angular, 3 for polarised, and
# 4 for average polarised
N_TYPE = 4
# The theta and phi angles simulated
THETA_ANGLE = ['00','25','53','90']
PHI_ANGLE = ['60']
# The element index of the excited atom all the elements in the system
# always the last element in the system, so if system contains H, C, Ag, C:exc
# it will be 4
ATOM = '4'
# Whether you want to output the individual atom contributions
ATOM_contribute = False

##### SETUP ALL LIST AND VARIABLES NEEDED #######################

# Create a list of all the atoms
NUMBERS = list(range(NUM_START, NUM_END + 1))

# Set up a list of all the directories all the data is in C48/, C49/... C57/
FOLDERS = []
for n in NUMBERS:
    FOLDERS.append(ELEMENT + str(n) + '/')

# Create variable with a string for the delta file which will be read
FILE_NAME = '/' + MOLECULE + '_' + METAL + '_' + ATOM + '_1_1_1_deltas.dat'

# Get the length of the deltas file
BANDS = np.loadtxt(ELEMENT + str(NUMBERS[0]) + '/t' + THETA_ANGLE[0] + '_p' + PHI_ANGLE[0] + FILE_NAME)
BANDS_NUM = len(BANDS)

# Create arrays with sizes of the system to use
PEAKS = np.zeros([len(NUMBERS),BANDS_NUM])
I = np.zeros([len(NUMBERS),BANDS_NUM])

###########################################################
def main():
    # Loop over both the theta and phi angles and the atom directories
    for t_a in THETA_ANGLE:
        for p_a in PHI_ANGLE:
            for i,direc in enumerate(FOLDERS):
                # Load the data from the MolPDOS file
                NEX_DATA = np.loadtxt(direc + 't' + t_a + '_p' + p_a + FILE_NAME)
                X, Y = NEX_DATA[:,0], NEX_DATA[:,N_TYPE]
                PEAKS[i,:] = X
                I[i,:] = Y
            # Write out all of the data into a delta peaks file
            NEX_DEL_FILE = open(MOLECULE + '_' + METAL + '_deltas_t' + t_a + '_p' + p_a + '.txt','w')
            NEX_DEL_FILE.write('#   <x in eV>     Intensity\n')
            for p,i in zip(PEAKS.flatten(), I.flatten()):
                NEX_DEL_FILE.write('{0:16.8f}    {1:16.8f}\n'.format(p,i))
            NEX_DEL_FILE.close()
            # Apply the pseudovoigt broadening to the data
            X, Y = dos_binning(PEAKS.flatten(), BROADENING=BROAD_1, MIX_1=MIX_1, MIX_2=MIX_2, START=X_START, STOP=X_STOP,
                COEFFS = I.flatten(), BROADENING_2=BROAD_2, EWID_1=EWID_1, EWID_2=EWID_2)
            # Write out spectrum into a text file
            NEX_SPEC_FILE = open(MOLECULE + '_' + METAL + '_spectrum_t' + t_a + '_p' + p_a + '.txt', 'w')
            for (X_I, Y_I) in zip(X,Y):
                NEX_DATA = str(X_I) + ' ' + str(Y_I) + '\n'
                NEX_SPEC_FILE.write(NEX_DATA)
            NEX_SPEC_FILE.close()
            # Calculates the individual atom contributions of the NEXAFS spectra if selected above
            if ATOM_contribute == True:
                XS = []
                YS = []
                for z in range(len(NUMBERS)):
                    X_TMP, Y_TMP = dos_binning(PEAKS[z,:], BROADENING=BROAD_1, MIX_1=MIX_1, MIX_2=MIX_2, START=X_START, STOP=X_STOP,
                            COEFFS = I[z,:], BROADENING_2=BROAD_2, EWID_1=EWID_1, EWID_2=EWID_2)
                    XS.append(X_TMP)
                    YS.append(Y_TMP)

                    ATOM_FILE = open(MOLECULE + '_' + METAL + '_' + ELEMENT + str(z) + '_t' + t_a + '_p' + p_a + '.txt', 'w')
                    for (XS_I, YS_I) in zip(X_TMP, Y_TMP):
                        ATOM_DATA = str(XS_I) + ' ' + str(YS_I) + '\n'
                        ATOM_FILE.write(ATOM_DATA)
                    ATOM_FILE.close()

############################################################
def gaussian(X, X_MEAN, BROADENING):

    GAUSSIAN_VAL = np.sqrt((4*np.log(2))/(np.pi*(BROADENING**2)))* np.exp(-((4*np.log(2))/(BROADENING**2))*(X-X_MEAN)**2);
    return GAUSSIAN_VAL

def lorentzian(X, X_MEAN, BROADENING):

    LORENTZIAN_VAL = (1/(2*np.pi))* (BROADENING)/(((BROADENING/2)**2)+(X-X_MEAN)**2);
    return LORENTZIAN_VAL

def PseudoVoigt(X, X_MEAN, BROADENING, MIXING):
    """ 
    Combines gaussian and lorentzian schemes together
    """
    return (1-MIXING)*gaussian(X, X_MEAN, BROADENING)+MIXING*lorentzian(X, X_MEAN, BROADENING)

def dos_binning(EIGENVALUES,BROADENING=0.75, BIN_WIDTH=0.01, MIX_1=0., MIX_2 = None,
        COEFFS=None,START=0.0, STOP=10.0, BROADENING_2 = None, EWID_1 = 10.0, EWID_2 = 20.0):
    """ 
    performs binning for a given set of EIGENVALUES and 
    optionally weight COEFFS.
    """
    if BROADENING_2 is None:
        BROADENING_2 = BROADENING
    if COEFFS is None:
        COEFFS = np.ones(len(EIGENVALUES))

    LOWEST_E = START
    HIGHEST_E = STOP
    NUM_BINS = int((HIGHEST_E-LOWEST_E)/BIN_WIDTH)
    X_AXIS = np.zeros([NUM_BINS])
    data = np.zeros([NUM_BINS])
    # setting up x-axis
    for i in range(NUM_BINS):
        X_AXIS[i] = LOWEST_E + i * BIN_WIDTH
    # get DOS
    SIGMA = np.zeros((len(EIGENVALUES)))
    MIXING = np.zeros((len(EIGENVALUES)))

    for ei,e in enumerate(EIGENVALUES):
        if e<=(EWID_1):
            SIGMA[ei]=BROADENING
            MIXING[ei]=MIX_1
        elif e>(EWID_2):
            SIGMA[ei]=BROADENING_2
            MIXING[ei]=MIX_2
        else:
            SIGMA[ei]=BROADENING + ((BROADENING_2-BROADENING)/(EWID_2-EWID_1))*(e-EWID_1)
            MIXING[ei]=(MIX_1 + ((MIX_2-MIX_1)/(EWID_2-EWID_1))*(e-EWID_1))

    for i in range(NUM_BINS):
        PSEUDO_VOIGT_VEC = np.zeros((len(EIGENVALUES)))
        PSEUDO_VOIGT_VEC=PseudoVoigt(X_AXIS[i],EIGENVALUES,SIGMA,MIXING)*COEFFS
        data[i]= np.sum(PSEUDO_VOIGT_VEC)
    return X_AXIS, data

main()
