import numpy as np

######BROADENING PARAMETERS################################

xstart = 275. #Start Value
xstop = 330. #End Value
broad1 = 0.75 #Broadening value for first section
broad2 = 2.0 #Broadening value for last section
firstpeak = 290.0 
ewid1 = firstpeak+5.0 #Set range to linearly move from braod1 to broad2
ewid2 = firstpeak+15.0
mix1 = 0.2 #First G/L mix ratio
mix2 = 0.8 #Last G/L mix ratio

######SYSTEM PARAMTERS######################################

molecule = 'azulene' #Name of molecule
metal = 'Ag' #Surface in system
element = 'C'
num_start = 48 #First index number of the atom directories
num_end = 57 #Last index number of the atom directories
n_type = 4 #1 for Total NEXAFS, 2 for angular, 3 for polarised, 4 for average polarised
angle = ['t25','t53','t90'] #Incidence angles
atom = '4' #The number of the excited atom in the list of elements in the system
MO_start = 17 #First MO state to project out
MO_end = 28 #Last MO state to project

######SETUP ALL LIST AND VARIABLES#############################

#Create list of all the atoms and MO states
numbers = list(range(num_start,num_end+1))
MO = list(map(str, range(MO_start,MO_end+1)))

#Create list of all the folders all the data is in C48/, C49/... C57/
folders = []
for n in numbers:
    folders.append(element+str(n)+'/')

#Create variable with a string og the delta file which will be read
filename = '/'+molecule+'_'+metal+'_'+atom+'_1_1_1_deltas.dat'

#Get the number of kpoints used in calculation in order to correct the MO projected state
kpts = []
with open(element+str(numbers[0])+'/'+molecule+'_'+metal+'.bands', 'r') as bands:
    for line in bands:
        if 'Number of k-points' in line:
            for word in line.split():
                try:
                    kpts.append(float(word))
                except ValueError:
                    pass

#Get the length of the deltas file
bands = np.loadtxt(element+str(numbers[0])+'/'+angle[0]+'/'+filename)
bands_num = len(bands)

#Read .param file to see if calculation is spin_polarised and set up required settings
with open (element+str(numbers[0])+'/'+molecule+'_'+metal+'.param') as param:
    if 'SPIN_POLARIZED: TRUE' in param.read():
        spin = True
        spin_val = list(map(str,range(1,3)))
        spin_num = bands_num/2
    else:
        spin = False
        spin_val = list(map(str, range(1,2)))
        spin_num = bands_num

#Create arrays with sized of the system to use
peaks = np.zeros([len(numbers),int(spin_num)])
I = np.zeros([len(numbers), int(spin_num)])


###############################################
def main():
#Loop over spin, MOs and angles in the indivdual directories
    for s in spin_val:
        for m in MO:
            for a in angle:
                for i,direc in enumerate(folders):
#Load the data from the MolPDOS calculation
                    nex_data = np.loadtxt(direc+a+filename)
                    x, y = nex_data[:,0], nex_data[:,n_type]
#If spin polarized is on then split the data in half for each spin, x1, y1 and x2, y2
                    if spin == True:
                        x1, y1 = x[:int(spin_num)], y[:int(spin_num)]
                        x2, y2 = x[int(spin_num):], y[int(spin_num):]
                        spindict = {
                                'spin1x' : x1,
                                'spin2y' : y1,
                                'spin2x' : x2,
                                'spin2y' : y2
                                }
#If not spin polarized load all data as x, y
                    else:
                        spindict = {
                                'spin1x' : x,
                                'spin1y' : y,
                                }
#Load the indivdual MO data from the MolPDOS calculation and add the k-point scaling and add
#multiply the intesity of the MO with the overall spectrum
                    nex_data2 = np.loadtxt(direc+a+'/'+molecule+'_'+metal+'_'+m+'_spin'+s+'_deltas.dat')
                    nex_data2*= kpts
                    peaks[i,:] = spindict['spin'+s+'x']
                    I[i,:] = spindict['spin'+s+'y']*nex_data2[:,1]
#Write out all of the MO data into a delta file
                mo_del_file = open(molecule+'_'+metal+'_MO'+m+'_deltas_'+a+'_spin'+s+'.txt','w')
                mo_del_file.write('#   <x in eV>     Intensity\n')
                for p,i in zip(peaks.flatten(), I.flatten()):
                    mo_del_file.write('{0:16.8f}    {1:16.8f}\n'.format(p,i))
                mo_del_file.close()
#Apply the broadening
                x, y = dos_binning(peaks.flatten(), broadening=broad1, mix1=mix1, mix2=mix2, start=xstart, stop=xstop,
                        coeffs = I.flatten(), broadening2=broad2, ewid1=ewid1, ewid2=ewid2)
#Write out MO peak into a text file
                mo_file = open(molecule+'_'+metal+'_MO'+m+'_'+a+'_spin'+s+'.txt', 'w')
                for (xi, yi) in zip(x,y):
                    mo_data = str(xi) + ' ' + str(yi) + '\n'
                    mo_file.write(mo_data)
                mo_file.close()

##############################################################
def gaussian(x, x_mean, broadening):

    gaussian_val = np.sqrt((4*np.log(2))/(np.pi*(broadening**2)))* np.exp(-((4*np.log(2))/(broadening**2))*(x-x_mean)**2);
    return gaussian_val

def lorentzian(x, x_mean, broadening):

    lorentzian_val = (1/(2*np.pi))* (broadening)/(((broadening/2)**2)+(x-x_mean)**2);
    return lorentzian_val

def PseudoVoigt(x, x_mean, broadening, mixing):
    """ 
    Combines gaussian and lorentzian schemes together
    """
    return (1-mixing)*gaussian(x, x_mean, broadening)+mixing*lorentzian(x, x_mean, broadening)

def dos_binning(eigenvalues,broadening=0.75, bin_width=0.01, mix1=0., mix2 = None,
        coeffs=None,start=0.0, stop=10.0, broadening2 = None, ewid1 = 10.0, ewid2 = 20.0):
    """ 
    performs binning for a given set of eigenvalues and 
    optionally weight coeffs.
    """
    if broadening2 is None:
        broadening2 = broadening
    if coeffs is None:
        coeffs = np.ones(len(eigenvalues))
    lowest_e = start
    highest_e = stop
    num_bins = int((highest_e-lowest_e)/bin_width)
    x_axis = np.zeros([num_bins])
    data = np.zeros([num_bins])
    #setting up x-axis
    for i in range(num_bins):
        x_axis[i] = lowest_e + i * bin_width
    #get DOS
    sigma=np.zeros((len(eigenvalues)))
    mixing=np.zeros((len(eigenvalues)))

    for ei,e in enumerate(eigenvalues):
        if e<=(ewid1):
            sigma[ei]=broadening
            mixing[ei]=mix1
        elif e>(ewid2):
            sigma[ei]=broadening2
            mixing[ei]=mix2
        else:
            sigma[ei]=broadening + ((broadening2-broadening)/(ewid2-ewid1))*(e-ewid1)
            mixing[ei]=(mix1 + ((mix2-mix1)/(ewid2-ewid1))*(e-ewid1))
    for i in range(num_bins):
        pseudovoigt_vec = np.zeros((len(eigenvalues)))
        pseudovoigt_vec=PseudoVoigt(x_axis[i],eigenvalues,sigma,mixing)*coeffs
        data[i]= np.sum(pseudovoigt_vec)
    return x_axis, data
