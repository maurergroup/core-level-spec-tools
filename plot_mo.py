import numpy as np

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

######BROADENING PARAMETERS################################

xstart = 285.
xstop = 310.
broad1 = 0.75
broad2 = 2.0
firstpeak = 290.0 
ewid1 = firstpeak+5.0
ewid2 = firstpeak+15.0
mix1 = 0.2
mix2 = 0.8

######SYSTEM_PARAMTERS######################################

pol_method = 4 #1 for Total NEXAFS, 2 for angular, 3 for polarised, 4 for average polarised
MO = list(map(str, range(17,29))) #List of MO numbers you want 17-29 gives 17, 18 ,19... 28
angle = ['t25','t53','t90'] #Incidence angles
molecule = 'azulene' #Name of molecule
metal = 'Ag' #Surface in system
elem = 'C'
numbers = list(range(48,58)) #Set range to correspond to the directories C48, C49... C57
atom = '4' #The number of the exicted atom in the list of elements in the system, always the last, if system contains H, C, Ag, C:exc, it will be 4

######SETUP ALL LIST AND VARIABLES#############################

#Create list of all the folders all the data is in C48/, C49/... C57/
folders = []
for n in numbers:
    folders.append(elem+str(n)+'/')

#Create variable with a string og the delta file which will be read
filename = '/'+molecule+'_'+metal+'_'+atom+'_1_1_1_deltas.dat'

#Get the number of kpoints used in calculation in order to correct the MO projected state
kpts = []
with open(elem+str(numbers[0])+'/'+molecule+'_'+metal+'.bands', 'r') as bands:
    for line in bands:
        if 'Number of k-points' in line:
            for word in line.split():
                try:
                    kpts.append(float(word))
                except ValueError:
                    pass

#Get the length of the deltas file
bands = np.loadtxt(elem+str(numbers[0])+'/'+angle[0]+'/'+filename)
bands_num = len(bands)

#Read .param file to see if calculation is spin_polarised and set up required settings
with open (elem+str(numbers[0])+'/'+molecule+'_'+metal+'.param') as param:
    if 'SPIN_POLARIZED: TRUE' in param.read():
        spin = True
        spin_val = list(map(str,range(1,3)))
        spin_num - bands_num/2
    else:
        spin = False
        spin_val = list(map(str, range(1,2)))
        spin_num = bands_num

#Create arrays with sized of the system to use
peaks = np.zeros([len(numbers),int(spin_num)])
I = np.zeros([len(numbers), int(spin_num)])


###############################################

for s in spin_val:
    for m in MO:
        for a in angle:
            for i,direc in enumerate(folders):
            
                data = np.loadtxt(direc+a+filename)
                x, y = data[:,0], data[:,pol_method]

                if spin == True:
                    x1, y1 = x[:int(spin_num)], y[:int(spin_num)]
                    x2, y2 = x[int(spin_num):], y[int(spin_num):]
                    spindict = {
                            'spin1x' : x1,
                            'spin2y' : y1,
                            'spin2x' : x2,
                            'spin2y' : y2
                            }
                else:
                    spindict = {
                            'spin1x' : x,
                            'spin1y' : y,
                            }

                data2 = np.loadtxt(direc+a+'/'+molecule+'_'+metal+'_'+m+'_spin'+s+'_deltas.dat')
                data2*= kpts
                peaks[i,:] = spindict['spin'+s+'x']
                I[i,:] = spindict['spin'+s+'y']*data2[:,1]
        
            fileout = open(molecule+'_'+metal+'_MO'+m+'_deltas_'+a+'_spin'+s+'.txt','w')
            fileout.write('#   <x in eV>     Intensity\n')
            for p,i in zip(peaks.flatten(), I.flatten()):
                fileout.write('{0:16.8f}    {1:16.8f}\n'.format(p,i))
            fileout.close()
            
            x, y = dos_binning(peaks.flatten(), broadening=broad1, mix1=mix1, mix2=mix2, start=xstart, stop=xstop,
                    coeffs = I.flatten(), broadening2=broad2, ewid1=ewid1, ewid2=ewid2)
        
            datafile = open(molecule+'_'+metal+'_MO'+m+'_'+a+'spin'+s+'.txt', 'w')
            for (xi, yi) in zip(x,y):
                asd = str(xi) + ' ' + str(yi) + '\n'
                datafile.write(asd)
            datafile.close()

