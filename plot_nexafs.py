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

#########################################
#Broadening Paramters
xstart = 285.
xstop = 310.
broad1 = 0.75
broad2 = 2.0
firstpeak = 290.0 
ewid1 = firstpeak+5.0
ewid2 = firstpeak+15.0
mix1 = 0.2
mix2 = 0.8
########################################

#File parameter
angle = ['t25','t53','t90']
molecule = 'benzene'
metal = 'Cu'

filename = '/'+molecule+'_'+metal+'_4_1_1_1_deltas.dat'

peaks = np.zeros([6,38048]) #Empty array of No. of carbons and No. of states
I = np.zeros([6,38048])

for a in angle:
    for i,direc in enumerate(['C1/','C2/','C3/','C4/','C5/','C6/']):

        data = np.loadtxt(direc+a+filename)
        x, y = data[:,0], data[:,3]
        peaks[i,:] = x
        I[i,:] = y

    fileout = open(molecule+'_'+metal+'_deltas_'+a+'.txt','w')
    fileout.write('#   <x in eV>     Intensity\n')
    for p,i in zip(peaks.flatten(), I.flatten()):
        fileout.write('{0:16.8f}    {1:16.8f}\n'.format(p,i))
    fileout.close()

    x, y = dos_binning(peaks.flatten(), broadening=broad1, mix1=mix1, mix2=mix2, start=xstart, stop=xstop,
        coeffs = I.flatten(), broadening2=broad2, ewid1=ewid1, ewid2=ewid2)

    datafile = open(molecule+'_'+metal+'_spectrum_'+a+'.txt', 'w')
    for (xi, yi) in zip(x,y):
        asd = str(xi) + ' ' + str(yi) + '\n'
        datafile.write(asd)
    datafile.close()
