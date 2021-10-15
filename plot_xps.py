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

###################################
xstart = 285.
xstop = 305.
broad1 = 0.7 
broad2 = 0.7 
firstpeak = 285.0
ewid1 = firstpeak+1.0
ewid2 = firstpeak+2.0
mix1 = 0.3
mix2 = 0.3
########################################

#Set what element you have calculated XPS for
element = 'C'
#Read in the XPS peaks in generated with python script
data = np.loadtxt(element+'_XPS_peaks.txt')
print(data)

#Apply the broadening
x, y = dos_binning(data, broadening=broad1, mix1=mix1, mix2=mix2, start=xstart, stop=xstop,
                coeffs = None, broadening2=broad2, ewid1=ewid1, ewid2=ewid2)

#Write out the spectrum to a text file
fileout = open(element+'_XPS_spectrum.txt', 'w')
for (xi, yi) in zip(x,y):
    dat = str(xi) + ' ' + str(yi) + '\n'
    fileout.write(dat)
fileout.close()

#To get the indivdual atom peaks uncomment the assert 0
assert 0

xs = []
ys = []

for z in range(10):
    peak = []
    peak.append(data[z])
    x_tmp, y_tmp = dos_binning(peak, broadening=broad1, mix1=mix1, mix2=mix2, start=xstart, stop=xstop,
                coeffs = None, broadening2=broad2, ewid1=ewid1, ewid2=ewid2)
    xs.append(x_tmp)
    ys.append(y_tmp)

    txtfile = open(element+'_XPS_spectrum_'+element+str(z)+'.txt','w')
    for (xsz, ysz) in zip(x_tmp, y_tmp):
        txt = str(xsz) + ' ' + str(ysz) + '\n'
        txtfile.write(txt)
    txtfile.close()
