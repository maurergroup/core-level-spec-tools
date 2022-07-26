import numpy as np

######BROADENING PARAMETERS###################################

xstart = 285. #Start Value
xstop = 305. #End Value
broad = 0.7 #Broadening value for first section 
mix = 0.3 #G/L mix ratio

##############################################################

atom_contribute = False

def main():
#Set what element you have calculated XPS for
    element = 'C'
#Read in the XPS peaks in generated with python script
    xps_data = np.loadtxt(element+'_XPS_peaks.txt')
    print(xps_data)

#Apply the broadening
    x, y = dos_binning(xps_data, broadening=broad, mix=mix, start=xstart, stop=xstop,
                coeffs = None)

#Write out the spectrum to a text file
    spec_file = open(element+'_XPS_spectrum.txt', 'w')
    for (xi, yi) in zip(x,y):
        spec_dat = str(xi) + ' ' + str(yi) + '\n'
        spec_file.write(spec_dat)
    spec_file.close()

    if atom_contribute == True:
        xs = []
        ys = []

        for z in range(len(xps_data)):
            peak = []
            peak.append(xps_data[z])
            x_tmp, y_tmp = dos_binning(peak, broadening=broad, mix=mix, start=xstart, stop=xstop,
                        coeffs = None,)
            xs.append(x_tmp)
            ys.append(y_tmp)

            atom_file = open(element+'_XPS_spectrum_'+element+str(z)+'.txt','w')
            for (xsz, ysz) in zip(x_tmp, y_tmp):
                atom_data = str(xsz) + ' ' + str(ysz) + '\n'
                atom_file.write(atom_data)
            atom_file.close()
        else:
            quit()

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

def dos_binning(eigenvalues, broadening=0.75, bin_width=0.01, mix=0.,
        coeffs=None, start=0.0, stop=10.0):
    """ 
    performs binning for a given set of eigenvalues and 
    optionally weight coeffs.
    """
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
    sigma=broadening
    mixing=mix

    for i in range(num_bins):
        pseudovoigt_vec = np.zeros((len(eigenvalues)))
        pseudovoigt_vec=PseudoVoigt(x_axis[i],eigenvalues,sigma,mixing)*coeffs
        data[i]= np.sum(pseudovoigt_vec)
    return x_axis, data

main()