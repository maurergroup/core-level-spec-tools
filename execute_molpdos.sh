#!/bin/bash

##### INPUT PARAMETERS ##########################################

# Set the range of numbers of the atom directories
Array=($(seq 48 1 57))

# Set both the theta and phi angles you want to simulate
ThetaArray=("00" "25" "53" "90")
PhiArray=("60")

# Settings of the system being investigated
molecule="azulene"
metal="Ag"
element="C"

###### GET XPS BINDING ENERGIES#####################################

# Set the system name
system="${molecule}_${metal}"

# Read the XPS binding energies and store in an array
readarray -t XPSArray < ../XPS/${element}_XPS_peaks.txt

# Search through all the atom directories and add the XPS energy to each of the .molpdos file
for index in ${!Array[@]}; do
    cd ${element}${Array[$index]}/
    sed -i "s/nexafs_xshift        :  [0-9]*.[0-9]*/nexafs_xshift        :  ${XPSArray[$index]}/g" ${system}.molpdos >> ${system}.molpdos
    cd ../
done

###### RUN MOLPDOS ######################################################

# Enter all of the directories and run MolPDOS for each angle stated and printing out the progress
for atom in ${Array[@]}; do
  echo $element$atom
  for theta in ${ThetaArray[@]}; do
    for phi in ${PhiArray[@]}; do
      cd ${element}$atom/
      echo $element$atom t$theta p$phi
      mkdir t${theta}_p$phi
      sed -i "s/nexafs_phi           :   [0-9]*/nexafs_phi           :   $phi/g" ${system}.molpdos >> ${system}.molpdos
      sed -i "s/nexafs_theta         :   [0-9]*/nexafs_theta         :   $theta/g" ${system}.molpdos >> ${system}.molpdos
      MolPDOS $system
      mv *.dat t${theta}_p$phi/
      echo $element$atom t$theta p$phi "done"
      cd ../
    done
  done
  echo $element$atom "done"
done
