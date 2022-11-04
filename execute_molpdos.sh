#!/bin/bash

# Set the range of numbers of the atom directories
Array=($(seq 48 1 57))

# Set the angles you want to simulate
AngleArray=("00")

# Change to the MOLECULE and METAL of the system for filename
molecule="azulene"
metal="Ag"
element="C"

system="${molecule}_${metal}"

# Read the XPS binding energies and store in an array
readarray -t ArrayX < ../XPS/${element}_XPS_peaks.txt

# Search through all the directories and add the XPS energy to each of the .molpdos file
for index in ${!Array[@]}; do
    cd ${element}${Array[$index]}/
    sed -i "s/nexafs_xshift        :  [0-9]*.[0-9]*/nexafs_xshift        :  ${ArrayX[$index]}/g" ${system}.molpdos >> ${system}.molpdos
    cd ../
done

# Enter all of the directories and run MolPDOS for each ANGLE stated and printing out the progress
for atom in ${Array[@]}; do
  echo $atom
  for angle in ${AngleArray[@]}; do
      cd ${element}$atom/
      echo $atom $angle
      mkdir t$angle
      sed -i "s/nexafs_theta         :   [0-9]*/nexafs_theta         :   $angle/g" ${system}.molpdos >> ${system}.molpdos
      MolPDOS $system
      mv *.dat t$angle/
      echo $atom $angle "done"
      cd ../
  done
  echo $atom "done"
done
