#!/bin/bash

##### INPUT PARAMETERS ##########################################

# Set both the theta and phi angles to simulate
PhiArray=("$1")
ThetaArray=("$2")

# Settings of the system being investigated
adsorbate="$3"
surface="$4"
element="$5"

# Get the atom directories
Array=( "$(ls -d "${element}"[0-9]{1,4}$)" )

###### GET XPS BINDING ENERGIES #####################################

# Set the system name
system="${adsorbate}_${surface}"

# Read the XPS binding energies and store in an array
readarray -t XPSArray <../XPS/"${element}"_XPS_peaks.txt

# Search through all the atom directories and add the XPS energy to each of the .molpdos file
for i in "${!Array[@]}"; do (
  cd "${Array[i]}" || exit
  sed -i "s/nexafs_xshift        :  [0-9]*.[0-9]*/nexafs_xshift        :  ${XPSArray[i]}/g" "${system}".molpdos
) done

###### RUN MOLPDOS ######################################################

# Enter all of the directories and run MolPDOS for each angle stated and printing out the progress
for atom in "${Array[@]}"; do (
  echo "$atom"
  for theta in "${ThetaArray[@]}"; do (
    for phi in "${PhiArray[@]}"; do (
      cd "$atom" || continue
      echo "$atom" t"$theta" p"$phi"
      mkdir t"${theta}"_p"$phi"
      sed -i "s/nexafs_phi           :   [0-9]*/nexafs_phi           :   $phi/g" "${system}".molpdos
      sed -i "s/nexafs_theta         :   [0-9]*/nexafs_theta         :   $theta/g" "${system}".molpdos
      MolPDOS "$system"
      mv ./*.dat t"${theta}"_p"$phi"
      echo "$atom" t"$theta" p"$phi" "done"
    ) done
  ) done
  echo "$atom" "done"
) done
