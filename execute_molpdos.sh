#!/bin/bash

#Set the range of of numbers of the element directories
#seq(first step last) so seq 0 1 9) is 0, 1, 2... 9
#Set the angles wanted to run MolPDOS for
declare -a Array=($(seq 0 1 9))
declare -a AngleArray=("00")

#Change to the molecule and metal of the system for filename
molecule="azulene"
metal="gas"
element="C"

system="${molecule}_${metal}"

#Read the XPS binding energies and store in an array
readarray -t ArrayX < ../XPS/${element}_XPS_peaks.txt

#Search through all the directories and add the XPS energy to each of the .molpdos file
for index in ${!Array[@]}; do
    cd ${element}${Array[$index]}/
    sed -i "s/nexafs_xshift        :  [0-9]*.[0-9]*/nexafs_xshift        :  ${ArrayX[$index]}/g" ${system}.molpdos >> ${system}.molpdos
    cd ../
done

#Enter all of the directories and run MolPDOS for each angle stated
for number in ${Array[@]}; do
#  fileName=$(basename $file)
  echo $number
#  var=$(pwd)
#  echo $var
  for angle in ${AngleArray[@]}; do
      cd ${element}$number/
#      rm -r t$angle
      echo $number $angle
#      var=$(pwd)
#      echo $var
      mkdir t$angle
      sed -i "s/nexafs_theta         :   [0-9]*/nexafs_theta         :   $angle/g" ${system}.molpdos >> ${system}.molpdos
      MolPDOS $system
      mv *.dat t$angle/
      echo $number $angle "done"
      cd ../
  done
  echo $number "done"
done
