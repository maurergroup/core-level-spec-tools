#!/bin/bash

declare -a Array=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9")
declare -a AngleArray=("00")

molecule="azulene"
metal="gas"

system="${molecule}_${metal}"

readarray -t Array1 < ../XPS/XPS_peaks.txt

for index in ${!Array[@]}; do
    cd C${Array[$index]}/
    sed -i "s/nexafs_xshift        :  [0-9]*.[0-9]*/nexafs_xshift        :  ${Array1[$index]}/g" ${system}.molpdos >> ${system}.molpdos
    cd ../
done

for atom in ${Array[@]}; do
#  fileName=$(basename $file)
  echo $atom
#  var=$(pwd)
#  echo $var
  for angle in ${AngleArray[@]}; do
      cd C$atom/
#      rm -r t$angle
      echo $atom $angle
#      var=$(pwd)
#      echo $var
      mkdir t$angle
      sed -i "s/nexafs_theta         :   [0-9]*/nexafs_theta         :   $angle/g" ${system}.molpdos >> ${system}.molpdos
      MolPDOS $system
      mv *.dat t$angle/
      echo $atom $angle "done"
      cd ../
  done
  echo $atom "done"
done
