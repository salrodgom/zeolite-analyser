#!/bin/bash
min_cell_size=10.0
function make_INPUT {
echo "SimulationType                MonteCarlo
NumberOfCycles                0
NumberOfInitializationCycles  0
PrintEvery                    0

Forcefield                    CrystalGenerator
RemoveAtomNumberCodeFromLabel yes

Framework 0
FrameworkName    STRUCTURE
UnitCells        Ua Ub Uc
" > INPUT
}
# Download: Force-Brute from europe.iza-structure.org/IZA-SC
# wait 30 min.
#for a in "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N" "O" "P" "Q" "R" "S" "T" "U" "V" "W" "X" "Y" "Z" ; do
# for b in "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N" "O" "P" "Q" "R" "S" "T" "U" "V" "W" "X" "Y" "Z" ; do
#  for c in "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N" "O" "P" "Q" "R" "S" "T" "U" "V" "W" "X" "Y" "Z" ; do
#   ABC=${a}${b}${c}
#   echo $ABC
#   wget https://europe.iza-structure.org/IZA-SC/cif/${ABC}.cif
#   sleep 0.1
#  done
# done
#done
# Back the original CIF Files:
# Transform to P1 with RASPA:
loc=$(pwd)
back="Original_Structures"
make_INPUT
if [ ! -d ReferenceDataBase ] ; then mkdir ReferenceDataBase ; fi
if [ ! -d ${back} ] ; then mkdir ${back} ; fi
for file in *.cif ; do
 name=$(echo $file | sed "s/\.cif//g")
 abc=$(echo $name | sed 's/_/ /g' | awk '{print $1}')
 sed "s/STRUCTURE/${name}/g"  INPUT > simulation.input
#_audit_creation_method             'generated by GULP'
#_symmetry_space_group_name_H-M     '(unknown)       '
#_symmetry_Int_Tables_number        1
#loop_
#_symmetry_equiv_pos_site_id
#_symmetry_equiv_pos_as_xyz
#_cell_length_a                       16.1703
#_cell_length_b                       16.1703
#_cell_length_c                       16.1703
#_cell_angle_alpha                  109.4712
#_cell_angle_beta                   109.4712
#_cell_angle_gamma                  109.4712
 sed -i -e 's/(unknown)/P 1/g' -e '/_symmetry_equiv_pos_site_id/d' -e '/_symmetry_equiv_pos_as_xyz/d' $file
#
 a=$(grep "_cell_length_a" $file | awk '{print $2*1.02}' | sed 's/(//g' | sed 's/)//g')
 b=$(grep "_cell_length_b" $file | awk '{print $2*1.02}' | sed 's/(//g' | sed 's/)//g')
 c=$(grep "_cell_length_c" $file | awk '{print $2*1.02}' | sed 's/(//g' | sed 's/)//g')
 u_a=1  ; while [ $(echo "${u_a} * $a < ${min_cell_size} " | bc -lq) == 1 ] ; do let u_a+=1 ; done
 u_b=1  ; while [ $(echo "${u_b} * $b < ${min_cell_size}" | bc -lq) == 1 ] ; do let u_b+=1 ; done
 u_c=1  ; while [ $(echo "${u_c} * $c < ${min_cell_size}" | bc -lq) == 1 ] ; do let u_c+=1 ; done
 echo $abc: $a $b $c : $u_a $u_b $u_c : 
 sed -i "s/UnitCells        Ua Ub Uc/UnitCells        $u_a $u_b $u_c/g" simulation.input
 mkdir tmp
 cd tmp
  cp ../simulation.input ../$file .
  simulate simulation.input
  cp -rf Movies/System_0/Framework_0_final_${u_a}_${u_b}_${u_c}_P1.cif ../ReferenceDataBase/${name}.cif
 cd ${loc}
 rm -rf tmp
 mv $file ${back}/.
done
rm -rf INPUT simulation.input
exit 0
