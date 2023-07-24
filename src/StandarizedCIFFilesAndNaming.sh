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
# Transform to P1 with RASPA:
loc=$(pwd)
back_IZA="Original_Structures"
make_INPUT
if [ ! -d SuperCell_P1 ] ; then mkdir SuperCell_P1 ; fi
if [ ! -d ${back_IZA} ] ; then mkdir ${back_IZA} ; fi
for file in *_${1}.cif ; do
 name=$(echo $file | sed "s/\_${1}.cif//g")
 abc=$(echo $name | sed 's/_/ /g' | awk '{print $1}')
 sed "s/STRUCTURE/${name}_${1}/g"  INPUT > simulation.input
 sed -i -e 's/(unknown)/P 1/g' -e '/_symmetry_equiv_pos_site_id/d' -e '/_symmetry_equiv_pos_as_xyz/d' $file
#
 a=$(grep "_cell_length_a" $file | awk '{print $2}' | sed 's/(//g' | sed 's/)//g')
 b=$(grep "_cell_length_b" $file | awk '{print $2}' | sed 's/(//g' | sed 's/)//g')
 c=$(grep "_cell_length_c" $file | awk '{print $2}' | sed 's/(//g' | sed 's/)//g')
 u_a=1  ; while [ $(echo "${u_a} * $a < ${min_cell_size} " | bc -lq) == 1 ] ; do let u_a+=1 ; done
 u_b=1  ; while [ $(echo "${u_b} * $b < ${min_cell_size}" | bc -lq) == 1 ] ; do let u_b+=1 ; done
 u_c=1  ; while [ $(echo "${u_c} * $c < ${min_cell_size}" | bc -lq) == 1 ] ; do let u_c+=1 ; done
 echo $abc: $a $b $c : $u_a $u_b $u_c : 
 sed -i "s/UnitCells        Ua Ub Uc/UnitCells        $u_a $u_b $u_c/g" simulation.input
 mkdir tmp
 cd tmp
  cp ../simulation.input ../$file .
  simulate simulation.input
  cp -rf Movies/System_0/Framework_0_final_${u_a}_${u_b}_${u_c}_P1.cif ../SuperCell_P1/${abc}.cif
 cd ${loc}
 rm -rf tmp
 mv $file ${back_IZA}/.
done
rm -rf INPUT simulation.input
exit 0
