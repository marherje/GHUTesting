#!/bin/bash

rm *.txt
#grep -r "SM:" ../DataChi2_headerA.txt > Model_SMA.txt
#sed -i "s/SM: //g" Model_SMA.txt

#grep -r "SM:" ../DataChi2_headerB_typocorrection_couplingsA.txt > Model_SMB.txt
#sed -i "s/SM: //g" Model_SMB.txt

grep -r "SM:" ../NewData.txt > Model_SMA.txt
sed -i "s/SM: //g" Model_SMA.txt

cp Model_SMA.txt Model_SMB.txt

for model in 1 2 3 4 5 6 7 8
do
    grep -r "Model "${model}":" ../NewData.txt > Model_${model}.txt
    sed -i "s/Model "${model}": //g" Model_${model}.txt
done

mv Model_1.txt Model_Am.txt
mv Model_2.txt Model_Ap.txt
mv Model_3.txt Model_Bm.txt
mv Model_4.txt Model_Bp.txt
mv Model_5.txt Model_Cm.txt
mv Model_6.txt Model_Cp.txt
mv Model_7.txt Model_A2.txt
mv Model_8.txt Model_A3.txt
