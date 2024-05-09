#!/bin/bash
# run a parameter study based on materialA
rm -rf study
mkdir -p study
echo '(cd ~/Code/Calibration/ && git update)' >> study/run.sh
echo '(cd ~/Code/Calibration/build && cmake .)' >> study/run.sh
echo '(cd ~/Code/Calibration/build/Drivers/Calibration && make -j)' >> study/run.sh
chmod +x study/run.sh

for Exe0 in CalibrationHeap CalibrationDrum CalibrationShearCell; do
rc0=0.5
mus0=0.5
mur0=0.1
bo0=0.0
Exe='./'$Exe0' -species LinearViscoelasticFrictionReversibleAdhesiveSpecies -density 1500 -collisionTime 0.000067997 -torsionFriction 0 -normalStress 1000 800 600 400 200 -output -restitutionCoefficient '$rc0' -slidingFriction '$mus0' -rollingFriction '$mur0' -bondNumber '$bo' -psd logNormal volume diameter 15e-4 '

for var in '0e-4' '1e-4' '2e-4' '4e-4' '8e-4'; do
  ExeParam=$Exe$var' -param _P'$var
  echo $ExeParam >> runDoE_B2.sh
  ExeOrig='~/Code/Calibration/build/Drivers/Calibration/'$Exe0
  sed "s|REPLACE_EXE_ORIG|$ExeOrig|g;s|REPLACE_EXE_PARAM|$ExeParam|g" slurmScript.sh > study/$Exe0$var
  echo 'sbatch -J '$Exe0$var' '$Exe0$var >> study/run.sh
  echo 'sleep 1' >> study/run.sh
done
done

scp -r study $msm3:
