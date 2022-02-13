mkdir PHE
cd PHE
cp ../package/PHE/type.raw .
cp ../input.PHE .
cp ../package/PHE/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf PHE
