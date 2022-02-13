mkdir LEU
cd LEU
cp ../package/LEU/type.raw .
cp ../input.LEU .
cp ../package/LEU/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf LEU
