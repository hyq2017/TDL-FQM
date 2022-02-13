mkdir GLN
cd GLN
cp ../package/GLN/type.raw .
cp ../input.GLN .
cp ../package/GLN/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf GLN
