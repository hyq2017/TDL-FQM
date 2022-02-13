mkdir GLY-ALA
cd GLY-ALA
cp ../package/GLY-ALA/type.raw .
cp ../input.GLY-ALA .
cp ../package/GLY-ALA/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf GLY-ALA
