mkdir ALA-GLY
cd ALA-GLY
cp ../package/ALA-GLY/type.raw .
cp ../input.ALA-GLY .
cp ../package/ALA-GLY/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf ALA-GLY
