mkdir ALA-ALA
cd ALA-ALA
cp ../package/ALA-ALA/type.raw .
cp ../input.ALA-ALA .
cp ../package/ALA-ALA/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf ALA-ALA
