mkdir ALA
cd ALA
cp ../package/ALA/type.raw .
cp ../input.ALA .
cp ../package/ALA/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf ALA
