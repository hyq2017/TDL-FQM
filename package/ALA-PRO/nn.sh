mkdir ALA-PRO
cd ALA-PRO
cp ../package/ALA-PRO/type.raw .
cp ../input.ALA-PRO .
cp ../package/ALA-PRO/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf ALA-PRO
