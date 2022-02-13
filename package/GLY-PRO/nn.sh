mkdir GLY-PRO
cd GLY-PRO
cp ../package/GLY-PRO/type.raw .
cp ../input.GLY-PRO .
cp ../package/GLY-PRO/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf GLY-PRO
