mkdir PRO-GLY
cd PRO-GLY
cp ../package/PRO-GLY/type.raw .
cp ../input.PRO-GLY .
cp ../package/PRO-GLY/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf PRO-GLY
