mkdir PRO-ALA
cd PRO-ALA
cp ../package/PRO-ALA/type.raw .
cp ../input.PRO-ALA .
cp ../package/PRO-ALA/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf PRO-ALA
