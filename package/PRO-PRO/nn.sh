mkdir PRO-PRO
cd PRO-PRO
cp ../package/PRO-PRO/type.raw .
cp ../input.PRO-PRO .
cp ../package/PRO-PRO/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf PRO-PRO
