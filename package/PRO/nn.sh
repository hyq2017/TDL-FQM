mkdir PRO
cd PRO
cp ../package/PRO/type.raw .
cp ../input.PRO .
cp ../package/PRO/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf PRO
