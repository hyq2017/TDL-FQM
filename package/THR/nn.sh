mkdir THR
cd THR
cp ../package/THR/type.raw .
cp ../input.THR .
cp ../package/THR/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf THR
