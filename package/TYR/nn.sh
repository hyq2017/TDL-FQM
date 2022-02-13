mkdir TYR
cd TYR
cp ../package/TYR/type.raw .
cp ../input.TYR .
cp ../package/TYR/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf TYR
