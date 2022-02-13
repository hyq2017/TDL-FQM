mkdir CYS
cd CYS
cp ../package/CYS/type.raw .
cp ../input.CYS .
cp ../package/CYS/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf CYS
