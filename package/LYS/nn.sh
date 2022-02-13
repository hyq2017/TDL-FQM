mkdir LYS
cd LYS
cp ../package/LYS/type.raw .
cp ../input.LYS .
cp ../package/LYS/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf LYS
