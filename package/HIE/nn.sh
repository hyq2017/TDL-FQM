mkdir HIE
cd HIE
cp ../package/HIE/type.raw .
cp ../input.HIE .
cp ../package/HIE/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf HIE
