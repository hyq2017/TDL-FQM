mkdir GLU
cd GLU
cp ../package/GLU/type.raw .
cp ../input.GLU .
cp ../package/GLU/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf GLU
