mkdir ARG
cd ARG
cp ../package/ARG/type.raw .
cp ../input.ARG .
cp ../package/ARG/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf ARG
