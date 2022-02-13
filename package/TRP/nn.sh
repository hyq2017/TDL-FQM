mkdir TRP
cd TRP
cp ../package/TRP/type.raw .
cp ../input.TRP .
cp ../package/TRP/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf TRP
