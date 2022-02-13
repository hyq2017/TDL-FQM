mkdir MET
cd MET
cp ../package/MET/type.raw .
cp ../input.MET .
cp ../package/MET/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf MET
