mkdir VAL
cd VAL
cp ../package/VAL/type.raw .
cp ../input.VAL .
cp ../package/VAL/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf VAL
