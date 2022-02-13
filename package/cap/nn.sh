mkdir cap
cd cap
cp ../package/cap/type.raw .
cp ../input.cap .
cp ../package/cap/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf cap
