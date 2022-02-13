mkdir ASP
cd ASP
cp ../package/ASP/type.raw .
cp ../input.ASP .
cp ../package/ASP/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf ASP
