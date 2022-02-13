mkdir ILE
cd ILE
cp ../package/ILE/type.raw .
cp ../input.ILE .
cp ../package/ILE/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf ILE
