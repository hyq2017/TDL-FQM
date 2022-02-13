mkdir ASN
cd ASN
cp ../package/ASN/type.raw .
cp ../input.ASN .
cp ../package/ASN/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf ASN
