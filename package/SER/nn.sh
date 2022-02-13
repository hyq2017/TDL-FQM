mkdir SER
cd SER
cp ../package/SER/type.raw .
cp ../input.SER .
cp ../package/SER/energy_force.py .
python energy_force.py
mv *.log ../
cd ..
rm -rf SER
