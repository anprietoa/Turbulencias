g++ ../src/Paper.cpp 
./a.out 1 velMAG.dat 
./a.out 0 velNM.dat

python ../src/Paper.py

mv velMAG.dat ../outp/
mv velNM.dat ../outp/
mv Corr.png ../outp/

rm a.out
