g++ -g -O3 ../src/Disco_Fcentral.cpp
./a.out 1 > fcen.g
./a.out 0 > fcen.dat

gnuplot fcen.g

mv fcen.gif ../outp/
mv fcen.dat ../outp/

rm fcen.g a.out
