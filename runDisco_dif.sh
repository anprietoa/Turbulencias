g++ -g -O3 ../src/Disco_difusion.cpp
./a.out 1 > dif.g
./a.out 0 > dif.dat

gnuplot dif.g

mv dif.gif ../outp/
mv dif.dat ../outp/

rm dif.g a.out
