g++ -g -O3 ../src/Disco_Magnus.cpp
./a.out 1 > mag.g
./a.out 0 > mag.dat

gnuplot mag.g

mv mag.gif ../outp/
mv mag.dat ../outp/

rm mag.g a.out
