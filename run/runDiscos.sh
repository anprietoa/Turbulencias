g++ -g -O3 ../src/Discos.cpp
./a.out 1 > mdisc.g
./a.out 0 > mdisc.dat

gnuplot mdisc.g

mv mdisc.gif ../outp/
mv mdisc.dat ../outp/

rm mdisc.g a.out