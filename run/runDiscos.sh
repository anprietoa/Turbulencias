g++ -g -O3 ../src/Discos.cpp
./a.out 2 > mdisc.g

gnuplot mdisc.g

mv mdisc.gif ../outp/
mv mdisc.dat ../outp/

rm mdisc.g a.out