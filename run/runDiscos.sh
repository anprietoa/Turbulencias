g++ -g -O3 ../src/Discos.cpp
./a.out 2 > mdisc.g

mv mdisc.dat ../outp/

gnuplot mdisc.g
mv mdisc.gif ../outp/

rm mdisc.g a.out