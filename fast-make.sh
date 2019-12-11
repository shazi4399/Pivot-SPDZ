make clean
make lib
#make mpir
make -j 8 Setup.x
make -j 8 Player-Online.x
make -j 8 Fake-Offline.x
make -j 8 semi-party.x
make -j 8 vfl-logistic-func-test.x
make -j 8 vfl-decision-tree-test.x
./Setup.x 3 128 128
Scripts/setup-online.sh 3 128 128
