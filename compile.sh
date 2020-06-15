g++ -std=c++0x -D LINUX "NetworkBatch.cpp" -o "./bin/net.out" -lgsl -lgslcblas -static
g++ -std=c++0x -D LINUX "NeuronBatch.cpp" -o "./bin/sn.out" -lgsl -lgslcblas -static
