g++ -std=c++0x -D WINDOWS "NetworkBatch.cpp" -o "bin\net.exe" -lgsl -lgslcblas -static
g++ -std=c++0x -D WINDOWS "NeuronBatch.cpp" -o "bin\sn.exe" -lgsl -lgslcblas -static
