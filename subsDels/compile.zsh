g++ -std=c++0x -O3 -I /home/apankov/chiaPetGraph/bamtools/include -L /home/apankov/chiaPetGraph/bamtools/lib -Wl,-rpath,/home/apankov/chiaPetGraph/bamtools/lib -lbamtools -o subsDels main.cpp genomeFile.cpp