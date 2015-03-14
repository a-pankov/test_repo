//
//  subset_name.cpp
//  subset by read name
//
//  Created by Alex on 1/21/14.
//  Copyright (c) 2014 Alex. All rights reserved.
//

#include <iostream>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <tr1/unordered_map>

#include <ctime>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

struct passwd *pw = getpwuid(getuid());
std::string homedir = pw->pw_dir;


using namespace BamTools;


void create_sam_index( std::string fname){
    std::cerr << "\tIndexing " + fname << "...\n";
    std::string sys_cmd = "/home/apankov/bin/samtools index " + fname;
    const char * c = sys_cmd.c_str();
    system(c);
}



int main(int argc, const char * argv[])
{
    clock_t t1 = std::clock();
    time_t tm = std::time(NULL);
    std::srand(static_cast<unsigned int>( tm ));
    std::cerr << std::asctime(std::localtime(&tm)) << "\n";
    
    using std::string;
    string inputFileName, fileName_subset, prefix, read_name;
    
    if ( argc == 2 )
    {
        inputFileName  = argv[1] ;
    }
    else {
        std::cerr << "Wrong number of arguments!\n\tname_subset <bam file> <read name> <output prefix> " << std::endl;
        exit(EXIT_FAILURE);
    }
    
    //  using namespace BamTools;
    
    //    Open files for reading and writing
    
    BamReader reader;
    if (!reader.Open(inputFileName) ) {
        std::cerr << "Cant open " + inputFileName +"\n"   ;
        exit(EXIT_FAILURE);
    }
    else{
        std::cerr << "Reading in BAM file:" + inputFileName + "\n";
    }
  
    
  BamAlignment  al;
  unsigned long long total(0);
  unsigned short found(0);
  std::tr1::unordered_map<int32_t, unsigned long long> length_counts_map ;
  
    while ( reader.GetNextAlignment(al) && ++total) {
      ++length_counts_map[al.Length];
    }    
       
    for(auto it = length_counts_map.begin(); it != length_counts_map.end(); ++it){
      std::cout << it->first << '\t' << it->second << '\n';
    }

    reader.Close();
    
    std::cerr << "\nTotal count: " << total << "\nFound: " << found << "reads.";
  
    std::cerr << "\nRunning Time in secs: " << (float) (std::clock() - t1)/CLOCKS_PER_SEC << '\n' << std::flush;
    
    return 0;
}

