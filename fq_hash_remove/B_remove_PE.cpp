//
//  main.cpp
//  hash_remove_PE
//
//  Created by Alex on 11/12/13.
//  Copyright (c) 2013 Alex. All rights reserved.
//
#include "FastQStatus.h"
#include "FastQFile.h"
#include <iostream>
#include "gzstream.h"
//#include <zlib.h>
#include <fstream>
//#include <sstream>
//
//#include <boost/iostreams/filtering_streambuf.hpp>
//#include <boost/iostreams/copy.hpp>
//#include <boost/iostreams/filter/gzip.hpp>
//#include <boost/iostreams/stream.hpp>


static int usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:   fq_hash_remove_PE <first mate fastq/fastq.gz file> <second mate fastq/fastq.gz file> <Output directory and prefix; e.g. ./test_Reads> \n");
  fprintf(stderr, "\n");
  return 1;
}

int main(int argc, char ** argv)
{
  
  std::string filename_1, filename_2, prefix, of1, of2;
  
  if (argc != 4) return usage();
  else {
    filename_1 = argv[1];
    filename_2 = argv[2];
    prefix = argv[3];
  }
  of1 = prefix + "_1.fq.gz";
  of2 = prefix + "_2.fq.gz";

  FastQFile fastQFile_1(0), fastQFile_2(0);
   
  fastQFile_1.disableMessages();
  fastQFile_2.disableMessages();
  fastQFile_1.disableSeqIDCheck();
  fastQFile_2.disableSeqIDCheck();
  
  // Open the fastqfile with the default UNKNOWN space type which will determine the
  // base type from the first character in the sequence.
  if(fastQFile_1.openFile(filename_1.c_str()) != FastQStatus::FASTQ_SUCCESS)
  {
    std::cerr << "Failed to open file: " << filename_1 << "\n" ;
    exit(EXIT_FAILURE);
  }
  if(fastQFile_2.openFile(filename_2.c_str()) != FastQStatus::FASTQ_SUCCESS)
  {
    std::cerr << "Failed to open file: " << filename_2 << "\n" ;
    exit(EXIT_FAILURE);
  }
  
  /*****************************************************
   Writing to Gzip
   ******************************************************/
  
  ogzstream rout_1( of1.c_str());

  ogzstream rout_2( of2.c_str());

  /*****************************************************
   ******************************************************/

  
  unsigned long total_valid(0), total_all(0);
  // Keep reading the file until there are no more fastq sequences to process.
  while (fastQFile_1.keepReadingFile() && fastQFile_2.keepReadingFile() && ++total_all)
  {
        // Read one sequence. This call will read all the lines for
        // one sequence.
        /////////////////////////////////////////////////////////////////
        // NOTE: It is up to you if you want to process only for success:
        //    if(readFastQSequence() == FASTQ_SUCCESS)
        // or for FASTQ_SUCCESS and FASTQ_INVALID:
        //    if(readFastQSequence() != FASTQ_FAILURE)
        // Do NOT try to process on a FASTQ_FAILURE
        /////////////////////////////////////////////////////////////////
    if(fastQFile_1.readFastQSequence() == FastQStatus::FASTQ_SUCCESS && fastQFile_2.readFastQSequence() == FastQStatus::FASTQ_SUCCESS){
          // Both sequences are valid.
      int start_1(0), end_1( fastQFile_1.myQualityString.Length()), start_2(0), end_2( fastQFile_2.myQualityString.Length());
      for(;start_1 != end_1 && (fastQFile_1.myQualityString[start_1] == 'B'); ++start_1){}
      if(start_1 == end_1) continue;
      else {
        
        for(;start_2 != end_2 && (fastQFile_2.myQualityString[start_2] == 'B'); ++start_2){}
        if(start_2 == end_2) continue;
        else {
          --end_2;
          for(;end_2 != start_2 && (fastQFile_2.myQualityString[end_2] == 'B'); --end_2){}
          rout_2 << fastQFile_2.mySequenceIdLine << '\n' << fastQFile_2.myRawSequence.Mid(start_2, end_2) << '\n' << fastQFile_2.myPlusLine << '\n' << 
fastQFile_2.myQualityString.Mid(start_2, end_2) <<'\n';
//rout_2.flush();
          
          --end_1;
          for(;end_1 != start_1 && (fastQFile_1.myQualityString[end_1] == 'B'); --end_1){}
          rout_1 << fastQFile_1.mySequenceIdLine << '\n' << fastQFile_1.myRawSequence.Mid(start_1, end_1) << '\n' << fastQFile_1.myPlusLine << '\n' << 
fastQFile_1.myQualityString.Mid(start_1, end_1) << '\n';
// rout_1.flush();
          
          ++total_valid;

        }
      }
    }
  }
    // Finished processing all of the sequences in the file.
    // Close the input file.
  rout_1.close();
  rout_2.close();
  fastQFile_1.closeFile();
  fastQFile_2.closeFile();
  std::cout << "Clipped # ends from " << filename_1 << " and " << filename_2 << "\nTotal number of alignments:\t" << total_all << "\nTotal number of kept alignments:\t" << total_valid << '\n';
  return 0; // It is up to you to determine your return.
}
