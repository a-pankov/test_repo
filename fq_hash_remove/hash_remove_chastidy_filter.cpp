#include "FastQStatus.h"
#include "FastQFile.h"
#include "StringBasics.h"
#include <iostream>

static int usage()
{
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   fq_hash_remove <fastq/fastq.gz file> \n");
        fprintf(stderr, "\n");
        return 1;
}
 
int main(int argc, char ** argv)
{

   FastQFile fastQFile;
   fastQFile.disableMessages();
   fastQFile.disableSeqIDCheck();
  

   String filename;

   if (argc == 1) return usage();
   else filename = argv[1];

   // Open the fastqfile with the default UNKNOWN space type which will determine the 
   // base type from the first character in the sequence.
   if(fastQFile.openFile(filename) != FastQStatus::FASTQ_SUCCESS)
   {
      // Failed to open the specified file.
      // Report the error and exit (handled by error).
//      error("Failed to open file: %s", filename.c_str());
      std::cerr << "Failed to open file: " << filename << "\n" ;
      exit(EXIT_FAILURE);
   }
  unsigned long total_valid(0), total_all(0), total_chastidy_filter(0);
   // Keep reading the file until there are no more fastq sequences to process.
   while (fastQFile.keepReadingFile() && ++total_all)
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
      if(fastQFile.readFastQSequence() == FastQStatus::FASTQ_SUCCESS)
      {
         vector< String > * split_vec(fastQFile.mySequenceIdLine.Split(':'));
         String chastidy_filt = (*split_vec)[7];
//         std::cout << (chastidy_filt == "Y") << '\t';
         if( chastidy_filt == "Y" && ++total_chastidy_filter ) continue;


         // The sequence is valid.
         int start(0), end( fastQFile.myQualityString.Length());
         for(;start != end && (fastQFile.myQualityString[start] == '#'); ++start){}
         if(start == end) continue;
         else {
            --end;
            for(;end != start && (fastQFile.myQualityString[end] == '#'); --end){}
            std::cout << fastQFile.mySequenceIdLine << '\n' << fastQFile.myRawSequence.Mid(start, end) << '\n' << fastQFile.myPlusLine << '\n' << fastQFile.myQualityString.Mid(start, end) << '\n';
            ++total_valid;
         }

//         <Your Processing Here>
         // For example if you want to print the lines of the sequence:
         
//         printf("The Sequence ID Line is: %s", fastQFile.mySequenceIdLine.c_str());
//         printf("The Sequence ID is: %s", fastQFile.mySequenceIdentifier.c_str());
//         printf("The Sequence Line is: %s", fastQFile.myRawSequence.c_str());
//         printf("The Plus Line is: %s", fastQFile.myPlusLine.c_str());
//         printf("The Quality String Line is: %s", fastQFile.myQualityString.c_str());
      }
   }
   // Finished processing all of the sequences in the file.
   // Close the input file.
   fastQFile.closeFile();
   std::cerr << "Clipped # ends from " << argv[1] << "\nTotal number of reads:\t" << total_all << "\nTotal number of kept reads:\t" << total_valid << "\nTotal number of chastidy filtered reads:\t" << total_chastidy_filter << '\n';
   return 0; // It is up to you to determine your return.
 }
