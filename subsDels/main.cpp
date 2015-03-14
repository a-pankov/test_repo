//
//  main.cpp
//  myPileup
//
//  Created by Alex Pankov on 9/19/13.
//  Copyright (c) 2013 Alex Pankov. All rights reserved.
//
//
///*************************************************************************///
//  CHROMOSOME COORDINATES are 0-BASED!!!!!!!!!!!!!
//
//  Output File Format - Tab-Delimited ('\t') with the following columns:
//  1. Alignment ID
//  2. Chromosome (string)
//  3. Coordinate along the chromosome (0-Based)
//  4. Strand (string)
//  5. Reference Base
//  6. Substitution Base or "DEL" for deletions
//  7. Phred Quality of base in the query read or Phred Qualities of flanking bases for deletions
//  8. Position along the read
//
///*************************************************************************///
//

#include <iostream>
#include <sstream>
#include <fstream>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <vector>
#include <tr1/unordered_map>

#include <ctime>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include "genomeFile.h"
#include "bedFile.h"
#include "lineFileUtilities.h"

struct passwd *pw = getpwuid(getuid());
std::string homedir = pw->pw_dir;

void CoverageBam( std::string, std::ofstream & );
void GetBamBlocks(const BamAlignment &bam,
                  const string &chrom,
                  bedVector &bedBlocks, std::ofstream & fout);

void AddCoverage(int start, int end, vector<DEPTH> & _currChromCoverage);
void AddBlockedCoverage(const vector<BED> &bedBlocks, vector<DEPTH> & _currChromCoverage, int _currChromSize);
void ReportChromCoverage(const vector<DEPTH> &chromCov, const int &chromSize, const string &chrom);
void printSub( const BamAlignment &bam, const string &chrom, int currPosition , int readPos, char sub, std::ofstream &fout);





int main(int argc, const char * argv[])
{
  
  clock_t t1 = std::clock();
  time_t tm = std::time(NULL);
  std::srand(static_cast<unsigned int>( tm ));
  std::cerr << std::asctime(std::localtime(&tm)) << "\n";
  
  using std::string;
  string inputFileName, outputFileName;
  
  if ( argc == 3 )
  {
    inputFileName = argv[1];
    outputFileName = argv[2];
    
  }
  else {
    std::cerr << "Wrong number of arguments, BAM file and Output file name must follow command invocation." << std::endl;
//    inputFileName = "/Users/alexp/data/blelloch/RBP_seq_2013/tests/../alignments/R-Pos_accepted_hits_forwardStrand.bam";
    exit(EXIT_FAILURE);
  }
  
  std::ofstream fout;
  fout.open(outputFileName.c_str());
  
  if (fout.is_open()) {
    CoverageBam(inputFileName, fout);
    fout.close();
  } else
  {
    std::cerr << "Cannot open " << outputFileName << ".\n";
    exit(EXIT_FAILURE);
  }
  

  std::cerr << "\nRunning Time in secs: " << (float) (std::clock() - t1)/CLOCKS_PER_SEC << '\n';
  
}





void CoverageBam(std::string bamFile, std::ofstream & fout) {
  
  using namespace BamTools;
  using std::string;
  
  string _currChromName = "";
  int _currChromSize = 0 ;
  std::vector<DEPTH> _currChromCoverage{};
  
  
// open the BAM file
  BamReader reader;
  if (!reader.Open(bamFile)) {
    std::cerr << "Failed to open BAM file " << bamFile << std::endl;
    exit(1);
  }
  
// get header & reference information
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();
  
// load the BAM header references into a BEDTools "genome file"
  GenomeFile * _genome = new GenomeFile(refs);
  
// convert each aligned BAM entry to BED
// and compute coverage on B
  BamAlignment bam;
//  int cnt(0);
  while (reader.GetNextAlignment(bam) /*&& (++cnt) < 100*/ ) {
// skip if the read is unaligned
    if (bam.IsMapped() == false)
      continue;
    
// extract the chrom, start and end from the BAM alignment
    string chrom(refs.at(bam.RefID).RefName);
    
// are we on a new chromosome?
    if ( chrom != _currChromName ){
      if (_currChromName.length() > 0) {
        ReportChromCoverage(_currChromCoverage, _currChromSize, _currChromName);
      }
      
// empty the previous chromosome and reserve new
      std::vector<DEPTH>().swap(_currChromCoverage);
      
      _currChromName = chrom;
      
// get the current chrom size and allocate space
      _currChromSize = _genome->getChromSize(_currChromName);
      
      if (_currChromSize >= 0)
        _currChromCoverage.resize(_currChromSize);
      else {
        cerr << "Input error: Chromosome " << _currChromName << " found in your input file but not in your genome file." << endl;
        exit(1);
      }

    }
    
// add coverage accordingly.
    bedVector bedBlocks;
    GetBamBlocks(bam, refs.at(bam.RefID).RefName, bedBlocks, fout);
    AddBlockedCoverage(bedBlocks, _currChromCoverage, _currChromSize);
  }
// close the BAM
  reader.Close();
  ReportChromCoverage(_currChromCoverage, _currChromSize, _currChromName);
}

void GetBamBlocks(const BamAlignment &bam,
                  const string &chrom,
                  bedVector &bedBlocks, std::ofstream & fout)
{
  vector<int> starts;
  vector<int> lengths;
  starts.push_back(0);
  
  string strand, md;
  char sub;
  bam.IsReverseStrand() ? strand = "-" : strand = "+";
  CHRPOS currPosition = bam.Position;
  int blockLength  = 0, nm, readPos(0), mdInt, cigMatched(0), loc(0);
  
  bam.GetTag("NM",  nm);

//  if (nm > 0) {
//    bam.GetTag("MD",  md);
//    stringstream ss;
//    ss << md;
//    int number;
//    ss >> number;
//    std::cout << number << '\t';
////    std::cout << md << '\t';
//    while (!ss.eof()) {
//      char c = ss.peek();
//      if (isdigit(c)) {
////        int number;
//        ss >> number;
//        std::cout << '\t' << number << '\t';
//      } else {
//        ss >> c;
//        std::cout << c;
//      }
//    }
//    
//    std::cout << '\n';
//  } else {
//    
//  }
  
  
//  Rip through the CIGAR ops and figure out if there is more
//  than one block for this alignment
  stringstream ss;
  if (nm > 0) {
    bam.GetTag("MD",  md);
    ss << md;
    ss >> mdInt;
    readPos += mdInt;
  }
  vector<CigarOp>::const_iterator cigItr = bam.CigarData.begin();
  vector<CigarOp>::const_iterator cigEnd = bam.CigarData.end();
  for (; cigItr != cigEnd; ++cigItr) {
    switch (cigItr->Type) {
      case ('M') :
      {
        blockLength += cigItr->Length;
        cigMatched += cigItr->Length;
        if( nm >0){
          while (readPos < cigMatched) {
            loc += mdInt;
            ss >> sub;
            printSub(bam, chrom, currPosition + loc , readPos++, sub, fout);
            loc++;
            ss >> mdInt;
            readPos += mdInt;
          }
          loc = 0;
        }
      }
      break;
      case ('I') : 
        readPos += cigItr->Length;
        break;
      case ('S') : break;
      case ('D') :
      {
        bedBlocks.push_back( BED(chrom, currPosition, currPosition + blockLength,
                                 bam.Name, ToString(bam.MapQuality), strand) );
        
        fout  /* << bam.Name << '\t' */ << chrom << '\t' << currPosition + blockLength << '\t' << strand << '\t';
        currPosition += cigItr->Length + blockLength;
        int cnt = cigItr->Length;
        ss >> sub;
        if (sub == '^') {
          ss >> sub;
        }
        while (--cnt) {
          fout << sub;
          ss >> sub;
        }
        fout << sub;
        if (isdigit(ss.peek())) {
          ss >> mdInt;
          fout <<  "\tDEL\t" << bam.Qualities.at(readPos-1) << bam.Qualities.at(readPos) << "\t" << readPos << '\n';  // if out of bounds, check conditions
          readPos += mdInt;
        }
        blockLength = 0;
      }
      break;
      case ('P') : break;
      case ('N') :
      {
        bedBlocks.push_back( BED(chrom, currPosition, currPosition + blockLength,
                                 bam.Name, ToString(bam.MapQuality), strand) );
        currPosition += cigItr->Length + blockLength;
        blockLength = 0;
      }
      break;
      case ('H') : break;                             // for 'H' - do nothing, move to next op
      default    :
      {
        printf("ERROR: Invalid Cigar op type\n");   // shouldn't get here
        exit(1);
      }
      break;
    }
  }
  bedBlocks.push_back( BED(chrom, currPosition, currPosition + blockLength,
                           bam.Name, ToString(bam.MapQuality), strand) );
}


void AddCoverage(int start, int end, vector<DEPTH> & _currChromCoverage, int _currChromSize) {
  // process the first line for this chromosome.
  // make sure the coordinates fit within the chrom
  if (start < _currChromSize)
    _currChromCoverage[start].starts++;
  if (end < _currChromSize)
    _currChromCoverage[end].ends++;
  else
    _currChromCoverage[_currChromSize-1].ends++;
}


void AddBlockedCoverage(const vector<BED> &bedBlocks, vector<DEPTH> & _currChromCoverage, int _currChromSize) {
  vector<BED>::const_iterator bedItr = bedBlocks.begin();
  vector<BED>::const_iterator bedEnd = bedBlocks.end();
  for (; bedItr != bedEnd; ++bedItr) {
    // the end - 1 must be done because BamAncillary::getBamBlocks
    // returns ends uncorrected for the genomeCoverageBed data structure.
    // ugly, but necessary.
    AddCoverage(bedItr->start, bedItr->end - 1, _currChromCoverage, _currChromSize);
  }
}

void printSub( const BamAlignment &bam, const string &chrom, int currPosition , int readPos, char sub, std::ofstream & fout){
  fout /* << bam.Name << '\t' */  << chrom << '\t' << currPosition <<'\t' << (bam.IsReverseStrand() ? "-" : "+") << '\t' << sub << '\t' << bam.QueryBases.at(readPos) << '\t' << bam.Qualities.at(readPos) << '\t' << readPos << '\n';
}

void ReportChromCoverage(const vector<DEPTH> &chromCov, const int &chromSize, const string &chrom) {
  
  int depth = 0; // initialize the depth
//  int offset = 0; // offset = (_eachBaseZeroBased)?0:1;
  for (int pos = 0; pos < chromSize; pos++) {
    
    depth += chromCov[pos].starts;
    // report the depth for this position.
    if (depth>0)
      std::cout << chrom << "\t" << pos << "\t" << depth<< std::endl;
    depth = depth - chromCov[pos].ends;
  }
}
