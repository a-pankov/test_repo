//
//  bedFile.h
//  myPileup
//
//  Created by Alex Pankov on 10/2/13.
//  Copyright (c) 2013 Alex Pankov. All rights reserved.
//

#ifndef myPileup_bedFile_h
#define myPileup_bedFile_h

typedef uint32_t CHRPOS;

struct BED {
  
// Regular BED fields
  string chrom;
  CHRPOS start;
  CHRPOS end;
  string name;
  string score;
  string strand;
  // all of the original fields in the record
  vector<string> fields;
  // indices of the "other" fields
  vector<uint16_t> other_idxs;
  // is this a zero length feature: i.e., start == end
  bool   zeroLength;
  
public:
  // constructors
  
  // Null
  BED()
  : chrom(""),
  start(0),
  end(0),
  name(""),
  score(""),
  strand(""),
  fields(),
  other_idxs(),
  zeroLength(false)
  {}
  
  // BED3
  BED(string chrom, CHRPOS start, CHRPOS end)
  : chrom(chrom),
  start(start),
  end(end),
  name(""),
  score(""),
  strand(""),
  fields(),
  other_idxs(),
  zeroLength(false)
  {}
  
  // BED4
  BED(string chrom, CHRPOS start, CHRPOS end, string strand)
  : chrom(chrom),
  start(start),
  end(end),
  name(""),
  score(""),
  strand(strand),
  fields(),
  other_idxs(),
  zeroLength(false)
  {}
  
  // BED6
  BED(string chrom, CHRPOS start, CHRPOS end, string name,
      string score, string strand)
  : chrom(chrom),
  start(start),
  end(end),
  name(name),
  score(score),
  strand(strand),
  fields(),
  other_idxs(),
  zeroLength(false)
  {}
  
  // BEDALL
  BED(string chrom, CHRPOS start, CHRPOS end, string name,
      string score, string strand, vector<string> fields,
      vector<uint16_t> other_idxs)
  : chrom(chrom),
  start(start),
  end(end),
  name(name),
  score(score),
  strand(strand),
  fields(fields),
  other_idxs(other_idxs),
  zeroLength(false)
  {}
  
  int size() {
    return end-start;
  }
  
}; // BED

typedef vector<BED>    bedVector;

struct DEPTH {
  uint32_t starts;
  uint32_t ends;
};



#endif
