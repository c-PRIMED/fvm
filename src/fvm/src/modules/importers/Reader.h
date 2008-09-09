#ifndef _READER_H_
#define _READER_H_

#include <iostream>

#include "RLogInterface.h"

using namespace std;

class Reader
{
public:
  Reader(const string& fileName);
  virtual ~Reader();
  void resetFilePtr();
  string readLine();
protected:
  const string _fileName; 
  FILE *_fp;
};
  
#endif
