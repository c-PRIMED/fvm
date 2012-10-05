// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _READER_H_
#define _READER_H_

#include <iostream>
#include <stdio.h>

#include "RLogInterface.h"

using namespace std;

class Reader
{
public:
  Reader(const string& fileName);
  virtual ~Reader();
  void resetFilePtr();
  string readLine();
  void close();
  
protected:
  const string _fileName; 
  FILE *_fp;
};
  
#endif
