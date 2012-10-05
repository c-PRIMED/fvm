// This file os part of FVM
// Copyright (c) 2012 FVM Authors
// See LICENSE file for terms.

#ifndef _CEXCEPTION_H_
#define _CEXCEPTION_H_

#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;


class CException : public runtime_error
{
public:
  CException(const string& what);
  virtual ~CException() throw() {}
};

#endif
