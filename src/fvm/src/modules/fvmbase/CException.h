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
