#ifndef _SCHEMEREADER_H_
#define _SCHEMEREADER_H_

#include <iostream>
#include "Reader.h"

using namespace std;

class SchemeReader : public Reader
{
public:
  SchemeReader(const string& fileName);
  virtual ~SchemeReader();

protected:
  int getNextSection();
  void closeSection();
  int closeSectionBinary(const int currentId);

  int readInt(const bool isBinary);
  void skipInt(const int n, const bool isBinary);
  
  char getNextChar();
  int moveToListOpen();
  void moveToListClose();
  void moveToListCloseBinary();
  void readHeader(int& i1, int& i2, int& i3, int& i4, int& i5);

  int readListLength();
  void readList(char *buffer);
};
  
#endif
