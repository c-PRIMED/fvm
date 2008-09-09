#include "Reader.h"


Reader::Reader(const string& fileName) :
  _fileName(fileName),
  _fp(0)
{
  resetFilePtr();
}

Reader::~Reader()
{}



void
Reader::resetFilePtr()
{
  if (_fp)
    fclose(_fp);
  _fp = fopen(_fileName.c_str(),"rb");
}

string
Reader::readLine()
{
  char buf[256];
  fgets(buf,256,_fp);
  string s(buf);
  return s;
}
    
