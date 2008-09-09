#include "SchemeReader.h"

SchemeReader::SchemeReader(const string& fileName) :
  Reader(fileName)
{}

SchemeReader::~SchemeReader()
{}

char
SchemeReader::getNextChar()
{
  char c;

  while ((c = getc(_fp)) != EOF)
    if (isprint(c) && !isspace(c))
      break;
  return(c);
}


/* move _fp just past next opening paren */
int
SchemeReader::moveToListOpen()
{
  char c;

  while ((c = getNextChar()))
    if ((char) c == '(')
      return 0;
    else if (c == EOF) 
      return EOF;

  return EOF;
}

/* Move _fp just past closing paren of current list */
void
SchemeReader::moveToListClose()
{
  char c;

  while ((c = getc(_fp)) != EOF)
    switch (c)
      {
      case '(':
	moveToListClose();
	break;
      case ')':
	return;
      }
}

int
SchemeReader::readListLength()
{
  char c;
  int parenLevel=0;
  int buffSize=0;

  
  while ((c = getc(_fp)) != EOF)
  {
      
      if (c == '(')
        parenLevel++;

      if (parenLevel > 0)
      {
          buffSize++;
      }
      
      if (c == ')')
        parenLevel--;
      
      if (parenLevel == 0  && buffSize > 0)
        return buffSize;
  }

  throw CException("EOF reached while reading list");
}

void
SchemeReader::readList(char *buffer)
{
  char c;
  int parenLevel=0;
  int buffSize=0;
  
  while ((c = getc(_fp)) != EOF)
  {
      if (c == '(')
        parenLevel++;

      if (parenLevel > 0)
      {
          buffer[buffSize] = c;
          buffSize++;
      }
      
      if (c == ')')
        parenLevel--;
      
      if (parenLevel == 0  && buffSize > 0)
        return;
  }

  throw CException("EOF reached while reading list");
}

void
SchemeReader::moveToListCloseBinary()
{
  char c;
  int id;
  moveToListOpen();
  do
    {
      while ((c = getc(_fp)) != '\n')
	if (EOF == c)
	  return;
    }
  while (1 != fscanf(_fp,"End of Binary Section %6d",&id));
}


int
SchemeReader::getNextSection()
{
  if (moveToListOpen() == EOF)
    return EOF;
  
  int id;
  if (fscanf(_fp,"%d", &id) != 1)
  {
      cerr << "error reading id "<< endl;
  }

  return id;
}


void
SchemeReader::closeSection()
{
  moveToListClose();
}

int
SchemeReader::closeSectionBinary(const int currentId)
{
  //  char c;
  int id;
  //  moveToListOpen();
  //moveToListClose();
  //moveToListOpen();
  do
  {
      getc(_fp);
//       while ((c = getc(_fp)) != '\n')
// 	if (EOF == c)
//         {
//             cerr << "error closing binary section: expected " << currentId
//                  << " , found EOF " << endl;

//             return EOF;
//         }
    }
  while (1 != fscanf(_fp,"End of Binary Section %d",&id));

  if (currentId != id)
    cerr << "error closing binary section: expected " << currentId
         << " , found " << id << endl;
  return id;
}

void
SchemeReader::readHeader(int& i1, int& i2, int& i3, int& i4, int& i5)
{
  moveToListOpen();

  if (fscanf(_fp, "%x%x%x%x",&i1,&i2,&i3,&i4) != 4)
  {
      cerr << "error reading header" << endl;
      //      readerError("Error reading the mesh header", DBG_HERE);
  }

  /* read shape, not always availabe */
  if (fscanf(_fp,"%x",&i5) != 1)
    i5 = -1;
  
  moveToListClose();
}

int
SchemeReader::readInt(const bool isBinary)
{
  int i;
  if (isBinary ?
      (1 != fread(&i, sizeof(int), 1, _fp)) :
      (1 != fscanf(_fp,"%x",&i)))
    cerr << "Error reading int " << endl;
  return i;
}

void
SchemeReader::skipInt(const int count, const bool isBinary)
{
  int i;
  for (int n=0; n<count; n++)
  {
      if (isBinary ?
          (1 != fread(&i, sizeof(int), 1, _fp)) :
          (1 != fscanf(_fp,"%x",&i)))
        cerr << "Error reading int " << endl;
  }
}
