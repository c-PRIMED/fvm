#ifndef _MMREADER_H_
#define _MMREADER_H_

#include "Reader.h"
#include "MultiField.h"
#include "MultiFieldMatrix.h"
#include "CRMatrix.h"
#include "Field.h"

using namespace std;

class MMReader : public Reader
{
public:

  typedef CRMatrix<double,double,double> MatrixType;
  
  MMReader(const string& matrixFileName, const string& rhsFileName);
  virtual ~MMReader();

  DEFINE_TYPENAME("MMReader");

  shared_ptr<LinearSystem> getLS();

private:
  const string _rhsFileName;
  Field _field;
  shared_ptr<StorageSite> _site;
  shared_ptr<CRConnectivity> _cm;
  
  void readHeader(int& nRows, int& nCols,
                  int& nNonZeroes, bool& isSymmetric);

};

#endif
