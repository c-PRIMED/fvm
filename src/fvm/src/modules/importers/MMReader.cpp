#include "MMReader.h"
#include "Array.h"
#include "CRConnectivity.h"
#include "LinearSystem.h"


MMReader::MMReader(const string& matrixFileName, const string& rhsFileName) :
  Reader(matrixFileName),
  _rhsFileName(rhsFileName),
  _field("test"),
  _site(),
  _cm()
{}

MMReader::~MMReader()
{}

void
MMReader::readHeader(int& nRows, int& nCols,
                     int& nNonZeroes, bool& isSymmetric)
{
  if (_fp==0)
    throw CException("cannot open matrix file ");

  char buf[32];

  fscanf(_fp,"%s",buf);

  string h1(buf);
  if (h1 == "%%")
  {
      fscanf(_fp,"%s",buf);
      string h2(buf);
      if (h2 != "MatrixMarket")
        throw CException("not a MatrixMarket file");
  }
  else if (h1 != "%%MatrixMarket")
    throw CException("not a MatrixMarket file");    
  
  fscanf(_fp,"%s",buf);
  string matrixType(buf);

  fscanf(_fp,"%s",buf);
  string coordinate(buf);
  
  fscanf(_fp,"%s",buf);
  string ftype(buf);

  fscanf(_fp,"%s",buf);
  string symm(buf);
  
  isSymmetric = false;
  

  if (matrixType != "matrix")
    throw CException("not a MatrixMarket file");

  if (coordinate != "coordinate")
    throw CException("not a sparse matrix");

  if (ftype != "real")
    throw CException("not a real matrix");

  if (symm == "symmetric")
    isSymmetric = true;
  else if (symm != "general")
    throw CException("not symmetric or general matrix");

  fscanf(_fp, "%d %d %d", &nRows, &nCols, &nNonZeroes);

  if (nRows != nCols)
    throw CException("not a square matrix");
}

shared_ptr<LinearSystem>
MMReader::getLS()
{
  int nRows, nCols, nNonZeroes;
  bool isSymmetric;

  readHeader(nRows,nCols,nNonZeroes,isSymmetric);
  
  cout << " nRow " << nRows << " x " << nCols
       << " matrix with " << nNonZeroes << " entries "
       << endl;

  cout << "reading sizes" << endl;

  _site = shared_ptr<StorageSite>(new StorageSite(nRows,0));

  _cm = shared_ptr<CRConnectivity>(new CRConnectivity(*_site,*_site));
  
  
  _cm->initCount();

  for(int nr=0; nr<nNonZeroes; nr++)
  {
      int i, j;
      double c;
      fscanf(_fp,"%d %d %lf",&i,&j,&c);
      i-=1;
      j-=1;
      if (i!=j)
      {
          _cm->addCount(i,1);
          if (isSymmetric)
            _cm->addCount(j,1);
      }
  }
  
  _cm->finishCount();

  shared_ptr<MatrixType> m(new MatrixType(*_cm));

  MultiField::ArrayIndex rowI(&_field,_site.get());

  shared_ptr<Array<double> > xPtr(new Array<double>(nRows));

  xPtr->zero();

  shared_ptr<LinearSystem> ls(new LinearSystem());
  ls->getX().addArray(rowI,xPtr);
  ls->getMatrix().addMatrix(rowI,rowI,m);

  ls->initAssembly();

  Array<double>& diag = m->getDiag();
  Array<double>& offDiag = m->getOffDiag();

  // rewind file 
  cout << "rewinding and reading coeffs" << endl;

  resetFilePtr();
  
  readHeader(nRows,nCols,nNonZeroes,isSymmetric);

  for(int nr=0; nr<nNonZeroes; nr++)
  {
      int i, j;
      double c;
      fscanf(_fp,"%d %d %lf",&i,&j,&c);
      i-=1;
      j-=1;
      if (i!=j)
      {
          int pos = _cm->add(i,j);
          offDiag[pos] = c;
          
          if (isSymmetric)
          {
              pos = _cm->add(j,i);
              offDiag[pos] = c;
          }
      }
      else
        diag[i] = c;
  }
  
  _cm->finishAdd();
  cout << "finished reading " << endl;
  close();


  //ls->getB().addArray(rowI,bPtr);
  
                     
  Array<double>& b = dynamic_cast<Array<double>&>(ls->getB()[rowI]);
  
  FILE *bFile = fopen(_rhsFileName.c_str(), "rb");

  for(int i=0; i<nRows; i++)
  {
      double r;
      fscanf(bFile,"%lf",&r);
      b[i]=-r;
  }
  
  ls->initSolve();
  return ls;
}
