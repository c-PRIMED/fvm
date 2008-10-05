#ifndef _NUMTYPE_H_
#define _NUMTYPE_H_

#include <string>
#include <math.h>
#include <stdlib.h>

using namespace std;

template<class T>
struct NumTypeTraits : public T
{};

template<>
struct   NumTypeTraits<bool>
{
  typedef bool T_Scalar;
  typedef bool T_BuiltIn;

  static string getTypeName() {return "bool";}
  static int getDimension() {return 0;}
  static void getShape(int *shp) {}
  static int getDataSize()  {return  sizeof(bool);}
  static bool getZero() {return false;}
  static bool getUnity() {return true;}
  static bool getNegativeUnity() {return false;}
  
  static bool sqrt(const bool x) {return false;}

  static void write(FILE* fp, const bool x) {fprintf(fp,"%d",(int)x);}
  static void accumulateOneNorm(bool& sum, const bool& v) {}
  static void accumulateDotProduct(bool& sum, const bool& v0, const bool& v1)
  {}

  static void safeDivide(bool& x, const bool& y) {}
};

template<>
struct   NumTypeTraits<int>
{
  typedef int T_Scalar;
  typedef int T_BuiltIn;

  static string getTypeName() {return "int";}
  static int getDimension() {return 0;}
  static void getShape(int *shp) {}
  static int getDataSize()  {return  sizeof(int);}
  static int getZero() {return 0;}
  static int getUnity() {return 1;}
  static int getNegativeUnity() {return -1;}
  static int sqrt(const int x) {return (int)::sqrt((double)x);}
  static void write(FILE* fp, const int x) {fprintf(fp,"%d",x);}
  static void accumulateOneNorm(int& sum, const int& v) { sum += abs(v);}
  static void accumulateDotProduct(int& sum, const int& v0, const int& v1)
  { sum += v0*v1;}

  static void safeDivide(int& x, const int& y) {if (y!=0) x/=y;}
};

template<>
struct   NumTypeTraits<double>
{
  typedef double T_Scalar;
  typedef double T_BuiltIn;
  static string getTypeName() {return "double";}
  static int getDimension() {return 0;}
  static void getShape(int *shp) {}
  static int getDataSize()  {return  sizeof(double);}
  static double getZero() {return 0;}
  static double getUnity() {return 1.;}
  static double getNegativeUnity() {return -1.;}
  static double sqrt(const double x) {return ::sqrt(x);}
  static void write(FILE* fp, const double x) {fprintf(fp,"%lf",x);}
  static double doubleMeasure(const double x) {return fabs(x);}

  static void accumulateOneNorm(double& sum, const double& v) { sum += fabs(v);}
  static void accumulateDotProduct(double& sum, const double& v0, const double& v1)
  { sum += v0*v1;}

  static void safeDivide(double& x, const double& y) {if (y!=0) x/=y;}

};

template<>
struct NumTypeTraits<float>
{ 
  typedef float T_Scalar;
  typedef float T_BuiltIn;

  static string getTypeName() {return "float";}
  static int getDimension() {return 0;}
  static void getShape(int *shp) {}
  static int getDataSize()  {return  sizeof(double);}
  static float getZero() {return 0;}
  static float getUnity() {return 1.;}
  static float getNegativeUnity() {return -1.;}
  static float sqrt(const float x) {return ::sqrt(x);}
  static void write(FILE* fp, const float x) {fprintf(fp,"%f",x);}
  static double doubleMeasure(const float x) {return fabs(double(x));}
  static void accumulateOneNorm(float& sum, const float& v) { sum += fabs(v);}
  static void accumulateDotProduct(float& sum, const float& v0, const float& v1)
  { sum += v0*v1;}
  static void safeDivide(float& x, const float& y) {if (y!=0) x/=y;}
};

#endif
