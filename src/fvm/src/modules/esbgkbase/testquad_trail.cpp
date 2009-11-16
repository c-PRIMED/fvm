#include <iostream>
#include <math.h>
#include "Quad_trial.h"
using namespace std;
int main(){
  quadrature myquad;
  int N1=8,N2=12,N3=13;
  int n_int=4,nphi_int=4;
  myquad.set_size(N1,N2,N3);
  // myquad.set_abscii_cartesian(1.0); //set cartesian type
  //myquad.set_abscii(2,1,1,4,4); //set spherical type all 3 in one function
 myquad.set_absci1(N1);myquad.set_absci2(1,n_int);myquad.set_absci3(1,nphi_int);
cout << myquad.getmember(4)<< '\n';
cout <<myquad.getmember(6)<< '\n';
  double density=0.5;
  double temperature=2.0;
  double x_velocity=0.1;
  double y_velocity=0.01;
  return 0;
}
