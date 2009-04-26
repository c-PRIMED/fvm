#ifndef _MATRIXOPERATION_H_
#define _MATRIXOPERATION_H_



template <class T> 
class matrix{
 public:

  matrix() {};
  ~matrix() {};

//source matrix b
//invert matrix a
void Invert4x4(T b[4][4], T a[4][4])
  {
    int indxc[4], indxr[4], ipiv[4];
    int i, icol, irow, j, ir, ic;
    T big, dum, pivinv, temp, bb;
    ipiv[0] = -1;
    ipiv[1] = -1;
    ipiv[2] = -1;
    ipiv[3] = -1;
    a[0][0] = b[0][0];
    a[1][0] = b[1][0];
    a[2][0] = b[2][0];
    a[3][0] = b[3][0];
    a[0][1] = b[0][1];
    a[1][1] = b[1][1];
    a[2][1] = b[2][1];
    a[3][1] = b[3][1];
    a[0][2] = b[0][2];
    a[1][2] = b[1][2];
    a[2][2] = b[2][2];
    a[3][2] = b[3][2];
    a[0][3] = b[0][3];
    a[1][3] = b[1][3];
    a[2][3] = b[2][3];
    a[3][3] = b[3][3];
    for (i = 0; i < 4; i++) {
      big = 0.0f;
      for (j = 0; j < 4; j++) {
        if (ipiv[j] != 0) {
            if (ipiv[0] == -1) {
                if ((bb = (float) fabs(a[j][0])) > big) {
                    big = bb;
                    irow = j;
                    icol = 0;
		}
	    } else if (ipiv[0] > 0) {
	      return;
	    }
	    if (ipiv[1] == -1) {
	      if ((bb = (float) fabs((float) a[j][1])) > big) {
		big = bb;
		irow = j;
		icol = 1;
	      }
	    } else if (ipiv[1] > 0) {
	      return;
	    }
	    if (ipiv[2] == -1) {
    if ((bb = (float) fabs((float) a[j][2])) > big) {
        big = bb;
        irow = j;
        icol = 2;
    }
} else if (ipiv[2] > 0) {
    return;
}
if (ipiv[3] == -1) {
    if ((bb = (float) fabs((float) a[j][3])) > big) {
        big = bb;
        irow = j;
        icol = 3;
    }
} else if (ipiv[3] > 0) {
    return;
        }
    }
}
++(ipiv[icol]);
if (irow != icol) {
    temp = a[irow][0];
    a[irow][0] = a[icol][0];
    a[icol][0] = temp;
    temp = a[irow][1];
    a[irow][1] = a[icol][1];
    a[icol][1] = temp;
    temp = a[irow][2];
    a[irow][2] = a[icol][2];
    a[icol][2] = temp;
    temp = a[irow][3];
    a[irow][3] = a[icol][3];
    a[icol][3] = temp;
}
indxr[i] = irow;
indxc[i] = icol;
if (a[icol][icol] == 0.0) {
  throw CException("4x4 singular matrix!" );
    return;
}
pivinv = 1.0f / a[icol][icol];
a[icol][icol] = 1.0f;
a[icol][0] *= pivinv;
a[icol][1] *= pivinv;
a[icol][2] *= pivinv;
a[icol][3] *= pivinv;
    if (icol != 0) {
        dum = a[0][icol];
        a[0][icol] = 0.0f;
        a[0][0] -= a[icol][0] * dum;
        a[0][1] -= a[icol][1] * dum;
        a[0][2] -= a[icol][2] * dum;
        a[0][3] -= a[icol][3] * dum;
    }
    if (icol != 1) {
        dum = a[1][icol];
        a[1][icol] = 0.0f;
        a[1][0] -= a[icol][0] * dum;
        a[1][1] -= a[icol][1] * dum;
        a[1][2] -= a[icol][2] * dum;
        a[1][3] -= a[icol][3] * dum;
    }
    if (icol != 2) {
        dum = a[2][icol];
        a[2][icol] = 0.0f;
        a[2][0] -= a[icol][0] * dum;
        a[2][1] -= a[icol][1] * dum;
        a[2][2] -= a[icol][2] * dum;
        a[2][3] -= a[icol][3] * dum;
    }
    if (icol != 3) {
        dum = a[3][icol];
        a[3][icol] = 0.0f;
        a[3][0] -= a[icol][0] * dum;
        a[3][1] -= a[icol][1] * dum;
        a[3][2] -= a[icol][2] * dum;
        a[3][3] -= a[icol][3] * dum;
    }
}
if (indxr[3] != indxc[3]) {
    ir = indxr[3];
    ic = indxc[3];
    temp = a[0][ir];
    a[0][ir] = a[0][ic];
    a[0][ic] = temp;
    temp = a[1][ir];
    a[1][ir] = a[1][ic];
    a[1][ic] = temp;
    temp = a[2][ir];
    a[2][ir] = a[2][ic];
    a[2][ic] = temp;
    temp = a[3][ir];
    a[3][ir] = a[3][ic];
    a[3][ic] = temp;
}
if (indxr[2] != indxc[2]) {
    ir = indxr[2];
    ic = indxc[2];
    temp = a[0][ir];
    a[0][ir] = a[0][ic];
    a[0][ic] = temp;
    temp = a[1][ir];
    a[1][ir] = a[1][ic];
    a[1][ic] = temp;

      temp = a[2][ir];
      a[2][ir] = a[2][ic];
      a[2][ic] = temp;
      temp = a[3][ir];
      a[3][ir] = a[3][ic];
      a[3][ic] = temp;
  }
  if (indxr[1] != indxc[1]) {
      ir = indxr[1];
      ic = indxc[1];
      temp = a[0][ir];
      a[0][ir] = a[0][ic];
      a[0][ic] = temp;
      temp = a[1][ir];
      a[1][ir] = a[1][ic];
      a[1][ic] = temp;
      temp = a[2][ir];
      a[2][ir] = a[2][ic];
      a[2][ic] = temp;
      temp = a[3][ir];
      a[3][ir] = a[3][ic];
      a[3][ic] = temp;
  }
  if (indxr[0] != indxc[0]) {
      ir = indxr[0];
      ic = indxc[0];
      temp = a[0][ir];
      a[0][ir] = a[0][ic];
      a[0][ic] = temp;
      temp = a[1][ir];
      a[1][ir] = a[1][ic];
      a[1][ic] = temp;
      temp = a[2][ir];
      a[2][ir] = a[2][ic];
      a[2][ic] = temp;
      temp = a[3][ir];
      a[3][ir] = a[3][ic];
      a[3][ic] = temp;
  }
}


//source matrix b
//invert matrix a
void Invert3x3(T b[3][3], T a[3][3])
{
  const double det = b[0][0]*b[1][1]*b[2][2]+
    b[0][1]*b[1][2]*b[2][0]+
    b[0][2]*b[2][0]*b[2][1]-
    b[0][0]*b[1][2]*b[2][1]-
    b[0][1]*b[1][0]*b[2][2]-
    b[0][2]*b[1][1]*b[2][0];

  //  if(det == 0)
  //throw CException ("sigular matrix");


  a[0][0] = b[1][1]*b[2][2]-b[1][2]*b[2][1];
  
  a[0][1] = b[0][2]*b[2][1]-b[0][1]*b[2][2];

  a[0][2] = b[0][1]*b[1][2]-b[0][2]*b[1][1];
  
  a[1][0] = b[1][2]*b[2][0]-b[1][0]*b[2][2];

  a[1][1] = b[0][0]*b[2][2]-b[0][2]*b[2][0];
  
  a[1][2] = b[0][2]*b[1][0]-b[0][0]*b[1][2];
  
  a[2][0] = b[1][0]*b[2][1]-b[1][1]*b[2][0];
  
  a[2][1] = b[0][1]*b[2][0]-b[0][0]*b[2][1];

  a[2][2] = b[0][0]*b[1][1]-b[0][1]*b[1][0];

  for(int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      a[i][j] /= det;
    }
  }

}




};


#endif
