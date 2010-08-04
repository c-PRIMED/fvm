#ifndef _MATRIXOPERATION_H_
#define _MATRIXOPERATION_H_



template <class T> 
class matrix{
 public:

  matrix() {};
  ~matrix() {};

//source matrix b
//invert matrix a
void Inverse4x4(T b[4][4], T a[4][4])
  {
    int indxc[4], indxr[4], ipiv[4];
    int i, icol, irow, j, ir, ic;
    T big, dum, pivinv, temp, bb;
    ipiv[0] = -1;
    ipiv[1] = -1;
    ipiv[2] = -1;
    ipiv[3] = -1;

    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
	a[i][j]=b[i][j];
      }
    }
    
    for (i = 0; i < 4; i++) {
      big = 0.0;
      for (j = 0; j < 4; j++) {
        if (ipiv[j] != 0) {
            if (ipiv[0] == -1) {
                if ((bb = fabs(a[j][0])) > big) {

                    big = bb;
                    irow = j;
                    icol = 0;
		}
	    } else if (ipiv[0] > 0) {
	      return;
	    }
	    if (ipiv[1] == -1) {

	      if ((bb =  fabs( a[j][1])) > big) {

		big = bb;
		irow = j;
		icol = 1;
	      }
	    } else if (ipiv[1] > 0) {
	      return;
	    }
	    if (ipiv[2] == -1) {

    if ((bb =  fabs( a[j][2])) > big) {

        big = bb;
        irow = j;
        icol = 2;
    }
} else if (ipiv[2] > 0) {
    return;
}
if (ipiv[3] == -1) {

    if ((bb =  fabs(a[j][3])) > big) {

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
pivinv = 1.0 / a[icol][icol];
a[icol][icol] = 1.0;
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
//inversed matrix a
void Inverse3x3(T b[3][3], T a[3][3])
{
  const T det = b[0][0]*b[1][1]*b[2][2]+
    b[0][1]*b[1][2]*b[2][0]+
    b[0][2]*b[2][0]*b[2][1]-
    b[0][0]*b[1][2]*b[2][1]-
    b[0][1]*b[1][0]*b[2][2]-
    b[0][2]*b[1][1]*b[2][0];

  if(det == 0)
    throw CException ("sigular matrix");


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


//source matrix b
//inversed matrix a
void Inverse6x6(T b[6][6], T a[6][6])
 {
   const int size = 6; 
   const double det = detMatrix6x6(b, size);
   if(det == 0)
     {
       for (int i=0; i<size; i++){
	 for (int j=0; j<size; j++){
	   cout<<b[i][j]<<"  ";
	 }
	 cout<<endl;
       }
       throw CException ("sigular matrix");
     }
   T tmp[6][6], fac[6][6], trans[6][6];
   int p, q, m, n, i, j;
   
   for(q=0; q<size; q++)
     {
       for(p=0; p<size; p++)
	 {
	   m=0;
	   n=0;
	   for(i=0; i<size; i++)
	     {
	       for(j=0; j<size; j++)
		 {
		   tmp[i][j]=0;
		   if(i!=q&&j!=p)
		     {
		       tmp[m][n]=b[i][j];
		       if(n<(size-2))
			 n++;
		       else
			 {
			   n=0;
			   m++;
			 }
		     }
		 }
	     }
	   fac[q][p]=pow(-1,q+p)*detMatrix6x6(tmp,size-1);
	 }
     }
  
   for(i=0; i<size; i++)
     {
       for(j=0; j<size; j++)
	 {
	   trans[i][j]=fac[j][i];
	 }
     }

   for(i=0; i<size; i++)
     {
       for(j=0; j<size; j++)
	 {
	   a[i][j]=trans[i][j]/det;
	 }
     }
 }


//find the determinant of a matrix
double detMatrix6x6 (T matrix[6][6], int size)
{
   if (size > 6)
     throw CException("matrix size is too large!" );
   
   double s = 1.0;
   double det = 0.0;
   double b[6][6];
   int m,n;
   if(size==1){
     det = matrix[0][0];
     return det;
   }
   else{
     for(int k=0; k<size; k++){
       m = 0;
       n = 0;
       for(int i=0; i<size; i++){
	 for(int j=0; j<size; j++){
	   b[i][j]=0.0;
	   if(i!=0&&j!=k){
	     b[m][n]=matrix[i][j];
	     if(n<(size-2))
	       n++;
	     else{
	       n=0;
	       m++;
	     }
	   }
	 }
       }
       det = det + s*(matrix[0][k]*detMatrix6x6(b, size-1));
       s=-1*s;
     }
   }
   return det;
}


//source matrix b
//inversed matrix a
void Inverse10x10(T b[10][10], T a[10][10])
 {
   const int size = 10; 
   const double det = detMatrix10x10(b, size);
   if(det == 0)
     {
       for (int i=0; i<size; i++){
	 for (int j=0; j<size; j++){
	   cout<<b[i][j]<<"  ";
	 }
	 cout<<endl;
       }
       throw CException ("sigular matrix");
     }
   T tmp[10][10], fac[10][10], trans[10][10];
   int p, q, m, n, i, j;
   
   for(q=0; q<size; q++)
     {
       for(p=0; p<size; p++)
	 {
	   m=0;
	   n=0;
	   for(i=0; i<size; i++)
	     {
	       for(j=0; j<size; j++)
		 {
		   tmp[i][j]=0;
		   if(i!=q&&j!=p)
		     {
		       tmp[m][n]=b[i][j];
		       if(n<(size-2))
			 n++;
		       else
			 {
			   n=0;
			   m++;
			 }
		     }
		 }
	     }
	   fac[q][p]=pow(-1,q+p)*detMatrix10x10(tmp,size-1);
	 }
     }
  
   for(i=0; i<size; i++)
     {
       for(j=0; j<size; j++)
	 {
	   trans[i][j]=fac[j][i];
	 }
     }

   for(i=0; i<size; i++)
     {
       for(j=0; j<size; j++)
	 {
	   a[i][j]=trans[i][j]/det;
	 }
     }
 }


//find the determinant of a matrix
double detMatrix10x10 (T matrix[10][10], int size)
{
   if (size > 10)
     throw CException("matrix size is too large!" );
   
   double s = 1.0;
   double det = 0.0;
   double b[10][10];
   int m,n;
   if(size==1){
     det = matrix[0][0];
     return det;
   }
   else{
     for(int k=0; k<size; k++){
       m = 0;
       n = 0;
       for(int i=0; i<size; i++){
	 for(int j=0; j<size; j++){
	   b[i][j]=0.0;
	   if(i!=0&&j!=k){
	     b[m][n]=matrix[i][j];
	     if(n<(size-2))
	       n++;
	     else{
	       n=0;
	       m++;
	     }
	   }
	 }
       }
       det = det + s*(matrix[0][k]*detMatrix10x10(b, size-1));
       s=-1*s;
     }
   }
   return det;
}

//find the determinant of a matrix
double detMatrix4x4 (T matrix[4][4], int size)
{
   if (size > 4)
     throw CException("matrix size is too large!" );
   
   double s = 1.0;
   double det = 0.0;
   double b[4][4];
   int m,n;
   if(size==1){
     det = matrix[0][0];
     return det;
   }
   else{
     for(int k=0; k<size; k++){
       m = 0;
       n = 0;
       for(int i=0; i<size; i++){
	 for(int j=0; j<size; j++){
	   b[i][j]=0.0;
	   if(i!=0&&j!=k){
	     b[m][n]=matrix[i][j];
	     if(n<(size-2))
	       n++;
	     else{
	       n=0;
	       m++;
	     }
	   }
	 }
       }
       det = det + s*(matrix[0][k]*detMatrix4x4(b, size-1));
       s=-1*s;
     }
   }
   return det;
}


void InverseGauss (T a[10][10], T x[10][10])
{
  const int size = 10; 
  int indx[10];
  int i, j, k, itmp;
  double c1, pi, pi1, pj;
  T b[10][10], c[10];

  for (i = 0; i < size; i++)  {
    for (j = 0; j < size; j++)   {
      b[i][j] = 0;
    }
  }
  for (i = 0; i < size; i++)  {
    b[i][i] = 1;
  }
  /* Initialize the index */
  for (i = 0; i < size; i++)
  {
    indx[i] = i;
  }
  /* Find the rescaling factors, one from each row */ 
  for (i = 0; i < size; i++)
  {
    c1 = 0;
    for (j = 0; j < size; j++)
    {
      if (fabs(a[i][j]) > c1) c1 = fabs(a[i][j]);
    }
    c[i] = c1;
  }
  /* Search the pivoting (largest) element from each column */ 
  for (j = 0; j < size-1; j++)
  {
    pi1 = 0;
    for (i = j; i < size; i++)
    {
      pi = fabs(a[indx[i]][j])/c[indx[i]];
      if (pi > pi1)
      {
        pi1 = pi;
        k = i;
      }
    }
    /* Interchange the rows via indx[] to record pivoting order */
    itmp = indx[j];
    indx[j] = indx[k];
    indx[k] = itmp;
    for (i = j+1; i < size; i++)
    {
      pj = a[indx[i]][j]/a[indx[j]][j];
      /* Record pivoting ratios below the diagonal */
      a[indx[i]][j] = pj;
      /* Modify other elements accordingly */
      for (k = j+1; k < size; k++)
      {
        a[indx[i]][k] = a[indx[i]][k]-pj*a[indx[j]][k];
      }
    }
  }
  for (i = 0; i < size-1; i++)
  {
    for (j = i+1; j < size; j++)
    {
      for (k = 0; k < size; k++)
      {
        b[indx[j]][k] = b[indx[j]][k]-a[indx[j]][i]*b[indx[i]][k];
      }
    }
  }

  for (i = 0; i < size; i++)
  {
    x[size-1][i] = b[indx[size-1]][i]/a[indx[size-1]][size-1];
    for (j = size-2; j >= 0; j = j-1)
    {
      x[j][i] = b[indx[j]][i];
      for (k = j+1; k < size; k++)
      {
        x[j][i] = x[j][i]-a[indx[j]][k]*x[k][i];
      }
      x[j][i] = x[j][i]/a[indx[j]][j];
    }
  }
}





};





#endif

