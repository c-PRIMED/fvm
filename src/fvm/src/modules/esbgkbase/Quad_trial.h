
class quadrature {
  int N1, N2, N3, N123;                 //N123=N1*N2*N3
  double *absci1,*absci2,*absci3;       //abscissa in 3 dimensions
  double *wts1, *wts2, *wts3;           //weights in 3 dimensions
  double *cx, *cy, *cz, *dcxyz;        // set here to be used by ESBGK solver
  public:
  void set_size (int,int,int);          // sets memory space given N1,N2,N3
  void set_abscii_cartesian(double T2); //sets cartesian type velocity coordinates
  void set_absci1(int option_ur);       // sets quadrature for velocity magnitude(0,2,8,16,80,160)
  void set_absci2(int option_theta, int n_int); //sets angle theta (0 for uniform, 1 for 3/8th rule: 3*n_int points)
  void set_absci3(int option_phi, int nphi_int);//sets angle phi (0 for uniform, 1 for 3/8th rule: 3*nphi_int+1 points)
  void set_cxyz(int option_velgrid);           // sets cx, cy, cz and dcxyz for N1*N2*N3 directions
double getmember(int x);
  // int jmax () {return (N1*N2*N3);}
};
double quadrature::getmember(int x){
return absci1[x];
}
void quadrature::set_size (int a, int b, int c) {
  int j1;
  N1 = a;N2 = b;N3 = c;
  N123= a*b*c;
  absci1=new double[N1];absci2=new double[N2];absci3=new double[N3];
  wts1=new double[N1]; wts2=new double[N2]; wts3=new double[N3];
  cx=new double[N123];cy=new double[N123];cz=new double[N123];
  dcxyz=new double[N123];
}
void quadrature::set_abscii_cartesian(double T2){
  double clim = 5.5;
  double cxmin,cxmax,cymin,cymax,czmin,czmax;
  double dcx,dcy,dcz;
  int j1,j2,j3,j;
 
  //velocity space
  cxmin=-clim*sqrt(0.5*T2);cxmax=clim*sqrt(0.5*T2);
  cymin=-clim*sqrt(0.5*T2);cymax=clim*sqrt(0.5*T2);
  czmin=-clim*sqrt(0.5*T2);czmax=clim*sqrt(0.5*T2);
  dcx=(cxmax-cxmin)/(N1-1.0);
  dcy=(cymax-cymin)/(N2-1.0);  
  dcz=(czmax-czmin)/(N3-1.0);
  
  absci1[0]=cxmin;absci2[0]=cymin; absci3[0]=czmin;
  for  (j3=1;j3<N3;j3++){
    absci3[j3]=absci3[j3-1]+dcz;
    wts3[j3]=dcz;
  }
  for (j2=1;j2<N2;j2++){
    absci2[j2]=absci2[j2-1]+dcy;
    wts2[j2]=dcy;
  }
  for  (j1=1;j1<N1;j1++){
    absci1[j1]=absci1[j1-1]+dcx;
    wts1[j1]=dcx;
  }

}

//void quadrature::set_abscii(int option_ur, int option_theta, int option_phi int n_int, int nphi_int){
void quadrature::set_absci1(int option_ur){
  //Spherical type coordinates
  switch(option_ur){
    int j1;
  case 0:   // constant difference for Ur
    {double dh1=sqrt(3.0)*3.889/(N1); //neglect ur=.0; 
      for (j1=0;j1<N1;j1++){
	absci1[j1]=(j1+1)*dh1;
	wts1[j1]=dh1*pow(absci1[j1],2.0);
      }
    }break; 
  case 2:{
      absci1[0]=0.7539869079898871;
      absci1[1]=0.1734055298879163E+1;
      wts1[0]=0.2738413467266824; 
      wts1[1]=0.1692721159996965;
      for (j1=0;j1<N1;j1++){
	wts1[j1]=wts1[j1]*exp(pow(absci1[j1],2.0));
      }
    }break;
  case 4:            //4-gaussian for y^2*exp(-y^2)
    { absci1[0]=0.4238628193900528;
      absci1[1]=0.1014332104566760E+1;
      absci1[2]=0.1742437375162050E+1;
      absci1[3]=0.2639813333635586E+1;
      wts1[0]=0.7649092266787873E-1;
      wts1[1]= 0.2435439494642453;
      wts1[2]= 0.1162953035510695;
      wts1[3]= 0.6783287043185401E-2;
      for (j1=0;j1<N1;j1++){
	wts1[j1]=wts1[j1]*exp(pow(absci1[j1],2.0));
      }
    }break; 
  case 8:   //8-gaussian for y^2*exp(-y^2)
    { absci1[0]=0.1990000637984294;
      absci1[1]=0.5059526450205794;
      absci1[2]=0.9041682182040568;
      absci1[3]=0.1372615723971598E+1;
      absci1[4]=0.1900969572329702E+1;
      absci1[5]=0.2490479841967435E+1;
      absci1[6]=0.3158780677105240E+1;
      absci1[7]=0.3966720403265353E+1;
      
      wts1[0]=0.9599144336400067E-2;
      wts1[1]=0.7072944976303661E-1;
      wts1[2]=0.157366887003943;
      wts1[3]=0.1429322724003870;
      wts1[4]=0.5431444004253597E-1;
      wts1[5]=0.7835224153141577E-2;
      wts1[6]=0.3338952597020048E-3;
      wts1[7]=0.2149767232664775E-5;
      for(j1=0;j1<N1;j1++){
	wts1[j1]=wts1[j1]*exp(pow(absci1[j1],2.0));
      }
    }break;
  case 80://gaussian for y*exp(-y^2)
    { absci1[0]=0.1218127678061463;
      absci1[1]=0.3882449491473571;
      absci1[2]=0.7651497067658092;
      absci1[3]=0.1224690624761160E+1;
      absci1[4]=0.1751398297664409E+1;
      absci1[5]=0.2343383197810315E+1;
      absci1[6]=0.3016608849956826E+1;
      absci1[7]=0.3831371300820741E+1;
      wts1[0]=0.23978773117765308E-1;
      wts1[1]=0.1092506819189940;
      wts1[2]=0.1797622678433810;
      wts1[3]=0.1351751653621029;
      wts1[4]=0.4552181928573556E-1;
      wts1[5]=0.6064921853788935E-2;
      wts1[6]=0.2448536436477049E-3;
      wts1[7]=0.1516914696753451E-5;
      for (j1=0;j1<N1;j1++){
	wts1[j1]=wts1[j1]*exp(pow(absci1[j1],2.0))*absci1[j1];
      }
    }break;
  case 16: //weights & absci for y^2 exp(-y^2)
    { absci1[0]= 0.8174913389984520E-1;
      absci1[1]=0.2154761962759740;
      absci1[2]=0.4003530517087630;
      absci1[3]=0.6298538771405607;
      absci1[4]=0.8976124329697087;
      absci1[5]=0.1198149260240153E+1;
      absci1[6]=0.1527188184719962E+1;
      absci1[7]=0.1881745606015598E+1;
      absci1[8]=0.2260132964654088E+1;
      absci1[9]=0.2661980315279350E+1;
      absci1[10]=0.3088376381635592E+1;
      absci1[11]=0.3542256017930265E+1;
      absci1[12]=0.4029312272760483E+1;
      absci1[13]=0.4560203031431090E+1;
      absci1[14]=0.5156826768007481E+1;
      absci1[15]=0.58781144889155572E+1;
      
      wts1[0]=0.7050727473210895E-3;
      wts1[1]=0.7107111654073120E-2;
      wts1[2]=0.2844188515941899E-1;
      wts1[3]=0.6660235171398239E-1;
      wts1[4]=0.1025785712747278;
      wts1[5]=0.1077502032531791;
      wts1[6]=0.7747156370638879E-1;
      wts1[7]=0.3763106373385135E-1;
      wts1[8]=0.1204873635560290E-1;
      wts1[9]=0.2453208613776865E-2;
      wts1[10]=0.3020309847850189E-3;
      wts1[11]=0.2092121075871870E-4;
      wts1[12]=0.7314637349679360E-6;
      wts1[13]=0.1080646863902574E-7;
      wts1[14]=0.4828081616137754E-10;
      wts1[15]=0.2840126937112534E-13;
      for (j1=0;j1<N1;j1++){
	wts1[j1]=wts1[j1]*exp(pow(absci1[j1],2.0));
      }
    }break;
  case 160:    //Gaussian with N=16  y*exp(-y^2)
    { absci1[0]=0.477579959543737674E-1;
      absci1[1]=0.1575643611266753;
      absci1[2]=0.3236556568455920;
      absci1[3]=0.5391473546675038;
      absci1[4]=0.7970053979972014;
      absci1[5]=0.1090958307363892E+1;
      absci1[6]=0.14159759700714936E+1;
      absci1[7]=0.1768437030466615E+1;
      absci1[8]=0.2146149962010079E+1;
      absci1[9]=0.2548365652625752E+1;
      absci1[10]=0.2975896592510777E+1;
      absci1[11]=0.3431483868308089E+1;
      absci1[12]=0.3920694119664905E+1;
      absci1[13]=0.4454120573510955E+1;
      absci1[14]=0.5053674269642785E+1;
      absci1[15]=0.5778478847939104E+1;
      
      wts1[0]=0.3795307814831678E-2;
      wts1[1]=0.2136808301992049E-1;
      wts1[2]=0.5595857089379011E-1;
      wts1[3]=0.9587168277747507E-1;
      wts1[4]=0.1169082070371872;
      wts1[5]=0.1029363012162623;
      wts1[6]=0.6468246716393942E-1;
      wts1[7]=0.2831911613754905E-1;
      wts1[8]=0.8362647991652432E-2;
      wts1[9]=0.1597736202726321E-2;
      wts1[10]=0.1870134647150351E-3;
      wts1[11]=0.1243935496206526E-4;
      wts1[12]=0.4208466925294357E-6;
      wts1[13]=0.6051847030054333E-8;
      wts1[14]=0.2643406562982473E-10;
      wts1[15]=0.1524594098604790E-13;
      for (j1=0;j1<N1;j1++){
	wts1[j1]=wts1[j1]*exp(pow(absci1[j1],2.0))*absci1[j1];
      }
    }
  }  
  }

void quadrature::set_absci2(int option_theta, int n_int){
  switch(option_theta){
   
  case 0:   //constant difference for Theta
    { double pi=3.14159;
      double dh2=2*pi/(N2); //number of intervals=no. of ordinates
      for (int j2=0;j2<N2;j2++){
	absci2[j2]= dh2*j2;
	wts2[j2]=dh2;
      }
    }break;
  case 1:  //three-eights for theta
    { double pi=3.14159;
      double h=pi;                     //the range of quadrature from [-h,h];
      double dh=2*h/n_int;             //anglular increment
      for (int i=0;i<n_int;i++){
	absci2[3*i]=-h+dh/3.0*(3.0*i);       
	absci2[3*i+1]=-h+dh/3.0*(3.0*i+1);   
	absci2[3*i+2]=-h+dh/3.0*(3.0*i+2);   
	if (i == 0){
	  wts2[i]=dh/4.0;  //here 0 -pi and pi are the same number
	}
	else {
	  wts2[3*i]=dh/4.0;
	}
	wts2[3*i+1]=dh/8.0*3.0;
	wts2[3*i+2]=dh/8.0*3.0;
      }
      //absci2(N2)=h;wts2(N2)=dh/8.0;
    }break;
  }
  }
void quadrature::set_absci3(int option_phi, int nphi_int){
  switch(option_phi){
   
  case 0: //constant difference for Phi
    {
double pi=3.14159;
      double dh3=pi/(N3-1); //number of intervals < number of ordinates
      for (int j3=0;j3<N3;j3++){
	absci3[j3]= dh3*j3;
	wts3[j3]=dh3*sin(absci3[j3]);
      }
    }break;
  case 1: //three-eights for phi
    { double pi=3.14159;
      double dh_phi=pi/nphi_int;  //anglular increment
      for (int i=0;i<nphi_int;i++){
	absci3[3*i]=dh_phi/3.0*(3.0*i);      
	absci3[3*i+1]=dh_phi/3.0*(3.0*i+1);  
	absci3[3*i+2]=dh_phi/3.0*(3.0*i+2);  
	if (i == 0){  
	  wts3[i]=dh_phi/8.0*sin(absci3[i]);
	}
	else{
	  wts3[3*i]=dh_phi/4.0*sin(absci3[3*i]);
	}      
	wts3[3*i+1]=dh_phi/8.0*3.0*sin(absci3[3*i+1]);
	wts3[3*i+2]=dh_phi/8.0*3.0*sin(absci3[3*i+2]);
      }
      absci3[N3-1]=pi;
      wts3[N3-1]=dh_phi/8.0*sin(absci3[N3-1]);
    }break;
  }
}
void quadrature::set_cxyz(int option_velgrid){
  
  switch(option_velgrid){
  case 1:{ //cartesian_type
int j,j1,j2,j3;
 j=0;
  for(j3=0;j3<N3;j3++){
    for (j2=0;j2<N2;j2++){
      for (j1=0;j1<N1;j1++){
	cx[j]=absci1[j1];	
	cy[j]=absci2[j2];
	cz[j]=absci3[j3];
	dcxyz[j]=wts1[j1]*wts2[j2]*wts3[j3];
	j=j+1;
      }}}
  }break;
  case 2:{  //spherical type
    int j,j1,j2,j3;
    j=0;
  for (j1=0;j1<N1;j1++){ //Ur
    for (j2=0;j2<N2;j2++){ //theta
      for (j3=0;j3<N3;j3++){ //phi
	cx[j]=absci1[j1]*cos(absci2[j2])*sin(absci3[j3]); //cx=Ur*cos(theta)*sin(phi)
	cy[j]=absci1[j1]*sin(absci2[j2])*sin(absci3[j3]); //cy=Ur*sin(theta)*sin(phi)
	cz[j]=absci1[j1]*cos(absci3[j3]);                 //Ur*cos(phi)
	dcxyz[j]=wts1[j1]*wts2[j2]*wts3[j3];
	j=j+1;       
      }}}
  }break;
  }
}
/*
int quadrature::distribution_function(int ix,int iy, int j){
  double pi=3.14159;
  for (int j=0;j<N123;j++){
    // f[ix][iy][j]=Rho[ix][iy]/pow((pi*Temp[ix][iy]),1.5)* exp(-((pow(cx[j]-xVel[ix][iy]),2.0)+pow((cy[j]-yVel[ix][iy]),2.0)+pow(cz[j],2.0))/Temp[ix][iy]);
    //f[1][1][j]=1.0;
    f[j][1][1]=density/pow((pi*temperature),1.5)* exp(-(pow((cx[j]-x_velocity),2.0)+pow((cy[j]-y_velocity),2.0)+pow(cz[j],2.0))/temperature);
  }
}

int main_old()
{
  int N1,N2,N3,N123; 
  int n_int,nphi_int;
  int option_velgrid =2, option_ur=8,option_theta=1,option_phi=1,solver=2;
  N123=N1*N2*N3;
  double f[N123][1][1];
  // distribution_function(N123,cx,cy,cz,density,temperature,x_velocity,y_velocity,f);
  double Rho[1][1],Temp[1][1],xVel[1][1],yVel[1][1],xTemp[1][1];
  
  double density=0.5;
  double temperature=2.0;
  double x_velocity=0.1;
  double y_velocity=0.01;
Rho[1][1]=density;
  Temp[1][1]=temperature;
  xVel[1][1]=x_velocity;
  yVel[1][1]=y_velocity;
  // Macroparameters(N123,cx,cy,cz,Rho,Temp,xVel,yVel,xTemp,f,dcxyz);
    cout << "Rho "<< Rho[1][1]<<'\n';
cout << "Temp "<< Temp[1][1]<<'\n';
  cout << "xVel "<< xVel[1][1]<<'\n'; 
  cout << "yVel "<< yVel[1][1]<<'\n';
  cout << "xTemp " << xTemp[1][1]<< '\n';
  
return 0;
  
}
*/



