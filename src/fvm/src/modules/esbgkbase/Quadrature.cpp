#include <math.h>
 
void velocity_quadrature(){

int j3,j2,j1,j;
int N1,N2,N3;
 double cxmin,cymin,czmin,cxmax,cymax,czmax;
  double clim=5.5;
if(option_velgrid == 1){   //Cartesian type velocity discretization
   
N1=14;N2=14;N3=14;
    //velocity space(cx(j3,j2,j1),cy(j3,j2,j1),cz(j3,j2,j1))
    cxmin=-clim*sqrt(0.5*T2);cxmax=clim*sqrt(0.5*T2);
    cymin=-clim*sqrt(0.5*T2);cymax=clim*sqrt(0.5*T2);
    czmin=-clim*sqrt(0.5*T2);czmax=clim*sqrt(0.5*T2);
    
    dcx=(cxmax-cxmin)/(N1-1.0);
    dcy=(cymax-cymin)/(N2-1.0);  //old
    dcz=(czmax-czmin)/(N3-1.0);
    
    cz0[0]=czmin;cy0[0]=cymin;cx0[0]=cxmin;
    for  (j3=1;j3<N3;j3++){
        cz0[j3]=cz0[j3-1]+dcz;
    }
    for (j2=1;j2<N2;j2++){
        cy0[j2]=cy0[j2-1]+dcy;
    }
    for  (j1=1;j1<N1;j1++){
        cx0[j1]=cx0[j1-1]+dcx;
    }
    //cz0(N3-1)=czmax;cy0(N2-1)=cymax;cx0(N1-1)=cxmax ;
   
    //volume of cell in velocity space
    j=0;
    for(j3=0;j3<N3;j3++){
      for (j2=0;j2<N2;j2++){
	for (j1=1;j1<N1;j1++){
	  cx[j]=cx0[j1];
	  cy[j]=cy0[j2];
	  cz[j]=cz0[j3];
	  dcxyz[j]=dcx[j1]*dcy[j2]*dcz[j3];
	  j=j+1;
	}
      }
    }
    // delete cx0,cy0,cz0,dcx,dcy,dcz,cxmin,cymin,czmin,cxmax,cymax,czmax,clim
 }
 else if (option_velgrid == 2){
   
   //Spherical type coordinates
   if(option_ur == 0){
     // constant difference for Ur
     dh1=sqrt(3.0)*3.889/(N1-1);  //dh1=sqrt(3.0)*cxmax *sqrt(0.5*T2)/(N1-1) //(cxmax^2+cymax^2+czmax^2)
     for (j1=0;j1<N1;j1++){
       absci1[j1]=j1*dh1;
       wts1[j1]=dh1*pow(absci1[j1],2.0);
     }
   }
   else if(option_ur==2){
     //dh1=1.0;
     absci1[]={0.7539869079898871,0.1734055298879163E+1};
     wts1[]={0.2738413467266824; 0.1692721159996965};
     for (j1=0;j1<N1;j1++){
       //absci1(j1)=absci1(j1)*sqrt(T2)
       wts1[j1]=wts1[j1]*exp(pow(absci1[j1],2.0));
     }
   }
   else if(option_ur == 4){               //4-gaussian for y^2*exp(-y^2)
     //dh1=1.0;                                //dummy used in dcxyz
     absci1[]={0.4238628193900528,
	     0.1014332104566760E+1,
	     0.1742437375162050E+1,
	     0.2639813333635586E+1};
     wts1[]={0.7649092266787873E-1,
	   0.2435439494642453,
	   0.1162953035510695,
	   0.6783287043185401E-2};
     for (j1=0;j1<N1;j1++){
       //absci1(j1)=absci1(j1)*sqrt(T2)
       wts1[j1]=wts1[j1]*exp(pow(absci1[j1],2.0));
     } 
   }
   else if(option_ur == 8){                //8-gaussian for y^2*exp(-y^2)
     //dh1=1.0;                                  //dummy used in dcxyz
     absci1[]={0.1990000637984294,
	     0.5059526450205794,
	     0.9041682182040568,
	     0.1372615723971598E+1,
	     0.1900969572329702E+1,
	     0.2490479841967435E+1,
	     0.3158780677105240E+1,
	     0.3966720403265353E+1};
     
     wts1[]={0.9599144336400067E-2,
	   0.7072944976303661E-1,
	   0.157366887003943 ,
	   0.1429322724003870,
	   0.5431444004253597E-1,
	   0.7835224153141577E-2,
	   0.3338952597020048E-3,
	   0.2149767232664775E-5};
     for(j1=0;j1<N1;j1++){
       wts1[j1]=wts1[j1]*exp(pow(absci1[j1],2.0));
     }
   }
   else if(option_ur == 80){                   //gaussian for y*exp(-y^2)
     absci1[0]=0.1218127678061463;
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
     
   }
   else if(option_ur == 16){   //weights & absci for y^2 exp(-y^2) from Shizgal
     
     absci1[]={0.8174913389984520E-1,
	       0.2154761962759740,
	       0.4003530517087630,
	       0.6298538771405607,
	       0.8976124329697087,
	       0.1198149260240153E+1,
	       0.1527188184719962E+1,
	       0.1881745606015598E+1,
	       0.2260132964654088E+1,
	       0.2661980315279350E+1,
	       0.3088376381635592E+1,
	       0.3542256017930265E+1,
	       0.4029312272760483E+1,
	       0.4560203031431090E+1,
	       0.5156826768007481E+1,
	       0.58781144889155572E+1};
     
     wts1[]={0.7050727473210895E-3,
	     0.7107111654073120E-2,
	     0.2844188515941899E-1,
	     0.6660235171398239E-1,
	     0.1025785712747278,
	     0.1077502032531791,
	     0.7747156370638879E-1,
	     0.3763106373385135E-1,
	     0.1204873635560290E-1,
	     0.2453208613776865E-2,
	     0.3020309847850189E-3,
	     0.2092121075871870E-4,
	     0.7314637349679360E-6,
	     0.1080646863902574E-7,
	     0.4828081616137754E-10,
	     0.2840126937112534E-13};
        
     for (j1=0;j1<N1;j1++){
       wts1[j1]=wts1[j1]*exp(pow(absci1[j1],2.0));
     }
   }
   else if(option_ur == 160){   //Gaussian with N=16  y*exp(-y^2)
      
     
     absci1[0]=0.477579959543737674E-1;
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

   if(option_theta == 0){ //constant difference for Theta
     dh2=2*pi/(N2); //number of intervals=no. of ordinates
     for (j2=0;j2<N2;j2++){
       absci2[j2]= dh2*j2;
       //absci2(j2)=(j2-1)*dh2+0.5*dh2
	  wts2[j2]=dh2;
     }
   }
   else if(option_theta == 1){
     //three-eights for theta
     h=pi;                     //the range of quadrature from [-h,h];
     // np=3*n_int;             //number of ordinates for integration
     dh=2*h/n_int;             //anglular increment
     
     for (int i=0;i<n_int;i++){
       absci2[3*i]=-h+dh/3.0*(3.0*i);       //abs[3*i]  =-h+dh/3*(3*i);
       absci2[3*i+1]=-h+dh/3.0*(3.0*i+1);   //abs[3*i+1]=-h+dh/3*(3*i+1);
       absci2[3*i+2]=-h+dh/3.0*(3.0*i+2);   //abs[3*i+2]=-h+dh/3*(3*i+2);
       if (i == 0){
	 wts2[i]=dh/4.0;  //here 0 -pi and pi are the same number
       }
       else {
	 wts2[3*i]=dh/4.0;
       }
       wts2[3*i+1]=dh/8.0*3.0;
       wts2[3*i+2]=dh/8.0*3.0;
     }
	 //absci2(N2)=h;
	 //wts2(N2)=dh/8.0;
   }
    
   if(option_phi == 0){ //constant difference for Phi
     dh3=pi/(N3-1); //number of intervals < number of ordinates
     for (j3=0;j3<N3;j3++){
       absci3[j3]= dh3*j3;
       wts3[j3]=dh3*sin(absci3[j3]);
     }
     else if(option_phi == 1){ //three-eights for phi
       h_phi=pi;               //the range of quadrature from [0,h];
       np_phi=3*nphi_int+1;    //number of ordinates for integration
       dh_phi=h_phi/nphi_int;  //anglular increment
       
       for (int i=0;i<nphi_int;i++){
	 absci3[3*i]=dh_phi/3.0*(3.0*i);       //abs[3*i]  =-h+dh/3*(3*i);
	 absci3[3*i+1]=dh_phi/3.0*(3.0*i+1);   //abs[3*i+1]=-h+dh/3*(3*i+1);
	 absci3[3*i+2]=dh_phi/3.0*(3.0*i+2);   //abs[3*i+2]=-h+dh/3*(3*i+2);
	 if (i == 0){   //i==0
	   wts3[i]=dh_phi/8.0*sin(absci3[i]);
	 }
	 else{
	   wts3[3*i]=dh_phi/4.0*sin(absci3[3*i]);
	 }      
	 wts3[3*i+1]=dh_phi/8.0*3.0*sin(absci3[3*i+1]);
	 wts3[3*i+2]=dh_phi/8.0*3.0*sin(absci3[3*i+2]);
       }
       absci3[3*nphi_int]=h_phi*sin(absci3[3*nphi_int]);
       wts3[3*nphi_int]=dh_phi/8.0*sin(absci3[3*nphi_int]);
     }
   }

   //weights for calculating moments
   j=0;
   for (j1=0;j1<N1;j1++){ //Ur
     for (j2=0;j2<N2;j2++){ //theta
       for (j3=0;j3<N3;j3++){ //phi
	 cx[j]=absci1[j1]*cos(absci2[j2])*sin(absci3[j3]); //cx=Ur*cos(theta)*sin(phi)
	 cy[j]=absci1[j1]*sin(absci2[j2])*sin(absci3[j3]); //cy=Ur*sin(theta)*sin(phi)
	 cz[j]=absci1[j1]*cos(absci3[j3]);                 //Ur*cos(phi)
	 //constant diff ur,theta,phi
	 //dcxyz(j3,j2,j1)=(absci1(j1))^2.0*sin(absci3(j3))*dh1*dh2*dh3 %Ur^2.0 sin(phi)d(theta)d(phi)
	 dcxyz[j]=wts1[j1]*wts2[j2]*wts3[j3];
	 j=j+1;       
       }
     }
   }
   
 }
}
    
