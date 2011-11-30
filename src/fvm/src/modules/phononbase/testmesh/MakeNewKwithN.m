clear all

load LA_K.txt
load LA_freq.txt
load LA_vg.txt
load TA_K.txt
load TA_freq.txt
load TA_vg.txt
load LO_K.txt
load LO_freq.txt
load LO_vg.txt
load TO_K.txt
load TO_freq.txt
load TO_vg.txt

%The lattice constant is adjusted to make the specific heat 
%the correct value.  This just changes the volume of the BZ
a=.275e-9;  %lattice constant (m)
Kmax=pi/a;
ntheta=4;
nphi=16;
nK=10;
npol=6;
dtheta=pi/ntheta/2;
dphi=2*pi/nphi;

%relaxation time stuff
Tlat=300;
A=1.32e-45;
B=1.73e-19;
Bl=2e-24;  % from Asheghi et al. JAP 91 (8), pp. 5079 (2002)
C=137.39;

%LA_K=LA_K(2:nK+1);
%LA_freq=LA_freq(2:nK+1);
%LA_vg=LA_vg(2:nK+1);

[newLA_K,newLA_freq]=transform(LA_K,LA_freq,nK);
[newTA_K,newTA_freq]=transform(TA_K,TA_freq,nK);
[newLO_K,newLO_freq]=transform(LO_K,LO_freq,nK);
[newTO_K,newTO_freq]=transform(TO_K,TO_freq,nK);


% plot(TA_K,TA_freq,'o');
% hold on
% plot(newTA_K,newTA_freq,'*')

[newLA_K,newLA_vg]=transform(LA_K,LA_vg,nK);
[newTA_K,newTA_vg]=transform(TA_K,TA_vg,nK);
[newLO_K,newLO_vg]=transform(LO_K,LO_vg,nK);
[newTO_K,newTO_vg]=transform(TO_K,TO_vg,nK);

fp = fopen('SiliconIsotropicN300K.txt','w');

fprintf(fp,'%u\n',npol);
fprintf(fp,'%u\n',nK);
fprintf(fp,'%u\n',nphi*ntheta);

dK=newTA_K(1)*2*Kmax;

totalvolume=0.;

for k=1:1:nK
    for theta=1:1:ntheta
        th=(theta-1)*dtheta+dtheta/2;
        for phi=1:1:nphi
            
            weight=(Kmax*newTA_K(k))^2*sin(th)*dK*dtheta*dphi/8/pi^3*2;
            totalvolume=totalvolume+weight;
            ph=(phi-1)*dphi+dphi/2;
            
            x=sin(th)*sin(ph);
            y=sin(th)*cos(ph);
            z=cos(th);
            
            kx=newTA_K(k)*Kmax*x;
            ky=newTA_K(k)*Kmax*y;
            kz=newTA_K(k)*Kmax*z;
            
            vx_TA1=newTA_vg(k)*x;
            vy_TA1=newTA_vg(k)*y;
            vz_TA1=newTA_vg(k)*z;
            
            freq_TA1=newTA_freq(k);
            tau_TA1=(A*freq_TA1^4)+(B*Tlat*freq_TA1^2*exp(-C/Tlat));
            tau_TA1n=(Bl*freq_TA1^2*Tlat^3)^(-1);
            tau_TA2n=tau_TA1n;
            tau_TA1=tau_TA1^-1;
            tau_TA2=tau_TA1;
            
            vx_LA=newLA_vg(k)*x;
            vy_LA=newLA_vg(k)*y;
            vz_LA=newLA_vg(k)*z;
            
            freq_LA=newLA_freq(k);
            tau_LA=(A*freq_LA^4)+(B*Tlat*freq_LA^2*exp(-C/Tlat));
            tau_LAn=(Bl*freq_LA^2*Tlat^3)^(-1);
            tau_LA=tau_LA^-1;
            
            vx_TA2=newTA_vg(k)*x;
            vy_TA2=newTA_vg(k)*y;
            vz_TA2=newTA_vg(k)*z;
            
            vx_TO1=newTO_vg(k)*x;
            vy_TO1=newTO_vg(k)*y;
            vz_TO1=newTO_vg(k)*z;
            
            freq_TO1=newTO_freq(k);
            tau_TO1=(A*freq_TO1^4)+(B*Tlat*freq_TO1^2*exp(-C/Tlat));
            tau_TO1n=(Bl*freq_TO1^2*Tlat^3)^(-1);
            tau_TO2n=tau_TO1n; 
            tau_TO1=tau_TO1^-1;
            tau_TO2=tau_TO1;
            
            vx_TO2=newTO_vg(k)*x;
            vy_TO2=newTO_vg(k)*y;
            vz_TO2=newTO_vg(k)*z;
            
            vx_LO=newLO_vg(k)*x;
            vy_LO=newLO_vg(k)*y;
            vz_LO=newLO_vg(k)*z;
            
            freq_LO=newLO_freq(k);
            tau_LO=(A*freq_LO^4)+(B*Tlat*freq_LO^2*exp(-C/Tlat));
            tau_LOn=(Bl*freq_LA^2*Tlat^3)^(-1);
            tau_LO=tau_LO^-1;
            
            %TA1 is first
            fprintf(fp,'%f\n',1.0);
            fprintf(fp,'%f\n',weight);
            fprintf(fp,'%f\n',freq_TA1);
            fprintf(fp,'%f\n',kx);
            fprintf(fp,'%f\n',ky);
            fprintf(fp,'%f\n',kz);
            fprintf(fp,'%f\n',vx_TA1);
            fprintf(fp,'%f\n',vy_TA1);
            fprintf(fp,'%f\n',vz_TA1);
            fprintf(fp,'%.15e\n',tau_TA1);
            fprintf(fp,'%.15e\n',tau_TA1n);
            
            %TA2 is second
            fprintf(fp,'%f\n',2.0);
            fprintf(fp,'%f\n',weight);
            fprintf(fp,'%f\n',freq_TA1);
            fprintf(fp,'%f\n',kx);
            fprintf(fp,'%f\n',ky);
            fprintf(fp,'%f\n',kz);
            fprintf(fp,'%f\n',vx_TA1);
            fprintf(fp,'%f\n',vy_TA1);
            fprintf(fp,'%f\n',vz_TA1);
            fprintf(fp,'%.15e\n',tau_TA1);
            fprintf(fp,'%.15e\n',tau_TA1n);
            
            %LA is third
            fprintf(fp,'%f\n',3.0);
            fprintf(fp,'%f\n',weight);
            fprintf(fp,'%f\n',freq_LA);
            fprintf(fp,'%f\n',kx);
            fprintf(fp,'%f\n',ky);
            fprintf(fp,'%f\n',kz);
            fprintf(fp,'%f\n',vx_LA);
            fprintf(fp,'%f\n',vy_LA);
            fprintf(fp,'%f\n',vz_LA);
            fprintf(fp,'%.15e\n',tau_LA);
            fprintf(fp,'%.15e\n',tau_LAn);
            
            %TO1 is fourth
            fprintf(fp,'%f\n',4.0);
            fprintf(fp,'%f\n',weight);
            fprintf(fp,'%f\n',freq_TO1);
            fprintf(fp,'%f\n',kx);
            fprintf(fp,'%f\n',ky);
            fprintf(fp,'%f\n',kz);
            fprintf(fp,'%f\n',vx_TO1);
            fprintf(fp,'%f\n',vy_TO1);
            fprintf(fp,'%f\n',vz_TO1);
            fprintf(fp,'%.15e\n',tau_TO1);
            fprintf(fp,'%.15e\n',tau_TO1n);
            
            %TO2 is fifth
            fprintf(fp,'%f\n',5.0);
            fprintf(fp,'%f\n',weight);
            fprintf(fp,'%f\n',freq_TO1);
            fprintf(fp,'%f\n',kx);
            fprintf(fp,'%f\n',ky);
            fprintf(fp,'%f\n',kz);
            fprintf(fp,'%f\n',vx_TO1);
            fprintf(fp,'%f\n',vy_TO1);
            fprintf(fp,'%f\n',vz_TO1);
            fprintf(fp,'%.15e\n',tau_TO1);
            fprintf(fp,'%.15e\n',tau_TO1n);
            
            %LO is sixth
            fprintf(fp,'%f\n',6.0);
            fprintf(fp,'%f\n',weight);
            fprintf(fp,'%f\n',freq_LO);
            fprintf(fp,'%f\n',kx);
            fprintf(fp,'%f\n',ky);
            fprintf(fp,'%f\n',kz);
            fprintf(fp,'%f\n',vx_LO);
            fprintf(fp,'%f\n',vy_LO);
            fprintf(fp,'%f\n',vz_LO);
            fprintf(fp,'%.15e\n',tau_LO);
            fprintf(fp,'%.15e\n',tau_LO);
        end
    end
end

totalvolume

fclose(fp);
