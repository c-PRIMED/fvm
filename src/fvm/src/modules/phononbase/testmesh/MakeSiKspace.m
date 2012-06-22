clear all

load SiKpts.txt;
load SiLAfreq.txt;
load SiLOfreq.txt;
load SiTAfreq.txt;
load SiTOfreq.txt;
load SiLAvg.txt;
load SiLOvg.txt;
load SiTAvg.txt;
load SiTOvg.txt;

ntheta=8;
nphi=4*ntheta;
nK=8;
npol=4;
dtheta=pi/ntheta/2;   %only does the top hemisphere
dphi=2*pi/nphi;

thstr=int2str(ntheta);
phstr=int2str(nphi);
nkstr=int2str(nK);
filename=['SiIsoEDIP_',thstr,'x',phstr,'x',nkstr,'.txt'];

olddk=SiKpts(1)*2;
Kmax=olddk*length(SiKpts);
dK=Kmax/nK;

[newK,newSiLAfreq]=InterpolateW(SiKpts,SiLAfreq,nK);
[newK,newSiLOfreq]=InterpolateW(SiKpts,SiLOfreq,nK);
[newK,newSiTAfreq]=InterpolateW(SiKpts,SiTAfreq,nK);
[newK,newSiTOfreq]=InterpolateW(SiKpts,SiTOfreq,nK);

[newSiLAvg]=InterpolateV(SiKpts,SiLAfreq,SiLAvg,newK,newSiLAfreq);
[newSiLOvg]=InterpolateV(SiKpts,SiLOfreq,SiLOvg,newK,newSiLOfreq);
[newSiTAvg]=InterpolateV(SiKpts,SiTAfreq,SiTAvg,newK,newSiTAfreq);
[newSiTOvg]=InterpolateV(SiKpts,SiTOfreq,SiTOvg,newK,newSiTOfreq);

[newSiLAtau]=InterpolateTau(SiKpts,SiLAfreq,SiLAvg,newK,newSiLAfreq,newSiLAvg);
[newSiLOtau]=InterpolateTau(SiKpts,SiLOfreq,SiLOvg,newK,newSiLOfreq,newSiLOvg);
[newSiTAtau]=InterpolateTau(SiKpts,SiTAfreq,SiTAvg,newK,newSiTAfreq,newSiTAvg);
[newSiTOtau]=InterpolateTau(SiKpts,SiTOfreq,SiTOvg,newK,newSiTOfreq,newSiTOvg);

fp = fopen(filename,'w');

fprintf(fp,'%u\n',npol);
fprintf(fp,'%u\n',nK);
fprintf(fp,'%u\n',nphi*ntheta);

totalvolume=0.;
count=0;

for k=1:1:nK
    for theta=1:1:ntheta
        th=(theta-1)*dtheta+dtheta/2;
        for phi=1:1:nphi
            
	    kweight=(newK(k))^2*dK+dK^3/12.;
            weight=kweight*2*sin(th)*sin(dtheta/2)*dphi/8/pi^3*2;
            totalvolume=totalvolume+weight;
            ph=(phi-1)*dphi+dphi/2;
            
            x=sin(th)*sin(ph);
            y=sin(th)*cos(ph);
            z=cos(th);
            
            kx=newK(k)*x;
            ky=newK(k)*y;
            kz=newK(k)*z;
            
            vx_LA=newSiLAvg(k)*x;
            vy_LA=newSiLAvg(k)*y;
            vz_LA=newSiLAvg(k)*z;
            
            vx_TA=newSiTAvg(k)*x;
            vy_TA=newSiTAvg(k)*y;
            vz_TA=newSiTAvg(k)*z;
            
            vx_TO=newSiTOvg(k)*x;
            vy_TO=newSiTOvg(k)*y;
            vz_TO=newSiTOvg(k)*z;
            
            vx_LO=newSiLOvg(k)*x;
            vy_LO=newSiLOvg(k)*y;
            vz_LO=newSiLOvg(k)*z;

            freq_TA=newSiTAfreq(k);
            freq_LA=newSiLAfreq(k);
            freq_TO=newSiTOfreq(k);
            freq_LO=newSiLOfreq(k);

            tau_TA=newSiTAtau(k);
            tau_LA=newSiLAtau(k);
            tau_TO=newSiTOtau(k);
            tau_LO=newSiLOtau(k);
            
            %TA is first
            fprintf(fp,'%f\n',1.0);
            fprintf(fp,'%f\n',weight);
            fprintf(fp,'%f\n',freq_TA);
            fprintf(fp,'%f\n',kx);
            fprintf(fp,'%f\n',ky);
            fprintf(fp,'%f\n',kz);
            fprintf(fp,'%f\n',vx_TA);
            fprintf(fp,'%f\n',vy_TA);
            fprintf(fp,'%f\n',vz_TA);
            fprintf(fp,'%.15e\n',tau_TA);
            fprintf(fp,'%.15e\n',tau_TA);
            
            %LA is second
            fprintf(fp,'%f\n',2.0);
            fprintf(fp,'%f\n',weight);
            fprintf(fp,'%f\n',freq_LA);
            fprintf(fp,'%f\n',kx);
            fprintf(fp,'%f\n',ky);
            fprintf(fp,'%f\n',kz);
            fprintf(fp,'%f\n',vx_LA);
            fprintf(fp,'%f\n',vy_LA);
            fprintf(fp,'%f\n',vz_LA);
            fprintf(fp,'%.15e\n',tau_LA);
            fprintf(fp,'%.15e\n',tau_LA);
            
            %TO is third
            fprintf(fp,'%f\n',3.0);
            fprintf(fp,'%f\n',weight);
            fprintf(fp,'%f\n',freq_TO);
            fprintf(fp,'%f\n',kx);
            fprintf(fp,'%f\n',ky);
            fprintf(fp,'%f\n',kz);
            fprintf(fp,'%f\n',vx_TO);
            fprintf(fp,'%f\n',vy_TO);
            fprintf(fp,'%f\n',vz_TO);
            fprintf(fp,'%.15e\n',tau_TO);
            fprintf(fp,'%.15e\n',tau_TO);
            
            %LO is fourth
            fprintf(fp,'%f\n',4.0);
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

            count=count+1;

        end
    end
end

	totalvolume/(4/3*pi*Kmax^3/8/pi^3)    %should come out to one
	count*npol
fclose(fp);
