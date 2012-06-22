function [ outKvec, outW ] = InterpolateW( inKvec, inW, desiredPanels )
%InterpolateW -- makes a smaller dispersion

  inPts=length(inKvec);
dk=inKvec(1)*2;
Kmax=dk*(inPts);
newdk=Kmax/desiredPanels;
outKvec(1:desiredPanels)=0;
outW(1:desiredPanels)=0;
kb=8.617343e-5;  %eV
hbar=6.582119e-16;  %eV s
T=300;

for k=1:1:desiredPanels
	klow=newdk*(k-1)/Kmax;
        khigh=newdk*k/Kmax;
        coarseKvol=4/3*pi*(khigh^3-klow^3);
        coarseCpTot=0;
        aveW=0;

        for kk=1:1:inPts                 
                 oldKmid=inKvec(kk)/Kmax;
                 oldKlow=oldKmid-dk/2/Kmax;
                 oldKhigh=oldKmid+dk/2/Kmax;
                 x=hbar*inW(kk)/kb/T;
                 cp=x^2*exp(x)/(exp(x)-1)^2;
                 if oldKlow>=klow && oldKhigh<khigh
                   coarseCpTot=coarseCpTot+4/3*pi*((oldKhigh)^3-(oldKlow)^3)*cp;
                   aveW=aveW+4/3*pi*((oldKhigh)^3-(oldKlow)^3)*x;
                 end
		 if oldKlow<klow && oldKhigh>klow
                   coarseCpTot=coarseCpTot+4/3*pi*((oldKhigh)^3-(klow)^3)*cp;
                   aveW=aveW+4/3*pi*((oldKhigh)^3-(klow)^3)*x;
		 end
		 if oldKlow<khigh && oldKhigh>khigh
                   coarseCpTot=coarseCpTot+4/3*pi*((khigh)^3-(oldKlow)^3)*cp;
                   aveW=aveW+4/3*pi*((khigh)^3-(oldKlow)^3)*x;
		 end
        end

	aveW=aveW/coarseKvol;
        x=findX(aveW,coarseCpTot/coarseKvol);
        outW(k)=x*kb*T/hbar;
        outKvec(k)=newdk*(k-.5);
        %outKvec(k)=(((khigh)^3+(klow))/2)^(1/3)*Kmax;
    
end

end
