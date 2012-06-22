function [ newTau ] = InterpolateTau( oldKvec, oldW, oldV, newKvec, newW, newV )
%InterpolateTau -- makes a smaller dispersion

  inPts=length(oldKvec);
desiredPanels=length(newKvec);
dk=oldKvec(1)*2;
Kmax=dk*(inPts);
newdk=Kmax/desiredPanels;
outV(1:desiredPanels)=0;
kb=8.617343e-5;  %eV
hbar=6.582119e-16;  %eV s
T=300;

A=1.32e-45;
B=1.73e-19;
Bl=2e-24;  % from Asheghi et al. JAP 91 (8), pp. 5079 (2002)
C=137.39;

for k=1:1:desiredPanels
	klow=newdk*(k-1)/Kmax;
        khigh=newdk*k/Kmax;
        coarsex=hbar*newW(k)/kb/T;
        coarseKvol=4/3*pi*(khigh^3-klow^3);
        fineCpVTot=0;
        coarseCpV2=coarsex^2*exp(coarsex)/(exp(coarsex)-1)^2/kb*newV(k)^2;

        for kk=1:1:inPts                 
                 oldKmid=oldKvec(kk)/Kmax;
                 oldKlow=oldKmid-dk/2/Kmax;
                 oldKhigh=oldKmid+dk/2/Kmax;
                 x=hbar*oldW(kk)/kb/T;
                 cp=x^2*exp(x)/(exp(x)-1)^2/kb;
                 v=oldV(kk);
                 tau=A*oldW(kk)^4+B*T*oldW(kk)^2*exp(-C/T);
                 tau=tau^-1;
                 if oldKlow>=klow && oldKhigh<khigh
                   fineCpVTot=fineCpVTot+4/3*pi*((oldKhigh)^3-(oldKlow)^3)*cp*v^2*tau;
                 end
		 if oldKlow<klow && oldKhigh>klow
             	   fineCpVTot=fineCpVTot+4/3*pi*((oldKhigh)^3-(klow)^3)*cp*v^2*tau;
		 end
		 if oldKlow<khigh && oldKhigh>khigh
                   fineCpVTot=fineCpVTot+4/3*pi*((khigh)^3-(oldKlow)^3)*cp*v^2*tau;
		 end
        end

        newTau(k)=fineCpVTot/coarseKvol/coarseCpV2;
    
end

end
