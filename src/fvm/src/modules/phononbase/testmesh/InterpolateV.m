function [ outV ] = InterpolateV( oldKvec, oldW, oldV, newKvec, newW  )
%InterpolateV -- makes a smaller dispersion

  inPts=length(oldKvec);
desiredPanels=length(newKvec);
dk=oldKvec(1)*2;
Kmax=dk*(inPts);
newdk=Kmax/desiredPanels;
outV(1:desiredPanels)=0;
kb=8.617343e-5;  %eV
hbar=6.582119e-16;  %eV s
T=300;

for k=1:1:desiredPanels
	klow=newdk*(k-1)/Kmax;
        khigh=newdk*k/Kmax;
        coarsex=hbar*newW(k)/kb/T;
        coarseKvol=4/3*pi*(khigh^3-klow^3);
        fineCpVTot=0;
        coarseCp=coarsex^2*exp(coarsex)/(exp(coarsex)-1)^2/kb;

        for kk=1:1:inPts                 
                 oldKmid=oldKvec(kk)/Kmax;
                 oldKlow=oldKmid-dk/2/Kmax;
                 oldKhigh=oldKmid+dk/2/Kmax;
                 x=hbar*oldW(kk)/kb/T;
                 cp=x^2*exp(x)/(exp(x)-1)^2/kb;
                 v=oldV(kk);
                 if oldKlow>=klow && oldKhigh<khigh
                   fineCpVTot=fineCpVTot+4/3*pi*((oldKhigh)^3-(oldKlow)^3)*cp*v;
                 end
		 if oldKlow<klow && oldKhigh>klow
                   fineCpVTot=fineCpVTot+4/3*pi*((oldKhigh)^3-(klow)^3)*cp*v;
		 end
		 if oldKlow<khigh && oldKhigh>khigh
                   fineCpVTot=fineCpVTot+4/3*pi*((khigh)^3-(oldKlow)^3)*cp*v;
		 end
        end

        outV(k)=fineCpVTot/coarseKvol/coarseCp;
    
end

end
