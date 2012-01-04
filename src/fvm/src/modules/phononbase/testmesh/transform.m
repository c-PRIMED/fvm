function [ outKvec, outVal ] = transform( inKvec, inVal, desiredPoints )
%transform -- this takes an input function inVal=f(inKvec)
%             and changes the number of points in the curve


inPoints=length(inKvec);  %Kvec will always go from 0-1
dkOut=1.0/desiredPoints;

outKvec=dkOut/2.0:dkOut:(1.0-dkOut/2.0);
outVal(1:desiredPoints)=0.0;

for i=1:desiredPoints
    k=outKvec(i);
    for j=1:inPoints
        if k<inKvec(j)
            x2=inKvec(j);
            x1=inKvec(j-1);
            y2=inVal(j);
            y1=inVal(j-1);
            break;
        end
    end
    outVal(i)=findVal(x1,x2,y1,y2,k);
end

end

