function [ x ] = findX( guess, cpTot )
%findX -- calculates proper x=hbar*w/kb/T 

deltaX=1;
iters=0;

while deltaX>1e-10 && iters<100
     fx=guess^2*exp(guess)/(exp(guess)-1)^2;
     fpx=(guess*exp(guess)*(2+guess)*(exp(guess)-1)-2*guess^2*exp(2*guess))/(exp(guess)-1)^3;
     deltaX=(fx-cpTot)/fpx;
     newguess=guess-deltaX;
     deltaX=abs(deltaX)/guess;
     guess=newguess;
     iters=iters+1;
end

x=guess;

end
