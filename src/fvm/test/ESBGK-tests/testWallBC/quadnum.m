% quadrature numbering
clc; clear all; close all;
N1=8;N2=8;N3=4;

abx=linspace(-3,3,N1);
aby=linspace(-3,3,N2);
abz=linspace(-3,3,N3);
j=1;tab=[];
for j1=1:N1
    for j2=1:N2
        for j3=1:N3
            cx(j)=abx(j1);
            cy(j)=aby(j2);
            cz(j)=abz(j3);
            
            jn=j3+(j2-1)*(N3)+(j1-1)*(N2)*(N3);
            tab=[tab;[j,j1,j2,j3]];
            j=j+1;
            wts(j)=j;
        end
    end
end

%incident
ux=3;uy=3;uz=-3;
dcx=(abx(N1)-abx(1))/(N1-1);
dcy=(aby(N2)-aby(1))/(N2-1);
dcz=(abz(N3)-abz(1))/(N3-1);

ix=(ux-abx(1))/dcx+1;
iy=(uy-aby(1))/dcy+1;
iz=(uz-abz(1))/dcy+1;
ixm1=ix-1;
iym1=iy-1;
izm1=iz-1;
ixyz=iz+(iy-1)*(N3)+(ix-1)*(N2*N3);
jn=izm1+iym1*N3+ixm1*N2*N3;
[cx(jn+1),cy(jn+1),cz(jn+1)]

N123=N1*N2*N3;
index=linspace(1,N123,N123);
[index',cx',cy',cz',tab];
%reflected
vx=-1;vy=3;vz=3;