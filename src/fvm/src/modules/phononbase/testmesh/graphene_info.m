clear all
load output1.mat
fp = fopen(['graphene_data' num2str(temp) '.txt'],'w');

ntheta=180;
nphi=1;
nq=25;
npol=6;
index=0;
final_data=zeros(ntheta*nq*npol,8);

fprintf(fp,'%u\n',npol);
fprintf(fp,'%u\n',nq);
fprintf(fp,'%u\n',nphi*ntheta);

% kmax = dq(:,1)*25; %kmax(theta)=dq(theta)*nq (I don't actually need kmax)
for theta_=1:ntheta
    for q_=1:nq
        for pol_=1:npol
            index=index+1;
            final_data(index,1)=pol_; %Polarization
            final_data(index,2)=qqq(q_,theta_)*dq(theta_)*dtheta*pi/180/c0/4/pi^2; %Weights
            final_data(index,3)=freq(theta_,q_,pol_); %Frequency
            %I am assuming 25 divisions of dq are equally spaced
            final_data(index,4)=(dq(theta_)/2+(q_-1)*dq(theta_))*cos(theta_*pi/180); %kx
            final_data(index,5)=(dq(theta_)/2+(q_-1)*dq(theta_))*sin(theta_*pi/180); %ky
            final_data(index,6)=vgx(theta_,q_,pol_); %vx
            final_data(index,7)=vgy(theta_,q_,pol_); %vy
            final_data(index,8)=reltime(theta_,q_,pol_); %Reltime (U+N combined)
%             fprintf(fp, '%f %f %f %f %f %f %.10f %.15f\n',final_data(index,:));
            
            %%%%%%%%%% Begin write data %%%%%%%%%%%%%%%%%%%%%
            fprintf(fp,'%f\n', pol_); %Polarization
            fprintf(fp,'%f\n',final_data(index,2)); %Weight
            fprintf(fp,'%f\n',final_data(index,3)); %Frequency
            fprintf(fp,'%f\n',final_data(index,4)); %Kx
            fprintf(fp,'%f\n',final_data(index,5)); %Ky
            fprintf(fp,'%f\n',0); %Kz
            fprintf(fp,'%f\n',final_data(index,6)); %vgx
            fprintf(fp,'%f\n',final_data(index,7)); %vgy
            fprintf(fp,'%.15e\n',0); %vgz ??????
            fprintf(fp,'%.15e\n',final_data(index,8)); %Reltime
            fprintf(fp,'%.15e\n',0); %Reltime with N
            %%%%%%%%% End write data %%%%%%%%%%%%%%%%%%%%%
        end
    end
end
fclose(fp);


