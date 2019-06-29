function [LHS_U ] = ADI_U( 	Dt,Re,Delx,Dely, Nx, Ny);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

temp1 = zeros((Nx+1).*Ny,1);
temp2 =temp1;
temp3 = temp1;
temp1(1:Ny) = 0;
temp1((Ny+1):(Ny.*Nx)) = (-Dt/(2.*Re))./Delx;
temp1((Ny.*Nx+1):(Ny.*(Nx+1))) = (-Dt./(Re))./Delx;
% temp1=temp1./(Dely.*Delx);

temp1=reshape(temp1,[Ny,(Nx+1)]);
temp1=temp1';
temp1=reshape(temp1,[],1);



temp2(1:Ny) = 1;
temp2((Ny+1):(Ny.*Nx)) = Delx+(Dt./(1.*Re))./Delx;
temp2((Ny.*Nx+1):(Ny.*(Nx+1))) = Delx+(Dt./(1.*Re))./Delx; 
% temp2 = temp2./(Dely.*Delx);

temp2=reshape(temp2,[Ny,Nx+1]);
temp2=temp2';
temp2=reshape(temp2,[],1);



temp3(1:Ny) = 0;
temp3((Ny+1):(Ny.*Nx)) = (-Dt./(2.*Re))./Delx;
temp3((Ny.*Nx+1):(Ny.*(Nx+1))) = 0; 
% temp3 = temp3./(Dely.*Delx);

temp3=reshape(temp3,[Ny,Nx+1]);
temp3=temp3';
temp3=reshape(temp3,[],1);

LHS_U(:,1:3) =[ temp3 temp2 temp1 ];


temp1(1:(Nx+1)) = 0;
temp1((Nx+2):(Nx+1).*Ny) = (-Dt./(2.*Re))./Dely;
%temp1 = temp1./(Delx.*Dely);

 temp1=reshape(temp1,[(Nx+1),Ny]);
 temp1=temp1';
 temp1=reshape(temp1,[],1);
 temp1(1:Ny) = 0;   
temp2(1:(Nx+1)) =  Dely+(Dt./(2.*Re))./Dely;
temp2((Nx+2):((Nx+1).*(Ny-1))) =  Dely+(Dt./(1.*Re))./Dely;
temp2((Nx+1).*(Ny-1)+1:(Nx+1).*Ny) = Dely+(Dt./(2.*Re))./Dely;
% temp2 = temp2./(Dely.*Delx);

 temp2=reshape(temp2,[(Nx+1),(Ny)]);
 temp2=temp2';
 temp2=reshape(temp2,[],1);
 temp2(1:Ny) =1;

temp3(1:((Nx+1).*(Ny-1))) = (-Dt./(2.*Re))./Dely;
temp3(((Nx+1).*(Ny-1)+1):(Nx+1).*Ny) = 0;
% temp3 = temp3./(Delx.*Dely);


temp3=reshape(temp3,[(Nx+1),(Ny)]);
temp3=temp3';
temp3=reshape(temp3,[],1);
temp3(1:Ny)=0;

LHS_U(:,4:6) = [ temp3 temp2 temp1 ];

end

