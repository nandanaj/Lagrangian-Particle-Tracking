function [ LHS_V  ] = ADI_V(Dt,Re,Delx,Dely, Nx, Ny );
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
temp1 = zeros((Nx).*(Ny+1),1);
temp2 =temp1;
temp3 = temp1;



temp1(1:(Ny+1)) = 0;
temp1((Ny+2):((Ny+1).*(Nx-1))) = (-Dt/(2.*Re))./Delx;
temp1(((Ny+1).*(Nx-1)+1):((Ny+1).*(Nx))) = (-Dt./(2.*Re))./Delx;


% temp1 = temp1./Delx;

temp1=reshape(temp1,[(Ny+1),(Nx)]);
temp1=temp1';
temp1=reshape(temp1,[],1);


temp2(1:(Ny+1)) = Delx+(3.*Dt./(2.*Re))./Delx;
temp2((Ny+2):((Ny+1).*(Nx-1))) = Delx+(Dt./(1.*Re))./Delx;
temp2(((Ny+1).*(Nx-1)+1):((Ny+1).*(Nx))) = Delx+(Dt./(2.*Re))./Delx; 


% 
temp2=(reshape(temp2,[(Ny+1),(Nx)]));
temp2=temp2';
temp2=reshape(temp2,[],1);



temp3((1:(Ny+1)))= (-Dt./(2.*Re))./Delx;
temp3((Ny+2):((Ny+1).*(Nx-1))) = (-Dt./(2.*Re))./Delx;
temp3(((Ny+1).*(Nx-1)+1):((Ny+1).*(Nx))) = 0; 

% 
temp3=(reshape(temp3,[(Ny+1),(Nx)]));
temp3=temp3';
temp3=reshape(temp3,[],1);


temp1(1:Nx) = 0;
temp2(1:Nx) = 1;
temp3(1:Nx) = 0;

temp1(Ny.*Nx+1:(Ny+1).*Nx) =0;
temp2(Ny.*Nx+1:(Ny+1).*Nx) = 1;
temp3(Ny.*Nx+1:(Ny+1).*Nx) = 0;

% temp1((Nx-1).*(Ny+1)+1:(Nx).*(Ny+1)) =0;
% temp2((Nx-1).*(Ny+1)+1:Nx.*(Ny+1)) = 1;
% temp3((Nx-1).*(Ny+1)+1:Nx.*(Ny+1)) = 0;


LHS_V(:,1:3) = [temp3 temp2 temp1];


temp1(1:(Nx)) = 0;
temp1((Nx+1):(Nx+1).*Ny-Nx) = (-Dt./(2.*Re))./Dely;
temp1((Nx+1).*Ny-Nx+1:(Ny+1).*Nx) = 0;


 temp1=reshape(temp1,[(Nx),Ny+1]);
 temp1=temp1';
 temp1=reshape(temp1,[],1);

temp2(1:(Nx)) =  1;
temp2((Nx+1):((Nx+1).*(Ny)-Nx)) =  Dely+(Dt./(1.*Re))./Dely;
temp2((Nx+1).*Ny-Nx+1:(Ny+1).*Nx) = 1;


 temp2=reshape(temp2,[(Nx),(Ny+1)]);
 temp2=temp2';
 temp2=reshape(temp2,[],1);


temp3(1:(Nx)) = 0;
temp3((Nx+1):(Nx+1).*Ny-Nx) = (-Dt./(2.*Re))./Dely;
temp3((Nx+1).*Ny-Nx+1:(Ny+1).*Nx) = 0;


temp3=reshape(temp3,[(Nx),(Ny+1)]);
temp3=temp3';
temp3=reshape(temp3,[],1);


LHS_V(:,4:6) = [temp3 temp2 temp1];



end

