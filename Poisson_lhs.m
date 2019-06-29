function [ M ] = Poisson_lhs( Nx,Ny,Delx,Dely );
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
M(1:Ny,3) = 0;
M((Ny+1):Nx.*Ny,3) = Dely./Delx;
M(Nx.*Ny-Ny+1:Nx.*Ny,4) = 0;
M(1:Nx.*Ny-Ny,4) = Dely./Delx;
M(:,1) = Delx./Dely;
M(:,2) =Delx./Dely;
for i=1:1:Nx
M(i.*Ny,1) = 0;
M((i-1).*Ny+1,2) = 0;
end
M(:,5) = -2.*(Dely./Delx+Delx./Dely);
M(Nx.*Ny-Ny+1:Nx.*Ny,5) = -2.*(Dely./Delx) - Delx./Dely;
M(1:Ny,5) = -2.*(Dely./Delx) - Delx./Dely;

for i=1:1:Nx
M(i.*Ny,5) = -(Dely./Delx) - 2.*Delx./Dely;
M((i-1).*Ny+1,5) = -(Dely./Delx) -2.* Delx./Dely;
end
M(1,5) = -(Dely./Delx) - Delx./Dely;
M(Ny,5) = -(Dely./Delx) - Delx./Dely;
M(Nx.*Ny,5) = -(Dely./Delx) - Delx./Dely;
M(Nx.*Ny-Ny+1,5) = -(Dely./Delx) - Delx./Dely;

end

