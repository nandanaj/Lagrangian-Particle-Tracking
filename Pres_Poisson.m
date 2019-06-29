function [ Unew Vnew ] = Pres_Poisson( Unew, Vnew, Nx,Ny, Delx,Dely )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
for i=1:1:Nx
for j=1:1:Ny
M(Ny.*(i-1)+j,6) = (1./Dt).*((Unew(Ny.*(i)+j) - Unew(Ny.*(i-1)+j)).*Dely + (Vnew((Ny+1).*(i-1)+j+1) - Vnew((Ny+1).*(i-1)+j)).*Delx);
end
end



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
M(Nx.*Ny-Ny+1:Nx.*Ny,5) = -2.*(Dely./Delx) - Delx./Dely
M(1:Ny,5) = -2.*(Dely./Delx) - Delx./Dely;

for i=1:1:Nx
M(i.*Ny,5) = -(Dely./Delx) - 2.*Delx./Dely;
M((i-1).*Ny+1,5) = -(Dely./Delx) -2.* Delx./Dely;
end
M(1,5) = -(Dely./Delx) - Delx./Dely;
M(Ny,5) = -(Dely./Delx) - Delx./Dely;
M(Nx.*Ny,5) = -(Dely./Delx) - Delx./Dely;
M(Nx.*Ny-Ny+1,5) = -(Dely./Delx) - Delx./Dely;
omega =1.7;
Pnew = sorSolver(M,Nx,Ny,omega);

for i=1:1:(Nx+1)
for j=1:1:Ny
if i==1||i==(Nx+1)
Unew(Ny.*(i-1)+j) = Unew(Ny.*(i-1)+j);
else
Unew(Ny.*(i-1)+j) = Unew(Ny.*(i-1)+j) -(Dt./Delx).*(Pnew(Ny.*(i)+j) -Pnew(Ny.*(i-1)+j));
end
end
end

for i=1:1:(Nx)
for j=1:1:(Ny+1)
if i==1||i==(Nx+1)
Vnew(Ny.*(i-1)+j) = Vnew(Ny.*(i-1)+j);
else
Vnew(Ny.*(i-1)+j) = Vnew(Ny.*(i-1)+j) -(Dt./Dely).*(Pnew(Ny.*(i-1)+j+1) -Pnew(Ny.*(i-1)+j));
end
end
end

end

