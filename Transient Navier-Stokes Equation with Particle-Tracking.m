clc
clear all
close all

%% Initial Variable Definitions
Nx = 128;
Ny =64;
Lx =2;
Ly=1;
Re =100;
Dt = .002;
Delx = Lx./Nx;
Dely = Ly./Ny;
tTot = 0.002;
omega = 1;
y = 0.5*Dely:Dely:1-0.5*Dely;
x = 0.5*Delx:Delx:2-0.5*Delx;

%% initialization of constant parameters - particle flow
st = .1;             % stokes number
Np = 4;            % particle flow rate
%end_time = 50;
%% initialisation of partile position and velocity. Prticle velocity has been initialized as a Gaussian distribution
y_random = 0.1*rand(Np,1)+0.45;
for m = 1:Np
    xp(m) = 0;
    yp(m) = y_random(m);%0.45 + m*0.1/Np;
        up(m) = exp(-(yp(m)-0.5)^2/0.05^2);
            vp(m) = 0;
end

%% ADI LHS definition
LHS_U = ADI_U(Dt,Re,Delx,Dely,Nx,Ny);
LHS_V = ADI_V(Dt,Re,Delx,Dely,Nx,Ny);

P_LHS = zeros( Nx.*Ny,6);
P_LHS(:,1:5) = Poisson_lhs(Nx,Ny,Delx,Dely);
hunew = zeros((Nx+1),Ny);
hvnew = zeros((Nx),(Ny+1));

huold = hunew;
hvold =hvnew;

U = zeros((Nx+1).*Ny,1);
V = zeros((Ny+1).*Nx,1);
P = zeros(Nx.*Ny,1);

%%U velocity initialized with a Gaussian Profile for the left-side boundary
for j=1:1:(Nx+1).*Ny

    
    % v taken care during computation
    % P taken care during computation
    if(((j-0.5)*Dely< .45) || ((j-0.5)*Dely > .55))
        U(j,1) = 0;
    else
        U(j,1) = exp(-((j-0.5)*Dely-0.5)^2/0.05^2);
    
end
end
time1= 0;
time2 = 0;
time3 = 0;


tic;
for t=1:1:4000.*tTot/Dt

%before beginning of t loop

% for i = 2:Nx
%     for j = 2:Ny
%                 huold(i,j) = 0;
%                 hvold(i,j) = 0;
%     end
% end
%
% for i = 2:Nx
%     hu(i,1,1) = 0;
% end
% for j = 2:Ny
%     hv(1,j,1) = 0;
% end

% after beginning of time loop
u = reshape(U,[(Ny),(Nx+1)])';
v =reshape(V,[(Ny+1),Nx])';
% left boundary

% right boundary
% for i = 1:Ny
%     u(Nx+1,j) = (4*u(Nx,j) - u(Nx-1,j))/3;
%     p(Nx+1,j) = p(Nx,j);
%     v(Nx+1,j) = v(Nx,j);
% end
% top boundary

%     u(i,Ny+1) = u(i,Ny);
%     p(i,Ny+1) = p(i,Ny);


% bottom boundary

% u velocity taken care during computation
% P taken care during computation

%Convective terms

% hu is the convective term for u momentum
% hv is the convective term for v momentum
%% Convective term for u
for i = 2:(Nx+1)
    for j = 1:Ny
        if i==(Nx+1)
            hunew(i,j) = 0;
        else
            hunew(i,j) = (u(i+1,j) + u(i,j))^2 - (u(i,j) + u(i-1,j))^2;
        
        hunew(i,j) = hunew(i,j)*Dely;
        if(j == 1)
            hunew(i,j) = hunew(i,j) + ((u(i,j+1)+u(i,j))*(v(i-1,j+1)+v(i,j+1)) - 0) * Delx;
        elseif (j==Ny)
            hunew(i,j) = hunew(i,j) + (0 - (u(i,j)+u(i,j-1))*(v(i,j)+v(i-1,j))) * Delx;
        else
            hunew(i,j) = hunew(i,j) + ((u(i,j+1)+u(i,j))*(v(i-1,j+1)+v(i,j+1)) - (u(i,j)+u(i,j-1))*(v(i,j)+v(i-1,j))) * Delx;
        end
        end
        hunew(i,j) = 0.25*hunew(i,j);
    end
end
%% Convective terms for v

for i = 1:Nx
    for j = 2:(Ny)
        if(i==1)
            hvnew(i,j) = (u(i+1,j) + u(i+1,j-1))*(v(i+1,j) + v(i,j));
        elseif(i==Nx)
            hvnew(i,j) = -(u(i,j)+u(i,j-1))*(v(i,j)+v(i-1,j));
        else
            hvnew(i,j) = (u(i+1,j) + u(i+1,j-1))*(v(i+1,j) + v(i,j)) - (u(i,j)+u(i,j-1))*(v(i,j)+v(i-1,j));
        end
        hvnew(i,j) = hvnew(i,j)*Dely;
        hvnew(i,j) = hvnew(i,j) + ((v(i,j+1)+v(i,j))^2 - (v(i,j)+v(i,j-1))^2)*Delx;
        hvnew(i,j) = hvnew(i,j)*0.25;
    end
end
%% Diffusive terms for u
for i = 2:(Nx+1)
    for j = 1:Ny
        % diff_u is the diffusion term for u momentum
        if i==(Nx+1)
            diff_u =  (- 2.*(u(i,j) + u(i-1,j))*Dely/Delx);
        else
            diff_u = ((u(i+1,j) - 2*u(i,j) + u(i-1,j))*Dely/Delx);
        end
        if(j==1)
            diff_u = diff_u + ((u(i,j+1) - u(i,j))*Delx/Dely);
        elseif (j==Ny)
            diff_u = diff_u + ((u(i,j-1) - u(i,j))*Delx/Dely);
        else
            diff_u = diff_u + ((u(i,j+1) - 2*u(i,j) + u(i,j-1))*Delx/Dely);
        end
        diff_u = .5*diff_u /Re;
        
        src_u(i,j) = (-1.5*hunew(i,j) + 0.5*huold(i,j) + diff_u).*(Dt);
    end
end
src_u(1,1:Ny) = 0;
huold = hunew;

% Linear solver for u

src_u=reshape(src_u,[],1);
% temp = [LHS_U(:,1:3) src_u];
temp=LHS_U(:,1:3);
temp(:,4)=src_u;
src_u = TDMA(temp);

src_u=reshape(src_u,[Nx+1,Ny]);
src_u=reshape(src_u',[],1);
% temp = [LHS_U(:,4:6) src_u];
temp=LHS_U(:,4:6);
temp(:,4)=src_u;

src_u = reshape(src_u,[(Ny),(Nx+1)])';

Ucorr = TDMA(temp);

U = U+Ucorr;

% 

% after linear solver for u
for i = 1:Nx
    for j = 2:Ny
        % diff_v is the diffusion term for v momentum
        if(i == 1)
            diff_v = ((v(i+1,j) - 3*v(i,j))*Dely/Delx);
        elseif i==(Nx)
            diff_v = (( -1*v(i,j) + v(i-1,j))*Dely/Delx);
        else
            diff_v = ((v(i+1,j) - 2*v(i,j) + v(i-1,j))*Dely/Delx);
        end
        diff_v = diff_v + ((v(i,j+1) - 2*v(i,j) + v(i,j-1))*Delx/Dely);
        diff_v = .5*diff_v /Re;
        src_v(i,j) = (-1.5*hvnew(i,j) + 0.5*hvold(i,j) + diff_v).*Dt;
    end
end
hvold = hvnew;
src_v(1:Nx,1) = 0;
src_v(1:Nx,(Ny+1)) = 0;
% %after linear solver for v
% %% Linear solver for v
src_v=reshape(src_v,[],1);
%  temp = [LHS_V(:,1:3) src_v];
temp=LHS_V(:,1:3);
temp(:,4)=src_v;
src_v =TDMA(temp);
%  
%  
%  
% 
% 
src_v=reshape(src_v,[Nx,Ny+1]);

src_v=reshape(src_v',[],1);
% 
% 
% 
% 
% 
% 
% temp = [LHS_V(:,4:6) src_v];
temp=LHS_V(:,4:6);
temp(:,4)=src_v;
src_v = reshape(src_v,[(Ny+1),(Nx)])';


Vcorr = TDMA(temp);

V = V+Vcorr;
time1 = time1 + toc;
%% Presure Poisson


if (mean(U(Ny.*Nx+1:(Nx+1).*Ny,1))<10^-4)
Factor = mean(U(1:Ny,1)) - mean(U(Ny.*Nx+1:(Nx+1).*Ny,1));
U(Ny.*Nx+1:(Nx+1).*Ny,1) = Factor + U(Ny.*Nx-Ny+1:Nx.*Ny,1);
else
    Factor = mean(U(1:Ny,1))./mean(U(Ny.*Nx+1:(Nx+1).*Ny,1));
    U(Ny.*Nx+1:(Nx+1).*Ny,1) = Factor .* U(Ny.*Nx-Ny+1:Nx.*Ny,1);
end
for i=1:1:Nx
for j=1:1:Ny
P_LHS(Ny.*(i-1)+j,6) = ((1./Dt.*1).*((U(Ny.*(i)+j) - U(Ny.*(i-1)+j)).*Dely + (V((Ny+1).*(i-1)+j+1) - V((Ny+1).*(i-1)+j)).*Delx));
end
end
P = poisSolver(P_LHS,Nx,Ny,omega,P,t);
%% % Trial linsolve
% K = zeros(Nx*Ny);
% for i=1:1:Nx
%     for j=1:1:Ny
%         if (i>1)
%              K(Ny.*(i-1)+j,Ny.*(i-2)+j) = P_LHS(Ny.*(i-1)+j,4);
%         end
%         if(i<(Nx))
%              K(Ny.*(i-1)+j,Ny.*(i)+j) = P_LHS(Ny.*(i-1)+j,3);
%         end
%         if(j>1)
%             K(Ny.*(i-1)+j,Ny.*(i-1)+j-1) = P_LHS(Ny.*(i-1)+j,2);
%         end
%         if(j<Ny)
%             K(Ny.*(i-1)+j,Ny.*(i-1)+j+1) = P_LHS(Ny.*(i-1)+j,1);
%         end
%         K(Ny.*(i-1)+j,Ny.*(i-1)+j) = P_LHS(Ny.*(i-1)+j,5);
%     end
% end
% 
% Ptrial = K\P_LHS(:,6);%linsolve(K,P_LHS(:,6));

%% Final velocity update

for i=1:1:(Nx+1)
for j=1:1:Ny
if i==1||i==(Nx+1)
U(Ny.*(i-1)+j) = U(Ny.*(i-1)+j);
else
U(Ny.*(i-1)+j) = U(Ny.*(i-1)+j) -(Dt.*1./Delx).*(P(Ny.*(i-1)+j) -P(Ny.*(i-2)+j));
% U(Ny.*(i-1)+j) = U(Ny.*(i-1)+j) -(Dt./Delx).*(Ptrial(Ny.*(i-1)+j) -Ptrial(Ny.*(i-2)+j));
end
end
end

for i=1:1:(Nx)
for j=1:1:(Ny+1)
if j==1||j==(Ny+1)
V((Ny+1).*(i-1)+j) = V((Ny+1).*(i-1)+j);
else
V((Ny+1).*(i-1)+j) = V((Ny+1).*(i-1)+j) -(Dt.*1./Dely).*(P(Ny.*(i-1)+j) -P(Ny.*(i-1)+j-1));
%  V((Ny+1).*(i-1)+j) = V((Ny+1).*(i-1)+j) -(Dt./Dely).*(Ptrial(Ny.*(i-1)+j) -Ptrial(Ny.*(i-1)+j-1));
end
end
end

U = reshape(U,Ny,Nx+1)';
V = reshape(V,Ny+1,Nx)';
time2 = time2+toc;
tic;

%for t = 1:end_time
    count = 0;
    loop_count = length(xp);            % loop_count gives the number of particles
    
for m = 1:loop_count
    
    xp1 = xp(m);
    yp1 = yp(m);
    up1 = up(m);
    vp1 = vp(m);

%2nd order RK method
    for r = 1:2

        if((xp1 >= 0) && (xp1 <= Lx) && (yp1 >= 0) && (yp1 <= Ly))                
        if(xp1<0.5*Delx)
            i = 0;
            x_factor = xp1/Delx;
            else
                i = floor((xp1-0.5*Delx)/Delx) + 1;
                x_factor = (xp1 - x(i))/Delx;
        end
            
                j = floor(yp1/Dely) + 1;
                y_factor = (yp1 - (y(j)-0.5*Dely))/Dely;

        if(i == 0)
            vnw = 0;
        else if(j==Ny)
                vnw = 0;
            else
                vnw = V(i,j+1);
            end
        end
        if(i == Nx)
            if(j == Ny)
                vne = 0;
            else
            vne = V(i,j+1);
            end
        else if(j == Ny)
                vne = 0;
            else
            vne = V(i+1,j+1);
            end
        end
        if(i == 0)
            vsw = 0;
        else
            vsw = V(i,j);
        end
        if(i == Nx)
            vse = V(i,j);
        else
        vse = V(i+1,j);
        end
                vfp = (1-y_factor)*((1-x_factor)*vsw + x_factor*vse);
                vfp = vfp + y_factor*((1-x_factor)*vnw + x_factor*vne);
 
        if(yp1<0.5*Dely)
            j = 0;
            y_factor = yp1/Dely;
            else
                j = floor((yp1-0.5*Dely)/Dely) + 1;
                y_factor = (yp1 - y(j))/Delx;
        end

            i = floor((xp1)/Delx) + 1;
            x_factor = (xp1 - (x(i)-0.5*Delx))/Delx;
        
       if(j == Ny)
           unw = U(i,j);
       else
           unw = U(i,j+1);
       end
       if(j == Ny)
           if(i == Nx)
               une = U(i,j);
           else
               U(i+1,j);
           end
       else if(i == Nx)
               une = U(i,j+1);
           else         
               une = U(i+1,j+1);
           end
       end
       if(j == 0)
           usw = U(i,1);
       else
           usw = U(i,j);
       end
       if(i == Nx)
           if(j == 0)
               use = U(Nx,1);
           else
           use = U(Nx,j);
           end
       else
           if(j == 0)
               use = U(i+1,1);
           else
               use = U(i+1,j);
           end
       end
            ufp = (1-y_factor)*((1-x_factor)*usw + x_factor*use);
            ufp = ufp + y_factor*((1-x_factor)*unw + x_factor*une);
     
% simultaneous 2rk

k(1,r) = up1;
k(2,r) = vp1;
k(3,r) = (ufp - up1)/st;
k(4,r) = (vfp - vp1)/st;
if(r ==1)
    xp1 = xp1 + Dt*k(1,r);
    yp1 = yp1 + Dt*k(2,r);
    up1 = up1 + Dt*k(3,r);
    vp1 = vp1 + Dt*k(4,r);
end
        else 
        xp1 = -1;
        yp1 = -1;
        end
    
              
    end

% isolating particles out of domain
if(xp1 == -1)
    xp(m) = -1;
yp(m) = -1;
up(m) = -1;
vp(m) = -1;
count = count+1;
else
    % copying into temp location, partciles within the boundary
xp_new(m-count) = xp(m) + 0.5*Dt*(k(1,1)+k(1,2));
yp_new(m-count) = yp(m) + 0.5*Dt*(k(2,1)+k(2,2));
up_new(m-count) = up(m) + 0.5*Dt*(k(3,1)+k(3,2));
vp_new(m-count) = vp(m) + 0.5*Dt*(k(4,1)+k(4,2));
end
end
% reassigning the particles to the original variables
clearvars xp yp up vp;
xp = xp_new;
yp = yp_new;
up = up_new;
vp = vp_new;
% clearing the temp variables
clearvars xp_new yp_new up_new vp_new;
% scatter(xp,yp);
% axis([0 2 0 1]);
loop_count = length(xp);
% adding Np particles for the next time step
y_random = 0.1*rand(Np,1)+0.45;
    for m = 1:Np
    xp(loop_count+m) = 0;
    yp(loop_count+m) = y_random(m);%0.45 + m*0.1/Np;
        up(loop_count+m) = exp(-(yp(m)-0.5)^2/0.05^2);
            vp(loop_count+m) = 0;
    end
%end

% U=reshape(U,[],1);
U=reshape(U',[],1);

% V=reshape(V,[],1);
V=reshape(V',[],1);
time3 = time3+toc;
if(t==50)
    U1=U;
    V1=V;
    X1=xp';
    Y1=yp';
end
if(t==100)
    U2=U;
    V2=V;
    X2=xp';
    Y2=yp';
end
if(t==150)
    U3=U;
    V3=V;
    X3=xp';
    Y3=yp';
end
if(t==500)
    U4=U;
    V4=V;
    X4=xp';
    Y4=yp';
end
if(t==1000)
    U5=U;
    V5=V;
    X5=xp';
    Y5=yp';
    
end
if(t==2000)
    U6=U;
    V6=V;
    X6=xp';
    Y6=yp';
end
    

end

U = reshape(U,Ny,Nx+1);
V = reshape(V,Ny+1,Nx);

contourf(U,20)
hold on
grid on
% scatter(xp,yp);
% axis([0 2 0 1]);


