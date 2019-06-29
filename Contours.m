close all;
scale = 4;
linespec = 'black';
TIME = [0.1 0.2  1 2 4];
for i =1:Nx
    for j= 1: Ny
        Ucell(i,j) = U1(Ny*(i-1)+j) + U1(Ny*(i)+j);
        Ucell(i,j) = Ucell(i,j) * 0.5;
        Vcell(i,j) = V1(Ny*(i-1)+j) + V1(Ny*(i-1)+j+1);
        Vcell(i,j) = Vcell(i,j) * 0.5;
    end
end

q = 1;
figure(q);
scatter(X1,Y1,'.');
hold on
grid on;
quiver(x',y',Ucell',Vcell',scale,linespec);
axis([0 2 0 1]);
title(sprintf('Particle distribution in the flow field; Re = %d; St = %d; time = %d',Re,st,TIME(q)));
% str = '# of Particles: ';

q = q+1;

for i =1:Nx
    for j= 1: Ny
        Ucell(i,j) = U2(Ny*(i-1)+j) + U2(Ny*(i)+j);
        Ucell(i,j) = Ucell(i,j) * 0.5;
        Vcell(i,j) = V2(Ny*(i-1)+j) + V2(Ny*(i-1)+j+1);
        Vcell(i,j) = Vcell(i,j) * 0.5;
    end
end

figure(q);
scatter(X2,Y2,'.');
hold on
grid on;
quiver(x',y',Ucell',Vcell',scale,linespec);
axis([0 2 0 1]);
title(sprintf('Particle distribution in the flow field; Re = %d; St = %d; time = %d',Re,st,TIME(q)));

q = q+1;

for i =1:Nx
    for j= 1: Ny
        Ucell(i,j) = U3(Ny*(i-1)+j) + U3(Ny*(i)+j);
        Ucell(i,j) = Ucell(i,j) * 0.5;
        Vcell(i,j) = V3(Ny*(i-1)+j) + V3(Ny*(i-1)+j+1);
        Vcell(i,j) = Vcell(i,j) * 0.5;
    end
end
figure(q);
scatter(X3,Y3,'.');
hold on
grid on;
quiver(x',y',Ucell',Vcell',scale,linespec);
axis([0 2 0 1]);
title(sprintf('Particle distribution in the flow field; Re = %d; St = %d; time = %d',Re,st,TIME(q)));


q = q+1;

for i =1:Nx
    for j= 1: Ny
        Ucell(i,j) = U4(Ny*(i-1)+j) + U4(Ny*(i)+j);
        Ucell(i,j) = Ucell(i,j) * 0.5;
        Vcell(i,j) = V4(Ny*(i-1)+j) + V4(Ny*(i-1)+j+1);
        Vcell(i,j) = Vcell(i,j) * 0.5;
    end
end
figure(q);
scatter(X4,Y4,'.');
hold on
grid on;
quiver(x',y',Ucell',Vcell',scale,linespec);
axis([0 2 0 1]);
title(sprintf('Particle distribution in the flow field; Re = %d; St = %d; time = %d',Re,st,TIME(q)));


q = q + 1;

for i =1:Nx
    for j= 1: Ny
        Ucell(i,j) = U5(Ny*(i-1)+j) + U5(Ny*(i)+j);
        Ucell(i,j) = Ucell(i,j) * 0.5;
        Vcell(i,j) = V5(Ny*(i-1)+j) + V5(Ny*(i-1)+j+1);
        Vcell(i,j) = Vcell(i,j) * 0.5;
    end
end

figure(q);
scatter(X5,Y5,'.');
hold on
grid on;
quiver(x',y',Ucell',Vcell',scale,linespec);
axis([0 2 0 1]);
title(sprintf('Particle distribution in the flow field; Re = %d; St = %d; time = %d',Re,st,TIME(q)));


%%
pp=1;
for i =1:Nx
    for j= 1: Ny
        Ucell(i,j) = U1(Ny*(i-1)+j) + U1(Ny*(i)+j);
        Ucell(i,j) = Ucell(i,j) * 0.5;
        Vcell(i,j) = V1(Ny*(i-1)+j) + V1(Ny*(i-1)+j+1);
        Vcell(i,j) = Vcell(i,j) * 0.5;
    end
end
q=q+1;
figure(q);
scatter(X1,Y1,'.');
title(sprintf('Particle distribution in the flow field; Re = %d; St = %d; time = %d',Re,st,TIME(pp)));
xlabel('X');
ylabel('Y');
axis([0 2 0 1]);
grid on;
q=q+1;

figure(q);
quiver(x',y',Ucell',Vcell',scale,linespec);
title(sprintf('Flow Direction Re = %d; St = %d; time = %d',Re,st,TIME(pp)));
xlabel('X');
ylabel('Y');
axis([0 2 0 1]);

str = '# of Particles: ';

q = q+1;
pp=pp+1;

for i =1:Nx
    for j= 1: Ny
        Ucell(i,j) = U2(Ny*(i-1)+j) + U2(Ny*(i)+j);
        Ucell(i,j) = Ucell(i,j) * 0.5;
        Vcell(i,j) = V2(Ny*(i-1)+j) + V2(Ny*(i-1)+j+1);
        Vcell(i,j) = Vcell(i,j) * 0.5;
    end
end

figure(q);
scatter(X2,Y2,'.');
axis([0 2 0 1]);
title(sprintf('Particle distribution in the flow field; Re = %d; St = %d; time = %d',Re,st,TIME(pp)));
xlabel('X');
ylabel('Y');
grid on;
q=q+1;
figure(q);
quiver(x',y',Ucell',Vcell',scale,linespec);
title(sprintf('Flow Direction Re = %d; St = %d; time = %d',Re,st,TIME(pp)));
xlabel('X');
ylabel('Y');
axis([0 2 0 1]);

q = q+1;
pp=pp+1;
for i =1:Nx
    for j= 1: Ny
        Ucell(i,j) = U3(Ny*(i-1)+j) + U3(Ny*(i)+j);
        Ucell(i,j) = Ucell(i,j) * 0.5;
        Vcell(i,j) = V3(Ny*(i-1)+j) + V3(Ny*(i-1)+j+1);
        Vcell(i,j) = Vcell(i,j) * 0.5;
    end
end
figure(q);
scatter(X3,Y3,'.');
axis([0 2 0 1]);
title(sprintf('Particle distribution in the flow field; Re = %d; St = %d; time = %d',Re,st,TIME(pp)));
xlabel('X');
ylabel('Y');
grid on;
q=q+1;
figure(q);
quiver(x',y',Ucell',Vcell',scale,linespec);
title(sprintf('Flow Direction Re = %d; St = %d; time = %d',Re,st,TIME(pp)));
xlabel('X');
ylabel('Y');
axis([0 2 0 1]);

pp=pp+1;
q = q+1;

for i =1:Nx
    for j= 1: Ny
        Ucell(i,j) = U4(Ny*(i-1)+j) + U4(Ny*(i)+j);
        Ucell(i,j) = Ucell(i,j) * 0.5;
        Vcell(i,j) = V4(Ny*(i-1)+j) + V4(Ny*(i-1)+j+1);
        Vcell(i,j) = Vcell(i,j) * 0.5;
    end
end
figure(q);
scatter(X4,Y4,'.');
axis([0 2 0 1]);
title(sprintf('Particle distribution in the flow field; Re = %d; St = %d; time = %d',Re,st,TIME(pp)));
xlabel('X');
ylabel('Y');
grid on;
q=q+1;
figure(q);
quiver(x',y',Ucell',Vcell',scale,linespec);
title(sprintf('Flow Direction Re = %d; St = %d; time = %d',Re,st,TIME(pp)));
xlabel('X');
ylabel('Y');
axis([0 2 0 1]);

pp=pp+1;
q = q + 1;

for i =1:Nx
    for j= 1: Ny
        Ucell(i,j) = U5(Ny*(i-1)+j) + U5(Ny*(i)+j);
        Ucell(i,j) = Ucell(i,j) * 0.5;
        Vcell(i,j) = V5(Ny*(i-1)+j) + V5(Ny*(i-1)+j+1);
        Vcell(i,j) = Vcell(i,j) * 0.5;
    end
end

figure(q);
scatter(X5,Y5,'.');
axis([0 2 0 1]);
title(sprintf('Particle distribution in the flow field; Re = %d; St = %d; time = %d',Re,st,TIME(pp)));
xlabel('X');
ylabel('Y');
grid on;
q=q+1;
figure(q);
quiver(x',y',Ucell',Vcell',scale,linespec);
title(sprintf('Flow Direction Re = %d; St = %d; time = %d',Re,st,TIME(pp)));
xlabel('X');
ylabel('Y');
axis([0 2 0 1]);

%%
pp=0;
UU1=reshape(U1,Ny,Nx+1);
UU2=reshape(U2,Ny,Nx+1);
UU3=reshape(U3,Ny,Nx+1);
UU4=reshape(U4,Ny,Nx+1);
UU5=reshape(U5,Ny,Nx+1);
pp=pp+1;
q=q+1;
figure(q);
contourf(UU1,20)
grid on
colorbar
xlabel('X');
ylabel('Y');
title(sprintf('X-Velocity: Re=%d St= %d time=%d s',Re,st,TIME(pp)));

pp=pp+1;
q=q+1;
figure(q);
contourf(UU2,20)
grid on
colorbar
xlabel('X');
ylabel('Y');
title(sprintf('X-Velocity: Re=%d St= %d time=%d s',Re,st,TIME(pp)));

pp=pp+1;
q=q+1;
figure(q);
contourf(UU3,20)
grid on
colorbar
xlabel('X');
ylabel('Y');
title(sprintf('X-Velocity: Re=%d St= %d time=%d s',Re,st,TIME(pp)));

pp=pp+1;
q=q+1;
figure(q);
contourf(UU4,20)
grid on
colorbar
xlabel('X');
ylabel('Y');
title(sprintf('X-Velocity: Re=%d St= %d time=%d s',Re,st,TIME(pp)));

pp=pp+1;
q=q+1;
figure(q);
contourf(UU5,20)
grid on
colorbar
xlabel('X');
ylabel('Y');
title(sprintf('X-Velocity: Re=%d St= %d time=%d s',Re,st,TIME(pp)));


