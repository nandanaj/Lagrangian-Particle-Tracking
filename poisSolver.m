function [initPhi,L]=poisSolver(M,nx,ny,omega,initPhi,step)
an=M(:,1);
as=M(:,2);
aw=M(:,3);
ae=M(:,4);
ap=M(:,5);
source=M(:,6);
L2=1;
itr=1;
k=1;
% initPhi=zeros(nx*ny,1);
tempPhi=initPhi;
while((L2>1E-5)&&(itr<1000)) 
    L(itr)=0;
    
for i=1:nx*ny
    s=0;
    if(rem(i,ny)==0)   %North
        s=s-0;
    else
        s=s-an(i)*initPhi(i+1);
    end

    if(rem(i-1,ny)==0)   %South
        s=s-0;
    else
        s=s-as(i)*initPhi(i-1);
    end

    if(i<=nx*ny-ny)   %East
        s=s-ae(i)*initPhi(i+ny);
    else
        s=s-0;
    end
    
    if(i<=ny)   %West
        s=s-0;
    else
        s=s-aw(i)*initPhi(i-ny);
    end

    initPhi(i)=(omega*(source(i)+s)/ap(i))+(1-omega)*tempPhi(i);
    err(i,1)=tempPhi(i)-initPhi(i);
%     errsq(i,1)=err(i,1).^2;
%     L(itr)=L(itr)+errsq(i);

end
    L(itr)=(norm(err)./(nx*ny));
    L2=L(itr);
    error=max(abs(err));
    itr=itr+1; 
%     if itr==1000
%         break
%     end
    tempPhi=initPhi;
    fprintf('\nSOR: %d %d %d',itr-1,L2,step);
%     q=L(itr-1);
end
n=1;
end
