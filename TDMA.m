function [ y ] = TDMA( M )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
size = numel(M(:,1));
C = zeros(size,1);
A = zeros(size,1);
alfa=-M(:,1);
beta=-M(:,3);
D=M(:,2);
RHS=M(:,4);
for i=1:1:size
    if i==1
        C(i) = RHS(i)./D(i);
        A(i) = alfa(i)./D(i);
    else
        C(i) = (beta(i).*C(i-1)+RHS(i))./(D(i)-beta(i).*A(i-1));
        A(i) = alfa(i)./(D(i)-beta(i).*A(i-1));
    end
end

 for i=size:-1:1
     if i==size
         y(i,1)=C(i);
     else
         y(i,1) =A(i).*y(i+1)+C(i);
     end
 end
end