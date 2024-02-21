function [Ind]=insideVT(V,Th,Z)
nT=size(Th,1);
tol=-1e3*eps; % some small number;
npix=size(Z,1);
xvec=Z(:,1); yvec=Z(:,2); zvec=Z(:,3);
L=HQbary(V,Th,xvec,yvec,zvec);
Ind=zeros(npix,1);
for j=1:nT
   I=find(L((j-1)*4+1,:)>=tol & L((j-1)*4+2,:)>=tol...
     & L((j-1)*4+3,:)>=tol & L((j-1)*4+4,:)>=tol);
   Ind(I)=1;
end;
