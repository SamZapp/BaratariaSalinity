function HS=getextraHsea(HS,h,MASK,dx,N,M);

% LL=2000/dx;
% HSextra=MASK*0;
% %HSextra(end-100:end,:)=0.5*[0:100]'/100*ones(1,M);
% HSextra(end-LL:end,221-5:263+5)=1*exp(-[LL:-1:0]'*dx/100*0.5)*ones(1,263-221+1+2*5);
% HSextra(end,:)=0;
% %figure;imagesc(HSextra');pause

%G=MASK*0;G(end,221-5-10+10:263-20)=1;
G=MASK*0;G(end-1:end,285:315)=1;
GG=double(bwdist(G)*dx);
%HSextra=1*exp(-GG/300);%USED UNTIL JULY10th
%HSextra(GG>1500)=0;%USED UNTIL JULY10th
HSextra=1*exp(-GG/2500);%new test
HSextra(GG>12000)=0;%new test
HSextra(MASK==0)=0;


HSextra(:,301-32:end)=0;
HSextra(1:end-39,301-38:end)=0;
HSextra(1:end-45,301-43:end)=0;

%size(HSextra)
%size(HS)
%pause
%figure;imagesc(GG);pause
%figure;imagesc(HSextra);pause
HS=HS+min(HSextra,h*0.5);