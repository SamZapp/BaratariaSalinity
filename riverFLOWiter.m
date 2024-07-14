function [U,Ux,Uy,P,q1,qm1,qN,qmN,U1,Um1,UN,UmN]=flowBasin(A,manning,h,ho,dx,fTide,Q,Uo,directionQ,imposenocrossflux,FcrUR,DD1,DDN,DDUlimit);
facADVh=0.5;

%A(h<0.2)=0;
%h=max(0.1,h);
%Uo=A*0+1;
Q=abs(Q);
%A(A==11)=1; %the river front cells are normal cells
%to impose no cross flux in inlet


%Uo=0.1;
%manning=0*A+n_manning;

%manning(h<0.2)=10;
%manning(VEG==1)=1;%figure;imagesc(manning);pause

A(A==22)=1;  %this behaves as normal flow %but do not update A!
%consider the pond cells as A==0
A(A==3)=1; %the isoalted pond behaves as normal cell (btu different depth...) %but do not update A!


%csi=h.^(1/3)  ./(manning.^2.*Uo);
csi=max(1,h).^(1/3)  ./(manning.^2.*Uo);
Icsi=1./csi;

G=0*A;a=find(A~=0);NN=length(a);G(a)=[1:NN];
rhs=zeros(NN,1);  
for i=1:length(Q)
    rhs(G(A==10+i-1))=Q(i)/dx;
end


[N,M] = size(G);i=[];j=[];s=[];

%boundary conditions imposed water level
a=find(A==2);
i=[i;G(a)]; j=[j;G(a)]; s=[s;ones(size(a))];rhs(G(a))=0;%water level zero

S=0*G;
%exclude the NOLAND CELLS (A==0)
p = find(A==1 | (A>=10 & A<=19));[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]

%avoid to the the cells out of the domain (risk to make it periodic...)
if k==N;a=find(col+1<=M);end;if k==-N;a=find(col-1>0);end;if k==-1;a=find(row-1>0);end;if k==1;a=find(row+1<=N);end;

q=p+k;%the translated cell
a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells


if abs(k)==1
DDU=h.*DD1/9.81;
elseif abs(k)==N
DDU=h.*DDN/9.81;
end
DDUm=sign(DDU(p(a))+DDU(q(a))).*min(abs(DDU(p(a))),abs(DDU(q(a)))).*(1-exp(- facADVh*min(h(p(a)),h(q(a))) )).*(h(p(a))>1 & h(q(a))>1);
Icsim=(Icsi(p(a))+Icsi(q(a)))/2;
hm=(h(p(a))+h(q(a)))/2;

val=min(0,max(-DDUlimit*Icsim,DDUm));
DD=1./(Icsim+val)   .*hm.^2 /(dx^2);


if imposenocrossflux==1;
for iii=1:length(Q)  
    if abs(directionQ(iii))==1 & abs(k)==N;
        cross=find(A(p(a))==10+iii-1);DD(cross)=0;
    elseif abs(directionQ(iii))==2 & abs(k)==1;
        cross=find(A(p(a))==10+iii-1);DD(cross)=0;
    end
end
end


S(p(a))=S(p(a))+DD; %exit from that cell
i=[i;G(q(a))]; j=[j;G(p(a))]; s=[s;-DD]; %gain from the neigborh cell
end

%summary of the material that exits the cell
i=[i;G(p)]; j=[j;G(p)]; s=[s;S(p)];

ds2 = sparse(i,j,s);%solve the matrix inversion

p=ds2\rhs;
P=G;P(G>0)=full(p(G(G>0)));
P(A==2 | A==21)=0;  %need when swtinching q and p

U1=0*A;Um1=0*A;UN=0*A;UmN=0*A;
q1=0*A;qm1=0*A;qN=0*A;qmN=0*A;
Ux=0*A;Uy=0*A;

p = find(A==1 | A==2 | (A>=10 & A<=19));[row col]=ind2sub(size(A),p);
for k = [N -1 1 -N]
if k==N; a=find(col+1<=M);end;if k==-N;a=find(col-1>0);end;if k==-1;a=find(row-1>0);end;if k==1; a=find(row+1<=N);end;
q=p+k;%the translated cell
a=a(A(q(a))>0);%exlcude the translated cell that are NOLAND cells






if abs(k)==1
DDU=h.*DD1/9.81;
elseif abs(k)==N
DDU=h.*DDN/9.81;
end
DDUm=sign(DDU(p(a))+DDU(q(a))).*min(abs(DDU(p(a))),abs(DDU(q(a)))).*(1-exp(- facADVh*min(h(p(a)),h(q(a))) )).*(h(p(a))>1 & h(q(a))>1);
Icsim=(Icsi(p(a))+Icsi(q(a)))/2;
hm=(h(p(a))+h(q(a)))/2;

val=min(0,max(-DDUlimit*Icsim,DDUm));
DD=1./(Icsim+val)   .*hm.^2 /(dx^2) *dx;

if imposenocrossflux==1;
for iii=1:length(Q)  
    if abs(directionQ(iii))==1 & abs(k)==N;
        cross=find(A(p(a))==10+iii-1);DD(cross)=0;
    elseif abs(directionQ(iii))==2 & abs(k)==1;
        cross=find(A(p(a))==10+iii-1);DD(cross)=0;
    end
end
end

if (k==1); q1(p(a))=q1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD; 
elseif (k==-1); qm1(p(a))=qm1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD; 
elseif (k==N); qN(p(a))=qN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD;  
elseif (k==-N); qmN(p(a))=qmN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD; 
end


if (k==1); U1(p(a))=U1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD./hm; 
elseif (k==-1); Um1(p(a))=Um1(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD./hm; 
elseif (k==N); UN(p(a))=UN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD./hm;  
elseif (k==-N); UmN(p(a))=UmN(p(a))+sign(k)*(P(p(a))-P(q(a))).*DD./hm; 
end


end







for i=1:length(Q)  
    
if directionQ(i)==1 %up
qm1(A==10+i-1)=q1(A==10+i-1);

elseif directionQ(i)==2 %left
qmN(A==10+i-1)=qN(A==10+i-1);

elseif directionQ(i)==-2 %right
qN(A==10+i-1)=qmN(A==10+i-1);
end

end






%ONLY THIS WORKS FOR THE RIVER!!!!
Ux=(q1+qm1)/2./h; 
Uy=(qN+qmN)/2./h;


U=sqrt(Ux.^2+Uy.^2);
Uo=min(U,FcrUR*sqrt(9.81*h));

Ux=Uo.*Ux./U;Ux(U==0)=0;
Uy=Uo.*Uy./U;Uy(U==0)=0;
U=sqrt(Ux.^2+Uy.^2);






















% 
% %this will impose more uniform discharge over the oulet
% for i=1:length(Q)
%     
%     
%     if directionQ(i)==1 %up
%     Uy(A==10+i-1)=Q(i)./h(A==10+i-1); Ux(A==10+i-1)=0;  %this will impose more uniform discharge over the oulet
%     qN(A==10+i-1)=0;qmN(A==10+i-1)=0;%mouth cell
%     %q1(A==10+i-1)=Q(i);
%     %for k=[N -1 1 -N];a=find(A==10+i-1);[a,q]=excludeboundarycell(k,N,M,a);q=q(A(q(a))==1);qm1(q)=Q(i);end %cell in front of the mouth
% %     
%     elseif directionQ(i)==2 %left
%     Ux(A==10+i-1)=Q(i)./h(A==10+i-1); Uy(A==10+i-1)=0;  %this will impose more uniform discharge over the oulet
%     q1(A==10+i-1)=0;qm1(A==10+i-1)=0;%mouth cell
%     %qN(A==10+i-1)=Q(i);
%     %for k=[N -1 1 -N];a=find(A==10+i-1);[a,q]=excludeboundarycell(k,N,M,a);q=q(A(q(a))==1);qmN(q)=Q(i);end %cell in front of the mouth
%     
%     elseif directionQ(i)==-2 %right
%     Ux(A==10+i-1)=Q(i)./h(A==10+i-1); Uy(A==10+i-1)=0;  %this will impose more uniform discharge over the oulet
%     q1(A==10+i-1)=0;qm1(A==10+i-1)=0;%mouth cell
%     %qmN(A==10+i-1)=-Q(i);
%     %for k=[N -1 1 -N];a=find(A==10+i-1);[a,q]=excludeboundarycell(k,N,M,a);q=q(A(q(a))==1);qN(q)=-Q(i);end %cell in front of the mouth
%    
%     end
% end

% figure
% subplot(2,1,1);
% imagesc(q1);caxis([0 82]);colormap('jet')
% subplot(2,1,2);
% imagesc(qm1);caxis([0 82]);colormap('jet')
% pause




%q=U.*h;












% q1=A*0+0.1;
% qm1=A*0+0.1;
% U=A*0+0.1;

%figure;imagesc(qm1);pause
% 
% figure;imagesc(P);
% 
% Dx=2*P-[P(1,:); P(1:end-1,:)]-[P(2:end,:); P(end,:)];
% Dy=2*P-[P(:,1) P(:,1:end-1)]-[P(:,2:end) P(:,end)];
% 
% 
% figure;
% subplot(3,1,1);imagesc(Dx);xlim([0 10]);ylim([80 120]);caxis([-1 1]/50);colormap('jet')
% subplot(3,1,2);imagesc(Dy);xlim([0 10]);ylim([80 120]);caxis([-1 1]/50);colormap('jet')
% subplot(3,1,3);imagesc(Dx+Dy);xlim([0 10]);ylim([80 120]);caxis([-1 1]/50);colormap('jet')
% 
% 
% 




% 
%%%check for divergence free

% Dx=[U1(1,:); U1(1:end-1,:)]-[Um1(2:end,:); Um1(end,:)];
% Dy=[UN(:,1) UN(:,1:end-1)]-[UmN(:,2:end) UmN(:,end)];
% 
% figure;
% subplot(3,1,1);imagesc(Dx);caxis([-1 1]/5);colormap('jet')%;xlim([0 10]);ylim([80 120])
% subplot(3,1,2);imagesc(-Dy);caxis([-1 1]/5);colormap('jet')
% subplot(3,1,3);imagesc(Dx+Dy);caxis([-1 1]/5);colormap('jet')
% pause







% % 
% % % P1=[P(:,2:end) P(:,end)];P2=[P(:,1) P(:,1:end-1) ];P3=[P(2:end,:); P(end,:)];P4=[P(1,:); P(1:end-1,:)];
% % % D1=[D(:,2:end) D(:,end)];D2=[D(:,1) D(:,1:end-1) ];D3=[D(2:end,:); D(end,:)];D4=[D(1,:); D(1:end-1,:)];
% % % A1=[A(:,2:end) A(:,end)];A2=[A(:,1) A(:,1:end-1) ];A3=[A(2:end,:); A(end,:)];A4=[A(1,:); A(1:end-1,:)];
% % % 
% % % pp1=P(A1==0);pp2=P(A2==0);pp3=P(A3==0);pp4=P(A4==0);
% % % P1(A1==0)=pp1;P2(A2==0)=pp2;P3(A3==0)=pp3;P4(A4==0)=pp4;
% % % 
% % % pp1=D(A1==0);pp2=D(A2==0);pp3=D(A3==0);pp4=D(A4==0);
% % % D1(A1==0)=pp1;D2(A2==0)=pp2;D3(A3==0)=pp3;D4(A4==0)=pp4;
% % % 
% % % DD1=(D+D1)/2;DD2=(D+D2)/2;DD3=(D+D3)/2;DD4=(D+D4)/2;
% % % 
% % %   Ux=0.5*((P-P1).*DD1 + (P2-P).*DD2)./h*dx;
% % %   Uy=0.5*((P-P3).*DD3 + (P4-P).*DD4)./h*dx;
% % %   Ux(A==2 | A==10)=2*Ux(A==2 | A==10);
% % %   Uy(A==2 | A==10)=2*Uy(A==2 | A==10);
% % % %G=cat(3,(P-P1).*DD1,(P2-P).*DD2);Ux=nanmean(G,3);
% % % %G=cat(3,(P-P3).*DD3,(P4-P).*DD4);Uy=nanmean(G,3);
% % 


