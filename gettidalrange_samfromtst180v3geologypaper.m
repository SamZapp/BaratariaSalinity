function [Trange90]=gettidalrange(Ttide90,Trange90_o,Cb,Cv,dx,x,zG,A);
%updated gettidalrange 11/18/22: h_n_eq has a lower bound to account for
%sub-grid size channels, b has a very small lower bound (for numeric
%stability, 3 tidal cycles are run to ensure that there's sufficient time
%to flush water setup in microtidal environments, water setup is now
%set by average topography


%zG=flip(zG,1);
n_add=0.02%0.023%+0.006 %INCREASE MANNING JUST HERE
zG(A==0)=NaN;
zG=zG(end:-1:1,:);
nG=Cb+zG*0+n_add;
nG(zG>=0)=Cv;
nG(isnan(zG))=NaN;
[N,M]=size(zG);

Ttide90=Ttide90*24*3600;
%T=12*3600;%tidal period
%r=0.7;%tidal range
%n=0.02;%manning
%bc=100;%channel width
%b=bc;%5000;%total width (channel + marsh/mudalt)
%b=500;%total width (channel + marsh/mudalt)
%ho=2;

%dx=50;
%x=[0:dx:20*1000];
%x=[1:5];
%N=length(x);

dt=3600*0.5;
t=[0:dt:3*Ttide90];
Y=Trange90_o/2*cos(t*2*pi/Ttide90+pi/2);
eta=x'*Trange90_o/2*10^-5+Y(1);%initial water level slope, made it sloping out
%eta=x'*0+Y(1);%inital = flat. Make the transient longer

%eta=mean(zG,2); %problem here?


G=zeros(length(x),length(t));
%U=zeros(length(x),length(t));
%figure
for j=1:length(t);
ii=[];jj=[];AA=[];

    
    %h=ho+eta;%boundary condition
    eta2d=eta*ones(1,M);
    hG=max(0,eta2d-zG);%.*(eta2d>zextra);
    
    
    deta=([eta(1); eta(1:end-1)]-[eta(2:end); eta(end)])/2/dx;deta(1)=deta(1)*2;
    %b=dx*M;
    %b=dx*sum((eta2d>zG),2);b=max(500,b);%set minimum bound* removed min bound for adjustment on 11/9/22
    %b=bc+(eta+r/2)*100;
    %Do=h.^(5/3)*bc./b/n./max(10^-10,sqrt(abs(deta)));
%     b_SGCadj=1*2*(sum(hG==0,2)*dx)/50; %subgrid channel adjustment, assumes a 1m x2m chan every 50m
%     b=b+b_SGCadj; %%WRONG!!
    
    h_n_eq=(nansum((hG.^(5/3)./nG),2)*dx).^(3/5);
    
    A_R=2; %aspect ratio
    wc=2; %width of subgrid channel (make it a function of dx?)
    dc=-wc/A_R; %aspect ratio
    hc=max(0,eta2d-dc); hc=hc(:,1); %
    hc=max(hc,abs(dc)); %technically should comment this out, but needed to prevent h_n_eq from reaching zero
    chanfreq=200; %meters/channel
    %heq=(sum((hG.^(5/3)),2)*dx./b).^(3/5);
    h_n_eq_SGC=((hc.^(5/3)./(Cb+n_add).*(M*dx/chanfreq)*wc).^(3/5));
    %h_n_eq_SGC=((hc.^(5/3)./(Cb+n_add)).^(3/5))*(M*dx/chanfreq)*wc; %units: put in m
    h_n_eq=max(h_n_eq,h_n_eq_SGC);
    
    %Do=sum((hG.^(5/3)),2)*dx./b/n./max(10^-10,sqrt(abs(deta)));
%     b=b+(M*dx/chanfreq)*wc
    Do=h_n_eq.^(5/3)./max(10^-10,sqrt(abs(deta)));
    %Do=heq.^(5/3)./n./max(10^-10,sqrt(abs(deta)));
    
    %Do=(sum((hG.^(5/3)),2)*dx)./b/n./max(10^-10,sqrt(abs(deta)));
    D  =Do/dx^2*dt;
    %D=D*0+1;
 
    b=dx*sum((eta2d>zG),2);b=max(b,0.01);%b=max(1*(M*dx/chanfreq)*wc,b);%b=max(hc*(M*dx/chanfreq)*wc,b); %technically should be the commented part but can cause singular matrix
%OLD METHOD IMPLICIT< NEGELCTS SOME DERIVATIVE
%    value=ones(1,N)+2*[0 D(2:end)'];   
%    ii=[ii [1:N]];   jj=[jj [1:N]];   AA=[AA value];
% 
%    value=-D(2:end)';
%    ii=[ii [2:N]];   jj=[jj [[3:N] N]];   AA=[AA value];
%    
%    value=-D(2:end)';
%    ii=[ii [2:N]];   jj=[jj [[1:N-1]]];   AA=[AA value];
   
   
   
%NEW METHOD IMPLICIT   
   Dleft=(D(2:end)+D(1:end-1))'/2;
   Drigth=(D(2:end)+[D(3:end); D(end)])'/2;
    
   value=ones(1,N)+[0 Drigth+Dleft]./b';  
   ii=[ii [1:N]];   jj=[jj [1:N]];   AA=[AA value];

   value=-Dleft./b(2:end)';
   ii=[ii [2:N]];   jj=[jj [[1:N-1]]];   AA=[AA value];
   
   value=-Drigth./b(2:end)';
   ii=[ii [2:N]];   jj=[jj [[3:N] N]];   AA=[AA value];
   
   
   
   
   
   I = sparse(ii,jj,AA);
   

   eta(1)=Y(j);
   eta=I\eta;
    

   G(:,j)=eta;
end
 
tiderange=max(G(:,2*floor(length(t)/3):end)')-min(G(:,2*floor(length(t)/3):end)');
%tiderange=max(G(:,floor(length(t)/2):end)')-min(G(:,floor(length(t)/2):end)');
Trange90_1D=tiderange;
Trange90=repmat(flip(Trange90_1D'),1,length(zG(1,:))); %pad array to make 2d
%hightide=max(G(:,floor(length(t)/2):end)');
%lowtide=min(G(:,floor(length(t)/2):end)');

%added line below 4/28/22
%Trange90(Trange90<0.02)=0.02;%min limit on trange

% figure
% plot(b)
% 
% figure
% imagesc(eta2d)
% pause