function [N,M,dx,A,AW,Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3,Active,x,y,msl,SPCLcell,B,directionQ]=initializegeometry(P);


dx=30;
%a= imread('BU1932bathy','tiff');
a= imread('BU2007bathy','tiff');a=double(a);
[aa] = imread('BU1932','tiff');aa=double(aa);aa(aa==2)=0;
%[aa] = imread('BU2015','tiff');for i=1:1980;aa(1:478+floor(i*0.2),i)=aa(1:478+floor(i*0.2),i)-1700;end;aa=aa>6300;

%best guess for unknow bathymetry
%a(aa==1)=0.1;%min(a(aa==1),0.3);Marsh
a(aa==1)=0.3;%min(a(aa==1),0.3);marsh+0.1m
a(aa==0)=min(a(aa==0),-0.8);

%a=min(0.3,a);

z=a;

%figure;imagesc(z);pause
% figure
% imagesc(z);caxis([-2 1])
% axis([800 2000 1000 2500 ])
% cch = roipoly;save cch cch;pause
% pause
load cch;z(cch==1)=-10;


%GIWW BARTAT
%figure;imagesc(z);caxis([-5 1]);colormap('jet');axis([1 1200 500 1100]);BW_1 = roipoly;save BW_1 BW_1;pause
load BW_1;z(BW_1==1)=-5;


%figure;imagesc(z);caxis([-5 1]);colormap('jet');axis([1 1200 500 1100]);BW_2 = roipoly;save BW_2 BW_2;pause
%load BW_2;z(BW_1==1)=-5;

% %GIWW BARTAT-upper Rigolettes
% %figure;imagesc(z);caxis([-5 1]);colormap('jet');axis([2600 3600 1870 2500]);BW3 = roipoly;save BW3 BW3;pause
% load BW3;z(BW3==1)=-10;




% figure
% imagesc(z);caxis([-2 1])
% axis([800 2000 2000 2500 ])
% c8 = roipoly;save c8 c8;pause
% pause

%barrier islands
load c7;z(c7==1)=9990;
load c8;z(c8==1)=9990;


%port fourcnohn corner
load c6
z(c6==1)=9990;


load c1
load c2
load c3
z(c1==1 | c2==1)=9990;
z(end-200:end,1:400)=9990;
z(c3==1)=9992;


%figure;imagesc(z);caxis([-1 1]);pause

%cut the bottom
z=z(1:end-100,:);

%resamling
%z=z(1:10:end,1:10:end);dx=dx*10;
%z=z(1:2:end,1:2:end);dx=dx*2;
z=z(1:4:end,1:4:end);dx=dx*4;


[N,M]=size(z);
A=ones(N,M);


A(z==9990)=0;
A(z==9992)=2;

z=-z;
%z=z*0+5;

rivermouthfront=[];

% mouthW=1; %2
% A(1:2,1:M/2-1-mouthW)=0;%create a concreate wall on the sides
% A(1:2,M/2+1+mouthW:end)=0;%create a concreate wall on the sides
% A(1,M/2-mouthW:M/2+mouthW)=10;
% %these are the cells in front of the river mouth
% S=A*0;S(2,M/2-mouthW:M/2+mouthW)=1;rivermouthfront=find(S==1);clear S;

SPCLcell=struct;
SPCLcell.rivermouthfront=rivermouthfront;

%bathymetry
%initial profile. ELEVATION AT MSL
x=[0:N-1]*dx/1000;y=[0:M-1]*dx/1000;


z(A==0)=NaN;
z(A==2)=50;

%%%%%%%%%%%%%%%%%%

% 
% %%nearshoremask
% dM=120*250/dx;
% for i=M-dM:M
%     A(end-min(15000/dx,floor(0.75*(dM-(M-i)))):end,i)=2;
% end
% 
% %%New Orelans
% dM=170*250/dx;
% for i=M-dM:M
%     A(1:min(14000/dx,floor(1*(dM-(M-i)))),i)=0;
% end
% 
% %%Deltafarms
% dN=53000/dx;
% for i=32500/dx:dN
%     A(i,1:min(104000/dx,floor(0.5*(dN-(N-i)))))=0;
% end
% 
% %%Marina Near mid diversion
% for i=13500/dx:24000/dx
%     A(i,M-25500/dx:M)=0;
% end
% dN=33000/dx;
% for i=24000/dx:dN+1000/dx
%     A(i,end+12000/dx+(floor(2*(dN-(N-i)))):end)=0;
% end



A(:,1)=0;
A(:,end)=0;
A(1,:)=0;

% %Left side
% A(M/2-1:M/2+1,1)=10;
% z(A==10)=P.hmouth(1);

% 


% %Right side
% A(M/2-1:M/2+1,end)=11;
% z(A==11)=P.hmouth(2);
% %z(1:5,N/2-1:N/2+1)=P.hmouth(1);

%north side-DAVIS POND
A(1,96:96)=10;
z(A==10)=P.hmouth(1);
z(1:5,96:96)=P.hmouth(1);

%GIWW
A(255,1)=11;
z(A==10)=P.hmouth(2);
z(255,1:2)=P.hmouth(2);


%Mid barataria
%A(198+14-1,308-1)=12;
%A(211,307)=12;
A(207,302)=12;
A(207,301)=1;
z(207,302-5:302)=P.hmouth(3);
z(208,302-10:302-5)=P.hmouth(3);
z(209,302-15:302-10)=P.hmouth(3);
z(210,302-20:302-15)=P.hmouth(3);
z(211,302-25:302-20)=P.hmouth(3);
z(A==12)=P.hmouth(3);
%z(211-1:211+1,307-1:307+1)=P.hmouth(3);

%Northwest corner of Lake Slavador
A(158,1)=13;
%A(158:159,1)=13;
z(A==10)=P.hmouth(4);
z(158:159,1:2)=P.hmouth(4);
%figure;imagesc(z);pause


%%GIW ON THE EAST SIDE< FROM NOLA
%A(198+14-1,308-1)=12;
%A(211,307)=12;
A(85,236)=14;
A(85,234:235)=0;
z(86:89,236)=P.hmouth(5);
z(A==14)=P.hmouth(5);

z(102,196:215)=5;
z(100,215:221)=5;
z(97:100,221)=5;
z(92:97,224)=5;
z(92,224:236)=5;



%figure;imagesc(A);pause

rn=2*(rand(N,M)-0.5)*0.05;
z=z+rn;

msl=0;
zbedo=-z-msl;
Active=zbedo<2;

G=A;
G=[A(2:end,:); A(end,:)];
G(A==1 & G==2 & z<0)=999;G(G~=999)=0;
G=max(G,[G(2:end,:); G(end,:)]);
A(G==999)=0;
%figure;imagesc(G);pause
%figure;imagesc(A);pause

B=A*0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Yb,Y1,Y2,Y3,zb,zs,plyr,flyr1,flyr2,flyr3,flyrb1,flyrb2,flyrb3]=initializestratigraphy_3sediments(z,N,M,P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%wave boundary conditions
%the lateral boundaries for waves
%postivie is left boundary; negative is right boudndary; 
%1 is no-gradient in wave. THIS IS ALSO no-gradient in wave-indcued lateral
%sediment transport
%2 is zero wave height . Also implies no sand transport at inlet
AW=A*0;%the internal points
%right boundary
AW(:,end)=-1;
%left boundary
AW(:,1)=1;
%AW(1:50,1)=1;
%%%%%%%%%%%%%%%%%%%%%
