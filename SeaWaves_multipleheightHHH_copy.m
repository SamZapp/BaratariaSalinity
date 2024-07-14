function [Umi,Tpseai,HS,F,PWi,QsWslope_seai]=SeaWaves(A,hw,angle,hwSea_lim,range,wind,VEG,ndir,dx,z,msl,tempdeltaMSL,ws1,fTide,extraHseaimposed,addextrafetch,extrafetch);

Nh=10;


[N,M]=size(hw);
Lbasin=0;%1000/dx;
Fetchlim=0;%max(50,dx*2);%dx*2;%600;%dx*2*10;


MASK=0*A+1;MASK(hw<=hwSea_lim | A==0 | VEG==1)=0;

%The standard way
if addextrafetch==0
F=calculatefetch(MASK,ndir,dx,angle);
elseif addextrafetch==1;
F=calculatefetchWITHEXTRAS(MASK,ndir,dx,angle,extrafetch,Lbasin,MASK);
end

Fo=F;F(Fo<=Fetchlim)=0;

%diffuse the fetch field    
alphadiffusefetch=0.01;   %messo 10 for the VCR wave validation 10;%0;%%%QUESTO ERA 1 FINO AD APRILE 23 2018!!!!!
F=diffusefetchPERPEND(MASK,F,alphadiffusefetch,dx,angle); 
F(Fo<=Fetchlim | MASK==0)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=find(Fo>Fetchlim & hw>hwSea_lim & F>0 & MASK==1);%h>dlo & %a=find(Fo>dx*2);%h>dlo & %a=find(h>dlo);

PWi=0;
QsWslope_seai=0;

z=z-tempdeltaMSL;%temporary shift to change MSL at ever tide, trick!!!
for i=1:Nh

dHW=max(0,-z+msl+range/2);%water depth at MHW
h=max(0,dHW-range*(i-1)/(Nh-1));
D=h(a);Ff=F(a);
  

TP=0*h;HS=0*h;
[Hs,Tp]=YeV(Ff,wind,D);%[Hs,Tp]=YeV(Ff,wind,min(3,D));  %TRUCCO PER EVITARE LARGE WAVES IN CHANELS
HS(a)=Hs;TP(a)=Tp;TP(TP==0)=1;
 
if extraHseaimposed==1
HS=getextraHsea(HS,h,MASK,dx,N,M);
end

HS(h<hwSea_lim)=0;
TP(h<hwSea_lim)=2;

kwave=0*h;kk=wavek(1./TP(a),h(a));kwave(a)=kk;
kwave(kwave==0 | h<hwSea_lim)=1;

Um=pi*HS./(TP.*sinh(kwave.*h));
cg=(2*pi./kwave./TP)*0.5.*(1+2*kwave.*h./(sinh(2*kwave.*h)));
PW=cg*1030*9.8.*HS.^2/16;

Um(h<hwSea_lim)=0;%
PW(h<hwSea_lim)=0;

Um(MASK==0)=0;
PW(MASK==0)=0;
Umi(:,:,i)=Um;
Tpseai(:,:,i)=TP;
PWi=PWi+PW;

rhos=2650;ss=1.5728;
QsWslope_sea=WaveSedimentTransport(HS,h,kwave,rhos,N,M,TP,dx,ss,ws1,hwSea_lim,NaN);
QsWslope_sea(HS==0)=0;
QsWslope_seai=QsWslope_seai+QsWslope_sea;

end

PWi=PWi/Nh;
QsWslope_seai=QsWslope_seai/Nh;

