function awVAR=getvariableedgeerosion(A,aw,SALTC);

awVAR=A*0+aw;
awVAR(SALTC>0.5 & SALTC<5)=aw*10;
