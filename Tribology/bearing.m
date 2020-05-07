%This program is prepared by Suman Kumar                                                              
clc;%clear command window
clear all;%graphic clear
Nx=input(' Nodes in X-Direction   :');
Nz=input(' Nodes in Z-Direction   :');
rpm=input('Speed of Runner        :');
eta0=0.06;

%low viscous bearing SAE 30 40 50 
%unit of  viscocoty is Pa-s;
% R1=0.075;
% R2=0.0925;                  %Radius distance in m
% Rm=(R1+R2)/2                %mean Radius0.0837
U= (0.0088*rpm);        %sliding Velocity (in m/s)0.0837
L=0.0657*rpm;           %2piN/60
W=0.0175;
h2=3E-6;                    %Minimum film thickness
sh=3E-6;
h1=h2+sh;
delxb=1/(Nx-1);
delzb=2/(Nz-1);
%*******************************************************************

for i=1:1:Nx
    for k=1:1:Nz
        xb(i)=(i-1)*delxb;
        zb(k)=-1+(k-1)*delzb;
        pb(i,k)=0;
        hb(i,k)=1+(sh/h2)*[1-xb(i)];
    end
end

iter1=1;
error1=1;
 while(error1>=1e-6)
     pbdifsum=0;
     pbsum=0;
     for i=2:1:Nx-1
         for k=2:1:Nz-1 

term1= 3*(hb(i,k))^2;
term2= (pb(i+1,k)-pb(i-1,k))/(2*delxb);
term3= (hb(i+1,k)-hb(i-1,k))/(2*delxb);
term4= term1*term2*term3;
term5= ((hb(i,k))^3)/delxb^2;
term6= pb(i+1,k)+pb(i-1,k);
term7= (2*L/W)^2;
term8= (hb(i,k+1)-hb(i,k-1))/(2*delzb);
term9= (pb(i,k+1)-pb(i,k-1))/(2*delzb);
term10= ((hb(i,k))^3)/(delzb^2);
term11= pb(i,k+1)+pb(i,k-1);
term12=(3*(hb(i+1,k)-hb(i-1,k)))/(delxb);
term13= term4 + term5* term6;
term14= 2*term5;
term15= term7*term1*term8*term9;
term16= term7*term10*term11;
term17= 2*term7*term10;
term18= term13+term15+term16-term12;
term19= term14+term17;
term20= (term18/term19);
pp= term20;
dif=(pp-pb(i,k));
pb(i,k)=pp;
pbdifsum=pbdifsum+abs(dif);
pbsum=pbsum+pb(i,k);
       end
     end
    error1=pbdifsum/pbsum;
 end
 
 for i=1:1:Nx
    for k=1:1:Nz        
        p(i,k)= pb(i,k)*eta0*U*1/(h2*h2)%%%check;
    end
 end
 hold on
 [X,Z]=meshgrid(i,k);
 surfc(p);
 view(110,45);
 zlabel('width');
 ylabel('pressure');
 xlabel('Length');
%************************************************************************
%solution for Renoylds Equation

