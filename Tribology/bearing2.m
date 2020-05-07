clc
clear all
Nx=input('nodes in x direction    :');
Nz=input('nodes in z direction    :');
rpm=input('enter the speed of runner   :');
eta0=0.06; % unit of viscosity is Pa-s
R1=0.075; % unit of radius is m
R2=0.0925*3.5; % unit of radius is m
Rm=(R1+R2)/2; % mean radius in m
U=(2*pi*Rm*rpm)/60; % sliding velocity in m/s
L=45*(pi/180)*Rm; % pad length corresponding to Rm
W= (R2-R1); % wirth of the pad
h2=3E-5;   % mim film thickness 30 micrometer
sh=3E-5;  % Shoulder height of fixed pad 20 micrometer
h1=(h2+sh); % inlet film thickness 
delxb=1/(Nx-1);
delzb=2/(Nz-1);
%...........................................

for i=1:1:Nx
    for k=1:1:Nz
        xb(i)=(i-1)*delxb;
        zb(k)=(k-1)*delzb;
        pb(i,k)=0;
        hb(i,k)=1+(sh/h2)*(1-xb(i));
        taub=1/hb(i,k);
    end
end
iter1=1;
error1=1;
while(error1>=1e-6)
    
     pbdifsum=0;
        pbsum=0;
%solution of Reynolds equation
for i=2:1:Nx-1
    for k=2:1:Nz-1
       
Term1= 3*(hb(i,k))^2;
Term2= (pb(i+1,k)-pb(i-1,k))/(2*delxb);
Term3= (hb(i+1,k)*hb(i-1,k))/(2*delxb);
Term4= Term1*Term2*Term3;
Term5= ((hb(i,k))^3)/delxb^2;
Term6= pb(i+1,k)+pb(i-1,k);
Term7= (2*L/W)^2;
Term8= (hb(i,k+1)-hb(i,k-1))/(2*delzb);
Term9= (pb(i,k+1)-pb(i,k-1))/(2*delzb);
Term10= ((hb(i,k))^3)/(delzb^2);
Term11= pb(i,k+1)+pb(i,k-1);
Term12=(3*(hb(i+1,k)-hb(i-1,k)))/(delxb);
Term13= Term4 + Term5*Term6;
Term14= 2*Term5;
Term15= Term7*Term1*Term8*Term9;
Term16= Term7*Term10*Term11;
Term17= 2*Term7*Term10;
 Term18= Term13 + Term15 + Term16 - Term12;
Term19= Term14 + Term17;
Term20= (Term18/Term19);
pp= Term20;
dif=pp-pb(i,k);
pb(i,k)=pp;
pbdifsum=pbdifsum+abs(dif);
pbsum=pbsum+pb(i,k);
taub(i,k)=(0.5*hb(i,k)*(p(i+1,k)-pb(i-1,k))/(2*delxb))+1/hb(i,k);


    end
end
error1=pbdifsum/pbsum;
end
for i=1:1:Nx
    for k=1:1:Nz
    p(i,k)=pb(i,k)*eta0*U*L/(h2*h2);
    end
end
[x,z]=meshgrid(i,k);
y=p(x,z);
surfc(p)
view(110,40)
zlabel('pressure')
ylabel('width')
xlabel('length')