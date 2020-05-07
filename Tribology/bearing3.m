clear all
close all
for k=1:1:Nz
    Wsum1=p(i,k)+p(Nx,k)
    Tsum1=taub(i,j)+taub(Nx,j);
    Wsum2=0.0;
    Tsum2=0.0
    for i=2:2:Nx-1
        Wsum2=Wsum2+pb(i,3);
        Tsum2=Tsum2+pb(