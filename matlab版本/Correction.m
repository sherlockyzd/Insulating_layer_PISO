function [ u_new,v_new,p_new,src ] = Correction( u,apu,v,apv,p,pprime,NxU,NyU,NxV,NyV,NxP,NyP,WSOR,data,rhoold,rho)

dx = data.dx;
dy= data.dy;
% dt = data.dt;
% hx = data.hx;
% hy = data.hy;
% vx = data.vx;
% vy = data.vy;
% U_lid = data.U_lid;
% rho = data.rho;
% mu = data.mu;


u_new = u;
v_new = v;
p_new = p;
source = zeros(NxP,NyP);


for i = 2:NxP-1
    for j = 2:NyP-1
        p_new(i,j) = p(i,j) + WSOR.alpha_p*pprime(i,j);    %修正P
    end
end

for i = 3:NxU-1
    for j = 2:NyU-1
        uprime = (dy(i,j)/apu(i,j))*( pprime(i-1,j) - pprime(i,j) ); %修正u，uprime就是u一撇
        u_new(i,j) = u(i,j) + uprime;
    end
end

for i = 2:NxV-1
    for j = 3:NyV-1
        vprime = (dx(i,j)/apv(i,j))*( pprime(i,j-1) - pprime(i,j) );
        v_new(i,j) = v(i,j) + vprime;
    end
end

for i = 2:NxP-1
    for j = 2:NyP-1
        rhoe=0.5*(rho(i,j)+rho(i+1,j));
        rhow=0.5*(rho(i,j)+rho(i-1,j));
        rhon=0.5*(rho(i,j)+rho(i,j+1));
        rhos=0.5*(rho(i,j)+rho(i,j-1));
        
        um = dy(i,j)*( rhoe*u_new(i+1,j) - rhow*u_new(i,j) );
        vm = dx(i,j)*( rhon*v_new(i,j+1) - rhos*v_new(i,j) );
        source(i,j) = -um-vm;%+(rhoold(i,j)-rho(i,j))*dx(i,j)*dy(i,j)/dt;%
    end
end
src =sum(sum(source.^2))/(NxP*NyP);

end

