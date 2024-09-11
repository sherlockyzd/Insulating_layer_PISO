function [ u_new,v_new ] = Correction2( u0,u,apu,AEu,AWu,ANu,ASu,v0,v,apv,AEv,AWv,ANv,ASv,NxU,NyU,NxV,NyV)

% dx = data.dx;
% dy = data.dy;
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

up=u-u0;
vp=v-v0;
for i = 3:NxU-1
    for j = 2:NyU-1
        Be = -AEu(i,j)*up(i+1,j);
        Bw = -AWu(i,j)*up(i-1,j);
        Bn = -ANu(i,j)*up(i,j+1);
        Bs = -ASu(i,j)*up(i,j-1);
        u_new(i,j) =  u(i,j)+(Be+Bw+Bn+Bs)/apu(i,j);
    end
end

for i = 2:NxV-1
    for j = 3:NyV-1
        Be = -AEv(i,j)*vp(i+1,j);
        Bw = -AWv(i,j)*vp(i-1,j);
        Bn = -ANv(i,j)*vp(i,j+1);
        Bs = -ASv(i,j)*vp(i,j-1);
        v_new(i,j) =  v(i,j)+(Be+Bw+Bn+Bs)/apv(i,j);
    end
end

end

