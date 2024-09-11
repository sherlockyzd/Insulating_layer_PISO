function [ v_new , apv, AE , AW , AN , AS ] = V_Slover( u,v,p,Nx,Ny,data,WSOR,rhoold,rho,mu,T)

dx = data.dx;
dy = data.dy;
hx = data.hx;
hy = data.hy;
theta=data.theta;
rp=data.rp;
cita=data.cita;
T0=data.Tf;
% dt = data.dt;
% vx = data.vx;
% vy = data.vy;
% U_lid = data.U_lid;
% rho = data.rho;
% mu = data.mu;

MaxG = WSOR.MaxG;
MaxI = WSOR.MaxU;
MinI = WSOR.MinU;

AP = zeros(Nx,Ny);
AE = zeros(Nx,Ny);
AW = zeros(Nx,Ny);
AN = zeros(Nx,Ny);
AS = zeros(Nx,Ny);
SP = zeros(Nx,Ny);
apv = ones(Nx,Ny);
v_new = v;
v_old = v;

% A = zeros(Nx*Ny,Nx*Ny);
% B = zeros(Nx*Ny,1);


for i = 2:Nx-1
    for j = 3:Ny-1
        % ip = (j-1)*Nx+i;
        rhoe=0.25*(rho(i,j)+rho(i+1,j)+rho(i,j-1)+rho(i+1,j-1));
        rhow=0.25*(rho(i,j)+rho(i-1,j)+rho(i-1,j-1)+rho(i,j-1));
        rhon=rho(i,j);
        rhos=rho(i,j-1);

        mue=0.25*(mu(i,j)+mu(i+1,j)+mu(i,j-1)+mu(i+1,j-1));
        muw=0.25*(mu(i,j)+mu(i-1,j)+mu(i-1,j-1)+mu(i,j-1));
        mun=mu(i,j);
        mus=mu(i,j-1);

        
        % Pe=0.5*rhoe*( u(i+1,j)   + u(i+1,j-1)   )*dx(i,j)/mue;
        % Pw=0.5*rhow*( u(i,j) + u(i,j-1) )*dx(i,j)/muw;
        % Pn=0.5*rhon*( v_old(i,j+1)   + v_old(i,j)   )*dy(i,j)/mun;
        % Ps=0.5*rhos*( v_old(i,j-1)   + v_old(i,j)   )*dy(i,j)/mus;
        % De0 = (mue*dy(i,j)*hx(i,j));
        % Dw0 = (muw*dy(i,j)*hx(i,j));
        % Dn0 = (mun*dx(i,j)*hy(i,j));
        % Ds0 = (mus*dx(i,j)*hy(i,j));
        % De=De0*max(1-0.1*abs(Pe)^5,0);
        % Dw=Dw0*max(1-0.1*abs(Pw)^5,0);
        % Dn=Dn0*max(1-0.1*abs(Pn)^5,0);
        % Ds=Ds0*max(1-0.1*abs(Ps)^5,0);
        % De=De0*max(1-0.5*Pe,0);
        % Dw=Dw0*max(1-0.5*Pw,0);
        % Dn=Dn0*max(1-0.5*Pn,0);
        % Ds=Ds0*max(1-0.5*Ps,0);
        De = (mue*dy(i,j)*hx(i,j));
        Dw = (muw*dy(i,j)*hx(i,j));
        Dn = (mun*dx(i,j)*hy(i,j));
        Ds = (mus*dx(i,j)*hy(i,j));
        
        Fe = 0.5*rhoe*( u(i+1,j)   + u(i+1,j-1)   )*dy(i,j);
        Fw = 0.5*rhow*( u(i,j) + u(i,j-1) )*dy(i,j);
        Fn = 0.5*rhon*( v_old(i,j+1)   + v_old(i,j)   )*dx(i,j);
        Fs = 0.5*rhos*( v_old(i,j-1)   + v_old(i,j)   )*dx(i,j);

        % me=0.5*(abs(Fe)-Fe);
        % mw=0.5*(abs(Fw)+Fe);
        % mn=0.5*(abs(Fn)-Fn);
        % ms=0.5*(abs(Fs)+Fs);
        me=max(-Fe,0);
        mw=max( Fw,0);
        mn=max(-Fn,0);
        ms=max( Fs,0);
        cds_upw=0.5*(Fe*(v(i,j)+v(i+1,j))-Fw*(v(i,j)+v(i-1,j))+Fn*(v(i,j)+v(i,j+1))-Fs*(v(i,j)+v(i,j-1)))-...
        ((me+mw+mn+ms)*v(i,j)-(me*v(i+1,j)+mw*v(i-1,j)+mn*v(i,j+1)+ms*v(i,j-1)));

        Dp = (p(i,j-1)-p(i,j))*dx(i,j);
        Sc=mu(i,j)*dy(i,j)/rp(i,j)*(u(i,j-1)+u(i,j)-u(i+1,j)-u(i+1,j-1))...
            +rhoold(i,j)*(0.25*(u(i,j-1)+u(i,j)+u(i+1,j)+u(i+1,j-1)))^2*(theta(i+1,j)-theta(i,j))*dy(i,j)...
            -rhoold(i,j)*(-cita*(T(i,j)-T0))*(9.8)*cos(theta(i,j))*dx(i,j)*dy(i,j);
        Sp=-mu(i,j)*dy(i,j)*(theta(i+1,j)-theta(i,j))/rp(i,j);
        Spmax=max(-Sp,0);
        
        ae = me + De;
        aw = mw + Dw;
        an = mn + Dn;
        as = ms + Ds;
        if i==2   
            aw = 0;%because zeroGradient
        end
        if i==Nx-1
            ae = 0;%because zeroGradient
        end
        % ap0= rhoold(i,j)*dy(i,j)*dx(i,j)/dt;
        ap = ae + aw + an + as  + Spmax;%+ ap0;% + Fe - Fw + Fn - Fs;
        ap = ap/WSOR.alpha_v;

        AP(i,j) =  ap;
        AE(i,j) = -ae;
        AW(i,j) = -aw;
        AN(i,j) = -an;
        AS(i,j) = -as;
        SP(i,j) =  Dp+Sc-WSOR.beta*cds_upw+(Sp+Spmax)*v_old(i,j);%+ap0*v_old(i,j);%

        apv(i,j) = ap; %ap;
        
        % A(ip,ip)    = AP(i,j);
        % A(ip,ip+1)  = AE(i,j);
        % A(ip,ip-1)  = AW(i,j);
        % A(ip,ip+Nx) = AN(i,j);
        % A(ip,ip-Nx) = AS(i,j);
        % B(ip,1)     = SP(i,j);
        
    end
end

GError = 1.0;
iter = 1;
while ( GError > MaxG )

    for i = 2:Nx-1
        for j = 3:Ny-1
            Be = AE(i,j)*v_new(i+1,j);
            Bw = AW(i,j)*v_new(i-1,j);
            Bn = AN(i,j)*v_new(i,j+1);
            Bs = AS(i,j)*v_new(i,j-1);
            v_new(i,j) = (1-WSOR.alpha_v)*v_old(i,j)+...
                (SP(i,j)-Be-Bw-Bn-Bs)/AP(i,j);
            % v_new(i,j) = WSOR.alpha_v*v_new(i,j)+(1.0-WSOR.alpha_v)*v_old(i,j);
        end
    end
    
    GError = sqrt((sum(sum((v_new-v_old).^2.0)))/(Nx*Ny));
    
    v_old = v_new;
    
    iter = iter+1;
    
    if( iter > MaxI )
        break;
    end
    if( iter < MinI )
        GError = 1.0;
    end
    
end


end

