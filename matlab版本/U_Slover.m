function [ u_new , apu, AE , AW , AN , AS ] = U_Slover( u,v,p,Nx,Ny,data,WSOR,rhoold,rho,mu,T)

dx = data.dx;
dy = data.dy;
hx = data.hx;
hy = data.hy;
theta=data.theta;
rp=data.rp;
cita=data.cita;
T0=data.Tf;
% dt = data.dt;

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
apu = ones(Nx,Ny);
u_new = u;
u_old = u;

% A = zeros(Nx*Ny,Nx*Ny);
% B = zeros(Nx*Ny,1);



for i = 3:Nx-1
    for j = 2:Ny-1
        % ip = (j-1)*Nx+i;
        rhoe=rho(i,j);
        rhow=rho(i-1,j);
        rhon=0.25*(rho(i,j)+rho(i,j+1)+rho(i-1,j+1)+rho(i-1,j));
        rhos=0.25*(rho(i,j)+rho(i-1,j)+rho(i-1,j-1)+rho(i,j-1));
        mue=mu(i,j);
        muw=mu(i-1,j);
        mun=0.25*(mu(i,j)+mu(i,j+1)+mu(i-1,j+1)+mu(i-1,j));
        mus=0.25*(mu(i,j)+mu(i-1,j)+mu(i-1,j-1)+mu(i,j-1));
        
        
        % Pe=0.5*rhoe*( u_old(i+1,j)   + u_old(i,j)   )*dx(i,j)/mue;
        % Pw=0.5*rhow*( u_old(i-1,j)   + u_old(i,j)   )*dx(i,j)/muw;
        % Pn=0.5*rhon*( v(i,j+1)   + v(i-1,j+1) )*dy(i,j)/mun;
        % Ps=0.5*rhos*( v(i,j)     + v(i-1,j) )*dy(i,j)/mus;

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
        
        Fe = 0.5*rhoe*( u_old(i+1,j)   + u_old(i,j)   )*dy(i,j);
        Fw = 0.5*rhow*( u_old(i-1,j)   + u_old(i,j)   )*dy(i,j);
        Fn = 0.5*rhon*( v(i,j+1)   + v(i-1,j+1) )*dx(i,j);
        Fs = 0.5*rhos*( v(i,j)     + v(i-1,j) )*dx(i,j);

        % ae = -0.5*Fe + De;
        % aw =  0.5*Fw + Dw;
        % an = -0.5*Fn + Dn;
        % as =  0.5*Fs + Ds;
        % ap = ae + aw + an + as + (Fe - Fw) + (Fn - Fs);
        % me=0.5*(abs(Fe)-Fe);
        % mw=0.5*(abs(Fw)+Fe);
        % mn=0.5*(abs(Fn)-Fn);
        % ms=0.5*(abs(Fs)+Fs);
        me=max(-Fe,0);
        mw=max( Fw,0);
        mn=max(-Fn,0);
        ms=max( Fs,0);
        cds_upw=0.5*(Fe*(u(i,j)+u(i+1,j))-Fw*(u(i,j)+u(i-1,j))+Fn*(u(i,j)+u(i,j+1))-Fs*(u(i,j)+u(i,j-1)))...
            -((me+mw+mn+ms)*u(i,j)-(me*u(i+1,j)+mw*u(i-1,j)+mn*u(i,j+1)+ms*u(i,j-1)));

        Dp = (p(i-1,j)-p(i,j))*dy(i,j);
        Sc=mu(i,j)*dy(i,j)/rp(i,j)*(v(i,j+1)+v(i,j)-v(i-1,j)-v(i-1,j+1))...
            +rhoold(i,j)*(-cita*(T(i,j)-T0))*(9.8)*sin(theta(i,j))*dx(i,j)*dy(i,j);
        Sp=-rhoold(i,j)*0.25*(v(i,j+1)+v(i,j)+v(i-1,j)+v(i-1,j+1))*(theta(i+1,j)-theta(i,j))*dy(i,j)...
            -mu(i,j)*dy(i,j)*(theta(i+1,j)-theta(i,j))/rp(i,j);
        Spmax=max(-Sp,0);
        
        ae = me + De;
        aw = mw + Dw;
        an = mn + Dn;
        as = ms + Ds;
        if j==2   
            as = 0; %because south area=0
        end
        % if j==Ny-1
        %     an = mn + 2*Dn;
        % end
        % ap0= rhoold(i,j)*dy*dx/dt;
        ap = ae + aw + an + as + Spmax ;%+ ap0;% + Fe - Fw + Fn - Fs;
        ap = ap/WSOR.alpha_u;

        
        AP(i,j) =  ap;
        AE(i,j) = -ae;
        AW(i,j) = -aw;
        AN(i,j) = -an;
        AS(i,j) = -as;
        SP(i,j) =  Dp+Sc-WSOR.beta*cds_upw+(Sp+Spmax)*u_old(i,j);%+ap0*u_old(i,j);%
        apu(i,j) = ap;%ap;

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
    for i = 3:Nx-1
        for j = 2:Ny-1
            Be = AE(i,j)*u_new(i+1,j);
            Bw = AW(i,j)*u_new(i-1,j);
            Bn = AN(i,j)*u_new(i,j+1);
            Bs = AS(i,j)*u_new(i,j-1);
            u_new(i,j) = (1-WSOR.alpha_u)*u_old(i,j)+...
                (SP(i,j)-Be-Bw-Bn-Bs)/AP(i,j);
            % u_new(i,j) = WSOR.alpha_u*u_new(i,j)+(1.0-WSOR.alpha_u)*u_old(i,j);
        end
    end
    
    GError = sqrt((sum(sum((u_new-u_old).^2.0)))/(Nx*Ny));
    
    u_old = u_new;
    
    iter = iter+1;
    
    if( iter > MaxI )
        break;
    end
    if( iter < MinI )
        GError = 1.0;
    end
    
end


end

