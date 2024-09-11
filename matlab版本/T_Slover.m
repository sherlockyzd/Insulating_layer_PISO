function [ T_new,GError,iter ] = T_Slover( u,v,T0,Nx,Ny,data,rho,kappa,WSOR)

dx = data.dx;
dy = data.dy;
hx = data.hx;
hy = data.hy;

% U_lid = data.U_lid;
% rho = data.rho;
% kappa = data.kappa;
Nu=size(u);

MaxG = WSOR.MaxTG;
MaxI = WSOR.MaxT;
MinI = WSOR.MinT;

AP = zeros(Nx,Ny);
AE = zeros(Nx,Ny);
AW = zeros(Nx,Ny);
AN = zeros(Nx,Ny);
AS = zeros(Nx,Ny);
SP = zeros(Nx,Ny);
% apu = ones(Nx,Ny);
% u_new = u;
% u_old = u;
T_new = T0;
T_old = T0;

% A = zeros(Nx*Ny,Nx*Ny);
% B = zeros(Nx*Ny,1);




for i = 2:Nx-1
    for j = 2:Ny-1
      
        if i<Nu(1) && j<Nu(2)
            rhoe=0.5*(rho(i,j)+rho(i+1,j));
            rhow=0.5*(rho(i,j)+rho(i-1,j));
            rhon=0.5*(rho(i,j)+rho(i,j+1));
            rhos=0.5*(rho(i,j)+rho(i,j-1));
            Fe = rhoe*u(i+1,j)*dy(i,j);
            Fw = rhow*u(i,j)*dy(i,j);
            Fn = rhon*v(i,j+1)*dx(i,j);
            Fs = rhos*v(i,j)*dx(i,j);
            % ae = -0.5*Fe + De;
            % aw =  0.5*Fw + Dw;
            % an = -0.5*Fn + Dn;
            % as =  0.5*Fs + Ds;
            % ap = ae + aw + an + as + (Fe - Fw) + (Fn - Fs);
            % me=0.5*(abs(Fe)-Fe);
            % mw=0.5*(abs(Fw)+Fe);
            % mn=0.5*(abs(Fn)-Fn);
            % ms=0.5*(abs(Fs)+Fs);
            % cds_upw=0.5*(Fe*(T_old(i,j)+T_old(i+1,j))-Fw*(T_old(i,j)+T_old(i-1,j))+Fn*(T_old(i,j)+T_old(i,j+1))-Fs*(T_old(i,j)+T_old(i,j-1)))...
            %     -((me+mw+mn+ms)*T_old(i,j)-(me*T_old(i+1,j)+mw*T_old(i-1,j)+mn*T_old(i,j+1)+ms*T_old(i,j-1)));
            me=max(-Fe,0);
            mw=max( Fw,0);
            mn=max(-Fn,0);
            ms=max( Fs,0);
        else
            me=0;
            mw=0;
            mn=0;
            ms=0;
        end

        
        % De0 = (kappa*dy(i,j)*hx(i,j));
        % Dw0 = (kappa*dy(i,j)*hx(i,j));
        % Dn0 = (kappa*dx(i,j)*hy(i,j));
        % Ds0 = (kappa*dx(i,j)*hy(i,j));
        % Pe=rhoe*u(i+1,j)*dx(i,j)/kappa;
        % Pw=rhow*u(i,j)*dx(i,j)/kappa;
        % Pn=rhon*v(i,j+1)*dy(i,j)/kappa;
        % Ps=rhos*v(i,j)*dy(i,j)/kappa;
        % De=De0*max(1-0.1*abs(Pe)^5,0);
        % Dw=Dw0*max(1-0.1*abs(Pw)^5,0);
        % Dn=Dn0*max(1-0.1*abs(Pn)^5,0);
        % Ds=Ds0*max(1-0.1*abs(Ps)^5,0);
        % De=De0*max(1-0.5*Pe,0);
        % Dw=Dw0*max(1-0.5*Pw,0);
        % Dn=Dn0*max(1-0.5*Pn,0);
        % Ds=Ds0*max(1-0.5*Ps,0);
        De = (kappa(i,j)*dy(i,j)*hx(i,j));
        Dw = (kappa(i,j)*dy(i,j)*hx(i,j));
        Dn = (kappa(i,j)*dx(i,j)*hy(i,j));
        Ds = (kappa(i,j)*dx(i,j)*hy(i,j));
        
        ae = me + De;
        aw = mw + Dw;
        an = mn + Dn;
        as = ms + Ds;

        if j==Ny-1 %北边界
            an = 0;
            % Sp=-dx(i,j)/(1/data.alpha+dy(i,j)/data.lammda);
            Spmax=dx(i,j)/(1/data.alpha+dy(i,j)/data.lammda);
            Sc=Spmax*data.Tf;
        else
            Sc=0;
            Spmax=0;
        end

        if j==2   %南边界
            as = 0;%because south area=0
            aw = 0;%because zeroGradient
            ae = 0;%because zeroGradient
        end

        if i==2     %西边界
            aw = 0;%because zeroGradient
        end

        if i==Nx-1  %东边界
            ae = 0;%because zeroGradient
        end


        ap = ae + aw + an + as + Spmax ;%+ ap0;% + Fe - Fw + Fn - Fs;
        ap=ap/WSOR.alpha_u;

        Dp = 0;

        AP(i,j) =  ap;
        AE(i,j) = -ae;
        AW(i,j) = -aw;
        AN(i,j) = -an;
        AS(i,j) = -as;
        % SP(i,j) =  Dp;%-WSOR.beta*cds_upw;
        SP(i,j) =  Dp+Sc;%+(Sp+Spmax)*T_old(i,j);%-WSOR.beta*cds_upw+ap0*u_old(i,j);%
        % apu(i,j) = ap;%ap;

        % A(ip,ip)    = AP(i,j);
        % A(ip,ip+1)  = AE(i,j);
        % A(ip,ip-1)  = AW(i,j);
        % A(ip,ip+Nx) = AN(i,j);
        % A(ip,ip-Nx) = AS(i,j);
        % B(ip,1)     = SP(i,j);
    end
end
% % update coefficients along boundaries
% AE(Nx-1, :)  = 0.0d0;
% AW(2,  :)  = 0.0d0;
% AN(:, Ny-1)  = 0.0d0;
% AS(:, 2 )  = 0.0d0;
% AP=-(AE + AW + AN + AS);
% %setting reference cell pressure to initialized pressure
% AP(2, 2) = 1.0d40 ;
% SP(2, 2) =0;




GError = 1.0;
iter = 1;
while ( GError > MaxG )
    for i = 2:Nx-1
        for j = 2:Ny-1
            Be = AE(i,j)*T_new(i+1,j);
            Bw = AW(i,j)*T_new(i-1,j);
            Bn = AN(i,j)*T_new(i,j+1);
            Bs = AS(i,j)*T_new(i,j-1);
            T_new(i,j) = (1-WSOR.alpha_u)*T_old(i,j)+...
                (SP(i,j)-Be-Bw-Bn-Bs)/AP(i,j);
            % u_new(i,j) = WSOR.alpha_u*u_new(i,j)+(1.0-WSOR.alpha_u)*u_old(i,j);
        end
    end
    
    GError = (sum(sum((T_new-T_old).^2.0)))/(Nx*Ny);
    
    T_old = T_new;
    
    iter = iter+1;
    
    if( iter > MaxI )
        break;
    end
    if( iter < MinI )
        GError = 1.0;
    end
    
end


end

% function b=sig(a)
%     if a<=0
%         b=0;
%     else
%         b=a;
%     end
% end