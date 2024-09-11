function [ pprime ] = P_Slover( us,apu,vs,apv,Nx,Ny,data,WSOR,rhoold,rho)

SType = WSOR.PType;

dx = data.dx;
dy = data.dy;
% dt = data.dt;
% hx = data.hx;
% hy = data.hy;
dx2 = data.dx2;
dy2 = data.dy2;
% rho = data.rho;

MaxG = WSOR.MaxPG;
MaxI = WSOR.MaxP;
MinI = WSOR.MinP;

% AP = zeros(Nx,Ny);
AE = zeros(Nx,Ny);
AW = zeros(Nx,Ny);
AN = zeros(Nx,Ny);
AS = zeros(Nx,Ny);
SP = zeros(Nx,Ny);
pprime  = zeros(Nx,Ny);
% p_new = p;

pprime_old = pprime;



for i = 2:Nx-1
    for j = 2:Ny-1
        
        rhoe=0.5*(rho(i,j)+rho(i+1,j));
        rhow=0.5*(rho(i,j)+rho(i-1,j));
        rhon=0.5*(rho(i,j)+rho(i,j+1));
        rhos=0.5*(rho(i,j)+rho(i,j-1));
        
        um = dy(i,j)*( rhoe*us(i+1,j) - rhow*us(i,j) );
        vm = dx(i,j)*( rhon*vs(i,j+1) - rhos*vs(i,j) );
        
        ae = rhoe*dy2(i,j)/apu(i+1,j);
        aw = rhow*dy2(i,j)/apu(i,j);
        an = rhon*dx2(i,j)/apv(i,j+1);
        as = rhos*dx2(i,j)/apv(i,j);

        AE(i,j) = -ae;
        AW(i,j) = -aw;
        AN(i,j) = -an;
        AS(i,j) = -as;
        SP(i,j) = -um-vm;%+(rhoold(i,j)-rho(i,j))*dx(i,j)*dy(i,j)/dt;%
       
    end
end

% % update coefficients along boundaries
% AE(Nx-1, :)  = 0.0d0;
% AW(2,  :)  = 0.0d0;
% AN(:, Ny-1)  = 0.0d0;
% AS(:, 2 )  = 0.0d0;
% update coefficients along boundaries
AE(Nx-1, :)  = 0.0d0;
AW(2,  :)  = 0.0d0;
AN(:, Ny-1)  = 0.0d0;
AS(:, 2 )  = 0.0d0;
AP=-(AE + AW + AN + AS);
%setting reference cell pressure to initialized pressure
AP(2, 2)=1.0d40 ;
AE(2, 2)=0.0d40 ;
AW(2, 2)=0.0d40 ;
AN(2, 2)=0.0d40 ;
AS(2, 2)=0.0d40 ;
SP(2, 2) =0;

% AP(2,2) = 1.0;
% A(Nx+2,Nx+2) = AP(2,2);
% B(Nx+2,1) = 0.0;
% 
if ( SType == 1 )
    % A = zeros(Nx*Ny,Nx*Ny);
    X_sparse=zeros(Nx*Ny,1);
    Y_sparse=zeros(Nx*Ny,1);
    A_sparse=zeros(Nx*Ny,1);
    B = zeros(Nx*Ny,1);
    Isum=0;
    for i = 2:Nx-1
        for j = 2:Ny-1
            ip = (j-1)*Nx+i;
            % A(ip,ip)    = AP(i,j);
            % A(ip,ip+1)  = AE(i,j);
            % A(ip,ip-1)  = AW(i,j);
            % A(ip,ip+Nx) = AN(i,j);
            % A(ip,ip-Nx) = AS(i,j);
            Isum=Isum+1;
            X_sparse(Isum)=ip;
            Y_sparse(Isum)=ip;
            A_sparse(Isum)=AP(i,j);

            Isum=Isum+1;
            X_sparse(Isum)=ip;
            Y_sparse(Isum)=ip+1;
            A_sparse(Isum)=AE(i,j);

            Isum=Isum+1;
            X_sparse(Isum)=ip;
            Y_sparse(Isum)=ip-1;
            A_sparse(Isum)=AW(i,j);

            Isum=Isum+1;
            X_sparse(Isum)=ip;
            Y_sparse(Isum)=ip+Nx;
            A_sparse(Isum)=AN(i,j);

            Isum=Isum+1;
            X_sparse(Isum)=ip;
            Y_sparse(Isum)=ip-Nx;
            A_sparse(Isum)=AS(i,j);


            B(ip,1)     = SP(i,j);
        end
    end
    i = 1;
    for j = 1:Ny
        ip = (j-1)*Nx+i;
        % A(ip,ip) = 1.0;
        % A(ip,ip+1) = -1.0;
        Isum=Isum+1;
        X_sparse(Isum)=ip;
        Y_sparse(Isum)=ip;
        A_sparse(Isum)=1.0;

        Isum=Isum+1;
        X_sparse(Isum)=ip;
        Y_sparse(Isum)=ip+1;
        A_sparse(Isum)=-1.0;


        B(ip,1) = 0.0;
    end
    i = Nx;
    for j = 1:Ny
        ip = (j-1)*Nx+i;
        % A(ip,ip) = 1.0;
        % A(ip,ip-1) = -1.0;
        Isum=Isum+1;
        X_sparse(Isum)=ip;
        Y_sparse(Isum)=ip;
        A_sparse(Isum)=1.0;

        Isum=Isum+1;
        X_sparse(Isum)=ip;
        Y_sparse(Isum)=ip-1;
        A_sparse(Isum)=-1.0;

        B(ip,1) = 0.0;
    end

    j = 1;
    for i = 2:Nx-1
        ip = (j-1)*Nx+i;
        % A(ip,ip) = 1.0;
        % A(ip,ip+Nx) = -1.0;
        Isum=Isum+1;
        X_sparse(Isum)=ip;
        Y_sparse(Isum)=ip;
        A_sparse(Isum)=1.0;

        Isum=Isum+1;
        X_sparse(Isum)=ip;
        Y_sparse(Isum)=ip+Nx;
        A_sparse(Isum)=-1.0;

        B(ip,1) = 0.0;
    end
    j = Ny;
    for i = 2:Nx-1
        ip = (j-1)*Nx+i;
        % A(ip,ip) = 1.0;
        % A(ip,ip-Nx) = -1.0;
        Isum=Isum+1;
        X_sparse(Isum)=ip;
        Y_sparse(Isum)=ip;
        A_sparse(Isum)=1.0;

        Isum=Isum+1;
        X_sparse(Isum)=ip;
        Y_sparse(Isum)=ip-Nx;
        A_sparse(Isum)=-1.0;

        B(ip,1) = 0.0;
    end

    % Xp = A\B;
    % tol=1e-8;
    % maxit=100;

    A=sparse(X_sparse(1:Isum),Y_sparse(1:Isum),A_sparse(1:Isum),Nx*Ny,Nx*Ny);%第一二列代表位置，第三列代表元素
    % L = ichol(A);
    % [Xp,~,relres,iter]=bicg(A,B,MaxG,MaxI);
    % setup = struct('type','ilutp','droptol',1e-6);
    % [L,U] = ilu(A,setup);
    % [Xp,~,relres,iter]=bicg(A,B,MaxG,MaxI,L,U);
    [Xp]=pcg(A,B,MaxG,MaxI);
    % Xp=lsqr(A,B,MaxG,MaxI);

    for i = 1:Nx
        for j = 1:Ny
            ip = (j-1)*Nx+i;
            pprime(i,j) = Xp(ip);
        end
    end

else
    % % update coefficients along boundaries
    % AE(Nx-1, :)  = 0.0d0;
    % AW(2,  :)  = 0.0d0;
    % AN(:, Ny-1)  = 0.0d0;
    % AS(:, 2 )  = 0.0d0;
    % AP=-(AE + AW + AN + AS);
    % %setting reference cell pressure to initialized pressure
    % AP(2, 2)=1.0d40 ;
    % SP(2, 2) =0;
        
    GError = 1.0;
    iter = 1;
    while ( GError > MaxG )
        for i = 2:Nx-1
            for j = 2:Ny-1
                Be = AE(i,j)*pprime(i+1,j);
                Bw = AW(i,j)*pprime(i-1,j);
                Bn = AN(i,j)*pprime(i,j+1);
                Bs = AS(i,j)*pprime(i,j-1);
                pprime(i,j) = (SP(i,j)-Be-Bw-Bn-Bs)/AP(i,j);
            end
        end
        
        GError = (sum(sum((pprime-pprime_old).^2.0)))/(Nx*Ny);
        pprime_old = pprime;
        
        iter = iter+1;
        
        if( iter > MaxI )
            break;
        end
        if( iter < MinI )
            GError = 1.0;
        end
        
    end
    
end

end

