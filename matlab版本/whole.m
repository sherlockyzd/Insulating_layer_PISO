%***************************************************************************************************
%*   Simulate lid-driven cavity flow , using SIMPLE algorithm, by presented code.
%*   I take no responsibilities for any errors in the code or damage thereby.
%*   Please notify me at sherlockyzd1994@gmail.com if the code is used in any type of application.
%***************************************************************************************************
%*   Developer   :  Yu Zhengdong(04-28-2024)
%***************************************************************************************************
%*   References  :
%*   An Introduction to Computational Fluid Dynamics: The Finite Volume Method 2nd Edition
%*   by H. Versteeg (Author), W. Malalasekera (Author)
%***************************************************************************************************
%*   Navier–Stokes equations
%***************************************************************************************************
clear,clc
close all
format compact
% format long
%%
%物理模型参数
Nr = 50;      %r方向划分的数目
Ntetta = 50;  %theta方向划分的数目
R1=0.17;     %流体的半径
R = 0.21;   %保温层的外半径
rho0 = 870.0;   %初始密度
rholay = 50; %隔热层密度
% Re=10000;         %雷诺数
% mu0 = rho0*Vi*Lx/Re;
% mu0 =rho0*sqrt(9.8*R)*Lx/Re;    %粘度的公式，可以自己给
mu0 = 10*exp(3.332-0.0382*50);    %采用拟合公式计算粘度
% data.U_lid = Vi;
% data.kappa=0.001;   % 等于 导热系数k/cp比热容
lammda = 0.1667;
cp = 2000;  %内部流体导热系数
kappa0 = lammda/cp;   % 等于导热系数k/cp比热容
lammdalay =0.03;%隔热层导热系数
cplay = 1500; %隔热层比热容
kappalay = lammdalay/cplay;   % 隔热层导热系数k/cp比热容
data.lammda=lammdalay; % 北边界为隔热层
data.alpha = 25;
cita=0.67;      %温度系数
data.cita = cita;
data.Tf = 20;      %动量方程体积力用到的参数,-g*rou*beta*(T-Tf)*sintheta
T_init=50.d0;   %初始化为50度
%%
%数值模拟参数设置
MaxI =1200;  % max iteration
MaxG= 1e-6;   %收敛误差
IW = 10;       % write frequency每5步屏幕输出一次
CFL=1.0;      %克朗数，目前代码用不上
WSOR.alpha_u = 0.5;     % u方向的松弛系数
WSOR.alpha_v = 0.5;     % v方向的松弛系数
WSOR.alpha_p = 0.01;   % p方向的松弛系数
WSOR.beta=0.0;% convection mixture coffecient对流混合系数，中心差分和迎风格式的混合系数
WSOR.MaxG = 1e-7;    % U方程迭代的收敛残差，在Usolver中17行
WSOR.MaxU = 100;     % U方程迭代的最大步数
WSOR.MinU = 20;      % U方程迭代的最小步数
WSOR.PType=2;        % Type=2表示高斯赛德尔求解
WSOR.MaxPG = 1e-7;   % P方程迭代的收敛残差，在Psolver中
WSOR.MaxP = 1000;    % P方程迭代的最大步数
WSOR.MinP = 800;     % P方程迭代的最小步数
WSOR.MaxTG = 1e-7;
WSOR.MaxT = 3000;
WSOR.MinT = 800;
%%
Lx = pi*1;  %seta的范围
Ly = R;     %r的范围
%网格坐标设置
M = Ntetta;     %中心控制体的个数
N = floor(Nr*R1/R);
MT = Ntetta;     %中心控制体的个数
NT = Nr;
NxU = M+2;    %加了两个虚拟边界网格
NyU = N+2;    %后文用到创建U的矩阵
NxV = M+2;
NyV = N+2;
NxP = M+2;
NyP = N+2;
NxT = MT+2;
NyT = NT+2;
rp=Ly/NT*[0,0.5:1:NyT-1];     %所有的rp的集合
theta=Lx/MT*[-0.5,0.5:1:NxT-1]';   %theta的范围
dx = repmat(rp*Lx/MT,NxT,1);%+zeros(NxP,NyP);   %就是rp*delta_theta
dy = Ly/NT+zeros(NxT,NyT);  %就是delta_r
hx = 1.0./dx;    %为了计算方便
hy = 1.0./dy;    %为了计算方便，算好放那儿
dx2 = dx.*dx;    %为了计算方便，算好放那儿
dy2 = dy.*dy;    %为了计算方便，算好放那儿
rpmatrix=repmat(rp,NxT,1);
thetamatrix=repmat(theta,1,NyT);
data.rp = rpmatrix;    % data.rp是一个结构体，为了方便减少形参，如Usolver中3~6行
data.theta =thetamatrix; %repmat命令是复制平铺矩阵
data.dx = dx;
data.dy = dy;
data.hx = hx;
data.hy = hy;
data.dx2 = dx2;
data.dy2 = dy2;
m = createMeshRadial2D(Nr, Ntetta, R, Lx);
[X, Y] = ndgrid(m.cellcenters.x, m.cellcenters.y);
[Xf,Yf]=ndgrid(m.facecenters.x, m.facecenters.y);
xp = X.*sin(Y);   %曲线坐标转换为笛卡尔坐标
yp = X.*cos(Y);   %曲线坐标转换为笛卡尔坐标
xfp = Xf.*sin(Yf);   %曲线坐标转换为笛卡尔坐标
yfp = Xf.*cos(Yf);   %曲线坐标转换为笛卡尔坐标
% dt=1e-5;
% data.dt = dt;
%initialize primary variables初始化最初变量
u = zeros(NxU,NyU);
v = zeros(NxV,NyV);
p = zeros(NxP,NyP);
T0 = T_init+zeros(NxT,NyT);    %初始化为50度
% data.rho = rho0-cita*(T-20);
% data.mu = mu0*data.rho./rho0;
% rho = rho0-cita*(T-20);
% mu = mu0*rho./rho0;
rho = rho0+zeros(NxT,NyT);  %初始化整个流体的密度
rho(:,N+2:NT+2) = rholay;
mu = mu0+zeros(NxU,NyU);    %初始化整个流体的粘度
kappa = kappa0+zeros(NxT,NyT);    %初始化整个流体的粘度
kappa(:,N+2:NT+2) = kappalay;

% set bc - east boundary
u(M+2, :) = 0.0d0;    %d0是双精度，小数点后面8个0，M+2是最右边的行
% v(M+2, :) = 0.0d0;      zeroGradient
% T0(M+2, :) = 20.0d0;    zeroGradient

% set bc - west boundary
u(2, :)    = 0.0d0;    %因交错网格的影响，左边界为第2行
% v(1, :)    = 0.0d0;    zeroGradient
% T0(1, :)    = 20.0d0;  zeroGradient

% set bc - north boundary
u(:, N+2) = 0.0d0;
v(:, N+2) = 0.0d0;
T0(:,N+2:NT+2) = data.Tf;

% set bc - south boundary
% u(:, 1)    = 0.0d0;    AS=0;
v(:, 2)    = 0.0d0;
% T0(:, 1)    = 20.0d0;  AS=0;

T=T0;

%%
Time=0.0;
% Told = T;
rhoold=rho;
i = 1;
k=1;

zz=1;
shuju(zz,1)=Time;
shuju(zz,2)=mean(T(2:end-1,28));
MaxTime=10;
dt=1e-2;
data.dt = dt;
%%
while (  Time < MaxTime  )
    j=1;
    global_ErrorU=1.0;
    [ ustar , apu , AEu , AWu , ANu , ASu] = U_Slover( u,v,p,NxU,NyU,data,WSOR,rhoold,rho,mu,T);
    [ vstar , apv , AEv , AWv , ANv , ASv] = V_Slover( u,v,p,NxV,NyV,data,WSOR,rhoold,rho,mu,T);
    while (  global_ErrorU > MaxG  )  %总误差
        [pprime ] = P_Slover( ustar,apu,vstar,apv,NxP,NyP,data,WSOR,rhoold,rho); %求解压力泊松方程，即连续性方程%pprime就是p一撇
        [ustarstar,vstarstar,pstar,src]   = Correction(  ustar,apu,vstar,apv,p,pprime,NxU,NyU,NxV,NyV,NxP,NyP,WSOR,data,rhoold,rho); %修正u和v和p
        [ustarstarprime,vstarstarprime] = Correction2( ustar,ustarstar,apu, AEu , AWu , ANu , ASu,...
            vstar,vstarstar, apv , AEv , AWv , ANv , ASv,NxU,NyU,NxV,NyV);
        [pprimeprime ] = P_Slover( ustarstarprime,apu,vstarstarprime,apv,NxP,NyP,data,WSOR,rhoold,rho);
        [ustarstarstar,vstarstarstar,pstarstar,src] = Correction( ustarstarprime,apu,vstarstarprime,apv,pstar,pprimeprime,NxU,NyU,NxV,NyV,NxP,NyP,WSOR,data,rhoold,rho);
        UError = (sum(sum((ustarstarstar-ustar).^2.0)))/(NxU*NyU);   %二范数，U方向的误差
        VError = (sum(sum((vstarstarstar-vstar).^2.0)))/(NxV*NyV);   %定义的v方向的误差
        PError = (sum(sum((pstarstar-p).^2.0)))/(NxP*NyP);   %定义的p方向的误差
        SaveErrors(i,1) = i;   % i是迭代次数，存下这些数便于后续画图
        SaveErrors(i,2) = UError;
        SaveErrors(i,3) = VError;
        SaveErrors(i,4) = PError;
        SaveErrors(i,5) = src;
        global_ErrorU=UError+VError+PError;%;+src;
        % if global_ErrorU<1  %当动量方程的误差足够小时，小于0.01，才计算温度方程
        % WSOR.PType=1;
        % end
        i = i+1;
        j = j+1;
        ustar = ustarstarstar;
        vstar = vstarstarstar;
        p     = pstarstar;
        % if global_ErrorU<0.1
        %     WSOR.PType=1;
        % else
        %     WSOR.PType=2;
        % end
        if ( rem(j,IW) == 0.0 )   % rem是整除取余数的语句，即i除以IW余数为0，能整除5
            % Writer( u,v,p,xp,yp,data,Lx,Ly,XL,YL,NL,M,N,i );
            % fprintf('Time: %4d , Dt:%15.12f\n',Time,dt);
            fprintf('Time: %4d,Iter: %4d , Xvelocity Residual: %15.12f \n',Time,j,UError);
            fprintf('Time: %4d,Iter: %4d , Yvelocity Residual: %15.12f \n',Time,j,VError);
            fprintf('Time: %4d,Iter: %4d , Pressure  Residual: %15.12f \n',Time,j,PError);
            % fprintf('Iter: %4d , netMass   Residual: %15.12f \n',i,src);
            fprintf('Time: %4d,Iter: %4d , globalSum Residual: %15.12f \n',Time,j,global_ErrorU);
            fprintf('------------------------------------------\n');
        end
        
        if( j > MaxI )%达到最大迭代步数，就跳出
            break;
        end
        
    end
    
    u = ustarstarstar;
    v = vstarstarstar;
    p     = pstarstar;
    [ T,GError,iter ] = T_Slover( u,v,T0,NxT,NyT,data,rho,kappa,WSOR);
    TError = (sum(sum((T-T0).^2.0)))/(NxP*NyP);
    fprintf('Time: %4d , Temprature  Residual: %15.12f \n',Time,TError);
    % SaveErrors(i,5) = TError;
    T0=T;     %给下一步迭代时候的值赋为T0=50
    k=k+1;
    if ( rem(k,10) == 0.0 )
        zz=zz+1;
        shuju(zz,1)=Time;
        shuju(zz,2)=mean(T(2:end-1,28));
    end
    U_magmax=sqrt(max(max((u.^2+v.^2))));
    Time=Time+dt;
end
%%
% [ s ] = Writer( unew,vnew,p,xp,yp,data,Lx,Ly,XL,YL,NL,M,N,i );
% xpp=linspace(0,Lx,NxP);
% xp(1) = 0.0; xp(1:NxP-2) = 0.5*dx:dx:Lx; xp(NxP) = Lx;
% yp(1) = 0.0; yp(2:NyP-1) = 0.5*dy:dy:Ly; yp(NyP) = Ly;
close all
usita_interp=0.5*(u(2:end-1,2:end-1)+u(3:end,2:end-1)); %interp表示插值，把面心的u值插值到体心上去
vr_interp=0.5*(v(2:end-1,2:end-1)+v(2:end-1,3:end));  %把面心的v值插值到体心上去
u_interp=usita_interp.*cos((Y(1:N,1:M))')+vr_interp.*sin((Y(1:N,1:M))');   %把theta和r方向的速度转化到x方向速度
v_interp=-usita_interp.*sin((Y(1:N,1:M))')+vr_interp.*cos((Y(1:N,1:M))');  %把theta和r方向的速度转化到y方向速度
p_interp=p(2:end-1,2:end-1);  %把虚拟边界上的值去除
T_interp=T(2:end-1,2:end-1);  %把虚拟边界上的值去除
format short

figure(1);hold on;grid on;
semilogy(SaveErrors(2:end,1),SaveErrors(2:end,5),'b','linewidth',2);
title('Convergency of Mass','FontSize',14,'FontWeight','bold');
xlabel('Time','FontSize',14,'FontWeight','bold');
ylabel('Log(Mass)','FontSize',14,'FontWeight','bold');

figure(2);hold on;grid on;
semilogy(SaveErrors(2:end,1),SaveErrors(2:end,2),'r','linewidth',2);
title('Convergency of Xvelocity','FontSize',14,'FontWeight','bold');
xlabel('Time','FontSize',14,'FontWeight','bold');
ylabel('Log(Norm2(Error))','FontSize',14,'FontWeight','bold');

figure(3);hold on;grid on;
semilogy(SaveErrors(2:end,1),SaveErrors(2:end,3),'g','linewidth',2);
title('Convergency of Yvelocity','FontSize',14,'FontWeight','bold');
xlabel('Time','FontSize',14,'FontWeight','bold');
ylabel('Log(Norm2(Error))','FontSize',14,'FontWeight','bold');


figure(4);
pcolor(xp(1:N,1:M), yp(1:N,1:M),u_interp');   %把u_interp转置
grid off;
set(gca,'FontSize',10,'FontName','Times New Roman');
axis equal;   %X轴和Y轴的比例尺相等
% xlim([0,1]);
% ylim([0,1]);
title('Xvelocity','FontSize',14,'FontWeight','bold');
xlabel('X','FontSize',14,'FontWeight','bold');
ylabel('Y','FontSize',14,'FontWeight','bold');
shading interp;box off;%去除上边框和右边框
set(gca,'DataAspectRatio',[1 1 1],'linewidth',1.2)
set(gca,'TickDir','out','TickLength',[0.005,0.005]);
set(gca,'XLim',[0,R],'XTick',(0:0.1:1)*R);
set(gca,'xTickLabel',num2str(get(gca,'xTick')','%.2f'));
set(gca,'YLim',[-R,R],'YTick',(-1:0.1:1)*R);
set(gca,'yTickLabel',num2str(get(gca,'yTick')','%.2f'));

% caxis([0,1]);
colormap(jet);
ch=colorbar;title(ch,'\fontsize{15}u');
% set(ch,'YTick',0:0.1:1);
set(ch,'TickDir','out','TickLength',0.01);
% set(ch,'yTickLabel',num2str(get(ch,'yTick')','%.1f'));

% axpos=get(gca,'Position');
% ch.Position(3)=1*ch.Position(3);
% set(gca,'Position',axpos);
% title('\fontsize{10}(a)','position',[-1.2,2]);
% print(gcf, '-dtiff','-r1500','veloc_x.tif');



figure(5);
% contour(v_interp',20);
pcolor(xp(1:N,1:M), yp(1:N,1:M),v_interp');  %把v_interp转置了
grid off;
set(gca,'FontSize',10,'FontName','Times New Roman');
axis equal;
% xlim([0,1]);
% ylim([0,1]);
title('Yvelocity','FontSize',14,'FontWeight','bold');
xlabel('X','FontSize',14,'FontWeight','bold');
ylabel('Y','FontSize',14,'FontWeight','bold');
shading interp;box off;
set(gca,'DataAspectRatio',[1 1 1],'linewidth',1.2)
set(gca,'TickDir','out','TickLength',[0.005,0.005]);
set(gca,'XLim',[0,R],'XTick',(0:0.1:1)*R);
set(gca,'xTickLabel',num2str(get(gca,'xTick')','%.2f'));
set(gca,'YLim',[-R,R],'YTick',(-1:0.1:1)*R);
set(gca,'yTickLabel',num2str(get(gca,'yTick')','%.2f'));

% caxis([-0.5,0.5]);
colormap(jet);
ch=colorbar;title(ch,'\fontsize{15}v');
% set(ch,'YTick',0:0.1:1);
set(ch,'TickDir','out','TickLength',0.01);
% set(ch,'yTickLabel',num2str(get(ch,'yTick')','%.1f'));

% axpos=get(gca,'Position');
% ch.Position(3)=1*ch.Position(3);
% set(gca,'Position',axpos);
% title('\fontsize{10}(a)','position',[-1.2,2]);
% print(gcf, '-dtiff','-r1500','veloc_y.tif');



figure(6);
% contourf(p',25);
pcolor(xp(1:N,1:M), yp(1:N,1:M),p_interp');
grid off;
set(gca,'FontSize',10,'FontName','Times New Roman');
axis equal;
% xlim([0,1]);
% ylim([0,1]);
title('Pressure','FontSize',14,'FontWeight','bold');
xlabel('X','FontSize',14,'FontWeight','bold');
ylabel('Y','FontSize',14,'FontWeight','bold');
shading interp;box off;
set(gca,'DataAspectRatio',[1 1 1],'linewidth',1.2)
set(gca,'TickDir','out','TickLength',[0.005,0.005]);
set(gca,'XLim',[0,R],'XTick',(0:0.1:1)*R);
set(gca,'xTickLabel',num2str(get(gca,'xTick')','%.2f'));
set(gca,'YLim',[-R,R],'YTick',(-1:0.1:1)*R);
set(gca,'yTickLabel',num2str(get(gca,'yTick')','%.2f'));

% clim([-1,1]);
% caxis([-1,1]);
colormap(jet);
ch=colorbar;title(ch,'\fontsize{15}p');
% set(ch,'YTick',0:0.1:1);
set(ch,'TickDir','out','TickLength',0.01);
% set(ch,'yTickLabel',num2str(get(ch,'yTick')','%.1f'));

% axpos=get(gca,'Position');
% ch.Position(3)=1*ch.Position(3);
% set(gca,'Position',axpos);
% title('\fontsize{10}(a)','position',[-1.2,2]);
% print(gcf, '-dtiff','-r1500','pressure.tif');

format short


figure(7);hold on;grid on;
% contourf(T',25);
pcolor(xp,yp,T_interp');
grid off;
set(gca,'FontSize',10,'FontName','Times New Roman');
axis equal;
% xlim([0,1]);
% ylim([0,1]);
% title('窑炉温度(对称)','FontSize',14,'FontWeight','bold');
xlabel('X','FontSize',14,'FontWeight','bold');
ylabel('Y','FontSize',14,'FontWeight','bold');
shading interp;box off;
set(gca,'DataAspectRatio',[1 1 1],'linewidth',1.2)
set(gca,'TickDir','out','TickLength',[0.005,0.005]);
set(gca,'XLim',[0,R],'XTick',(0:0.1:1)*R);
set(gca,'xTickLabel',num2str(get(gca,'xTick')','%.2f'));
set(gca,'YLim',[-R,R],'YTick',(-1:0.1:1)*R);
set(gca,'yTickLabel',num2str(get(gca,'yTick')','%.2f'));

axis off
% caxis([20,47]);
colormap(jet);
ch=colorbar;title(ch,'\fontsize{15}T');


% axpos=get(gca,'Position');
% ch.Position(3)=1*ch.Position(3);
% set(gca,'Position',axpos);
% title('\fontsize{10}(a)','position',[-1.2,2]);
% print(gcf, '-dtiff','-r1500','Tempreture.tif');


figure(8);
% plot(xp(1:N,1:M), yp(1:N,1:M), 'og','Markersize',1),hold on;
% plot(xp(N+1:NT,1:MT), yp(N+1:NT,1:MT), 'or','Markersize',1),hold on;
plot(xfp(1:NyP-1,1:NxP-1), yfp(1:NyP-1,1:NxP-1), '-k', (xfp(1:NyP-1,1:NxP-1))', (yfp(1:NyP-1,1:NxP-1))', '-k'),hold on;
plot(xfp(NyP-1:NyT-1,1:NxT-1), yfp(NyP-1:NyT-1,1:NxT-1), '-m', (xfp(NyP-1:NyT-1,1:NxT-1))', (yfp(NyP-1:NyT-1,1:NxT-1))', '-m')
axis equal
axis off

figure(9);
quiver(xp(1:N,1:M), yp(1:N,1:M),u_interp',v_interp'),hold on;    %quiver画出向量箭头
% sx=xp(:,1);sy=yp(:,1);
% hlines = streamline(xp,yp,u_interp',v_interp');
% set(hlines,'LineWidth',2,'Color','r')
grid off;
set(gca,'FontSize',10,'FontName','Times New Roman');
axis equal;
% xlim([0,1]);
% ylim([0,1]);
title('Velocity','FontSize',14,'FontWeight','bold');
xlabel('X','FontSize',14,'FontWeight','bold');
ylabel('Y','FontSize',14,'FontWeight','bold');
shading interp;box off;
set(gca,'DataAspectRatio',[1 1 1],'linewidth',1.2)
set(gca,'TickDir','out','TickLength',[0.005,0.005]);
set(gca,'XLim',[0,R],'XTick',(0:0.1:1)*R);
set(gca,'xTickLabel',num2str(get(gca,'xTick')','%.2f'));
set(gca,'YLim',[-R,R],'YTick',(-1:0.1:1)*R);
set(gca,'yTickLabel',num2str(get(gca,'yTick')','%.2f'));

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

function MS = createMesh2D(varargin)
% MeshStructure = createMesh2D(Nx, Ny, Lx, Ly)
% MeshStructure = createMesh2D(facelocationX, facelocationY)
% creates a uniform 2D mesh:
% Nx is the number of cells in x (horizontal) direction
% Ny is the number of cells in y (vertical) direction
% Lx is the domain length in x direction
% Ly is the domain length in y direction
%
% SYNOPSIS:
%   MeshStructure = createMesh2D(Nx, Ny, Lx, Ly)
%
% PARAMETERS:
%   Nx: number of cells in the x direction
%   Lx: domain length in x direction
%   Ny: number of cells in the y direction
%   Ly: domain length in y direction
%
% RETURNS:
%   MeshStructure.
%                 dimensions=2 (2D problem)
%                 numbering: shows the indexes of cellsn from left to right
%                 and top to bottom
%                 cellsize: x and y elements of the cell size =[Lx/Nx,
%                 Ly/Ny]
%                 cellcenters.x: location of each cell in the x direction
%                 cellcenters.y: location of each cell in the y direction
%                 facecenters.x: location of interface between cells in the
%                 x direction
%                 facecenters.y: location of interface between cells in the
%                 y direction
%                 numberofcells: [Nx, Ny]
%
%
% EXAMPLE:
%   Nx = 5;
%   Ny = 7;
%   Lx = 10;
%   Ly = 20;
%   m = createMesh2D(Nx, Ny, Lx, Ly);
%   [X, Y] = ndgrid(m.cellcenters.x, m.cellcenters.y);
%   [Xf,Yf]=ndgrid(m.facecenters.x, m.facecenters.y);
%   plot(X, Y, 'or', ...
%        Xf, Yf, '-b', Xf', Yf', '-b');
%
% SEE ALSO:
%     createMesh1D, createMesh3D, createMeshCylindrical1D, ...
%     createMeshCylindrical2D, createCellVariable, createFaceVariable

% Written by Ali A. Eftekhari
% See the license file

if nargin==4
    % uniform 1D mesh
    Nx=varargin{1};
    Ny=varargin{2};
    Width=varargin{3};
    Height=varargin{4};
    % cell size is dx
    dx = Width/Nx;
    dy = Height/Ny;
    G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
    CellSize.x= dx*ones(Nx+2,1);
    CellSize.y= dy*ones(Ny+2,1);
    CellSize.z= [0.0];
    CellLocation.x= [1:Nx]'*dx-dx/2;
    CellLocation.y= [1:Ny]'*dy-dy/2;
    CellLocation.z= [0.0];
    FaceLocation.x= [0:Nx]'*dx;
    FaceLocation.y= [0:Ny]'*dy;
    FaceLocation.z= [0.0];
elseif nargin==2
    % nonuniform 1D mesh
    facelocationX=varargin{1};
    facelocationY=varargin{2};
    facelocationX=facelocationX(:);
    facelocationY=facelocationY(:);
    Nx = length(facelocationX)-1;
    Ny = length(facelocationY)-1;
    G=reshape(1:(Nx+2)*(Ny+2), Nx+2, Ny+2);
    CellSize.x= [facelocationX(2)-facelocationX(1); ...
        facelocationX(2:end)-facelocationX(1:end-1); ...
        facelocationX(end)-facelocationX(end-1)];
    CellSize.y= [facelocationY(2)-facelocationY(1); ...
        facelocationY(2:end)-facelocationY(1:end-1); ...
        facelocationY(end)-facelocationY(end-1)];
    CellSize.z= [0.0];
    CellLocation.x= 0.5*(facelocationX(2:end)+facelocationX(1:end-1));
    CellLocation.y= 0.5*(facelocationY(2:end)+facelocationY(1:end-1));
    CellLocation.z= [0.0];
    FaceLocation.x= facelocationX;
    FaceLocation.y= facelocationY;
    FaceLocation.z= [0.0];
end
c=G([1,end], [1,end]);
MS=MeshStructure(2, ...
    [Nx,Ny], ...
    CellSize, ...
    CellLocation, ...
    FaceLocation, ...
    c(:), ...
    [1]);
end

function MS = createMeshRadial2D(varargin)
% MeshStructure = createMeshRadial2D(Nr, Ntetta, Lr, Tetta)
% MeshStructure = createMeshRadial2D(facelocationR, facelocationTetta)
% creates a uniform 2D mesh on a Radia coordinate:
% imagine a top view slice of pie
% Nr is the number of cells in r (radial) direction
% Ntetta is the number of cells in tetta direction
% Lr is the domain length in r direction
% Tetta is the domain length in tetta direction
%
% SYNOPSIS:
%   MeshStructure = createMeshRadial2D(Nr, Ntetta, Lr, Tetta)
%
% PARAMETERS:
%   Nr: number of cells in the x direction
%   Lr: domain length in x direction
%   Ntetta: number of cells in the y direction
%   Tetta: domain length in y direction 0<Tetta<=2*pi
%
% RETURNS:
%   MeshStructure.
%                 dimensions=2.8 (2D problem, radial coordinate)
%                 numbering: shows the indexes of cellsn from left to right
%                 and top to bottom
%                 cellsize: r and y elements of the cell size =[Lr/Nr,
%                 Ly/Ntetta]
%                 cellcenters.x: location of each cell in the r direction
%                 cellcenters.y: location of each cell in the y direction
%                 facecenters.x: location of interface between cells in the
%                 r direction
%                 facecenters.y: location of interface between cells in the
%                 y direction
%                 numberofcells: [Nr, Ntetta]
%
%
% EXAMPLE:
%   Nr = 5;
%   Ntetta = 7;
%   R = 10;
%   Ly = 2*pi;
%   m = createMeshRadial2D(Nx, Ntetta, Lx, Ly);
%   [X, Y] = ndgrid(m.cellcenters.x, m.cellcenters.y);
%   [Xf,Yf]=ndgrid(m.facecenters.x, m.facecenters.y);
%   plot(X, Y, 'or', ...
%        Xf, Yf, '-b', Xf', Yf', '-b');
%
% SEE ALSO:
%     createMesh1D, createMesh3D, createMeshCylindrical1D, ...
%     createMeshCylindrical2D, createCellVariable, createFaceVariable

% Written by Ali A. Eftekhari
% See the license file

if nargin==4
    % uniform 1D mesh
    Nx=varargin{1};
    Ny=varargin{2};
    Width=varargin{3};
    Tetta=varargin{4};
    if Tetta>2*pi
        warning('Tetta is higher than 2*pi. It is scaled to 2*pi');
        Tetta = 2*pi;
    end
    MS=createMesh2D(Nx, Ny, Width, Tetta);
elseif nargin==2
    % nonuniform 1D mesh
    facelocationX=varargin{1};
    facelocationY=varargin{2};
    if facelocationY(end)>2*pi
        facelocationY = facelocationY/facelocationY(end)*2.0*pi;
        warning('The domain size adjusted to match a maximum of 2*pi.')
    end
    MS=createMesh2D(facelocationX, facelocationY);
end
MS.dimension=2.8;
end
classdef MeshStructure
    %MeshStructure class
    % contains information about the domain and the mesh size
    
    properties
        dimension
        dims
        cellsize
        cellcenters
        facecenters
        corners
        edges
    end
    
    methods
        function meshVar = MeshStructure(dimension, dims, cellsize, ...
                cellcenters, facecenters, corners, edges)
            if nargin>0
                meshVar.dimension = dimension;
                meshVar.dims = dims;
                meshVar.cellsize = cellsize;
                meshVar.cellcenters = cellcenters;
                meshVar.facecenters = facecenters;
                meshVar.corners= corners;
                meshVar.edges= edges;
            end
        end
    end
end


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
% endfunction [ u_new , apu, AE , AW , AN , AS ] = U_Slover( u,v,p,Nx,Ny,data,WSOR,rhoold,rho,mu,T)

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


