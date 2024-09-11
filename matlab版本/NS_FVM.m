
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



% caxis([0,1]);
% colormap(jet);
% ch=colorbar;title(ch,'\fontsize{15}u');
% % set(ch,'YTick',0:0.1:1);
% set(ch,'TickDir','out','TickLength',0.01);
% set(ch,'yTickLabel',num2str(get(ch,'yTick')','%.1f'));
