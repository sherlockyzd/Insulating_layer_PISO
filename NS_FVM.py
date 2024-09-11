# ***************************************************************************************************
# *   Simulate lid-driven cavity flow , using PISO algorithm, by presented code.
# *   I take no responsibilities for any errors in the code or damage thereby.
# *   Please notify me at sherlockyzd1994@gmail.com if the code is used in any type of application.
# ***************************************************************************************************
# *   Developer   :  Yu Zhengdong(04-28-2024)
# ***************************************************************************************************
# *   References  : 
# *   An Introduction to Computational Fluid Dynamics: The Finite Volume Method 2nd Edition
# *   by H. Versteeg (Author), W. Malalasekera (Author)
# ***************************************************************************************************
# *   Navier–Stokes equations  
# ***************************************************************************************************
import numpy as np
import math
import matplotlib.pyplot as plt
import MeshFVM as msh
import solver as slv

# 物理模型参数
Nr = 50      # r方向划分的数目
Ntetta = 50  # theta方向划分的数目
R1 = 0.45    # 流体的半径
R = 0.5     # 保温层的外半径
rho0 = 800.0 # 初始密度
rholay = 50  # 隔热层密度
mu0 = 10 * np.exp(3.332 - 0.0382 * 50)  # 采用拟合公式计算粘度
lammda = 0.1667
cp = 2000    # 内部流体导热系数
kappa0 = lammda / cp   # 等于导热系数 k / cp 比热容
lammdalay = 0.003       # 隔热层导热系数
cplay = 3000           # 隔热层比热容
kappalay = lammdalay / cplay   # 隔热层导热系数
cita=0.01
# 传热系数，热对流系数，辐射系数，流体温度，环境温度
T_init = 50.0    # 初始化为 50 度

# 数值模拟参数设置
IW = 20      # write frequency 每5步屏幕输出一次
CFL = 1.0    # 克朗数，目前代码用不上

Lx = np.pi
Ly = R

# 网格坐标设置
M = Ntetta
N = int(Nr * R1 / R)
MT = Ntetta
NT = Nr

NxU = M + 2
NyU = N + 2
NxV = M + 2
NyV = N + 2
NxP = M + 2
NyP = N + 2
NxT = MT + 2
NyT = NT + 2

rp = Ly / NT * np.array([0.5 + i for i in range(NyT)])
theta = Lx / MT * np.array([-0.5 + i for i in range(NxT)])
dx = np.tile(rp * Lx / MT, (NxT, 1))
dy = Ly / NT * np.ones((NxT, NyT))
hx = 1.0 / dx
hy = 1.0 / dy
dx2 = dx ** 2
dy2 = dy ** 2
rpmatrix = np.tile(rp, (NxT, 1))
thetamatrix = np.tile(theta, (NyT, 1)).T


# 创建径向网格
m = msh.MeshGenerator().create_radial_mesh_2d(Nr, Ntetta, R, Lx)
# 生成网格
X, Y = np.meshgrid(m.cellcenters['x'], m.cellcenters['y'], indexing='ij')  # 'ij' indexing for ndgrid equivalent in MATLAB
Xf, Yf = np.meshgrid(m.facecenters['x'], m.facecenters['y'], indexing='ij')

# 曲线坐标转换为笛卡尔坐标
xp = X * np.sin(Y)
yp = X * np.cos(Y)
xfp = Xf * np.sin(Yf)
yfp = Xf * np.cos(Yf)


# 绘制网格
plt.figure(1,figsize=(8, 8))
# 绘制网格线
plt.plot(xfp[:N+1,:], yfp[:N+1,:], 'r-', lw=0.5)  # 红色的网格线
plt.plot(xfp[:N+1,:].T, yfp[:N+1,:].T, 'r-', lw=0.5)  # 转置来画另一组网格线
plt.plot(xfp[N:,:], yfp[N:,:], 'k-', lw=0.5)  # 红色的网格线
plt.plot(xfp[N:,:].T, yfp[N:,:].T, 'k-', lw=0.5)  # 转置来画另一组网格线
# 绘制面中心的空心圆圈
plt.plot(xp, yp, 'bo', markerfacecolor='none', markersize=1)  # 蓝色的空心圆圈
# 设置绘图比例相同
plt.gca().set_aspect('equal', adjustable='box')
plt.title('Radial Mesh in Cartesian Coordinates')
plt.xlabel('X')
plt.ylabel('Y')
# 显示图形
# plt.show()
plt.pause(0.2)
plt.savefig('RadialMesh.png')



# 初始化变量
u = np.zeros((NxU, NyU))
v = np.zeros((NxV, NyV))
p = np.zeros((NxP, NyP))
T0 = T_init * np.ones((NxT, NyT))

rho = rho0 * np.ones((NxT, NyT))
rho[:, N + 2:NT + 2] = rholay
mu = mu0 * np.ones((NxU, NyU))
kappa = kappa0 * np.ones((NxT, NyT))
kappa[:, N + 2:NT + 2] = kappalay


# data = {
#     'lammda': lammdalay,'alpha': 25,'cita': cita,'Tf': 20,'Tenv': 5,'rp': rpmatrix,
#     'theta': thetamatrix,'dx': dx,'dy': dy,'hx': hx,'hy': hy,'dx2': dx2,'dy2': dy2
# }

Tf=20
Tenv=5
# 初始化边界条件
#北边界
T0[:, N + 1:NT] = Tf
T0[:, NT + 1] = Tenv
u[:, N + 1] = 0.0
v[:, N + 1] = 0.0
p[:, N + 1] = p[:, N]
#南边界
T0[:, 0] = T0[:,1]
u[:, 0] = u[:,1]
v[:, 0] = v[:,1]
p[:, 0] = p[:,1]
#东边界
T0[MT + 1, :] = T0[MT, :]
u[M + 1, :] = 0.0
v[M + 1, :] = v[M, :]
p[M + 1, :] = p[M, :]
#西边界
T0[0, :] = T0[1, :]
u[0, :] = 0.0
v[0, :] = v[1, :]
p[0, :] = p[1, :]

T = T0
Time = 0.0
rhoold = rho
k = 0
j=0
zz = 1
shuju = np.zeros((100, 2))  # 创建用于存储的数组
shuju[zz - 1, 0] = Time
shuju[zz - 1, 1] = np.mean(T[1:-1, 28])

MaxTime = 3600
dt = 60.0
# data['dt'] = dt
SaveErrors = []  # 创建一个二维数组
##
alpha_u=0.5
alpha_v= 0.5
alpha_p=0.001
alpha_T=0.9
beta_T=0.5
alpha=25
beta= 0.0
MaxI = 800  # max iteration
MaxG = 1e-6  # 收敛误差
MaxU=500
MinU=100
PType=2
MaxPG=1e-5
MaxP=1000
MinP=500
MaxTG=1e-5
MaxT=10000
MinT=1000

# WSOR = {
#     'alpha_u': 0.5, 'alpha_v': 0.5, 'alpha_p': 0.001,
#     'beta': 0.0, 'MaxG': MaxG, 'MaxU': 100, 'MinU': 20,
#     'PType': 2, 'MaxPG': 1e-7, 'MaxP': 1000, 'MinP': 500,
#     'MaxTG': 1e-7, 'MaxT': 3000, 'MinT': 800
# }
# 迭代求解
while Time < MaxTime:
    i = 0
    global_ErrorU = 1.0
    # ustar, apu, AEu, AWu, ANu, ASu = slv.U_Solver0(u, v, p, NxU, NyU, data, WSOR, rhoold, rho, mu, T)
    # vstar, apv, AEv, AWv, ANv, ASv = slv.V_Solver0(u, v, p, NxV, NyV, data, WSOR, rhoold, rho, mu, T)
    ustar, apu, AEu, AWu, ANu, ASu = slv.U_Solver(u, v, p, NxU, NyU, dx, dy, hx, hy, thetamatrix, rpmatrix, cita, Tf, dt, MaxG, MaxU, MinU, alpha_u, beta, rhoold, rho, mu, T)
                                     
    vstar, apv, AEv, AWv, ANv, ASv = slv.V_Solver(u, v, p, NxV, NyV, dx, dy, hx, hy, thetamatrix, rpmatrix, cita, Tf, dt, MaxG, MaxU, MinU, alpha_v, beta, rhoold, rho, mu, T)
    while global_ErrorU > MaxG:
        # pprime = slv.P_Solver(ustar, apu, vstar, apv, NxP, NyP, data, WSOR, rhoold, rho)
        pprime=slv.P_Solver(ustar, apu, vstar, apv, NxP, NyP, dx, dy, dx2, dy2, dt, MaxPG, MaxP, MinP, rhoold, rho)
        ustarstar, vstarstar, pstar, src = slv.Correction(ustar, apu, vstar, apv, p, pprime, NxU, NyU, NxV, NyV, NxP, NyP,rhoold, rho,alpha_p,dx,dy,dt)
        ustarstarprime, vstarstarprime = slv.Correction2(ustar, ustarstar, apu, AEu, AWu, ANu, ASu, vstar, vstarstar, apv, AEv, AWv, ANv, ASv, NxU, NyU, NxV, NyV)
        pprimeprime = slv.P_Solver(ustarstarprime, apu, vstarstarprime, apv, NxP, NyP, dx, dy, dx2, dy2, dt, MaxPG, MaxP, MinP, rhoold, rho)
        ustarstarstar, vstarstarstar, pstarstar, src = slv.Correction(ustarstarprime, apu, vstarstarprime, apv, pstar, pprimeprime, NxU, NyU, NxV, NyV, NxP, NyP,rhoold, rho,alpha_p,dx,dy,dt)
        
        UError = np.sqrt(np.sum((ustarstarstar - ustarstar) ** 2) / (NxU * NyU))
        VError = np.sqrt(np.sum((vstarstarstar - vstarstar) ** 2) / (NxV * NyV))
        PError = np.sqrt(np.sum((pstarstar - pstar) ** 2) / (NxP * NyP))
        
        global_ErrorU = UError + VError + PError #+src
        # 在每次迭代时追加数据
        SaveErrors.append([k, UError, VError, PError, src])
        if i % IW == 0:
            print(f"Time: {Time:.4f}, Iter: {i}, Xvelocity Residual: {UError:.12f},Yvelocity Residual: {VError:.12f}, Pressure Residual: {PError:.12f}, Global Error: {global_ErrorU:.12f}")
            T, GError, iter = slv.T_Solver(ustarstarstar, vstarstarstar, T0, NxT, NyT, rho, kappa, dx, dy, hx, hy, dt,MaxTG, MaxT, MinT,beta_T,alpha_T,alpha,lammdalay,Tenv)
            TError =  np.sqrt(np.sum((T - T0) ** 2) / (NxP * NyP))
            T0 = T
            print(f"Time: {Time:.4f}, Temperature Residual: {TError:.12f}")
        
        if i > MaxI:
            break
        
        i += 1
        j += 1
        k += 1
        ustar = ustarstarstar
        vstar = vstarstarstar
        p = pstarstar

    u = ustarstarstar
    v = vstarstarstar
    p = pstarstar
    T, GError, iter = slv.T_Solver(u, v, T0, NxT, NyT, rho, kappa, dx, dy, hx, hy, dt,MaxTG, MaxT, MinT,beta_T,alpha_T,alpha,lammdalay,Tenv)
    TError =  np.sqrt(np.sum((T - T0) ** 2) / (NxP * NyP))
    
    print(f"Time: {Time:.4f}, Temperature Residual: {TError:.12f}")
    T0 = T
    k += 1
    if k % 10 == 0:
        zz += 1
        shuju[zz - 1, 0] = Time
        shuju[zz - 1, 1] = np.mean(T[1:-1, 28])
    
    Time += dt
    # 温度分布
    plt.figure(8)
    # plt.pcolormesh(xfp, yfp, T_interp.T, shading='auto', cmap='inferno')
    # beautify_plot('Temperature', 'X', 'Y', 'T')
    # plt.pause(0.5)
    T_interp = T[0:-1, 0:-1]
    Tmax=math.ceil(max(T_interp.T.flatten()))
    Tmin=math.floor(min(T_interp.T.flatten()))
    levels = np.linspace(Tmin, Tmax, 50)  # 自定义等高线值
    # levels_log = np.logspace(np.log10(Tenv+1e-3), np.log10(Tmax), 20)  # 10 个对数分布的级别
    levels_low = np.linspace(Tenv, Tenv + (Tmax - Tenv) * 0.9, 5)  # 低温区少等高线
    levels_high = np.geomspace(Tenv + (Tmax - Tenv) * 0.92, Tmax, 10)  # 高温区多等高线
    levels_log = np.concatenate([levels_low, levels_high])
    # plt.contourf(X, Y, T_new, levels=levels, cmap='coolwarm')
    # 绘制填充等高线图
    # 绘制黑色的等高线
    cs = plt.contour(xfp, yfp, T_interp.T,levels=levels_log, colors='black', linewidths=0.8)  # 绘制黑色的等高线
    # plt.clabel(cs, inline=True, fontsize=8)  # 添加等高线标签
    plt.contourf(xfp, yfp, T_interp.T, levels=levels, cmap='coolwarm')  # 选择颜色映射和等高线级别
    plt.colorbar(label='Temperature (T)')  # 添加颜色条显示温度范围
    # beautify_plot('Temperature Contour Plot', 'X', 'Y', 'T','coolwarm')
    plt.title('Temperature Contour Plot')
    plt.xlabel('X', fontsize=14, fontweight='bold')
    plt.ylabel('Y', fontsize=14, fontweight='bold')
    plt.grid(False)
    plt.axis('equal')
    plt.tight_layout()
    # plt.pause(1.5)
    # plt.show()
    plt.savefig('Temperature.png')
    plt.pause(2.5)
    plt.close(8)

# Plotting functions
# 根据需要使用 matplotlib 进行绘图
# 绘图部分

SaveErrors = np.array(SaveErrors)
# Figure 1: 收敛性图（质量残差）
plt.figure(2)
plt.grid(True)
plt.semilogy(SaveErrors[:, 0], SaveErrors[:, 4], 'b', linewidth=1)
plt.title('Convergence of Mass', fontsize=14, fontweight='bold')
plt.xlabel('Time', fontsize=14, fontweight='bold')
plt.ylabel('Log(Mass)', fontsize=14, fontweight='bold')
plt.pause(0.5)
plt.savefig('MassResidual.png')

# Figure 2: 收敛性图（X方向速度残差）
plt.figure(3)
plt.grid(True)
plt.semilogy(SaveErrors[1:, 0], SaveErrors[1:, 1], 'r', linewidth=1)
plt.title('Convergence of X velocity', fontsize=14, fontweight='bold')
plt.xlabel('Time', fontsize=14, fontweight='bold')
plt.ylabel('Log(Norm2(Error))', fontsize=14, fontweight='bold')
plt.pause(0.5)
plt.savefig('XvelocityResidual.png')

# Figure 3: 收敛性图（Y方向速度残差）
plt.figure(4)
plt.grid(True)
plt.semilogy(SaveErrors[1:, 0], SaveErrors[1:, 2], 'g', linewidth=1)
plt.title('Convergence of Y velocity', fontsize=14, fontweight='bold')
plt.xlabel('Time', fontsize=14, fontweight='bold')
plt.ylabel('Log(Norm2(Error))', fontsize=14, fontweight='bold')
plt.pause(0.5)
# 格式化输出为短格式
# np.set_printoptions(precision=4)
plt.pause(0.5)
plt.savefig('YvelocityResidual.png')

# 关闭所有图形窗口
plt.close('all')


# usita_interp表示插值，把面心的u值插值到体心上去
usita_interp = 0.5 * (u[0:-1, 0:-1] + u[1:, 0:-1])
# vr_interp表示插值，把面心的v值插值到体心上去
vr_interp = 0.5 * (v[0:-1, 0:-1] + v[0:-1, 1:])
# 把theta和r方向的速度转化到x方向速度
u_interp = usita_interp * np.cos(Yf[0:N+1, 0:M+1].T) + vr_interp * np.sin(Yf[0:N+1, 0:M+1].T)
# 把theta和r方向的速度转化到y方向速度
v_interp = -usita_interp * np.sin(Yf[0:N+1, 0:M+1].T) + vr_interp * np.cos(Yf[0:N+1, 0:M+1].T)
# 把虚拟边界上的值去除
p_interp = p[0:-1, 0:-1]
T_interp = T[0:-1, 0:-1]


# 美化作图函数
def beautify_plot(title, xlabel, ylabel, cbar_label, cmap='viridis'):
    plt.title(title, fontsize=16, fontweight='bold')
    plt.xlabel(xlabel, fontsize=14, fontweight='bold')
    plt.ylabel(ylabel, fontsize=14, fontweight='bold')
    plt.colorbar(label=cbar_label)
    plt.grid(False)
    plt.axis('equal')
    plt.tight_layout()

# X方向速度分布
plt.figure(5)
# 绘制网格线
plt.plot(xfp[:N+1, :M+1], yfp[:N+1, :M+1], 'r-', lw=0.1)  # 红色的网格线
plt.plot(xfp.T[:M+1, :N+1], yfp.T[:M+1, :N+1], 'r-', lw=0.1)  # 转置来画另一组网格线
plt.pcolormesh(xfp[:N+1, :M+1], yfp[:N+1, :M+1], (u_interp.T), shading='gouraud', cmap='plasma')
beautify_plot('X Velocity', 'X', 'Y', 'u')
plt.pause(0.5)
plt.savefig('XVelocity.png')


# Y方向速度分布
plt.figure(6)
plt.plot(xfp[:N+1, :M+1], yfp[:N+1, :M+1], 'r-', lw=0.1)  # 红色的网格线
plt.plot(xfp.T[:M+1, :N+1], yfp.T[:M+1, :N+1], 'r-', lw=0.1)  # 转置来画另一组网格线
plt.pcolormesh(xfp[:N+1, :M+1], yfp[:N+1, :M+1], v_interp.T, shading='auto', cmap='plasma')
beautify_plot('Y Velocity', 'X', 'Y', 'v')
plt.pause(0.5)
plt.savefig('YVelocity.png')

# 压力分布
plt.figure(7)
plt.plot(xfp[:N+1, :M+1], yfp[:N+1, :M+1], 'r-', lw=0.1)  # 红色的网格线
plt.plot(xfp.T[:M+1, :N+1], yfp.T[:M+1, :N+1], 'r-', lw=0.1)  # 转置来画另一组网格线
plt.pcolormesh(xfp[:N+1, :M+1], yfp[:N+1, :M+1], p_interp.T, shading='auto', cmap='coolwarm')
beautify_plot('Pressure', 'X', 'Y', 'p')
plt.pause(0.5)
plt.savefig('Pressure.png')

# # 温度分布
# plt.figure(8)
# # plt.pcolormesh(xfp, yfp, T_interp.T, shading='auto', cmap='inferno')
# # beautify_plot('Temperature', 'X', 'Y', 'T')
# # plt.pause(0.5)
# import math
# Tmax=math.ceil(max(T_interp.T.flatten()))
# # Tmin=math.floor(min(T_interp.T.flatten()))
# levels = np.arange(Tenv, Tmax, 0.1)  # 自定义等高线值
# # plt.contourf(X, Y, T_new, levels=levels, cmap='coolwarm')
# # 绘制填充等高线图
# # 绘制黑色的等高线
# # cs = plt.contour(xfp, yfp, T_interp.T, colors='black')  # 绘制黑色的等高线
# # plt.clabel(cs, inline=True, fontsize=8)  # 添加等高线标签
# plt.contourf(xfp, yfp, T_interp.T, levels=levels, cmap='coolwarm')  # 选择颜色映射和等高线级别
# plt.colorbar(label='Temperature (T)')  # 添加颜色条显示温度范围
# # beautify_plot('Temperature Contour Plot', 'X', 'Y', 'T','coolwarm')
# plt.title('Temperature Contour Plot')
# plt.xlabel('X', fontsize=14, fontweight='bold')
# plt.ylabel('Y', fontsize=14, fontweight='bold')
# plt.grid(False)
# plt.axis('equal')
# plt.tight_layout()


# plt.savefig('Temperature.png')

# 创建图像并设置图像大小和分辨率
plt.figure(9,figsize=(10, 8), dpi=100)
# 下采样数据，步长为2以减少矢量密度
step = 1
xp_sample = xfp[:N+1:step, :M+1:step]
yp_sample = yfp[:N+1:step, :M+1:step]
u_sample = u_interp[::step, ::step].T
v_sample = v_interp[::step, ::step].T
# 计算速度大小用于颜色映射，并手动设置颜色映射范围
magnitude = np.sqrt(u_sample**2 + v_sample**2)  # 速度大小
vmin = np.min(magnitude)
vmax = np.max(magnitude)
# 矢量场的参数调整
plt.quiver(xp_sample, yp_sample, u_sample, v_sample, magnitude, 
           cmap='plasma', scale=20, width=0.003, headwidth=4, headlength=6, alpha=1.0, clim=(vmin, vmax))
# 添加标题和标签，并设置字体大小和加粗效果
plt.title('Velocity Field', fontsize=18, fontweight='bold')
plt.xlabel('X', fontsize=14, fontweight='bold')
plt.ylabel('Y', fontsize=14, fontweight='bold')
# 隐藏网格并保持比例一致
plt.grid(False)
plt.axis('equal')
# 添加颜色条，并设置其标签的字体大小
cbar = plt.colorbar(label='Velocity Magnitude', shrink=0.8)
# cbar.set_clim(vmin, vmax)  # 设置颜色条的范围
cbar.set_label('Velocity Magnitude', fontsize=12)
# 自动调整布局，避免标签和图像重叠
plt.tight_layout()
plt.pause(1.5)
plt.savefig('VelocityQuiver.png')

plt.close('all')



pass
