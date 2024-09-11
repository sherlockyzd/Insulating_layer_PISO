import numpy as np
# from numba import jit
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import gmres,cg,spilu,LinearOperator,lsmr  # Conjugate Gradient solver
from numba import njit, prange


# @njit(parallel=True)
def U_Solver(u, v, p, Nx, Ny, dx, dy, hx, hy, theta, rp, cita, T0, dt, MaxG, MaxI, MinI, alpha_u, beta, rhoold, rho, mu, T):
    AP = np.zeros((Nx, Ny))
    AE = np.zeros((Nx, Ny))
    AW = np.zeros((Nx, Ny))
    AN = np.zeros((Nx, Ny))
    AS = np.zeros((Nx, Ny))
    SP = np.zeros((Nx, Ny))
    apu = np.ones((Nx, Ny))
    
    u_new = np.copy(u)
    u_old = np.copy(u)
    
    # 并行化循环
    for i in prange(1, Nx-1):  # 使用 prange 代替 range 以实现并行化
        for j in range(1, Ny-1):
            rhoe = rho[i, j]
            rhow = rho[i-1, j]
            rhon = 0.25 * (rho[i, j] + rho[i, j+1] + rho[i-1, j+1] + rho[i-1, j])
            rhos = 0.25 * (rho[i, j] + rho[i-1, j] + rho[i-1, j-1] + rho[i, j-1])
            mue = mu[i, j]
            muw = mu[i-1, j]
            mun = 0.25 * (mu[i, j] + mu[i, j+1] + mu[i-1, j+1] + mu[i-1, j])
            mus = 0.25 * (mu[i, j] + mu[i-1, j] + mu[i-1, j-1] + mu[i, j-1])
            
            De = mue * dy[i, j] * hx[i, j]
            Dw = muw * dy[i, j] * hx[i, j]
            Dn = mun * dx[i, j] * hy[i, j]
            Ds = mus * dx[i, j] * hy[i, j]
            
            Fe = 0.5 * rhoe * (u_old[i+1, j] + u_old[i, j]) * dy[i, j]
            Fw = 0.5 * rhow * (u_old[i-1, j] + u_old[i, j]) * dy[i, j]
            Fn = 0.5 * rhon * (v[i, j+1] + v[i-1, j+1]) * dx[i, j]
            Fs = 0.5 * rhos * (v[i, j] + v[i-1, j]) * dx[i, j]
            
            me = max(-Fe, 0)
            mw = max(Fw, 0)
            mn = max(-Fn, 0)
            ms = max(Fs, 0)
            
            cds_upw = (0.5 * (Fe * (u[i, j] + u[i+1, j]) - Fw * (u[i, j] + u[i-1, j]) + 
                              Fn * (u[i, j] + u[i, j+1]) - Fs * (u[i, j] + u[i, j-1])) - 
                       ((me + mw + mn + ms) * u[i, j] - 
                        (me * u[i+1, j] + mw * u[i-1, j] + mn * u[i, j+1] + ms * u[i, j-1])))
            Dp = (p[i-1, j] - p[i, j]) * dy[i, j]

            Sc = (mu[i, j] * dy[i, j] / rp[i, j] * (v[i, j+1] + v[i, j] - v[i-1, j] - v[i-1, j+1]) +
                  rhoold[i, j] * (1-cita * (T[i, j] - T0)) * 9.8 * np.sin(theta[i, j]) * dx[i, j] * dy[i, j])
            
            Sp = (-rhoold[i, j] * 0.25 * (v[i, j+1] + v[i, j] + v[i-1, j] + v[i-1, j+1]) * 
                  (theta[i+1, j] - theta[i, j]) * dy[i, j] -
                  mu[i, j] * dy[i, j] * (theta[i+1, j] - theta[i, j]) / rp[i, j])
            Spmax = max(-Sp, 0)
            
            ae = me + De
            aw = mw + Dw
            an = mn + Dn
            as_ = ms + Ds
            if j == 1:   
                as_ = 0  # South boundary condition
            
            ap0 = rhoold[i, j] * dy[i, j] * dx[i, j] / dt
            ap = ae + aw + an + as_ + Spmax + ap0
            ap /= alpha_u
            
            AP[i, j] = ap
            AE[i, j] = -ae
            AW[i, j] = -aw
            AN[i, j] = -an
            AS[i, j] = -as_
            SP[i, j] = Dp + Sc - beta * cds_upw + (Sp + Spmax) * u_old[i, j] + ap0 * u_old[i, j]
            apu[i, j] = ap
    
    GError = 1.0
    iter = 1
    while GError > MaxG:
        for i in prange(1, Nx-1):  # 同样并行化外层循环
            for j in range(1, Ny-1):
                Be = AE[i, j] * u_new[i+1, j]
                Bw = AW[i, j] * u_new[i-1, j]
                Bn = AN[i, j] * u_new[i, j+1]
                Bs = AS[i, j] * u_new[i, j-1]
                u_new[i, j] = (1 - alpha_u) * u_old[i, j] + (SP[i, j] - Be - Bw - Bn - Bs) / AP[i, j]
        
        GError = np.sqrt(np.sum((u_new - u_old)**2) / (Nx * Ny))

        u_new[:, Ny-1] = 0.0
        u_new[:, 0] = u_new[:,1]
        u_new[Nx-1, :] = 0.0
        u_new[0, :] = 0.0

        u_old[:] = u_new
        
        iter += 1
        if iter > MaxI:
            break
        if iter < MinI:
            GError = 1.0

    return u_new, apu, AE, AW, AN, AS


# @njit(parallel=True)
def V_Solver(u, v, p, Nx, Ny, dx, dy, hx, hy, theta, rp, cita, T0, dt, MaxG, MaxI, MinI, alpha_v, beta, rhoold, rho, mu, T):
    AP = np.zeros((Nx, Ny))
    AE = np.zeros((Nx, Ny))
    AW = np.zeros((Nx, Ny))
    AN = np.zeros((Nx, Ny))
    AS = np.zeros((Nx, Ny))
    SP = np.zeros((Nx, Ny))
    apv = np.ones((Nx, Ny))
    
    v_new = np.copy(v)
    v_old = np.copy(v)

    # 并行化主要循环
    for i in prange(1, Nx-1):  # 并行化外层循环
        for j in range(1, Ny-1):
            # 计算密度和粘度的均值
            rhoe = 0.25 * (rho[i, j] + rho[i+1, j] + rho[i, j-1] + rho[i+1, j-1])
            rhow = 0.25 * (rho[i, j] + rho[i-1, j] + rho[i-1, j-1] + rho[i, j-1])
            rhon = rho[i, j]
            rhos = rho[i, j-1]
            
            mue = 0.25 * (mu[i, j] + mu[i+1, j] + mu[i, j-1] + mu[i+1, j-1])
            muw = 0.25 * (mu[i, j] + mu[i-1, j] + mu[i-1, j-1] + mu[i, j-1])
            mun = mu[i, j]
            mus = mu[i, j-1]
            
            # 计算离散化系数
            De = mue * dy[i, j] * hx[i, j]
            Dw = muw * dy[i, j] * hx[i, j]
            Dn = mun * dx[i, j] * hy[i, j]
            Ds = mus * dx[i, j] * hy[i, j]
            
            Fe = 0.5 * rhoe * (u[i+1, j] + u[i+1, j-1]) * dy[i, j]
            Fw = 0.5 * rhow * (u[i, j] + u[i, j-1]) * dy[i, j]
            Fn = 0.5 * rhon * (v_old[i, j+1] + v_old[i, j]) * dx[i, j]
            Fs = 0.5 * rhos * (v_old[i, j-1] + v_old[i, j]) * dx[i, j]
            
            me = max(-Fe, 0)
            mw = max(Fw, 0)
            mn = max(-Fn, 0)
            ms = max(Fs, 0)
            
            cds_upw = (0.5 * (Fe * (v[i, j] + v[i+1, j]) - Fw * (v[i, j] + v[i-1, j]) +
                              Fn * (v[i, j] + v[i, j+1]) - Fs * (v[i, j] + v[i, j-1])) - 
                       ((me + mw + mn + ms) * v[i, j] - (me * v[i+1, j] + mw * v[i-1, j] +
                                                        mn * v[i, j+1] + ms * v[i, j-1])))
            
            Dp = (p[i, j-1] - p[i, j]) * dx[i, j]
            Sc = (mu[i, j] * dy[i, j] / rp[i, j] * (u[i, j-1] + u[i, j] - u[i+1, j] - u[i+1, j-1]) +
                  rhoold[i, j] * (0.25 * (u[i, j-1] + u[i, j] + u[i+1, j] + u[i+1, j-1]))**2 *
                  (theta[i+1, j] - theta[i, j]) * dy[i, j] -
                  rhoold[i, j] * (1-cita * (T[i, j] - T0)) * 9.8 * np.cos(theta[i, j]) * dx[i, j] * dy[i, j])
            Sp = -mu[i, j] * dy[i, j] * (theta[i+1, j] - theta[i, j]) / rp[i, j]
            Spmax = max(-Sp, 0)
            
            ae = me + De
            aw = mw + Dw
            an = mn + Dn
            as_ = ms + Ds
            
            # 边界条件
            if i == 1:   
                aw = 0  # 西侧零梯度
            if i == Nx-2:
                ae = 0  # 东侧零梯度
            
            ap0 = rhoold[i, j] * dy[i, j] * dx[i, j] / dt
            ap = ae + aw + an + as_ + Spmax + ap0
            ap /= alpha_v
            
            AP[i, j] = ap
            AE[i, j] = -ae
            AW[i, j] = -aw
            AN[i, j] = -an
            AS[i, j] = -as_
            SP[i, j] = Dp + Sc - beta * cds_upw + (Sp + Spmax) * v_old[i, j] + ap0 * v_old[i, j]
            apv[i, j] = ap
    
    GError = 1.0
    iter = 1
    while GError > MaxG:
        for i in prange(1, Nx-1):  # 并行化循环
            for j in range(1, Ny-1):
                Be = AE[i, j] * v_new[i+1, j]
                Bw = AW[i, j] * v_new[i-1, j]
                Bn = AN[i, j] * v_new[i, j+1]
                Bs = AS[i, j] * v_new[i, j-1]
                v_new[i, j] = (1 - alpha_v) * v_old[i, j] + (SP[i, j] - Be - Bw - Bn - Bs) / AP[i, j]
        
        GError = np.sqrt(np.sum((v_new - v_old)**2) / (Nx * Ny))

        # 更新边界条件
        v_new[:, Ny-1] = 0.0  # 北边界
        v_new[:, 0] = 0.0     # 南边界
        v_new[Nx-1, :] = v_new[Nx-2, :]  # 东边界
        v_new[0, :] = v_new[1, :]        # 西边界

        v_old[:] = v_new  # 避免频繁复制，直接更新

        iter += 1
        if iter > MaxI:
            break
        if iter < MinI:
            GError = 1.0

    return v_new, apv, AE, AW, AN, AS


# @njit(parallel=True)
def T_Solver(u, v, T0, Nx, Ny, rho, kappa, dx, dy, hx, hy, dt, MaxG, MaxI, MinI,beta,alpha_T,alpha,lammda,Tenv): 
    Nu = u.shape
    AP = np.zeros((Nx, Ny))
    AE = np.zeros((Nx, Ny))
    AW = np.zeros((Nx, Ny))
    AN = np.zeros((Nx, Ny))
    AS = np.zeros((Nx, Ny))
    SP = np.zeros((Nx, Ny))
    T_new = np.copy(T0)
    T_old = np.copy(T0)
    
    # 并行化主循环
    for i in prange(1, Nx-1):  # 使用 prange 并行化外层循环
        for j in range(1, Ny-1):
            if i < Nu[0]-1 and j < Nu[1]-1:
                rhoe = 0.5 * (rho[i, j] + rho[i+1, j])
                rhow = 0.5 * (rho[i, j] + rho[i-1, j])
                rhon = 0.5 * (rho[i, j] + rho[i, j+1])
                rhos = 0.5 * (rho[i, j] + rho[i, j-1])
                Fe = rhoe * u[i+1, j] * dy[i, j]
                Fw = rhow * u[i, j] * dy[i, j]
                Fn = rhon * v[i, j+1] * dx[i, j]
                Fs = rhos * v[i, j] * dx[i, j]
                
                me = max(-Fe, 0)
                mw = max(Fw, 0)
                mn = max(-Fn, 0)
                ms = max(Fs, 0)
                # cds_upw=0.5*(Fe*(T_old(i,j)+T_old(i+1,j))-Fw*(T_old(i,j)+T_old(i-1,j))+Fn*(T_old(i,j)+T_old(i,j+1))-Fs*(T_old(i,j)+T_old(i,j-1)))-((me+mw+mn+ms)*T_old(i,j)-(me*T_old(i+1,j)+mw*T_old(i-1,j)+mn*T_old(i,j+1)+ms*T_old(i,j-1)))
                # me, mw, mn, ms 是边界通量/质量流量，Fe, Fw, Fn, Fs 是系数
                # cds_upw = 0.5 * (Fe * (T_old[i, j] + T_old[i+1, j]) - Fw * (T_old[i, j] + T_old[i-1, j])
                #                 + Fn * (T_old[i, j] + T_old[i, j+1]) - Fs * (T_old[i, j] + T_old[i, j-1])) \
                #         - ((me + mw + mn + ms) * T_old[i, j] - (me * T_old[i+1, j] + mw * T_old[i-1, j]
                #                                                 + mn * T_old[i, j+1] + ms * T_old[i, j-1]))
            else:
                me, mw, mn, ms,cds_upw = 0, 0, 0, 0,0
            
            De = (kappa[i, j] * dy[i, j] * hx[i, j])
            Dw = (kappa[i, j] * dy[i, j] * hx[i, j])
            Dn = (kappa[i, j] * dx[i, j] * hy[i, j])
            Ds = (kappa[i, j] * dx[i, j] * hy[i, j])
            
            ae = me + De
            aw = mw + Dw
            an = mn + Dn
            as_ = ms + Ds
            
            if j == Ny-2:  # North boundary
                an = 0
                Spmax = dx[i, j] / (1 / alpha + dy[i, j] / lammda)
                Sc = Spmax * Tenv
            else:
                Sc = 0
                Spmax = 0
            
            if j == 1:  # South boundary
                as_ = 0
                aw = 0
                ae = 0
            
            if i == 1:  # West boundary
                aw = 0
            
            if i == Nx-2:  # East boundary
                ae = 0
            
            ap0 = rho[i, j] * dy[i, j] * dx[i, j] / dt
            ap = ae + aw + an + as_ + Spmax + ap0
            ap /= alpha_T
            
            Dp = 0
            
            AP[i, j] = ap
            AE[i, j] = -ae
            AW[i, j] = -aw
            AN[i, j] = -an
            AS[i, j] = -as_
            SP[i, j] = Dp + Sc + ap0 * T_old[i, j] #- beta*cds_upw
    
    GError = 1.0
    iter = 1
    APinv=1.0/AP
    while GError > MaxG:
        # for i in prange(1, Nx-1):  # 并行化外层循环
        #     for j in range(1, Ny-1):
        #         Be = AE[i, j] * T_new[i+1, j]
        #         Bw = AW[i, j] * T_new[i-1, j]
        #         Bn = AN[i, j] * T_new[i, j+1]
        #         Bs = AS[i, j] * T_new[i, j-1]
        #         T_new[i, j] = (1 - alpha_u) * T_old[i, j] + (SP[i, j] - Be - Bw - Bn - Bs) / AP[i, j]
        
        # GError = np.sum((T_new - T_old)**2) / (Nx * Ny)

        # 更新内部网格点的温度矩阵（不包括边界）
        Be = AE[1:-1, 1:-1] * T_new[2:, 1:-1]  # i+1
        Bw = AW[1:-1, 1:-1] * T_new[:-2, 1:-1] # i-1
        Bn = AN[1:-1, 1:-1] * T_new[1:-1, 2:]  # j+1
        Bs = AS[1:-1, 1:-1] * T_new[1:-1, :-2] # j-1
        # 矢量化的温度更新公式
        f1=(SP[1:-1, 1:-1] - Be - Bw - Bn - Bs)* APinv[1:-1, 1:-1]
        T_new[1:-1, 1:-1] = (1 - alpha_T) * T_old[1:-1, 1:-1] +  f1

        Be = AE[1:-1, 1:-1] * T_new[2:, 1:-1]  # i+1
        f2=(SP[1:-1, 1:-1] - Be - Bw - Bn - Bs)* APinv[1:-1, 1:-1]
        T_new[1:-1, 1:-1] = (1 - alpha_T) * T_new[1:-1, 1:-1] +  f2

        Bw = AW[1:-1, 1:-1] * T_new[:-2, 1:-1] # i-1
        f3=(SP[1:-1, 1:-1] - Be - Bw - Bn - Bs)* APinv[1:-1, 1:-1]
        T_new[1:-1, 1:-1] = (1 - alpha_T) * T_new[1:-1, 1:-1] +  f3

        Bn = AN[1:-1, 1:-1] * T_new[1:-1, 2:]  # j+1
        f4=(SP[1:-1, 1:-1] - Be - Bw - Bn - Bs)* APinv[1:-1, 1:-1]
        T_new[1:-1, 1:-1] = (1 - alpha_T) * T_new[1:-1, 1:-1] +  f4

        Bs = AS[1:-1, 1:-1] * T_new[1:-1, :-2] # j-1
        f5=(SP[1:-1, 1:-1] - Be - Bw - Bn - Bs)* APinv[1:-1, 1:-1]
        T_new[1:-1, 1:-1] = (1 - alpha_T) * T_new[1:-1, 1:-1] +  0.2*(f1+f2+f3+f4+f5)
        
        # 更新边界条件
        # T_new[:, Ny-1] = Tenv           # 北边界
        T_new[:, 0] = T_new[:, 1]       # 南边界
        T_new[Nx-1, :] = T_new[Nx-2, :] # 东边界
        T_new[0, :] = T_new[1, :]       # 西边界

        # 计算全局误差（矢量化计算）
        GError = np.sum((T_new - T_old) ** 2) / (Nx * Ny)
        # 避免频繁的内存复制，直接更新
        T_old[:] = T_new
        
        iter += 1
        if iter > MaxI:
            break
        if iter < MinI:
            GError = 1.0
    
    return T_new, GError, iter


# @njit(parallel=True)
def P_Solver(us, apu, vs, apv, Nx, Ny, dx, dy, dx2, dy2, dt, MaxG, MaxI, MinI, rhoold, rho):
    # 初始化系数矩阵
    AE = np.zeros((Nx, Ny))
    AW = np.zeros((Nx, Ny))
    AN = np.zeros((Nx, Ny))
    AS = np.zeros((Nx, Ny))
    SP = np.zeros((Nx, Ny))
    pprime = np.zeros((Nx, Ny))
    pprime_old = np.copy(pprime)

    # 边界条件初始化
    pprime[0, :] = pprime[1, :]
    pprime[-1, :] = pprime[-2, :]
    pprime[:, 0] = pprime[:, 1]
    pprime[:, -1] = pprime[:, -2]

    # 并行化主要循环
    for i in prange(1, Nx-1):  # 并行化外层循环
        for j in range(1, Ny-1):
            # 计算系数
            rhoe = 0.5 * (rho[i, j] + rho[i+1, j])
            rhow = 0.5 * (rho[i, j] + rho[i-1, j])
            rhon = 0.5 * (rho[i, j] + rho[i, j+1])
            rhos = 0.5 * (rho[i, j] + rho[i, j-1])

            um = dy[i, j] * (rhoe * us[i+1, j] - rhow * us[i, j])
            vm = dx[i, j] * (rhon * vs[i, j+1] - rhos * vs[i, j])

            ae = rhoe * dy2[i, j] / apu[i+1, j]
            aw = rhow * dy2[i, j] / apu[i, j]
            an = rhon * dx2[i, j] / apv[i, j+1]
            as_ = rhos * dx2[i, j] / apv[i, j]

            # 设置离散化矩阵
            AE[i, j] = -ae
            AW[i, j] = -aw
            AN[i, j] = -an
            AS[i, j] = -as_
            SP[i, j] = -um - vm + (rhoold[i, j] - rho[i, j]) * dx[i, j] * dy[i, j] / dt

    # 边界条件处理
    AE[Nx-2, :] = 0.0
    AW[1, :] = 0.0
    AN[:, Ny-2] = 0.0
    AS[:, 1] = 0.0

    AP = -(AE + AW + AN + AS)
    AP[1, 1] = 1.0e40  # 防止奇异性
    SP[1, 1] = 0.0

    # 初始化误差和迭代次数
    GError = 1.0
    iter_ = 1
    inv_AP = 1.0 / AP  # 预先计算 AP 的倒数用于加速

    while GError > MaxG:
        # 使用向量化和切片操作来替代循环
        Be = AE[1:-1, 1:-1] * pprime[2:, 1:-1]
        Bw = AW[1:-1, 1:-1] * pprime[:-2, 1:-1]
        Bn = AN[1:-1, 1:-1] * pprime[1:-1, 2:]
        Bs = AS[1:-1, 1:-1] * pprime[1:-1, :-2]
        pprime[1:-1, 1:-1] = (SP[1:-1, 1:-1] - Be - Bw - Bn - Bs) * inv_AP[1:-1, 1:-1]

        # 处理边界条件
        pprime[0, :] = pprime[1, :]
        pprime[-1, :] = pprime[-2, :]
        pprime[:, 0] = pprime[:, 1]
        pprime[:, -1] = pprime[:, -2]
        # 使用向量化和切片操作来替代循环
        Be = AE[1:-1, 1:-1] * pprime[2:, 1:-1]
        pprime[1:-1, 1:-1] = (SP[1:-1, 1:-1] - Be - Bw - Bn - Bs) * inv_AP[1:-1, 1:-1]
        Bw = AW[1:-1, 1:-1] * pprime[:-2, 1:-1]
        pprime[1:-1, 1:-1] = (SP[1:-1, 1:-1] - Be - Bw - Bn - Bs) * inv_AP[1:-1, 1:-1]
        Bn = AN[1:-1, 1:-1] * pprime[1:-1, 2:]
        pprime[1:-1, 1:-1] = (SP[1:-1, 1:-1] - Be - Bw - Bn - Bs) * inv_AP[1:-1, 1:-1]
        Bs = AS[1:-1, 1:-1] * pprime[1:-1, :-2]
        pprime[1:-1, 1:-1] = (SP[1:-1, 1:-1] - Be - Bw - Bn - Bs) * inv_AP[1:-1, 1:-1]

        # 处理边界条件
        pprime[0, :] = pprime[1, :]
        pprime[-1, :] = pprime[-2, :]
        pprime[:, 0] = pprime[:, 1]
        pprime[:, -1] = pprime[:, -2]

        # 计算全局误差
        GError = np.sum((pprime - pprime_old)**2) / (Nx * Ny)
        pprime_old[:] = pprime  # 避免复制，直接覆盖

        iter_ += 1
        if iter_ > MaxI:
            break
        if iter_ < MinI:
            GError = 1.0

    return pprime


# @jit(nopython=True)
def U_Solver0(u, v, p, Nx, Ny, data, WSOR, rhoold, rho, mu, T):
    dx = data['dx']
    dy = data['dy']
    hx = data['hx']
    hy = data['hy']
    theta = data['theta']
    rp = data['rp']
    cita = data['cita']
    T0 = data['Tf']
    dt = data['dt']
    
    MaxG = WSOR['MaxG']
    MaxI = WSOR['MaxU']
    MinI = WSOR['MinU']
    
    AP = np.zeros((Nx, Ny))
    AE = np.zeros((Nx, Ny))
    AW = np.zeros((Nx, Ny))
    AN = np.zeros((Nx, Ny))
    AS = np.zeros((Nx, Ny))
    SP = np.zeros((Nx, Ny))
    apu = np.ones((Nx, Ny))
    u_new = np.copy(u)
    u_old = np.copy(u)
    
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
            rhoe = rho[i, j]
            rhow = rho[i-1, j]
            rhon = 0.25 * (rho[i, j] + rho[i, j+1] + rho[i-1, j+1] + rho[i-1, j])
            rhos = 0.25 * (rho[i, j] + rho[i-1, j] + rho[i-1, j-1] + rho[i, j-1])
            mue = mu[i, j]
            muw = mu[i-1, j]
            mun = 0.25 * (mu[i, j] + mu[i, j+1] + mu[i-1, j+1] + mu[i-1, j])
            mus = 0.25 * (mu[i, j] + mu[i-1, j] + mu[i-1, j-1] + mu[i, j-1])
            
            De = (mue * dy[i, j] * hx[i, j])
            Dw = (muw * dy[i, j] * hx[i, j])
            Dn = (mun * dx[i, j] * hy[i, j])
            Ds = (mus * dx[i, j] * hy[i, j])
            
            Fe = 0.5 * rhoe * (u_old[i+1, j] + u_old[i, j]) * dy[i, j]
            Fw = 0.5 * rhow * (u_old[i-1, j] + u_old[i, j]) * dy[i, j]
            Fn = 0.5 * rhon * (v[i, j+1] + v[i-1, j+1]) * dx[i, j]
            Fs = 0.5 * rhos * (v[i, j] + v[i-1, j]) * dx[i, j]
            
            me = max(-Fe, 0)
            mw = max(Fw, 0)
            mn = max(-Fn, 0)
            ms = max(Fs, 0)
            cds_upw = (0.5 * (Fe * (u[i, j] + u[i+1, j]) - Fw * (u[i, j] + u[i-1, j]) + 
                              Fn * (u[i, j] + u[i, j+1]) - Fs * (u[i, j] + u[i, j-1])) -
                       ((me + mw + mn + ms) * u[i, j] - (me * u[i+1, j] + mw * u[i-1, j] + 
                                                        mn * u[i, j+1] + ms * u[i, j-1])))
            
            Dp = (p[i-1, j] - p[i, j]) * dy[i, j]
            Sc = (mu[i, j] * dy[i, j] / rp[i, j] * (v[i, j+1] + v[i, j] - v[i-1, j] - v[i-1, j+1]) +
                  rhoold[i, j] * (1-cita * (T[i, j] - T0)) * (9.8) * np.sin(theta[i, j]) * dx[i, j] * dy[i, j])
            Sp = (-rhoold[i, j] * 0.25 * (v[i, j+1] + v[i, j] + v[i-1, j] + v[i-1, j+1]) * 
                  (theta[i+1, j] - theta[i, j]) * dy[i, j] -
                  mu[i, j] * dy[i, j] * (theta[i+1, j] - theta[i, j]) / rp[i, j])
            Spmax = max(-Sp, 0)
            
            ae = me + De
            aw = mw + Dw
            an = mn + Dn
            as_ = ms + Ds
            if j == 1:   
                as_ = 0  # South boundary condition
            
            ap0 = rhoold[i, j] * dy[i, j] * dx[i, j] / dt
            ap = ae + aw + an + as_ + Spmax + ap0
            ap /= WSOR['alpha_u']
            
            AP[i, j] = ap
            AE[i, j] = -ae
            AW[i, j] = -aw
            AN[i, j] = -an
            AS[i, j] = -as_
            SP[i, j] = Dp + Sc - WSOR['beta'] * cds_upw + (Sp + Spmax) * u_old[i, j] + ap0 * u_old[i, j]
            apu[i, j] = ap
    
    GError = 1.0
    iter = 1
    while GError > MaxG:
        for i in range(1, Nx-1):
            for j in range(1, Ny-1):
                Be = AE[i, j] * u_new[i+1, j]
                Bw = AW[i, j] * u_new[i-1, j]
                Bn = AN[i, j] * u_new[i, j+1]
                Bs = AS[i, j] * u_new[i, j-1]
                u_new[i, j] = (1 - WSOR['alpha_u']) * u_old[i, j] + (SP[i, j] - Be - Bw - Bn - Bs) / AP[i, j]
        
        GError = np.sqrt(np.sum((u_new - u_old)**2) / (Nx * Ny))
        #北边界
        u_new[:, Ny-1] = 0.0
        #南边界
        u_new[:, 0] = u_new[:,1]
        #东边界
        u_new[Nx-1, :] = 0.0
        #西边界
        u_new[0, :] = 0.0

        u_old = np.copy(u_new)
        
        iter += 1
        if iter > MaxI:
            break
        if iter < MinI:
            GError = 1.0

    return u_new, apu, AE, AW, AN, AS

# @jit(nopython=True)
def V_Solver0(u, v, p, Nx, Ny, data, WSOR, rhoold, rho, mu, T):
    dx = data['dx']
    dy = data['dy']
    hx = data['hx']
    hy = data['hy']
    theta = data['theta']
    rp = data['rp']
    cita = data['cita']
    T0 = data['Tf']
    dt = data['dt']
    
    MaxG = WSOR['MaxG']
    MaxI = WSOR['MaxU']
    MinI = WSOR['MinU']
    
    AP = np.zeros((Nx, Ny))
    AE = np.zeros((Nx, Ny))
    AW = np.zeros((Nx, Ny))
    AN = np.zeros((Nx, Ny))
    AS = np.zeros((Nx, Ny))
    SP = np.zeros((Nx, Ny))
    apv = np.ones((Nx, Ny))
    v_new = np.copy(v)
    v_old = np.copy(v)
    
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
            rhoe = 0.25 * (rho[i, j] + rho[i+1, j] + rho[i, j-1] + rho[i+1, j-1])
            rhow = 0.25 * (rho[i, j] + rho[i-1, j] + rho[i-1, j-1] + rho[i, j-1])
            rhon = rho[i, j]
            rhos = rho[i, j-1]
            mue = 0.25 * (mu[i, j] + mu[i+1, j] + mu[i, j-1] + mu[i+1, j-1])
            muw = 0.25 * (mu[i, j] + mu[i-1, j] + mu[i-1, j-1] + mu[i, j-1])
            mun = mu[i, j]
            mus = mu[i, j-1]
            
            De = (mue * dy[i, j] * hx[i, j])
            Dw = (muw * dy[i, j] * hx[i, j])
            Dn = (mun * dx[i, j] * hy[i, j])
            Ds = (mus * dx[i, j] * hy[i, j])
            
            Fe = 0.5 * rhoe * (u[i+1, j] + u[i+1, j-1]) * dy[i, j]
            Fw = 0.5 * rhow * (u[i, j] + u[i, j-1]) * dy[i, j]
            Fn = 0.5 * rhon * (v_old[i, j+1] + v_old[i, j]) * dx[i, j]
            Fs = 0.5 * rhos * (v_old[i, j-1] + v_old[i, j]) * dx[i, j]
            
            me = max(-Fe, 0)
            mw = max(Fw, 0)
            mn = max(-Fn, 0)
            ms = max(Fs, 0)
            cds_upw = (0.5 * (Fe * (v[i, j] + v[i+1, j]) - Fw * (v[i, j] + v[i-1, j]) +
                              Fn * (v[i, j] + v[i, j+1]) - Fs * (v[i, j] + v[i, j-1])) -
                       ((me + mw + mn + ms) * v[i, j] - (me * v[i+1, j] + mw * v[i-1, j] +
                                                        mn * v[i, j+1] + ms * v[i, j-1])))
            
            Dp = (p[i, j-1] - p[i, j]) * dx[i, j]
            Sc = (mu[i, j] * dy[i, j] / rp[i, j] * (u[i, j-1] + u[i, j] - u[i+1, j] - u[i+1, j-1]) +
                  rhoold[i, j] * (0.25 * (u[i, j-1] + u[i, j] + u[i+1, j] + u[i+1, j-1]))**2 *
                  (theta[i+1, j] - theta[i, j]) * dy[i, j] -
                  rhoold[i, j] * (1-cita * (T[i, j] - T0)) * (9.8) * np.cos(theta[i, j]) * dx[i, j] * dy[i, j])
            Sp = -mu[i, j] * dy[i, j] * (theta[i+1, j] - theta[i, j]) / rp[i, j]
            Spmax = max(-Sp, 0)
            
            ae = me + De
            aw = mw + Dw
            an = mn + Dn
            as_ = ms + Ds
            if i == 1:   
                aw = 0  # Zero gradient
            if i == Nx-2:
                ae = 0  # Zero gradient
            
            ap0 = rhoold[i, j] * dy[i, j] * dx[i, j] / dt
            ap = ae + aw + an + as_ + Spmax + ap0
            ap /= WSOR['alpha_v']
            
            AP[i, j] = ap
            AE[i, j] = -ae
            AW[i, j] = -aw
            AN[i, j] = -an
            AS[i, j] = -as_
            SP[i, j] = Dp + Sc - WSOR['beta'] * cds_upw + (Sp + Spmax) * v_old[i, j] + ap0 * v_old[i, j]
            apv[i, j] = ap
    
    GError = 1.0
    iter = 1
    while GError > MaxG:
        for i in range(1, Nx-1):
            for j in range(1, Ny-1):
                Be = AE[i, j] * v_new[i+1, j]
                Bw = AW[i, j] * v_new[i-1, j]
                Bn = AN[i, j] * v_new[i, j+1]
                Bs = AS[i, j] * v_new[i, j-1]
                v_new[i, j] = (1 - WSOR['alpha_v']) * v_old[i, j] + (SP[i, j] - Be - Bw - Bn - Bs) / AP[i, j]
        
        GError = np.sqrt(np.sum((v_new - v_old)**2) / (Nx * Ny))
        #北边界
        v_new[:, Ny-1] = 0.0
        #南边界
        v_new[:, 0] = 0.0
        #东边界
        v_new[Nx-1, :] = v_new[Nx-2, :]
        #西边界
        v_new[0, :] = v_new[1, :]

        v_old = np.copy(v_new)
        
        iter += 1
        if iter > MaxI:
            break
        if iter < MinI:
            GError = 1.0
    
    return v_new, apv, AE, AW, AN, AS

# @jit(nopython=True)
def T_Solver0(u, v, T0, Nx, Ny, data, rho, kappa, WSOR):
    dx = data['dx']
    dy = data['dy']
    hx = data['hx']
    hy = data['hy']
    dt = data['dt']
    Nu = u.shape
    
    MaxG = WSOR['MaxTG']
    MaxI = WSOR['MaxT']
    MinI = WSOR['MinT']
    
    AP = np.zeros((Nx, Ny))
    AE = np.zeros((Nx, Ny))
    AW = np.zeros((Nx, Ny))
    AN = np.zeros((Nx, Ny))
    AS = np.zeros((Nx, Ny))
    SP = np.zeros((Nx, Ny))
    
    T_new = np.copy(T0)
    T_old = np.copy(T0)
    
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
            if i < Nu[0]-1 and j < Nu[1]-1:
                rhoe = 0.5 * (rho[i, j] + rho[i+1, j])
                rhow = 0.5 * (rho[i, j] + rho[i-1, j])
                rhon = 0.5 * (rho[i, j] + rho[i, j+1])
                rhos = 0.5 * (rho[i, j] + rho[i, j-1])
                Fe = rhoe * u[i+1, j] * dy[i, j]
                Fw = rhow * u[i, j] * dy[i, j]
                Fn = rhon * v[i, j+1] * dx[i, j]
                Fs = rhos * v[i, j] * dx[i, j]
                
                me = max(-Fe, 0)
                mw = max(Fw, 0)
                mn = max(-Fn, 0)
                ms = max(Fs, 0)
            else:
                me, mw, mn, ms = 0, 0, 0, 0
            
            De = (kappa[i, j] * dy[i, j] * hx[i, j])
            Dw = (kappa[i, j] * dy[i, j] * hx[i, j])
            Dn = (kappa[i, j] * dx[i, j] * hy[i, j])
            Ds = (kappa[i, j] * dx[i, j] * hy[i, j])
            
            ae = me + De
            aw = mw + Dw
            an = mn + Dn
            as_ = ms + Ds
            
            if j == Ny-2:  # North boundary
                an = 0
                Spmax = dx[i, j] / (1 / data['alpha'] + dy[i, j] / data['lammda'])
                Sc = Spmax * data['Tenv']
            else:
                Sc = 0
                Spmax = 0
            
            if j == 1:  # South boundary
                as_ = 0
                aw = 0
                ae = 0
            
            if i == 1:  # West boundary
                aw = 0
            
            if i == Nx-2:  # East boundary
                ae = 0
            
            ap0 = rho[i, j] * dy[i, j] * dx[i, j] / dt
            ap = ae + aw + an + as_ + Spmax + ap0
            ap /= WSOR['alpha_u']
            
            Dp = 0
            
            AP[i, j] = ap
            AE[i, j] = -ae
            AW[i, j] = -aw
            AN[i, j] = -an
            AS[i, j] = -as_
            SP[i, j] = Dp + Sc + ap0 * T_old[i, j]
    
    GError = 1.0
    iter = 1
    while GError > MaxG:
        for i in range(1, Nx-1):
            for j in range(1, Ny-1):
                Be = AE[i, j] * T_new[i+1, j]
                Bw = AW[i, j] * T_new[i-1, j]
                Bn = AN[i, j] * T_new[i, j+1]
                Bs = AS[i, j] * T_new[i, j-1]
                T_new[i, j] = (1 - WSOR['alpha_u']) * T_old[i, j] + (SP[i, j] - Be - Bw - Bn - Bs) / AP[i, j]
        
        GError = np.sum((T_new - T_old)**2) / (Nx * Ny)
        #北边界
        T_new[:, Ny-1] = data['Tenv']
        #南边界
        T_new[:, 0] = T_new[:, 1]
        #东边界
        T_new[Nx-1, :] = T_new[Nx-2, :]
        #西边界
        T_new[0, :] = T_new[1, :]
        
        T_old = np.copy(T_new)
        
        iter += 1
        if iter > MaxI:
            break
        if iter < MinI:
            GError = 1.0
    
    return T_new, GError, iter


# def P_Solver(us, apu, vs, apv, Nx, Ny, data, WSOR, rhoold, rho):
#     SType = WSOR['PType']

#     dx = data['dx']
#     dy = data['dy']
#     dt = data['dt']
#     dx2 = data['dx2']
#     dy2 = data['dy2']

#     MaxG = WSOR['MaxPG']
#     MaxI = WSOR['MaxP']
#     MinI = WSOR['MinP']

#     AE = np.zeros((Nx, Ny))
#     AW = np.zeros((Nx, Ny))
#     AN = np.zeros((Nx, Ny))
#     AS = np.zeros((Nx, Ny))
#     SP = np.zeros((Nx, Ny))
#     pprime = np.zeros((Nx, Ny))
#     pprime_old = np.copy(pprime)

#     # Coefficient calculations
#     for i in range(1, Nx-1):
#         for j in range(1, Ny-1):
#             rhoe = 0.5 * (rho[i, j] + rho[i+1, j])
#             rhow = 0.5 * (rho[i, j] + rho[i-1, j])
#             rhon = 0.5 * (rho[i, j] + rho[i, j+1])
#             rhos = 0.5 * (rho[i, j] + rho[i, j-1])

#             um = dy[i, j] * (rhoe * us[i+1, j] - rhow * us[i, j])
#             vm = dx[i, j] * (rhon * vs[i, j+1] - rhos * vs[i, j])

#             ae = rhoe * dy2[i, j] / apu[i+1, j]
#             aw = rhow * dy2[i, j] / apu[i, j]
#             an = rhon * dx2[i, j] / apv[i, j+1]
#             as_ = rhos * dx2[i, j] / apv[i, j]

#             AE[i, j] = -ae
#             AW[i, j] = -aw
#             AN[i, j] = -an
#             AS[i, j] = -as_
#             SP[i, j] = -um - vm + (rhoold[i, j] - rho[i, j]) * dx[i, j] * dy[i, j] / dt

#     # Boundary conditions
#     AE[Nx-2, :] = 0.0
#     AW[1, :] = 0.0
#     AN[:, Ny-2] = 0.0
#     AS[:, 1] = 0.0

#     AP = -(AE + AW + AN + AS)
#     AP[1, 1] = 1.0e40
#     AE[1, 1] = 0.0
#     AW[1, 1] = 0.0
#     AN[1, 1] = 0.0
#     AS[1, 1] = 0.0
#     SP[1, 1] = 0.0

#     if SType == 1:
#         # Sparse matrix setup for Conjugate Gradient method
#         X_sparse = []
#         Y_sparse = []
#         A_sparse = []
#         B = np.zeros((Nx * Ny, 1))

#         Isum = 0
#         for i in range(1, Nx-1):
#             for j in range(1, Ny-1):
#                 ip = (j-1) * Nx + i
#                 X_sparse.append(ip)
#                 Y_sparse.append(ip)
#                 A_sparse.append(AP[i, j])

#                 X_sparse.append(ip)
#                 Y_sparse.append(ip+1)
#                 A_sparse.append(AE[i, j])

#                 X_sparse.append(ip)
#                 Y_sparse.append(ip-1)
#                 A_sparse.append(AW[i, j])

#                 X_sparse.append(ip)
#                 Y_sparse.append(ip+Nx)
#                 A_sparse.append(AN[i, j])

#                 X_sparse.append(ip)
#                 Y_sparse.append(ip-Nx)
#                 A_sparse.append(AS[i, j])

#                 B[ip, 0] = SP[i, j]

#         # Sparse matrix in COO format
#         A = coo_matrix((A_sparse, (X_sparse, Y_sparse)), shape=(Nx * Ny, Nx * Ny))

#         # Conjugate Gradient solver
#         Xp, info = cg(A, B.flatten(), tol=MaxG, maxiter=MaxI)

#         if info != 0:
#             print(f"Conjugate gradient did not converge: info = {info}")

#         # Reshape solution back into pprime
#         for i in range(Nx):
#             for j in range(Ny):
#                 ip = (j-1) * Nx + i
#                 pprime[i, j] = Xp[ip]

#     else:
#         GError = 1.0
#         iter_ = 1
#         while GError > MaxG:
#             for i in range(1, Nx-1):
#                 for j in range(1, Ny-1):
#                     Be = AE[i, j] * pprime[i+1, j]
#                     Bw = AW[i, j] * pprime[i-1, j]
#                     Bn = AN[i, j] * pprime[i, j+1]
#                     Bs = AS[i, j] * pprime[i, j-1]
#                     pprime[i, j] = (SP[i, j] - Be - Bw - Bn - Bs) / AP[i, j]

#             GError = np.sum((pprime - pprime_old)**2) / (Nx * Ny)
#             pprime_old = np.copy(pprime)

#             iter_ += 1
#             if iter_ > MaxI:
#                 break
#             if iter_ < MinI:
#                 GError = 1.0

#     return pprime


def P_Solver0(us, apu, vs, apv, Nx, Ny, data, WSOR, rhoold, rho):
    SType = WSOR['PType']

    dx = data['dx']
    dy = data['dy']
    dt = data['dt']
    dx2 = data['dx2']
    dy2 = data['dy2']

    MaxG = WSOR['MaxPG']
    MaxI = WSOR['MaxP']
    MinI = WSOR['MinP']

    # 使用 NumPy 创建矩阵
    AE = np.zeros((Nx, Ny))
    AW = np.zeros((Nx, Ny))
    AN = np.zeros((Nx, Ny))
    AS = np.zeros((Nx, Ny))
    SP = np.zeros((Nx, Ny))
    pprime = np.zeros((Nx, Ny))
    pprime_old = np.copy(pprime)
    pprime[0,:] = pprime[1,:]
    pprime[-1,:] = pprime[-2,:]
    pprime[:,0] = pprime[:,1]
    pprime[:,-1] = pprime[:,-2]

    # # 使用向量化操作来计算矩阵系数，代替手动循环
    # rho_e = 0.5 * (rho[:-1, :] + rho[1:, :])[:,:Ny]
    # rho_w = 0.5 * (rho[1:, :] + rho[:-1, :])[:,:Ny]
    # rho_n = 0.5 * (rho[:, :-1] + rho[:, 1:])[:,1:Ny]
    # rho_s = 0.5 * (rho[:, 1:] + rho[:, :-1])[:,1:Ny]

    # um = dy[1:,:Ny] * (rho_e * us[1:, :] - rho_w * us[:-1, :])
    # vm = dx[:,1:Ny] * (rho_n * vs[:, 1:] - rho_s * vs[:, :-1])

    # ae = rho_e * dy2[1:, :] / apu[1:, :]
    # aw = rho_w * dy2[:-1, :] / apu[:-1, :]
    # an = rho_n * dx2[:, 1:] / apv[:, 1:]
    # as_ = rho_s * dx2[:, :-1] / apv[:, :-1]
    # AE[1:-1, 1:-1] = -ae
    # AW[1:-1, 1:-1] = -aw
    # AN[1:-1, 1:-1] = -an
    # AS[1:-1, 1:-1] = -as_
    # SP[1:-1, 1:-1] = -um[1:-1, 1:-1] - vm[1:-1, 1:-1] + (rhoold[1:-1, 1:-1] - rho[1:-1, 1:-1]) * dx[1:-1, 1:-1] * dy[1:-1, 1:-1] / dt
    for i in range(1, Nx-1):
        for j in range(1, Ny-1):
            rhoe = 0.5 * (rho[i, j] + rho[i+1, j])
            rhow = 0.5 * (rho[i, j] + rho[i-1, j])
            rhon = 0.5 * (rho[i, j] + rho[i, j+1])
            rhos = 0.5 * (rho[i, j] + rho[i, j-1])

            um = dy[i, j] * (rhoe * us[i+1, j] - rhow * us[i, j])
            vm = dx[i, j] * (rhon * vs[i, j+1] - rhos * vs[i, j])

            ae = rhoe * dy2[i, j] / apu[i+1, j]
            aw = rhow * dy2[i, j] / apu[i, j]
            an = rhon * dx2[i, j] / apv[i, j+1]
            as_ = rhos * dx2[i, j] / apv[i, j]

            AE[i, j] = -ae
            AW[i, j] = -aw
            AN[i, j] = -an
            AS[i, j] = -as_
            SP[i, j] = -um - vm + (rhoold[i, j] - rho[i, j]) * dx[i, j] * dy[i, j] / dt


    # 边界条件处理
    AE[Nx-2, :] = 0.0
    AW[1, :] = 0.0
    AN[:, Ny-2] = 0.0
    AS[:, 1] = 0.0

    AP = -(AE + AW + AN + AS)
    AP[1, 1] = 1.0e40
    SP[1, 1] = 0.0

    if SType == 1:
        # 使用预分配数组替换 append 操作，节约内存和时间
        num_entries = 5 * (Nx - 2) * (Ny - 2)
        X_sparse = np.zeros(num_entries, dtype=int)
        Y_sparse = np.zeros(num_entries, dtype=int)
        A_sparse = np.zeros(num_entries)
        B = np.zeros((Nx - 2) * (Ny - 2))

        idx = 0
        for i in range(1, Nx-1):
            for j in range(1, Ny-1):
                ip = (j-1) * (Nx-2) + (i-1)  # 修正索引，排除边界
                # 确保索引号在有效范围内
                X_sparse[idx] = ip
                Y_sparse[idx] = ip
                A_sparse[idx] = AP[i, j]
                idx += 1

                if i < Nx-2:
                    X_sparse[idx] = ip
                    Y_sparse[idx] = ip+1
                    A_sparse[idx] = AE[i, j]
                    idx += 1

                if i > 1:
                    X_sparse[idx] = ip
                    Y_sparse[idx] = ip-1
                    A_sparse[idx] = AW[i, j]
                    idx += 1

                if j < Ny-2:
                    X_sparse[idx] = ip
                    Y_sparse[idx] = ip+(Nx-2)
                    A_sparse[idx] = AN[i, j]
                    idx += 1

                if j > 1:
                    X_sparse[idx] = ip
                    Y_sparse[idx] = ip-(Nx-2)
                    A_sparse[idx] = AS[i, j]
                    idx += 1

                B[ip] = SP[i, j]

        # 使用 scipy 的稀疏矩阵构建函数，只使用有效条目
        A = coo_matrix((A_sparse[:idx], (X_sparse[:idx], Y_sparse[:idx])), shape=((Nx-2) * (Ny-2), (Nx-2) * (Ny-2)))

        # 预处理加速求解器（如不完全LU分解）
        try:
            ilu = spilu(A.tocsc())
            M = LinearOperator(A.shape, ilu.solve)

            # 使用预处理器的共轭梯度求解器
            Xp, info = cg(A, B, tol=MaxG, maxiter=MaxI, M=M)

            if info != 0:
                print(f"Conjugate gradient did not converge: info = {info}")
                raise RuntimeError("CG did not converge")

        except RuntimeError as e:
            print(f"LU decomposition or CG failed: {e}")
            # 使用 gmres 方法处理奇异矩阵，调整参数以加速收敛
            # tol = 1e-6  # 容差
            # maxiter = 1000  # 最大迭代次数
            # 使用合适的初始条件来加速 lsmr 方法的收敛
            x0 = np.zeros((Nx-2) * (Ny-2))  # 初始条件可以是零向量或其他合适的初始猜测
            # 使用 lsmr 方法处理奇异矩阵，调整参数以加速收敛
            damp = 1e-8  # 正则化参数，防止奇异矩阵导致的数值不稳定
            atol = 1e-6  # 绝对容差
            btol = 1e-6  # 相对容差
            conlim = 1e8  # 条件数限制
            Xp = lsmr(A, B, damp=damp, atol=atol, btol=btol, conlim=conlim, x0=x0)[0]
            # Xp, info = gmres(A, B, tol=tol, restart=50, maxiter=maxiter)
            # 使用 lsmr 方法处理奇异矩阵
            # Xp = lsmr(A, B)[0]
            # if info != 0:
            #     print(f"GMRES did not converge: info = {info}")
            #     raise RuntimeError("GMRES did not converge")

        # 将解重新赋值回 pprime
        # for i in range(Nx):
        #     for j in range(Ny):
        #         ip = (j) * Nx + i
        #         pprime[i, j] = Xp[ip]
        for i in range(1, Nx-1):
            for j in range(1, Ny-1):
                ip = (j-1) * (Nx-2) + (i-1)
                pprime[i, j] = Xp[ip]
        # 使用NumPy的切片操作来处理边界条件
        pprime[0, :] = pprime[1, :]
        pprime[-1, :] = pprime[-2, :]
        pprime[:, 0] = pprime[:, 1]
        pprime[:, -1] = pprime[:, -2]


    else:
        GError = 1.0
        iter_ = 1
        # 预先计算常量项
        inv_AP = 1.0 / AP
        while GError > MaxG:
            # for i in range(1, Nx-1):
            #     for j in range(1, Ny-1):
            #         Be = AE[i, j] * pprime[i+1, j]
            #         Bw = AW[i, j] * pprime[i-1, j]
            #         Bn = AN[i, j] * pprime[i, j+1]
            #         Bs = AS[i, j] * pprime[i, j-1]
            #         pprime[i, j] = (SP[i, j] - Be - Bw - Bn - Bs) / AP[i, j]
            # 使用NumPy的向量化操作来替代循环
            Be = AE[1:-1, 1:-1] * pprime[2:, 1:-1]
            Bw = AW[1:-1, 1:-1] * pprime[:-2, 1:-1]
            Bn = AN[1:-1, 1:-1] * pprime[1:-1, 2:]
            Bs = AS[1:-1, 1:-1] * pprime[1:-1, :-2]
            pprime[1:-1, 1:-1] = (SP[1:-1, 1:-1] - Be - Bw - Bn - Bs) * inv_AP[1:-1, 1:-1]
            Be = AE[1:-1, 1:-1] * pprime[2:, 1:-1]
            pprime[1:-1, 1:-1] = (SP[1:-1, 1:-1] - Be - Bw - Bn - Bs) * inv_AP[1:-1, 1:-1]
            Bw = AW[1:-1, 1:-1] * pprime[:-2, 1:-1]
            pprime[1:-1, 1:-1] = (SP[1:-1, 1:-1] - Be - Bw - Bn - Bs) * inv_AP[1:-1, 1:-1]
            Bn = AN[1:-1, 1:-1] * pprime[1:-1, 2:]
            pprime[1:-1, 1:-1] = (SP[1:-1, 1:-1] - Be - Bw - Bn - Bs) * inv_AP[1:-1, 1:-1]
            Bs = AS[1:-1, 1:-1] * pprime[1:-1, :-2]
            pprime[1:-1, 1:-1] = (SP[1:-1, 1:-1] - Be - Bw - Bn - Bs) * inv_AP[1:-1, 1:-1]

            # 使用NumPy的切片操作来处理边界条件
            pprime[0, :] = pprime[1, :]
            pprime[-1, :] = pprime[-2, :]
            pprime[:, 0] = pprime[:, 1]
            pprime[:, -1] = pprime[:, -2]
            GError = np.sum((pprime - pprime_old)**2) / (Nx * Ny)
            pprime_old = np.copy(pprime)

            iter_ += 1
            if iter_ > MaxI:
                break
            if iter_ < MinI:
                GError = 1.0

    return pprime


# @njit(parallel=True)
def Correction(u, apu, v, apv, p, pprime, NxU, NyU, NxV, NyV, NxP, NyP, rhoold, rho, alpha_p, dx, dy, dt):
    # 预分配数组
    u_new = np.copy(u)
    v_new = np.copy(v)
    p_new = np.copy(p)
    source = np.zeros((NxP, NyP))

    # 并行修正压力
    for i in prange(1, NxP-1):  # 并行化
        for j in range(1, NyP-1):
            p_new[i, j] = p[i, j] + alpha_p * pprime[i, j]

    # 并行修正速度u
    for i in prange(1, NxU-1):  # 并行化
        for j in range(1, NyU-1):
            uprime = (dy[i, j] / apu[i, j]) * (pprime[i-1, j] - pprime[i, j])
            u_new[i, j] = u[i, j] + uprime

    # 并行修正速度v
    for i in prange(1, NxV-1):  # 并行化
        for j in range(1, NyV-1):
            vprime = (dx[i, j] / apv[i, j]) * (pprime[i, j-1] - pprime[i, j])
            v_new[i, j] = v[i, j] + vprime

    # 并行计算 source 项
    for i in prange(1, NxP-1):  # 并行化
        for j in range(1, NyP-1):
            rhoe = 0.5 * (rho[i, j] + rho[i+1, j])
            rhow = 0.5 * (rho[i, j] + rho[i-1, j])
            rhon = 0.5 * (rho[i, j] + rho[i, j+1])
            rhos = 0.5 * (rho[i, j] + rho[i, j-1])

            um = dy[i, j] * (rhoe * u_new[i+1, j] - rhow * u_new[i, j])
            vm = dx[i, j] * (rhon * v_new[i, j+1] - rhos * v_new[i, j])
            source[i, j] = -um - vm + (rhoold[i, j] - rho[i, j]) * dx[i, j] * dy[i, j] / dt

    # 计算 source 项的总和
    src = np.sqrt(np.sum(source**2) / (NxP * NyP))

    return u_new, v_new, p_new, src

# @jit(nopython=True)
def Correction0(u, apu, v, apv, p, pprime, NxU, NyU, NxV, NyV, NxP, NyP, rhoold, rho,alpha_p,dx,dy,dt):

    u_new = np.copy(u)
    v_new = np.copy(v)
    p_new = np.copy(p)
    source = np.zeros((NxP, NyP))

    # 修正压力
    for i in range(1, NxP-1):
        for j in range(1, NyP-1):
            p_new[i, j] = p[i, j] + alpha_p * pprime[i, j]

    # 修正速度u
    for i in range(1, NxU-1):
        for j in range(1, NyU-1):
            uprime = (dy[i, j] / apu[i, j]) * (pprime[i-1, j] - pprime[i, j])
            u_new[i, j] = u[i, j] + uprime

    # 修正速度v
    for i in range(1, NxV-1):
        for j in range(1, NyV-1):
            vprime = (dx[i, j] / apv[i, j]) * (pprime[i, j-1] - pprime[i, j])
            v_new[i, j] = v[i, j] + vprime

    # 计算 source 项
    for i in range(1, NxP-1):
        for j in range(1, NyP-1):
            rhoe = 0.5 * (rho[i, j] + rho[i+1, j])
            rhow = 0.5 * (rho[i, j] + rho[i-1, j])
            rhon = 0.5 * (rho[i, j] + rho[i, j+1])
            rhos = 0.5 * (rho[i, j] + rho[i, j-1])

            um = dy[i, j] * (rhoe * u_new[i+1, j] - rhow * u_new[i, j])
            vm = dx[i, j] * (rhon * v_new[i, j+1] - rhos * v_new[i, j])
            source[i, j] = -um - vm + (rhoold[i, j] - rho[i, j]) * dx[i, j] * dy[i, j] / dt

    src = np.sum(source**2) / (NxP * NyP)

    return u_new, v_new, p_new, src

# @njit(parallel=True)
def Correction2(u0, u, apu, AEu, AWu, ANu, ASu, v0, v, apv, AEv, AWv, ANv, ASv, NxU, NyU, NxV, NyV):
    u_new = np.copy(u)
    v_new = np.copy(v)

    up = u - u0
    vp = v - v0

    # 修正u方向速度
    for i in prange(1, NxU-1):
        for j in range(1, NyU-1):
            Be = -AEu[i, j] * up[i+1, j]
            Bw = -AWu[i, j] * up[i-1, j]
            Bn = -ANu[i, j] * up[i, j+1]
            Bs = -ASu[i, j] * up[i, j-1]
            u_new[i, j] = u[i, j] + (Be + Bw + Bn + Bs) / apu[i, j]

    # 修正v方向速度
    for i in prange(1, NxV-1):
        for j in range(1, NyV-1):
            Be = -AEv[i, j] * vp[i+1, j]
            Bw = -AWv[i, j] * vp[i-1, j]
            Bn = -ANv[i, j] * vp[i, j+1]
            Bs = -ASv[i, j] * vp[i, j-1]
            v_new[i, j] = v[i, j] + (Be + Bw + Bn + Bs) / apv[i, j]

    return u_new, v_new
