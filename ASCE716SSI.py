# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 17:41:39 2020

@author: sebas
"""
def detalpha(R):
    """
    

    Parameters
    ----------
    R : TYPE
        DESCRIPTION.

    Returns
    -------
    alpha : TYPE
        DESCRIPTION.

    """
    if R <= 3:
        alpha = 0.70
    elif R > 3 and R < 6:
        alpha = 0.50 + R/15
    elif R >= 6:
        alpha = 0.90
        
    # Returning value of alpha.
    return alpha
    
def detWgor(Mast):
    """
    

    Parameters
    ----------
    Mast : TYPE
        DESCRIPTION.

    Returns
    -------
    Wgor : TYPE
        DESCRIPTION.

    """
    Wgor = Mast*9.81
    
    # Returning value of Wgor
    return Wgor

def dethgor(hn):
    """
    

    Parameters
    ----------
    hn : TYPE
        DESCRIPTION.

    Returns
    -------
    hgor : TYPE
        DESCRIPTION.

    """
    hgor = 0.70*hn
    # Returning value of hgor.
    return hgor

def detkgor(Mast,T):
    """
    

    Parameters
    ----------
    Mast : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.

    Returns
    -------
    kgor : TYPE
        DESCRIPTION.

    """
    import numpy as np
    
    # Previous definitions
    # --------------------
    Wgor = detWgor(Mast)
    kgor = 4*np.pi*(Wgor/(9.81*T**2))
    # Returning value of kgor
    return kgor

def detGsoil(Ss,gamma,Vs30):
    """
    

    Parameters
    ----------
    Ss : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    Vs30 : TYPE
        DESCRIPTION.

    Returns
    -------
    G : TYPE
        DESCRIPTION.
    Go : TYPE
        DESCRIPTION.
    Vs : TYPE
        DESCRIPTION.
    G_Go : TYPE
        DESCRIPTION.
    Vs_Vso : TYPE
        DESCRIPTION.
    SDS : TYPE
        DESCRIPTION.

    """
    import ASCE716Seismic as ASCE716
    # Previous calculations
    # ---------------------
    SiteClass = ASCE716.detSiteClass(Vs30)
    SDS = ASCE716.detSDS(Ss,Vs30)
    Vso = Vs30
    Go = gamma*Vso**2/9.81
    
    # Determining relations from tables 19.3-1 y 19.3-2 of ASCE7-16
    # --------------------------------------------------------------
    if SDS/2.5 == 0 and SiteClass in ('A','B','C','D','E'):
        Vs_Vso = 1.0
        G_Go = 1.0
    if SiteClass in ('A'):
        Vs_Vso = 1.0
        G_Go = 1.0
    elif SiteClass in ('B'):
        if SDS/2.5 <= 0.1 and SDS/2.5 > 0:
            Vs_Vso = 1.0
            G_Go = 1.0
        elif SDS/2.5 > 0.1 and SDS/2.5 <= 0.40:
            Vs_Vso = 1.00 + (SDS/2.5 - 0.10)*(0.97 - 1.00)/(0.40 - 0.10)
            G_Go = 1.00 + (SDS/2.5 - 0.10)*(0.95 - 1.00)/(0.40 - 0.10)
        elif SDS/2.5 > 0.40 and SDS/2.5 < 0.80:
            Vs_Vso = 0.97 + (SDS/2.5 - 0.40)*(0.95 - 0.97)/(0.8 - 0.40)
            G_Go = 0.95 + (SDS/2.5 - 0.40)*(0.90 - 0.95)/(0.80 - 0.40)
        elif SDS/2.5 >= 0.80:
            Vs_Vso, G_Go = 0.95, 0.90
    elif SiteClass in ('C'):
        if SDS/2.5 <= 0.1 and SDS/2.5 > 0:
            Vs_Vso = 0.97
            G_Go = 0.95
        elif SDS/2.5 > 0.1 and SDS/2.5 <= 0.40:
            Vs_Vso = 0.97 + (SDS/2.5 - 0.10)*(0.87 - 0.97)/(0.40 - 0.10)
            G_Go = 0.95 + (SDS/2.5 - 0.10)*(0.75 - 0.95)/(0.40 - 0.10)
        elif SDS/2.5 > 0.40 and SDS/2.5 < 0.80:
            Vs_Vso = 0.87 + (SDS/2.5 - 0.40)*(0.77 - 0.87)/(0.8 - 0.40)
            G_Go = 0.75 + (SDS/2.5 - 0.40)*(0.60 - 0.75)/(0.80 - 0.40)
        elif SDS/2.5 >= 0.80:
            Vs_Vso, G_Go = 0.77, 0.60
    elif SiteClass in ('D'):
        if SDS/2.5 <= 0.1 and SDS/2.5 > 0:
            Vs_Vso = 0.95
            G_Go = 0.90
        elif SDS/2.5 > 0.1 and SDS/2.5 <= 0.40:
            Vs_Vso = 0.95 + (SDS/2.5 - 0.10)*(0.71 - 0.95)/(0.40 - 0.10)
            G_Go = 0.90 + (SDS/2.5 - 0.10)*(0.50 - 0.90)/(0.40 - 0.10)
        elif SDS/2.5 > 0.40 and SDS/2.5 < 0.80:
            Vs_Vso = 0.71 + (SDS/2.5 - 0.40)*(0.32 - 0.71)/(0.8 - 0.40)
            G_Go = 0.50 + (SDS/2.5 - 0.40)*(0.22 - 0.50)/(0.80 - 0.40)
        elif SDS/2.5 >= 0.80:
            Vs_Vso, G_Go = 0.32, 0.22
    elif SiteClass in ('E'):
        if SDS/2.5 <= 0.1 and SDS/2.5 > 0:
            Vs_Vso = 0.77
            G_Go = 0.60
        elif SDS/2.5 > 0.1 and SDS/2.5 <= 0.40:
            Vs_Vso = 0.77 + (SDS/2.5 - 0.10)*(0.22 - 0.77)/(0.40 - 0.10)
            G_Go = 0.60 + (SDS/2.5 - 0.10)*(0.17 - 0.60)/(0.40 - 0.10)
        elif SDS/2.5 > 0.40 and SDS/2.5 < 0.80:
            Vs_Vso = 0.22 + (SDS/2.5 - 0.40)*(0.13 - 0.22)/(0.8 - 0.40)
            G_Go = 0.17 + (SDS/2.5 - 0.40)*(0.10 - 0.17)/(0.80 - 0.40)
        elif SDS/2.5 >= 0.80:
            Vs_Vso, G_Go = 0.13, 0.10
    elif SiteClass in ('F'):
            print('More detailed soil testing required!')
    G = Go*G_Go
    Vs = Vso*Vs_Vso
    
    # Returnig value for Gsoil and other parameters.
    return G,Go,Vs,G_Go,Vs_Vso,SDS

def detKiisur(Ss,Vs30,gamma,nu,B,L,D=0.,omega=0.):
    """
    

    Parameters
    ----------
    Ss : TYPE
        DESCRIPTION.
    Vs30 : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    nu : TYPE
        DESCRIPTION.
    B : TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.

    Returns
    -------
    Kzsur : TYPE
        DESCRIPTION.
    Kysur : TYPE
        DESCRIPTION.
    Kxsur : TYPE
        DESCRIPTION.
    Kzzsur : TYPE
        DESCRIPTION.
    Kyysur : TYPE
        DESCRIPTION.
    Kxxsur : TYPE
        DESCRIPTION.

    """
    # Previous calculations
    # ----------------------
    (G,Go,Vs,G_Go,Vs_Vso,SDS) = detGsoil(Ss,gamma,Vs30)
    
    # Using Pais abd Kausel 1988 expresions for determining 
    # Static Stiffness of rigid footings at  ground surface.
    # =======================================================================
    # Translation along Z axis.
    # -------------------------
    Kzsur = G*(0.50*B)/(1 - nu)*(3.1*(L/B)**0.75 + 1.60)
    # Translation along Y axis.
    # -------------------------
    Kysur = G*(0.50*B)/(2 - nu)*(6.80*(L/B)**0.65 + 0.80*(L/B) + 1.60)
    # Translation along X axis.
    # -------------------------
    Kxsur = G*(0.50*B)/(2 - nu)*(6.80*(L/B)**0.65 + 2.4)
    # Torsion about Z axis.
    # ---------------------
    Kzzsur = G*(0.50*B)**3*(4.25*(L/B)**2.45 + 4.06)
    # Rocking about Y axis.
    # ---------------------
    Kyysur = G*(0.50*B)**3/(1 - nu)*(3.73*(L/B)**2.4 + 0.27)
    # Rocking about X axis.
    # ---------------------
    Kxxsur = G*(0.50*B)**3/(1 - nu)*(3.2*(L/B) + 0.80)
    
    if D != 0.:
        # Pais and Kausel 1988  Embedment correction factors for Static Stiffness
        # and Rigid Footings to determine the value of eta_i.
        # =====================================================================
        # Translation along Z axis.
        # -------------------------
        eta_z = (1 + (0.25 + 0.25/(L/B))*(D/(0.5*B))**0.80)
        # Translation along Y axis.
        # -------------------------
        eta_y = (1 + (0.33 + 1.34/(1 + (L/B)))*(D/(0.5*B))**0.8)
        # Translation along X axis.
        # -------------------------
        eta_x = eta_y
        # Torsion about Z axis.
        # ---------------------
        eta_zz = (1 + (1.3 + 1.32/(L/B))*(D/(0.5*B))**0.9)
        # Rocking about Y axis.
        # ---------------------
        eta_yy = (1. +  D/(0.5*B) + (1.60/(0.35 + (L/B)**4))*(D/(0.5*B))**2)
        # Rocking about X axis.
        # ---------------------
        eta_xx = (1. +  D/(0.5*B) + (1.60/(0.35 + (L/B)))*(D/(0.5*B))**2)
    else:
        eta_z, eta_y, eta_x, eta_zz, eta_xx, eta_yy = 1., 1., 1., 1., 1., 1.
    
    if omega != 0.:
        # Dynamic Stiffness modifiers for rigid footings as per Pais and Kausel
        # ======================================================================
        a_o = omega*(0.5*B)/Vs
        # Translation along Z axis.
        # -------------------------
        alpha_z = 1. - ((0.4 + 0.2/(L/B))*a_o**2/((10/(1+3*((L/B) - 1))) + a_o**2))
        # Translation along Y axis.
        # -------------------------
        alpha_y = 1.
        # Translation along X axis.
        # -------------------------
        alpha_x = 1. 
        # Torsion about Z axis.
        # ---------------------
        alpha_zz = 1. - ((0.33 - 0.03*((L/B) - 1.)**0.5)*a_o**2/((0.80/(1. + 0.33*((L/B) - 1))) + a_o**2))
        # Rocking about Y axis.
        # ---------------------
        alpha_yy = 1. - (0.55*a_o**2/((0.6 + 1.4/(L/B)**3) + a_o**2))
        # Rocking about X axis.
        # ---------------------
        alpha_xx = 1. - ((0.55 + 0.01*((L/B) - 1.)**0.5)*a_o**2/((2.4 - 0.4/(L/B)**3) + a_o**2))
    else:
        alpha_z, alpha_y, alpha_x, alpha_zz, alpha_yy, alpha_xx = 1., 1., 1., 1., 1., 1.
        
    # Stiffness parameters modified for embedment and dynamic uses.
    # =============================================================
    kz = Kzsur*eta_z*alpha_z
    ky = Kysur*eta_y*alpha_y
    kx = Kxsur*eta_x*alpha_x
    kzz = Kzzsur*eta_zz*alpha_zz
    kyy = Kyysur*eta_yy*alpha_yy
    kxx = Kxxsur*eta_xx*alpha_xx
    
    # Returning values of stiffness.
    return kz, ky, kx, kzz, kyy, kxx

def detTgor(Ss,Vs30,gamma,nu,B,L,hn,Mast,T,direction,D=0,omega=0):
    """
    

    Parameters
    ----------
    Ss : TYPE
        DESCRIPTION.
    Vs30 : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    nu : TYPE
        DESCRIPTION.
    B : TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.
    hn : TYPE
        DESCRIPTION.
    Mast : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    direction : TYPE
        DESCRIPTION.

    Returns
    -------
    Tgor : TYPE
        DESCRIPTION.

    """
    # Previous definitions
    # --------------------
    hgor = dethgor(hn)
    kgor = detkgor(Mast,T)
    (_,ky,kx,_,kyy,kxx) = \
        detKiisur(Ss,Vs30,gamma,nu,B,L,D,omega)
        
    if direction in ('x','X'):
        Ky = kx
        Ktheta = kyy
    elif direction in ('z','Z'):
        Ky = ky
        Ktheta = kxx
    
    # Calculation of Tgor value
    # -------------------------
    Tgor = T*(1 + kgor/Ky*(1 + Ky*hgor**2/Ktheta))**(1/2)
    
    # Returning value of  Tgor
    return Tgor

def detTg_T_eff(Ss,R,Omega_o,Vs30,gamma,nu,B,L,hn,Mast,T,direction,Tg=0,D=0,omega=0):
    """
    

    Parameters
    ----------
    Ss : TYPE
        DESCRIPTION.
    R : TYPE
        DESCRIPTION.
    Omega_o : TYPE
        DESCRIPTION.
    Vs30 : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    nu : TYPE
        DESCRIPTION.
    B : TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.
    hn : TYPE
        DESCRIPTION.
    Mast : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    direction : TYPE
        DESCRIPTION.

    Returns
    -------
    Tg_T_eff : TYPE
        DESCRIPTION.

    """
    # Previous calculations
    # ---------------------
    mu = R/Omega_o
    if Tg == 0.:
        Tg = detTgor(Ss,Vs30,gamma,nu,B,L,hn,Mast,T,direction,D,omega)
    # Calculation of Tgor value.
    # --------------------------
    Tg_T_eff = (1 + 1/mu*((Tg/T)**2 - 1))**0.5
    
    # Returning Tgor value
    return Tg_T_eff

def detbeta_s(Ss,Vs30):
    """
    

    Parameters
    ----------
    Ss : TYPE
        DESCRIPTION.
    Vs30 : TYPE
        DESCRIPTION.

    Returns
    -------
    beta_s : TYPE
        DESCRIPTION.

    """
    import ASCE716Seismic as ASCE716
    # Previos calculations
    # ---------------------
    SDS = ASCE716.detSDS(Ss, Vs30)
    SiteClass = ASCE716.detSiteClass(Vs30)
    
    # Calculationof beta_s
    # --------------------
    # Otabined from Table 19.3-3 of ASCE7-16.
    if SDS/2.5 == 0 and SiteClass in ('C','D','E'):
        beta_s = 0.01
    elif SiteClass in ('C'):
        if SDS/2.5 > 0. and SDS/2.5 <= 0.1:
            beta_s = 0.01
        elif SDS/2.5 > 0.1 and SDS/2.5 <= 0.40:
            beta_s = 0.01 + (SDS/2.5 - 0.10)*(0.03 - 0.01)/(0.40 - 0.10)
        elif SDS/2.5 > 0.40 and SDS/2.5 < 0.80:
            beta_s = 0.03 + (SDS/2.5 - 0.40)*(0.05 - 0.03)/(0.8 - 0.40)
        elif SDS/2.5 >= 0.80:
            beta_s = 00.05
    elif SiteClass in ('D'):
        if SDS/2.5 > 0. and SDS/2.5 <= 0.1:
            beta_s =  0.01 + (SDS/2.5 - 0.0)*(0.02 - 0.01)/(0.1 - 0.0)
        elif SDS/2.5 > 0.1 and SDS/2.5 <= 0.40:
            beta_s = 0.02 + (SDS/2.5 - 0.10)*(0.07 - 0.02)/(0.40 - 0.10)
        elif SDS/2.5 > 0.40 and SDS/2.5 < 0.80:
            beta_s = 0.07 + (SDS/2.5 - 0.40)*(0.15 - 0.07)/(0.8 - 0.40)
        elif SDS/2.5 >= 0.80:
            beta_s = 0.15
    elif SiteClass in ('E'):
        if SDS/2.5 > 0. and SDS/2.5 <= 0.1:
            beta_s =  0.01 + (SDS/2.5 - 0.0)*(0.05 - 0.01)/(0.1 - 0.0)
        elif SDS/2.5 > 0.1 and SDS/2.5 <= 0.40:
            beta_s = 0.05 + (SDS/2.5 - 0.10)*(0.2 - 0.05)/(0.40 - 0.10)
        elif SDS/2.5 > 0.40 and SDS/2.5 < 0.80:
            beta_s = 0.20 + (SDS/2.5 - 0.40)*(0.38 - 0.20)/(0.8 - 0.40)
        elif SDS/2.5 >= 0.80:
            beta_s = 0.38

    # Returning value of beta_s
    return beta_s

def detbeta_rd(Ss,B,L,hn,T,Mast,Vs30,gamma,nu,direction,Tg=0.,D=0,omega=0):
    """
    

    Parameters
    ----------
    Ss : TYPE
        DESCRIPTION.
    B : TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.
    hn : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    Mast : TYPE
        DESCRIPTION.
    Vs30 : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    nu : TYPE
        DESCRIPTION.
    direction : TYPE
        DESCRIPTION.

    Returns
    -------
    beta_rd : TYPE
        DESCRIPTION.

    """
    import numpy as np
    # Previous calculations 
    # ----------------------
    (G,_,Vs,_,_,_) = detGsoil(Ss,gamma,Vs30)
    Kxx = G*(0.5*B)**3/(1 - nu)*(3.2*(L/B) + 0.80)                                              # (19.3-9)
    if Tg == 0:
        Tg = detTgor(Ss,Vs30,gamma,nu,B,L,hn,Mast,T,direction,D,omega)
    a_o = 2*np.pi*(0.5*B)/(Tg*Vs)                                                               # (19.3-11)
    alpha_xx = 1.0 - ((0.55 + 0.01*((L/B) - 1.0)**0.5)*a_o**2/((2.4 - 0.4/(L/B)**3) + a_o**2))  # (19.3-14)
    hgor = dethgor(hn)
    Txx = 2*np.pi*((Mast*hgor**2)/(alpha_xx*Kxx))**0.5                                          # (19.3-7)
    Ky = G*(0.5*B)/(2 - nu)*(6.8*(L/B)**0.65 + 0.80*(L/B) + 1.60)                               # (19.3-8)
    Ty = 2*np.pi*(Mast/Ky)**0.5                                                                 # (19.3-6)
    beta_y = (4*(L/B)/(Ky/(G*0.5*B)))*(a_o/2)                                                   # (19.3-10)
    PHI = (2*(1 - nu)/(1 - 2*nu))**0.5                                                          # (19.3-13)
    if PHI > 2.5:
        PHI = 2.5
    beta_xx = ((4*PHI/3)*(L/B)*a_o**2/((Kxx/(G*(0.5*B)**3))*((2.2 - 0.4/(L/B)**3) + a_o**2)))*(a_o/(2*alpha_xx))   # (19.3-12)
    # Calculation of beta_rd
    # -----------------------
    # Calculation according to equations 19.3-5 to 19.3-14 of ASCE 7-16.
    beta_rd = 1/(Tg/Ty)**2*beta_y + 1/(Tg/Txx)**2*beta_xx   #  (19.3-5)
    
    # returning value of beta_rd
    return beta_rd

def detbeta_f(Ss,B,L,hn,T,Mast,Vs30,gamma,nu,direction,Tg=0,D=0,omega=0):
    """
    

    Parameters
    ----------
    Ss : TYPE
        DESCRIPTION.
    B : TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.
    hn : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    Mast : TYPE
        DESCRIPTION.
    Vs30 : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    nu : TYPE
        DESCRIPTION.
    direction : TYPE
        DESCRIPTION.

    Returns
    -------
    beta_f : TYPE
        DESCRIPTION.

    """
    # Previous calculations
    # ----------------------
    if Tg == 0:
        Tg = detTgor(Ss,Vs30,gamma,nu,B,L,hn,Mast,T,direction,D,omega)
    beta_s = detbeta_s(Ss,Vs30)
    beta_rd = detbeta_rd(Ss,B,L,hn,T,Mast,Vs30,gamma,nu,direction,Tg,D,omega)
    # Calculation of beta_f
    # ---------------------
    beta_f = (((Tg/T)**2 - 1)/(Tg/T)**2)*beta_s + beta_rd
    
    # Returning value of beta_f
    return beta_f

def detbeta_o(Ss,R,Omega_o,Vs30,gamma,nu,B,L,hn,Mast,T,direction,beta=0.05,Tg=0,D=0,omega=0):
    """
    

    Parameters
    ----------
    Ss : TYPE
        DESCRIPTION.
    R : TYPE
        DESCRIPTION.
    Omega_o : TYPE
        DESCRIPTION.
    Vs30 : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    nu : TYPE
        DESCRIPTION.
    B : TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.
    hn : TYPE
        DESCRIPTION.
    Mast : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    direction : TYPE
        DESCRIPTION.
    beta : TYPE, optional
        DESCRIPTION. The default is 0.05.

    Returns
    -------
    beta_o : TYPE
        DESCRIPTION.

    """
    # Previous definitions.
    # ----------------------
    Tg_T_eff = detTg_T_eff(Ss,R,Omega_o,Vs30,gamma,nu,B,L,hn,Mast,T,direction,Tg,D,omega)
    beta_f = detbeta_f(Ss,B,L,hn,T,Mast,Vs30,gamma,nu,direction,Tg,D,omega)
    # Determining value of beta_o
    # ----------------------------
    beta_o = beta_f + beta/(Tg_T_eff)**2
    if beta_o > 0.20:
        beta_o = 0.20
    elif beta_o <= 0.20:
        beta_o = beta_o
        
    # Returning vaule of beta_o
    return beta_o

def detB_SSI(Ss,R,Omega_o,Vs30,gamma,nu,B,L,hn,Mast,T,direction,beta=0.05,Tg=0,D=0,omega=0):
    """
    

    Parameters
    ----------
    Ss : TYPE
        DESCRIPTION.
    R : TYPE
        DESCRIPTION.
    Omega_o : TYPE
        DESCRIPTION.
    Vs30 : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    nu : TYPE
        DESCRIPTION.
    B : TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.
    hn : TYPE
        DESCRIPTION.
    Mast : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    direction : TYPE
        DESCRIPTION.
    beta : TYPE, optional
        DESCRIPTION. The default is 0.05.

    Returns
    -------
    B_SSI : TYPE
        DESCRIPTION.

    """
    import numpy as np
    
    # Previous calculations
    # ---------------------
    beta_o = detbeta_o(Ss,R,Omega_o,Vs30,gamma,nu,B,L,hn,Mast,T,direction,beta,Tg,D,omega)
    B_SSI = 4*(5.6 - np.log(100*beta_o))
    
    # Retiurning value of B_SSI
    return B_SSI

def detDeltaV(Ss,S1,R,Omega_o,Ie,Vs30,gamma,nu,B,L,hn,Mast,T,direction,beta=0.05,Tg=0,D=0,omega=0):
    """
    

    Parameters
    ----------
    Ss : TYPE
        DESCRIPTION.
    S1 : TYPE
        DESCRIPTION.
    R : TYPE
        DESCRIPTION.
    Omega_o : TYPE
        DESCRIPTION.
    Ie : TYPE
        DESCRIPTION.
    Vs30 : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    nu : TYPE
        DESCRIPTION.
    B : TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.
    hn : TYPE
        DESCRIPTION.
    Mast : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    direction : TYPE
        DESCRIPTION.
    beta : TYPE, optional
        DESCRIPTION. The default is 0.05.

    Returns
    -------
    DeltaVD : TYPE
        DESCRIPTION.
    DeltaVF : TYPE
        DESCRIPTION.

    """
    import ASCE716Seismic as ASCE716
    # Previous definitions
    # -------------------
    Wgor = detWgor(Mast)
    if Tg == 0:
        Tg = detTgor(Ss,Vs30,gamma,nu,B,L,hn,Mast,T,direction,D,omega)
    B_SSI = detB_SSI(Ss,R,Omega_o,Vs30,gamma,nu,B,L,hn,Mast,T,direction,beta,Tg,D,omega)
    # Values of Cs with the fundamental fixed-base period
    # ---------------------------------------------------
    CsD = ASCE716.detCsDrifts(Ss, S1, Vs30, T, R, Ie)
    CsF = ASCE716.detCs(Ss, S1, Vs30, T, R, Ie)
    # Values of Csg, with the modified fundamental period.
    # ----------------------------------------------------
    CsDg = ASCE716.detCsDrifts(Ss, S1, Vs30, Tg, R, Ie)
    CsFg = ASCE716.detCs(Ss, S1, Vs30, Tg, R, Ie)
    
    DeltaVD = (CsD - CsDg/B_SSI)*Wgor
    DeltaVF = (CsF - CsFg/B_SSI)*Wgor
    
    # Returning DeltaV value.
    return DeltaVD, DeltaVF

def detVgor(VD,VF,Ss,S1,R,Omega_o,Ie,Vs30,gamma,nu,B,L,hn,Mast,T,direction,beta=0.05,Tg=0,D=0,omega=0):
    """
    

    Parameters
    ----------
    VD : TYPE
        DESCRIPTION.
    VF : TYPE
        DESCRIPTION.
    Ss : TYPE
        DESCRIPTION.
    S1 : TYPE
        DESCRIPTION.
    R : TYPE
        DESCRIPTION.
    Omega_o : TYPE
        DESCRIPTION.
    Ie : TYPE
        DESCRIPTION.
    Vs30 : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    nu : TYPE
        DESCRIPTION.
    B : TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.
    hn : TYPE
        DESCRIPTION.
    Mast : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    direction : TYPE
        DESCRIPTION.
    beta : TYPE, optional
        DESCRIPTION. The default is 0.05.

    Returns
    -------
    VgorD : TYPE
        DESCRIPTION.
    VgorF : TYPE
        DESCRIPTION.
    deltaVD : TYPE
        DESCRIPTION.
    deltaVF : TYPE
        DESCRIPTION.

    """
    # Previous calculations
    # ---------------------
    alpha = detalpha(R)
    (deltaVD,deltaVF) = detDeltaV(Ss,S1,R,Omega_o,Ie,Vs30,gamma,nu,B,L,hn,Mast,T,direction,beta,Tg,D,omega)
    # Reduced base shear calcultaion.
    # -------------------------------
    VgorD = VD - deltaVD
    VgorF = VF - deltaVF
    # Limit verification
    # ------------------
    if VgorD < alpha*VD:
        VgorD = alpha*VD
    else:
        VgorD = VgorD
        
    if VgorF < alpha*VF:
        VgorF = alpha*VF
    else:
        VgorF = VgorF
        
    # Returning value of reduced shear and deltaV
    return VgorD, VgorF, deltaVD, deltaVF
