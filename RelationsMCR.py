# -*- coding: utf-8 -*-
"""
Created on Wed May  6 18:35:09 2020

@author: sebas
"""

"""
    Definitions to generate Moment-Curvature-Rotation relations for 
    rectangular reinforced concrete sections.
"""

"""
    First, definitions of stress-strain behavior of reinforcing steel 
    and concrete.
    Uses relations from:
        Steel - (Bai & Au, 2011) and (Mohele) Seismic Design of Reinforced Concete
        Buildings.
        Concrete Compression - (Mohele) Seismic Design of Reinforced Concete
        Buildings, Hognestad 2015, (Kaklauskas & Ghaboussi, 2001) and Mander 1984.
        Concrete Tension - (Kaklauskas & Ghaboussi, 2001)
"""

def Resultantfct(length, b, eps_t, eps_cr, fcr, alpha1=1.0, alpha2=5.0):
    """
    

    Parameters
    ----------
    length : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    eps_t : TYPE
        DESCRIPTION.
    eps_cr : TYPE
        DESCRIPTION.
    fcr : TYPE
        DESCRIPTION.
    alpha1 : TYPE, optional
        DESCRIPTION. The default is 1.0.
    alpha2 : TYPE, optional
        DESCRIPTION. The default is 5.0.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    from scipy.integrate import quad

    # Defining function to be integrated. Stress related to discretized length
    # ------------------------------------------------------------------------
    def fct(cti):
        if eps_t <= eps_cr:
            eps_ti = eps_t*cti/length
            if eps_ti <= eps_cr:
                return fcr*eps_ti/eps_cr
            elif eps_ti > eps_cr and eps_ti <= alpha2*eps_cr:
                return alpha1*fcr*(alpha2*eps_cr - eps_ti)/(alpha2*eps_cr - eps_cr)
            else:
                return 0.0
        else:
            return 0.0
    
    # Calculating resultant force from stress distribution along tensile length.
    # --------------------------------------------------------------------------
    Res = quad(fct,0,length)                 # [res,err] resultant area in length.
    Tc = Res[0]*b                            # [MN] resultant force from volume.
    
    # Defining function to calculate area moments to locate resultant force.
    # ----------------------------------------------------------------------
    def fctmoment(cti):
        if eps_t <= eps_cr:
            eps_ti = eps_t*cti/length
            if eps_ti <= eps_cr:
                return cti*fcr*eps_ti/eps_cr
            elif eps_ti > eps_cr and eps_ti <= alpha2*eps_cr:
                return cti*alpha1*fcr*(alpha2*eps_cr - eps_ti)/(alpha2*eps_cr - eps_cr)
            else:
                return 0.0
        else:
            return 0.0
    # Calculating location of resultant with respect to neutral axis.
    # ---------------------------------------------------------------
    Res1 = quad(fctmoment,0,length)         # [res,err] resultant moment in length.
    if Res[0] == 0:
        yt = 0
    else:
        yt = Res1[0]/Res[0]                 # [m] distance from neutral axis to resultant.
    
    return Tc, yt


def Resultantfcc(length, b, eps_c, eps_0, eps_cult, fpc):
    """
    

    Parameters
    ----------
    length : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    eps_c : TYPE
        DESCRIPTION.
    eps_0 : TYPE
        DESCRIPTION.
    eps_cult : TYPE
        DESCRIPTION.
    fpc : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    """
    from scipy.integrate import quad
    
    # Defining function to be integrated. Stress related to discretized length
    # ------------------------------------------------------------------------
    def fcc(cci):
        eps_ci = eps_c*cci/length
        if eps_ci < eps_0:
            return fpc*(2*eps_ci/eps_0 - (eps_ci/eps_0)**2)
        elif eps_ci >= eps_0 and eps_ci <= eps_cult:
            return fpc*(1 - 0.15*(eps_ci - eps_0)/(eps_cult - eps_0))
    
    # Calculating resultant force from stress distribution along tensile length.
    # --------------------------------------------------------------------------
    Res = quad(fcc,0.00001,length)                 # [res,err] resultant area in length.
    Cc = Res[0]*b                            # [MN] resultant force from volume.
    
    # Defining function to calculate area moments to locate resultant force.
    # ----------------------------------------------------------------------
    def fccmoment(cci):
        eps_ci = eps_c*cci/length
        if eps_ci < eps_0:
            return cci*fpc*(2*eps_ci/eps_0 - (eps_ci/eps_0)**2)
        elif eps_ci >= eps_0 and eps_ci <= eps_cult:
            return cci*fpc*(1 - 0.15*(eps_ci - eps_0)/(eps_cult - eps_0))
    
    # Calculating location of resultant with respect to neutral axis.
    # ---------------------------------------------------------------
    Res1 = quad(fccmoment,0,length)         # [res,err] resultant moment in length.
    yc = Res1[0]/Res[0]                     # [m] distance from neutral axis to resultant.
    
    return Cc, yc

def deteps_si(eps_c,c,di):
    """
    

    Parameters
    ----------
    eps_c : TYPE
        DESCRIPTION.
    c : TYPE
        DESCRIPTION.
    di : TYPE
        DESCRIPTION.

    Returns
    -------
    eps_si : TYPE
        DESCRIPTION.

    """
    eps_si = eps_c*(c-di)/c        # [m/m] strain of steel reinforcement layer.
    
    return eps_si

def Resultantfst(Asi,eps_si,Es,Esh,fy,fsu,eps_sh,eps_su,eps_y):
    """
    

    Parameters
    ----------
    Asi : TYPE
        DESCRIPTION.
    eps_si : TYPE
        DESCRIPTION.
    Es : TYPE
        DESCRIPTION.
    Esh : TYPE
        DESCRIPTION.
    fy : TYPE
        DESCRIPTION.
    fsu : TYPE
        DESCRIPTION.
    eps_sh : TYPE
        DESCRIPTION.
    eps_su : TYPE
        DESCRIPTION.
    eps_y : TYPE
        DESCRIPTION.

    Returns
    -------
    Tsi : TYPE
        DESCRIPTION.

    """
    import numpy as np
    signo = np.sign(eps_si)                 # [adim] compresssion traction indicator.
    eps_si = abs(eps_si)                    # [m/m] absolute value of strain for calculations.
    
    # Deterination of steel stress depending on strain.
    # -------------------------------------------------
    if eps_si <= eps_y:
        fsi = eps_si*Es
    elif eps_si > eps_y and eps_si < eps_sh:
        fsi = fy
    elif eps_si >= eps_sh:
        P = Esh*(eps_su - eps_sh)/(fsu - fy)
        fsi = fsu + (fy - fsu)*((eps_su - eps_si)/(eps_su - eps_sh))**P
        
    # Calculation of foce
    # --------------------
    Tsi = signo*fsi*Asi                     # [MN] resultant force of steel layer.
    
    return Tsi


"""
    This definition script, calculates the moment-curvature-rotation
    relationship in a single layer reinforced concrete rectangular beam for My.
"""

def MyMCR(b,h,l,Asi,di,fpc,eps_cu,Es,fy,fsu,
          db = 0.025, eps_sh = 0.02,eps_su = 0.10,eps_cult = 0.01,fract = 0.0):
    """
    

    Parameters
    ----------
    b : TYPE
        DESCRIPTION.
    h : TYPE
        DESCRIPTION.
    l : TYPE
        DESCRIPTION.
    Asi : TYPE
        DESCRIPTION.
    di : TYPE
        DESCRIPTION.
    fpc : TYPE
        DESCRIPTION.
    eps_cu : TYPE
        DESCRIPTION.
    Es : TYPE
        DESCRIPTION.
    fy : TYPE
        DESCRIPTION.
    fsu : TYPE
        DESCRIPTION.
    db : TYPE, optional
        DESCRIPTION. The default is 0.025.
    eps_sh : TYPE, optional
        DESCRIPTION. The default is 0.02.
    eps_su : TYPE, optional
        DESCRIPTION. The default is 0.10.
    eps_cult : TYPE, optional
        DESCRIPTION. The default is 0.01.
    fract : TYPE, optional
        DESCRIPTION. The default is 0.0.

    Returns
    -------
    My : TYPE
        DESCRIPTION.
    curvature : TYPE
        DESCRIPTION.
    rotation : TYPE
        DESCRIPTION.
    lp : TYPE
        DESCRIPTION.

    """
    
    import unitskN_m_C as units
    
    # MECHANICAL PROPERTIES
    # ======================
    # Concrete
    # --------
    fcr = (0.62*(fpc/units.MPa)**0.5)*units.MPa         # [kPa] modulus of rupture of tension concrete.
    Ec = (4700*(fpc/units.MPa)**0.5)*units.MPa          # [kPa] modulus of elasticity of concrete.
    eps_0 = 2*fpc/Ec                                    # [m/m] strain at concrete peak stength.
    eps_cr = fcr/Ec                                     # [m/m] concrete strain at tensile fracture.
    
    # Steel
    # ------
    eps_y = fy/Es                                       # [m/m] yielding strain of steel.
    Esh = 0.05*Es                                       # [kPa] modulus at the begnining of strenght hardening region.
    
    # Initial Axial Loads
    # -------------------
    Po = 0.80*(0.85*fpc*(b*h - sum(Asi)) + sum(Asi)*fy)                # [kN] Force applied to rectangular reinforced section.
    P = fract*Po
    
    eps_step = 1e-7
    eps_cc = eps_step
    diff = -1e3
    
    while diff < 0.0 and eps_cc <= eps_cult:
        c = eps_cc/(eps_cc + eps_y)*di[0]
        eps_ct = eps_cc*(h - c)/c
        Ps = []                                                        # [kN] forces from steel layers.
        # Induced forces on concrete.
        # ----------------------------
        [Pcc,ycc] = Resultantfcc(c, b, eps_cc, eps_0, eps_cult, fpc)
        [Pct,yct] = Resultantfct(h-c, b, eps_ct, eps_cr, fcr)
        # Induced forces on reinforcing steel layers.
        # -------------------------------------------
        for j in range(len(Asi)):
            Ps.append(Resultantfst(Asi[j],deteps_si(eps_cc,c,di[j]),Es,Esh,
                                      fy,fsu,eps_sh,eps_su,eps_y))
        Pcalc = Pcc - Pct + sum(Ps) - P
        diff = Pcalc
        eps_cc += eps_step
    
    # Moment components calculation
    # -----------------------------
    Msi = []
    for j in range(len(Asi)):
        Msi.append(Resultantfst(Asi[j],deteps_si(eps_cc,c,di[j]),
                                   Es,Esh,fy,fsu,eps_sh,eps_su,eps_y)*(c - di[j]))
    # Storing Moment and curvature
    # -----------------------------
    My = Pcc*ycc + Pct*yct + sum(Msi)
    curvature = eps_ct/(h - c)
    
    # Plastic length calculation
    # --------------------------
    lp = 0.08*l + 0.022*db*(fy/units.MPa)          # [m] plastic length at the end of elements.
    rotation = curvature*lp                        # [rad] rotation along plastic hinge length.
    
    return My, curvature, rotation, lp, c