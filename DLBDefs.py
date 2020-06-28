#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 10:29:27 2019

@author: jsbaquero
"""

"""
This file contains all definitions used to analyse and design double layered 
reinforced concrete beam elements. Units kN, m, kPa, sec, m/s

Created by: Juan Sebastian Baquero M.
Date : december 2019
"""

def round5m(x):
    """
    Function that rounds values to the nearest minor multiple of 0.05. Usefull for
    rounding dimensions of concrete section elements measured in meters [m].

    Parameters
    ----------
    x : [m] dimension to round.

    Returns
    -------
    xround : [m] rounded dimension.

    """
    dec = int(x*10)
    cent = round(10*(x*10 - dec))
    if cent < 5:
        cent = 0
    elif cent >= 5:
        cent = 5
    xround = dec/10 + cent/100
    
    return xround

def round5M(x):
    """
    Function that rounds values to the nearest major multiple of 0.05. Usefull for
    rounding dimensions of concrete section elements measured in meters [m].

    Parameters
    ----------
    x : [m] dimension to round.

    Returns
    -------
    xround : [m] rounded dimension.

    """
    dec = round(x*10)
    cent = round(10*(x*10 - int(x*10)))
    if cent < 5:
        cent = 5
    else:
        cent = 0
    xround = dec/10 + cent/100
    return xround

def Rn_funcrho(fprimac,fy,rho):
    """
    Function to calculate Rn [MPa] according to mechanical characteristics
    and to rho calculated. FOr doubly reinforced concrete sections, the rho
    value corresponds to rhomax.

    Parameters
    ----------
    fprimac : [kPa] concrete compressive strength.
    fy : [kPa] yielding strength of reinforcement steel.
    rho : [m2/m2] ratio of reinforcement area against concrete section area.

    Returns
    -------
    Rn : [kPa] resistance factor for a one layer reinfroced concrete erctangular
               section.

    """
    
    # Rn calculated in MPa
    # --------------------
    Rn = fy*rho*(1 - (rho*fy)/(1.7*fprimac))
    
    return Rn

def beta1conc(fprimac):
    """
    Functino to calculate Beta1 factor accordin to ACI318-14 for concrete.

    Parameters
    ----------
    fprimac : [kPa] compresive stength of concrete.

    Returns
    -------
    beta1 : [adim] concrete strength depending factor.

    """
    if fprimac >= 17e3 and fprimac <= 28e3:
        beta1 = 0.85
    elif fprimac > 28e3 and fprimac < 55e3:
        beta1 = 0.85 - 0.05 * (fprimac - 28e3)/7e3
    elif fprimac >= 55e3:
        beta1 = 0.65
        
    return beta1

def c_2layers(fprimac,fy,As,Aprimas,beta1,b):
    """
    Function to calculate distance to neutral axis measured from farest compressed
    fiber on a rectangualr section. This function asumes that f's = fy.

    Parameters
    ----------
    fprimac : [kPa] concrete compressive strength.
    fy : [kPa] yielding strength of reinforcement steel.
    As : [m2] total reinforcement area of section.
    Aprimas : [m2] reinforcement area of the compression steel layer.
    beta1 : [adim] concrete depending factor according to ACI318-14.
    b : [m] section width.

    Returns
    -------
    c : [m] distance from farest compressed fiber to de neutral axis on a 
            rectangular section.

    """
    
    c = abs(As - Aprimas)*fy/(0.85*fprimac*beta1*b)
    
    return c

def epsilon_t(d,c):
    """
    Function to calculate tension strain on the tractioned reinforcement layer.

    Parameters
    ----------
    d : [m] distance from farest compressed concrete fiber to tractioned
            reinforcement layer. Effective depth of section.
    c : [m] distance from farest compressed concrete fiber to neutral axis.

    Returns
    -------
    e_t : [m/m] tension strain on tractioned reinforcement layer.

    """
    e_t = 0.003*(d - c)/c
    
    return e_t

def epsilon_ty(Es,fy):
    """
    Function to calculate yielding strain on a tractioned reinforcement layer.

    Parameters
    ----------
    Es : [kPa] steel modulus of elasticity.
    fy : [kPa] yielding strength of reinforcement steel.

    Returns
    -------
    epsilon_yield : [m/m] yielding strain of tractioned layer.

    """
    epsilon_yield = fy/Es
    return epsilon_yield

def epsilonprima_s(dprima,c):
    """
    Function to calculate strain on a compressed layer of reinforcement steel.

    Parameters
    ----------
    dprima : [m] distance from farest compressed concrete fiber to the first
                 compressed layer of steel reinforcement.
    c : [m] distance from farest compressed concrete fiber to neutral axis.

    Returns
    -------
    ept : [m/m] strain of compressed reinforcement steel layer.

    """
    ept = 0.003*(c - dprima)/c
    if ept > 0.0:
        ept = 0.003*(c - dprima)/c
    elif ept <= 0.0:
        ept = 0.0
    return ept

def Mn_DLayer(As1,As2,d,dprima,a,fy):
    """
    Function to calculate the nominal moment, Mn, of a rectangular doubly
    reinforced beam section.

    Parameters
    ----------
    As1 : [m2] reinforcement steel area related with the nominal moment Mn1.
               Steel + concrete forces.
    As2 : [m2] reinforcement steel area related with the nominal moment Mn2.
               Steel forces onlu.
    d : [m] effective section depth.
    dprima : [m] distance to the compressed reinforcement steel layer.
    a : [m] depth of the Whitney rectangle.
    fy : [kPa] yielding strength of reinforcement steel.

    Returns
    -------
    Mn : [kN-m] nominal moment of section.

    """
    
    Mn = As1*fy*(d - a/2) + As2*fy*(d - dprima)
    
    return Mn

def phi_AxFlx(et,ety,transtype):
    """
    Function to calculate phi reduction factor. This factor transforms nominal
    section forces to design formces.

    Parameters
    ----------
    et : [m/m] strain of tractioned reinforcement steel layer.
    ety : [m/m] yielding strain for reinforcement steel.
    transtype : STRING indicating the type of transversal reinforcement.
                It could be:
                    transtype = 'spiral'
                      or
                    transtype = 'other'

    Returns
    -------
    phi : [adim] design factor.

    """
    if transtype in ('spiral'):
        if et <= ety:
            phi = 0.75
        elif et > ety and et < 0.005:
            phi = 0.75 + 0.15*((et - ety)/(0.005 - ety))
        elif et >= 0.005:
            phi = 0.90
    elif transtype in ('other'):
        if et <= ety:
            phi = 0.65
        elif et > ety and et < 0.005:
            phi = 0.65 + 0.25*(et - ety)/(0.005 - ety)
        elif et >= 0.005:
            phi = 0.90
    return phi

def rho_gen(As,b,d):
    """
    Function to calculate the general amount of steel on a rectangular section.

    Parameters
    ----------
    As : [m2] total area of steel reinforcement on a section.
    b : [m] section width.
    d : [m] section effective depth.

    Returns
    -------
    quant : [m2/m2] ratio of reinforcement area against total section area.

    """
    quant = As/(b*d)
    return quant

def rho_max1Layer(fprimac,fy,beta1,etlim):
    """
    Function to calculate maximum amount of reinforcement area (area ratio) for
    a single reinforced rectangular concrete beam section.

    Parameters
    ----------
    fprimac : [kPa] concrete compressive strength.
    fy : [kPa] yielding strength of reinforcement steel.
    beta1 : [adim] Beta1 factor deppending on concrete strength accordin to
                   ACI318-14.
    etlim : [m/m] limit strain for ductile section consideration.

    Returns
    -------
    rho_max : [m2/m2] maximum amount of steel per section area required for a
                      singly reinforced beam section.

    """
    
    rho_max = (0.85*fprimac*beta1/fy)*(0.003/(0.003 + etlim))
    
    return rho_max
    

def procAnalysisMn(fprimac,fy,Es,As,Aprimas,b,d,dprima,transtype):
    """
    

    Parameters
    ----------
    fprimac : TYPE
        DESCRIPTION.
    fy : TYPE
        DESCRIPTION.
    Es : TYPE
        DESCRIPTION.
    As : TYPE
        DESCRIPTION.
    Aprimas : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    d : TYPE
        DESCRIPTION.
    dprima : TYPE
        DESCRIPTION.
    transtype : TYPE
        DESCRIPTION.

    Returns
    -------
    Mn : TYPE
        DESCRIPTION.
    phi : TYPE
        DESCRIPTION.
    phiMn : TYPE
        DESCRIPTION.

    """
    
    # Direct calculation of 'c' value for time advantage
    # ===================================================
    beta1 = beta1conc(fprimac)
    c1 = (-(0.003*Es*Aprimas-As*fy) + \
          ((0.003*Es*Aprimas-As*fy)**2 - \
           4*(0.85*fprimac*beta1*b)*(-0.003*Es*Aprimas*dprima))**(0.5))\
        /(2*(0.85*fprimac*beta1*b))
    c2 = (-(0.003*Es*Aprimas-As*fy) - \
          ((0.003*Es*Aprimas-As*fy)**2 - \
           4*(0.85*fprimac*beta1*b)*(-0.003*Es*Aprimas*dprima))**(0.5))\
        /(2*(0.85*fprimac*beta1*b))
    c = []
    if c1 >= 0.:
        c.append(c1)
    if c2 >= 0.:
        c.append(c2)
    c = min(c)
    # a = beta1*c
    
    # Calculation of fs and fps according to its strain state.
    # ========================================================
    (fs,fsr) = detfs(fy,Es,c,d)
    (fps,fsr) = detfs(fy,Es,c,dprima)
    
    # Nominal Moment calculation using generalized expresion
    # =======================================================
    Mn = As*fs*(c - d) + Aprimas*fps*(c - dprima) + \
        0.85*fprimac*beta1*b*c**2*(1 - beta1/2)
    
    # Phi calculation according to etreme tractinoed bar layer
    # ========================================================
    ets = epsilon_t(d,c)
    ety = epsilon_ty(Es,fy)
    phi = phi_AxFlx(ets,ety,transtype)
    
    phiMn = phi*Mn
    
    return Mn, phi, phiMn

def procDesignDLayer(fprimac,fy,Es,etlim,b,h,Mu,transtype,Mutype):
    """
    Function used to design doubly reinforced concrete beam sections. In this
    function, the input parameters are used to first determine whether a singly
    reinforced beam is sufficient or not to withstand the demand moment, Mu.
    It forces the section to become a doubly reinforced one by decreasing its
    geometry dimensions --> Not True! Modified... I will revise all later.

    Parameters
    ----------
    fprimac : [kPa] concrete compressive strength.
    fy : [kPa] yielding strength of reinforcement steel.
    Es : [kPa] steel modulus of elasticity.
    etlim : [m/m] limit strain for ductile section consideration.
    b : [m] section width.
    h : [m] section height.
    Mu : [kN-m] demand factored moment.
    transtype : STRING indicating the type of transversal reinforcement.
                It could be:
                    transtype = 'spiral'
                      or
                    transtype = 'other'
    Mutype : STRING indicating the sign of moment demand to considerations
             related to effective width considerations.
             It could be:
                    Mutype = 'positive'
                       or
                    Mutype = 'negative'

    Returns
    -------
    breq : [m] required section width.
    hreq : [m] required section height.
    dreq : [m] calculated effective hsection eight.
    dprima : [m] calculated height to compressive steel layer.
    Asreq : [m2] required tensile reinforcement steel area.
    Aprimas : [m2] required compressive reinforcement steel area.
    dsgnchck : binary integer defining whether the design process has been
               accomplished (1) or not (0).
    phi : [adim] nominal to design factor.
    Mn : [kN-m] nominal moment of section.
    phiMn : [kN-m] design moment of section.
    etcalc : [m/m] strain state of tensile reinforcement steel layer.
    rhocalc : [m2/m2] calculated amount of steel.

    """
    
    
    # Initial definitions for the function
    # ====================================
    # import DLBDefs
    hreq = h
    breq = b
    recprima = 0.0625
    d = hreq - recprima
    dprima = recprima
    phi = 0.90                     # assumed initial value of phi.
    beta1 = beta1conc(fprimac)
    ety = epsilon_ty(Es, fy)
    rhomax = rho_max1Layer(fprimac, fy, beta1, etlim)
    Rn = Rn_funcrho(fprimac, fy, rhomax)

    # Calculations for design process
    # ================================
    Mn = Mu*1.02/phi    # Mu augmented in 2% to accomplish design.
    As1 = rhomax*breq*d
    Mu1 = Rn*phi*breq*d**2
    Mn1 = Mu1/phi
    Mn2 = Mn - Mn1

    deltahreq = 0.01  # [m] decreasing value for the section heigth, if necessary.
    
    # In the case in which Mn1 >= Mn, the magnitude of Mn2 would be <= to zero,
    # thus, there would not be necesary a doble layered reinforced beam. In
    # such case, the section geometry could be optimized by decresing its section
    # until Mn2 be significant and so the second layer becomes necessary.
    
    if Mutype in ('negative'):
        
        if Mn2 <= 0.0:
            hiter = hreq
            while Mn2 <= 0.0:
                hiter = hiter - deltahreq
                hreq = round5M(hiter)
                breq = round5M(hreq/2)
                dreq = hreq - recprima
                As1 = rhomax*breq*dreq
                Mu1 = Rn*phi*breq*dreq**2
                Mn1 = Mu1/phi
                Mn2 = Mn - Mn1
        elif Mn2 > 0.0:
            hreq = round5M(hreq)
            breq = round5M(hreq/2)
            dreq = hreq - recprima
            As1 = rhomax*breq*dreq
            Mu1 = Rn*phi*breq*dreq**2
            Mn1 = Mu1/phi
            Mn2 = Mn - Mn1
            if Mn2 <= 0.0:
                hiter = hreq
                while Mn2 <= 0.0:
                    hiter = hiter - deltahreq
                    hreq = round5M(hiter)
                    breq = round5M(hreq/2)
                    dreq = hreq - recprima
                    As1 = rhomax*breq*dreq
                    Mu1 = Rn*phi*breq*dreq**2
                    Mn1 = Mu1/phi
                    Mn2 = Mn - Mn1
        # In the special case in which the section width becomes smaller than 0.25[m],
        # this dimension shall be stablished as 0.25[m] because is the minimum value
        # permited by ACI318-14. Besides, when Mn2 results less than 0.0[MN-m], the
        # calculation values will be asumed as 0.0 in order to establish a single
        # layered reinforced concrete section.
        if breq <= 0.25:
            breq = 0.25
            hreq = round5M(breq*1.5)
            dreq = hreq - recprima
            As1 = rhomax*breq*dreq
            Mu1 = Rn*phi*breq*dreq**2
            Mn1 = Mu1/phi
            Mn2 = Mn - Mn1
            if Mn2 <= 0.0:
                Mn2 = 0.0
                
        # Let's first suppose that fprimas == fy.
        # =======================================
        c = (As1*fy)/(0.85*fprimac*breq*beta1)
        eprimas = epsilonprima_s(dprima, c)
        
        if eprimas >= ety:
            Aprimas = Mn2/(fy*(dreq - dprima))
            As2 = Aprimas
            Asreq = As1 + As2
            etcalc = (dreq - c)/c*0.003
            phi = phi_AxFlx(etcalc, ety, transtype)
            Mn = (Asreq - As2)*fy*(dreq - beta1*c/2) + As2*fy*(dreq - dprima)
            phiMn = phi*Mn
        elif eprimas < ety:
            fprimas = eprimas*Es
            Aprimas = Mn2/(fprimas*(dreq - dprima))
            As2 = Aprimas*fprimas/fy
            Asreq = As1 + As2
            etcalc = (dreq - c)/c*0.003
            phi = phi_AxFlx(etcalc, ety, transtype)
            Mn = As1*fy*(dreq - beta1*c/2) + Aprimas*fprimas*(dreq - dprima)
            phiMn = phi*Mn
        
    elif Mutype in ('positive'):
        # There will be assumed that always produces a single layer reinforced
        # beam due to the dimensions selected from negative moment demand design
        # and the only proceure to execute will be to reduce As until Mn_design
        # reaches Mn_required.
        
        if Mn2 <= 0.0:
            Mn2 = 0.0
            Aprimas = 0.0
            dreq = d
            c = (As1*fy)/(0.85*fprimac*breq*beta1)
            eprimas = epsilonprima_s(dprima, c)
            Asreq = As1
            etcalc = (dreq - c)/c*0.003
            phi = phi_AxFlx(etcalc, ety, transtype)
            Mn = (Asreq)*fy*(dreq - beta1*c/2)
            phiMn = phi*Mn
            while phiMn > 1.02*Mu:
                Asreq = Asreq - 0.00005
                c = (Asreq*fy)/(0.85*fprimac*breq*beta1)
                eprimas = epsilonprima_s(dprima, c)
                etcalc = (dreq - c)/c*0.003
                phi = phi_AxFlx(etcalc, ety, transtype)
                Mn = (Asreq)*fy*(dreq - beta1*c/2)
                phiMn = phi*Mn
        
    rhocalc = (Asreq + Aprimas)/(breq*hreq)
    
    # Checking of design accomplishment.
    # ==================================
    if phiMn >= Mu:
        dsgnchck = 1
    elif phiMn < Mu:
        dsgnchck = 0
            
    return breq, hreq, dreq, dprima, Asreq, Aprimas, dsgnchck, phi, Mn, phiMn, etcalc, rhocalc


def procDesignDLayerMod(fprimac,fy,Es,etlim,b,h,Mu,transtype,Mutype):
    """
    Function used to design doubly reinforced concrete beam sections. In this
    function, the input parameters are used to first determine whether a singly
    reinforced beam is sufficient or not to withstand the demand moment, Mu.
    It forces the section to become a doubly reinforced one by decreasing its
    geometry dimensions --> Not True! Modified... I will revise all later.

    Parameters
    ----------
    fprimac : [kPa] concrete compressive strength.
    fy : [kPa] yielding strength of reinforcement steel.
    Es : [kPa] steel modulus of elasticity.
    etlim : [m/m] limit strain for ductile section consideration.
    b : [m][LIST] section width containing both bw and bweff calculated according
                  to location of beam element into the diaphragm.
                  b = [bw, bweff]
    h : [m] section height.
    Mu : [kN-m] demand factored moment.
    transtype : STRING indicating the type of transversal reinforcement.
                It could be:
                    transtype = 'spiral'
                      or
                    transtype = 'other'
    Mutype : STRING indicating the sign of moment demand to considerations
             related to effective width considerations.
             It could be:
                    Mutype = 'positive'
                       or
                    Mutype = 'negative'

    Returns
    -------
    breq : [m] required section width.
    hreq : [m] required section height.
    dreq : [m] calculated effective hsection eight.
    dprima : [m] calculated height to compressive steel layer.
    Asreq : [m2] required tensile reinforcement steel area.
    Aprimas : [m2] required compressive reinforcement steel area.
    dsgnchck : binary integer defining whether the design process has been
               accomplished (1) or not (0).
    phi : [adim] nominal to design factor.
    Mn : [kN-m] nominal moment of section.
    phiMn : [kN-m] design moment of section.
    etcalc : [m/m] strain state of tensile reinforcement steel layer.

    """
    # Initial definitions for the function
    # ====================================
    # import DLBDefs
    hreq = h
    breq = b[1]
    recprima = 0.0625
    d = hreq - recprima
    dprima = recprima
    phi = 0.90   # valor inicial de phi asumido.
    beta1 = beta1conc(fprimac)
    ety = epsilon_ty(Es, fy)
    rhomax = rho_max1Layer(fprimac, fy, beta1, etlim)
    Rn = Rn_funcrho(fprimac,fy,rhomax)

    # Calculations for design process
    # ================================
    Mn = Mu*1.02/phi    # Mu augmented in 2% to accomplish design.
    As1 = rhomax*breq*d
    Mu1 = Rn*phi*breq*d**2
    Mn1 = Mu1/phi
    Mn2 = Mn - Mn1

    dreq = d
    # In the case in which Mn1 >= Mn, the magnitude of Mn2 would be <= to zero,
    # thus, there would not be necesary a doble layered reinforced beam. In
    # such case, the section geometry could be optimized by decresing its section
    # until Mn2 be significant and so the second layer becomes necessary.
    #  - (11/feb/2020) But for seismic design what controls the dimensions on
    # elements would be the drift ratio magnitud and the stiffness required 
    # for the structure to behave properly.
    # So, in this modified function I use the dimensions and let the beam be a
    # one layer beam if it is so.
    
    if Mutype in ('negative','positive'):
        
        if Mn2 > 0.0:
            
            c = (As1*fy)/(0.85*fprimac*breq*beta1)
            eprimas = epsilonprima_s(dprima, c)
            if eprimas >= ety:
                Aprimas = Mn2/(fy*(dreq - dprima))
                As2 = Aprimas
                Asreq = As1 + As2
                etcalc = (dreq - c)/c*0.003
                phi = phi_AxFlx(etcalc,ety,transtype)
                Mn = (Asreq - As2)*fy*(dreq - beta1*c/2) + As2*fy*(dreq - dprima)
                phiMn = phi*Mn
            elif eprimas < ety:
                fprimas = eprimas*Es
                Aprimas = Mn2/(fprimas*(dreq - dprima))
                As2 = Aprimas*fprimas/fy
                Asreq = As1 + As2
                etcalc = (dreq - c)/c*0.003
                phi = phi_AxFlx(etcalc, ety, transtype)
                Mn = As1*fy*(dreq - beta1*c/2) + Aprimas*fprimas*(dreq - dprima)
                phiMn = phi*Mn
                
        elif Mn2 <= 0.0:
            
            Mn2 = 0.0
            Aprimas = 0.0
            dreq = d
            c = (As1*fy)/(0.85*fprimac*breq*beta1)
            eprimas = epsilonprima_s(dprima, c)
            Asreq = As1
            etcalc = (dreq - c)/c*0.003
            phi = phi_AxFlx(etcalc, ety, transtype)
            Mn = (Asreq)*fy*(dreq - beta1*c/2)
            phiMn = phi*Mn
            while phiMn > 1.02*Mu:
                Asreq = Asreq - 0.00005
                c = (Asreq*fy)/(0.85*fprimac*breq*beta1)
                eprimas = epsilonprima_s(dprima, c)
                etcalc = (dreq - c)/c*0.003
                phi = phi_AxFlx(etcalc, ety, transtype)
                Mn = (Asreq)*fy*(dreq - beta1*c/2)
                phiMn = phi*Mn
    
    import unitskN_m_C as units

    Asreq = max(0.25*(fprimac*(1/units.MPa))**(0.5)*min(2*b[0],b[1])*dreq/(fy*(1/units.MPa)),\
                1.4/(fy*(1/units.MPa))*min(2*b[0],b[1])*dreq,Asreq)
    rhocalc = (Asreq + Aprimas)/(breq*hreq)
    
    # Checking of design accomplishment.
    # ==================================
    if phiMn >= Mu:
        dsgnchck = 1
    elif phiMn < Mu:
        dsgnchck = 0
            
    return dreq, dprima, Asreq, Aprimas, dsgnchck, phi, Mn, phiMn, etcalc, rhocalc

def SpecialBeamDesign(fprimac,fy,Es,etlim,b,h,Muipos,Muineg,Mumpos,Mumneg,Mujpos,\
                      Mujneg,transtype,lper,lpar,beamloc,hslab=0.25):
    """
    

    Parameters
    ----------
    fprimac : [kPa] TYPE
        DESCRIPTION.
    fy : [kPa] TYPE
        DESCRIPTION.
    Es : [kPa] TYPE
        DESCRIPTION.
    etlim : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    h : TYPE
        DESCRIPTION.
    Muipos : TYPE
        DESCRIPTION.
    Muineg : TYPE
        DESCRIPTION.
    Mumpos : TYPE
        DESCRIPTION.
    Mumneg : TYPE
        DESCRIPTION.
    Mujpos : TYPE
        DESCRIPTION.
    Mujneg : TYPE
        DESCRIPTION.
    transtype : TYPE
        DESCRIPTION.
    lper : TYPE
        DESCRIPTION.
    lpar : TYPE
        DESCRIPTION.
    beamloc : TYPE
        DESCRIPTION.
    hslab : TYPE, optional
        DESCRIPTION. The default is 0.25.

    Returns
    -------
    Asi : TYPE
        DESCRIPTION.
    Apsi : TYPE
        DESCRIPTION.
    Asm : TYPE
        DESCRIPTION.
    Apsm : TYPE
        DESCRIPTION.
    Asj : TYPE
        DESCRIPTION.
    Apsj : TYPE
        DESCRIPTION.

    """
    # DESIGN FOR THE SIX MOMENT VALUES
    #   ==================================
    """
        The design for the six moment values is excecuted separately in order to 
        determine the definitive dimensions of the beam first, and the reinforcement
        of the center section and both i and j ends depending on its Seismic Design
        Category (SDC).
    """
    # Design for negative moments
    # ============================
    # Node i
    # -------
    (dimin,dpriimin,Asimin,Aprimasimin,dsgnchckimin1,\
     phiimin,Mnimin,phiMnimin,etcalcimin,rhocalcimin) = \
        procDesignDLayerMod(fprimac,fy,Es,etlim,[b,b],h,Muineg,transtype,'negative')
    # Node j
    # -------
    (djmin,dprijmin,Asjmin,Aprimasjmin,dsgnchckjmin1,\
     phijmin,Mnjmin,phiMnjmin,etcalcjmin,rhocalcjmin) = \
        procDesignDLayerMod(fprimac,fy,Es,etlim,[b,b],h,Mujneg,transtype,'negative')
    # Node m
    # -------
    (dmmin,dprimmin,Asmmin,Aprimasmmin,dsgnchckmmin1,\
     phimmin,Mnmmin,phiMnmmin,etcalcmmin,rhocalcmmin) = \
        procDesignDLayerMod(fprimac,fy,Es,etlim,[b,b],h,Mumneg,transtype,'negative')
    
    """
        Mostly, the negative moment on a beam end section controls the design and
        there the dimensions determined for a beam. Assuming this fact, from negative
        moment demand designs, the general dimensions are establidhed and used to
        execute design for positive moments using effective width including slab
        contribution.
    """
    
    bbeam = b    # redefining the defitive width.
    hbeam = h    # redefining the definitive height.
    
    # DETERMINATION OF EFFECTIVE WIDTH FOR POSITIVE MOMENT ON BEAMS
    # =============================================================
    """
        According to ACI318, when on positive moments actions over beams, there
        is a need to determine the effective width to calculate de nominal Moment,
        Mn of the section considering an area effective area of slab. The length used
        as effective wxtra width depends on the location of the beam; being interior
        or border beam.
        
    """
    sw = lper - bbeam       # [m] free distance between webs of adjacent beams.
    ln = lpar - 0.40    # [m] free spam length parallel to beam axis assuming 
                        #     minimum pressumed column size.
    if beamloc in ('edge'):
        bf1 = 6*hslab
        bf2 = sw/2
        bf3 = ln/12
        bf = min(bf1,bf2,bf3)
        bweff = bbeam + bf
    elif beamloc in ('interior'):
        bf1 = 8*hslab
        bf2 = sw/2
        bf3 = ln/8
        bf = min(bf1,bf2,bf3)
        bweff = bbeam + 2*bf
    
    # Design for positive moments
    # ============================
    # Node i
    # ------
    (dimax,dpriimax,Asimax,Aprimasimax,dsgnchckimax1,\
     phiimax,Mnimax,phiMnimax,etcalcimax,rhocalcimax) = \
        procDesignDLayerMod(fprimac,fy,Es,etlim,[b,bweff],hbeam,Muipos,transtype,'positive')
    # Node j
    # ------
    (djmax,dprijmax,Asjmax,Aprimasjmax,dsgnchckjmax1,\
     phijmax,Mnjmax,phiMnjmax,etcalcjmax,rhocalcjmax) = \
        procDesignDLayerMod(fprimac,fy,Es,etlim,[b,bweff],hbeam,Mujpos,transtype,'positive')
    # Node m
    # -------
    (dmmax,dprimmax,Asmmax,Aprimasmmax,dsgnchckmmax1,\
     phimmax,Mnmmax,phiMnmmax,etcalcmmax,rhocalcmmax) = \
        procDesignDLayerMod(fprimac,fy,Es,etlim,[b,bweff],hbeam,Mumpos,transtype,'positive')
    
    # DEFINING DEFITIVE REINFORCEMENT OF BEAM
    #   =========================================
    """
        The maximum amount of steel reinforcement for both top and bottom layers
        are defined according to the area needed by moment demand in each of the
        three sections for negative and positive moments.
    """
    Asi = max(Asimax, Aprimasimin)
    Apsi = max(Aprimasimax, Asimin)
    # Asm = max(Asmmax, Aprimasmmin)
    # Apsm = max(Aprimasmmax, Asmmin)
    Asj = max(Asjmax,Aprimasjmin)
    Apsj = max(Aprimasjmax, Asjmin)
    
    return Asi, Apsi, Asj, Apsj
    # return Asi, Apsi, Asm, Apsm, Asj, Apsj

def BeamEndsMn(fprimac,fy,Es,b,h,Asi,Apsi,Asj,Apsj,transtype):
    """
    

    Parameters
    ----------
    fprimac : TYPE
        DESCRIPTION.
    fy : TYPE
        DESCRIPTION.
    Es : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    h : TYPE
        DESCRIPTION.
    Asi : TYPE
        DESCRIPTION.
    Apsi : TYPE
        DESCRIPTION.
    Asj : TYPE
        DESCRIPTION.
    Apsj : TYPE
        DESCRIPTION.
    transtype : TYPE
        DESCRIPTION.

    Returns
    -------
    Mnipos : TYPE
        DESCRIPTION.
    Mnineg : TYPE
        DESCRIPTION.
    Mnjpos : TYPE
        DESCRIPTION.
    Mnjneg : TYPE
        DESCRIPTION.

    """
    # Previos calculations an data
    # -----------------------------
    dprima = 0.0625
    d = h - dprima  
    # Nominal Moment for node i
    # =========================
    (Mnipos,_,_) = procAnalysisMn(fprimac,fy,Es,Asi,Apsi,b,d,dprima,transtype)
    (Mnineg,_,_) = procAnalysisMn(fprimac,fy,Es,Apsi,Asi,b,d,dprima,transtype)
    # Nominal Moment for node j
    # =========================
    (Mnjpos,_,_) = procAnalysisMn(fprimac,fy,Es,Asj,Apsj,b,d,dprima,transtype)
    (Mnjneg,_,_) = procAnalysisMn(fprimac,fy,Es,Apsj,Asj,b,d,dprima,transtype)
    
    return Mnipos, Mnineg, Mnjpos, Mnjneg

def detfs(fy,Es,c,d):
    """
    

    Parameters
    ----------
    fy : [kPa] TYPE
        DESCRIPTION.
    Es : [kPa] TYPE
        DESCRIPTION.
    c : TYPE
        DESCRIPTION.
    d : TYPE
        DESCRIPTION.

    Returns
    -------
    fsy : [kPa] TYPE
        DESCRIPTION.
    fsr : [kPa] TYPE
        DESCRIPTION.

    """
    # calculates the state of stress of a steel reinforcement layer. Used
    # for column analysis of nominal moments.
    import numpy as np
    fsr = 0.003*Es*(c-d)/c;    # [kN/m2] real stress on the steel layer.
    signo = np.sign(fsr)
    if abs(fsr) > fy:
        fsy = signo*fy
    else:
        fsy = fsr
    
    return fsy, fsr

def PnCol(fprimac,fy,Es,b,d,As,nlayers,c):
    """
    

    Parameters
    ----------
    fprimac : [kPa] TYPE
        DESCRIPTION.
    fy : [kPa] TYPE
        DESCRIPTION.
    Es : [kPa] TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    d : TYPE
        DESCRIPTION.
    As : TYPE
        DESCRIPTION.
    nlayers : TYPE
        DESCRIPTION.
    c : TYPE
        DESCRIPTION.

    Returns
    -------
    Pn : [kN] TYPE
        DESCRIPTION.

    """
    # Determines nominal axial strength of a column, given a value of c.
    beta1 = beta1conc(fprimac)
    Pc = 0.85*fprimac*beta1*c*b    # [kN] contribution of concrete
    Ps = 0.0                       # [kN] contribution of steel layers.
    for i in range(nlayers):
        (fsy,fsr) = detfs(fy,Es,c,d[i])
        Ps += As[i]*fsy
    Pn = (Pc + Ps)
    
    return Pn

def MnCol(fprimac,fy,Es,b,h,d,As,nlayers,c):
    """
    

    Parameters
    ----------
    fprimac : [kPa] TYPE
        DESCRIPTION.
    fy : [kPa] TYPE
        DESCRIPTION.
    Es : [kPa] TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    h : TYPE
        DESCRIPTION.
    d : TYPE
        DESCRIPTION.
    As : TYPE
        DESCRIPTION.
    nlayers : TYPE
        DESCRIPTION.
    c : TYPE
        DESCRIPTION.

    Returns
    -------
    Mn : [kN-m] TYPE
        DESCRIPTION.

    """
    # Calulates the nominal moment of a column given a value of c and column
    # geometry.
    beta1 = beta1conc(fprimac)
    Mc = 0.85*fprimac*beta1*c*b*0.5*(h-c*beta1)    # [kN-m] Nominal moment given just by concrete.
    Ms = 0.0                                       # [kN-m] nominal comtribution of steel.
    for i in range(nlayers):
        (fsy,fsr) = detfs(fy,Es,c,d[i])
        Ms += As[i]*fsy*(0.5*h - d[i])
    Mn = Mc + Ms
    
    return Mn

def ColAnalysisMn(fprimac,fy,Es,ety,Pu,b,h,d,As,nlayers,transtype):
    """
    

    Parameters
    ----------
    fprimac : [kPa] TYPE
        DESCRIPTION.
    fy : [kPa] TYPE
        DESCRIPTION.
    Es : [kPa] TYPE
        DESCRIPTION.
    ety : TYPE
        DESCRIPTION.
    Pu : [kN] TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    h : TYPE
        DESCRIPTION.
    d : TYPE
        DESCRIPTION.
    As : TYPE
        DESCRIPTION.
    nlayers : TYPE
        DESCRIPTION.
    transtype : TYPE
        DESCRIPTION.

    Returns
    -------
    Mn : [kN-m] TYPE
        DESCRIPTION.
    Pn : [kN] TYPE
        DESCRIPTION.
    phi : TYPE
        DESCRIPTION.

    """
    # Calculates a Mn according to an applied Pu on the column.
    
    # Previous calculations
    # ----------------------
    beta1 = beta1conc(fprimac)
    cmin = 0.001     # [m] minumum value of c.
    cmax = h/beta1   # [m] maximum value of c.
    
    # Iterative process to find Pn-Mn that matches Pu.
    # ------------------------------------------------
    citer = cmin
    while citer < cmax:
        Mn = MnCol(fprimac,fy,Es,b,h,d,As,nlayers,citer)   # [kN-m]
        Pn = PnCol(fprimac,fy,Es,b,d,As,nlayers,citer)     # [kN]
        et = epsilon_t(d[-1],citer)
        phi = phi_AxFlx(et,ety,transtype)
        phiPn = phi*Pn
        if phiPn >= Pu:
            break
        citer += 0.001
    
    return Mn, Pn, phi