# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 13:16:51 2020

@author: sebas
"""

"""
    Procedures to calculate seismic forces according to ASCE7-16.
"""

def detFa(Ss,Vs30):
    """
    Function to determine Fa factor according to Site Class in Table 11.4-1
    of ASCE7-16. Useful to determine Spectral Response acceleration parameters
    adjusted to soil characteristics.

    Parameters
    ----------
    Ss : [g] MCE_R ground motion parameter for 0.20s Spectral response Acceleration. 
    SiteClass : [STRING] Site Class according to ASCE7-16. Could be:
                        'A' , 'B' , 'C' , 'D' , 'E' or 'F'

    Returns
    -------
    Fa : [adim] Site coeficient.

    """
    
    import numpy as np
    
    if Ss <= 0.25:
        Fa = np.exp(-0.727*np.log(Vs30/760) - 0.2298*(np.exp(-0.00638*(min(Vs30,760) - 360)) - \
                                                      np.exp(-0.00638*400))*np.log(((0.25/2.30) + 0.10)/0.10))
    elif Ss >= 1.50:
        Fa = np.exp(-0.727*np.log(Vs30/760) - 0.2298*(np.exp(-0.00638*(min(Vs30,760) - 360)) - \
                                                      np.exp(-0.00638*400))*np.log(((1.50/2.30) + 0.10)/0.10))
    elif Ss > 0.25 and Ss < 1.50:
        Fa = np.exp(-0.727*np.log(Vs30/760) - 0.2298*(np.exp(-0.00638*(min(Vs30,760) - 360)) - \
                                            np.exp(-0.00638*400))*np.log(((Ss/2.30) + 0.10)/0.10))
    
    sigmaFa = 0.67   # standar deviation of Fa values.
    
    if Vs30 < 180 and Ss < 1.0:
        Fa = np.exp(np.log(Fa) + 0.5*sigmaFa)
    elif Vs30 < 180 and Ss >= 1.0 and Ss < 1.50:
        Fa = np.exp(-0.727*np.log(489/760) - 0.2298*(np.exp(-0.00638*(min(489,760) - 360)) - \
                                                     np.exp(-0.00638*400))*np.log(((Ss/2.30) + 0.10)/0.10))
            # that of SiteClass C with Ss >=1.0 and < 1.50
    elif Vs30 < 180 and Ss >= 1.50:
        Fa = np.exp(-0.727*np.log(489/760) - 0.2298*(np.exp(-0.00638*(min(489,760) - 360)) - \
                                                     np.exp(-0.00638*400))*np.log(((1.50/2.30) + 0.10)/0.10))
    
        
    return Fa

def detFv(S1,Vs30):
    """
    Function to determine Fv factor according to Site Class in Table 11.4-2
    of ASCE7-16. Useful to determine Spectral Response acceleration parameters
    adjusted to soil characteristics.

    Parameters
    ----------
    S1 : [g] MCE_R ground motion parameter for 1.0s Spectral response Acceleration.
    SiteClass : [STRING] Site Class according to ASCE7-16. Could be:
                        'A' , 'B' , 'C' , 'D' , 'E' or 'F'

    Returns
    -------
    Fv : [adim] Site coeficient.

    """
    import numpy as np
    
    if S1 <= 0.10:
        Fv = np.exp(-1.03*np.log(Vs30/760) - 0.118*(np.exp(-0.00756*(min(Vs30,760) - 360)) - \
                                                    np.exp(-0.00756*400)*np.log(((0.10/0.7) + 0.1)/0.1)))
    elif S1 >= 0.60:
        Fv = np.exp(-1.03*np.log(Vs30/760) - 0.118*(np.exp(-0.00756*(min(Vs30,760) - 360)) - \
                                                    np.exp(-0.00756*400)*np.log(((0.60/0.7) + 0.1)/0.1)))
    elif S1 > 0.1 and S1 < 0.6:
        Fv = np.exp(-1.03*np.log(Vs30/760) - 0.118*(np.exp(-0.00756*(min(Vs30,760) - 360)) - \
                                                    np.exp(-0.00756*400)*np.log(((S1/0.7) + 0.1)/0.1)))
    
    sigmaFv = 0.58   # standar deviation of Fa values.
    
    if Vs30 < 180:
        Fv = np.exp(np.log(Fv) + 0.5*sigmaFv)
        
    return Fv

def detSMS(Ss,Vs30):
    """
    Function to determine Spectral response acceleration parameters for short
    periods, SMS. It is calculated acording to equation 11.4-1 on ASCE7-16.

    Parameters
    ----------
    Ss : [g] MCE_R ground motion parameter for 0.20s Spectral response Acceleration. 
    SiteClass : [STRING] Site Class according to ASCE7-10. Could be:
                        'A' , 'B' , 'C' , 'D' , 'E' or 'F'

    Returns
    -------
    SMS : [g] Spectral response Acceleration parameter adjusted to site class.

    """
    Fa = detFa(Ss,Vs30)
    SMS = Fa*Ss
    
    return SMS

def detSM1(S1,Vs30):
    '''
    Function to determine Spectral response acceleration parameters for 1.0s
    period, SM1. It is calculated acording to equation 11.4-2 on ASCE7-10.

    Parameters
    ----------
    S1 : [g] MCE_R ground motion parameter for 1.0s Spectral response Acceleration.
    SiteClass : [STRING] Site Class according to ASCE7-10. Could be:
                        'A' , 'B' , 'C' , 'D' , 'E' or 'F'

    Returns
    -------
    SM1 : [g] Spectral response Acceleration parameter adjusted to site class.

    '''
    Fv = detFv(S1,Vs30)
    SM1 = Fv*S1
    
    return SM1

def detSDS(Ss,Vs30):
    """
    Function to determine DEIGN Spectral response acceleration parameters for 
    short periods, SMS. It is calculated acording to equation 11.4-3 on ASCE7-10.

    Parameters
    ----------
    Ss : [g] MCE_R ground motion parameter for 0.20s Spectral response Acceleration. 
    SiteClass : [STRING] Site Class according to ASCE7-10. Could be:
                        'A' , 'B' , 'C' , 'D' , 'E' or 'F'

    Returns
    -------
    SDS : [g] DESIGN Spectral response Acceleration parameter.

    """
    Fa = detFa(Ss,Vs30)
    SMS = Fa*Ss
    SDS = 2/3*SMS
    
    return SDS
    
def detSD1(S1,Vs30):
    """
    Function to determine DEIGN Spectral response acceleration parameters for 
    1.0s period, SM1. It is calculated acording to equation 11.4-4 on ASCE7-10.

    Parameters
    ----------
    S1 : [g] MCE_R ground motion parameter for 1.0s Spectral response Acceleration.
    SiteClass : [STRING] Site Class according to ASCE7-10. Could be:
                        'A' , 'B' , 'C' , 'D' , 'E' or 'F'

    Returns
    -------
    SD1 : [g] DESIGN Spectral response Acceleration parameter.

    """
    Fv = detFv(S1,Vs30)
    SM1 = Fv*S1
    SD1 = 2/3*SM1
    
    return SD1
    
def detSiteClass(Vs30):
    # SiteClass definition according to Chapter 20 of ASCE7-16.
    if Vs30 > 1524:
        SiteClass = 'A'
    elif Vs30 > 762 and Vs30 <= 1524:
        SiteClass = 'B'
    elif Vs30 > 366 and Vs30 <= 762:
        SiteClass = 'C'
    elif Vs30 > 180 and Vs30 <= 366:
        SiteClass = 'D'
    elif Vs30 <= 180:
        SiteClass = 'E'
        
    return SiteClass

def Sa(Ss,S1,Vs30,T,TL=8.0):
    """
    Function to calculate Design Response Spectrum Acceleration, Sa according
    to procedures described in Section 11.4.5 of ASCE7-10.

    Parameters
    ----------
    Ss : [g] MCE_R ground motion parameter for 0.20s Spectral response Acceleration. 
    S1 : [g] MCE_R ground motion parameter for 1.0s Spectral response Acceleration.
    SiteClass : [STRING] Site Class according to ASCE7-10. Could be:
                        'A' , 'B' , 'C' , 'D' , 'E' or 'F'
    T : [s] Calculated or approximated period of vibration of the structure.
    TL : [s] Long period value according to Figures 22-12 through 22-16. 
             The default is 8.0.

    Returns
    -------
    Sa : [g] Design Response Spectrum Acceleration.
    Sa_MCER : [g] Response Spectrum Acceleration for Maximum Considered Earthquake.
    To : [s] Limit period for the begining of constant acceleration zone.
    Ts : [s] Limit period for the end of constant acceleration zone.

    """
    # Previous calculations concerning to limiting periods.
    # ======================================================
    SD1 = detSD1(S1,Vs30)
    SDS = detSDS(Ss,Vs30)
    To = 0.20*SD1/SDS
    Ts = SD1/SDS
    # Sa calculation
    # --------------
    if T < To:
        Sa = SDS*(0.40 + 0.60*T/To)
    elif T >= To and T <= Ts:
        Sa = SDS
    elif T > Ts and T <= TL:
        Sa = SD1/T
    elif T > TL:
        Sa = SD1*TL/T**2
    Sa_MCER = Sa*1.50
    
    return Sa, Sa_MCER, To, Ts
    
def detSDCat(RiskCat,Ss,S1,Vs30):
    """
    Function to determine the Seismic Design Category (SDC)  according to
    Sectionn 11.6 of ASCE7-16. It makes use of Tables 11.6-1 and 11.6-2 along
    with Risk Category data from Table 1.5-1.

    Parameters
    ----------
    RiskCat : [STRING] Risk Category of buildings according to Table 1.5-1
                       of ASCE7-10. It could be:
                           'I' , 'II', 'III' or 'IV'
    Ss : [g] MCE_R ground motion parameter for 0.20s Spectral response Acceleration.
    S1 : [g] MCE_R ground motion parameter for 1.0s Spectral response Acceleration.
    SiteClass :[STRING] Site Class according to ASCE7-10. Could be:
                        'A' , 'B' , 'C' , 'D' , 'E' or 'F'

    Returns
    -------
    SDC : [STRING] Seismic Design Category used to determine the type of 
                   structural systems in seismic design.

    """
    # Previous calculations related to SDS and SD1
    # --------------------------------------------
    SDS = detSDS(Ss,Vs30)
    SD1 = detSD1(S1,Vs30)
    Categories = ['None','A','B','C','D','E','F']
    # Previous consierations
    # ----------------------
    if RiskCat in ('I','II','III') and S1 >= 0.75:
        SDC1 = 5
    elif RiskCat in ('IV') and S1 >= 0.75:
        SDC1 = 6
    elif RiskCat in ('I','II','III','IV') and S1 < 0.75:
        SDC1 = 0
        
    # According to Table 11.6-1 
    # --------------------------
    if RiskCat in ('I','II','III'):
        if SDS < 0.167:
            SDC2 = 1
        elif SDS >= 0.167 and SDS < 0.33:
            SDC2 = 2
        elif SDS >= 0.33 and SDS < 0.50:
            SDC2 = 3
        elif SDS >= 0.50:
            SDC2 = 4
    elif RiskCat in ('IV'):
        if SDS < 0.167:
            SDC2 = 1
        elif SDS >= 0.167 and SDS < 0.33:
            SDC2 = 3
        elif SDS >= 0.33 and SDS < 0.50:
            SDC2 = 4
        elif SDS >= 0.50:
            SDC2 = 4
            
    # According to Table 11.6-2 
    # --------------------------
    if RiskCat in ('I','II','III'):
        if SD1 < 0.067:
            SDC3 = 1
        elif SD1 >= 0.067 and SD1 < 0.133:
            SDC3 = 2
        elif SD1 >= 0.133 and SD1 < 0.20:
            SDC3 = 3
        elif SD1 >= 0.20:
            SDC3 = 4
    elif RiskCat in ('IV'):
        if SD1 < 0.067:
            SDC3 = 1
        elif SD1 >= 0.067 and SD1 < 0.133:
            SDC3 = 3
        elif SD1 >= 0.133 and SD1 < 0.20:
            SDC3 = 4
        elif SD1 >= 0.20:
            SDC3 = 4
    
    # Selecting the appropriate Seismic Design Category (SDC)
    # -------------------------------------------------------
    SDC = Categories[max(SDC1,SDC2,SDC3)]
    
    return SDC
    
def RedFact(SDC):
    """
    Redundancy factor, rho, determined according to Section 12.3.4 of ASCE7-16.

    Parameters
    ----------
    SDC : [STRING] Seismic Design Category used to determine the type of 
                   structural systems in seismic design.
                   It could be:
                       'A','B','C','D','E' or 'F'

    Returns
    -------
    redunrho : [adim] redundancy factor according to SDC.

    """
    if SDC in ('A','B','C'):
        redunrho = 1.0
    elif SDC in ('D','E','F'):
        redunrho = 1.30
        
    return redunrho

def detCu(S1,Vs30):
    """
    Function that determines the COefficient for upper limit on calculated
    period, Cu. It is obtained interpolating data from Table 12.8-1 on ASCE7-16.

    Parameters
    ----------
    S1 : [g] MCE_R ground motion parameter for 1.0s Spectral response Acceleration.
    SiteClass : [STRING] Site Class according to ASCE7-16. Could be:
                        'A' , 'B' , 'C' , 'D' , 'E' or 'F'

    Returns
    -------
    Cu : [adim] Coefficient for upper limit on calculated period.
    SD1 : [g] DESIGN Spectral response Acceleration parameter.

    """
    SD1 = detSD1(S1,Vs30)
    if SD1 >= 0.30:
        Cu = 1.40
    elif SD1 >= 0.20 and SD1 < 0.30:
        Cu = 1.50 + (SD1 - 0.20)*(1.40 - 1.50)/(0.30 - 0.20)
    elif SD1 >= 0.15 and SD1 < 0.20:
        Cu = 1.60 + (SD1 - 0.15)*(1.50 - 1.60)/(0.20 - 0.15)
    elif SD1 >= 0.10 and SD1 < 0.15:
        Cu = 1.70 + (SD1 - 0.10)*(1.60 - 1.70)/(0.15 - 0.10)
    elif SD1 < 0.10:
        Cu = 1.70
        
    return  Cu, SD1

def detTa(hn,StrType):
    """
    Function to calculate the approximate period of vibration of the building
    according to Section 12.8.2 of ASCE7-16.

    Parameters
    ----------
    hn : [m] Building heigth measured from the base to the roof.
    StrType : [INTEGER] it is an integer defining the structural system type
              according to Table 12.8-2 on ASCE7-10. The values are:
                  [0] For Steel moment resisting frames.
                  [1] For concrete moment resisting frames.
                  [2] Steel excentrically braced  frames.
                  [3] Steel buckling-restrained braced frames.
                  [4] All other structural systems.

    Returns
    -------
    Ta : [s] approximated structural fundamental period of vibration.

    """
    # Data for structural type
    # ------------------------
    Ctst = [0.0724,0.0466,0.0731,0.0731,0.0488]
    xst = [0.8,0.9,0.75,0.75,0.75]
    # Calculation of approximate period.
    # ----------------------------------
    Ct = Ctst[StrType]
    x = xst[StrType]
    Ta = Ct*hn**x
    
    return Ta
    
def detCs(Ss,S1,Vs30,T,R,Ie,TL=8.0):
    """
    Function to determine the Seismic Response Coefficient according to Section
    12.8.1.1 on ASCE7-16

    Parameters
    ----------
    Ss : [g] MCE_R ground motion parameter for 0.20s Spectral response Acceleration.
    S1 : [g] MCE_R ground motion parameter for 1.0s Spectral response Acceleration.
    SiteClass : [STRING] Site Class according to ASCE7-10. Could be:
                        'A' , 'B' , 'C' , 'D' , 'E' or 'F'
    
    T : [s] Calculated or approximated period of vibration of the structure.
    R : [adim] Response modification coefficient according to Table 12.2-1.
    Ie : [adim] Importance factor according to Table 1.5-2. It could be:
                      For Risk Category       Seismic Importance Factor, Ie
                      -----------------        ---------------------------
                             I                           1.00
                             II                          1.00
                             III                         1.25
                             IV                          1.50
    TL : [s] Long period value according to Figures 22-12 through 22-16. 
             The default is 8.0.

    Returns
    -------
    Cs : [adim] Seismic Response Coefficient.

    """
    SDS = detSDS(Ss,Vs30)
    SD1 = detSD1(S1,Vs30)
    SiteClass = detSiteClass(Vs30)
    Ta = T
    (_,_,_,Ts) = Sa(Ss,S1,Vs30,T)

    Cscalc = SDS/(R/Ie)                  # According to 12.8-2    
    if SiteClass in ('D') and S1 >= 0.20:   # According to 11.4.8
        if Ta <= 1.50*Ts:
            Csmax = SDS/(R/Ie)
        elif Ta > 1.50*Ts and Ta <= TL:
            Csmax = 1.5*SD1/(Ta*R/Ie)            # According to 12.8-3
        elif Ta > TL:
            Csmax = 1.5*SD1*TL/((Ta**2)*R/Ie)    # According to 12.8-4
    else:
        Cscalc = SDS/(R/Ie)                  # According to 12.8-2
        if Ta <= TL:
            Csmax = SD1/(Ta*R/Ie)            # According to 12.8-3
        elif Ta > TL:
            Csmax = SD1*TL/((Ta**2)*R/Ie)    # According to 12.8-4
        
    Csmin1 = max(0.044*SDS*Ie,0.01)      # According to 12.8-5
    
    if S1 >= 0.60:
        Csmin2 = 0.50*S1/(R/Ie)          # According to 12.8-6
    elif S1 < 0.60:
        Csmin2 = 0.00
    Csmin = max(Csmin1,Csmin2)
    
    if Csmax > Csmin:
        if Cscalc >= Csmax:
           Cs = Csmax
        elif Cscalc <= Csmin:
           Cs = Csmin
        else:
           Cs = Cscalc
    elif Csmax < Csmin:
        Cs = max(Csmin,Csmax)
       
    return Cs

def detCsDrifts(Ss,S1,Vs30,T,R,Ie,TL=8.0):
    """
    Function to determine the Seismic Response Coefficient according to Section
    12.8.1.1 on ASCE7-16

    Parameters
    ----------
    Ss : [g] MCE_R ground motion parameter for 0.20s Spectral response Acceleration.
    S1 : [g] MCE_R ground motion parameter for 1.0s Spectral response Acceleration.
    SiteClass : [STRING] Site Class according to ASCE7-10. Could be:
                        'A' , 'B' , 'C' , 'D' , 'E' or 'F'
    
    T : [s] Calculated or approximated period of vibration of the structure.
    R : [adim] Response modification coefficient according to Table 12.2-1.
    Ie : [adim] Importance factor according to Table 1.5-2. It could be:
                      For Risk Category       Seismic Importance Factor, Ie
                      -----------------        ---------------------------
                             I                           1.00
                             II                          1.00
                             III                         1.25
                             IV                          1.50
    TL : [s] Long period value according to Figures 22-12 through 22-16. 
             The default is 8.0.

    Returns
    -------
    Cs : [adim] Seismic Response Coefficient.

    """
    SDS = detSDS(Ss,Vs30)
    SD1 = detSD1(S1,Vs30)
    SiteClass = detSiteClass(Vs30)
    Ta = T
    (_,_,_,Ts) = Sa(Ss,S1,Vs30,T)
    
    Cscalc = SDS/(R/Ie)                  # According to 12.8-2
    if SiteClass in ('D') and S1 >= 0.20:   # According to 11.4.8
        if Ta <= 1.50*Ts:
            Csmax = SDS/(R/Ie)
        elif Ta > 1.50*Ts and Ta <= TL:
            Csmax = 1.5*SD1/(Ta*R/Ie)            # According to 12.8-3
        elif Ta > TL:
            Csmax = 1.5*SD1*TL/((Ta**2)*R/Ie)    # According to 12.8-4
    else:
        Cscalc = SDS/(R/Ie)                  # According to 12.8-2
        if Ta <= TL:
            Csmax = SD1/(Ta*R/Ie)            # According to 12.8-3
        elif Ta > TL:
            Csmax = SD1*TL/((Ta**2)*R/Ie)    # According to 12.8-4
        
    # Csmin1 = max(0.044*SDS*Ie,0.01)      # According to 12.8-5
    Csmin1 = 0.0                         # Suppresing the usage of 12.8-5 by 12.8.6.1
    
    if S1 >= 0.60:
        Csmin2 = 0.50*S1/(R/Ie)          # According to 12.8-6
    elif S1 < 0.60:
        Csmin2 = 0.00
    Csmin = max(Csmin1,Csmin2)
    
    if Csmax > Csmin:
        if Cscalc >= Csmax:
           Cs = Csmax
        elif Cscalc <= Csmin:
           Cs = Csmin
        else:
           Cs = Cscalc
    elif Csmax < Csmin:
        Cs = max(Csmin,Csmax)
               
    return Cs

def ELFP(Ss,S1,Vs30,Wx,R,Ie,hx,StrType,Tcalc=0.0,TL=8.0):
    """
    Function to determine seismic forces for static seismic analysis using
    procedures stablished in Section 12.8 of ASCE7-16, by the Equivalent Lateral
    Force Procedure.
Treal
    Parameters
    ----------
    Ss : [g] MCE_R ground motion parameter for 0.20s Spectral response Acceleration.
    S1 : [g] MCE_R ground motion parameter for 1.0s Spectral response Acceleration.
    SiteClass : [STRING] Site Class according to ASCE7-10. Could be:
                        'A' , 'B' , 'C' , 'D' , 'E' or 'F'
    Wx : [kN][LIST] is a list of 1xN elements in which one of them relates to
                    the weight of each floor above the base. The element from
                    the list, Wx[0,0] must be the weight of the first floor. 
                    The Wx[0,N-1]  shall be the weight of the roof.
                    N is the total number of stories above the base level.
                    Please notice that the N-1 index indicates the last position
                    of the list according to Python indexing notation.
    R : [adim] Response modification coefficient according to Table 12.2-1.
    Ie : [adim] Importance factor according to Table 1.5-2. It could be:
                      For Risk Category       Seismic Importance Factor, Ie
                      -----------------        ---------------------------
                             I                           1.00
                             II                          1.00
                             III                         1.25
                             IV                          1.50
    hx : [m][LIST] is a list of 1xN elemnts in which each element relates to
                   the height of each floor level above the base. The element
                   from the list, indexed as hx[0,0] must be the height of the
                   first floor level. The hx[0,N-1] shall be the height from the
                   base to the roof level.
                   N is the total number of stories above the base level.
                   Please notice that the N-1 index indicates the last position
                   of the list according to Python indexing notation.
    StrType : [INTEGER] it is an integer defining the structural system type
              according to Table 12.8-2 on ASCE7-10. The values are:
                  [0] For Steel moment resisting frames.
                  [1] For concrete moment resisting frames.
                  [2] Steel excentrically braced  frames.
                  [3] Steel buckling-restrained braced frames.
                  [4] All other structural systems.
    Tcalc : [s] Real period of the structure in direction of analysis, previously
                calculated.
                The default is 0.0.
    TL : [s] Long period value according to Figures 22-12 through 22-16. 
             The default is 8.0.

    Returns
    -------
    SDS : [g] DESIGN Spectral response Acceleration parameter for 0.20s.
    SD1 : [g] DESIGN Spectral response Acceleration parameter for 1.00s.
    Ta : [s] approximated structural fundamental period of vibration.
    Cu : [adim] Coefficient for upper limit on calculated period.
    Tmax : [s] maximum limited fundamental period of vibration.
    Cs :[adim] Seismic Response Coefficient.
    CsSa : [adim] Seismic Response Coefficient from Spectrum of Section 11.4.5.
    V : [kN] Base seismic shear value for seismic design.
    Cvx : [adim][LIST] Vertical distribution factor list of 1xN elements, where
                       N is the total number of stories above the base level.
                       The element indexed as Cvx[0,0] is the factor related to
                       the first story and Cvx[0,N-1] is the one of the roof.
    Fx : [kN][LIST] Vertical distribution o seismic forces; list of 1xN elements, 
                    where N is the total number of stories above the base level.
                    The element indexed as Fx[0,0] is the force related to
                    the first story and Fx[0,N-1] is the one of the roof.
    Vx : [kN][LIST] Vertical distribution o seismic story shear; list of 1xN elements, 
                    where N is the total number of stories above the base level.
                    The element indexed as Vx[0,0] is the force related to
                    the N story and Vx[0,N-1] is the one of the BASE.
                    Notice that Vx[0,N-1] must equals the Value of V of the 
                    base shear calculated above.
    To : [s] Limit period for the begining of constant acceleration zone.
    Ts : [s] Limit period for the end of constant acceleration zone.

    """
    # Previosu calculations related to spectral data and structural period.
    # ---------------------------------------------------------------------
    SDS = detSDS(Ss,Vs30)
    SD1 = detSD1(S1,Vs30)
    # Approximate period of vibration of the structure.
    # Here, if no period is previously calculated and entered as an argument,
    # i.e Tcalc = 0.0, the function works with the calculated approximate value
    # of period. If Tcalc is introduced, it works with the real value.
    # Ta is a good start for determining the forces for a first execution.
    Ta = detTa(hx[-1],StrType)
    if Tcalc == 0.0:
        T = Ta
    elif Tcalc != 0.0:
        T = Tcalc
    # Maximum period of vibration for the structure
    # ---------------------------------------------
    Cu = detCu(S1,Vs30)
    Tmax = Ta*Cu[0]
    # Period used for determining lateral forces for strength design and drifts
    # for dimensioning. Considering 12.8.6.2 of ASCE7-16.
    # ------------------------------------------------------------------------
    if T > Tmax:
        Tforces = Tmax
        Tdrifts = T
    elif T < Ta:
        Tforces = Ta
        Tdrifts = Ta
    else:
        Tforces = T
        Tdrifts = T
    
    # Cs Calculation
    # --------------
    CsF = detCs(Ss,S1,Vs30,Tforces,R,Ie,TL=8.0)  # for strength design Forces determining.
    CsD = detCsDrifts(Ss,S1,Vs30,Tdrifts,R,Ie,TL=8.0)    # for drift analysis force determining.

    # Base shear calculation.
    # =======================
    VF = CsF*sum(Wx)
    VD = CsD*sum(Wx)
    
    # Vertical distribution of base shear.
    # ------------------------------------
    # Calculation of k
    # -----------------
    # For FORCES
    # ----------
    if Tforces <= 0.50:
        kF = 1.0
    elif Tforces > 0.50 and Tforces < 2.50:
        kF = 0.75 + 0.50*Tforces
    elif Tforces >= 2.50:
        kF = 2.0
    # For DRIFTS
    # ----------
    if Tdrifts <= 0.50:
        kD = 1.0
    elif Tdrifts > 0.50 and Tdrifts < 2.50:
        kD = 0.75 + 0.50*Tdrifts
    elif Tdrifts >= 2.50:
        kD = 2.0
    
    # For STRENGTH design determinations
    # ----------------------------------
    num = []
    for i in range(len(hx)):    # from first floor to top.
        deno = Wx[i]*hx[i]**kF
        num.append(deno)
    Cvx = []
    for i in range(len(num)):   # from first floor to top.
        Cvx.append(num[i]/sum(num))
    FxF = []
    for i in range(len(Cvx)):
        FxF.append(Cvx[i]*VF)
    VxF = []
    VxF.append(FxF[-1])
    cont = -1
    for i in range(len(FxF)-1):
        cont -= 1
        VxF.append(VxF[i] + FxF[cont])
        
    # For DRIFT EVALUATION determinations
    # ----------------------------------
    num = []
    for i in range(len(hx)):
        deno = Wx[i]*hx[i]**kD
        num.append(deno)
    Cvx = []
    for i in range(len(num)):
        Cvx.append(num[i]/sum(num))
    FxD = []
    for i in range(len(Cvx)):
        FxD.append(Cvx[i]*VD)
    VxD = []
    VxD.append(FxD[-1])
    cont = -1
    for i in range(len(FxD)-1):
        cont -= 1
        VxD.append(VxD[i] + FxD[cont])
    
    # Just for controling the value of V and Vx to be propperly calculated.
    # =====================================================================
    # difV = Vx[-1] - V
    # if difV <= 0.01:
    #     print('V = Vx[-1] verified!')
    # else:
    #     print('Vx not propperly calculated')
    
    return SDS, SD1, Ta, Tforces, Tdrifts, Cu[0], Tmax, CsF, CsD, VF, VD, FxF, FxD, VxF, VxD

def detMo(Fx,hx):
    """
    Function to calculate the Overturning Moment at the base of the structure
    determined by multiplying the applied seismic forces of each floor by its
    height from the base and then summing all of them.

    Parameters
    ----------
    Fx : [kN][LIST] Vertical distribution o seismic forces; list of 1xN elements, 
                    where N is the total number of stories above the base level.
                    The element indexed as Fx[0,0] is the force related to
                    the first story and Fx[0,N-1] is the one of the roof.
    hx : [m][LIST] is a list of 1xN elemnts in which each element relates to
                   the height of each floor level above the base. The element
                   from the list, indexed as hx[0,0] must be the height of the
                   first flor level. The hx[0,N-1] shall be the height from the
                   base to the roof level.
                   N is the total number of stories above the base level.
                   Please notice that the N-1 index indicates the last position
                   of the list according to Python indexing notation.

    Returns
    -------
    Mo : [kN-m] the overturning moment at the base using the unmodified 
                seismic forces and not including the reduction permitted in the 
                design of the foundation.

    """
    M = []
    for i in range(len(Fx)):
        M.append(Fx[i]*hx[i])
    Mo = sum(M)
    return Mo

def detdelta_a(SDC,RiskCat):
    """
    

    Parameters
    ----------
    SDC : TYPE
        DESCRIPTION.
    RiskCat : TYPE
        DESCRIPTION.

    Returns
    -------
    alldelta : TYPE
        DESCRIPTION.

    """
    redunrho = RedFact(SDC)
    if RiskCat in ('I','II'):
        alldelta = 0.02/redunrho
    elif RiskCat in ('III'):
        alldelta = 0.015/redunrho
    elif RiskCat in ('IV'):
        alldelta = 0.010/redunrho
    
    return alldelta

def vertdistV(V,T,hx,Wx):
    # Vertical distribution of V deppending on period of vibration according to
    # ASCE 7-16.
    # Vertical distribution of base shear.
    # ------------------------------------
    # Calculation of k
    # -----------------
    if T <= 0.50:
        k = 1.0
    elif T > 0.50 and T < 2.50:
        k = 0.75 + 0.50*T
    elif T >= 2.50:
        k = 2.0
    
    # For STRENGTH design determinations
    # ----------------------------------
    num = []
    for i in range(len(hx)):
        deno = Wx[i]*hx[i]**k
        num.append(deno)
    Cvx = []
    for i in range(len(num)):
        Cvx.append(num[i]/sum(num))
    Fx = []
    for i in range(len(Cvx)):
        Fx.append(Cvx[i]*V)
    Vx = []
    Vx.append(Fx[-1])
    cont = -1
    for i in range(len(Fx)-1):
        cont -= 1
        Vx.append(Vx[i] + Fx[cont])
        
    return Fx, Vx