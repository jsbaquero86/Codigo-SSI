# -*- coding: utf-8 -*-
"""
Created on Thu May 14 14:57:04 2020

@author: sebas
"""

"""
    Definitions for functions according to ASCE41-17.
"""


def round5M(x):
    """
    Function that rounds values to the nearest major detltiple of 0.05. Usefull for
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

def detBPOE(RiskCat):
    # List of Structural Performance Levels
    # --------------------------------------
    SPL = ['S-1','S-2','S-3','S-4','S-5']
    """
    S-1: Immediate Occupancy Structural Performance Level.
    S-2: Damage Control Structural Performance Level.
    S-3: Life Safety Structural Performance Level.
    S-4: Limited Safety Structural Performance Level.
    S-5: Collapse Prevention Structural Performance Level.
    """
    # List of Seismic Hazard Level
    # ----------------------------
    SHL = ['BSE-1E','BSE-2E']
    
    if RiskCat in ('I','II'):
        BPOE = dict([(SHL[0],SPL[2]),(SHL[1],SPL[4])])
    elif RiskCat in ('III'):
        BPOE = dict([(SHL[0],SPL[1]),(SHL[1],SPL[3])])
    elif RiskCat in ('IV'):
        BPOE = dict([(SHL[0],SPL[0]),(SHL[1],SPL[2])])
        
    return BPOE


def detBPON(RiskCat):
    # List of Structural Performance Levels
    # --------------------------------------
    SPL = ['S-1','S-2','S-3','S-4','S-5']
    """
    S-1: Immediate Occupancy Structural Performance Level.
    S-2: Damage Control Structural Performance Level.
    S-3: Life Safety Structural Performance Level.
    S-4: Limited Safety Structural Performance Level.
    S-5: Collapse Prevention Structural Performance Level.
    """
    # List of Seismic Hazard Level
    # ----------------------------
    SHL = ['BSE-1N','BSE-2N']
    
    if RiskCat in ('I','II'):
        BPON = dict([(SHL[0],SPL[2]),(SHL[1],SPL[4])])
    elif RiskCat in ('III'):
        BPON = dict([(SHL[0],SPL[1]),(SHL[1],SPL[3])])
    elif RiskCat in ('IV'):
        BPON = dict([(SHL[0],SPL[0]),(SHL[1],SPL[2])])
        
    return BPON


def detFa(Ss,Vs30):

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



def detSs(lon,lat,curvetype='RT',PAi=1/2475):
    # Determination of Ss values (T = 0.2s) for a determined location and 
    # for an specific Annual Frequency of Exceedance.
    # curvetype could be 'RT' per Risk-Targeted MCE or 'PEF' for probability
    # of exceedance of the a determined hazard in a frequency as was used
    # in previous issues of ASCE7.
    
    import pandas as pd
    import numpy as np
    
    rlon = round(round5M(lon),2)
    rlat = round(round5M(lat),2)
    
    # Obtaining Data Frame with Pandas
    # --------------------------------
    if curvetype in ('RT'):
        df = pd.read_csv('2010_ASCE-7-US-Ss-S1-PGA-0p05.csv')
        dfarray = df.to_numpy()
        
        # Determining geographical location index within the dataFrame.
        # -------------------------------------------------------------
        cocount = 0
        while rlon != dfarray[cocount,1] or rlat != dfarray[cocount,0]:
            cocount +=1
        
        Ss = dfarray[cocount,2]/100
        
    elif curvetype in ('PEF'):
        df = pd.read_csv('US_2018_nshmp_hazard_curves_BC_Vs30_760_0p2s.csv')
        dfarray = df.to_numpy()
        accels = dfarray[0,2:]
        
        # Determining geographical location index within the dataFrame.
        # -------------------------------------------------------------
        cocount = 0
        while rlon != dfarray[cocount,0] or rlat != dfarray[cocount,1]:
            cocount +=1
        # Defining vector of Hazard Curve for the specific location.
        PAs = dfarray[cocount,2:]
        
        # Determining indexing for PAi within the hazard curve values
        # ----------------------------------------------------------
        pacount = 0
        while PAi < PAs[pacount]:
            pacount += 1
        
        # Values of annual frequency of exceedance for interpolation
        # ----------------------------------------------------------
        PA2 = PAs[pacount]
        PA1 = PAs[pacount-1]
        # Values of Ground Motion for interpolation
        # -----------------------------------------
        GM2 = accels[pacount]
        GM1 = accels[pacount-1]
        
        # Logaritmic interpolation of PAi
        # --------------------------------
        X = 10**GM1 + (10**PAi - 10**PA1)*(10**GM2 - 10**GM1)/(10**PA2 - 10**PA1)
        Ss = np.log10(X)
    
    return Ss


def detS1(lon,lat,curvetype='RT',PAi=1/2475):
    # Determination of Ss values (T = 1.0s) for a determined location and 
    # for an specific Annual Frequency of Exceedance.
    # curvetype could be 'RT' per Risk-Targeted MCE or 'PEF' for probability
    # of exceedance of the a determined hazard in a frequency as was used
    # in previous issues of ASCE7.
    
    import pandas as pd
    import numpy as np
    
    rlon = round(round5M(lon),2)
    rlat = round(round5M(lat),2)
    
    # Obtaining Data Frame with Pandas
    # --------------------------------
    if curvetype in ('RT'):
        df = pd.read_csv('2010_ASCE-7-US-Ss-S1-PGA-0p05.csv')
        dfarray = df.to_numpy()
        
        # Determining geographical location index within the dataFrame.
        # -------------------------------------------------------------
        cocount = 0
        while rlon != dfarray[cocount,1] or rlat != dfarray[cocount,0]:
            cocount +=1
        
        S1 = dfarray[cocount,3]/100
        
    elif curvetype in ('PEF'):
        df = pd.read_csv('US_2018_nshmp_hazard_curves_BC_Vs30_760_1s.csv')
        dfarray = df.to_numpy()
        accels = dfarray[0,2:]
        
        # Determining geographical location index within the dataFrame.
        # -------------------------------------------------------------
        cocount = 0
        while rlon != dfarray[cocount,0] or rlat != dfarray[cocount,1]:
            cocount +=1
        # Defining vector of Hazard Curve for the specific location.
        PAs = dfarray[cocount,2:]
        
        # Determining indexing for PAi within the hazard curve values
        # ----------------------------------------------------------
        pacount = 0
        while PAi < PAs[pacount]:
            pacount += 1
        
        # Values of annual frequency of exceedance for interpolation
        # ----------------------------------------------------------
        PA2 = PAs[pacount]
        PA1 = PAs[pacount-1]
        # Values of Ground Motion for interpolation
        # -----------------------------------------
        GM2 = accels[pacount]
        GM1 = accels[pacount-1]
        
        # Logaritmic interpolation of PAi
        # --------------------------------
        X = 10**GM1 + (10**PAi - 10**PA1)*(10**GM2 - 10**GM1)/(10**PA2 - 10**PA1)
        S1 = np.log10(X)
    
    return S1



def detSxi(lon,lat,Vs30,SHL):
    # Determination of design spectral acceleration parameters Sxs and Sx1
    # according to Chapter 2 of ASCE41-17. It depends on the Seismic Hazard
    # Level nitroduced as input argument.
    
    # DETERMINATION OF RISK-TAGETED SPECTRAL COEFFICIENTS
    # ==================================================
    SsRT = detSs(lon,lat,curvetype='RT')
    S1RT = detS1(lon,lat,curvetype='RT')
    
    # DETERMINATION OF SITE COEFFICIENTS FOR ADJUSTING PURPOSES.
    # ==========================================================
    Fa = detFa(SsRT,Vs30) 
    Fv = detFv(S1RT,Vs30)
    
    if SHL in ('BSE-2N'):       # SHL: Seismic Hazard Level 
        Sxs = 2/3*Fa*SsRT
        Sx1 = 2/3*Fv*S1RT
    elif SHL in ('BSE-1N'):
        Sxs = 4/9*Fa*SsRT
        Sx1 = 4/9*Fv*S1RT
    elif SHL in ('BSE-2E'):
        # Determining Uniform Hazard Response Spectrum coefficients for
        # 5% of probability of exceedance in 50 years.
        # --------------------------------------------
        Ss_5_50 = detSs(lon,lat,curvetype='PEF',PAi=0.00102534)
        S1_5_50 = detS1(lon,lat,curvetype='PEF',PAi=0.00102534)
        # Determining Design coefficients
        # -------------------------------
        Sxs = Fa*Ss_5_50
        Sx1 = Fv*S1_5_50
        if Sxs > 2/3*Fa*SsRT:
            Sxs = 2/3*Fa*SsRT
        if Sx1 > 2/3*Fv*S1RT:
            Sx1 = 2/3*Fv*S1RT
    elif SHL in ('BSE-1E'):
        # Determining Uniform Hazard Response Spectrum coefficients for
        # 20% of probability of exceedance in 50 years.
        # --------------------------------------------
        Ss_20_50 = detSs(lon,lat,curvetype='PEF',PAi=0.00445293)
        S1_20_50 = detS1(lon,lat,curvetype='PEF',PAi=0.00445293)
        # Determining Design coefficients
        # -------------------------------
        Sxs = Fa*Ss_20_50
        Sx1 = Fv*S1_20_50
        if Sxs > 4/9*Fa*SsRT:
            Sxs = 4/9*Fa*SsRT
        if Sx1 > 4/9*Fv*S1RT:
            Sx1 = 4/9*Fv*S1RT
        
    return Sxs, Sx1


def detGHRS(Sxs,Sx1,T,beta=0.05,TL=8.,Tf=6.0,Npts=500):
    # Function definition to construct the General Horizontal Response Spectrum
    # (GHRS).
    import numpy as np
    # Previous calculations.
    # ----------------------
    T0 = 0.20*(Sx1/Sxs)
    Ts = Sx1/Sxs
    B1 = 4/(5.6 - np.log(100*beta))
    Sds = Sxs/B1
    Sd1 = Sx1/B1
    PGA = 0.40*Sxs
    
    # Specific Sa value, depending on fundamental period of vibration.
    # ----------------------------------------------------------------
    if T > 0 and T < T0:
        Sa = Sxs*((5/B1 - 2)*T/Ts + 0.40)
    elif T >= T0 and T <= Ts:
        Sa = Sds
    elif T > Ts and T <= TL:
        Sa = Sd1/T
    elif T > TL:
        Sa = Sd1*TL/T**2
        
    # Data series for spectrum plotting.
    # ----------------------------------
    SaData = []
    periods = [Tf*x/Npts for x in range(Npts)]
    periods = periods + [T0,Ts]
    if Tf > TL:
        periods += [TL]
    periods.sort()
    for i in range(len(periods)):
        if periods[i] >= 0 and periods[i] < T0:
            SaData.append(Sxs*((5/B1 - 2)*periods[i]/Ts + 0.40))
        elif periods[i] >= T0 and periods[i] <= Ts:
            SaData.append(Sds)
        elif periods[i] > Ts and periods[i] <= TL:
            SaData.append(Sd1/periods[i])
        elif periods[i] > TL:
            SaData.append(Sd1*TL/periods[i]**2)
            
    return Sa, Sds, Sd1, PGA, [T,T0,Ts], periods, SaData


def detSeismLevel(lon,lat,Vs30):
    # Function definition to determine Seimicity Level according to ASCE41-17
    
    # Spectral acceleration parameters determination
    # -----------------------------------------------
    SsRT = detSs(lon,lat,curvetype='RT')
    S1RT = detS1(lon,lat,curvetype='RT')
    
    # Site Class factors
    # -------------------
    Fa = detFa(SsRT, Vs30)
    Fv = detFv(S1RT, Vs30)
    
    # Design spectral acceleration parameters
    # ---------------------------------------
    SDS = 2/3*SsRT*Fa
    SD1 = 2/3*S1RT*Fv
    
    # Defining numeric value for Seimicity Level
    # ------------------------------------------
    """
    In this case, we are going to set numerical values to each level, as follows:
        [0] - Very Low
        [1] - Low
        [2] - Moderate
        [3] - High
    """
    num = [0,1,2,3]
    SLevel = ['Very Low','Low','Moderate','High']
    SLdict = dict(zip(num,SLevel))
    
    # Considering SDS
    # ---------------
    if SDS < 0.167:
        SL_SDS = 0
    elif SDS >= 0.167 and SDS < 0.33:
        SL_SDS = 1
    elif SDS >= 0.33 and SDS < 0.50:
        SL_SDS = 2
    elif SDS >= 0.50:
        SL_SDS = 3
        
    # Considering SD1
    # ----------------
    if SD1 < 0.067:
        SL_SD1 = 0
    elif SD1 >= 0.067 and SD1 < 0.133:
        SL_SD1 = 1
    elif SD1 >= 0.133 and SD1 < 0.20:
        SL_SD1 = 2
    elif SD1 >= 0.20:
        SL_SD1 = 3
    
    SLindex = max(SL_SDS,SL_SD1)
    SL = SLdict[SLindex]
    
    return SLindex, SL


def detkappa(LOK,PO,CA):
    # For determining the knowledge factor, kappa, to account for the
    # uncertainty in the collection of data the Selected Performance Objective,
    # analysis procedure and data collection process shall be used.
    # The following nomenclature is required as input.
    # LOK (Level of Knowledge):
        # - 'M' for minumum.
        # - 'U' for usual.
        # - 'C' for comprehensive.
    # PO (Performance Objective):
        # - 'S-3' for Life safety or lower.
        # - 'S-2' for damage control or lower.
        # - 'S-1' for Immediate Occupancy or lower.
    # CA (Condition Assessment):
        # - 'V' for visual assessment.
        # - 'Cp' for comprehensive assessment.
    
    if LOK in ('M'):
        if PO in ('S-3','S-4','S-5','S-6'):
            if CA in ('V'):
                kappa = 0.90
            elif CA in ('Cp'):
                kappa = 0.75
        else:
            kappa = 0.75
            print('No such PO permitted for the level of knowledge.')
    elif LOK in ('U','C'):
        kappa = 1.0
        
    return kappa



def detCm(NStr,BldTypCm,Te):
    # Effective mass factor to account for higher modal mass participation 
    # effects obtained from Table 7-4.
    
    # Cm factor calculation according to Table 7-4 and controled by Te.
    # -----------------------------------------------------------------
    if BldTypCm in ('C-MF', 'C-SHW', 'C-PS', 'S-MF', 'S-CBF', 'S-EBF', 'Other') and NStr <= 2:
        Cm = 1.0
    elif BldTypCm in ('C-MF', 'S-CBF', 'S-CBF', 'S-EBF') and NStr >= 3:
        Cm = 0.90
    elif BldTypCm in ('C-SHW', 'C-PS') and NStr >= 3:
        Cm = 0.80
    elif BldTypCm in ('Other') and NStr >= 3:
        Cm = 1.0
    if Te > 1.0:
        Cm = 1.0
    
    return Cm



def mu_strength(lon,lat,Vs30,SHL,T,W,R,Ie,Cm,beta=0.05,TL=8.,Tf=6.0,Npts=500,Vy=0,Te=0):
    # Ratio Ratio of elastic strength demand to yield strength coefficient 
    # calculated in accordance with Eq. (7-31) with the elastic base shear 
    # capacity substituted for shear yield strength, Vy.
    # For this definition, the Te value is the effective fundamental period of
    # the structure calculated as in Eq.. (7-27) of ASCE41-17. For the first
    # approximation of NSP response, the value used for Te is the real dynamicaly
    # calculated period, T; and then calculated again with the obtained Te value.
    
    # Definition of whether to use T or Te for calculations.
    if Te != 0:
        T = Te
    
    import ASCE716Seismic as ASCE716
    
    
    # Determination ratio of elastic strength demand
    # mu_str = Sa/(Vy/W)*Cm
    # According to ASCE 7-16 - Section C12.1.1:
        # The result is that structures typically have a much higher
        # lateral strength than that specified as the minimum by the standard, 
        # and the first significant yielding of structures may occur at lateral 
        # load levels that are 30% to 100% higher than the prescribed design 
        # seismic forces. This assumption is used when an estimate of the parameter
        # is required for the first run of NSP. After that point, one can already
        # define a value of Vy from the idealized curve of Section 7.4.3.2.4.
    # Determination of sesimic spectral acceleration parameters.
    (Sxs,Sx1) = detSxi(lon,lat,Vs30,SHL)
    # Determination of Sa value according to SHL (Seismic Hazard Level).
    (Sa,_,_,_,_,_,_) = detGHRS(Sxs,Sx1,T,beta,TL,Tf,Npts)
    if Vy == 0:
        # Determination of Seismic Coefficient value, Cs, to Estimate Vy.
        Cs = ASCE716.detCs(Sxs,Sx1,Vs30,T,R,Ie,TL)
        mu_str = Sa/(1.65*Cs)*Cm
    else:
        mu_str = Sa/(Vy/W)*Cm
    
    return mu_str



def det_afactor(Vs30):
    # Determination of a factor which is a Site Class factor.
    
    import ASCE716Seismic as ASCE716
    
    # Site Class determination
    SiteClass = ASCE716.detSiteClass(Vs30)
    # a factor determination according to site Class. Section 7.4.3.3.2 ASCE41-13
    if SiteClass in ('A', 'B'):
        a = 130
    elif SiteClass in ('C'):
        a = 90
    elif SiteClass in ('D', 'E', 'F'):
        a = 60
        
    return a



def detCo(NStr,BldTypCo):
    # Determination of  Co value according to Table 7-5 of ASCE4117.
    # Co factor dtermination from Table 7-5 ASCE41-13.
    # ------------------------------------------------
    if BldTypCo in ('Sh-T', 'Sh-U', 'O-A') and NStr == 1:
        Co = 1.0
    elif BldTypCo in ('Sh-T', 'O-A') and NStr == 2:
        Co = 1.20
    elif BldTypCo in ('Sh-U') and NStr == 2:
        Co = 1.15
    elif BldTypCo in ('Sh-T', 'Sh-U') and NStr == 3:
        Co = 1.20
    elif BldTypCo in ('O-A') and NStr == 3:
        Co = 1.30
    elif BldTypCo in ('Sh-T') and NStr == 4:
        Co = 1.25
    elif BldTypCo in ('Sh-U') and NStr == 4:
        Co = 1.20
    elif BldTypCo in ('O-A') and NStr == 4:
        Co = 1.35
    elif BldTypCo in ('Sh-T') and NStr >= 5:
        Co = 1.30
    elif BldTypCo in ('Sh-U') and NStr >= 5:
        Co = 1.20
    elif BldTypCo in ('O-A') and NStr >= 5 and NStr < 10:
        Co = 1.40
    elif BldTypCo in ('O-A') and NStr > 10:
        Co = 1.50

    return Co



def detC1(lon,lat,Vs30,SHL,T,W,R,Ie,Cm,beta=0.05,TL=8.,Tf=6.0,Npts=500,Vy=0,Te=0):
    # To determine C1 modification factor according to 7.4.1.3.1 pf ASCE41-17.
    # Modification factor to relate expected maximum inelastic displacements 
    # to displacements calculated for linear elastic response.
    # For this definition, the Te value is the effective fundamental period of
    # the structure calculated as in Eq.. (7-27) of ASCE41-17. For the first
    # approximation of NSP response, the value used for Te is the real dynamicaly
    # calculated period, T; and then calculated again with the obtained Te value.
    
    # Definition of whether to use T or Te for calculations.
    if Te != 0:
        T = Te
    
    if T < 0.20:
        # Determination of mu_str parameter
        mu_str = mu_strength(lon,lat,Vs30,SHL,T,W,R,Ie,Cm,beta,TL,Tf,Npts,Vy,T)
        # Determination of a_factor.
        a = det_afactor(Vs30)
        # Calculation of C1 factor.
        C1_T = 1 + (mu_str - 1)/(a*T**2)
        
        # Same calculations but with T = 0.2
        mu_str = mu_strength(lon,lat,Vs30,SHL,T,W,R,Ie,Cm,beta,TL,Tf,Npts,Vy,T)
        C1_02 = 1 + (mu_str - 1)/(a*0.20**2)
        C1 = min(C1_T,C1_02)
        
    elif T >= 0.20 and T <= 1.0:
        # Determination of mu_str parameter
        mu_str = mu_strength(lon,lat,Vs30,SHL,T,W,R,Ie,Cm,beta,TL,Tf,Npts,Vy,T)
        # Determination of a_factor.
        a = det_afactor(Vs30)
        # Calculation of C1 factor.
        C1 = 1 + (mu_str - 1)/(a*T**2)
    elif T > 1.0:
        C1 = 1.0
    
    return C1



def detC2(lon,lat,Vs30,SHL,T,W,R,Ie,Cm,beta=0.05,TL=8.,Tf=6.0,Npts=500,Vy=0,Te=0):
    # Modification factor to represent the effect of pinched hysteresis shape, 
    # cyclic stiffness degradation, and strength deterioration on maximum 
    # displacement response.
    # For this definition, the Te value is the effective fundamental period of
    # the structure calculated as in Eq.. (7-27) of ASCE41-17. For the first
    # approximation of NSP response, the value used for Te is the real dynamicaly
    # calculated period, T; and then calculated again with the obtained Te value.
    
    # Definition of whether to use T or Te for calculations.
    if Te != 0:
        T = Te
    
    if T <= 0.70:
        # Determination of mu_str ratio.
        mu_str = mu_strength(lon,lat,Vs30,SHL,T,W,R,Ie,Cm,beta,TL,Tf,Npts,Vy,T)
        # Calculation of C2 factor.
        C2 = 1 + 1/800*((mu_str - 1)/T)**2
    else:
        C2 = 1.0
    
    return C2



def LSPforces(lon,lat,Vs30,SHL,NStr,BldTypCm,T,Wx,hx,R,Ie,beta=0.05,TL=8.,Tf=6.0,Npts=500,Vy=0,Te=0):
    # Linear static procedure determined as per ASCE41-17.
    # For this definition, the Te value is the effective fundamental period of
    # the structure calculated as in Eq.. (7-27) of ASCE41-17. For the first
    # approximation of NSP response, the value used for Te is the real dynamicaly
    # calculated period, T; and then calculated again with the obtained Te value.
    
    import numpy as np
    
    # Definition of whether to use T or Te for calculations.
    if Te != 0:
        T = Te
        
    # Period determination
    # ====================
    # The period determination is already obtained with the Analytical Method
    # or Method 1 according to ASCE41.
    
    W = sum(Wx)
    # Modification factors.
    # ---------------------
    Cm = detCm(NStr,BldTypCm,T)
    C1 = detC1(lon,lat,Vs30,SHL,T,W,R,Ie,Cm,beta,TL,Tf,Npts,Vy,T)
    C2 = detC2(lon,lat,Vs30,SHL,T,W,R,Ie,Cm,beta,TL,Tf,Npts,Vy,T)
    
    
    # Determination of spectral acceleration, Sa
    # ------------------------------------------
    (Sxs,Sx1) = detSxi(lon,lat,Vs30,SHL)
    (Sa,_,_,_,_,_,_) = detGHRS(Sxs,Sx1,T,beta,TL,Tf,Npts)
    
    # Pseudo Seismic force for LSP
    # ============================
    V = C1*C2*Cm*Sa*W
    
    # Vertical distribution of of Seismic Forces
    # ===========================================
    Cvx = np.zeros(len(hx))
    # Calculation of k
    if T <= 0.50:
        kF = 1.0
    elif T > 0.50 and T < 2.50:
        kF = 0.75 + 0.50*T
    elif T >= 2.50:
        kF = 2.0
    # Determination of numerator vector
    numerator = []
    for i in range(len(hx)):
        numerator.append(Wx[i]*hx[i]**kF)
    # Calculation of Cvx
    for i in range(len(hx)):
        Cvx[i] = numerator[i]/sum(numerator)
    # Calculation of vertical distribution of shear base.
    # ---------------------------------------------------
    Fx = V*Cvx
    
    # Determination of Shear Distribution for each floor.
    # ===================================================
    Vx = np.zeros(len(Fx))
    for i in range(len(hx)):
        Vx[i] = sum(Fx[i:])
    
    return V, Fx, Vx



def detCMCR(NbaysX,NbaysZ,NStr,XbayL,ZbayL,StoryH,XBDims,ZBDims,ColDims,e_slb=0.25,flex='no'):
    
    import numpy as np
    
    # Previous calculations.
    # ----------------------
    NPX = NbaysX+1                   # [u] number of nodes in X direction.
    NPZ = NbaysZ+1                   # [u] number of nodes in Z direction.
    PPF = NPX*NPZ                    # [u] number of nodes/columns in floor plan.
    NbeamsX = NbaysX*NPZ             # [u] number of beams in X direction at a level.
    NbeamsZ = NbaysZ*NPX             # [u] number of beams in Z direction ar a level.
    B = ZbayL*NbaysZ                 # [m] short plan dimension of building.
    L = XbayL*NbaysX                 # [m] long plan dimension of building.
    
    
    # =================
    # CENTER OF MASSES
    # =================
    if flex in ('Y','YES','Yes','yES','yes','y'):
        CM_levelcoord = np.zeros((NStr+1,2))    # pre-definition of level CM coordinates.
        CM_cumulcoord = np.zeros((NStr+1,2))    # pre-definition of cumulative CM coordinates.
    else:
        CM_levelcoord = np.zeros((NStr,2))    # pre-definition of level CM coordinates.
        CM_cumulcoord = np.zeros((NStr,2))    # pre-definition of cumulative CM coordinates.
    
    # Volumes from area elements (slabs)
    # ---------------------------------
    # Assuming so far, regular rectangular buildings, slabs have the same
    # dimensions as building plan, and is common to all levels (storeys).
    A_slb = B*L                       # [m2] area of plan slab.
    Vslb = A_slb*e_slb                # [m3] volum of slab component.
    Coord_slb= [L/2,B/2]              # [m][LIST][X,Z] coordinates of CM slab component.
    
    # Beams in X direction for each floor
    # -----------------------------------
    # Pre-definition of matrices for locating values of volume and coordinates.
    XbeamsVols = np.zeros((NStr,NbeamsX))
    XbeamsXCoord = np.zeros((NStr,NbeamsX))
    XbeamsZCoord = np.zeros((NStr,NbeamsX))
    num = 0
    
    for i in range(NStr):
        for j in range(NbeamsX):
            
            # Volume calculation
            XbeamsVols[i,j] = XBDims[num,0]*XBDims[num,1]*XbayL
            # Coordinates determination
            if (j+1)%NbaysX == 0:
                iNode = (NbaysX - 1)*XbayL
                jNode = NbaysX*XbayL
            else:
                iNode = ((j+1)%NbaysX-1)*XbayL
                jNode = ((j+1)%NbaysX)*XbayL
            XbeamsXCoord[i,j] = 0.5*(iNode + jNode)
            XbeamsZCoord[i,j] = np.floor(j/NbaysX)*ZbayL
            num += 1
            
    # Beams in Z direction for each floor
    # -----------------------------------
    # Pre-definition of matrices for allocating values of volume and coordinates.
    ZbeamsVols = np.zeros((NStr,NbeamsZ))
    ZbeamsXCoord = np.zeros((NStr,NbeamsZ))
    ZbeamsZCoord = np.zeros((NStr,NbeamsZ))
    num = 0
    
    for i in range(NStr):
        for j in range(NbeamsZ):
            
            # Volume calculation
            ZbeamsVols[i,j] = ZBDims[num,0]*ZBDims[num,1]*ZbayL
            # Coordinates determination
            if (j+1)%NbaysZ == 0:
                iNode = (NbaysZ - 1)*ZbayL
                jNode = NbaysZ*ZbayL
            else:
                iNode = ((j+1)%NbaysZ - 1)*ZbayL
                jNode = ((j+1)%NbaysZ)*ZbayL
            ZbeamsXCoord[i,j] = np.floor(j/NbaysZ)*XbayL
            ZbeamsZCoord[i,j] = 0.5*(iNode + jNode)
            num += 1
            
    # Beams in Columns for each floor
    # -----------------------------------
    # Pre-definition of matrices for allocating values of volume and coordinates.
    ColsVols = np.zeros((NStr,PPF))
    ColsXCoord = np.zeros((NStr,PPF))
    ColsZCoord = np.zeros((NStr,PPF))
    num = 0
    
    for i in range(NStr):
        colindex = 0
        for j in range(NPZ):
            for k in range(NPX):
            
                # Volume calculation
                if i == (NStr-1):
                    ColsVols[i,colindex] = ColDims[num,0]*ColDims[num,1]*StoryH/2
                else:
                    ColsVols[i,colindex] = ColDims[num,0]*ColDims[num,1]*StoryH
                
                # Coordinates determination
                ColsXCoord[i,colindex] = k*XbayL
                ColsZCoord[i,colindex] = j*ZbayL
            
                colindex += 1
        num += 1
        
    # DETERMINATION OF CENTER OF MASS
    # ================================
    if flex in ('Y','YES','Yes','yES','yes','y'):
        CM_levelcoord[0,0] = sum(ColsVols[0,:]*ColsXCoord[0,:])/sum(ColsVols[0,:])
        CM_levelcoord[0,1] = sum(ColsVols[0,:]*ColsZCoord[0,:])/sum(ColsVols[0,:])
        for i in range(NStr):
            # X Coordinates
            CM_Xbeams = sum(XbeamsVols[i,:]*XbeamsXCoord[i,:])/sum(XbeamsVols[i,:])
            CM_Zbeams = sum(ZbeamsVols[i,:]*ZbeamsXCoord[i,:])/sum(ZbeamsVols[i,:])
            CM_cols = sum(ColsVols[i,:]*ColsXCoord[i,:])/sum(ColsVols[i,:])
            CM_levelcoord[i+1,0] = (CM_Xbeams*sum(XbeamsVols[i,:]) + CM_Zbeams*sum(ZbeamsVols[i,:])\
                                   + CM_cols*sum(ColsVols[i,:]) + Coord_slb[0]*Vslb)/(sum(XbeamsVols[i,:]) + \
                                                                   sum(ZbeamsVols[i,:]) + \
                                                                       sum(ColsVols[i,:]) + Vslb)
            # Z coordinates
            CM_Xbeams = sum(XbeamsVols[i,:]*XbeamsZCoord[i,:])/sum(XbeamsVols[i,:])
            CM_Zbeams = sum(ZbeamsVols[i,:]*ZbeamsZCoord[i,:])/sum(ZbeamsVols[i,:])
            CM_cols = sum(ColsVols[i,:]*ColsZCoord[i,:])/sum(ColsVols[i,:])
            CM_levelcoord[i+1,1] = (CM_Xbeams*sum(XbeamsVols[i,:]) + CM_Zbeams*sum(ZbeamsVols[i,:])\
                                   + CM_cols*sum(ColsVols[i,:]) + Coord_slb[1]*Vslb)/(sum(XbeamsVols[i,:]) + \
                                                                   sum(ZbeamsVols[i,:]) + \
                                                                       sum(ColsVols[i,:]) + Vslb)
    else:
        for i in range(NStr):
            # X Coordinates
            CM_Xbeams = sum(XbeamsVols[i,:]*XbeamsXCoord[i,:])/sum(XbeamsVols[i,:])
            CM_Zbeams = sum(ZbeamsVols[i,:]*ZbeamsXCoord[i,:])/sum(ZbeamsVols[i,:])
            CM_cols = sum(ColsVols[i,:]*ColsXCoord[i,:])/sum(ColsVols[i,:])
            CM_levelcoord[i,0] = (CM_Xbeams*sum(XbeamsVols[i,:]) + CM_Zbeams*sum(ZbeamsVols[i,:])\
                                   + CM_cols*sum(ColsVols[i,:]) + Coord_slb[0]*Vslb)/(sum(XbeamsVols[i,:]) + \
                                                                   sum(ZbeamsVols[i,:]) + \
                                                                       sum(ColsVols[i,:]) + Vslb)
            # Z coordinates
            CM_Xbeams = sum(XbeamsVols[i,:]*XbeamsZCoord[i,:])/sum(XbeamsVols[i,:])
            CM_Zbeams = sum(ZbeamsVols[i,:]*ZbeamsZCoord[i,:])/sum(ZbeamsVols[i,:])
            CM_cols = sum(ColsVols[i,:]*ColsZCoord[i,:])/sum(ColsVols[i,:])
            CM_levelcoord[i,1] = (CM_Xbeams*sum(XbeamsVols[i,:]) + CM_Zbeams*sum(ZbeamsVols[i,:])\
                                   + CM_cols*sum(ColsVols[i,:]) + Coord_slb[1]*Vslb)/(sum(XbeamsVols[i,:]) + \
                                                                   sum(ZbeamsVols[i,:]) + \
                                                                       sum(ColsVols[i,:]) + Vslb)
        
    # DETERMINATION OF THE CUMULATIVE POSITION OF CENTER OF MASS
    # =========================================================
    if flex in ('Y','YES','Yes','yES','yes','y'):
        floorVols = np.zeros((NStr+1,1))
        floorVols[0,0] = sum(ColsVols[i,:])/2
        for i in range(NStr):
            floorVols[i+1,0] = sum(XbeamsVols[i,:]) + sum(ZbeamsVols[i,:]) + sum(ColsVols[i,:]) + Vslb
        for i in range(NStr+1):
            CM_cumulcoord[i,0] = sum(floorVols[i:,0]*CM_levelcoord[i:,0])/sum(floorVols[i:,0])
            CM_cumulcoord[i,1] = sum(floorVols[i:,0]*CM_levelcoord[i:,1])/sum(floorVols[i:,0])
    else:
        floorVols = np.zeros((NStr,1))
        for i in range(NStr):
            floorVols[i,0] = sum(XbeamsVols[i,:]) + sum(ZbeamsVols[i,:]) + sum(ColsVols[i,:]) + Vslb
        for i in range(NStr):
            CM_cumulcoord[i,0] = sum(floorVols[i:,0]*CM_levelcoord[i:,0])/sum(floorVols[i:,0])
            CM_cumulcoord[i,1] = sum(floorVols[i:,0]*CM_levelcoord[i:,1])/sum(floorVols[i:,0])
    
    # LOCATION OF CM FOR ACCIDENTAL TORSIONAL MOMENT GENERATION
    # =========================================================
    CM_levelcoord_acc = np.zeros((len(CM_levelcoord),2))
    CM_levelcoord_acc[:,0] = CM_levelcoord[:,0]+0.05*L
    CM_levelcoord_acc[:,1] = CM_levelcoord[:,1]+0.05*B
            
    
    # ================================        
    # CENTER OF RIGIDITY DETERMINATION
    # ================================
    # Pre-allocation of CR coordinates
    if flex in ('Y','YES','Yes','yES','yes','y'):
        CR_levelcoord = np.zeros((NStr+1,2))
    else:
        CR_levelcoord = np.zeros((NStr,2))
    
    # Pre-allocation of stiffnesses of vertical elements in each floor.
    kz_cols = np.zeros((NStr,PPF))
    kx_cols = np.zeros((NStr,PPF))
    
    # Stiffnessess definition for each column in each floor
    num = 0
    for i in range(NStr):
        for j in range(PPF):
            # X Direction Stiffness
            kx_cols[i,j] = ColDims[num,0]*ColDims[num,1]**3
            kz_cols[i,j] = ColDims[num,1]*ColDims[num,0]**3
            num += 1
            
    # Already have coordinates of columns from CM determination
    # **********************************************************
    
    # Determination of CR coordinates for each floor
    if flex in ('Y','YES','Yes','yES','yes','y'):
        for i in range(NStr):
            CR_levelcoord[i+1,0] = sum(kz_cols[i,:]*ColsXCoord[i,:])/sum(kz_cols[i,:])
            CR_levelcoord[i+1,1] = sum(kx_cols[i,:]*ColsZCoord[i,:])/sum(kx_cols[i,:])
        CR_levelcoord[0,0] = sum(kz_cols[0,:]*ColsXCoord[0,:])/sum(kz_cols[0,:])
        CR_levelcoord[0,1] = sum(kx_cols[0,:]*ColsZCoord[0,:])/sum(kx_cols[0,:])
    else:
        for i in range(NStr):
            CR_levelcoord[i,0] = sum(kz_cols[i,:]*ColsXCoord[i,:])/sum(kz_cols[i,:])
            CR_levelcoord[i,1] = sum(kx_cols[i,:]*ColsZCoord[i,:])/sum(kx_cols[i,:])
    
    return CM_levelcoord, CM_cumulcoord, CM_levelcoord_acc, CR_levelcoord



def detMt(Vx_z,CM_cumulcoord,CR_levelcoord,B,L):
    
    # DETERMINATION OF ACTUAL AND ACCIDENTAL TORSIONAL MOMENT FOR DEFINITION OF ETA
    # ==============================================================================
    import numpy as np
    
    # Actual torsional moment at each story
    # -------------------------------------
    act_Mt = np.zeros((len(Vx_z),2))        # Pre.allocating matrix with [act_Mtx,act_Mtz]
    for i in range(len(Vx_z)):
        # Actual Torsional Moment X - Shear Force X with excentricity in Z.
        act_Mt[i,0] = Vx_z[i,0]*abs(CM_cumulcoord[i,1] - CR_levelcoord[i,1])
        # Actual Torsional Moment Z - Shear Force Z with excentricity in X.
        act_Mt[i,1] = Vx_z[i,1]*abs(CM_cumulcoord[i,0] - CR_levelcoord[i,0])
        
    # Accidental torsional moment at each story.
    # ------------------------------------------
    acc_Mt = np.zeros((len(Vx_z),2))     # Pre-allocating matrix with [acc_Mtx,acc_Mtz]
    for i in range(len(Vx_z)):
        # Accidental torsional moment X - Shear force X with acidental excentricity in Z.
        acc_Mt[i,0] = Vx_z[i,0]*0.05*B
        # Accidental torsional moment Z - Shear force Z with acidental excentricity in X.
        acc_Mt[i,1] = Vx_z[i,1]*0.05*L
        
    return act_Mt, acc_Mt


def eta_factor(ndm,ndf,NbaysX,NbaysZ,XbayL,ZbayL,StoryH,NStr,flex,fpc,Ec,fy,E0,bsteel,
               ColTransfType,XBDims,ZBDims,ColDims,BeamEffFact,ColEffFact,EMs,gamma_conc,g,
               Qsd,Ql,Qlr,lon,lat,Vs30,R,Ie,hx,StrType):
    
    # Torsional amplification multiplier
    # -----------------------------------
    # To determine eta factor for modifications due to torsional effects.
    # A lateral analysis should be made to determine this factor acording to
    # the deflections.
    # eta_actual-torsion is calculated
    # eta_(actual+accidental)-torsion is calculated
    
    import numpy as np
    import openseespy.opensees as ops
    import OPSDefsMOD as OPSDMOD
    import OPSDefsAN as OPSDAN
    import ASCE716Seismic as ASCE716
    
    L = NbaysX*XbayL                  # [m] large side of plan dimensions.
    B = NbaysZ*ZbayL                  # [m] short side of plan dimensions.
    Ss = detSs(lon, lat)              # [g] acceleration parameter.
    S1 = detS1(lon, lat)              # [g] acceleration parameter.
    
    # Pre-generation of eta_factor matrices 
    # -------------------------------------
    # Actual moment only.
    eta_act = np.zeros((len(hx),2))
    # Actual plus accidental.
    eta_actacc = np.zeros((len(hx),2))
    
    # *************************************************************************************
    # LATERAL EQUIVALENT LOAD CASE - X & Z DIRECTIONS -- ACTUAL TORSIONAL MOMENT
    # *************************************************************************************
    dirs = ['x','z']
    for i in range(len(dirs)):
        ops.wipe()
        OPSDMOD.ModelGen(ndm,ndf)
        OPSDMOD.NodeGen(NbaysX,NbaysZ,XbayL,ZbayL,StoryH,NStr,flex)
        OPSDMOD.MastNodeGen(NbaysX,NbaysZ,XbayL,ZbayL,StoryH,NStr,flex)
        OPSDMOD.SPConstGen(NbaysX,NbaysZ,flex)
        OPSDMOD.MPConstGen(NbaysX,NbaysZ,NStr,flex)
        OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
        OPSDMOD.GeomTransGen(ColTransfType)
        OPSDMOD.ElementGen(NbaysX,NbaysZ,XbayL,ZbayL,NStr,StoryH,XBDims,ZBDims,\
                                   ColDims,BeamEffFact,ColEffFact,Ec,fy,EMs)
        [Wx,MassInputMatr] = \
            OPSDMOD.LumpedMassGen(NbaysX,NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr,StoryH,Qsd,flex)
        
        # GRAVITY ANALYSIS EXCECUTED ACCOUNTING FOR LIVE, DEAD AND SUPERDEAD LOADS.
        # ========================================================================
        OPSDMOD.DeadLoadGen(NbaysX,NbaysZ,NStr,XBDims,ZBDims,ColDims,gamma_conc)
        OPSDMOD.SuperDeadLoadGen(NbaysX,NbaysZ,NStr,XbayL,ZbayL,Qsd)
        OPSDMOD.LiveLoadGen(NbaysX,NbaysZ,NStr,XbayL,ZbayL,Ql,Qlr)
        
        # MODAL ANALYSIS FOR DETERMINING FUNDAMENTAL PERIOD AND ITS DIRECTION.
        # ==========================================================================
        (T,Tmaxreal,Mast) = OPSDAN.ModalAnalysis(NStr,B,L,Wx,flex)
            
        # DETERMINATION OF LATERAL LOADS APPLIED TO EACH FLOOR, FxD.
        # =========================================================
        # Using functions from ASCE710Seismic Module.
        (SDS,SD1,Ta,Tforces,Tdrifts,Cu,Tmax,CsF,CsD,VF,VD,FxFX,FxDX,VxFX,VxDX) = \
            ASCE716.ELFP(Ss,S1,Vs30,Wx,R,Ie,hx,StrType,T[0,i])
        
                
        # EQUIVALENT LATERAL FORCE PROCEDURE - FORCES APPLICATION
        # =======================================================
        ops.loadConst('-time',0.0)   # Mantaining the gravity loads constant along the
                                     # ELF application.
        
        OPSDMOD.ELFPForceGen(dirs[i],FxDX,flex)
        
        # EQUIVALENT LATERAL FORCE CASE ANALYSIS - STATIC
        # ================================================
        ops.wipeAnalysis()
        ops.system('BandGeneral')
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.test('NormDispIncr', 1.0e-6, 100)
        ops.algorithm('Newton')
        ops.integrator('LoadControl', 1)
        ops.analysis('Static')
        ops.analyze(1)
        
        # Getting answers from that bi... building.
        # -----------------------------------------
        if dirs[i] in ['x','X']:
            for j in range(len(hx)):
                d_n1 = ops.nodeDisp(int(NbaysX+1+1e4+1e5*(j+1)),1)
                d_n2 = ops.nodeDisp(int((NbaysX+1)*(NbaysZ+1)+1e4+1e5*(j+1)),1)
                d_max = max(d_n1,d_n2)
                d_avg = (d_n1 + d_n2)/2
                eta_act[j,i] = d_max/d_avg
        elif dirs[i] in ['z','Z']:
            for j in range(len(hx)):
                d_n1 = ops.nodeDisp(int(NbaysX+1+1e4+1e5*(j+1)),3)
                d_n2 = ops.nodeDisp(int(1+1e4+1e5*(j+1)),3)
                d_max = max(d_n1,d_n2)
                d_avg = (d_n1 + d_n2)/2
                eta_act[j,i] = d_max/d_avg
    
    # *************************************************************************************
    # LATERAL EQUIVALENT LOAD CASE - X & Z DIRECTIONS -- ACTUAL + ACCIDENTAL TORSIONAL MOMENT
    # *************************************************************************************
    
    (_,_,CM_levelcoord_acc,_) = \
        detCMCR(NbaysX,NbaysZ,NStr,XbayL,ZbayL,StoryH,XBDims,ZBDims,ColDims,0.25,flex)
    
    for i in range(len(dirs)):
        ops.wipe()
        OPSDMOD.ModelGen(ndm,ndf)
        OPSDMOD.NodeGen(NbaysX,NbaysZ,XbayL,ZbayL,StoryH,NStr,flex)
        OPSDMOD.MastNodeGen(NbaysX,NbaysZ,XbayL,ZbayL,StoryH,NStr,flex,CM_levelcoord_acc)
        OPSDMOD.SPConstGen(NbaysX,NbaysZ,flex)
        OPSDMOD.MPConstGen(NbaysX,NbaysZ,NStr,flex)
        OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
        OPSDMOD.GeomTransGen(ColTransfType)
        OPSDMOD.ElementGen(NbaysX,NbaysZ,XbayL,ZbayL,NStr,StoryH,XBDims,ZBDims,\
                                   ColDims,BeamEffFact,ColEffFact,Ec,fy,EMs)
        [Wx,MassInputMatr] = \
            OPSDMOD.LumpedMassGen(NbaysX,NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr,StoryH,Qsd,flex)
        
        # GRAVITY ANALYSIS EXCECUTED ACCOUNTING FOR LIVE, DEAD AND SUPERDEAD LOADS.
        # ========================================================================
        OPSDMOD.DeadLoadGen(NbaysX,NbaysZ,NStr,XBDims,ZBDims,ColDims,gamma_conc)
        OPSDMOD.SuperDeadLoadGen(NbaysX,NbaysZ,NStr,XbayL,ZbayL,Qsd)
        OPSDMOD.LiveLoadGen(NbaysX,NbaysZ,NStr,XbayL,ZbayL,Ql,Qlr)
        
        # MODAL ANALYSIS FOR DETERMINING FUNDAMENTAL PERIOD AND ITS DIRECTION.
        # ==========================================================================
        (T,Tmaxreal,Mast) = OPSDAN.ModalAnalysis(NStr,B,L,Wx,flex)
            
        # DETERMINATION OF LATERAL LOADS APPLIED TO EACH FLOOR, FxD.
        # =========================================================
        # Using functions from ASCE710Seismic Module.
        (SDS,SD1,Ta,Tforces,Tdrifts,Cu,Tmax,CsF,CsD,VF,VD,FxFX,FxDX,VxFX,VxDX) = \
            ASCE716.ELFP(Ss,S1,Vs30,Wx,R,Ie,hx,StrType,T[0,i])
        
                
        # EQUIVALENT LATERAL FORCE PROCEDURE - FORCES APPLICATION
        # =======================================================
        ops.loadConst('-time',0.0)   # Mantaining the gravity loads constant along the
                                     # ELF application.
        
        OPSDMOD.ELFPForceGen(dirs[i],FxDX,flex)
        
        # EQUIVALENT LATERAL FORCE CASE ANALYSIS - STATIC
        # ================================================
        ops.wipeAnalysis()
        ops.system('BandGeneral')
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.test('NormDispIncr', 1.0e-6, 100)
        ops.algorithm('Newton')
        ops.integrator('LoadControl', 1)
        ops.analysis('Static')
        ops.analyze(1)
        
        # Getting answers from that bi... building.
        # -----------------------------------------
        if dirs[i] in ['x','X']:
            for j in range(len(hx)):
                d_n1 = ops.nodeDisp(int(NbaysX+1+1e4+1e5*(j+1)),1)
                d_n2 = ops.nodeDisp(int((NbaysX+1)*(NbaysZ+1)+1e4+1e5*(j+1)),1)
                d_max = max(d_n1,d_n2)
                d_avg = (d_n1 + d_n2)/2
                eta_actacc[j,i] = d_max/d_avg
        elif dirs[i] in ['z','Z']:
            for j in range(len(hx)):
                d_n1 = ops.nodeDisp(int(NbaysX+1+1e4+1e5*(j+1)),3)
                d_n2 = ops.nodeDisp(int(1+1e4+1e5*(j+1)),3)
                d_max = max(d_n1,d_n2)
                d_avg = (d_n1 + d_n2)/2
                eta_actacc[j,i] = d_max/d_avg
    
    # =============================================================
    # Determination of Ax amplification factor for Linear analyses
    # =============================================================
    Ax = np.zeros((len(hx),2))       # pre-definition of Ax matrix.
    for j in range(len(dirs)):
        for i in range(len(hx)):
            ax = (eta_actacc[i,j]/1.20)**2
            if ax <= 3 and ax >= 1:
                Ax[i,j] = ax
            elif ax < 1.0:
                Ax[i,j] = 1
            elif ax > 3.0:
                Ax[i,j] = 3.0
    
    return eta_act, eta_actacc, Ax



def TargetDisp(lon,lat,Vs30,SHL,NStr,BldTypCo,T,W,R,Ie,Cm,
               beta=0.05,TL=8.,Tf=6.0,Npts=500,Vy=0,Te=0):
    # Target displacement calculation according to ASCE4117.
    
    from math import pi
    
    g = 9.81                    # [m/s**2] acceleration of gravity.
    
    # Definition of whether to use T or Te for calculations.
    if Te != 0:
        T = Te
        
    # Determination of modification factors
    Co = detCo(NStr,BldTypCo)
    C1 = detC1(lon,lat,Vs30,SHL,T,W,R,Ie,Cm,beta,TL,Tf,Npts,Vy,T)
    C2 = detC2(lon,lat,Vs30,SHL,T,W,R,Ie,Cm,beta,TL,Tf,Npts,Vy,T)
    
    # Determination of spectral acceleration, Sa
    (Sxs,Sx1) = detSxi(lon,lat,Vs30,SHL)
    (Sa,_,_,_,_,_,_) = detGHRS(Sxs,Sx1,T,beta,TL,Tf,Npts)
    
    # Target Displacement calculation 
    # ================================
    dtarget = Co*C1*C2*Sa*T**2/(4*pi**2)*g

    return dtarget



def mu_maxima(lon,lat,Vs30,SHL,NStr,BldTypCo,T,W,R,Ie,Cm,Delta_y,Delta_Vbmax,alpha_PDelta,
               alpha2,beta=0.05,TL=8.,Tf=6.0,Npts=500,Vy=0,Te=0):
    
    # For this definition, the Te value is the effective fundamental period of
    # the structure calculated as in Eq.. (7-27) of ASCE41-17. For the first
    # approximation of NSP response, the value used for Te is the real dynamicaly
    # calculated period, T; and then calculated again with the obtained Te value.
    
    import numpy as np
    
    # Definition of whether to use T or Te for calculations.
    if Te != 0:
        T = Te
    
    # Determination of Seismic spectral parameters
    (Sxs,Sx1) = detSxi(lon,lat,Vs30,SHL)
    # Lambda factor of Near-Field effect
    if Sx1 >= 0.60 and SHL in ('BSE-2N'):
        lamb = 0.80
    elif Sx1 < 0.60 and SHL in ('BSE-2N'):
        lamb = 0.20
    else:
        lamb = 0
        
    # Target displacement calculation
    delta_t = TargetDisp(lon,lat,Vs30,SHL,NStr,BldTypCo,T,W,R,Ie,Cm,beta,TL,Tf,Npts,Vy,T)
    # Leser between the target displacement and the displacement corresponding
    # to maximum base shear according to idealized capacity curve.
    Delta_d = min(delta_t,Delta_Vbmax)
    h = 1 + 0.15*np.log(T)
    # Effective negative  post-yield slope ratio
    alpha_e = alpha_PDelta + lamb(alpha2 - alpha_PDelta)
    # Determination of the maximum strength ratio, mu_max.
    # ----------------------------------------------------
    mu_max = Delta_d/Delta_y + abs(alpha_e)**(-h)/4
    
    return mu_max



def IFDC(dtg,capacity_points):
    # Generation of a idealized force-displacement curve for NSP of analysis.
    
    import numpy as np
    
    # Previous calculations
    # ---------------------
    # Delta-d displacement according to ASCE4117 got from curve data.
    # Verify if target displacement is array of results from PO analysis, if True
    # the shear force is taken from the result array, if False it is calculated
    # by interpolation.
    if dtg in capacity_points[:,0]:
        indexVdtg = 0
        while capacity_points[indexVdtg,0] != dtg:
            indexVdtg += 1
        V_dtg = capacity_points[indexVdtg,1]
    else:
        indexVdtg = 0
        while capacity_points[indexVdtg,0] < dtg:
            indexVdtg += 1
        # Having the index of the value greater than dtg now the value
        # of V_dtg can be interpolated for the dtg value.
        V_dtg = capacity_points[indexVdtg-1,1] + (dtg - capacity_points[indexVdtg-1,0])*\
            (capacity_points[indexVdtg,1] - capacity_points[indexVdtg-1,1])/\
                (capacity_points[indexVdtg,0] - capacity_points[indexVdtg-1,0])
    
    # Maximum shear force obtained from the pushover analysis.
    maxV = max(capacity_points[:,1])
    # Index of maximum value of shear in  capactiy curve.
    index = 0
    while maxV != capacity_points[index,1]:
        index += 1
    maxV_disp = capacity_points[index,0]
    
    # Definition of point coordinates [Delta_d,V_d]
    # ---------------------------------------------
    if maxV_disp < dtg:
        Delta_d, V_d, ind = maxV_disp, maxV, index+1
    else:
        Delta_d, V_d, ind = dtg, V_dtg, indexVdtg
    
    id_area = 0            # [kN-m] area under idealized the capacity curve.
    real_area = np.trapz(capacity_points[:ind,1],capacity_points[:ind,0])
                           # [kN-m] area under real the capacity curve.
    alpha1 = 0.01          # [adim] slope ratio for second segment of capacity curve.
    
    varalpha = []
    diffarea = []
    
    count = 1
    while id_area < real_area and alpha1 > 0:
        Vi = capacity_points[count,1]         # [kN] shear force in step i
        Deltai = capacity_points[count,0]     # [m] displacement of roof in step i
        V_yi = 0.60**(-1)*Vi                  # [kN] assumed yield shear force
        Kei = Vi/Deltai                       # [kN/m] effective stiffnesss of building at 60%
        Delta_yi = (Deltai/Vi)*V_yi           # [m] calculated yielding deformation at top of building
        alpha1 = (V_d - V_yi)/(Delta_d - Delta_yi)*Kei**(-1)   # [adim] post yield slope ratio
        varalpha.append(alpha1)
        id_area = 0.5*Delta_yi*V_yi + 0.5*(V_yi + V_d)*(Delta_d - Delta_yi)
        diffarea.append(real_area - id_area)
        count += 1
        
    Delta_y = Delta_yi          # [m] definitive yielding deformation of idealized curve
    V_y = V_yi                  # [kN] definitive yielding shear of building.
    Ke = Kei                    # [kN/m] definitive bulding lateral stiffness.
    
    # Calculating an approximation of alpha2 ratio.
    # ---------------------------------------------
    alpha2 = (capacity_points[-1,1] - V_d)/\
        (capacity_points[-1,0] - Delta_d)*Ke**(-1)  # [adim] negative slope ratio
        
    return Delta_y, V_y, Delta_d, V_d, Ke, alpha1, alpha2