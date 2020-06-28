# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 09:38:49 2020

@author: sebas
"""

"""
Script generated with all algorithms used in dimensioning and reinforce
gathered in a single document and optimized (I hope!)
"""
# Importing needed modules
# ========================
import os
from timeit import default_timer as timer
import unitskN_m_C as units
import DesignProcedures as dsgnpcs

startime = timer()
directory = os.getcwd()

# IMPORTANT DATA FOR GENERATING THE MODEL.
# ========================================
# Model Definition Data
# ---------------------
ndm = 3
ndf = 6
# Geometric definitions
# ---------------------
NbaysX = [4,6,8,10]                               # [INTEGER] number of bays in X direction.
NbaysZ = 4                                 # [INTEGER] number of bays in X direction.
XbayL = 5.0                                # [m] length of bays in X direction.
ZbayL = 5.0                                # [m] length of bays in X direction.
StoryH = 3.0                               # [m] story height - uniform for all levels.
NStr = [5,8,10,12,15,18,20]                                 # [INTEGER] number of levels above base.
B = NbaysZ*ZbayL                           # [m] plan width of building. Constant! always < L.
XbeamOff = 0
ZbeamOff = 0
ColOff = 0

# Material definitions
# --------------------
fpc = 35.*units.MPa                        # [MPa] concrete compressive strength (absolut value).
fy = 420.*units.MPa                        # [MPa] steel yield strength in tension.
E0 = 200.*units.GPa                        # [GPa] initial elastic tangent of steel.
etlim = 0.0052                             # [m/m] limit strain for tensile reinforcing steel.
bsteel = 0.01                              # [ratio] strain-hardening ratio.
Ec = (4700*(fpc/units.MPa)**0.5)*units.MPa # [MPa] Concrete Modulus of Elasticity
nuconc = 0.2                               # [(m/m)/(m/m)] Poisson ratio
gamma_conc = 2.4*units.tonf/units.m3       # [tonf/m3] weight per unit volume of concrete.
rec = 0.04*units.m                         # [m] typical concrete cover for beam-column sections.
cover = 0.0625*units.m                     # [m] cover to centroid of longitudinal reinforcement.
transtype = 'other'                        # [STRING] type of transversal reinforcement. 'spiral' or 'other'.

# Data related to soil characteristics
# ------------------------------------
gamma_soil = 18.                # [kN/m3] average unit wight of a SC soil (clayey sand)
Vso = [100.,155.,200.,266.,300.,400.,489.]                               # [m/s] site class definition in means of shear wave velocity 
                                           #       according to ASCE7-16 commentaries of section 11.
nu_soil = 0.25                  # [(m/m)/(m/m)] Poisson ratio for soil.
omega_soil = 0.                 # [rad/s] frequency of soil taken as the one the structure has.
D = 0.0                         # [m] foundation depth.
Re = 0.3                        # [adim] factor indicating the fraction of length considered for kz_xx,yy

# Inertia reduction factor as by ACI318-14.
# -----------------------------------------
BeamEffFact = 0.35
ColEffFact = 0.70
# Other Useful data.
# -------------------
g = 9.807                                  # [m/s2] gravity acceleration in SI units.
ColTransfType = 'PDelta'                   # Could be 'PDelta', 'Linear' or 'Corotational'
N = 5                                      # [int(u)] number of nodes distributed along the element for integration
                                           #          when using elastic elements and a Lobatto integration method.
                                           #          always odd number!

# DISTRIBUTED LOADS TO INCLUDE IN MODAL ANALYSIS AS DYNAMIC WEIGHT (Dead + SuperDead)
# ===================================================================================
Qsd = 7.15                                 # [kN/m2] distributed dead load per unit area over the slab.
Ql = 2.00                                  # [kN/m2] distributed live load per unit area over the slab.
Qlr = 0.70                                 # [kN/m2] distributed live load per unit area over the roof slab.

# SEISMIC DATA
# ------------
Ss = 1.50                                  # [g] ground motion parameter for 0.20s spectral response acc.
S1 = 0.70                                  # [g] ground motion parameter for 1.00s spectral response acc.
R = 8.0                                    # [adim] response modification coefficient.
Omega_o = 3.                               # [adim] overstrengh factor.
Cd = 5.5                                   # [adim] deflection amplification factor.
Ie = 1.00                                  # [adim] importance factor.
StrType = 1                                # [INTEGER] 1 - corresponding to concrete moment resisting frames.
RiskCat = 'II'                             # [STRING] risk category definition according to ASCE7-16.



# ANALYSIS TYPE CONSIDERATIONS
# ===============================
dirs = ['X','Z']                           # Directions of analysis.
DsgnType = 'conv_dsgn'                     # [STRING] could be 'conv_dsgn' or 'SSI_dsgn'



# DESIGN PROCEDURES RECALL FOR COMPLETE BUILDING DESIGN
# ======================================================

if DsgnType in ('conv_dsgn'):
    [ColDims,Colreinf,XBDims,XBreinf,ZBDims,ZBreinf,Di_h,FxD,FxF,T,Ta,Tmax,VxD,VxF,Da,hxplot] = \
        dsgnpcs.RegularDesign(ndm,ndf,NbaysX,NbaysZ,XbayL,ZbayL,StoryH,NStr,B,XbeamOff,ZbeamOff,ColOff,
                   fpc,fy,E0,etlim,bsteel,Ec,nuconc,gamma_conc,rec,cover,transtype,BeamEffFact,
                   ColEffFact,g,ColTransfType,N,Qsd,Ql,Qlr,Ss,S1,R,Cd,Ie,Vso,StrType,RiskCat,
                   dirs,directory)
elif DsgnType in ('SSI_dsgn'):
    [ColDims,Colreinf,XBDims,XBreinf,ZBDims,ZBreinf,DiSSI_h,FxD,FxDSSI,YReaction,\
        FxF,FxFSSI,T,Ta,Tmax,VgorD,VgorF,VxD,VxF,DeltaVF,DeltaVD,DeltaVF,Da,hxplot] = \
        dsgnpcs.SSIDesign(ndm,ndf,NbaysX,NbaysZ,XbayL,ZbayL,StoryH,NStr,B,XbeamOff,ZbeamOff,ColOff,
                   fpc,fy,E0,etlim,bsteel,Ec,nuconc,gamma_conc,rec,cover,transtype,BeamEffFact,
                   ColEffFact,g,ColTransfType,N,gamma_soil,nu_soil,omega_soil,D,Re,Qsd,Ql,Qlr,
                   Ss,S1,R,Cd,Omega_o,Ie,Vso,StrType,RiskCat,dirs,directory)
            
elapsed = timer() - startime
print(f"Total elapsed time is {round(elapsed,0)} [s] - {round(elapsed/60,2)} [m]")