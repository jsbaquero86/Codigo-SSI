# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 17:08:35 2020

@author: sebas
"""

# Importing needed modules
# ========================
import os
from timeit import default_timer as timer
import unitskN_m_C as units
import matplotlib.pyplot as plt
import pickle
# import ASCE716Seismic as ASCE716
import openseespy.opensees as ops
import OPSDefsMOD as OPSDMOD
import OPSDefsAN as OPSDAN
import AnalysisProceduresASCE41 as anprcdrs
# import ASCE4117

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
NbaysX = [6]                               # [INTEGER] number of bays in X direction.
NbaysZ = 4                                 # [INTEGER] number of bays in X direction.
XbayL = 5.0                                # [m] length of bays in X direction.
ZbayL = 5.0                                # [m] length of bays in X direction.
StoryH = 3.0                               # [m] story height - uniform for all levels.
NStr = [12]                                 # [INTEGER] number of levels above base.
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
Vso = [266.]                               # [m/s] site class definition in means of shear wave velocity 
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
# The following acceleration parameters are obtained in the inner procces of analyses.
# (Sxs,Sx1) = detSxi(lon,lat,Vs30,SHL)       # [g] Sxs, ground motion parameter for 0.20s spectral response acc.
#                                            # [g] Sx1, ground motion parameter for 1.00s spectral response acc.
lon = -120                                 # [degrees] longitud - geographical location.
lat = 35                                 # [degrees] latitude - geographical location.
SHL = ''
R = 8.0                                    # [adim] response modification coefficient.
Omega_o = 3.                               # [adim] overstrengh factor.
Cd = 5.5                                   # [adim] deflection amplification factor.
Ie = 1.00                                  # [adim] importance factor.
StrType = 1                                # [INTEGER] 1 - corresponding to concrete moment resisting frames.
RiskCat = 'II'                             # [STRING] risk category definition according to ASCE7-16.
# SiteClass = ASCE716.detSiteClass(Vso[0])   # [STRING] site class determination according to ASCE7.


# ANALYSIS TYPE CONSIDERATIONS
# ===============================
dirs = ['X']                           # Directions of analysis.
DsgnType = 'conv_dsgn'                     # [STRING] could be 'conv_dsgn' or 'SSI_dsgn'
EMs = 'plastic-3'                            # [STRING] determines the element's mechanic
                                           #          characteristics: 'elastic'
                                           #                           'plastic-1'
                                           #                           'plastic-2'
                                           #                           'plastic-3'
BldTypCo = 'Sh-T'                          # [STRING] building type classification to obtain Co.
BldTypCm = 'C-MF'                          # [STRING] building type classification to obtain Co.

# PROCEDURES RECALL FOR ANALYSIS OF BUILDING/S
# =============================================

(results,dtg,tgfac,dispana) = anprcdrs.NSProcedure(ndm,ndf,NbaysX,NbaysZ,NStr,XbayL,ZbayL,StoryH,lon,lat,SHL,Vso,
                gamma_soil,nu_soil,Re,fpc,Ec,gamma_conc,fy,E0,bsteel,
                ColTransfType,BeamEffFact,ColEffFact,g,Qsd,Ql,Qlr,
                EMs,R,Ie,StrType,BldTypCo,BldTypCm,DsgnType,dirs,directory)



#%%


counter = 0
total = len(NbaysX)*len(NStr)*len(Vso)

L = NbaysX[0]*XbayL                        # [m] plan length of building. Always > B.

workpath = directory + '\\RegularDesign\\' + str(NbaysX[0])+'BayX'+str(NbaysZ)+\
    'BayZ'+str(NStr[0])+'FLRS'+str(Vso[0])+'.pickle'

# workpath = directory + '/RegularDesign/' +  str(NbaysX[0]) + 'BayX' + str(NbaysZ) + \
#     'BayZ'+str(NStr[0])+'FLRS'+str(Vso[0])+'.pickle'
    
with open(workpath,'rb') as f:
    ColDims,Colreinf,XBDims,XBreinf,ZBDims,ZBreinf,Di_h,FxD,FxF,T,Ta,Tmax,VxD,VxF,Da,hxplot = pickle.load(f)
    
# =================================================================================
ops.wipe()
# Model Generation
# ================
OPSDMOD.ModelGen(ndm,ndf)
# Basic nodes generation
# ======================
OPSDMOD.NodeGen(NbaysX[0],NbaysZ,XbayL,ZbayL,StoryH,NStr[0],flex)
# Generation of master nodes
# ===========================
OPSDMOD.MastNodeGen(NbaysX[0],NbaysZ,XbayL,ZbayL,StoryH,NStr[0],flex)
# Single Point constraint of base nodes - fixity of structure.
# =============================================================
OPSDMOD.SPConstGen(NbaysX[0],NbaysZ,flex)
# Multi-Point constraint - Rigid Diaphragm assignment.
# =====================================================
OPSDMOD.MPConstGen(NbaysX[0],NbaysZ,NStr[0],flex)
# Material generation for sections and elements.
# ==============================================
# Only oncrete and reinforcement steel.
OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,b)


# # ****************************************************************************
# # ****************************************************************************
# # ELASTIC ELEMENTS
# # ==================
# # Generation of elastic element sections 
# # ---------------------------------------
# OPSDMOD.ElConcSectGen(NbaysX[0],NbaysZ,NStr[0],Ec,nuconc,XBDims,ZBDims,ColDims,BeamEffFact,ColEffFact)
# # Generation of beam section integrators.
# # ---------------------------------------
# OPSDMOD.LobattoIntegratorGen(NbaysX[0],NbaysZ,NStr[0],N)

# if flex in ('Y','YES','Yes','yES','yes','y'):
#     # Interface elements generation for foundation flexibility considerations.
#     # =========================================================================
#     # Materials generation: stiffness constants accounting for soil flexibility.
#     # ---------------------------------------------------------------------------
#     OPSDMOD.FoundFlexMaterials(NbaysX[0],NbaysZ,XbayL,ZbayL,Ss,Vso[0],gamma_soil,nu_soil,B,L,Re,D,omega_soil,'lat')
#     # Zero-Length elements creation for connecting base nodes.
#     OPSDMOD.FoundFlexZLElements(NbaysX[0],NbaysZ,XbayL,ZbayL,B,L,Re)



# ****************************************************************************
# ****************************************************************************
# PLASTIC ELEMENTS
# ==================
# Generation of PLASTIC element sections 
# ---------------------------------------
OPSDMOD.PlConcSectGen(NbaysX[0],NbaysZ,Ec,nuconc,NStr[0],XBDims,ZBDims,ColDims,BeamEffFact,\
                    ColEffFact,cover,XBreinf,ZBreinf,Colreinf)
# Generation of beam-column section integrators.
# ---------------------------------------
OPSDMOD.RadauIntegratorGen(NbaysX[0],NbaysZ,NStr[0],XbayL,ZbayL,StoryH,fy)
if flex in ('Y','YES','Yes','yES','yes','y'):
    # Interface elements generation for foundation flexibility considerations.
    # =========================================================================
    # Materials generation: stiffness constants accounting for soil flexibility.
    # ---------------------------------------------------------------------------
    OPSDMOD.FoundFlexMaterials(NbaysX[0],NbaysZ,XbayL,ZbayL,Ss,Vso[0],gamma_soil,nu_soil,B,L,Re,D,omega_soil,'lat')
    # Zero-Length elements creation for connecting base nodes.
    OPSDMOD.FoundFlexZLElements(NbaysX[0],NbaysZ,XbayL,ZbayL,B,L,Re)

# ****************************************************************************
# ****************************************************************************
# Generation of Geometric Transformation objects for elements
# ===========================================================
# OPSDMOD.GeomTransGen(ColTransfType,XBD=XBDims[0,:],ZBD=ZBDims[0,:],ColD=ColDims[0,:])    # Rigid offsets on both columns and beams ends.
OPSDMOD.GeomTransGen(ColTransfType,XBD=XBDims[0,:],ZBD=ZBDims[0,:],ColD=[0,0])          # Offset on columns only.
# OPSDMOD.GeomTransGen(ColTransfType,XBD=[0,0],ZBD=[0,0],ColD=ColDims[0,:])               # Offset on beams only.
# OPSDMOD.GeomTransGen(ColTransfType,XBD=[0,0],ZBD=[0,0],ColD=[0,0])    # None rigid offsets..

# Element object generation, Columns and Beams in X and Z global directions.
# --------------------------------------------------------------------------
OPSDMOD.forceBeamColumnGen(NbaysX[0],NbaysZ,NStr[0],XbayL,ZbayL,XBDims,ZBDims,ColDims,gamma_conc=0,g=g,Qsd=0)
# OPSDMOD.forceBeamColumnGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,XBDims,ZBDims,ColDims,gamma_conc,g,Qsd)
# Seems like massdensity of columns are not considered from elements different 
# from horizontal ones to construct the consistent mass matrix. The Periods
# obtained using this mass matrix construction option, I consider, not adequate.
# Better to use my function to define lumped Masses from elemnts and loads below. (JSBM)
# Lumped mass determination by "manual" (element-oriented) calculation of masses.
# -------------------------------------------------------------------------------
# (Not to be used while using forceBeamColumnGen with non-zero values of gamma_conc and Qsd)
[Wx,MassInputMatr] = \
    OPSDMOD.LumpedMassGen(NbaysX[0],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[0],StoryH,Qsd,flex)
W = sum(Wx)
# ****************************************************************************
# ****************************************************************************
# GRAVITY ANALYSIS EXCECUTED ACCOUNING FOR LIVE, DEAD AND SUPERDEAD LOADS.
# ========================================================================
OPSDMOD.DeadLoadGen(NbaysX[0],NbaysZ,NStr[0],XBDims,ZBDims,ColDims,gamma_conc)
OPSDMOD.SuperDeadLoadGen(NbaysX[0],NbaysZ,NStr[0],XbayL,ZbayL,Qsd)
OPSDMOD.LiveLoadGen(NbaysX[0],NbaysZ,NStr[0],XbayL,ZbayL,Ql,Qlr)

# GRAVITY-LOADS-CASE ANALYSIS.
# ============================
ops.system('ProfileSPD')
ops.constraints('Transformation')
ops.numberer('RCM')
ops.test('NormDispIncr',1.0e-4,100)
ops.algorithm('KrylovNewton')
ops.integrator('LoadControl', 1)
ops.analysis('Static')
ops.analyze(1)
# Vertical reactions Calculation for verification.
# ------------------------------------------------
ops.reactions()
YReact = 0
for i in range((NbaysX[0]+1)*(NbaysZ+1)):
    if flex in ('Y','YES','Yes','yES','yes','y'):
        YReact += ops.nodeReaction(99000+i+1,2)
    else:
        YReact += ops.nodeReaction(i+1,2)
print('>>> Gravity loads applied...')

# ==========================================================================
# MODAL ANALYSIS FOR DETERMINING FUNDAMENTAL PERIOD AND ITS DIRECTION.
# ==========================================================================
(Tver,Tmaxver,Mast) = OPSDAN.ModalAnalysis(NStr[0],B,L,Wx,flex)

import ASCE716SSI as SSI
if flex in ('Y','YES','Yes','yES','yes','y'):
    # [m][LIST] with level height starting from first level or from base if foundation flexibility is included.
    hx = [0.0001]
else:
    hx = []

for i in range(NStr[0]):
    hx.append((i+1)*StoryH)
TgSSI710 = SSI.detTgor(Ss,Vso[0],gamma_soil,nu_soil,B,L,hx[-1],Mast,T[0,0],'X',D,omega_soil)

#%% ==================
# PUSHOVER ANALYSIS
# ==================
time_o = timer()

BldTypCo = 'Sh-T'
BldTypCm = 'C-MF'

dtarget = OPSDAN.POTargetDisp(Ss,S1,Tver[0,0],W,Vso[0],NStr[0],BldTypCo,BldTypCm)


# Aplication of forces to the model.
# ----------------------------------
ops.loadConst('-time',0.0)
MdeShape = OPSDMOD.POForceGen(NStr[0],FxF[:,0],dirs[0],flex)
# Execution of the analysis and output results.
# ---------------------------------------------
(results,dtg) = \
    OPSDAN.POAnalysisProc(Ss,S1,Vso[0],NbaysX[0],NbaysZ,NStr[0],StoryH,\
                        BldTypCo,BldTypCm,Tver[0,0],W,dirs[0])  # [Disp,Force]

ttime = timer()- time_o
print(f'Elapsed Time {round(ttime/60,2)} [m]')

plt.figure()
plt.plot(results[:,0],results[:,1])
plt.grid()

#%% ======================
# TIME HISTORY ANALYSIS
# =======================
# Including Rayleigh damping to the model
# ---------------------------------------
ops.loadConst('-time',0.0)
lambi = 0.05        # [ratio] damping ratio on the mode i
lambj = 0.05        # [ratio] damping ratio on mode j.
# OPSD.RayleighDamping(T[0,0],T[0,1],lambi,lambj)
OPSDAN.RayleighDamping(T[0,0],0.1*T[0,0],lambi,lambj)

# Timer to measure the process time lapsing.
# ------------------------------------------
startime = timer()

# Transient analysis
# ------------------
dt_rec = 0.01       # [s] record time interval according to accelerogram records.
NPts = 4896         # [u] number of discrete acceleration points.
FilePath = 'C:\\Users\\sebas\\Documents\\MEGAsync\\DOCTORADO_UPC\\INVESTIGACION\\'+\
    'INTERACCION SUELO-ESTRUCTURA\\TOOLBOX\\SP_17-14_ACI318-14\\PY_TOOLBX\\'+\
        'CompleteBuildingDesign\\Manta2016North.txt'

(THdisp,t_plot,floors) = OPSDAN.\
    THAnalysisProc(FilePath,dt_rec,NStr[0],dirs[0],gammaNw=1/2,betaNw=1/4,dt_an=dt_rec,fract=0.5)

totaltime = timer() - startime
print(f'Total elapsed time for T-H analysis is {round(totaltime/60,2)} [min].') 

plt.figure()
plt.plot(t_plot,THdisp['Str-5'])
plt.grid()

#%%
import ReadRecord
import matplotlib.pyplot as plt


# Set the gravity loads to be constant & reset the time in the domain
ops.loadConst('-time', 0.0)

# ----------------------------------------------------
# End of Model Generation & Initial Gravity Analysis
# ----------------------------------------------------

# Define nodal mass in terms of axial load on columns
g = 9.81


# Set some parameters
record = 'elCentro'

# Permform the conversion from SMD record to OpenSees record
dt, nPts = ReadRecord.ReadRecord(record+'.at2', record+'.dat')

# Set time series to be passed to uniform excitation
ops.timeSeries('Path', 7, '-filePath', record+'.dat', '-dt', dt, '-factor', g)

# Create UniformExcitation load pattern
#                         tag dir 
ops.pattern('UniformExcitation',  7,   1,  '-accel', 7)

# set the rayleigh damping factors for nodes & elements
ops.rayleigh(0.4274, 0.005827, 0.0, 0.0)

# Delete the old analysis and all it's component objects
ops.wipeAnalysis()

# Create the system of equation, a banded general storage scheme
ops.system('BandGeneral')

# Create the constraint handler, a plain handler as homogeneous boundary
ops.constraints('Transformation')

# Create the convergence test, the norm of the residual with a tolerance of 
# 1e-12 and a max number of iterations of 10
ops.test('NormDispIncr', 1.0e-4,  100 )

# Create the solution algorithm, a Newton-Raphson algorithm
ops.algorithm('Newton')

# Create the DOF numberer, the reverse Cuthill-McKee algorithm
ops.numberer('RCM')

# Create the integration scheme, the Newmark with alpha =0.5 and beta =.25
ops.integrator('Newmark',  0.5,  0.25 )

# Create the analysis object
ops.analysis('Transient')

# Perform an eigenvalue analysis
numEigen = 2
eigenValues = ops.eigen(numEigen)
print("eigen values at start of transient:",eigenValues)

# set some variables
tFinal = nPts*dt
tCurrent = ops.getTime()
ok = 0

time = [tCurrent]
u3 = [0.0]

# Perform the transient analysis
while ok == 0 and tCurrent < tFinal:
    
    ok = ops.analyze(1, .01)
    
    # if the analysis fails try initial tangent iteration
    if ok != 0:
        print("regular newton failed .. lets try an initail stiffness for this step")
        ops.test('NormDispIncr', 1.0e-12,  100, 0)
        ops.algorithm('ModifiedNewton', '-initial')
        ok =ops.analyze( 1, .01)
        if ok == 0:
            print("that worked .. back to regular newton")
        ops.test('NormDispIncr', 1.0e-12,  10 )
        ops.algorithm('Newton')
    
    tCurrent = ops.getTime()
    
    print(f'{tCurrent} [s] analyzed of {tFinal}')

    time.append(tCurrent)
    u3.append(ops.nodeDisp(2999,1))



# Perform an eigenvalue analysis
eigenValues = ops.eigen(numEigen)
print("eigen values at end of transient:",eigenValues)

results = open('results.out','a+')

if ok == 0:
    results.write('PASSED : RCFrameEarthquake.py\n');
    print("Passed!")
else:
    results.write('FAILED : RCFrameEarthquake.py\n');
    print("Failed!")

results.close()

plt.plot(time, u3)
plt.ylabel('Horizontal Displacement of node 3 (in)')
plt.xlabel('Time (s)')

plt.show()



print("==========================")

#%% TIME HISTORY ANALYSIS
# ==========================
import ReadRecord
import numpy as np

g = 9.81     # [m/s**2] acceleration of gravity in SI units.

# Establishing recorders for the output responses.
# -------------------------------------------------
ops.recorder('Node','-file','RoofDispX_Z.out','-time','-node',5999,'-dof',*[1,3],'disp')

# Including Rayleigh damping to the model
# ---------------------------------------
ops.loadConst('-time',0.0)

lambi = 0.05        # [ratio] damping ratio on the mode i
lambj = 0.05        # [ratio] damping ratio on mode j.
OPSDAN.RayleighDamping(T[0,0],T[0,1],lambi,lambj)

# Load record to analyse.
# -----------------------
# record = 'elCentro'
# # Transform the record to a readable format.
# # -------------------------------------------
# dt_record, NPts = ReadRecord.ReadRecord(record+'.at2', record+'.dat')   

record = np.loadtxt('Manta2016North.txt',delimiter='\n')
# record = record.tolist()
NPts = len(record)
dt_record = 0.01

# Set TimeSeries to be passed to uniform excitation
# --------------------------------------------------
# ops.timeSeries('Path',7,'-filePath',record+'.dat','-dt',dt_record,'-factor',g)
ops.timeSeries('Path',7,'-values',*record,'-dt',dt_record,'-factor',g)

# Creating a uniform excitation load pattern
# ------------------------------------------
ops.pattern('UniformExcitation',7,1,'-accel',7)

# ANALYSIS PARAMETERS
# ===================

# First, delete previously defined parameters
# -------------------------------------------
ops.wipeAnalysis()

# System of equations.
# --------------------
ops.system('BandGeneral')

# Create the constraints handler.
# -------------------------------
ops.constraints('Transformation')

# Convergency test
# -----------------
ops.test('NormDispIncr',1e-6,100)

# Solution algorithm
# ------------------
ops.algorithm('Newton')

# DOFs numberer.
# --------------
ops.numberer('RCM')

# Integration scheme for transient analysis.
# -----------------------------------------
ops.integrator('Newmark',0.5,0.25)

# Analysis object
# ----------------
ops.analysis('Transient')

# Perform the transient analysis
# ------------------------------
print('-------------------')
print('TH-analysis started!')
print(' ')
for i in range(NPts):
    ops.analyze(1,dt_record)
    recordtime = ops.getTime()
    print(f'Time-History analysis executing, {round(i/NPts*100,2)} % - {round(recordtime,2)} [s] of ground motion record...')
print(f'Time-History analysis executing, 100 % - {round(recordtime,2)} [s] of ground motion record...')

dispXroof = np.loadtxt('RoofDispX_Z.out')

plt.figure()
plt.plot(dispXroof[:,0],dispXroof[:,1])
plt.grid()