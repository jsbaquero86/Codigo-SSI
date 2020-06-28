# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 10:44:30 2020

@author: sebas
"""

"""
These definitions are written according to the procedures recommended by ASCE41
for the assessment of buildings. Both  Nonlinear static and dynamic analysis
procedures are defined as NSP and NDP, respectively.
"""

"""
    In order to determine which analysis procedure is suitable for use in the 
    assessment process, the following must be considered:
         - Component gravity loads depends on the type of analysis, wheter it
           is linear or nonlinear analysis --> Section 7.2.2.
         - Always the builing shall be modeled as 3D structure, unless it has
           rigid diaphragms and/or no torsional effects exceed limits of 7.2.3.2
           or are accounted appropriately.
         - Torsion shall be considered according to 7.2.3.2.
    
    The definitions for each analysis procedure are as follows.
"""


def NSProcedure(ndm,ndf,NbaysX,NbaysZ,NStr,XbayL,ZbayL,StoryH,lon,lat,SHL,Vso,
                gamma_soil,nu_soil,Re,fpc,Ec,gamma_conc,fy,E0,bsteel,
                ColTransfType,BeamEffFact,ColEffFact,g,Qsd,Ql,Qlr,
                EMs,R,Ie,StrType,BldTypCo,BldTypCm,DsgnType,dirs,directory):
    
    # The Nonlinear Static Procdure is executed in order to get responses
    # of a defined building.
    
    import openseespy.opensees as ops
    import OPSDefsMOD as OPSDMOD
    import OPSDefsAN as OPSDAN
    from timeit import default_timer as timer
    import ASCE716Seismic as ASCE716
    import ASCE4117
    import matplotlib.pyplot as plt
    import pickle
    import numpy as np
    
    for mm in range(len(NbaysX)):
        for nn in range(len(NStr)):
            for oo in range(len(Vso)):
                
                Teff = np.zeros((len(dirs)))     # [s][LIST] effective lateral period of building.
                
                for ii in range(len(dirs)):
                    
                    time_o = timer()
                    
                    # Some previous definitions
                    # ----------------------------
                    if DsgnType in ('conv_dsgn'):
                        flex = 'no'
                    elif DsgnType in ('ssi_dsgn'):
                        flex = 'yes'
                    
                    # SiteClass = ASCE716.detSiteClass(Vso[oo])
                    
                    # Unpicklin' some stored parameters from design process.
                    # ------------------------------------------------------
                    workpath = directory + '\\RegularDesign\\' + str(NbaysX[mm])+'BayX'+str(NbaysZ)+\
                        'BayZ'+str(NStr[nn])+'FLRS'+str(Vso[oo])+'.pickle'
                    
                    with open(workpath,'rb') as f:
                        ColDims,Colreinf,XBDims,XBreinf,ZBDims,ZBreinf,_,_,_,_,_,_,_,_,_,_ = pickle.load(f)
                    
                    # Calculation of height vector
                    # ------------------------------
                    if flex in ('Y','YES','Yes','yES','yes','y'):
                        # [m][LIST] with level height starting from first level or from base if foundation flexibility is included.
                        hx = [0.0001]
                    else:
                        hx = []
                    
                    for i in range(NStr[nn]):
                        hx.append((i+1)*StoryH)
                        
                    # Plan Dimensions of Building
                    B = NbaysZ*ZbayL                  # [m] short side of building plan.
                    L = NbaysX[mm]*XbayL              # [m] long side of building plan.
                    
                    # Determination of MCEr spectral acceleration parameters
                    # ------------------------------------------------------
                    (Sxs,Sx1) = ASCE4117.detSxi(lon,lat,Vso[oo],SHL)
                        
                    # MODELING OF THE STRUCTURE USING OPENSEESPY
                    # ===========================================
                    ops.wipe()
                    OPSDMOD.ModelGen(ndm,ndf)
                    OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex)
                    OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex,coords=0)
                    OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,flex)
                    OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],flex)
                    OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                    # OPSDMOD.GeomTransGen(ColTransfType,XBD=[min(XBDims[:,0]),min(XBDims[:,1])],\
                    #                      ZBD=[min(ZBDims[:,0]),min(ZBDims[:,1])],\
                    #                          ColD=[min(ColDims[:,0]),min(ColDims[:,1])])
                    OPSDMOD.GeomTransGen(ColTransfType,ColD=[min(ColDims[:,0]),min(ColDims[:,1])])
                    # OPSDMOD.GeomTransGen(ColTransfType)
                        
                    if flex in ('Y','YES','Yes','yES','yes','y'):
                        # Interface elements generation for foundation flexibility considerations.
                        # =========================================================================
                        # Materials generation: stiffness constants accounting for soil flexibility.
                        # ---------------------------------------------------------------------------
                        OPSDMOD.FoundFlexMaterials(NbaysX[mm],NbaysZ,XbayL,ZbayL,Sxs,Vso[oo],gamma_soil,nu_soil,B,L,Re,\
                                                   D=0,omega_soil=0,analtype='lat')
                        # Zero-Length elements creation for connecting base nodes.
                        OPSDMOD.FoundFlexZLElements(NbaysX[mm],NbaysZ,XbayL,ZbayL,B,L,Re)
                        
                    OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                               ColDims,BeamEffFact,ColEffFact,Ec,fy,EMs,\
                                                   XBreinf,ZBreinf,Colreinf,N=5,rec=0.0654,nuconc=0.2,dbar=0.025)
                    [Wx,MassInputMatr] = \
                        OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,flex)
                    W = sum(Wx)                # [kN] total weight of the building
                    
                    # GRAVITY LOADS APPLIED TO MODEL ACCORDINGO TO ASCE4117
                    # ======================================================
                    # According to ASCE4117 Section 7.2.2, equation (7-3), the combination
                    # of gravitational loads mus be as follows:
                        #       Qg = Qd + 0.25*Ql + Qs      (7-3)
                    OPSDMOD.DeadLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XBDims,ZBDims,ColDims,gamma_conc)
                    OPSDMOD.SuperDeadLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,Qsd)
                    OPSDMOD.LiveLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,0.25*Ql,0.25*Qlr)
                    
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
                    for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                        if flex in ('Y','YES','Yes','yES','yes','y'):
                            YReact += ops.nodeReaction(int(i+1),2)
                        else:
                            YReact += ops.nodeReaction(int(i+1+1e4),2)
                    
                    # ==========================================================================
                    # MODAL ANALYSIS FOR DETERMINING FUNDAMENTAL PERIOD AND ITS DIRECTION.
                    # ==========================================================================
                    (T,Tmaxver,Mast) = OPSDAN.ModalAnalysis(NStr[nn],B,L,Wx,flex)
                    print(T)
                    # ==================
                    # PUSHOVER ANALYSIS
                    # ==================
                    # Lateral Force used for analysis
                    # --------------------------------
                    # Using functions from ASCE716Seismic Module.
                    (_,_,_,_,_,_,_,_,_,_,_,FxF,_,_,_) = \
                        ASCE716.ELFP(Sxs,Sx1,Vso[oo],Wx,R,Ie,hx,StrType,T[0,ii])
                        
                    # Aplication of forces to the model.
                    # ----------------------------------
                    ops.loadConst('-time',0.0)
                    MdeShape = OPSDMOD.POForceGen(NStr[nn],FxF,dirs[ii],flex)
                    
                    # =============================================
                    # First Execution of the analysis and output results.
                    # =============================================
                    (results1,dtg1,tgfactor1) = \
                        OPSDAN.POAnalysisProc(lon,lat,Vso[oo],SHL,NbaysX[mm],NbaysZ,\
                                              NStr[nn],StoryH,R,Ie,BldTypCo,BldTypCm,\
                                                  T[0,ii],W,dirs[ii],flex,\
                                                      DispIncr=0,beta=0.05,TL=8.,Tf=6.0,Npts=500,Vy=0,Te=0)  # [Disp,Force]
                    
                    
                    # =====================================================
                    # Determining the first approximation of values of Vy 
                    # and calculation of effective fundamental period for NSP
                    # =====================================================
                    (Delta_y, V_y, Delta_d, V_d, Ke, alpha1, alpha2) = \
                        ASCE4117.IFDC(dtg1,results1)
                    
                    Ki = Mast[ii]*4*np.pi**2/T[0,ii]**2       # [kN/m] elastic lateral stiffness of the building.
                    Teff[ii] = T[0,ii]*(Ki/Ke)**0.5           # [s] effctive fundamental period of building.
                    
                    
                    # =============================================
                    # Second Execution of the analysis and output results.
                    # =============================================
                    ops.wipe()
                    OPSDMOD.ModelGen(ndm,ndf)
                    OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex)
                    OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex,coords=0)
                    OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,flex)
                    OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],flex)
                    OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                    # OPSDMOD.GeomTransGen(ColTransfType,XBD=[min(XBDims[:,0]),min(XBDims[:,1])],\
                    #                      ZBD=[min(ZBDims[:,0]),min(ZBDims[:,1])],\
                    #                          ColD=[min(ColDims[:,0]),min(ColDims[:,1])])
                    OPSDMOD.GeomTransGen(ColTransfType,ColD=[min(ColDims[:,0]),min(ColDims[:,1])])
                    # OPSDMOD.GeomTransGen(ColTransfType)
                        
                    if flex in ('Y','YES','Yes','yES','yes','y'):
                        # Interface elements generation for foundation flexibility considerations.
                        # =========================================================================
                        # Materials generation: stiffness constants accounting for soil flexibility.
                        # ---------------------------------------------------------------------------
                        OPSDMOD.FoundFlexMaterials(NbaysX[mm],NbaysZ,XbayL,ZbayL,Sxs,Vso[oo],gamma_soil,nu_soil,B,L,Re,\
                                                   D=0,omega_soil=0,analtype='lat')
                        # Zero-Length elements creation for connecting base nodes.
                        OPSDMOD.FoundFlexZLElements(NbaysX[mm],NbaysZ,XbayL,ZbayL,B,L,Re)
                        
                    OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                               ColDims,BeamEffFact,ColEffFact,Ec,fy,EMs,\
                                                   XBreinf,ZBreinf,Colreinf,N=5,rec=0.0654,nuconc=0.2,dbar=0.025)
                    [Wx,MassInputMatr] = \
                        OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,flex)
                    W = sum(Wx)                # [kN] total weight of the building
                    # GRAVITY LOADS APPLIED TO MODEL ACCORDINGO TO ASCE4117
                    # ======================================================
                    OPSDMOD.DeadLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XBDims,ZBDims,ColDims,gamma_conc)
                    OPSDMOD.SuperDeadLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,Qsd)
                    OPSDMOD.LiveLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,0.25*Ql,0.25*Qlr)
                    
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
                    # ==================
                    # PUSHOVER ANALYSIS
                    # ==================
                    # Lateral Force used for analysis
                    # --------------------------------
                    # Using functions from ASCE716Seismic Module.
                    (_,_,_,_,_,_,_,_,_,_,_,FxF,_,_,_) = \
                        ASCE716.ELFP(Sxs,Sx1,Vso[oo],Wx,R,Ie,hx,StrType,Teff[ii])
                        
                    # Aplication of forces to the model.
                    # ----------------------------------
                    ops.loadConst('-time',0.0)
                    MdeShape = OPSDMOD.POForceGen(NStr[nn],FxF,dirs[ii],flex)
                    
                    # =============================================
                    # First Execution of the analysis and output results.
                    # =============================================
                    (results2,dtg2,tgfactor2) = \
                        OPSDAN.POAnalysisProc(lon,lat,Vso[oo],SHL,NbaysX[mm],NbaysZ,\
                                              NStr[nn],StoryH,R,Ie,BldTypCo,BldTypCm,\
                                                  T[0,ii],W,dirs[ii],flex,\
                                                      DispIncr=0,beta=0.05,TL=8.,Tf=6.0,Npts=500,Vy=0,Te=Teff[ii])  # [Disp,Force].
                    
                    
                    # =====================================================
                    # Determining the "exact"  values of Vy 
                    # and calculation of effective fundamental period for NSP
                    # =====================================================
                    (Delta_y, V_y, Delta_d, V_d, Ke, alpha1, alpha2) = \
                        ASCE4117.IFDC(dtg1,results2)
                    
                    Ki = Mast[ii]*4*np.pi**2/T[0,ii]**2       # [kN/m] elastic lateral stiffness of the building.
                    Teff[ii] = T[0,ii]*(Ki/Ke)**0.5           # [s] effctive fundamental period of building.
                    
                    
                    
                    ttime = timer()- time_o
                    print(f'Elapsed Time {round(ttime/60,2)} [m]')
                    
                    plt.figure()
                    plt.plot(results2[:,0],results2[:,1])
                    plt.grid()
                    
                    print('The mode shape is:')
                    print(MdeShape)
                    
    return (results1,results2), (dtg1,dtg2), (tgfactor1,tgfactor2), (tgfactor1*dtg1,tgfactor2*dtg2)