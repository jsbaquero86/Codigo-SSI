# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 11:37:12 2020

@author: sebas
"""


def RegularDesign(ndm,ndf,NbaysX,NbaysZ,XbayL,ZbayL,StoryH,NStr,B,XbeamOff,ZbeamOff,ColOff,
                   fpc,fy,E0,etlim,bsteel,Ec,nuconc,gamma_conc,rec,cover,transtype,BeamEffFact,
                   ColEffFact,g,ColTransfType,N,Qsd,Ql,Qlr,Ss,S1,R,Cd,Ie,Vso,StrType,RiskCat,
                   dirs,directory,flex='No'):
    
    import numpy as np
    import matplotlib.pyplot as plt
    import pickle
    import DLBDefs
    import ASCE716Seismic as ASCE716
    import openseespy.opensees as ops
    import OPSDefsMOD as OPSDMOD
    import OPSDefsAN as OPSDAN

    total = len(NbaysX)*len(NStr)*len(Vso)      # [u.] total number of structures to be designed.
    counter = 0                                 # [INTEGER] counter for the evolution of design, one for each 
                                                #           designed building.
    for mm in range(len(NbaysX)):
        for nn in range(len(NStr)):
            for oo in range(len(Vso)):
                
                L = NbaysX[mm]*XbayL               # [m] plan length of building. L >= B always!
                
                if flex in ('Y','YES','Yes','yES','yes','y'):
                    # [m][LIST] with level height starting from first level or from base if foundation flexibility is included.
                    hx = [0.0001]
                else:
                    hx = []
                
                for i in range(NStr[nn]):
                    hx.append((i+1)*StoryH)
                    
                # Section Properties for element (as initial setup)
                # -------------------------------------------------
                XBDims = np.zeros((NbaysX[mm]*(NbaysZ + 1)*NStr[nn],2))        # for [b,h] dimensions. b is normal to local y-axis.
                ZBDims = np.zeros((NbaysZ*(NbaysX[mm] + 1)*NStr[nn],2))        # for [b,h] dimensions.
                ColDims = np.zeros(((NbaysX[mm] + 1)*(NbaysZ + 1)*NStr[nn],2)) # for [b,h] dimensions.
                for i in range(NbaysX[mm]*(NbaysZ + 1)*NStr[nn]):
                    XBDims[i,0] = 0.30
                    XBDims[i,1] = 0.45
                for i in range(NbaysZ*(NbaysX[mm] + 1)*NStr[nn]):
                    ZBDims[i,0] = 0.30
                    ZBDims[i,1] = 0.45
                for i in range((NbaysX[mm] + 1)*(NbaysZ + 1)*NStr[nn]):
                    ColDims[i,0] = 0.30
                    ColDims[i,1] = 0.30
                    
                # Seismic design Category needed for strenght design.
                # ==================================================
                SDC = \
                    ASCE716.detSDCat(RiskCat,Ss,S1,Vso[oo]) # [STRING] seismic design category of building.
                alldelta = ASCE716.detdelta_a(SDC,RiskCat)        # [m/m] allowable drift ratio
                
                # Data for ploting after analysis.
                # --------------------------------
                hxplot = []
                Da = []            # [m/m] allowable drift ratio.
                for i in range(NStr[nn]+1):
                    hxplot.append((i)*StoryH)
                    Da.append(alldelta)
                    
    # ****************************************************************************
    # ****************************************************************************
    # ****************************************************************************
    # ****************************************************************************
                    # DIMENSIONING OF BUILDING #
    # ****************************************************************************
    # ****************************************************************************
    # ****************************************************************************
    # ****************************************************************************                
                
                # ITERATIVE PROCESS TO DETERMINE ELEMNTS GEOMETRY.
                # =================================================
                maxDi_h = np.zeros((1,len(dirs)))            # [m/m] maximum inelastic Drift ratios in building height.
                niters = 1
                maxDi_h[0,0], maxDi_h[0,1] = 2*alldelta, 2*alldelta
                
                while maxDi_h[0,0] > alldelta or maxDi_h[0,1] > alldelta:
                    
                    if ColDims[0,0] > XbayL/2:
                        # print('Element dimensions too long for frame structure!')
                        break
                    
                    # Redefine sections of columns for both h and b dimensions depending on the
                    # plan aspect ratio of the structure.
                    if NbaysX[mm] != NbaysZ:
                        incbCol, inchCol  = 0.050, 0.025
                        incbB, inchB = 0.02, 0.025
                    elif NbaysX[mm] == NbaysZ:
                        incbCol, inchCol = 0.05, 0.05
                        incbB, inchB = 0.02, 0.025
                    
                    # Redefining elements sections.
                    # ----------------------------
                    for i in range(NbaysX[mm]*(NbaysZ + 1)*NStr[nn]):
                        XBDims[i,0] += incbB
                        XBDims[i,1] += inchB
                    for i in range(NbaysZ*(NbaysX[mm] + 1)*NStr[nn]):
                        ZBDims[i,0] += incbB
                        ZBDims[i,1] += inchB
                    for i in range((NbaysX[mm] + 1)*(NbaysZ + 1)*NStr[nn]):
                        ColDims[i,0] += incbCol
                        ColDims[i,1] += inchCol
                        
                    # Predefinition of variables.
                    # ---------------------------
                    FxF = np.zeros((NStr[nn],len(dirs)))              # [kN] variable list of Fx values for strenght design.
                    FxD = np.zeros((NStr[nn],len(dirs)))              # [kN] variable list of Fx values for drit analysis.
                    Di = np.zeros((NStr[nn]+1,len(dirs)))            # [m] inelastic Drift.
                    Di_h = np.zeros((NStr[nn]+1,len(dirs)))          # [m/m] inelastic Drift Ratio.
                    T = np.zeros((1,len(dirs)))  #   T = [Tx,Tz] # [s] two first periods of vibration. Assumed translational.
                    Tforces = np.zeros((1,len(dirs)))
                    Tdrifts = np.zeros((1,len(dirs)))
                    
                    for ii in range(len(dirs)):
                        
                        ops.wipe()
                        # Model Generation
                        # ================
                        OPSDMOD.ModelGen(ndm,ndf)
                        # Basic nodes generation
                        # ======================
                        # NbaysX,NbaysZ,XbayL,ZbayL,StoryH,NStr,CX=0.5,CZ=0.5,h_beam=0.6,flex='No'
                        OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex)
                        # Generation of master nodes
                        # ===========================
                        OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex)
                        # Single Point constraint of base nodes - fixity of structure.
                        # =============================================================
                        OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,flex)
                        # Multi-Point constraint - Rigid Diaphragm assignment.
                        # =====================================================
                        OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],flex)
                        # Material generation for sections and elements.
                        # ==============================================
                        # Only oncrete and reinforcement steel.
                        OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                        
                        # ****************************************************************************
                        # ****************************************************************************
                        # ELASTIC ELEMENTS
                        # ==================
                        # Generation of Geometric Transformation objects for elements
                        # ===========================================================
                        OPSDMOD.GeomTransGen(ColTransfType)
                        # Element object generation, Columns and Beams in X and Z global directions.
                        # --------------------------------------------------------------------------
                        OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                           ColDims,BeamEffFact,ColEffFact,Ec,fy)
                        # OPSD.forceBeamColumnGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,XBDims,ZBDims,ColDims,gamma_conc,g,Qsd)
                        # Seems like massdensity of columns are not considered from elements different 
                        # from horizontal ones to construct the consistent mass matrix. The Periods
                        # obtained using this mass matrix construction option, I consider, not adequate.
                        # Better to use my function to define lumped Masses from elemnts and loads below. (JSBM)
                        # Lumped mass determination by "manual" (element-oriented) calculation of masses.
                        # -------------------------------------------------------------------------------
                        # (Not to be used while using forceBeamColumnGen with non-zero values of gamma_conc and Qsd)
                        [Wx,MassInputMatr] = \
                            OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,flex)
                        
                        # ****************************************************************************
                        # ****************************************************************************
                        # GRAVITY ANALYSIS EXCECUTED ACCOUNING FOR LIVE, DEAD AND SUPERDEAD LOADS.
                        # ========================================================================
                        OPSDMOD.DeadLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XBDims,ZBDims,ColDims,gamma_conc)
                        OPSDMOD.SuperDeadLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,Qsd)
                        OPSDMOD.LiveLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,Ql,Qlr)
                        
                        
                        # GRAVITY-LOADS-CASE ANALYSIS.
                        # ============================
                        ops.system('BandGeneral')
                        ops.constraints('Transformation')
                        ops.numberer('RCM')
                        ops.test('NormDispIncr',1.0e-6,100)
                        ops.algorithm('Newton')
                        ops.integrator('LoadControl', 1)
                        ops.analysis('Static')
                        ops.analyze(1)
                        # Vertical reactions Calculation for verification.
                        # ------------------------------------------------
                        ops.reactions()
                        YReact = 0
                        for i in range(1,(NbaysX[mm]+1)*(NbaysZ+1)+1):
                            YReact += ops.nodeReaction(i+1e4,2)
                            
    
                        # ==========================================================================
                        # MODAL ANALYSIS FOR DETERMINING FUNDAMENTAL PERIOD AND ITS DIRECTION.
                        # ==========================================================================
                        (T,Tmaxreal,Mast) = OPSDAN.ModalAnalysis(NStr[nn],B,L,Wx,flex)
                            
                            
                        # ****************************************************************************
                        # ****************************************************************************
                        # DETERMINATION OF LATERAL LOADS APPLIED TO EACH FLOOR, FxD.
                        # =========================================================
                        # Using functions from ASCE710Seismic Module.
                        (SDS,SD1,Ta,Tforces,Tdrifts,Cu,Tmax,CsF,CsD,VF,VD,FxF[:,ii],FxD[:,ii],VxF,VxD) = \
                            ASCE716.ELFP(Ss,S1,Vso[oo],Wx,R,Ie,hx,StrType,T[0,ii])
                        
                        
                        # ****************************************************************************
                        # ****************************************************************************
                        # EQUIVALENT LATERAL FORCE PROCEDURE - FORCES APPLICATION
                        # =======================================================
                        # It is done according to the direction of analysis defined by the STRING object
                        # direction. The vector of Fx loads must be one with the loads per floor ordered
                        # from first floor to roof.
                        ops.loadConst('-time',0.0)   # Mantaining the gravity loads constant along the
                                                 # ELF application.
                        
                        OPSDMOD.ELFPForceGen(dirs[ii],FxD[:,ii],flex)    # Appliying lateral force for drift analysis.
                        
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
                        # Lateral reactions for verification of Base Shear.
                        # -------------------------------------------------
                        ops.reactions()
                        LatReact = 0
                        for i in range(1,(NbaysX[mm]+1)*(NbaysZ+1)+1):
                            if dirs[ii] in ('x','X'):
                                LatReact += ops.nodeReaction(i+1e4,1)
                            elif dirs[ii] in ('z','Z'):
                                LatReact += ops.nodeReaction(i+1e4,3)
                        
                        # ****************************************************************************
                        # ****************************************************************************
                        # STORY DRIFT CALCULATION
                        # ========================
                        dei = [0]          # [m] elastic displacement vector of Master Nodes.
                        deltaei = [0]      # [m] elastic drift.
                        if dirs[ii] in ('x','X'):
                            anDOF = 1
                        elif dirs[ii] in ('z','Z'):
                            anDOF = 3
                        for i in range(1,NStr[nn]+1):
                            dei.append(ops.nodeDisp(int(999+1e4+1e5*i),anDOF))
                        for i in range(1,NStr[nn]+1):
                            deltaei.append(dei[i] - dei[i-1])
                        for i in range(len(deltaei)):
                            Di[i,ii] = (deltaei[i]*Cd/Ie)
                            Di_h[i,ii] = Di[i,ii]/StoryH
                            
                        maxDi_h[0,ii] = max(Di_h[:,ii])
                        
                    niters += 1
                    # print(str(niters) + ' Iterations have been executed.')
                    # Ploting last results of Drift Ratio to compare previous results
                    # ================================================================
                    # Ploting first results of Drift Ratio to compare after iterations
                    # ================================================================
                    plt.plot(Di_h[:,0],hxplot,Di_h[:,1],hxplot)
                    
                plt.plot(Da,hxplot)
                plt.grid()
                
                # print('Dimensioning process finished.')
                
    # ****************************************************************************
    # ****************************************************************************
    # ****************************************************************************
    # ****************************************************************************
                    # REINFORCING PROCESS OF BUILDING #
    # ****************************************************************************
    # ****************************************************************************
    # ****************************************************************************
    # ****************************************************************************
                # print('-----------------------------------')
                # print('Reinfrcing process started')
                # MODEL GENERATION FOR THE ENTIRE BUILDING WITH ELASTIC ELEMENTS AND REDUCED INERTIA
                # ===================================================================================
                # Once for each load case
                # =================
                # STATIC LOAD CASES
                # =================
                
                # *************************************************************************************
                # DEAD LOAD CASE
                # *************************************************************************************
                ops.wipe()
                OPSDMOD.ModelGen(ndm,ndf)
                OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex)
                OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex)
                OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,flex)
                OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],flex)
                OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                OPSDMOD.GeomTransGen(ColTransfType)
                OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                           ColDims,BeamEffFact,ColEffFact,Ec,fy)
                [Wx,MassInputMatr] = \
                    OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,flex)
                    
                # DEAD-LOAD-CASE ANALYSIS.
                # ============================
                OPSDMOD.DeadLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XBDims,ZBDims,ColDims,gamma_conc)
                ops.system('BandGeneral')
                ops.constraints('Transformation')
                ops.numberer('RCM')
                ops.test('NormDispIncr',1.0e-6,100)
                ops.algorithm('Newton')
                ops.integrator('LoadControl', 1)
                ops.analysis('Static')
                ops.analyze(1)
                # Vertical reactions Calculation for verification.
                # ------------------------------------------------
                ops.reactions()
                YReactD = 0
                for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                    YReactD += ops.nodeReaction(i+1+1e4,2)
                    
                # Storing moments for beams and axial forces for columns
                # ======================================================
                # Previous Calculations
                # ----------------------
                NbeamsX = NbaysX[mm]*(NbaysZ + 1)  # number per floor
                NbeamsZ = NbaysZ*(NbaysX[mm] + 1)  # number per floor
                Ncols = (NbaysX[mm] + 1)*(NbaysZ + 1)
                # Beams on X
                # -----------
                MDXbeams = np.zeros((NbeamsX*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsX+1):
                        MDXbeams[loc,0] = ops.sectionForce(int(j+1e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mz,My,T
                        MDXbeams[loc,1] = ops.sectionForce(int(j+1e4+1e5*i),int((N-1)/2 + 1),2)
                        MDXbeams[loc,2] = ops.sectionForce(int(j+1e4+1e5*i),N,2)
                        loc += 1
                # Beams on Z
                # -----------
                MDZbeams = np.zeros((NbeamsZ*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsZ+1):
                        MDZbeams[loc,0] = ops.sectionForce(int(j+2e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mz,My,T
                        MDZbeams[loc,1] = ops.sectionForce(int(j+2e4+1e5*i),int((N-1)/2 + 1),2)
                        MDZbeams[loc,2] = ops.sectionForce(int(j+2e4+1e5*i),N,2)
                        loc += 1
                # Columns
                # -------
                PDcols = np.zeros((Ncols*NStr[nn],1))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,Ncols+1):
                        PDcols[loc,0] = ops.sectionForce(int(j+3e4+1e5*i),1,1)
                        loc += 1
                # print('Dead Load Case finished')        
                # *************************************************************************************
                # SUPER DEAD LOAD CASE
                # *************************************************************************************
                ops.wipe()
                OPSDMOD.ModelGen(ndm,ndf)
                OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex)
                OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex)
                OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,flex)
                OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],flex)
                OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                OPSDMOD.GeomTransGen(ColTransfType)
                OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                           ColDims,BeamEffFact,ColEffFact,Ec,fy)
                [Wx,MassInputMatr] = \
                    OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,flex)
                    
                # SUPERDEAD-LOAD-CASE ANALYSIS.
                # ============================
                OPSDMOD.SuperDeadLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,Qsd)
                ops.system('BandGeneral')
                ops.constraints('Transformation')
                ops.numberer('RCM')
                ops.test('NormDispIncr',1.0e-6,100)
                ops.algorithm('Newton')
                ops.integrator('LoadControl', 1)
                ops.analysis('Static')
                ops.analyze(1)
                # Vertical reactions Calculation for verification.
                # ------------------------------------------------
                ops.reactions()
                YReactSD = 0
                for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                    YReactSD += ops.nodeReaction(i+1+1e4,2)
                    
                # Storing moments for beams and axial forces for columns
                # ======================================================
                # Beams on X
                # -----------
                MSDXbeams = np.zeros((NbeamsX*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsX+1):
                        MSDXbeams[loc,0] = ops.sectionForce(int(j+1e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mz,My,T
                        MSDXbeams[loc,1] = ops.sectionForce(int(j+1e4+1e5*i),int((N-1)/2 + 1),2)
                        MSDXbeams[loc,2] = ops.sectionForce(int(j+1e4+1e5*i),N,2)
                        loc += 1
                # Beams on Z
                # -----------
                MSDZbeams = np.zeros((NbeamsZ*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsZ+1):
                        MSDZbeams[loc,0] = ops.sectionForce(int(j+2e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mz,My,T
                        MSDZbeams[loc,1] = ops.sectionForce(int(j+2e4+1e5*i),int((N-1)/2 + 1),2)
                        MSDZbeams[loc,2] = ops.sectionForce(int(j+2e4+1e5*i),N,2)
                        loc += 1
                # Columns
                # -------
                PSDcols = np.zeros((Ncols*NStr[nn],1))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,Ncols+1):
                        PSDcols[loc,0] = ops.sectionForce(int(j+3e4+1e5*i),1,1)
                        loc += 1
                # print('Super-Dead Load Case finished')
                # *************************************************************************************
                # LIVE LOAD CASE
                # *************************************************************************************
                ops.wipe()
                OPSDMOD.ModelGen(ndm,ndf)
                OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex)
                OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex)
                OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,flex)
                OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],flex)
                OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                OPSDMOD.GeomTransGen(ColTransfType)
                OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                           ColDims,BeamEffFact,ColEffFact,Ec,fy)
                [Wx,MassInputMatr] = \
                    OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,flex)
                    
                # LIVE-LOAD-CASE ANALYSIS.
                # ============================
                OPSDMOD.LiveLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,Ql,Qlr)
                ops.system('BandGeneral')
                ops.constraints('Transformation')
                ops.numberer('RCM')
                ops.test('NormDispIncr',1.0e-6,100)
                ops.algorithm('Newton')
                ops.integrator('LoadControl', 1)
                ops.analysis('Static')
                ops.analyze(1)
                # Vertical reactions Calculation for verification.
                # ------------------------------------------------
                ops.reactions()
                YReactL = 0
                for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                    YReactL += ops.nodeReaction(i+1+1e4,2)
                    
                # Storing moments for beams and axial forces for columns
                # ======================================================
                # Beams on X
                # -----------
                MLXbeams = np.zeros((NbeamsX*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsX+1):
                        MLXbeams[loc,0] = ops.sectionForce(int(j+1e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mz,My,T
                        MLXbeams[loc,1] = ops.sectionForce(int(j+1e4+1e5*i),int((N-1)/2 + 1),2)
                        MLXbeams[loc,2] = ops.sectionForce(int(j+1e4+1e5*i),N,2)
                        loc += 1
                # Beams on Z
                # -----------
                MLZbeams = np.zeros((NbeamsZ*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsZ+1):
                        MLZbeams[loc,0] = ops.sectionForce(int(j+2e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mz,My,T
                        MLZbeams[loc,1] = ops.sectionForce(int(j+2e4+1e5*i),int((N-1)/2 + 1),2)
                        MLZbeams[loc,2] = ops.sectionForce(int(j+2e4+1e5*i),N,2)
                        loc += 1
                # Columns
                # -------
                PLcols = np.zeros((Ncols*NStr[nn],1))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,Ncols+1):
                        PLcols[loc,0] = ops.sectionForce(int(j+3e4+1e5*i),1,1)
                        loc += 1
                        
                # print('Live Load Case finished')
                # *************************************************************************************
                # LATERAL EQUIVALENT LOAD CASE - X DIRECTION
                # *************************************************************************************
                ops.wipe()
                OPSDMOD.ModelGen(ndm,ndf)
                OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex)
                OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex)
                OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,flex)
                OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],flex)
                OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                OPSDMOD.GeomTransGen(ColTransfType)
                OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                           ColDims,BeamEffFact,ColEffFact,Ec,fy)
                [Wx,MassInputMatr] = \
                    OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,flex)
                
                # Lateral load determination, Fx - X DIRECTION
                # --------------------------------------------
                # Already calculated at the dimensioning stage.
                        
                # EQUIVALENT LATERAL FORCE PROCEDURE - FORCES APPLICATION
                # =======================================================
                ops.loadConst('-time',0.0)   # Mantaining the gravity loads constant along the
                                         # ELF application.
                
                OPSDMOD.ELFPForceGen(dirs[0],FxF[:,0],flex)
                
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
                # Lateral reactions for verification of Base Shear.
                # -------------------------------------------------
                ops.reactions()
                LatReactX = 0
                for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                    if dirs[0] in ('x','X'):
                        LatReactX += ops.nodeReaction(i+1+1e4,1)
                    elif dirs[0] in ('z','Z'):
                        LatReactX += ops.nodeReaction(i+1+1e4,3)
                
                # Storing moments for beams and axial forces for columns
                # ======================================================
                # Beams on X
                # -----------
                MELFXXbeams = np.zeros((NbeamsX*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsX+1):
                        MELFXXbeams[loc,0] = ops.sectionForce(int(j+1e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mx,My,T
                        MELFXXbeams[loc,1] = ops.sectionForce(int(j+1e4+1e5*i),int((N-1)/2 + 1),2)
                        MELFXXbeams[loc,2] = ops.sectionForce(int(j+1e4+1e5*i),N,2)
                        loc += 1
                # Beams on Z
                # -----------
                MELFXZbeams = np.zeros((NbeamsZ*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsZ+1):
                        MELFXZbeams[loc,0] = ops.sectionForce(int(j+2e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mx,My,T
                        MELFXZbeams[loc,1] = ops.sectionForce(int(j+2e4+1e5*i),int((N-1)/2 + 1),2)
                        MELFXZbeams[loc,2] = ops.sectionForce(int(j+2e4+1e5*i),N,2)
                        loc += 1
                # Columns
                # -------
                PELFXcols = np.zeros((Ncols*NStr[nn],1))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,Ncols+1):
                        PELFXcols[loc,0] = ops.sectionForce(int(j+3e4+1e5*i),1,1)
                        loc += 1
                
                # print('Lateral Load Case - X Direction finished')
                # *************************************************************************************
                # LATERAL EQUIVALENT LOAD CASE - Z DIRECTION
                # *************************************************************************************
                ops.wipe()
                OPSDMOD.ModelGen(ndm,ndf)
                OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex)
                OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex)
                OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,flex)
                OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],flex)
                OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                OPSDMOD.GeomTransGen(ColTransfType)
                OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                           ColDims,BeamEffFact,ColEffFact,Ec,fy)
                [Wx,MassInputMatr] = \
                    OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,flex)
                
                # Lateral load determination, Fx - Z DIRECTION
                # --------------------------------------------
                # Already calculated at the dimensioning stage.
                        
                # EQUIVALENT LATERAL FORCE PROCEDURE - FORCES APPLICATION
                # =======================================================
                ops.loadConst('-time',0.0)   # Mantaining the gravity loads constant along the
                                         # ELF application.
                
                OPSDMOD.ELFPForceGen(dirs[1],FxF[:,1],flex)
                
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
                # Lateral reactions for verification of Base Shear.
                # -------------------------------------------------
                ops.reactions()
                LatReactZ = 0
                for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                    if dirs[1] in ('x','X'):
                        LatReactZ += ops.nodeReaction(i+1+1e4,1)
                    elif dirs[1] in ('z','Z'):
                        LatReactZ += ops.nodeReaction(i+1+1e4,3)
                
                # Storing moments for beams and axial forces for columns
                # ======================================================
                # Beams on X
                # -----------
                MELFZXbeams = np.zeros((NbeamsX*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsX+1):
                        MELFZXbeams[loc,0] = ops.sectionForce(int(j+1e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mx,My,T
                        MELFZXbeams[loc,1] = ops.sectionForce(int(j+1e4+1e5*i),int((N-1)/2 + 1),2)
                        MELFZXbeams[loc,2] = ops.sectionForce(int(j+1e4+1e5*i),N,2)
                        loc += 1
                # Beams on Z
                # -----------
                MELFZZbeams = np.zeros((NbeamsZ*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsZ+1):
                        MELFZZbeams[loc,0] = ops.sectionForce(int(j+2e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mx,My,T
                        MELFZZbeams[loc,1] = ops.sectionForce(int(j+2e4+1e5*i),int((N-1)/2 + 1),2)
                        MELFZZbeams[loc,2] = ops.sectionForce(int(j+2e4+1e5*i),N,2)
                        loc += 1
                # Columns
                # -------
                PELFZcols = np.zeros((Ncols*NStr[nn],1))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,Ncols+1):
                        PELFZcols[loc,0] = ops.sectionForce(int(j+3e4+1e5*i),1,1)
                        loc += 1
                
                # print('Lateral Load Case - Z Direction finished')
                
                # Ploting first results of Drift Ratio to compare after iterations
                # # ================================================================
                # plt.figure()
                # plt.plot(Di_h[:,0],hxplot,Di_h[:,1],hxplot,Da,hxplot)
                # plt.grid()
                
                # ===========================================================================
                # LOAD COMBINATIONS FOR THE ANALYSED LOAD CASES 
                # ===========================================================================
                # There are the following load combinations:
                # U1 = 1.20*D + 1.00*L + 1.00*rho(EX+0.30EZ)
                # U2 = 1.20*D + 1.00*L + 1.00*rho(0.30EX+EZ)
                # U3 = 1.20*D + 1.00*L - 1.00*rho(Ex+0.30EZ)
                # U4 = 1.20*D + 1.00*L - 1.00*rho(0.30Ex+EZ)
                #       Where rho is the redundancy factor.
                rho = ASCE716.RedFact(SDC)
                # Combianations for Xbeams
                # ------------------------
                MuBX = np.zeros((NbeamsX*NStr[nn],6))   # [kN-m]
                for i in range(NbeamsX*NStr[nn]):
                    MuBX[i,0] = max(1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] + 1.00*rho*(MELFXXbeams[i,0]+0.30*MELFZXbeams[i,0]),
                                    1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] - 1.00*rho*(MELFXXbeams[i,0]+0.30*MELFZXbeams[i,0]),\
                                        1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] + 1.00*rho*(0.30*MELFXXbeams[i,0]+MELFZXbeams[i,0]),\
                                            1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] - 1.00*rho*(0.30*MELFXXbeams[i,0]+MELFZXbeams[i,0]))
                    MuBX[i,1] = min(1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] + 1.00*rho*(MELFXXbeams[i,0]+0.30*MELFZXbeams[i,0]),
                                    1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] - 1.00*rho*(MELFXXbeams[i,0]+0.30*MELFZXbeams[i,0]),\
                                        1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] + 1.00*rho*(0.30*MELFXXbeams[i,0]+MELFZXbeams[i,0]),\
                                            1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] - 1.00*rho*(0.30*MELFXXbeams[i,0]+MELFZXbeams[i,0]))
                    MuBX[i,2] = max(1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] + 1.00*rho*(MELFXXbeams[i,1]+0.30*MELFZXbeams[i,1]),
                                    1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] - 1.00*rho*(MELFXXbeams[i,1]+0.30*MELFZXbeams[i,1]),\
                                        1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] + 1.00*rho*(0.30*MELFXXbeams[i,1]+MELFZXbeams[i,1]),\
                                            1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] - 1.00*rho*(0.30*MELFXXbeams[i,1]+MELFZXbeams[i,1]))
                    MuBX[i,3] = min(1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] + 1.00*rho*(MELFXXbeams[i,1]+0.30*MELFZXbeams[i,1]),
                                    1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] - 1.00*rho*(MELFXXbeams[i,1]+0.30*MELFZXbeams[i,1]),\
                                        1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] + 1.00*rho*(0.30*MELFXXbeams[i,1]+MELFZXbeams[i,1]),\
                                            1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] - 1.00*rho*(0.30*MELFXXbeams[i,1]+MELFZXbeams[i,1]))
                    MuBX[i,4] = max(1.20*(MDXbeams[i,2]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] + 1.00*rho*(MELFXXbeams[i,2]+0.30*MELFZXbeams[i,2]),
                                    1.20*(MDXbeams[i,2]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] - 1.00*rho*(MELFXXbeams[i,2]+0.30*MELFZXbeams[i,2]),\
                                        1.20*(MDXbeams[i,2]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] + 1.00*rho*(0.30*MELFXXbeams[i,2]+MELFZXbeams[i,2]),\
                                            1.20*(MDXbeams[i,1]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] - 1.00*rho*(0.30*MELFXXbeams[i,2]+MELFZXbeams[i,2]))
                    MuBX[i,5] = min(1.20*(MDXbeams[i,2]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] + 1.00*rho*(MELFXXbeams[i,2]+0.30*MELFZXbeams[i,2]),
                                    1.20*(MDXbeams[i,2]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] - 1.00*rho*(MELFXXbeams[i,2]+0.30*MELFZXbeams[i,2]),\
                                        1.20*(MDXbeams[i,2]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] + 1.00*rho*(0.30*MELFXXbeams[i,2]+MELFZXbeams[i,2]),\
                                            1.20*(MDXbeams[i,1]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] - 1.00*rho*(0.30*MELFXXbeams[i,2]+MELFZXbeams[i,2]))
                        
                # Combianations for Zbeams
                # ------------------------
                MuBZ = np.zeros((NbeamsZ*NStr[nn],6))   # [kN-m]
                for i in range(NbeamsZ*NStr[nn]):
                    MuBZ[i,0] = max(1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] + 1.00*rho*(MELFXZbeams[i,0]+0.30*MELFZZbeams[i,0]),
                                    1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] - 1.00*rho*(MELFXZbeams[i,0]+0.30*MELFZZbeams[i,0]),\
                                        1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] + 1.00*rho*(0.30*MELFXZbeams[i,0]+MELFZZbeams[i,0]),\
                                            1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] - 1.00*rho*(0.30*MELFXZbeams[i,0]+MELFZZbeams[i,0]))
                    MuBZ[i,1] = min(1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] + 1.00*rho*(MELFXZbeams[i,0]+0.30*MELFZZbeams[i,0]),
                                    1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] - 1.00*rho*(MELFXZbeams[i,0]+0.30*MELFZZbeams[i,0]),\
                                        1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] + 1.00*rho*(0.30*MELFXZbeams[i,0]+MELFZZbeams[i,0]),\
                                            1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] - 1.00*rho*(0.30*MELFXZbeams[i,0]+MELFZZbeams[i,0]))
                    MuBZ[i,2] = max(1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] + 1.00*rho*(MELFXZbeams[i,1]+0.30*MELFZZbeams[i,1]),
                                    1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] - 1.00*rho*(MELFXZbeams[i,1]+0.30*MELFZZbeams[i,1]),\
                                        1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] + 1.00*rho*(0.30*MELFXZbeams[i,1]+MELFZZbeams[i,1]),\
                                            1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] - 1.00*rho*(0.30*MELFXZbeams[i,1]+MELFZZbeams[i,1]))
                    MuBZ[i,3] = min(1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] + 1.00*rho*(MELFXZbeams[i,1]+0.30*MELFZZbeams[i,1]),
                                    1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] - 1.00*rho*(MELFXZbeams[i,1]+0.30*MELFZZbeams[i,1]),\
                                        1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] + 1.00*rho*(0.30*MELFXZbeams[i,1]+MELFZZbeams[i,1]),\
                                            1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] - 1.00*rho*(0.30*MELFXZbeams[i,1]+MELFZZbeams[i,1]))
                    MuBZ[i,4] = max(1.20*(MDZbeams[i,2]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] + 1.00*rho*(MELFXZbeams[i,2]+0.30*MELFZZbeams[i,2]),
                                    1.20*(MDZbeams[i,2]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] - 1.00*rho*(MELFXZbeams[i,2]+0.30*MELFZZbeams[i,2]),\
                                        1.20*(MDZbeams[i,2]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] + 1.00*rho*(0.30*MELFXZbeams[i,2]+MELFZZbeams[i,2]),\
                                            1.20*(MDZbeams[i,1]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] - 1.00*rho*(0.30*MELFXZbeams[i,2]+MELFZZbeams[i,2]))
                    MuBZ[i,5] = min(1.20*(MDZbeams[i,2]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] + 1.00*rho*(MELFXZbeams[i,2]+0.30*MELFZZbeams[i,2]),
                                    1.20*(MDZbeams[i,2]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] - 1.00*rho*(MELFXZbeams[i,2]+0.30*MELFZZbeams[i,2]),\
                                        1.20*(MDZbeams[i,2]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] + 1.00*rho*(0.30*MELFXZbeams[i,2]+MELFZZbeams[i,2]),\
                                            1.20*(MDZbeams[i,1]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] - 1.00*rho*(0.30*MELFXZbeams[i,2]+MELFZZbeams[i,2]))
                        
                # Combianations for Columns
                # ------------------------
                PuCols = np.zeros((Ncols*NStr[nn],1))    # [kN] axial load in columns.
                for i in range(Ncols*NStr[nn]):
                    PuCols[i,0] = min(1.20*(PDcols[i,0]+PSDcols[i,0]) + 1.00*PLcols[i,0] + 1.00*rho*(PELFXcols[i,0]+0.30*PELFZcols[i,0]),
                                    1.20*(PDcols[i,0]+PSDcols[i,0]) + 1.00*PLcols[i,0] - 1.00*rho*(PELFXcols[i,0]+0.30*PELFZcols[i,0]),\
                                        1.20*(PDcols[i,0]+PSDcols[i,0]) + 1.00*PLcols[i,0] + 1.00*rho*(0.30*PELFXcols[i,0]+PELFZcols[i,0]),\
                                            1.20*(PDcols[i,0]+PSDcols[i,0]) + 1.00*PLcols[i,0] - 1.00*rho*(0.30*PELFXcols[i,0]+PELFZcols[i,0]))
                
                # print('Load Combinations done!')
                
                # ***********************************************************************************************************************************
                # ***********************************************************************************************************************************
                # REINFORCING OF BEAMS
                # ***********************************************************************************************************************************
                # ***********************************************************************************************************************************
                # BEAMS IN X DIRECTION
                # =====================
                
                # print('Starting Beam reinforcing....')
                
                MuBXDes = abs(MuBX)          # [kN-m] 
                                             # absolute value because the index into the object determines the sign within.
                XBreinf = np.zeros((NbaysX[mm]*(NbaysZ+1)*NStr[nn],4))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbaysZ+2):
                        for k in range(1,NbaysX[mm]+1):
                            if j == 1 or j == NbaysZ+1:
                                beamloc = 'edge'
                            else:
                                beamloc = 'interior'
                            XBreinf[loc,:] = \
                                DLBDefs.SpecialBeamDesign(abs(fpc),fy,E0,etlim,XBDims[loc,0],XBDims[loc,1],\
                                                          MuBXDes[loc,0],MuBXDes[loc,1],MuBXDes[loc,2],MuBXDes[loc,3],MuBXDes[loc,4],MuBXDes[loc,5],\
                                                              transtype,XbayL,ZbayL,beamloc,hslab=0.25)
                            loc += 1
                            
                # BEAMS IN Z DIRECTION
                # =====================
                MuBZDes = abs(MuBZ)     # [kN-m] . 
                                             # absolute value because the order into the object determines the sign within.
                ZBreinf = np.zeros((NbaysZ*(NbaysX[mm]+1)*NStr[nn],4))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbaysX[mm]+2):
                        for k in range(1,NbaysZ+1):
                            if j == 1 or j == NbaysX[mm]+1:
                                beamloc = 'edge'
                            else:
                                beamloc = 'interior'
                            ZBreinf[loc,:] = \
                                DLBDefs.SpecialBeamDesign(abs(fpc),fy,E0,etlim,ZBDims[loc,0],ZBDims[loc,1],\
                                                          MuBZDes[loc,0],MuBZDes[loc,1],MuBZDes[loc,2],MuBZDes[loc,3],MuBZDes[loc,4],MuBZDes[loc,5],\
                                                              transtype,XbayL,ZbayL,beamloc,hslab=0.25)
                            loc += 1
                
                # print('Beams reinforcing done!')
                # print('Calculating and alocating Nominal Moments of Beams')
                
                # NOMINAL MOMENTS FROM BEAMS ON NODES FOR COLUMN DESIGN.
                # ======================================================
                NAX = NbaysX[mm] + 1
                NAZ = NbaysZ + 1
                MnbXdir = np.zeros((Ncols*NStr[nn],1))   # [kN-m] nominal moments at nodes
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NAZ+1):
                        for k in range(1,NbaysX[mm]+1):
                            Nodei = k + (j-1)*NAX + (i-1)*Ncols - 1      # the order in which the list of columns is stored, not the node numbering.
                            Nodej = k + 1 + (j-1)*NAX + (i-1)*Ncols - 1
                            (Mnipos,Mnineg,Mnjpos,Mnjneg) = \
                                DLBDefs.BeamEndsMn(abs(fpc),fy,E0,XBDims[loc,0],XBDims[loc,1],\
                                                   XBreinf[loc,0],XBreinf[loc,1],XBreinf[loc,2],XBreinf[loc,3],transtype)
                            MnbXdir[Nodei,0] += Mnipos
                            MnbXdir[Nodej,0] += Mnjneg
                            loc += 1
                            # print(str(loc)+' X-Beams analized!')
                MnbZdir = np.zeros((Ncols*NStr[nn],1))   # [kN-m] nominal moments at nodes
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NAX+1):
                        for k in range(1,NbaysZ+1):
                            Nodei = j + (k-1)*NAX + (i-1)*Ncols - 1
                            Nodej = j + (k-1)*NAX + NAX + (i-1)*Ncols - 1
                            (Mnipos,Mnineg,Mnjpos,Mnjneg) = \
                                DLBDefs.BeamEndsMn(abs(fpc),fy,E0,XBDims[loc,0],XBDims[loc,1],\
                                                   ZBreinf[loc,0],ZBreinf[loc,1],ZBreinf[loc,2],ZBreinf[loc,3],transtype)
                            MnbZdir[Nodei,0] += Mnipos
                            MnbZdir[Nodej,0] += Mnjneg
                            loc += 1
                            # print(str(loc)+' Z-Beams analized!')
                
                # print('Nominal moments from beams stored in nodes')
                
                # ================================            
                # Design of column reinforcemwnt.
                # ================================
                # print('Starting Columns reinforcing....')
                loc = 0
                rhomincol = 0.01
                # rhomaxcol = 0.04
                Colreinf = np.zeros((Ncols*NStr[nn],10))
                for i in range(1,Ncols*NStr[nn]+1):
                    
                    # X Direction
                    # ------------
                    # bcolx = ColDims[loc,0]
                    hcolx = ColDims[loc,1]
                    # bcorex = bcolx - 2*cover
                    hcorex = hcolx - 2*cover
                    edx = hcorex/4
                    dx = [cover,cover+edx,cover+2*edx,cover+3*edx,cover+4*edx]
                    # Z Direction
                    # ------------
                    # bcolz = ColDims[loc,1]
                    hcolz = ColDims[loc,0]
                    # bcorez = bcolz - 2*cover
                    hcorez = hcolz - 2*cover
                    edz = hcorez/3          # because there are 4 layers in this direction!!
                    dz = [cover,cover+edz,cover+2*edz,cover+3*edz]
                    
                    
                    Ascol = rhomincol*ColDims[loc,0]*ColDims[loc,1]
                    Mnx = 0.0
                    Mnz = 0.0
                    while 2*Mnx < 1.20*MnbXdir[loc,0] or 2*Mnz < 1.20*MnbZdir[loc,0]:
                        Asz = [5/14*Ascol,2/14*Ascol,2/14*Ascol,5/14*Ascol]
                        Asx = [4/14*Ascol,2/14*Ascol,2/14*Ascol,2/14*Ascol,4/14*Ascol]
                        # X Direction
                        # ------------
                        (Mnx,Pn,phi) =\
                            DLBDefs.ColAnalysisMn(abs(fpc),fy,E0,etlim,abs(PuCols[loc,0]),ColDims[loc,1],ColDims[loc,0],\
                                                  dx,Asx,5,transtype)
                        # Z Direction
                        # ------------
                        (Mnz,Pn,phi) =\
                            DLBDefs.ColAnalysisMn(abs(fpc),fy,E0,etlim,abs(PuCols[loc,0]),ColDims[loc,0],ColDims[loc,1],\
                                                  dz,Asz,4,transtype)
                        Ascol += 0.0001
                        
                    Colreinf[loc,:5] = 4/14*Ascol, 4/14*Ascol, 2/14*Ascol, 2/14*Ascol, 2/14*Ascol
                    Colreinf[loc,5:] = 4/14*Ascol, 4/14*Ascol, 2/14*Ascol, 2/14*Ascol, 2/14*Ascol
                    loc += 1
                    # print(str(loc)+' Cols with reinforcement')
                
              
                counter += 1
                print(f"{round(counter/total*100,2)} %  --  {counter} designed out of {total}")
                varpath = directory+'\\RegularDesign\\'+\
                    str(NbaysX[mm])+'BayX'+str(NbaysZ)+'BayZ'+str(NStr[nn])+'FLRS'+str(Vso[oo])+'.pickle'
                with open(varpath,'wb') as f:
                    pickle.dump([ColDims,Colreinf,XBDims,XBreinf,ZBDims,ZBreinf,Di_h,FxD,FxF,T,Ta,Tmax,VxD,\
                                  VxF,Da,hxplot],f)
                        
    return ColDims,Colreinf,XBDims,XBreinf,ZBDims,ZBreinf,Di_h,FxD,FxF,T,Ta,Tmax,VxD,VxF,Da,hxplot



def SSIDesign(ndm,ndf,NbaysX,NbaysZ,XbayL,ZbayL,StoryH,NStr,B,XbeamOff,ZbeamOff,ColOff,
                   fpc,fy,E0,etlim,bsteel,Ec,nuconc,gamma_conc,rec,cover,transtype,BeamEffFact,
                   ColEffFact,g,ColTransfType,N,gamma_soil,nu_soil,omega_soil,D,Re,Qsd,Ql,Qlr,
                   Ss,S1,R,Cd,Omega_o,Ie,Vso,StrType,RiskCat,dirs,directory):
    
    import numpy as np
    import matplotlib.pyplot as plt
    import pickle
    import DLBDefs
    import ASCE716Seismic as ASCE716
    import ASCE716SSI as SSI716
    import openseespy.opensees as ops
    import OPSDefsMOD as OPSDMOD
    import OPSDefsAN as OPSDAN
    
    counter = 0
    total = len(NbaysX)*len(NStr)*len(Vso)
    
    for mm in range(len(NbaysX)):
        for nn in range(len(NStr)):
            for oo in range(len(Vso)):
                    
                L = NbaysX[mm]*XbayL               # [m] plan length of building. L >= B always!
                    
                # Section Properties for element (as initial setup)
                # -------------------------------------------------
                XBDims = np.zeros((NbaysX[mm]*(NbaysZ + 1)*NStr[nn],2))        # for [b,h] dimensions. b is normal to local y-axis.
                ZBDims = np.zeros((NbaysZ*(NbaysX[mm] + 1)*NStr[nn],2))        # for [b,h] dimensions.
                ColDims = np.zeros(((NbaysX[mm] + 1)*(NbaysZ + 1)*NStr[nn],2)) # for [b,h] dimensions.
                for i in range(NbaysX[mm]*(NbaysZ + 1)*NStr[nn]):
                    XBDims[i,0] = 0.30
                    XBDims[i,1] = 0.45
                for i in range(NbaysZ*(NbaysX[mm] + 1)*NStr[nn]):
                    ZBDims[i,0] = 0.30
                    ZBDims[i,1] = 0.45
                for i in range((NbaysX[mm] + 1)*(NbaysZ + 1)*NStr[nn]):
                    ColDims[i,0] = 0.30
                    ColDims[i,1] = 0.30
                    
                # Seismic design Category needed for strenght design.
                # ==================================================
                SDC = \
                    ASCE716.detSDCat(RiskCat,Ss,S1,Vso[oo]) # [STRING] seismic design category of building.
                alldelta = ASCE716.detdelta_a(SDC,RiskCat)        # [m/m] allowable drift ratio
                
                # Data for ploting after analysis.
                # --------------------------------
                hxplot = list(range(0,NStr[nn]*3+1,3))
                Da = alldelta*np.ones((len(hxplot)))
                    
    # ****************************************************************************
    # ****************************************************************************
    # ****************************************************************************
    # ****************************************************************************
                    # DIMENSIONING OF BUILDING #
    # ****************************************************************************
    # ****************************************************************************
    # ****************************************************************************
    # ****************************************************************************                
                
                # ITERATIVE PROCESS TO DETERMINE ELEMENTS GEOMETRY.
                # =================================================
                maxDi_h = np.zeros((1,len(dirs)))            # [m/m] maximum inelastic Drift ratios in building height.
                niters = 1
                maxDi_h[0,0], maxDi_h[0,1] = 2*alldelta, 2*alldelta
                
                while maxDi_h[0,0] > alldelta or maxDi_h[0,1] > alldelta:
                    
                    if ColDims[0,0] > XbayL/2:
                        # print('Element dimensions too long for frame structure!')
                        break
                    
                    # Redefine sections of columns for both h and b dimensions depending on the
                    # plan aspect ratio of the structure.
                    if NbaysX[mm] != NbaysZ:
                        incbCol, inchCol  = 0.050, 0.025
                        incbB, inchB = 0.02, 0.025
                    elif NbaysX[mm] == NbaysZ:
                        incbCol, inchCol = 0.05, 0.05
                        incbB, inchB = 0.02, 0.025
                    
                    # Redefining elements sections.
                    # ----------------------------
                    for i in range(NbaysX[mm]*(NbaysZ + 1)*NStr[nn]):
                        XBDims[i,0] += incbB
                        XBDims[i,1] += inchB
                    for i in range(NbaysZ*(NbaysX[mm] + 1)*NStr[nn]):
                        ZBDims[i,0] += incbB
                        ZBDims[i,1] += inchB
                    for i in range((NbaysX[mm] + 1)*(NbaysZ + 1)*NStr[nn]):
                        ColDims[i,0] += incbCol
                        ColDims[i,1] += inchCol
                        
                    # Predefinition of variables.
                    # ---------------------------
                    FxF = np.zeros((NStr[nn]+1,len(dirs)))              # [kN] variable list of Fx values for strenght design.
                    FxD = np.zeros((NStr[nn]+1,len(dirs)))              # [kN] variable list of Fx values for drit analysis.
                    # Di = np.zeros((NStr[nn]+1,len(dirs)))            # [m] inelastic Drift.
                    # Di_h = np.zeros((NStr[nn]+1,len(dirs)))          # [m/m] inelastic Drift Ratio.
                    T = np.zeros((1,len(dirs)))  #   T = [Tx,Tz] # [s] two first periods of vibration. Assumed translational.
                    Tforces = np.zeros((1,len(dirs)))
                    Tdrifts = np.zeros((1,len(dirs)))
                    
                    # SSI predefinitions
                    # ------------------
                    FxFSSI = np.zeros((NStr[nn]+1,len(dirs)))              # [kN] variable list of Fx values for strenght design.
                    FxDSSI = np.zeros((NStr[nn]+1,len(dirs)))              # [kN] variable list of Fx values for drit analysis considering foundation flex, hence the NStr+1.
                    VgorD = np.zeros((1,len(dirs)))                  # [kN] reduced base shear.
                    DeltaVD = np.zeros((1,len(dirs)))                # [kN] reduction for base shear accounting for SSI.
                    VgorF = np.zeros((1,len(dirs)))                  # [kN] reduced base shear.
                    DeltaVF = np.zeros((1,len(dirs)))                # [kN] reduction for base shear accounting for SSI.
                    # Tg = np.zeros((1,len(dirs)))                    # [s] modified period of the structure.
                    DiSSI = np.zeros((NStr[nn]+1,len(dirs)))            # [m] inelastic Drift.
                    DiSSI_h = np.zeros((NStr[nn]+1,len(dirs)))          # [m/m] inelastic Drift Ratio.
                    
                    flex = ['No','Yes']                    # [STRING] 'Yes' or 'No' to consider foundation flexibility at the foundation interface.
                    
                    T = np.zeros((len(flex),3))            # [s] first row of this matrix corresponds to fixed base periods, while the second to modified by foundation flexibility.
                    Tmaxreal = np.zeros((len(flex)))
                    Mast = np.zeros((len(flex),len(dirs)))
                    
                    for ii in range(len(dirs)):      # loop iterator for determining drift values in each direction.
                        
                        for jj in range(len(flex)):    # loop iterator for determining both fixed-base and flexible foundation fundamental periods.
                            
                            if flex[jj] in ('Y','YES','Yes','yES','yes','y'):
                                # [m][LIST] with level height starting from first level or from base if foundation flexibility is included.
                                hx = [0.0001]
                            else:
                                hx = []
                            for i in range(NStr[nn]):
                                hx.append((i+1)*StoryH)
                            
                            ops.wipe()
                            # Model Generation
                            # ================
                            OPSDMOD.ModelGen(ndm,ndf)
                            # Basic nodes generation
                            # ======================
                            OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex[jj])
                            # Generation of master nodes
                            # ===========================
                            OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],flex[jj])
                            # Single Point constraint of base nodes - fixity of structure.
                            # =============================================================
                            OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,flex[jj])
                            # Multi-Point constraint - Rigid Diaphragm assignment.
                            # =====================================================
                            OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],flex[jj])
                            # Material generation for sections and elements.
                            # ==============================================
                            # Only oncrete and reinforcement steel.
                            OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                            # ****************************************************************************
                            # ****************************************************************************
                            # ELASTIC ELEMENTS
                            # ==================
                            # Generation of elastic element sections 
                            # ---------------------------------------
                            if flex[jj] in ('Y','YES','Yes','yES','yes','y'):
                                # Materials generation: stiffness constants accounting for soil flexibility.
                                # ---------------------------------------------------------------------------
                                OPSDMOD.FoundFlexMaterials(NbaysX[mm],NbaysZ,XbayL,ZbayL,Ss,Vso[oo],gamma_soil,nu_soil,B,L,Re,D,omega_soil,'grav')
                                # Zero-Length elements creation for connecting base nodes.
                                OPSDMOD.FoundFlexZLElements(NbaysX[mm],NbaysZ,XbayL,ZbayL,B,L,Re)
                            
                            # ****************************************************************************
                            # ****************************************************************************
                            # Generation of Geometric Transformation objects for elements
                            # ===========================================================
                            OPSDMOD.GeomTransGen(ColTransfType)
                            # Element object generation, Columns and Beams in X and Z global directions.
                            # --------------------------------------------------------------------------
                            OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                               ColDims,BeamEffFact,ColEffFact,Ec,fy)
                            # OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                               # ColDims,BeamEffFact,ColEffFact,Ec,fy)
                            # Seems like massdensity of columns are not considered from elements different 
                            # from horizontal ones to construct the consistent mass matrix. The Periods
                            # obtained using this mass matrix construction option, I consider, not adequate.
                            # Better to use my function to define lumped Masses from elemnts and loads below. (JSBM)
                            # Lumped mass determination by "manual" (element-oriented) calculation of masses.
                            # -------------------------------------------------------------------------------
                            # (Not to be used while using forceBeamColumnGen with non-zero values of gamma_conc and Qsd)
                            [Wx,MassInputMatr] = \
                                OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,flex[jj])
                            
                            # ****************************************************************************
                            # ****************************************************************************
                            # GRAVITY ANALYSIS EXCECUTED ACCOUNING FOR LIVE, DEAD AND SUPERDEAD LOADS.
                            # ========================================================================
                            OPSDMOD.DeadLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XBDims,ZBDims,ColDims,gamma_conc)
                            OPSDMOD.SuperDeadLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,Qsd)
                            OPSDMOD.LiveLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,Ql,Qlr)
                            
                            # GRAVITY-LOADS-CASE ANALYSIS.
                            # ============================
                            ops.system('BandGeneral')
                            ops.constraints('Transformation')
                            ops.numberer('RCM')
                            ops.test('NormDispIncr',1.0e-5,100)
                            ops.algorithm('Newton')
                            ops.integrator('LoadControl', 1)
                            ops.analysis('Static')
                            ops.analyze(1)
                            # Vertical reactions Calculation for verification.
                            # ------------------------------------------------
                            ops.reactions()
                            YReact = 0
                            if flex in ('Y','YES','Yes','yES','yes','y'):
                                for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                                    YReact += ops.nodeReaction(int(i+1),2)
                            else:
                                for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                                    YReact += ops.nodeReaction(int(i+1+1e4),2)
                            
                            # ==========================================================================
                            # MODAL ANALYSIS FOR DETERMINING FUNDAMENTAL PERIOD AND ITS DIRECTION.
                            # ==========================================================================
                            (T[jj,:],Tmaxreal[jj],Mast[jj,:]) = OPSDAN.ModalAnalysis(NStr[nn],B,L,Wx,flex[jj])
                            
                        # ****************************************************************************
                        # ****************************************************************************
                        # DETERMINATION OF LATERAL LOADS APPLIED TO EACH FLOOR, FxD.
                        # =========================================================
                        # Using functions from ASCE710Seismic Module.
                        (SDS,SD1,Ta,Tforces,Tdrifts,Cu,Tmax,CsF,CsD,VF,VD,FxF[:,ii],FxD[:,ii],VxF,VxD) = \
                            ASCE716.ELFP(Ss,S1,Vso[oo],Wx,R,Ie,hx,StrType,T[0,ii])     # Using fixed-base period.
                        
                        # ****************************************************************************
                        # ****************************************************************************
                        # DEFINITION OF SSI PARAMETERS FOR DRIFT CALCULATION AND DEFINING DESIGN FORCES
                        # ****************************************************************************
                        # ****************************************************************************
                        hn = hx[-1]                  # [m] total building height.
                        (VgorD,VgorF,DeltaVD,DeltaVF) = SSI716.\
                            detVgor(VD,VF,Ss,S1,R,Omega_o,Ie,Vso[oo],gamma_soil,nu_soil,B,L,hn,\
                                    Mast[1,ii],T[0,ii],dirs[ii],beta=0.05,Tg=T[1,ii],D=0,omega=0)
                        
                        # Force Distribution on building height.
                        # --------------------------------------
                        (FxDSSI[:,ii],_) = ASCE716.vertdistV(VgorD,T[1,ii],hx,Wx)
                        (FxFSSI[:,ii],_) = ASCE716.vertdistV(VgorF,T[1,ii],hx,Wx)
                        
                        # ****************************************************************************
                        # ****************************************************************************
                        # EQUIVALENT LATERAL FORCE PROCEDURE - FORCES APPLICATION
                        # =======================================================
                        ops.wipe()
                        OPSDMOD.ModelGen(ndm,ndf)
                        OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],'yes')
                        OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],'yes')
                        OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,'yes')
                        OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],'yes')
                        OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                        # Interface elements generation for foundation flexibility considerations.
                        # =========================================================================
                        OPSDMOD.FoundFlexMaterials(NbaysX[mm],NbaysZ,XbayL,ZbayL,Ss,Vso[oo],gamma_soil,nu_soil,B,L,Re,D,omega_soil,'lat')
                        # Zero-Length elements creation for connecting base nodes.
                        OPSDMOD.FoundFlexZLElements(NbaysX[mm],NbaysZ,XbayL,ZbayL,B,L,Re)
                        # ------------------------------------------------------------------------
                        OPSDMOD.GeomTransGen(ColTransfType)
                        OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                               ColDims,BeamEffFact,ColEffFact,Ec,fy)
                        [Wx,MassInputMatr] = \
                            OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,'yes')
                        
                        OPSDMOD.DeadLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XBDims,ZBDims,ColDims,gamma_conc)
                        OPSDMOD.SuperDeadLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,Qsd)
                        OPSDMOD.LiveLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,Ql,Qlr)
                        
                        # GRAVITY-LOADS-CASE ANALYSIS.
                        # ============================
                        ops.system('BandGeneral')
                        ops.constraints('Transformation')
                        ops.numberer('RCM')
                        ops.test('NormDispIncr',1.0e-5,100)
                        ops.algorithm('Newton')
                        ops.integrator('LoadControl', 1)
                        ops.analysis('Static')
                        ops.analyze(1)
                        ops.reactions()
                        YReaction = 0
                        for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                            YReaction += ops.nodeReaction(int(i+1),2)
                        # It is done according to the direction of analysis defined by the STRING object
                        # direction. The vector of Fx loads must be one with the loads per floor ordered
                        # from first floor to roof.
                            
                        # *** NOTE: remember that the last analysis was made while flex = 'yes'. So,
                        #           is the correct gravitational analysis to use on lateral load analysis.
                            
                        ops.loadConst('-time',0.0)   # Mantaining the gravity loads constant along the
                                                     # ELF application.
                        
                        OPSDMOD.ELFPForceGen(dirs[ii],FxDSSI[:,ii],flex[1])    # Appliying lateral force for drift analysis.
                        
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
                        # Lateral reactions for verification of Base Shear.
                        # -------------------------------------------------
                        ops.reactions()
                        LatReact = 0
                        for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                            if dirs[ii] in ('x','X'):
                                LatReact += ops.nodeReaction(int(i+1),1)
                            elif dirs[ii] in ('z','Z'):
                                LatReact += ops.nodeReaction(int(i+1),3)
                                
                        # ****************************************************************************
                        # ****************************************************************************
                        # STORY DRIFT CALCULATION
                        # ========================
                        dei = []          # [m] elastic displacement vector of Master Nodes.
                        deltaei = [0]      # [m] elastic drift.
                        if dirs[ii] in ('x','X'):
                            anDOF = 1
                        elif dirs[ii] in ('z','Z'):
                            anDOF = 3
                        # With foundation flexibility considerations
                        # ------------------------------------------
                        istory = 0
                        for i in range(len(hx)):     # for elastic displacement in each master node including that at the base.
                            dei.append(ops.nodeDisp(int(999+1e4+1e5*istory),anDOF))
                            istory += 1
                        deinorm = np.zeros((len(dei)))
                        deinorm[:] = dei
                        deinorm = deinorm - dei[0]
                        for i in range(1,NStr[nn]+1):
                            deltaei.append(deinorm[i] - deinorm[i-1])
                        for i in range(len(deltaei)):
                            DiSSI[i,ii] = (deltaei[i]*Cd/Ie)
                            DiSSI_h[i,ii] = DiSSI[i,ii]/StoryH
                            
                        maxDi_h[0,ii] = max(DiSSI_h[:,ii])
                        
                    niters += 1
    
                    # print(str(niters) + ' Iterations have been executed.')
                    # Ploting last results of Drift Ratio to compare previous results
                    # ================================================================
                    # Ploting first results of Drift Ratio to compare after iterations
                    # ================================================================
                    plt.plot(DiSSI_h[:,0],hxplot,DiSSI_h[:,1],hxplot)
                    
                plt.plot(Da,hxplot)
                plt.grid()
                
    # ****************************************************************************
    # ****************************************************************************
                    # REINFORCING PROCESS OF BUILDING #
    # ****************************************************************************
    # ****************************************************************************
    
                # MODEL GENERATION FOR THE ENTIRE BUILDING WITH ELASTIC ELEMENTS AND REDUCED INERTIA
                # ===================================================================================
                # Once for each load case
                # =================
                # STATIC LOAD CASES
                # =================
                
                # *************************************************************************************
                # DEAD LOAD CASE
                # *************************************************************************************
                ops.wipe()
                OPSDMOD.ModelGen(ndm,ndf)
                OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],'yes')
                OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],'yes')
                OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,'yes')
                OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],'yes')
                OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                # Interface elements generation for foundation flexibility considerations.
                # =========================================================================
                OPSDMOD.FoundFlexMaterials(NbaysX[mm],NbaysZ,XbayL,ZbayL,Ss,Vso[oo],gamma_soil,nu_soil,B,L,Re,D,omega_soil,'grav')
                # Zero-Length elements creation for connecting base nodes.
                OPSDMOD.FoundFlexZLElements(NbaysX[mm],NbaysZ,XbayL,ZbayL,B,L,Re)
                
                OPSDMOD.GeomTransGen(ColTransfType)
                OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                               ColDims,BeamEffFact,ColEffFact,Ec,fy)
                [Wx,MassInputMatr] = \
                    OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,'yes')
                    
                # DEAD-LOAD-CASE ANALYSIS.
                # ============================
                OPSDMOD.DeadLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XBDims,ZBDims,ColDims,gamma_conc)
                ops.system('BandGeneral')
                ops.constraints('Transformation')
                ops.numberer('RCM')
                ops.test('NormDispIncr',1.0e-6,100)
                ops.algorithm('Newton')
                ops.integrator('LoadControl', 1)
                ops.analysis('Static')
                ops.analyze(1)
                # Vertical reactions Calculation for verification.
                # ------------------------------------------------
                ops.reactions()
                YReactD = 0
                for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                    YReactD += ops.nodeReaction(int(i+1),2)
                    
                # Storing moments for beams and axial forces for columns
                # ======================================================
                # Previous Calculations
                # ----------------------
                NbeamsX = NbaysX[mm]*(NbaysZ + 1)  # number per floor
                NbeamsZ = NbaysZ*(NbaysX[mm] + 1)  # number per floor
                Ncols = (NbaysX[mm] + 1)*(NbaysZ + 1)
                # Beams on X
                # -----------
                MDXbeams = np.zeros((NbeamsX*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsX+1):
                        MDXbeams[loc,0] = ops.sectionForce(int(j+1e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mz,My,T
                        MDXbeams[loc,1] = ops.sectionForce(int(j+1e4+1e5*i),int((N-1)/2 + 1),2)
                        MDXbeams[loc,2] = ops.sectionForce(int(j+1e4+1e5*i),N,2)
                        loc += 1
                # Beams on Z
                # -----------
                MDZbeams = np.zeros((NbeamsZ*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsZ+1):
                        MDZbeams[loc,0] = ops.sectionForce(int(j+2e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mz,My,T
                        MDZbeams[loc,1] = ops.sectionForce(int(j+2e4+1e5*i),int((N-1)/2 + 1),2)
                        MDZbeams[loc,2] = ops.sectionForce(int(j+2e4+1e5*i),N,2)
                        loc += 1
                # Columns
                # -------
                PDcols = np.zeros((Ncols*NStr[nn],1))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,Ncols+1):
                        PDcols[loc,0] = ops.sectionForce(int(j+3e4+1e5*i),1,1)
                        loc += 1
                # print('Dead Load Case finished')        
                # *************************************************************************************
                # SUPER DEAD LOAD CASE
                # *************************************************************************************
                ops.wipe()
                OPSDMOD.ModelGen(ndm,ndf)
                OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],'yes')
                OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],'yes')
                OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,'yes')
                OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],'yes')
                OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                # Interface elements generation for foundation flexibility considerations.
                # =========================================================================
                OPSDMOD.FoundFlexMaterials(NbaysX[mm],NbaysZ,XbayL,ZbayL,Ss,Vso[oo],gamma_soil,nu_soil,B,L,Re,D,omega_soil,'grav')
                # Zero-Length elements creation for connecting base nodes.
                OPSDMOD.FoundFlexZLElements(NbaysX[mm],NbaysZ,XbayL,ZbayL,B,L,Re)
                
                OPSDMOD.GeomTransGen(ColTransfType)
                OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                               ColDims,BeamEffFact,ColEffFact,Ec,fy)
                [Wx,MassInputMatr] = \
                    OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,'yes')
                    
                # SUPERDEAD-LOAD-CASE ANALYSIS.
                # ============================
                OPSDMOD.SuperDeadLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,Qsd)
                ops.system('BandGeneral')
                ops.constraints('Transformation')
                ops.numberer('RCM')
                ops.test('NormDispIncr',1.0e-5,100)
                ops.algorithm('Newton')
                ops.integrator('LoadControl', 1)
                ops.analysis('Static')
                ops.analyze(1)
                # Vertical reactions Calculation for verification.
                # ------------------------------------------------
                ops.reactions()
                YReactSD = 0
                for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                    YReactSD += ops.nodeReaction(int(i+1),2)
                    
                # Storing moments for beams and axial forces for columns
                # ======================================================
                # Beams on X
                # -----------
                MSDXbeams = np.zeros((NbeamsX*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsX+1):
                        MSDXbeams[loc,0] = ops.sectionForce(int(j+1e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mz,My,T
                        MSDXbeams[loc,1] = ops.sectionForce(int(j+1e4+1e5*i),int((N-1)/2 + 1),2)
                        MSDXbeams[loc,2] = ops.sectionForce(int(j+1e4+1e5*i),N,2)
                        loc += 1
                # Beams on Z
                # -----------
                MSDZbeams = np.zeros((NbeamsZ*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsZ+1):
                        MSDZbeams[loc,0] = ops.sectionForce(int(j+2e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mz,My,T
                        MSDZbeams[loc,1] = ops.sectionForce(int(j+2e4+1e5*i),int((N-1)/2 + 1),2)
                        MSDZbeams[loc,2] = ops.sectionForce(int(j+2e4+1e5*i),N,2)
                        loc += 1
                # Columns
                # -------
                PSDcols = np.zeros((Ncols*NStr[nn],1))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,Ncols+1):
                        PSDcols[loc,0] = ops.sectionForce(int(j+3e4+1e5*i),1,1)
                        loc += 1
                # *************************************************************************************
                # LIVE LOAD CASE
                # *************************************************************************************
                ops.wipe()
                OPSDMOD.ModelGen(ndm,ndf)
                OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],'yes')
                OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],'yes')
                OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,'yes')
                OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],'yes')
                OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                # Interface elements generation for foundation flexibility considerations.
                # =========================================================================
                OPSDMOD.FoundFlexMaterials(NbaysX[mm],NbaysZ,XbayL,ZbayL,Ss,Vso[oo],gamma_soil,nu_soil,B,L,Re,D,omega_soil,'grav')
                # Zero-Length elements creation for connecting base nodes.
                OPSDMOD.FoundFlexZLElements(NbaysX[mm],NbaysZ,XbayL,ZbayL,B,L,Re)
                
                OPSDMOD.GeomTransGen(ColTransfType)
                OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                               ColDims,BeamEffFact,ColEffFact,Ec,fy)
                [Wx,MassInputMatr] = \
                    OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,'yes')
                    
                # LIVE-LOAD-CASE ANALYSIS.
                # ============================
                OPSDMOD.LiveLoadGen(NbaysX[mm],NbaysZ,NStr[nn],XbayL,ZbayL,Ql,Qlr)
                ops.system('BandGeneral')
                ops.constraints('Transformation')
                ops.numberer('RCM')
                ops.test('NormDispIncr',1.0e-5,100)
                ops.algorithm('Newton')
                ops.integrator('LoadControl', 1)
                ops.analysis('Static')
                ops.analyze(1)
                # Vertical reactions Calculation for verification.
                # ------------------------------------------------
                ops.reactions()
                YReactL = 0
                for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                    YReactL += ops.nodeReaction(int(i+1),2)
                    
                # Storing moments for beams and axial forces for columns
                # ======================================================
                # Beams on X
                # -----------
                MLXbeams = np.zeros((NbeamsX*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsX+1):
                        MLXbeams[loc,0] = ops.sectionForce(int(j+1e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mz,My,T
                        MLXbeams[loc,1] = ops.sectionForce(int(j+1e4+1e5*i),int((N-1)/2 + 1),2)
                        MLXbeams[loc,2] = ops.sectionForce(int(j+1e4+1e5*i),N,2)
                        loc += 1
                # Beams on Z
                # -----------
                MLZbeams = np.zeros((NbeamsZ*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsZ+1):
                        MLZbeams[loc,0] = ops.sectionForce(int(j+2e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mz,My,T
                        MLZbeams[loc,1] = ops.sectionForce(int(j+2e4+1e5*i),int((N-1)/2 + 1),2)
                        MLZbeams[loc,2] = ops.sectionForce(int(j+2e4+1e5*i),N,2)
                        loc += 1
                # Columns
                # -------
                PLcols = np.zeros((Ncols*NStr[nn],1))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,Ncols+1):
                        PLcols[loc,0] = ops.sectionForce(int(j+3e4+1e5*i),1,1)
                        loc += 1
                        
                # print('Live Load Case finished')
                # *************************************************************************************
                # LATERAL EQUIVALENT LOAD CASE - X DIRECTION
                # *************************************************************************************
                ops.wipe()
                OPSDMOD.ModelGen(ndm,ndf)
                OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],'yes')
                OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],'yes')
                OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,'yes')
                OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],'yes')
                OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                # Interface elements generation for foundation flexibility considerations.
                # =========================================================================
                OPSDMOD.FoundFlexMaterials(NbaysX[mm],NbaysZ,XbayL,ZbayL,Ss,Vso[oo],gamma_soil,nu_soil,B,L,Re,D,omega_soil,'lat')
                # Zero-Length elements creation for connecting base nodes.
                OPSDMOD.FoundFlexZLElements(NbaysX[mm],NbaysZ,XbayL,ZbayL,B,L,Re)
                
                OPSDMOD.GeomTransGen(ColTransfType)
                OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                               ColDims,BeamEffFact,ColEffFact,Ec,fy)
                [Wx,MassInputMatr] = \
                    OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,'yes')
                
                # Lateral load determination, Fx - X DIRECTION
                # --------------------------------------------
                # Already calculated at the dimensioning stage.
                        
                # EQUIVALENT LATERAL FORCE PROCEDURE - FORCES APPLICATION
                # =======================================================
                ops.loadConst('-time',0.0)   # Mantaining the gravity loads constant along the
                                         # ELF application.
                
                # Force Distribution from SSI considerations previous calculations.
                # -----------------------------------------------------------------
                OPSDMOD.ELFPForceGen(dirs[0],FxFSSI[:,0],'yes')
                
                # EQUIVALENT LATERAL FORCE CASE ANALYSIS - STATIC
                # ================================================
                ops.wipeAnalysis()
                ops.system('BandGeneral')
                ops.constraints('Transformation')
                ops.numberer('RCM')
                ops.test('NormDispIncr', 1.0e-5, 100)
                ops.algorithm('Newton')
                ops.integrator('LoadControl', 1)
                ops.analysis('Static')
                ops.analyze(1)
                # Lateral reactions for verification of Base Shear.
                # -------------------------------------------------
                ops.reactions()
                LatReactX = 0
                for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                    if dirs[0] in ('x','X'):
                        LatReactX += ops.nodeReaction(int(i+1),1)
                    elif dirs[0] in ('z','Z'):
                        LatReactX += ops.nodeReaction(int(i+1),3)
                
                # Storing moments for beams and axial forces for columns
                # ======================================================
                # Beams on X
                # -----------
                MELFXXbeams = np.zeros((NbeamsX*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsX+1):
                        MELFXXbeams[loc,0] = ops.sectionForce(int(j+1e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mx,My,T
                        MELFXXbeams[loc,1] = ops.sectionForce(int(j+1e4+1e5*i),int((N-1)/2 + 1),2)
                        MELFXXbeams[loc,2] = ops.sectionForce(int(j+1e4+1e5*i),N,2)
                        loc += 1
                # Beams on Z
                # -----------
                MELFXZbeams = np.zeros((NbeamsZ*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsZ+1):
                        MELFXZbeams[loc,0] = ops.sectionForce(int(j+2e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mx,My,T
                        MELFXZbeams[loc,1] = ops.sectionForce(int(j+2e4+1e5*i),int((N-1)/2 + 1),2)
                        MELFXZbeams[loc,2] = ops.sectionForce(int(j+2e4+1e5*i),N,2)
                        loc += 1
                # Columns
                # -------
                PELFXcols = np.zeros((Ncols*NStr[nn],1))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,Ncols+1):
                        PELFXcols[loc,0] = ops.sectionForce(int(j+3e4+1e5*i),1,1)
                        loc += 1
                
                # print('Lateral Load Case - X Direction finished')
                # *************************************************************************************
                # LATERAL EQUIVALENT LOAD CASE - Z DIRECTION
                # *************************************************************************************
                ops.wipe()
                OPSDMOD.ModelGen(ndm,ndf)
                OPSDMOD.NodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],'yes')
                OPSDMOD.MastNodeGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,StoryH,NStr[nn],'yes')
                OPSDMOD.SPConstGen(NbaysX[mm],NbaysZ,'yes')
                OPSDMOD.MPConstGen(NbaysX[mm],NbaysZ,NStr[nn],'yes')
                OPSDMOD.MatGenRCB(fpc,Ec,fy,E0,bsteel)
                # Interface elements generation for foundation flexibility considerations.
                # =========================================================================
                OPSDMOD.FoundFlexMaterials(NbaysX[mm],NbaysZ,XbayL,ZbayL,Ss,Vso[oo],gamma_soil,nu_soil,B,L,Re,D,omega_soil,'lat')
                # Zero-Length elements creation for connecting base nodes.
                OPSDMOD.FoundFlexZLElements(NbaysX[mm],NbaysZ,XbayL,ZbayL,B,L,Re)
                
                OPSDMOD.GeomTransGen(ColTransfType)
                OPSDMOD.ElementGen(NbaysX[mm],NbaysZ,XbayL,ZbayL,NStr[nn],StoryH,XBDims,ZBDims,\
                                               ColDims,BeamEffFact,ColEffFact,Ec,fy)
                [Wx,MassInputMatr] = \
                    OPSDMOD.LumpedMassGen(NbaysX[mm],NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr[nn],StoryH,Qsd,'yes')
                
                # Lateral load determination, Fx - Z DIRECTION
                # --------------------------------------------
                # Already calculated at the dimensioning stage.
                        
                # EQUIVALENT LATERAL FORCE PROCEDURE - FORCES APPLICATION
                # =======================================================
                ops.loadConst('-time',0.0)   # Mantaining the gravity loads constant along the
                                         # ELF application.
                
                # Force Distribution from SSI considerations previous calculations.
                # -----------------------------------------------------------------            
                OPSDMOD.ELFPForceGen(dirs[1],FxFSSI[:,1],'yes')
                
                # EQUIVALENT LATERAL FORCE CASE ANALYSIS - STATIC
                # ================================================
                ops.wipeAnalysis()
                ops.system('BandGeneral')
                ops.constraints('Transformation')
                ops.numberer('RCM')
                ops.test('NormDispIncr', 1.0e-5, 100)
                ops.algorithm('Newton')
                ops.integrator('LoadControl', 1)
                ops.analysis('Static')
                ops.analyze(1)
                # Lateral reactions for verification of Base Shear.
                # -------------------------------------------------
                ops.reactions()
                LatReactZ = 0
                for i in range((NbaysX[mm]+1)*(NbaysZ+1)):
                    if dirs[1] in ('x','X'):
                        LatReactZ += ops.nodeReaction(int(i+1),1)
                    elif dirs[1] in ('z','Z'):
                        LatReactZ += ops.nodeReaction(int(i+1),3)
                
                # Storing moments for beams and axial forces for columns
                # ======================================================
                # Beams on X
                # -----------
                MELFZXbeams = np.zeros((NbeamsX*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsX+1):
                        MELFZXbeams[loc,0] = ops.sectionForce(int(j+1e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mx,My,T
                        MELFZXbeams[loc,1] = ops.sectionForce(int(j+1e4+1e5*i),int((N-1)/2 + 1),2)
                        MELFZXbeams[loc,2] = ops.sectionForce(int(j+1e4+1e5*i),N,2)
                        loc += 1
                # Beams on Z
                # -----------
                MELFZZbeams = np.zeros((NbeamsZ*NStr[nn],3))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbeamsZ+1):
                        MELFZZbeams[loc,0] = ops.sectionForce(int(j+2e4+1e5*i),1,2)    # Because the output DOFs for a section are P,Mx,My,T
                        MELFZZbeams[loc,1] = ops.sectionForce(int(j+2e4+1e5*i),int((N-1)/2 + 1),2)
                        MELFZZbeams[loc,2] = ops.sectionForce(int(j+2e4+1e5*i),N,2)
                        loc += 1
                # Columns
                # -------
                PELFZcols = np.zeros((Ncols*NStr[nn],1))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,Ncols+1):
                        PELFZcols[loc,0] = ops.sectionForce(int(j+3e4+1e5*i),1,1)
                        loc += 1
                
                # ===========================================================================
                # LOAD COMBINATIONS FOR THE ANALYSED LOAD CASES 
                # ===========================================================================
                # There are the following load combinations:
                # U1 = 1.20*D + 1.00*L + 1.00*rho(EX+0.30EZ)
                # U2 = 1.20*D + 1.00*L + 1.00*rho(0.30EX+EZ)
                # U3 = 1.20*D + 1.00*L - 1.00*rho(Ex+0.30EZ)
                # U4 = 1.20*D + 1.00*L - 1.00*rho(0.30Ex+EZ)
                #        Where rho is the redundancy factor.
                rho = ASCE716.RedFact(SDC)
                # Combianations for Xbeams
                # ------------------------
                MuBX = np.zeros((NbeamsX*NStr[nn],6))   # [kN-m]
                for i in range(NbeamsX*NStr[nn]):
                    MuBX[i,0] = max(1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] + 1.00*rho*(MELFXXbeams[i,0]+0.30*MELFZXbeams[i,0]),
                                    1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] - 1.00*rho*(MELFXXbeams[i,0]+0.30*MELFZXbeams[i,0]),\
                                        1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] + 1.00*rho*(0.30*MELFXXbeams[i,0]+MELFZXbeams[i,0]),\
                                            1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] - 1.00*rho*(0.30*MELFXXbeams[i,0]+MELFZXbeams[i,0]))
                    MuBX[i,1] = min(1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] + 1.00*rho*(MELFXXbeams[i,0]+0.30*MELFZXbeams[i,0]),
                                    1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] - 1.00*rho*(MELFXXbeams[i,0]+0.30*MELFZXbeams[i,0]),\
                                        1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] + 1.00*rho*(0.30*MELFXXbeams[i,0]+MELFZXbeams[i,0]),\
                                            1.20*(MDXbeams[i,0]+MSDXbeams[i,0]) + 1.00*MLXbeams[i,0] - 1.00*rho*(0.30*MELFXXbeams[i,0]+MELFZXbeams[i,0]))
                    MuBX[i,2] = max(1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] + 1.00*rho*(MELFXXbeams[i,1]+0.30*MELFZXbeams[i,1]),
                                    1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] - 1.00*rho*(MELFXXbeams[i,1]+0.30*MELFZXbeams[i,1]),\
                                        1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] + 1.00*rho*(0.30*MELFXXbeams[i,1]+MELFZXbeams[i,1]),\
                                            1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] - 1.00*rho*(0.30*MELFXXbeams[i,1]+MELFZXbeams[i,1]))
                    MuBX[i,3] = min(1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] + 1.00*rho*(MELFXXbeams[i,1]+0.30*MELFZXbeams[i,1]),
                                    1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] - 1.00*rho*(MELFXXbeams[i,1]+0.30*MELFZXbeams[i,1]),\
                                        1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] + 1.00*rho*(0.30*MELFXXbeams[i,1]+MELFZXbeams[i,1]),\
                                            1.20*(MDXbeams[i,1]+MSDXbeams[i,1]) + 1.00*MLXbeams[i,1] - 1.00*rho*(0.30*MELFXXbeams[i,1]+MELFZXbeams[i,1]))
                    MuBX[i,4] = max(1.20*(MDXbeams[i,2]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] + 1.00*rho*(MELFXXbeams[i,2]+0.30*MELFZXbeams[i,2]),
                                    1.20*(MDXbeams[i,2]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] - 1.00*rho*(MELFXXbeams[i,2]+0.30*MELFZXbeams[i,2]),\
                                        1.20*(MDXbeams[i,2]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] + 1.00*rho*(0.30*MELFXXbeams[i,2]+MELFZXbeams[i,2]),\
                                            1.20*(MDXbeams[i,1]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] - 1.00*rho*(0.30*MELFXXbeams[i,2]+MELFZXbeams[i,2]))
                    MuBX[i,5] = min(1.20*(MDXbeams[i,2]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] + 1.00*rho*(MELFXXbeams[i,2]+0.30*MELFZXbeams[i,2]),
                                    1.20*(MDXbeams[i,2]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] - 1.00*rho*(MELFXXbeams[i,2]+0.30*MELFZXbeams[i,2]),\
                                        1.20*(MDXbeams[i,2]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] + 1.00*rho*(0.30*MELFXXbeams[i,2]+MELFZXbeams[i,2]),\
                                            1.20*(MDXbeams[i,1]+MSDXbeams[i,2]) + 1.00*MLXbeams[i,2] - 1.00*rho*(0.30*MELFXXbeams[i,2]+MELFZXbeams[i,2]))
                        
                # Combianations for Zbeams
                # ------------------------
                MuBZ = np.zeros((NbeamsZ*NStr[nn],6))   # [kN-m]
                for i in range(NbeamsZ*NStr[nn]):
                    MuBZ[i,0] = max(1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] + 1.00*rho*(MELFXZbeams[i,0]+0.30*MELFZZbeams[i,0]),
                                    1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] - 1.00*rho*(MELFXZbeams[i,0]+0.30*MELFZZbeams[i,0]),\
                                        1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] + 1.00*rho*(0.30*MELFXZbeams[i,0]+MELFZZbeams[i,0]),\
                                            1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] - 1.00*rho*(0.30*MELFXZbeams[i,0]+MELFZZbeams[i,0]))
                    MuBZ[i,1] = min(1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] + 1.00*rho*(MELFXZbeams[i,0]+0.30*MELFZZbeams[i,0]),
                                    1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] - 1.00*rho*(MELFXZbeams[i,0]+0.30*MELFZZbeams[i,0]),\
                                        1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] + 1.00*rho*(0.30*MELFXZbeams[i,0]+MELFZZbeams[i,0]),\
                                            1.20*(MDZbeams[i,0]+MSDZbeams[i,0]) + 1.00*MLZbeams[i,0] - 1.00*rho*(0.30*MELFXZbeams[i,0]+MELFZZbeams[i,0]))
                    MuBZ[i,2] = max(1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] + 1.00*rho*(MELFXZbeams[i,1]+0.30*MELFZZbeams[i,1]),
                                    1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] - 1.00*rho*(MELFXZbeams[i,1]+0.30*MELFZZbeams[i,1]),\
                                        1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] + 1.00*rho*(0.30*MELFXZbeams[i,1]+MELFZZbeams[i,1]),\
                                            1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] - 1.00*rho*(0.30*MELFXZbeams[i,1]+MELFZZbeams[i,1]))
                    MuBZ[i,3] = min(1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] + 1.00*rho*(MELFXZbeams[i,1]+0.30*MELFZZbeams[i,1]),
                                    1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] - 1.00*rho*(MELFXZbeams[i,1]+0.30*MELFZZbeams[i,1]),\
                                        1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] + 1.00*rho*(0.30*MELFXZbeams[i,1]+MELFZZbeams[i,1]),\
                                            1.20*(MDZbeams[i,1]+MSDZbeams[i,1]) + 1.00*MLZbeams[i,1] - 1.00*rho*(0.30*MELFXZbeams[i,1]+MELFZZbeams[i,1]))
                    MuBZ[i,4] = max(1.20*(MDZbeams[i,2]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] + 1.00*rho*(MELFXZbeams[i,2]+0.30*MELFZZbeams[i,2]),
                                    1.20*(MDZbeams[i,2]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] - 1.00*rho*(MELFXZbeams[i,2]+0.30*MELFZZbeams[i,2]),\
                                        1.20*(MDZbeams[i,2]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] + 1.00*rho*(0.30*MELFXZbeams[i,2]+MELFZZbeams[i,2]),\
                                            1.20*(MDZbeams[i,1]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] - 1.00*rho*(0.30*MELFXZbeams[i,2]+MELFZZbeams[i,2]))
                    MuBZ[i,5] = min(1.20*(MDZbeams[i,2]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] + 1.00*rho*(MELFXZbeams[i,2]+0.30*MELFZZbeams[i,2]),
                                    1.20*(MDZbeams[i,2]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] - 1.00*rho*(MELFXZbeams[i,2]+0.30*MELFZZbeams[i,2]),\
                                        1.20*(MDZbeams[i,2]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] + 1.00*rho*(0.30*MELFXZbeams[i,2]+MELFZZbeams[i,2]),\
                                            1.20*(MDZbeams[i,1]+MSDZbeams[i,2]) + 1.00*MLZbeams[i,2] - 1.00*rho*(0.30*MELFXZbeams[i,2]+MELFZZbeams[i,2]))
                        
                # Combianations for Columns
                # ------------------------
                PuCols = np.zeros((Ncols*NStr[nn],1))    # [kN] axial load in columns.
                for i in range(Ncols*NStr[nn]):
                    PuCols[i,0] = min(1.20*(PDcols[i,0]+PSDcols[i,0]) + 1.00*PLcols[i,0] + 1.00*rho*(PELFXcols[i,0]+0.30*PELFZcols[i,0]),
                                    1.20*(PDcols[i,0]+PSDcols[i,0]) + 1.00*PLcols[i,0] - 1.00*rho*(PELFXcols[i,0]+0.30*PELFZcols[i,0]),\
                                        1.20*(PDcols[i,0]+PSDcols[i,0]) + 1.00*PLcols[i,0] + 1.00*rho*(0.30*PELFXcols[i,0]+PELFZcols[i,0]),\
                                            1.20*(PDcols[i,0]+PSDcols[i,0]) + 1.00*PLcols[i,0] - 1.00*rho*(0.30*PELFXcols[i,0]+PELFZcols[i,0]))
                
                
                # ***********************************************************************************************************************************
                # ***********************************************************************************************************************************
                # REINFORCING OF BEAMS
                # ***********************************************************************************************************************************
                # ***********************************************************************************************************************************
                # BEAMS IN X DIRECTION
                # =====================
                
                # print('Starting Beam reinforcing....')
                
                MuBXDes = abs(MuBX/1000)     # [MN-m] changing units for using in other function. 
                                             # absolute value because the index into the object determines the sign within.
                XBreinf = np.zeros((NbaysX[mm]*(NbaysZ+1)*NStr[nn],4))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbaysZ+2):
                        for k in range(1,NbaysX[mm]+1):
                            if j == 1 or j == NbaysZ+1:
                                beamloc = 'edge'
                            else:
                                beamloc = 'interior'
                            XBreinf[loc,:] = \
                                DLBDefs.SpecialBeamDesign(abs(fpc),fy,E0,etlim,XBDims[loc,0],XBDims[loc,1],\
                                                          MuBXDes[loc,0],MuBXDes[loc,1],MuBXDes[loc,2],MuBXDes[loc,3],MuBXDes[loc,4],MuBXDes[loc,5],\
                                                              transtype,XbayL,ZbayL,beamloc,hslab=0.25)
                            loc += 1
                            
                # BEAMS IN Z DIRECTION
                # =====================
                MuBZDes = abs(MuBZ/1000)     # [MN-m] changing units for using in other function. 
                                             # absolute value because the order into the object determines the sign within.
                ZBreinf = np.zeros((NbaysZ*(NbaysX[mm]+1)*NStr[nn],4))
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NbaysX[mm]+2):
                        for k in range(1,NbaysZ+1):
                            if j == 1 or j == NbaysX[mm]+1:
                                beamloc = 'edge'
                            else:
                                beamloc = 'interior'
                            ZBreinf[loc,:] = \
                                DLBDefs.SpecialBeamDesign(abs(fpc),fy,E0,etlim,ZBDims[loc,0],ZBDims[loc,1],\
                                                          MuBZDes[loc,0],MuBZDes[loc,1],MuBZDes[loc,2],MuBZDes[loc,3],MuBZDes[loc,4],MuBZDes[loc,5],\
                                                              transtype,XbayL,ZbayL,beamloc,hslab=0.25)
                            loc += 1
                
                # NOMINAL MOMENTS FROM BEAMS ON NODES FOR COLUMN DESIGN.
                # ======================================================
                NAX = NbaysX[mm] + 1
                NAZ = NbaysZ + 1
                MnbXdir = np.zeros((Ncols*NStr[nn],1))   # [MN-m] nominal moments at nodes
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NAZ+1):
                        for k in range(1,NbaysX[mm]+1):
                            Nodei = k + (j-1)*NAX + (i-1)*Ncols - 1      # the order in which the list of columns is stored, not the node numbering.
                            Nodej = k + 1 + (j-1)*NAX + (i-1)*Ncols - 1
                            (Mnipos,Mnineg,Mnjpos,Mnjneg) = \
                                DLBDefs.BeamEndsMn(abs(fpc),fy,E0,XBDims[loc,0],XBDims[loc,1],\
                                                   XBreinf[loc,0],XBreinf[loc,1],XBreinf[loc,2],XBreinf[loc,3],transtype)
                            MnbXdir[Nodei,0] += Mnipos
                            MnbXdir[Nodej,0] += Mnjneg
                            loc += 1
                            # print(str(loc)+' X-Beams analized!')
                MnbZdir = np.zeros((Ncols*NStr[nn],1))   # [MN-m] nominal moments at nodes
                loc = 0
                for i in range(1,NStr[nn]+1):
                    for j in range(1,NAX+1):
                        for k in range(1,NbaysZ+1):
                            Nodei = j + (k-1)*NAX + (i-1)*Ncols - 1
                            Nodej = j + (k-1)*NAX + NAX + (i-1)*Ncols - 1
                            (Mnipos,Mnineg,Mnjpos,Mnjneg) = \
                                DLBDefs.BeamEndsMn(abs(fpc),fy,E0,XBDims[loc,0],XBDims[loc,1],\
                                                   ZBreinf[loc,0],ZBreinf[loc,1],ZBreinf[loc,2],ZBreinf[loc,3],transtype)
                            MnbZdir[Nodei,0] += Mnipos
                            MnbZdir[Nodej,0] += Mnjneg
                            loc += 1
                            
                # ================================            
                # Design of column reinforcemwnt.
                # ================================
                # print('Starting Columns reinforcing....')
                loc = 0
                rhomincol = 0.01
                # rhomaxcol = 0.04
                Colreinf = np.zeros((Ncols*NStr[nn],10))
                for i in range(1,Ncols*NStr[nn]+1):
                    
                    # X Direction
                    # ------------
                    # bcolx = ColDims[loc,0]
                    hcolx = ColDims[loc,1]
                    # bcorex = bcolx - 2*cover
                    hcorex = hcolx - 2*cover
                    edx = hcorex/4
                    dx = [cover,cover+edx,cover+2*edx,cover+3*edx,cover+4*edx]
                    # Z Direction
                    # ------------
                    # bcolz = ColDims[loc,1]
                    hcolz = ColDims[loc,0]
                    # bcorez = bcolz - 2*cover
                    hcorez = hcolz - 2*cover
                    edz = hcorez/3          # because there are 4 layers in this direction!!
                    dz = [cover,cover+edz,cover+2*edz,cover+3*edz]
                    
                    
                    Ascol = rhomincol*ColDims[loc,0]*ColDims[loc,1]
                    Mnx = 0.0
                    Mnz = 0.0
                    while 2*Mnx < 1.20*MnbXdir[loc,0] or 2*Mnz < 1.20*MnbZdir[loc,0]:
                        Asz = [5/14*Ascol,2/14*Ascol,2/14*Ascol,5/14*Ascol]
                        Asx = [4/14*Ascol,2/14*Ascol,2/14*Ascol,2/14*Ascol,4/14*Ascol]
                        # X Direction
                        # ------------
                        (Mnx,Pn,phi) =\
                            DLBDefs.ColAnalysisMn(abs(fpc),fy,E0,etlim,abs(PuCols[loc,0]),ColDims[loc,1],ColDims[loc,0],\
                                                  dx,Asx,5,transtype)
                        # Z Direction
                        # ------------
                        (Mnz,Pn,phi) =\
                            DLBDefs.ColAnalysisMn(abs(fpc),fy,E0,etlim,abs(PuCols[loc,0]),ColDims[loc,0],ColDims[loc,1],\
                                                  dz,Asz,4,transtype)
                        Ascol += 0.0001
                        
                    Colreinf[loc,:5] = 4/14*Ascol, 4/14*Ascol, 2/14*Ascol, 2/14*Ascol, 2/14*Ascol
                    Colreinf[loc,5:] = 4/14*Ascol, 4/14*Ascol, 2/14*Ascol, 2/14*Ascol, 2/14*Ascol
                    loc += 1
                    # print(str(loc)+' Cols with reinforcement')
                
              
                counter += 1
                print(f"{round(counter/total*100,2)} %  --  {counter} designed out of {total}")
                varpath = directory+'\\SSIDesign\\SSI'+\
                    str(NbaysX[mm])+'BayX'+str(NbaysZ)+'BayZ'+str(NStr[nn])+'FLRS'+str(Vso[oo])+'.pickle'
                with open(varpath,'wb') as f:
                    pickle.dump([ColDims,Colreinf,XBDims,XBreinf,ZBDims,ZBreinf,DiSSI_h,\
                                  FxD,FxDSSI,FxF,FxFSSI,T,Ta,Tmax,VgorD,VgorF,VxD,VxF,\
                                      DeltaVF,DeltaVD,DeltaVF,Da,hxplot],f)
                        
    return ColDims,Colreinf,XBDims,XBreinf,ZBDims,ZBreinf,DiSSI_h,FxD,FxDSSI,YReaction,\
        FxF,FxFSSI,T,Ta,Tmax,VgorD,VgorF,VxD,VxF,DeltaVF,DeltaVD,DeltaVF,Da,hxplot