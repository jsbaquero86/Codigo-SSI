# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 12:54:41 2020

@author: sebas
"""

"""
    Definitions (Functions) for the automatic generation of the general model
    (nodes, elements, masses and gravity loads) are written. Almost all of the 
    definitions assumes that the building is completely regular.
    
"""
import openseespy.opensees as ops


def POTargetDisp(Ss, S1, Te, W, Vs30, NStr, BldTypCo, BldTypCm):
    # Calculates the target displacement according to ASCE41-17.

    import ASCE716Seismic as ASCE716
    from math import pi

    g = 9.81  # [m/s**2] acceleration of gravity in SI units.
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
    # a factor determination according to site Class. Section 7.4.3.3.2 ASCE41-13
    # ----------------------------------------------------------------------------
    SiteClass = ASCE716.detSiteClass(Vs30)
    if SiteClass in ('A', 'B'):
        a = 130
    elif SiteClass in ('C'):
        a = 90
    elif SiteClass in ('D', 'E', 'F'):
        a = 60
    # C1 factor determination from Equation 7-29 of ASCE41-13.
    # ---------------------------------------------------------
    (Sa, _, _, _) = ASCE716.Sa(Ss, S1, Vs30, Te, TL=8.0)
    Vy = 0.225*W          # Asumed to calculate Target Displacement
    mustr = Sa/(Vy/W)*Cm  # According to Eqn. 7-31
    C1 = 1 + (mustr - 1)/(a*Te**2)
    # --------------  For verification of max C1 value -----------------
    (Sa_02, _, _, _) = ASCE716.Sa(Ss, S1, Vs30, 0.2, TL=8.0)
    mustr_02 = Sa_02/(Vy/W)*Cm
    C1_02 = 1 + (mustr_02 - 1)/(a*0.20**2)
    if C1 > C1_02:
        C1 = C1_02
    if Te > 1.0:
        C1 = 1.0
    # C2 factor calculation according to eqn. 7-30 of ASCE41-13.
    # ----------------------------------------------------------
    C2 = 1 + 1/800*((mustr - 1)/Te)**2
    if Te > 0.70:
        C2 = 1.0

    dtarget = Co*C1*C2*Sa*Te**2/(4*pi**2)*g

    return dtarget


def POAnalysisProc(lon,lat,Vs30,SHL,NbaysX,NbaysZ,NStr,StoryH,R,Ie,BldTypCo,BldTypCm,T,W,
                   direction,flex = 'no',DispIncr=0,beta=0.05,TL=8.,Tf=6.0,Npts=500,Vy=0,Te=0):
    # This function generates all the necessary parameters and modifications
    # through the pushover analysis for minimizing the chances of convergence.
    # Also performs the analysis considering the step number entered by users.
    # Needs to introduce some data.

    # Previous definitions for analysis parameters.
    # ---------------------------------------------
    import numpy as np
    import ASCE4117
    
    # Definition of whether to use T or Te for calculations.
    if Te != 0:
        T = Te

    # Calculation of target displacement.
    # ----------------------------------
    # According to ASCE41-17, from 0 to 150% of targer displacement shall
    # be evaluated. Hence the 1.5 factor affecting POTargetDisp.
    Cm = ASCE4117.detCm(NStr, BldTypCm, T)
    
    targetfactor = 1.5
    TargetDisp = ASCE4117.TargetDisp(lon,lat,Vs30,SHL,NStr,BldTypCo,T,W,R,Ie,Cm,
               beta,TL,Tf,Npts,Vy,T)
    print(T)
    print(TargetDisp)
    print(TargetDisp*targetfactor)
    
    PPF = (NbaysX + 1)*(NbaysZ + 1)
    if DispIncr == 0:
        DispIncr = 0.001  # [m] increment of displacement
    CtrlNode = int(999+1e4+1e5*NStr)
    POdisp = []
    POforce = []

    if direction in ('x', 'X'):
        CtrlDOF = 1
    elif direction in ('z', 'Z'):
        CtrlDOF = 3

    ops.wipeAnalysis()
    # Initial Parameters of analysis.
    # --------------------------------
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')

    # Variable analysis parameters in case tha analysis fails in any step
    # due to unconvergence.
    # -------------------------------------------------------------------
    fract = 1/10  # fraction of displacement increment to use in case of analysis failure.
    iters = 100
    ops.integrator('DisplacementControl', CtrlNode, CtrlDOF, DispIncr)
    ops.test('NormDispIncr', 1.0e-4, iters)
    ops.algorithm('RaphsonNewton')
    ops.analysis('Static')

    # Storing initial state of structure.
    # ====================================
    # Displacement
    # -------------
    POdisp.append(ops.nodeDisp(CtrlNode, CtrlDOF))  # [m] initial displacement of control node.

    # Reactions
    # -------------
    ops.reactions()
    LatReact = 0
    if flex in ('Y','YES','Yes','yES','yes','y'):
        for j in range(PPF):
            LatReact += ops.nodeReaction(int(j+1), CtrlDOF)
    else:
        for j in range(PPF):
            LatReact += ops.nodeReaction(int(j+1+1e4), CtrlDOF)
    POforce.append(LatReact)  # [kN] initial base shear.

    # Predefining additional options of tester and algorithms to try if convergency
    # falures arrise during static nonlinear analysis.
    # ------------------------------------------------------------------------------
    options = [('ModifiedNewton', 'NormUnbalance'), ('ModifiedNewton', 'NormDispIncr'),
               ('ModifiedNewton', 'EnergyIncr'), \
               ('KrylovNewton', 'NormUnbalance'), ('KrylovNewton', 'NormDispIncr'), ('KrylovNewton', 'EnergyIncr'), \
               ('Broyden', 'NormUnbalance'), ('Broyden', 'NormDispIncr'), ('Broyden', 'EnergyIncr'), \
               ('NewtonLineSearch', 'NormUnbalance'), ('NewtonLineSearch', 'NormDispIncr'),
               ('NewtonLineSearch', 'EnergyIncr')]

    count = 0  # counter for locating current displacement, disp_current.
    killanal = 0  # [FLAG] flag to kill entire analysis because of convergency isues.

    print('-----------------------------------------------')
    print('  ***   Starting Static Nonlinear Analysis *** ')
    print('-----------------------------------------------')

    while POdisp[count] < targetfactor*TargetDisp:

        okPO = ops.analyze(1)

        if okPO != 0:  # First, I reduce analysis displacement increment to a predefined fraction.
            DispIncr_red = fract*DispIncr
            print(f' **Troubles** : No convergency - Trying reducing step interval!')
            print(f' **Notice** : Reducing displacement increment to a {round(fract*100, 2)} % - from {DispIncr} to {DispIncr_red} [m]')
            evolDispIncr_red = 0  # [s] evolution of reduced displacement increment.
            while evolDispIncr_red < DispIncr:
                okPO = ops.analyze(1)
                if okPO != 0:
                    print(f' **Troubles** : Trying to reach convergency with other algorithms and testers!')
                    select = 0  # counter for tester-algorithm combinations
                    while okPO != 0 and select < len(options):
                        ops.algorithm(options[select][0])
                        ops.test(options[select][1], 1.0e-4, iters)
                        okPO = ops.analyze(1)
                        select += 1
                    if okPO != 0 and select >= len(options):
                        print(f' **Error** : Not able to converge at {round(POdisp[count], 2)} [m] of total TargetDisp.')
                        killanal = 1
                        break
                    else:
                        ops.test('NormDispIncr', 1.0e-6, iters)
                        ops.algorithm('RaphsonNewton')
                        print(f' **Notice** : {options[select]} worked! Back to initial tester and algorithm.')
                        evolDispIncr_red += DispIncr_red
                        POdisp.append(ops.nodeDisp(CtrlNode, CtrlDOF))
                        count += 1
                        # Storing reactions from base nodes.
                        # ----------------------------------
                        ops.reactions()
                        LatReact = 0
                        if flex in ('Y','YES','Yes','yES','yes','y'):
                            for j in range(PPF):
                                LatReact += ops.nodeReaction(int(j+1), CtrlDOF)
                        else:
                            for j in range(PPF):
                                LatReact += ops.nodeReaction(int(j+1+1e4), CtrlDOF)
                        POforce.append(LatReact)
                        print(f'>>> Pushover analysis - {round(POdisp[count] / (targetfactor*TargetDisp)*100, 2)} % executed... ({round(POdisp[count], 2)} of {round(targetfactor*TargetDisp, 2)} [m])')
                else:
                    evolDispIncr_red += DispIncr_red
                    POdisp.append(ops.nodeDisp(CtrlNode, CtrlDOF))
                    count += 1
                    # Storing reactions from base nodes.
                    # ----------------------------------
                    ops.reactions()
                    LatReact = 0
                    if flex in ('Y','YES','Yes','yES','yes','y'):
                        for j in range(PPF):
                            LatReact += ops.nodeReaction(int(j+1), CtrlDOF)
                    else:
                        for j in range(PPF):
                            LatReact += ops.nodeReaction(int(j+1+1e4), CtrlDOF)
                    POforce.append(LatReact)
                    print(f'>>> Pushover analysis - {round(POdisp[count] / (targetfactor*TargetDisp)*100, 2)} % executed... ({round(POdisp[count], 2)} of {round(targetfactor*TargetDisp, 2)} [m])')

            if killanal != 0:
                break
            else:
                print(f' **Notice** : Going back to original displacement increment - from {DispIncr_red} to {DispIncr} [m]')
        else:
            POdisp.append(ops.nodeDisp(CtrlNode, CtrlDOF))
            count += 1
            # Storing reactions from base nodes.
            # ----------------------------------
            ops.reactions()
            LatReact = 0
            if flex in ('Y','YES','Yes','yES','yes','y'):
                for j in range(PPF):
                    LatReact += ops.nodeReaction(int(j+1), CtrlDOF)
            else:
                for j in range(PPF):
                    LatReact += ops.nodeReaction(int(j+1+1e4), CtrlDOF)
            POforce.append(LatReact)
            print(f'>>> Pushover analysis - {round(POdisp[count] / (targetfactor*TargetDisp)*100, 2)} % executed... ({round(POdisp[count], 2)} of {round(targetfactor*TargetDisp, 2)} [m])')

    if killanal == 0:
        print(f'>>> Pushover analysis 100% completed.')
    else:
        print(f' **Error** : Pushover analysis cancelled. No convergency at {POdisp[count]} [m] of total target Displacement!')

    POresults = np.zeros((len(POdisp), 2))
    for i in range(len(POdisp)):
        POresults[i, 0] = POdisp[i]
        POresults[i, 1] = -1.*POforce[i]

    return POresults, TargetDisp, targetfactor


def ModalAnalysis(NStr, B, L, Wx, flex='No'):
    # Returns modal parameters of the structure, including modal mass of the
    # fundamental mode.

    import numpy as np
    g = 9.81  # [m/s**2] acceleration of gravity.

    # ==========================================================================
    # MODAL ANALYSIS FOR DETERMINING FUNDAMENTAL PERIOD AND ITS DIRECTION.
    # ==========================================================================
    T = np.zeros((1, 3))
    w2 = ops.eigen('-genBandArpack', 5)
    for i in range(np.size(T, 1)):
        T[0, i] = 2*np.pi/(w2[i])**(0.5)

    # Determination of fundamental mode shape.
    # ----------------------------------------
    Mod1Disp = ops.nodeEigenvector(int(999+1e4+1e5*NStr), 1, 1)
    Mod2Disp = ops.nodeEigenvector(int(999+1e4+1e5*NStr), 2, 1)
    if abs(Mod1Disp) >= abs(Mod2Disp):
        TX = T[0, 0];
        fundamental = 1
        TZ = T[0, 1];
        nonfund = 2
    else:
        TX = T[0, 1];
        nonfund = 1
        TZ = T[0, 0];
        fundamental = 2
    T[0, 0], T[0, 1] = TX, TZ
    Tmaxreal = max(T[0, :2])

    Mast = np.zeros((2))

    # EFFECTIVE MASS OF FUNDAMENTAL MODE
    # ===================================
    # Mode 1 - EigenVector data extraction.
    # ------------------------------------
    Mod1 = np.zeros((len(Wx)*3, 1))
    if flex in ('Y', 'YES', 'Yes', 'yES', 'yes', 'y'):
        istory = 0
    else:
        istory = 1
    for i in range(len(Wx)):
        Mod1[i*3, 0] = ops.nodeEigenvector(int(999+1e4+1e5*istory), fundamental, 1)
        Mod1[i*3 + 1, 0] = ops.nodeEigenvector(int(999+1e4+1e5*istory), fundamental, 3)
        Mod1[i*3 + 2, 0] = ops.nodeEigenvector(int(999+1e4+1e5*istory), fundamental, 5)
        istory += 1

    # Mode Vector normalization.
    # --------------------------
    Mod1 = Mod1*(1/(np.sign(Mod1[0, 0])*max(abs(Mod1))))
    # Calculation of Modal effective mass on fundamental vibration mode.
    # ------------------------------------------------------------------
    MassM = np.zeros((len(Wx)*3, len(Wx)*3))  # Diagonal Mass matrix.
    for i in range(len(Wx)):
        MassM[i*3, i*3] = Wx[i]/g
        MassM[i*3 + 1, i*3 + 1] = Wx[i]/g
        MassM[i*3 + 2, i*3 + 2] = (Wx[i]/g)/12*(B**2 + L**2)

    Mast[0] = np.matmul(np.matmul(np.transpose(Mod1), MassM), Mod1)

    # EFFECTIVE MASS OF NONFUNDAMENTAL MODE
    # ===================================
    # Mode 2 - EigenVector data extraction.
    # ------------------------------------
    Mod2 = np.zeros((len(Wx)*3, 1))
    if flex in ('Y', 'YES', 'Yes', 'yES', 'yes', 'y'):
        istory = 0
    else:
        istory = 1
    for i in range(len(Wx)):
        Mod2[i*3, 0] = ops.nodeEigenvector(int(999+1e4+1e5*istory), nonfund, 1)
        Mod2[i*3 + 1, 0] = ops.nodeEigenvector(int(999+1e4+1e5*istory), nonfund, 3)
        Mod2[i*3 + 2, 0] = ops.nodeEigenvector(int(999+1e4+1e5*istory), nonfund, 5)
        istory += 1

    # Mode Vector normalization.
    # --------------------------
    Mod2 = Mod2*(1/(np.sign(Mod2[0, 0])*max(abs(Mod2))))
    # Calculation of Modal effective mass on fundamental vibration mode.
    # ------------------------------------------------------------------
    Mast[1] = np.matmul(np.matmul(np.transpose(Mod2), MassM), Mod2)

    return T, Tmaxreal, Mast


def RayleighDamping(Ti,Tj,lambi,lambj):
    # Rayleigh proportinoal damping definition.

    from math import pi

    # Calculation of frequencies
    # --------------------------
    wi = 2*pi/Ti
    wj = 2*pi/Tj

    # Calculation of proportion paramenters.
    # --------------------------------------
    alphaM = 2*wi*wj*(wi*lambj - wj*lambi)/(wi**2 - wj**2)
    betaK = 2*(wi*lambi - wj*lambj)/(wi**2 - wj**2)

    ops.rayleigh(alphaM, betaK, 0.0, 0.0)


def THAnalysisProc(FilePath, dt_rec, NStr, andirection, gammaNw=0.50, betaNw=0.25, dt_an=0, fract=0.25):
    # Definition to excecute a Time-History analysis with uniform exitation at
    # the base of the structure.
    import numpy as np
    g = 9.81  # [m/s**2] acceleration of gravity in SI.

    THtsTag = 6
    THpatTag = 6

    # TIME SERIES FOR UNIFORM EXCITATION
    # ----------------------------------
    record = np.loadtxt(FilePath, delimiter='\n')
    NPts = len(record)
    ops.timeSeries('Path', THtsTag, '-dt', dt_rec, '-values', *record, '-factor', g)

    # PATTERN FOR UNIFROM EXCITATION
    # ------------------------------
    if andirection in ('x', 'X'):
        direction = 1
    elif andirection in ('z', 'Z'):
        direction = 3

    # pattern('UniformExcitation', patternTag, dir, '-accel', accelSeriesTag)
    ops.pattern('UniformExcitation', THpatTag, direction, '-accel', THtsTag)

    # ANALYSIS PARAMETERS TO START ANALYSIS.
    # =======================================
    ops.wipeAnalysis()
    # Initial Parameters of analysis.
    # --------------------------------
    ops.constraints('Transformation')
    ops.numberer('RCM')
    ops.system('BandGeneral')
    # Variable analysis parameters in case tha analysis fails in any step
    # due to unconvergence.
    # -------------------------------------------------------------------
    # fract = 1/4   # fraction of time interval to modify in case of analysis failure.
    ops.integrator('Newmark', gammaNw, betaNw)
    iters = 100
    ops.test('NormDispIncr', 1.0e-4, iters)
    ops.algorithm('RaphsonNewton')
    ops.analysis('Transient')

    duration = NPts*dt_rec  # [s] duration of record according to Npts and dt.
    t_current = [ops.getTime()]  # [s] time in the domain. If loadConst() function
    #     has been called in general scope, it will be 0.0
    if dt_an == 0:
        dt_an = dt_rec  # [s] time interval of analysis. Initial the same as
        #     of the record.

    # Predefining dictionary to store displacement responses from floor Master Nodes.
    # -------------------------------------------------------------------------------
    floors = ['Str-' + str(xx) for xx in range(1, NStr + 1)]
    emptylist = [[0] for xx in range(NStr)]
    THdisp = dict(zip(floors, emptylist))

    # Predefining additional options of tester and algorithms to try if convergency
    # failures arrise during transient analysis.
    options = [('ModifiedNewton', 'NormUnbalance'), ('ModifiedNewton', 'NormDispIncr'),
               ('ModifiedNewton', 'EnergyIncr'), \
               ('KrylovNewton', 'NormUnbalance'), ('KrylovNewton', 'NormDispIncr'), ('KrylovNewton', 'EnergyIncr'), \
               ('Broyden', 'NormUnbalance'), ('Broyden', 'NormDispIncr'), ('Broyden', 'EnergyIncr'), \
               ('NewtonLineSearch', 'NormUnbalance'), ('NewtonLineSearch', 'NormDispIncr'),
               ('NewtonLineSearch', 'EnergyIncr')]

    count = 0  # counter for locating current time, t_current.
    killanal = 0  # [FLAG] flag to kill entire analysis because of convergency isues.

    print('--------------------------------------------')
    print('  ***  Starting Time - History Analysis *** ')
    print('--------------------------------------------')

    while t_current[count] < duration:

        okTH = ops.analyze(1, dt_an)

        if okTH != 0:  # First, I reduce analysis time increment to a predefined fraction.
            dt_red = fract*dt_rec
            print(f' **Troubles** : No convergency - Trying reducing time interval!')
            print(f' **Notice** : Reducing time interval to a {round(fract*100, 2)} % - from {dt_an} to {dt_red} [s]')
            evoldt_red = 0  # [s] evolution of reduced time increment.
            while evoldt_red < dt_an:
                okTH = ops.analyze(1, dt_red)
                if okTH != 0:
                    print(f' **Troubles** : Trying to reach convergency with other algorithms and testers!')
                    select = 0  # counter for tester-algorithm combinations
                    while okTH != 0 and select < len(options):
                        ops.algorithm(options[select][0])
                        ops.test(options[select][1], 1.0e-4, iters)
                        okTH = ops.analyze(1, dt_red)
                        select += 1
                    if okTH != 0 and select >= len(options):
                        print(f' **Error** : Not able to converge at {round(t_current[count], 2)} [s] of record time.')
                        killanal = 1
                        break
                    else:
                        ops.test('NormDispIncr', 1.0e-4, iters)
                        ops.algorithm('RaphsonNewton')
                        print(f' **Notice** : {options[select]} worked! Back to initial tester and algorithm.')
                        evoldt_red += dt_red
                        t_current.append(ops.getTime())
                        count += 1
                        # Storing displacements from master nodes in each floor.
                        # ------------------------------------------------------
                        for m in range(NStr):
                            THdisp[floors[m]].append(ops.nodeDisp(int(999+1e4+1e5*(m+1)), direction))
                        print(f'>>> Time-History analysis - {round(t_current[count]/duration*100, 2)} % excecuted... ({round(t_current[count], 2)} [s] of {round(duration, 2)} [s])')
                else:
                    evoldt_red += dt_red
                    t_current.append(ops.getTime())
                    count += 1
                    # Storing displacements from master nodes in each floor.
                    # ------------------------------------------------------
                    for m in range(NStr):
                        THdisp[floors[m]].append(ops.nodeDisp(int(999+1e4+1e5*(m+1)), direction))
                    print(f'>>> Time-History analysis - {round(t_current[count]/duration*100, 2)} % executed... ({round(t_current[count], 2)} [s] of {round(duration, 2)} [s])')

            if killanal != 0:
                break
            else:
                print(f' **Notice** : Going back to original time increment - from {dt_red} to {dt_an} [s]')
        else:
            t_current.append(ops.getTime())
            count += 1
            # Storing displacements from master nodes in each floor.
            # ------------------------------------------------------
            for m in range(NStr):
                THdisp[floors[m]].append(ops.nodeDisp(int(999+1e4+1e5*(m+1)), direction))
            print(
                f'>>> Time-History analysis - {round(t_current[count]/duration*100, 2)} % executed... ({round(t_current[count], 2)} [s] of {round(duration, 2)} [s])')

    if killanal == 0:
        print(f'>>> Time-History analysis 100% completed.')
    else:
        print(f' **Error** : Time-HIstory analysis cancelled. No convergency at {t_current[count]} [s] of record time!')

    return THdisp, t_current, floors
    
