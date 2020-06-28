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


def FiberSectGen(b,h,rec,As,Aprimas,As14,As24,As34,NbarsBot,NbarsTop,Nbarsh14,Nbarsh24,Nbarsh34,secTag,GJ=1.0e10):
    """
    asf

    Parameters
    ----------
    b : TYPE
        DESCRIPTION.
    h : TYPE
        DESCRIPTION.
    rec : TYPE
        DESCRIPTION.
    As : TYPE
        DESCRIPTION.
    Aprimas : TYPE
        DESCRIPTION.
    As14 : TYPE
        DESCRIPTION.
    As24 : TYPE
        DESCRIPTION.
    As34 : TYPE
        DESCRIPTION.
    NbarsBot : TYPE
        DESCRIPTION.
    NbarsTop : TYPE
        DESCRIPTION.
    Nbarsh14 : TYPE
        DESCRIPTION.
    Nbarsh24 : TYPE
        DESCRIPTION.
    Nbarsh34 : TYPE
        DESCRIPTION.
    secTag : TYPE
        DESCRIPTION.
    GJ : TYPE, optional
        DESCRIPTION. The default is 1.0e10.

    Returns
    -------
    None.

    """

    # Generates fiber section for doubly layered reinforced concrete rectangular
    # section.
    aStlFiberTop = Aprimas/NbarsTop
    aStlFiberBot = As/NbarsBot
    aStlFiber14 = As14/Nbarsh14
    aStlFiber24 = As24/Nbarsh24
    aStlFiber34 = As34/Nbarsh34
    # Previous calculations for geometry of section fibers.
    # -----------------------------------------------------
    coreb = b - 2*rec
    coreh = h - 2*rec
    # Fiber Section specification
    # ---------------------------
    ops.section('Fiber',secTag,'-GJ',GJ)
    # Concrete Core fibers
    # --------------------
    # patch('quad', matTag, numSubdivIJ, numSubdivJK, *crdsI, *crdsJ, *crdsK, *crdsL) --> Just for reference!
    ops.patch('quad',1,10,10,*[-coreh/2,-coreb/2],*[coreh/2,-coreb/2],*[coreh/2,coreb/2],*[-coreh/2, coreb/2])

    # Concrete cover fibers
    # ---------------------
    ops.patch('quad', 2, 10, 1,*[-h/2,coreb/2],*[h/2,coreb/2],*[h/2,b/2],*[-h/2,b/2])
    ops.patch('quad', 2, 10, 1,*[-h/2,-b/2],*[h/2,-b/2],*[h/2,-coreb/2],*[-h/2,-coreb/2])
    ops.patch('quad', 2, 1, 10,*[-h/2,-coreb/2],*[-coreh/2,-coreb/2],*[-coreh/2,coreb/2],*[-h/2,coreb/2])
    ops.patch('quad', 2, 1, 10,*[coreh/2,-coreb/2],*[h/2,-coreb/2],*[h/2,coreb/2],*[coreh/2,coreb/2])

    # Layers of reinforcement
    # -----------------------
    # Five layers are stablished to the reinforcement arrangement.
    # layer('straight', matTag, numFiber, areaFiber, *start, *end) --> Just for reference!
    ops.layer('straight', 3, NbarsTop, aStlFiberTop, *[coreh/2,coreb/2], *[coreh/2,-coreb/2])
    ops.layer('straight', 3, NbarsBot, aStlFiberBot, *[-coreh/2, coreb/2], *[-coreh/2, -coreb/2])
    ops.layer('straight', 3, Nbarsh14, aStlFiber14, *[-coreh/4, coreb/2], *[-coreh/4, -coreb/2])
    ops.layer('straight', 3, Nbarsh24, aStlFiber24, *[0., coreb/2], *[0., -coreb/2])
    ops.layer('straight', 3, Nbarsh34, aStlFiber34, *[coreh/4, coreb/2], *[coreh/4, -coreb/2])


def ModelGen(ndm, ndf):
    
    ops.model('basic', '-ndm', ndm, '-ndf', ndf)


def NodeGen(NbaysX,NbaysZ,XbayL,ZbayL,StoryH,NStr,flex='No'):

    import numpy as np
    NPX = NbaysX + 1          # number of nodes along X Global Axis.
    NPZ = NbaysZ + 1          # number of nodes along Z Global Axis.
    PPF = NPX*NPZ             # total number of points per flor.
    TNP = PPF*(NStr + 1)      # total number of points in the model.
    
    # Deatils of numbering will be described later, once all numbering tags are
    # coherent.
    
    # NODES RELATED TO BEAM-COLUMN JOINTS.
    # ====================================
    mat = np.zeros((TNP, 4))
    # Determination of integer tags for Nodes.
    # ----------------------------------------
    loc = 0
    for j in range(NStr + 1):
        pref = 1e5*j + 1e4
        for i in range(1, PPF + 1):
            mat[loc, 0] = i + pref
            loc = loc + 1
    # Determination of XCoord for each node.
    # ----------------------------------------
    loc = 0
    for i in range(NStr + 1):
        for j in range(NPZ):
            for k in range(NPX):
                mat[loc, 1] = k*XbayL
                loc = loc + 1
    # Determination of YCoord for each node.
    # ----------------------------------------
    loc = 0
    for i in range(NStr + 1):
        for j in range(PPF):
            mat[loc, 2] = i*StoryH
            loc = loc + 1
    # Determination of ZCoord for each node.
    # ----------------------------------------
    loc = 0
    for i in range(NStr + 1):
        for j in range(NPZ):
            for k in range(NPX):
                mat[loc, 3] = j*ZbayL
                loc = loc + 1
                
                
    # Node command application.
    # =========================
    for i in range(TNP):
        ops.node(int(mat[i, 0]), mat[i, 1], mat[i, 2], mat[i, 3])
        
        
    # If the foundation flexibility SSI consideration is applied, then, an ad-
    # ditional set of nodes are created with the same coordinates as for the nodes
    # at the base to be joined to the latter with zero-length elements.
    
    if flex in ('Y', 'YES', 'Yes', 'yES', 'yes', 'y'):
        # Creation of node numbering with 1e4 order for base designation.
        # ---------------------------------------------------------------
        bNssi = []
        for i in range(1, PPF + 1):
            bNssi.append(i)
        # Coordinates for these nodes.
        # ----------------------------
        # X Coordinates
        # -------------
        XCrdNssi = []
        for j in range(NPZ):
            for k in range(NPX):
                XCrdNssi.append(k*XbayL)
        # Z Coordinates
        # -------------
        ZCrdNssi = []
        for i in range(NPZ):
            for j in range(NPX):
                ZCrdNssi.append(i*ZbayL)
        # Opensees Command to Create the nodes.
        # --------------------------------------
        # matver = np.zeros((len(bNssi),4))
        for i in range(len(bNssi)):
            # matver[i,:] = int(bNssi[i]),XCrdNssi[i],0.0,ZCrdNssi[i]
            ops.node(int(bNssi[i]), XCrdNssi[i], 0.0, ZCrdNssi[i])
            

def MastNodeGen(NbaysX, NbaysZ, XbayL, ZbayL, StoryH, NStr, flex='No', coords=0):

    # Function to generate master nodes using geometric location for a regular
    # building with in-plan regularity. Also this function fixes the nodes in
    # such a way that it represents the only free DOFs of a rigid diaphragm.
    # MASTER NODE COORDINATES
    # =======================
    if type(coords) == int:  # Asumed coordinates as midpoint building plan sides.
        
        
        XCoordMN = NbaysX*XbayL/2.
        ZCoordMN = NbaysZ*ZbayL/2.
        for i in range(1, NStr + 1):
            ops.node(int(999 + 1e4 + 1e5*i), XCoordMN, i*StoryH, ZCoordMN)
        # MASTER NODES FIXITY
        # ===================
        fixity = [0, 1, 0, 1, 0, 1]
        for i in range(1, NStr + 1):
            ops.fix(int(999 + 1e4 + 1e5*i), *fixity)
    
        # If there are foundation flexibility considerations due to SSI effects,
        # we have to restrain the base nodes as well because the new fixed nodes
        # will be the base-duplicated ones.
        if flex in ('Y', 'YES', 'Yes', 'yES', 'yes', 'y'):
            ops.node(int(999+1e4), XCoordMN, 0.0, ZCoordMN)
            ops.fix(int(999+1e4), *fixity)
            
            
    else:     # Coordinates determined by ponderation of volumes Given Coordinates
        
        
        if flex in ('Y', 'YES', 'Yes', 'yES', 'yes', 'y'):
            for i in range(1, NStr + 1):
                ops.node(int(999 + 1e4 + 1e5*i), coords[i,0], i*StoryH, coords[i,1])
        else:
            for i in range(1, NStr + 1):
                ops.node(int(999 + 1e4 + 1e5*i), coords[i-1,0], i*StoryH, coords[i-1,1])
        # MASTER NODES FIXITY
        # ===================
        fixity = [0, 1, 0, 1, 0, 1]
        for i in range(1, NStr + 1):
            ops.fix(int(999 + 1e4 + 1e5*i), *fixity)
            
        # If there are foundation flexibility considerations due to SSI effects,
        # we have to restrain the base nodes as well because the new fixed nodes
        # will be the base-duplicated ones.
        if flex in ('Y', 'YES', 'Yes', 'yES', 'yes', 'y'):
            ops.node(int(999+1e4), coords[0,0], 0.0, coords[0,1])
            ops.fix(int(999+1e4), *fixity)


def MatGenRCB(fpc, Ec, fy, E0, b):
    
    
    import unitskN_m_C as units

    # Previous calculations and definitions
    # -------------------------------------
    ke = 0.85                                                             # [adim] assumed effective factor.
    rho_s = 0.0102                                                        # [m/m2] assumed confinement volumetric amount.
    fpco = 1.30*fpc                                                       # [kPa] unconfined concrete compressive strength.
    fpl = 0.50*ke*rho_s*fy                                                # [kPa] effective lateral confinement stress.
    fpcc = fpco*(-1.254 + 2.254*(1 + 7.94*fpl/fpco)**0.5 - 2* fpl/fpco)   # [kPa] confined strength of concrete.
    epsc0 = 2*fpco/Ec                                                     # [m/m] unconfined strain at peak stress value.
    epscc = epsc0*(1 + 5*(fpcc/fpco - 1))                                 # [m/m] confined strain at peak stress value.
    epssu = 0.12                                                          # [m/m] ultimate strain of reinforcement steel.
    epscu = 0.004 + 1.4*rho_s*fy*epssu/fpcc                               # [m/m] ultimate confined strain.
    fpcr = (0.63*(fpco/units.MPa)**0.5)*units.MPa                         # [kPa] rupture tensile modulus of concrete.
    epscr = fpcr/Ec                                                       # [m/m] concrete strain at rupture stress.
    alpha2 = 15                                                           # [adim] tensile concrete factor.
    Ets = abs(fpcr/(epscr*(alpha2 - 1)))                                  # [kPa] tension softening stiffness.
    epssp = 2*epsc0                                                       # [m/m] stress at spalling condition.
    
    # Compressive-strength properties.
    # ---------------------------------
    # Confinded concrete
    fpcConf = -1*fpcc               # [kPa] CONFINED concrete maximum stress.
    epsc0Conf = -1*epscc            # [m/m] concrete strain at maximum strength.
    fpcuConf = -0.10*fpcConf        # [kPa] concrete crushing strength.
    epsuConf = -epscu               # [m/m] concrete strain at crushin strength.
    lamb = 0.10                     # ratio between unloading slope at $epscu and initial slope.
    ftConf = fpcr                   # [kPa] tensile strength +tension
    EtsConf = Ets                   # tension softening stiffness
                   
    # Unconfined Concrete   
    fpcUnconf = -1*fpco               # [kPa] UNCONFINED concrete maximum stress.
    epsc0Unconf = -1*epsc0            # [m/m] concrete strain at maximum strength.
    fpcuUnconf = -0.10*fpcUnconf      # [kPa] concrete crushing strength.
    epsuUnconf = -1.5*epssp           # [m/m] concrete strain at crushin strength.
    ftUnconf = fpcr                   # [kPa] tensile strength +tension
    EtsUnconf = Ets                   # tension softening stiffness

    # CONCRETE UNIAXIAL MATERIAL GENERATION
    # =====================================
    # uniaxialMaterial('Concrete02', matTag, fpc, epsc0, fpcu, epsU, lambda, ft, Ets)
    ops.uniaxialMaterial('Concrete02', 1, fpcConf, epsc0Conf, fpcuConf, epsuConf, lamb, ftConf, EtsConf)              # concrete core
    ops.uniaxialMaterial('Concrete02', 2, fpcUnconf, epsc0Unconf, fpcuUnconf, epsuUnconf, lamb, ftUnconf, EtsUnconf)  # concrete cover

    # REINFORCING STEEL UNIAXIAL MATERIAL GENERATION
    # ===============================================
    R0 = 18
    cR1 = 0.925
    cR2 = 0.15
    # uniaxialMaterial('Steel02',matTag,Fy,E0,b,*params)
    ops.uniaxialMaterial('Steel02', 3, fy, E0, b, *[R0, cR1, cR2])


def FoundFlexMaterials(NbaysX,NbaysZ,XbayL,ZbayL,Ss,Vs30,gamma,nu,B,L,Re,D=0,omega=0,analtype='grav'):
    # If the foundation flexibility is considered, then the characteristics of
    # linear equivalent springs related to soil features shall be defined in this
    # function and incorporated at the model.
    # A one-spring stiffnes for both vertical and lateral flexibility are generated.
    # In the case of vertical flexibility, two uniaxial materials are generated
    # to account for rocking stiffness spring constants. These springs are meant
    # to be located at the perimeter of the base footprint within an effective length
    # meassured from the outer perimeter to the core of the building footprint area.

    import ASCE716SSI as SSI716
    import numpy as np
    # Previous calculations
    # ---------------------
    NAX = (NbaysX + 1)  # Number of Z axes along X direction.
    NAZ = (NbaysZ + 1)  # Number of X-axes along Z direction.
    PPF = NAX*NAZ

    # Unitary stiffness values for all the foundation.
    (kz, ky, kx, kzz, kyy, kxx) = SSI716.detKiisur(Ss, Vs30, gamma, nu, B, L, D, omega)  # [kN/m, kN-m/rad]
    # Regular vertical stiffness under every column.
    kzi = kz/PPF  # [kN/m]
    if analtype in ('grav'):
        # Regular lateral stiffness in Z global coordinates of the domain (y-direction of
        # literature axis convention).
        kyi = ky/PPF
        # Regular lateral stiffness in X global coordinates of the domain (X-direction of
        # literature axis convention).
        kxi = kx/PPF
    elif analtype in ('lat'):  # [kN/m]
        # Regular lateral stiffness in Z global coordinates of the domain (y-direction of
        # literature axis convention).
        kyi = ky*1e12/PPF
        # Regular lateral stiffness in X global coordinates of the domain (X-direction of
        # literature axis convention).
        kxi = kx*1e12/PPF

    # In accounting for rocking stiffness, stiffnesses for vertical springs are
    # calculated for the perimeter of foundation footprint equaling the moment
    # generated by these vertical springs and the rotational kxx or kyy.
    # The factors given at NIST GCR 12-917-21 are not correct.
    kz_xxtotRe = 12*kxx/(Re*B*(Re**2*B + 3*(1 - Re)*(Re*B + 2)))
    kz_yytotRe = 12*kyy/(Re*L*(Re**2*L + 3*(1 - Re)*(Re*L + 2)))
    kz_rockavgtotRe = (B *L**3*kz_yytotRe + B**3*L*kz_xxtotRe)/(B*L**3 + L*B**3)

    # The amount of coulumns within each ReB or ReL regions must be determined.
    # in order to define equivalent vertical stiffness for each column end.
    # -------------------------------------------------------------------------
    NAxesXReb = 1 + np.floor(Re*0.5*B/ZbayL)  # Number of axes in X direction within width Re*b; b = 0.5B
    NAxesZRel = 1 + np.floor(Re*0.5*L/XbayL)  # Number of axes in Z direction within width Re*l; l = 0.5L
    # Corner columns
    # --------------
    CornerCols = NAxesXReb*NAxesZRel  # number of columns at the corner of footprint plan included
    # within both Re*l and Re*b zones.
    NonCorner_yy = NAZ*NAxesZRel - 2*CornerCols  # Number of columns within Re*l not in corner zones.
    NonCorner_xx = NAX*NAxesXReb - 2*CornerCols  # Number of columns within Re*b not in corner zones.

    # Determining the unitary (per column base) spring stiffness constant.
    # --------------------------------------------------------------------
    kz_rockavg = kz_rockavgtotRe/(2*CornerCols)  # [kN/m] stiffness constant for corner zone columns.
    kz_yy = kz_yytotRe/NonCorner_yy  # [kN/m] stiffness constant for non corner zone columns, within Re*l.
    kz_xx = kz_xxtotRe/NonCorner_xx  # [kN/m] stiffness constant for non corner zone columns, within Re*b.

    # Creation of uniaxial materials with the stiffness created before.
    # =================================================================
    # The numbering of these materials starts from 101. In the case of corner
    # zones, this means, where columns bases corresponds to both areas within Re*L/2
    # or Re*B/2, a ponderate average is done to have a stiffness servin for
    # both directions.
    ops.uniaxialMaterial('Elastic', 101, kzi)            # Vertical regular spring for all columns.
    ops.uniaxialMaterial('Elastic', 102, kyi)            # Horizontal regular spring for all columns in Z-global direction.
    ops.uniaxialMaterial('Elastic', 103, kxi)            # Horizontal regular spring for all columns in X-global direction.
    ops.uniaxialMaterial('Elastic', 104, kz_xx)          # Spring for columns within Re*B/2 and no in corners accounting only for rocking.
    ops.uniaxialMaterial('Elastic', 105, kz_yy)          # Spring for columns within Re*L/2 and no in corners accounting only for rocking.
    ops.uniaxialMaterial('Elastic', 106, kz_rockavg)     # Spring for columns within corners accounting only for rocking.
    ops.uniaxialMaterial('Parallel', 107, *[101, 104])   # Vertical springs acounting for vertical support kzi and rocking around X-global axis.
    ops.uniaxialMaterial('Parallel', 108, *[101, 105])   # Vertical springs acounting for vertical support kzi and rocking around Z-global axis.
    ops.uniaxialMaterial('Parallel', 109, *[101, 106])   # Vertical corner zone springs acounting for vertical support kzi and rocking around both X- and Y-global axis.


def FoundFlexZLElements(NbaysX, NbaysZ, XbayL, ZbayL, B, L, Re):
    # This function generates Zero-Length elements connecting both nodes at the
    # same location set at the base of the structure.
    import numpy as np

    NAX = NbaysX + 1
    NAZ = NbaysZ + 1

    # The amount of coulumns within each ReB or ReL regions must be determined.
    # in order to define equivalent vertical stiffness for each column end.
    # -------------------------------------------------------------------------
    NAxesXReb = 1 + np.floor(Re*0.5*B/ZbayL)  # Number of axes in X direction within width Re*b; b = 0.5B
    NAxesZRel = 1 + np.floor(Re*0.5*L/XbayL)  # Number of axes in Z direction within width Re*l; l = 0.5L

    # Material lists
    # --------------
    matcorner = [103, 109, 102]
    matRel = [103, 108, 102]
    matReb = [103, 107, 102]
    matreg = [103, 101, 102]
    # Material directions
    # -------------------
    dirs = [1, 2, 3]

    loc = 1
    # element('zeroLength',eleTag,*eleNodes,'-mat',*matTags,'-dir',*dirs)
    for i in range(1, NAZ + 1):
        for j in range(1, NAX + 1):
            if i <= NAxesXReb or i > NAZ - NAxesXReb:  # Re*b zones for kz_xx equivalent stiffness.
                if j <= NAxesZRel or j > NAX - NAxesZRel:  # Defining corner zones.
                    ops.element('zeroLength',int(loc+4e4), *[int(loc),int(loc+1e4)], '-mat', *matcorner, '-dir', *dirs)
                    loc += 1
                else:
                    ops.element('zeroLength',int(loc+4e4), *[int(loc),int(loc+1e4)], '-mat', *matReb, '-dir', *dirs)
                    loc += 1
            elif i > NAxesXReb and i <= NAZ - NAxesXReb:
                if j <= NAxesZRel or j > NAX - NAxesZRel:
                    ops.element('zeroLength',int(loc+4e4), *[int(loc),int(loc+1e4)], '-mat', *matRel, '-dir', *dirs)
                    loc += 1
                elif j > NAxesZRel and j <= NAX - NAxesZRel:
                    ops.element('zeroLength',int(loc+4e4), *[int(loc),int(loc+1e4)], '-mat', *matreg, '-dir', *dirs)
                    loc += 1


def SPConstGen(NbaysX, NbaysZ, flex='No'):

    # To fix all nodes located at the base of the building. According to
    # numbering process executed on "NodeGen" function.

    # Numbers of base nodes.
    # ======================
    NPX = NbaysX + 1
    NPZ = NbaysZ + 1
    PPF = NPX*NPZ  # total number of nodes at base.
    constrVals = [1, 1, 1, 1, 1, 1]  # fully constrained nodes for fixed base.
    if flex in ('Y', 'YES', 'Yes', 'yES', 'yes', 'y'):
        for i in range(1, PPF + 1):
            ops.fix(int(i), *constrVals)
            ops.fix(int(i+1e4), *[0, 0, 0, 1, 1, 1])
    else:
        for i in range(1, PPF + 1):
            ops.fix(int(i+1e4), *constrVals)


def MPConstGen(NbaysX, NbaysZ, NStr, flex='No'):

    # To generate a Rigid Diaphragm MPCnostraint.
    NPX = NbaysX + 1  # number of nodes along X Global Axis.
    NPZ = NbaysZ + 1  # number of nodes along Z Global Axis.
    PPF = NPX*NPZ  # total number of points per flor.

    # Determination of integer tags for Nodes.
    # =======================================
    for j in range(1, NStr + 1):
        slaves = []                         # variable list taht will be filled with the floor nodes.
        master = int(999 + 1e4 + 1e5*j)     # Master node Tag for each floor.
        for i in range(1, PPF + 1):
            slaves.append(int(i + 1e4 + 1e5*j))    # slave nodes for beam-column joint nodes.

        ops.rigidDiaphragm(2, master, *slaves)

    # In case of accounting for flexibility at the base foundation.
    # ---------------------------------------------------------------
    # A rigid diaphragm is created for zero nodes as well. Nodes at the fixed-
    # base coordinates.
    if flex in ('Y', 'YES', 'Yes', 'yES', 'yes', 'y'):
        slavesBase = [int(number+1e4) for number in range(1,PPF+1)]
        ops.rigidDiaphragm(2, int(999+1e4), *slavesBase)

def GeomTransGen(ColTransfType, XBD=[0, 0], ZBD=[0, 0], ColD=[0, 0]):

    # The only input argument corresponds to whether the column shall be
    # use a geometric transformation other than 'Linear'.

    # BEAMS ALONG X GLOBAL AXIS.
    # ==========================
    vecxzBX = [0, 0, 1]
    dIBX = [0.5*ColD[1], 0, 0]
    dJBX = [-0.5*ColD[1], 0, 0]
    ops.geomTransf('Linear', 1, *vecxzBX, '-jntOffset', *dIBX, *dJBX)

    # BEAMS ALONG Z GLOBAL AXIS.
    # ==========================
    vecxzBZ = [-1, 0, 0]
    dIBZ = [0, 0, 0.5*ColD[0]]
    dJBZ = [0, 0, -0.5*ColD[0]]
    ops.geomTransf('Linear', 2, *vecxzBZ, '-jntOffset', *dIBZ, *dJBZ)

    # COLUMN ELEMENTS ALONG Z GLOBAL AXIS.
    # ==========================
    vecxzC = [0, 0, 1]
    dIC = [0, 0.5*max(XBD[1], ZBD[1]), 0]
    dJC = [0, -0.5*max(XBD[1], ZBD[1]), 0]
    ops.geomTransf(ColTransfType, 3, *vecxzC, '-jntOffset', *dIC, *dJC)


def LumpedMassGen(NbaysX,NbaysZ,XBDims,ZBDims,ColDims,gamma_conc,g,XbayL,ZbayL,NStr,StoryH,Qsd,flex='No'):
    # generates lumped masses to assign to each master node.
    # Includes SuperDead Loads converted to massed.

    import numpy as np
    NAX = NbaysX + 1
    NAZ = NbaysZ + 1
    PPF = (NbaysX + 1)*(NbaysZ + 1)
    Qsdmax = 0.50*Qsd*min([XbayL, ZbayL])  # distributed load per unit length

    # Mass of beam elements in X direction for each floor
    # ---------------------------------------------------
    loc = 0
    floormassXB = np.zeros((NStr, 1))
    for i in range(1, NStr + 1):
        fmassXB = []
        qSDmassXB = []
        for j in range(1, NAZ + 1):
            for k in range(1, NbaysX + 1):
                fmassXB.append(XBDims[loc, 0]*XBDims[loc, 1]*XbayL*gamma_conc/g)
                if j == 1 or j == NAZ:
                    if XbayL > ZbayL:
                        qdDistTRA = Qsdmax*(1 - 1/3*(min([XbayL, ZbayL])/max([XbayL, ZbayL]))**2)
                        qSDmassXB.append(qdDistTRA*XbayL/g)
                    else:
                        qdDistTRI = 2/3*Qsdmax
                        qSDmassXB.append(qdDistTRI*XbayL/g)
                else:
                    if XbayL > ZbayL:
                        qdDistTRA = Qsdmax*(1 - 1/3*(min([XbayL, ZbayL])/max([XbayL, ZbayL]))**2)
                        qSDmassXB.append(2*qdDistTRA*XbayL/g)
                    else:
                        qdDistTRI = 2/3*Qsdmax
                        qSDmassXB.append(2*qdDistTRI*XbayL/g)
        floormassXB[i - 1, 0] = sum(fmassXB) + sum(qSDmassXB)

    # Mass of beam elements in Z direction for each floor
    # ---------------------------------------------------
    loc = 0
    floormassZB = np.zeros((NStr, 1))
    for i in range(1, NStr + 1):
        fmassZB = []
        qSDmassZB = []
        for j in range(1, NAX + 1):
            for k in range(1, NbaysZ + 1):
                fmassZB.append(ZBDims[loc, 0]*ZBDims[loc, 1]*ZbayL*gamma_conc/g)
                if j == 1 or j == NAX:
                    if ZbayL > XbayL:
                        qdDistTRA = Qsdmax*(1 - 1/3*(min([XbayL, ZbayL])/max([XbayL, ZbayL]))**2)
                        qSDmassZB.append(qdDistTRA*ZbayL/g)
                    else:
                        qdDistTRI = 2/3*Qsdmax
                        qSDmassZB.append(qdDistTRI*ZbayL/g)
                else:
                    if ZbayL > XbayL:
                        qdDistTRA = Qsdmax*(1 - 1/3*(min([XbayL, ZbayL])/max([XbayL, ZbayL]))**2)
                        qSDmassZB.append(2*qdDistTRA*ZbayL/g)
                    else:
                        qdDistTRI = 2/3*Qsdmax
                        qSDmassZB.append(2*qdDistTRI*ZbayL/g)
        floormassZB[i - 1, 0] = sum(fmassZB) + sum(qSDmassZB)

    # Mass of column elements for each floor
    # ---------------------------------------
    loc = 0
    floormassCol = np.zeros((NStr, 1))
    for i in range(1, NStr + 1):
        fmassCol = []
        for j in range(1, PPF + 1):
            fmassCol.append(ColDims[loc, 0]*ColDims[loc, 1]*StoryH*gamma_conc/g)
        if i == NStr:
            floormassCol[i - 1, 0] = sum(fmassCol)/2
        else:
            floormassCol[i - 1, 0] = sum(fmassCol)

    # Total mass per floor and moment of inertia of floor
    # ---------------------------------------------------
    if flex in ('Y', 'YES', 'Yes', 'yES', 'yes', 'y'):
        # If SSI effects are considered, and flexibility of base is modeled
        # there will be another group of 3 DOFs corresponding to base elements
        # with a small mass and mass rotational inertia due to column elements
        # only. Hence, the vector of masses and weights will have NStr+1 rows.
        TMPF = np.zeros((NStr + 1, 1))
        MIF = np.zeros((NStr + 1, 1))
        TMPF[0, 0] = floormassCol[-1, 0]
        TWPF = [TMPF[0, 0] * g]
        MIF[0, 0] = TMPF[0, 0]/12*((NbaysX*XbayL)**2 + (NbaysZ*ZbayL)**2)
        for i in range(1, NStr + 1):
            TMPF[i, 0] = floormassXB[i - 1, 0] + floormassZB[i - 1, 0] + floormassCol[i - 1, 0]
            TWPF.append(TMPF[i, 0]*g)
            MIF[i, 0] = TMPF[i, 0]/12*((NbaysX*XbayL)**2 + (NbaysZ*ZbayL)**2)

        # Mass object definition for OpenSees
        # ------------------------------------
        matr = np.zeros((NStr + 1, 7))
        for i in range(NStr + 1):
            mx = TMPF[i]
            my = 1.0e-10
            mz = TMPF[i]
            ixx = 1.0e-10
            iyy = MIF[i]
            izz = 1.0e-10
            matr[i, 0], matr[i, 1], matr[i, 2], matr[i, 3], matr[i, 4], matr[i, 5], matr[i, 6], = \
                int(999 + 1e4 + 1e5*i), mx, my, mz, ixx, iyy, izz
            marg = []
            for j in range(6):
                marg.append(matr[i, j + 1])
            ops.mass(int(matr[i, 0]), *marg)

    else:
        TMPF = np.zeros((NStr, 1))
        MIF = np.zeros((NStr, 1))
        TWPF = []
        for i in range(1, NStr + 1):
            TMPF[i - 1, 0] = floormassXB[i - 1, 0] + floormassZB[i - 1, 0] + floormassCol[i - 1, 0]
            TWPF.append(TMPF[i - 1, 0] * g)
            MIF[i - 1, 0] = TMPF[i - 1, 0]/12*((NbaysX*XbayL)**2 + (NbaysZ*ZbayL)**2)

        # Mass object definition for OpenSees
        # ------------------------------------
        matr = np.zeros((NStr, 7))
        for i in range(1, NStr + 1):
            mx = TMPF[i - 1]
            my = 1.0e-10
            mz = TMPF[i - 1]
            ixx = 1.0e-10
            iyy = MIF[i - 1]
            izz = 1.0e-10
            matr[i - 1, 0], matr[i - 1, 1], matr[i - 1, 2], matr[i - 1, 3], matr[i - 1, 4], matr[i - 1, 5], matr[i - 1, 6], = \
                int(999 +1e4 + 1e5*i), mx, my, mz, ixx, iyy, izz
            marg = []
            for j in range(6):
                marg.append(matr[i - 1, j + 1])
            ops.mass(int(999+1e4+1e5*i), *marg)

    return TWPF, matr


def LiveLoadGen(NbaysX, NbaysZ, NStr, XbayL, ZbayL, Ql, Qlr):
    # Live load application to all elements in building. This function makes
    # diference between live load in common floors and at the roof.
    # LIVE LOADS OVER THE ELEMENTS
    # ============================
    import numpy as np
    NAX = NbaysX + 1
    NAZ = NbaysZ + 1
    Qlmax = 0.50*Ql*min([XbayL, ZbayL])  # distributed load per unit length
    Qlrmax = 0.50*Qlr*min([XbayL, ZbayL])  # distributed roof load per unit length

    LLtsTag = 3
    LLpatTag = 3

    # LOAD TIMESERIES GENERATION FOR LIVE LOAD
    # ========================================
    ops.timeSeries('Linear', LLtsTag)  # tsTag for Live Load will be always 3
    # Will always be Linear timeSeries in order to be applied variably
    # until a cFactor = 1.0 is reached. When the user wants to mantain the
    # load for a latter seismic or pushover analysis, the Openseespy
    # function loadConst() must be used.
    
    # LOAD PATTERN GENERATION FOR LIVE LOAD
    # =====================================
    ops.pattern('Plain', LLpatTag, LLtsTag)

    # Beams in X direction.
    # ---------------------
    loc = 0
    matrLLXB = np.zeros((NStr*NAZ*NbaysX, 1))   # including offset elements
    for i in range(1, NStr + 1):
        num = 1
        if i == NStr:
            for j in range(1, NAZ + 1):
                for k in range(1, NbaysX + 1):
                    if j == 1 or j == NAZ:
                        if XbayL > ZbayL:
                            qTRA = Qlrmax * (1 - 1/3*(min(XbayL, ZbayL)/max(XbayL, ZbayL))**2)
                            matrLLXB[loc, 0] = qTRA
                            ops.eleLoad('-ele', int(num+1e4+1e5*i),\
                                        '-type', '-beamUniform', -matrLLXB[loc, 0], 0.0)
                            loc += 1
                        else:
                            qTRI = Qlrmax*2/3
                            matrLLXB[loc, 0] = qTRI
                            ops.eleLoad('-ele', int(num+1e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLXB[loc, 0], 0.0)
                            loc += 1
                    else:
                        if XbayL > ZbayL:
                            qTRA = Qlrmax*(1 - 1/3*(min(XbayL, ZbayL)/max(XbayL, ZbayL))**2)
                            matrLLXB[loc, 0] = 2*qTRA
                            ops.eleLoad('-ele', int(num+1e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLXB[loc, 0], 0.0)
                            loc += 1
                        else:
                            qTRI = Qlrmax*2/3
                            matrLLXB[loc, 0] = 2*qTRI
                            ops.eleLoad('-ele', int(num+1e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLXB[loc, 0], 0.0)
                            loc += 1
                    num += 1
        else:
            for j in range(1, NAZ + 1):
                for k in range(1, NbaysX + 1):
                    if j == 1 or j == NAZ:
                        if XbayL > ZbayL:
                            qTRA = Qlmax*(1 - 1/3*(min(XbayL, ZbayL)/max(XbayL, ZbayL))**2)
                            matrLLXB[loc, 0] = qTRA
                            ops.eleLoad('-ele', int(num+1e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLXB[loc, 0], 0.0)
                            loc += 1
                        else:
                            qTRI = Qlmax*2/3
                            matrLLXB[loc, 0] = qTRI
                            ops.eleLoad('-ele', int(num+1e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLXB[loc, 0], 0.0)
                            loc += 1
                    else:
                        if XbayL > ZbayL:
                            qTRA = Qlmax*(1 - 1/3*(min(XbayL, ZbayL)/max(XbayL, ZbayL))**2)
                            matrLLXB[loc, 0] = 2*qTRA
                            ops.eleLoad('-ele', int(num+1e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLXB[loc, 0], 0.0)
                            loc += 1
                        else:
                            qTRI = Qlmax*2/3
                            matrLLXB[loc, 0] = 2*qTRI
                            ops.eleLoad('-ele', int(num+1e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLXB[loc, 0], 0.0)
                            loc += 1
                    num += 1
    
    # Beams in Z direction.
    # ---------------------
    loc = 0
    matrLLZB = np.zeros((NStr*NAX*NbaysZ, 1))
    for i in range(1, NStr + 1):
        num = 1
        if i == NStr:
            for j in range(1, NAX + 1):
                for k in range(1, NbaysZ + 1):
                    if j == 1 or j == NAX:
                        if ZbayL > XbayL:
                            qTRA = Qlrmax*(1 - 1/3*(min(XbayL, ZbayL)/max(XbayL, ZbayL))**2)
                            matrLLZB[loc, 0] = qTRA
                            ops.eleLoad('-ele', int(num+2e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLZB[loc, 0], 0.0)
                            loc += 1
                        else:
                            qTRI = Qlrmax*2/3
                            matrLLZB[loc, 0] = qTRI
                            ops.eleLoad('-ele', int(num+2e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLZB[loc, 0], 0.0)
                            loc += 1
                    else:
                        if ZbayL > XbayL:
                            qTRA = Qlrmax*(1 - 1/3*(min(XbayL, ZbayL)/max(XbayL, ZbayL))**2)
                            matrLLZB[loc, 0] = 2*qTRA
                            ops.eleLoad('-ele', int(num+2e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLZB[loc, 0], 0.0)
                            loc += 1
                        else:
                            qTRI = Qlrmax*2/3
                            matrLLZB[loc, 0] = 2*qTRI
                            ops.eleLoad('-ele', int(num+2e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLZB[loc, 0], 0.0)
                            loc += 1
                    num += 1
        else:
            for j in range(1, NAX + 1):
                for k in range(1, NbaysZ + 1):
                    if j == 1 or j == NAX:
                        if ZbayL > XbayL:
                            qTRA = Qlmax*(1 - 1/3*(min(XbayL, ZbayL)/max(XbayL, ZbayL))**2)
                            matrLLZB[loc, 0] = qTRA
                            ops.eleLoad('-ele', int(num+2e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLZB[loc, 0], 0.0)
                            loc += 1
                        else:
                            qTRI = Qlmax*2/3
                            matrLLZB[loc, 0] = qTRI
                            ops.eleLoad('-ele', int(num+2e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLZB[loc, 0], 0.0)
                            loc += 1
                    else:
                        if ZbayL > XbayL:
                            qTRA = Qlmax*(1 - 1/3*(min(XbayL, ZbayL)/max(XbayL, ZbayL))**2)
                            matrLLZB[loc, 0] = 2*qTRA
                            ops.eleLoad('-ele', int(num+2e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLZB[loc, 0], 0.0)
                            loc += 1
                        else:
                            qTRI = Qlmax*2/3
                            matrLLZB[loc, 0] = 2*qTRI
                            ops.eleLoad('-ele', int(num+2e4+1e5*i), \
                                        '-type', '-beamUniform', -matrLLZB[loc, 0], 0.0)
                            loc += 1
                    num += 1


def DeadLoadGen(NbaysX, NbaysZ, NStr, XBDims, ZBDims, ColDims, gamma_conc):
    # Generates only the loads due to self weight of elements. Taken from the
    # function definition of lumped masses generation.
    # DEAD LOADS FROM THE ELEMENTS
    # ============================
    import numpy as np
    NAX = NbaysX + 1
    NAZ = NbaysZ + 1

    DEADtsTag = 1
    DEADpatTag = 1

    # LOAD TIMESERIES GENERATION FOR DEAD LOAD FROM ELEMENTS
    # =======================================================
    ops.timeSeries('Linear', DEADtsTag)  # tsTag for DEAD Load will be always 1
    # Will always be Linear timeSeries in order to be applied variably
    # until a cFactor = 1.0 is reached. When the user wants to mantain the
    # load for a latter seismic or pushover analysis, the Openseespy
    # function loadConst() must be used.
    # LOAD PATTERN GENERATION FOR DEAD LOADS FROM ELEMENTS
    # =====================================================
    ops.pattern('Plain', DEADpatTag, DEADtsTag)

    # Beams in X direction.
    # ---------------------
    loc = 0
    matrDEADXB = np.zeros((NStr*NAZ*NbaysX, 1))
    for i in range(1, NStr + 1):
        num = 1
        for j in range(1, NAZ + 1):
            for k in range(1, NbaysX + 1):
                matrDEADXB[loc, 0] = XBDims[loc, 0]*XBDims[loc, 1]*gamma_conc
                ops.eleLoad('-ele', int(num+1e4+1e5*i), \
                            '-type', '-beamUniform', -matrDEADXB[loc, 0], 0.0)
                loc += 1
                num += 1

    # Beams in Z direction.
    # ---------------------
    loc = 0
    matrDEADZB = np.zeros((NStr*NAX*NbaysZ, 1))
    for i in range(1, NStr + 1):
        num = 1
        for j in range(1, NAX + 1):
            for k in range(1, NbaysZ + 1):
                matrDEADZB[loc, 0] = ZBDims[loc, 0]*ZBDims[loc, 1]*gamma_conc
                ops.eleLoad('-ele', int(num+2e4+1e5*i), \
                            '-type', '-beamUniform', -matrDEADZB[loc, 0], 0.0)
                loc += 1
                num += 1

    # Column elements.
    # ----------------
    loc = 0
    matrDEADCols = np.zeros((NStr*NAX*NAZ, 1))
    for i in range(1, NStr + 1):
        num = 1
        for j in range(1, NAX*NAZ + 1):
            matrDEADCols[loc, 0] = ColDims[loc, 0]*ColDims[loc, 1]*gamma_conc
            ops.eleLoad('-ele', int(num+3e4+1e5*i), \
                        '-type', '-beamUniform', 0.0, 0.0, -matrDEADCols[loc, 0])
            loc += 1
            num += 1


def SuperDeadLoadGen(NbaysX, NbaysZ, NStr, XbayL, ZbayL, Qsd):
    # Generates superDead patern distributed loads over elements.
    import numpy as np
    NAX = NbaysX + 1
    NAZ = NbaysZ + 1
    Qsdmax = 0.50*Qsd*min([XbayL, ZbayL])  # distributed load per unit length

    SDtsTag = 2
    SDpatTag = 2

    # LOAD TIMESERIES GENERATION FOR SUPER DEAD LOAD
    # ==============================================
    ops.timeSeries('Linear', SDtsTag)  # tsTag for SuperDead Load will be always 2
    # Will always be Linear timeSeries in order to be applied variably
    # until a cFactor = 1.0 is reached. When the user wants to mantain the
    # load for a latter seismic or pushover analysis, the Openseespy
    # function loadConst() must be used.
    
    # LOAD PATTERN GENERATION FOR SUPER DEAD LOAD
    # =====================================
    ops.pattern('Plain', SDpatTag, SDtsTag)

    # Beams in X direction.
    # ---------------------
    loc = 0
    matrSDXB = np.zeros((NStr*NAZ*NbaysX, 1))
    for i in range(1, NStr + 1):
        num = 1
        for j in range(1, NAZ + 1):
            for k in range(1, NbaysX + 1):
                if j == 1 or j == NAZ:
                    if XbayL > ZbayL:
                        qTRA = Qsdmax*(1 - 1/3*(min(XbayL, ZbayL)/max(XbayL, ZbayL))**2)
                        matrSDXB[loc, 0] = qTRA
                        ops.eleLoad('-ele', int(num+1e4+1e5*i), \
                                    '-type', '-beamUniform', -matrSDXB[loc, 0], 0.0)
                        loc += 1
                    else:
                        qTRI = Qsdmax*2/3
                        matrSDXB[loc, 0] = qTRI
                        ops.eleLoad('-ele', int(num+1e4+1e5*i), \
                                    '-type', '-beamUniform', -matrSDXB[loc, 0], 0.0)
                        loc += 1
                else:
                    if XbayL > ZbayL:
                        qTRA = Qsdmax*(1 - 1/3*(min(XbayL, ZbayL)/max(XbayL, ZbayL))**2)
                        matrSDXB[loc, 0] = 2*qTRA
                        ops.eleLoad('-ele', int(num+1e4+1e5*i), \
                                    '-type', '-beamUniform', -matrSDXB[loc, 0], 0.0)
                        loc += 1
                    else:
                        qTRI = Qsdmax*2/3
                        matrSDXB[loc, 0] = 2*qTRI
                        ops.eleLoad('-ele', int(num+1e4+1e5*i), \
                                    '-type', '-beamUniform', -matrSDXB[loc, 0], 0.0)
                        loc += 1
                num += 1

    # Beams in Z direction.
    # ---------------------
    loc = 0
    matrSDZB = np.zeros((NStr*NAX*NbaysZ, 1))
    for i in range(1, NStr + 1):
        num = 1
        for j in range(1, NAX + 1):
            for k in range(1, NbaysZ + 1):
                if j == 1 or j == NAX:
                    if ZbayL > XbayL:
                        qTRA = Qsdmax*(1 - 1/3*(min(XbayL, ZbayL)/max(XbayL, ZbayL))**2)
                        matrSDZB[loc, 0] = qTRA
                        ops.eleLoad('-ele', int(num+2e4+1e5*i), \
                                    '-type', '-beamUniform', -matrSDZB[loc, 0], 0.0)
                        loc += 1
                    else:
                        qTRI = Qsdmax*2/3
                        matrSDZB[loc, 0] = qTRI
                        ops.eleLoad('-ele', int(num+2e4+1e5*i), \
                                    '-type', '-beamUniform', -matrSDZB[loc, 0], 0.0)
                        loc += 1
                else:
                    if ZbayL > XbayL:
                        qTRA = Qsdmax*(1 - 1/3*(min(XbayL, ZbayL)/max(XbayL, ZbayL))**2)
                        matrSDZB[loc, 0] = 2*qTRA
                        ops.eleLoad('-ele', int(num+2e4+1e5*i), \
                                    '-type', '-beamUniform', -matrSDZB[loc, 0], 0.0)
                        loc += 1
                    else:
                        qTRI = Qsdmax*2/3
                        matrSDZB[loc, 0] = 2*qTRI
                        ops.eleLoad('-ele', int(num+2e4+1e5*i), \
                                    '-type', '-beamUniform', -matrSDZB[loc, 0], 0.0)
                        loc += 1
                num += 1


def ELFPForceGen(direction, Fx, flex='No'):
    # The applied loads into the list Fx must be in order from the first to floor
    # to the roof.

    import numpy as np

    ELFtsTag = 4
    ELFpatTag = 4

    # LOAD TIMESERIES GENERATION FOR QUIVALENT LATERAL FORCE
    # ======================================================
    ops.timeSeries('Linear', ELFtsTag)  # tsTag for Lateral Force Load will be always 4
    # Will always be Linear timeSeries in order to be applied variably
    # until a cFactor = 1.0 is reached. When the user wants to mantain the
    # load for a latter seismic or pushover analysis, the Openseespy
    # function loadConst() must be used.
    
    # LOAD PATTERN GENERATION FOR EQUIVALENT LATERAL FORCE
    # ====================================================
    ops.pattern('Plain', ELFpatTag, ELFtsTag)

    NLoads = len(Fx)
    MNLoads = np.zeros((NLoads, 7))
    if flex in ('Y', 'YES', 'Yes', 'yES', 'yes', 'y'):
        istory = 0
    else:
        istory = 1
    for i in range(NLoads):
        MNLoads[i, 0] = int(999+1e4+1e5*istory)
        istory += 1
        if direction in ('X', 'x'):
            MNLoads[i, 1] = Fx[i]
        elif direction in ('z', 'Z'):
            MNLoads[i, 3] = Fx[i]
        ops.load(int(MNLoads[i, 0]), *list(MNLoads[i, 1:]))


def POForceGen(NStr, Fx, direction, flex='No'):
    # This function executes a modal analysis in order to determine the load factors
    # used to determine the proportional magnitude of lateral forces according to
    # modal shape of the fundamental mode of vibration of the structure.

    import numpy as np

    # EIGEN ANALYSIS OF THE MODELED STRUCTURE
    # =========================================
    ops.eigen('-genBandArpack', 5)

    # Determination of fundamental mode shape.
    # ----------------------------------------
    Mod1Disp = ops.nodeEigenvector(int(999+1e4+1e5*NStr), 1, 1)
    Mod2Disp = ops.nodeEigenvector(int(999+1e4+1e5*NStr), 2, 1)
    if abs(Mod1Disp) >= abs(Mod2Disp):
        fundamentalDOF = 1
    else:
        fundamentalDOF = 3

    Vector = np.zeros((len(Fx), 1))

    if flex in ('Y', 'YES', 'Yes', 'yES', 'yes', 'y'):
        istory = 0
    else:
        istory = 1

    for i in range(len(Fx)):
        Vector[i, 0] = ops.nodeEigenvector(int(999+1e4+1e5*istory), 1, fundamentalDOF)
        istory += 1

    # Normalization of modal vector of fundamental period
    # --------------------------------------------------
    Vector = abs(Vector)
    maxVec = max(Vector)
    ModeShape = 1/maxVec*Vector

    POtsTag = 5
    POpatTag = 5

    # LOAD TIMESERIES GENERATION FOR PUSH-OVER LATERAL FORCE
    # ======================================================
    ops.timeSeries('Linear', POtsTag)  # tsTag for Push-over Force Load will be always 5
    # Will always be Linear timeSeries in order to be applied variably
    # until a cFactor = 1.0 is reached. When the user wants to mantain the
    # load for a latter seismic or pushover analysis, the Openseespy
    # function loadConst() must be used.
    # LOAD PATTERN GENERATION FOR PUSH-OVER LATERAL FORCE
    # ====================================================
    ops.pattern('Plain', POpatTag, POtsTag)

    # APPLICATION OF PROPORTIONAL FORCES TO MASTER NODES
    # ===================================================
    if flex in ('Y', 'YES', 'Yes', 'yES', 'yes', 'y'):
        istory = 0
    else:
        istory = 1
    MNLoads = np.zeros((len(Fx), 7))
    for i in range(len(Fx)):
        verif = []
        MNLoads[i, 0] = int(999+1e4+1e5*istory)
        istory += 1
        if direction in ('X', 'x'):
            MNLoads[i, 1] = ModeShape[i]*Fx[-1]
        elif direction in ('z', 'Z'):
            MNLoads[i, 3] = ModeShape[i]*Fx[-1]
        ops.load(int(MNLoads[i, 0]), *list(MNLoads[i, 1:]))
        verif.append(int(MNLoads[i, 0]))
        verif.append(list(MNLoads[i, 1:]))
        
    return ModeShape

def ElementGen(NbaysX,NbaysZ,XbayL,ZbayL,NStr,StoryH,XBDims,ZBDims,ColDims,BeamEffFact,ColEffFact,Ec,fy,
               EMs='elastic',XBreinf=0,ZBreinf=0,Colreinf=0,N=5,rec=0.0654,nuconc=0.2,dbar=0.025):
    # Defines the element along with the necesary combinations of integrator
    # adn sections depending on the type of element wanted. It could be elastic
    # or plastic with plastic concentrated hinges, distributed plasticity along
    # the hole length of element, or concentrated plasticity in element-end
    # regions.
    
    import openseespy.opensees as ops
    from OPSDefsMOD import FiberSectGen
    
    
    import numpy as np
    import unitskN_m_C as units
    
    # Previos Calculations
    # --------------------
    G = Ec/(2*(1 + nuconc))                  # [kPa] concrete shear modulus.
    NPX = NbaysX + 1                          # number of Z-axes along X-axis
    NPZ = NbaysZ + 1                          # number of X-axes along Z-axis.
    PPF = NPX*NPZ                             # total number of columns/joints 
                                              # per each floor.
    NbeamsX = NbaysX*NPZ                      # number of beams in X direction
                                              # per each floor.
    NbeamsZ = NbaysZ*NPX                      # number of beams in Z direction
                                              # per each floor.
    maxIter = 100
    tol = 1e-4
    
    # Elastic mechanical properties of element sections.
    # -------------------------------------------------
    # Beams in X-Direction
    xbelastic = np.zeros((NbeamsX*NStr,4))
    for i in range(NbeamsX*NStr):
        xbelastic[i,0] = XBDims[i,0]*XBDims[i,1]                          # [m2] section area.
        xbelastic[i,1] = BeamEffFact*XBDims[i, 0]*XBDims[i, 1]**3/12      # Iz = b*h^3/12
        xbelastic[i,2] = BeamEffFact*XBDims[i, 0]**3*XBDims[i, 1]/12      # Iy = b^3*h/12
        xbelastic[i,3] = xbelastic[i,1] + xbelastic[i,2]                  # J
    # Beams in Z-Direction
    zbelastic = np.zeros((NbeamsZ*NStr,4))
    for i in range(NbeamsZ*NStr):
        zbelastic[i,0] = ZBDims[i,0]*ZBDims[i,1]                          # [m2] section area.
        zbelastic[i,1] = BeamEffFact*ZBDims[i, 0]*ZBDims[i, 1]**3/12      # Iz = b*h^3/12
        zbelastic[i,2] = BeamEffFact*ZBDims[i, 0]**3*ZBDims[i, 1]/12      # Iy = b^3*h/12
        zbelastic[i,3] = zbelastic[i,1] + zbelastic[i,2]                  # J
    # Columns
    colelastic = np.zeros((PPF*NStr,4))
    for i in range(PPF*NStr):
        colelastic[i,0] = ColDims[i,0]*ColDims[i,1]                          # [m2] section area.
        colelastic[i,1] = ColEffFact*ColDims[i, 0]*ColDims[i, 1]**3/12       # Iz = b*h^3/12
        colelastic[i,2] = ColEffFact*ColDims[i, 0]**3*ColDims[i, 1]/12       # Iy = b^3*h/12
        colelastic[i,3] = colelastic[i,1] + colelastic[i,2]                  # J
    
    # Reinforcement cosniderations when using in analysis for design.
    # ===============================================================
    if type(XBreinf) == int:
        XBreinf = np.zeros((NbeamsX*NStr,4))
    if type(ZBreinf) == int:
        ZBreinf = np.zeros((NbeamsZ*NStr,4))
    if type(Colreinf) == int:
        Colreinf = np.zeros((PPF*NStr,10))
    
    # section('Elastic', secTag, E_mod, A, Iz, Iy, G_mod, Jxx)
    
    # beamIntegration('Lobatto', tag, secTag, N)
    # beamIntegration('HingeRadau', tag, secI, lpI, secJ, lpJ, secE)
    
    # element('forceBeamColumn', eleTag, *eleNodes, transfTag, integrationTag, 
    #                           '-iter', maxIter=10, tol=1e-12, '-mass', mass=0.0)
    
    
    # All elements are composed by the forceBeamColumn element opensees object
    # and are defined for 4 cases:
    # 1.- Elastic elements. 'elastic'
    # 2.- Plastic elements with distributed plasticity. 'plastic-1'
    # 3.- Plastic elements with concentrated  plasticity M-curv according to ASCE41-17.
    # 4.- Plastic elements with concentrated plasticity from fiber section.
    
    # IF ELEMENTS ARE MEANT TO BE ELASTIC.
    # ====================================
    if EMs in ('elastic'):
        
        # Beams in X-Direction
        # --------------------
        num = 0
        for i in range(1,NStr+1):
            for j in range(1,NbeamsX+1):
                
                # Section generation
                ops.section('Elastic',int(j+1e4+1e5*i),Ec,xbelastic[num,0],\
                            xbelastic[num,1],xbelastic[num,2],G,xbelastic[num,3])         # Central segment.
                        
                # Integrator Generation
                ops.beamIntegration('Lobatto', int(j+1e4+1e5*i), int(j+1e4+1e5*i), N)     # Central segment

                num += 1
        
        # Beams in Z-Direction
        # --------------------
        num = 0
        for i in range(1,NStr+1):
            for j in range(1,NbeamsZ+1):
                
                # Section generation
                ops.section('Elastic',int(j+2e4+1e5*i),Ec,zbelastic[num,0],\
                            zbelastic[num,1],zbelastic[num,2],G,zbelastic[num,3])          # Central segment.
                    
                # Integrator Generation
                ops.beamIntegration('Lobatto', int(j+2e4+1e5*i), int(j+2e4+1e5*i), N)      # Central segment.
                
                num += 1     
                    
        # Columns elements
        # ----------------
        num = 0
        for i in range(1,NStr+1):
            for j in range(1,PPF+1):
                
                # Section generation
                ops.section('Elastic',int(j+3e4+1e5*i),Ec,colelastic[num,0],\
                            colelastic[num,1],colelastic[num,2],G,colelastic[num,3])        # Central segment.
                    
                # Integrator Generation
                ops.beamIntegration('Lobatto', int(j+3e4+1e5*i), int(j+3e4+1e5*i), N)       # Central segment
                
                num += 1
                
                
    # IF ELEMENTS ARE MEANT TO BE MODELED WITH DISTRIBUTED PLASTICITY.
    # =============================================================
    elif EMs in ('plastic-1'):
        
        # Beams in X-Direction
        # --------------------
        num = 0
        for i in range(1,NStr+1):
            for j in range(1,NbeamsX+1):
                
                # Section generation
                # FiberSectGen(b,h,rec,As,Aprimas,As14,As24,As34,\
                    # NbarsBot,NbarsTop,Nbarsh14,Nbarsh24,Nbarsh34,secTag,GJ=1.0e10)
                FiberSectGen(XBDims[num,0],XBDims[num,1],rec,XBreinf[num,0],XBreinf[num,1],\
                                    0.,0.,0.,4,4,2,2,2,int(j+1e4+1e5*i),1e6)                     # Central segment.
            
                # Integrator Generation
                ops.beamIntegration('Lobatto', int(j+1e4+1e5*i), int(j+1e4+1e5*i), N)         # Central segment
                
                num += 1
        
        # Beams in Z-Direction
        # --------------------
        num = 0
        for i in range(1,NStr+1):
            for j in range(1,NbeamsZ+1):
                
                # Section generation
                # FiberSectGen(b,h,rec,As,Aprimas,As14,As24,As34,\
                    # NbarsBot,NbarsTop,Nbarsh14,Nbarsh24,Nbarsh34,secTag,GJ=1.0e10)
                FiberSectGen(ZBDims[num,0],ZBDims[num,1],rec,ZBreinf[num,0],ZBreinf[num,1],\
                                 0.,0.,0.,4,4,2,2,2,int(j+2e4+1e5*i),1e6)                        # central segment
                    
                # Integrator Generation
                ops.beamIntegration('Lobatto', int(j+2e4+1e5*i), int(j+2e4+1e5*i), N)         # Central segment
                
                num += 1     
                    
        # Columns elements
        # ----------------
        num = 0
        for i in range(1,NStr+1):
            for j in range(1,PPF+1):
                
                 # Section generation
                # FiberSectGen(b,h,rec,As,Aprimas,As14,As24,As34,\
                    # NbarsBot,NbarsTop,Nbarsh14,Nbarsh24,Nbarsh34,secTag,GJ=1.0e10)
                FiberSectGen(ColDims[num,0],ColDims[num,1],rec,Colreinf[num,0],Colreinf[num,1],\
                                 Colreinf[num,2],Colreinf[num,3],Colreinf[num,4],4,4,2,2,2,\
                                     int(j+3e4+1e5*i),1e6)                                          # Central Segment
                    
                # Integrator Generation
                ops.beamIntegration('Lobatto', int(j+3e4+1e5*i), int(j+3e4+1e5*i), N)               # Central segment
                    
                num += 1
                
                
    # IF ELEMENTS ARE MEANT TO BE PLASTIC WITH CONCENTRATED HINGE / MOMENT-CURV ASCE41-17.
    # ==================================================================================
    elif EMs in ('plastic-2'):
        
        # Beams in X-Direction
        # --------------------
        num = 0
        for i in range(1,NStr+1):
            for j in range(1,NbeamsX+1):
                
                # Section generation
                
                # Integrator Generation
                # ops.beamIntegration('HingeRadau', tag, secI, lpI, secJ, lpJ, secE)
                
                num += 1
        
        # Beams in Z-Direction
        # --------------------
        num = 0
        for i in range(1,NStr+1):
            for j in range(1,NbeamsZ+1):
                
                # Section generation
                
                # Integrator Generation
                # ops.beamIntegration('HingeRadau', tag, secI, lpI, secJ, lpJ, secE)
                
                num += 1     
                    
        # Columns elements
        # ----------------
        num = 0
        for i in range(1,NStr+1):
            for j in range(1,PPF+1):
                
                # Section generation
                
                # Integrator Generation
                # ops.beamIntegration('HingeRadau', tag, secI, lpI, secJ, lpJ, secE)
                
                num += 1
                
    # IF ELEMENTS ARE MEANT TO BE PLASTIC WITH DISTRIBUTED PLASTICITY AT ELEMENT-END REGIONS.
    # =======================================================================================
    elif EMs in ('plastic-3'):
        
        # Beams in X-Direction
        # --------------------
        num = 0
        for i in range(1,NStr+1):
            for j in range(1,NbeamsX+1):
                
                # Section generation
                ops.section('Elastic',int(j+1e4+1e5*i),Ec,xbelastic[num,0],\
                            xbelastic[num,1],xbelastic[num,2],G,xbelastic[num,3])                 # Central segment.
                # FiberSectGen(b,h,rec,As,Aprimas,As14,As24,As34,\
                    # NbarsBot,NbarsTop,Nbarsh14,Nbarsh24,Nbarsh34,secTag,GJ=1.0e10)
                FiberSectGen(XBDims[num,0],XBDims[num,1],rec,XBreinf[num,0],XBreinf[num,1],\
                                 0.,0.,0.,4,4,2,2,2,int(j+5e3+1e4+1e5*i),1e6)                     # Central segment i-end.
                FiberSectGen(XBDims[num,0],XBDims[num,1],rec,XBreinf[num,2],XBreinf[num,3],\
                                 0.,0.,0.,4,4,2,2,2,int(j+6e3+1e4+1e5*i),1e6)                       # Central segment j-end.
                        
                # Integrator Generation
                lp = 0.08*XbayL + 0.022*(dbar)*(fy/units.MPa)    # fy in [MPa]; [m] plastic length for beams in X direction.
                ops.beamIntegration('HingeRadau',int(j+1e4+1e5*i),int(j+5e3+1e4+1e5*i),lp,\
                                    int(j+6e3+1e4+1e5*i),lp,int(j+1e4+1e5*i))         # Central segment
                    
                num += 1
        
        # Beams in Z-Direction
        # --------------------
        num = 0
        for i in range(1,NStr+1):
            for j in range(1,NbeamsZ+1):
                
                # Section generation
                ops.section('Elastic',int(j+2e4+1e5*i),Ec,zbelastic[num,0],\
                            zbelastic[num,1],zbelastic[num,2],G,zbelastic[num,3])                   # Central segment.
                # FiberSectGen(b,h,rec,As,Aprimas,As14,As24,As34,\
                    # NbarsBot,NbarsTop,Nbarsh14,Nbarsh24,Nbarsh34,secTag,GJ=1.0e10)
                FiberSectGen(ZBDims[num,0],ZBDims[num,1],rec,ZBreinf[num,0],ZBreinf[num,1],\
                                 0.,0.,0.,4,4,2,2,2,int(j+5e3+2e4+1e5*i),1e6)                       # Central segment i-end.
                FiberSectGen(ZBDims[num,0],ZBDims[num,1],rec,ZBreinf[num,2],ZBreinf[num,3],\
                                 0.,0.,0.,4,4,2,2,2,int(j+6e3+2e4+1e5*i),1e6)                       # Central segment j-end.
                        
                # Integrator Generation
                lp = 0.08*ZbayL + 0.022*(dbar)*(fy/units.MPa)    # fy in [MPa]; [m] plastic length for beams in Z direction.
                ops.beamIntegration('HingeRadau',int(j+2e4+1e5*i),int(j+5e3+2e4+1e5*i),lp,\
                                    int(j+6e3+2e4+1e5*i),lp,int(j+2e4+1e5*i))                     # Central segment
                    
                num += 1     
                    
        # Columns elements
        # ----------------
        num = 0
        for i in range(1,NStr+1):
            for j in range(1,PPF+1):
                
                # Section generation
                ops.section('Elastic',int(j+3e4+1e5*i),Ec,colelastic[num,0],\
                            colelastic[num,1],colelastic[num,2],G,colelastic[num,3])                   # Central segment.
                # FiberSectGen(b,h,rec,As,Aprimas,As14,As24,As34,\
                    # NbarsBot,NbarsTop,Nbarsh14,Nbarsh24,Nbarsh34,secTag,GJ=1.0e10)
                FiberSectGen(ColDims[num,0],ColDims[num,1],rec,Colreinf[num,0],Colreinf[num,1],\
                                 Colreinf[num,2],Colreinf[num,3],Colreinf[num,4],4,4,2,2,2,\
                                     int(j+5e3+3e4+1e5*i),1e6)                                         # Central segment i-end.
                FiberSectGen(ColDims[num,0],ColDims[num,1],rec,Colreinf[num,5],Colreinf[num,6],\
                                 Colreinf[num,7],Colreinf[num,8],Colreinf[num,9],4,4,2,2,2,\
                                     int(j+6e3+3e4+1e5*i),1e6)                                         # Central segment j-end.
                    
                # Integrator Generation
                lp = 0.08*StoryH + 0.022*(dbar)*(fy/units.MPa)  # fy in [MPa]; [m] plastic length for columns.
                ops.beamIntegration('HingeRadau',int(j+3e4+1e5*i), int(j+5e3+3e4+1e5*i), lp, \
                                    int(j+6e3+3e4+1e5*i), lp, int(j+3e4+1e5*i))                   # Central segment
                    
                num += 1
                
    # =============================================================================================================================
    # =============================================================================================================================
    # =============================================================================================================================
    # GENERATION OF ELEMENT OBJECTS
    # =============================
    for i in range(1,NStr+1):
        # For beams in X direction
        # ------------------------
        for j in range(1,NbeamsX+1):
            # Element Generator
            if j%NbaysX == 0:
                iNode = NPX*int(j/NbaysX) - 1
                jNode = NPX*int(j/NbaysX)
            else:
                iNode = NPX*int(j/NbaysX) + j%NbaysX
                jNode = NPX*int(j/NbaysX) + j%NbaysX + 1
                
            ops.element('forceBeamColumn', int(j+1e4+1e5*i), *[int(iNode+1e4+1e5*i),(jNode+1e4+1e5*i)], \
                        1, int(j+1e4+1e5*i),'-iter', maxIter, tol)                                # central segment element
                
                
        # For Beams in Z Direction
        # ------------------------
        for j in range(1,NbeamsZ+1):
            # Element Generator
            
            if j%NbaysZ == 0:
                iNode = int(np.ceil(j/NbaysZ)) + ((j-1)%NbaysZ)*NPX
                jNode = int(np.ceil(j/NbaysZ)) + ((j-1)%NbaysZ)*NPX + NPX
            else:
                iNode = int(np.ceil(j/NbaysZ)) + (j%NbaysZ-1)*NPX
                jNode = int(np.ceil(j/NbaysZ)) + (j%NbaysZ-1)*NPX + NPX
                
            ops.element('forceBeamColumn', int(j+2e4+1e5*i), *[int(iNode+1e4+1e5*i),(jNode+1e4+1e5*i)], \
                        2, int(j+2e4+1e5*i),'-iter', maxIter, tol)                                # central segment element
           
            
        # For Columns
        # ------------
        for j in range(1,PPF+1):
            # Element Generator
            # -------------------------------------------------------------------------------------------------
            ops.element('forceBeamColumn', int(j+3e4+1e5*i), *[int(j+1e4+1e5*(i-1)),int(j+1e4+1e5*i)], \
                        3, int(j+3e4+1e5*i),'-iter', maxIter, tol)                                # central segment element