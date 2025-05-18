#-------------------------------------------------------------
#   SOLVES STEADY 1D HEAT CONDUCTION-CONVECTIN PROBLEM FOR:
#       - DIRICHLET BOUNDARY CONDITIONS
#       - CONSTANT CONDUCTIVITY, DENSITY AND SPECIFIC HEAT
#         AND A GIVEN VELOCITY
#   BY EMPLOYING:
#       - 1ST ORDER, 2ND ORDER OR A 1ST_2ND ORDER BLEND SCHEME
#         FOR CONVECTION APPROXIMATION
#       - TDMA OR GSM LINEAR EQUATION SOLVER
#-------------------------------------------------------------

import numpy as np
import math

#-------------------------------------------------------------
#   SUBROUTINE TDMA
#-------------------------------------------------------------
def TDMA():
    global T, AE, AW, B, AP, N

    #   CALCULATE NEW AP(i) AND B(i)
    # -------------------------------------------------------------
    for i in range(1,N):
        AP[i] = AP[i] - AW[i] * AE[i-1] / AP[i-1]
        B[i] = B[i] + AW[i] * B[i-1] / AP[i-1]

    #   CALCULATE VARIABLE VALUES - BACKWARD SUBSTITUTION
    # -------------------------------------------------------------
    T[N] = B[N-1] / AP[N-1]
    for i in range(N-1,0,-1):
        T[i] = (B[i-1] + AE[i-1] * T[i+1]) / AP[i-1]
    print(T)

#-------------------------------------------------------------
#   SUBROUTINE GSM
#-------------------------------------------------------------
def GSM():
    global T, AE, AW, B, AP, GSMTOL, N, ITGS

    # INITIALIZE FIELDS AND SET "BOUNDRY COEFFICIENTS" TO ZERO
    # -------------------------------------------------------------
    for i in range(1,N+1):
        T[i]=0

    AW[0]=0
    AE[N-1]=0

    #   CALCULATE RESIDUAL AND CHECK CONVERGENCE
    # -------------------------------------------------------------
    for ITGS in range(0,1000):
        RES=0
        for i in range(1,N+1):
            RES=RES+abs(AE[i - 1] * T[i + 1] + AW[i - 1] * T[i - 1] + B[i - 1] - AP[i - 1] * T[i])
        if ITGS==0:
            RESN=RES    # NORMALIZE THE RESIDUAL
            RES=RES/RESN
        if RES<GSMTOL:
            break       # CONVERGENCE ACHIEVED

        #   CALCULATE NEW SOLUTION
        # -------------------------------------------------------------
        for i in range(1, N + 1):
            T[i] = (AE[i - 1] * T[i + 1] + AW[i - 1] * T[i - 1] + B[i - 1]) / AP[i - 1]



#-------------------------------------------------------------
#   INITIALIZATION
#-------------------------------------------------------------

NCV=10
N=0
AE=np.zeros(NCV)
AW=np.zeros(NCV)
AP=np.zeros(NCV)
B=np.zeros(NCV)         # RIGHT HAND SIDE VECTOR
T=np.zeros(NCV+2)       # TEMPERATURE
ITGS=0                  # ITERATIONS OF GAUSS-SEIDEL
GSMTOL=0                # GSM TOLERANCE

TEX=np.zeros(NCV)       # EXACT TEMPERATURE
X=np.zeros(NCV+2)

#-------------------------------------------------------------
#   INPUT DATA
#-------------------------------------------------------------


#   DEFINE THE GRID
#-------------------------------------------------------------
N=10                    # NUMBER OF CELLS (INPUT NCV AS SAME)
X[0]=0                  # FIRST WALL POSITION
X[N+1]=10               # LAST WALL POSITION
DX=(X[N+1]-X[0])/N
S=1                     # WALL AREA
X[1]=X[0]+DX/2

for i in range(2,N+1):
    X[i]=X[i-1]+DX

#   BOUNDARY VALUES
#-------------------------------------------------------------

T[0]=10                # FIRST WALL BOUNDARY VALUE
T[N+1]=210              # LAST WALL BOUNDARY VALUE

#   READ CONDUCTIVITY, DENSITY AND SPECIFIC HEAT
#-------------------------------------------------------------

CON=1
DEN=1
SPH=1

#   READ VELOCITY
#-------------------------------------------------------------

VEL=0.1

#   CHOOSE CONVECTION SCHEME (1 - CD, 2 - UD, 3-CD_UD BLEND)
#   AND BLENDING FACTOR (1 - CD, 0 - UD)
#-------------------------------------------------------------

IC=2
GAMMA=1

#   CHOOSE SOLVER (1 - TDMA, 2 - GSM) AND GSM TOLERANCE
#-------------------------------------------------------------

IS=2                   # CHOOSE SOLVER
GSMTOL=0.000001         # GSM TOLERANCE

#   MAX NUMBER OF ITERATIONS (1 FOR LINEAR PROBLEM)
#-------------------------------------------------------------

MAXIT=10

#-------------------------------------------------------------
#   ITERATION LOOP (START)
#-------------------------------------------------------------

NIT = 0
while NIT<MAXIT:
    NIT=NIT+1
    #-------------------------------------------------------------
    #   ASSEMBLE COEFFICIENT MATRIX AND THE RIGHT HAND SIDE VECTOR
    #-------------------------------------------------------------

    #   DIFFUSION COEFFICIENTS
    #-------------------------------------------------------------

    for i in range(N):
        AE[i]=CON*S/DX
        AW[i]=CON*S/DX

    #   CONVECTION COEFFICIENTS
    # -------------------------------------------------------------

    FLUX = DEN * VEL * S

    if(IC==1):
        for i in range(N):
            AE[i] = AE[i] - 0.5 * SPH * FLUX
            AW[i] = AW[i] + 0.5 * SPH * FLUX

        AW[0] = 2 * AW[0]
        AE[N-1] = 2 * AE[N-1]
    else:
        AW[0] = 2 * AW[0]
        AE[N - 1] = 2 * AE[N - 1]

        for i in range(N):
            AE[i] = AE[i] - SPH * min(FLUX, 0)
            AW[i] = AW[i] - SPH * min(-FLUX, 0)

    #   CENTRAL COEFFICIENT AND RIGHT HAND SIDE VECTOR
    # -------------------------------------------------------------
    for i in range(N):
        AP[i]=AW[i]+AE[i]
        B[i]=0
    if(IC==3):
        for i in range(N):
            B[i]=B[i] + GAMMA * SPH * ((min(FLUX, 0) - 0.5 * FLUX) * (T[i+2] - T[i+1]) + (min(-FLUX,0) + 0.5 * FLUX)*(T[i]-T[i+1]))

    #   BOUNDARY CONDITIONS
    #-------------------------------------------------------------

    B[0] = B[0] + AW[0] * T[0]
    B[N-1] = B[N-1] + AE[N-1] * T[N+1]
    if(IC==3):
        B[0] = B[0] + GAMMA * SPH * 0.5 * FLUX * (T[0]-T[1])
        B[N-1] = B[N-1] - GAMMA * SPH * 0.5 * FLUX * (T[N+1] - T[N]) #PROBLEM OVDJE
    #-------------------------------------------------------------
    #   SOLVE LINEAR EQUATION SYSTEM
    #-------------------------------------------------------------

    if(IS==1):
        TDMA()
    elif(IS==2):
        GSM()

    #-------------------------------------------------------------
    #   ITERATION LOOP (END)
    #-------------------------------------------------------------

#-------------------------------------------------------------
#   PRINT RESULTS ON THE SCREEN AND ON THE OUTPUT FILE
#-------------------------------------------------------------
else:
    PE= VEL * (X[N+1] - X[0]) * SPH * DEN / CON
    print("\nPECLET NUMBER: PE = ",PE,"\n")
    if (IC==1):
        print("CDS USED FOR CONVECTION ")
    elif (IC==2):
        print("UDS USED FOR CONVECTION ")
    elif (IC==3):
        print("BLEND OF CDS AND UDS USED FOR CONVECTION: GAMMA = ",GAMMA)

    if(IS==1):
        print("\nTDMA SOLVER \n")
    elif(IS==2):
        print(f"\nGAUSS-SEIDEL SOLVER: {ITGS} ITERATIONS \n")

    print(f"{'I':<3} {'X':<5} {'T':<20} {'TEXACT':<20} {'ERROR':<20}")

    ERROR = 0
    for i in range(1,N+1):
        TEX[i-1]=T[0] + (math.exp(PE*(X[i]-X[0])/(X[N+1]-X[0]))-1)/(math.exp(PE)-1)*(T[N+1]-T[0])
        ERROR+=abs(TEX[i-1]-T[i])
        print(f"{i:<3} {X[i]:<5} {T[i]:<20} {TEX[i - 1]:<20} {TEX[i-1]-T[i]:<20}")

    ERROR=ERROR/N
    print("\nAVERAGE ERROR = ", ERROR)