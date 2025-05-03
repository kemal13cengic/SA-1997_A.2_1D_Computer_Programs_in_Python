##############################################

# WORK IN PROGRESS

# TESTED AND WORKING FOR LINEAR PROBLEMS (CONSTANT CONDUCTIVITY)


##############################################


#-------------------------------------------------------------
#   SOLVES STEADY 1D HEAT CONDUCTION PROBLEM FOR:
#       - DIRICHLET BOUNDARY CONDITIONS
#       - CONSTANT OR VARIABLE CONDUCTIVITY
#   BY EMPLOYING:
#       - TDMA OR GSM LINEAR EQUATION SOLVER
#-------------------------------------------------------------

import numpy as np
import math

np.set_printoptions(precision=30, threshold=np.inf)

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
ITGS=0                  #
GSMTOL=0                # GSM TOLERANCE

CON=np.zeros(NCV+2)     # CONDUCTIVITY
TEX=np.zeros(NCV)       #
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

T[0]=10                 # FIRST WALL BOUNDARY VALUE
T[N+1]=210              # LAST WALL BOUNDARY VALUE

#   CHOOSE SOLVER (1 - TDAMA, 2 - GSM) AND GSM TOLERANCE
#-------------------------------------------------------------

IS=1                   # CHOOSE SOLVER
GSMTOL=0.000001         # GSM TOLERANCE

#   MAX NUMBER OF ITERATIONS (1 FOR LINEAR PROBLEM)
#-------------------------------------------------------------

MAXIT=1

#   INITIALISE THE CONDUCTIVITY
#-------------------------------------------------------------
CON0=1
CON[:N+2]=CON0

#-------------------------------------------------------------
#   ITERATION LOOP (START)
#-------------------------------------------------------------

NIT = 0
while NIT<MAXIT:
    NIT=NIT+1
    #-------------------------------------------------------------
    #   ASSEMBLE COEFFICIENT MATRIX AND THE RIGHT HAND SIDE VECTOR
    #-------------------------------------------------------------

    for i in range(N):
        AE[i] = 0.5 * (CON[i+1] + CON[i+2]) * S / DX
        AW[i] = 0.5 * (CON[i+1] + CON[i]) * S / DX # PYTHON INDEX = FORTRAN INDEX - 1


    AW[0]=2*AW[0]
    AE[N-1]=2*AE[N-1]

    for i in range(N):
        AP[i]=AW[i]+AE[i]
        B[i]=0


    #   BOUNDARY CONDITIONS
    #-------------------------------------------------------------

    B[0] = B[0] + AW[0] * T[0]
    B[N-1] = B[N-1] + AE[N-1] * T[N+1]

    #-------------------------------------------------------------
    #   SOLVE LINEAR EQUATION SYSTEM
    #-------------------------------------------------------------

    if(IS==1):
        TDMA()
    elif(IS==2):
        GSM()

    #-------------------------------------------------------------
    #   UPDATE CONDUCTIVITY (FOR NONLINEAR PROBLEM)
    #-------------------------------------------------------------

    if(MAXIT>1):
        for i in range(1,N+1):
            CON[i]=CON0*T[i]

    #-------------------------------------------------------------
    #   ITERATION LOOP (END)
    #-------------------------------------------------------------

#-------------------------------------------------------------
#   PRINT RESULTS ON THE SCREEN AND ON THE OUTPUT FILE
#-------------------------------------------------------------
else:
    if(IS==1):
        print("\n TDMA SOLVER \n")
    if(IS==2):
        print(f"\n GAUSS-SEIDEL SOLVER: {ITGS} ITERATIONS \n")

print(f"{'I':<3} {'X':<5} {'T':<20} {'TEXACT':<20} {'ERROR':<20}")
ERROR =0
for i in range(1,N+1):
    if MAXIT==1:
        #   CASE 1 : CONSTANT CONDUCTIVITY
        TEX[i-1]=((T[i+1]-T[0])*(X[i]-X[0]))/(X[i+1]-X[0])+T[0] # FLOAT POINT ERROR PRESENT
    else:
        #   CASE 2 : CONDUCTIVITY PROPORTIONAL TO TEMPERATURE
        TEX[i-1]=math.sqrt(((T[i+1]**2-T[0]**2)*(X[i]-X[0]))/(X[i+1]-X[0])+T[0]**2)
    # TEX=[20,40,60,80,100,120,140,160,180,200]                   # UNCOMMENT FOR ORIGINAL INPUT (NO FLOAT POINT ERROR)
    ERROR=ERROR+abs(TEX[i-1]-T[i])
    print(f"{i:<3} {X[i]:<5.1f} {T[i]:<20.15f} {TEX[i-1]:<20.15f} {(T[i] - TEX[i-1]):<20.15e}")

    print(type(T[1]))








