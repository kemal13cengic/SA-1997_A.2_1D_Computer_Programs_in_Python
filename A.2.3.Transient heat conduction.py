#-------------------------------------------------------------
#   SOLVES UNSTEADY 1D HEAT CONDUCTION PROBLEM FOR:
#       - DIRICHLET BOUNDARY CONDITIONS
#       - CONSTANT CONDUCTIVITY
#   BY EMPLOYING:
#       - EXPLICIT, FULLY IMPLICIT AND CRANK-NICOLSON SCHEME
#       - TDMA LINEAR EQUATION SOLVER
#-------------------------------------------------------------

import numpy as np

NCV=10
AE=np.zeros(NCV)
AW=np.zeros(NCV)
AP=np.zeros(NCV)
B=np.zeros(NCV)
T=np.zeros(NCV+2)
T0=np.zeros(NCV+2)
X=np.zeros(NCV+2)

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

#-------------------------------------------------------------
#   INPUT DATA
#-------------------------------------------------------------


#   DEFINE THE GRID
#-------------------------------------------------------------
N=10                    # NUMBER OF CELLS
X[0]=0                  # FIRST WALL POSITION
X[N+1]=10               # LAST WALL POSITION
DX=(X[N+1]-X[0])/N
S=1                     # WALL AREA
V=S*DX                  # VOLUME OF CELL
X[1]=X[0]+DX/2

for i in range(2,N+1):
    X[i]=X[i-1]+DX

#   INITIAL AND BOUNDARY VALUES
#-------------------------------------------------------------

T[0]=100                 # FIRST WALL BOUNDARY VALUE
T[N+1]=210              # LAST WALL BOUNDARY VALUE
T[1:-1]=10              # INITIAL VALUE
T0[0]=T[0]
T0[N+1]=T[N+1]

CON=1                   # CONDUCTIVITY
DEN=1                   # DENSITY
SPH=1                   # SPECIFIC HEAT

DT=0.4                  # TIME STEP SIZE
NTSTEP=200              # NUMBER OF TIME STEPS
NPRF=10                 # FREQUENCY OF PRINT OUT

IT=3                    # TEMPORAL DIFFERENCING SCHEME (1-EX, 2-FI, 3-CN)

#-------------------------------------------------------------
#   TIME LOOP (START)
#-------------------------------------------------------------

TIME=0
NT=0
for i in range(1,NTSTEP+1):
    NT+=1
    TIME += DT
    T0[1:N+1] = T[1:N+1]


    #-------------------------------------------------------------
    #   ASSEMBLE COEFFICIENT MATRIX AND THE RIGHT HAND SIDE VECTOR
    #-------------------------------------------------------------

    #   DIFFUSION COEFFICIENTS
    #-------------------------------------------------------------

    for i in range(0,N):
        AE[i] = CON * S / DX
        AW[i] = CON * S / DX

    AW[0]=2*AW[0]
    AE[N-1]=2*AE[N-1]

    #   TRANSIENT TERM DISCRETISATION AND AP (AND T FOR EX SCHEME) CALCULATION
    #-------------------------------------------------------------

    if(IT==1):
        # EXPLICIT
        for i in range(0,N):
            B[i] = AE[i] * T0[i+2] + AW[i] * T0[i] + ((V/DT) * DEN * SPH - AE[i] - AW[i]) * T0[i+1]
            AP[i] = (V/DT) * DEN * SPH
            T[i+1] = B[i] / AP[i]
    elif(IT==2):
        # FULLY IMPLICIT
        for i in range(0,N):
            B[i] = (V/DT) * DEN * SPH * T0[i+1]
            AP[i] = (V/DT) * DEN * SPH + AE[i] + AW[i]
    elif(IT==3):
        # CRANK-NICOLSON
        for i in range(0,N):
            AE[i] = 0.5 * AE[i]
            AW[i] = 0.5 * AW[i]
            B[i] = AE[i] * T0[i+2] + AW[i] * T0[i] + ((V/DT) * DEN * SPH - AE[i] - AW[i]) * T0[i+1]
            AP[i] = (V/DT) * DEN * SPH + AE[i] + AW[i]

    #   BOUNDARY CONDITIONS
    #-------------------------------------------------------------

    B[0] = B[0] + AW[0] * T[0]
    B[N-1] = B[N-1] + AE[N-1] * T[N+1]

#-------------------------------------------------------------
#   SOLVE LINEAR EQATION SYSTEM (FOR F-I AND C-N SCHEME)
#-------------------------------------------------------------
    if IT!=1:
        TDMA()

#-------------------------------------------------------------
#   PRINT RESULTS ON THE SCREEN AND ON THE OUTPUT FILE
#-------------------------------------------------------------
    if NT==1:
        CO= CON * DT / (DEN * SPH * DX**2)
        print("\n COURANT NUMBER : CO = ",CO)
        print("\n TIME STEP : DT = ",DT)
        if   IT == 1:
           print("\n EXPLICIT TRANSIENT TERM DISCRETIZATION")
        elif IT == 2:
           print("\n IMPLCIT TRANSIENT TERM DISCRETIZATION")
        elif IT == 3:
           print("\n CRANK-NICOLSON TRANSIENT TERM DISCRETIZATION")

    if NT % NPRF == 0:
        print("\nTIMESTEP : ", NT, "TIME : T = ", TIME)
        print(f"{'I':<3} {'X':<5} {'T':<20}")
        for i in range(1, N + 1):
            print(f"{i:<3} {X[i]:<5} {T[i]:<20}")


