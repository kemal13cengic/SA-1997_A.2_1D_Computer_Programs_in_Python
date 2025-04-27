
##############################################

# WOOOOOOOOOOOOOOOOOORK IN PROOOOOOOOOOOGRESS

# PROVJERI CINJENICU DA u FORTRAN_ARRAY(i) = np.array(i-1)
# SVE MORAS NESTATI U TIME LOOP!!!

##############################################




#-------------------------------------------------------------
#   SOLVES UNSTEADY 1D HEAT CONDUCTION PROBLEM FOR:
#       - DIRICHLET BOUNDARY CONDITIONS
#       - CONSTANT CONDUCTIVITY
#   BY EMPLOYING:
#       - EXPLICIT, FULLY IMPLICIT AND CRANK-NICOLSON SCHEME
#       - TDMA LINEAR EQUATION SOLVER
#-------------------------------------------------------------

import numpy as np

NCV=100
AE=np.zeros(NCV)
AW=np.zeros(NCV)
AP=np.zeros(NCV)
B=np.zeros(NCV)
T=np.zeros(NCV+2)
T0=np.zeros(NCV+2)
X=np.zeros(NCV+2)

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

T[0]=10                 # FIRST WALL BOUNDARY VALUE
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

IT=1                    # TEMPORAL DIFFERENCING SCHEME (1-EX, 2-FI, 3-CN)

#-------------------------------------------------------------
#   TIME LOOP (START)
#-------------------------------------------------------------

TIME=0
for i in range(1,NTSTEP+1):
    TIME += DT
    T0[1:N+1] = T[1:N+1]


#-------------------------------------------------------------
#   ASSEMBLE COEFFICIENT MATRIX AND THE RIGHT HAND SIDE VECTOR
#-------------------------------------------------------------


#   DIFFUSION COEFFICIENTS
#-------------------------------------------------------------

for i in range(1,N):
    AE[i] = CON * S / DX
    AW[i] = CON * S / DX

AW[1]=2*AW[1]
AE[N]=2*AE[N]

#   TRANSIENT TERM DISCRETISATION AND AP (AND T FOR EX SCHEME) CALCULATION
#-------------------------------------------------------------

if(IT==1):
    # EXPLICIT
    for i in range(1,N):
        B[i] = AE[i] * T0[i+1] + AW[i] * T0[i-1] + ((V/DT) * DEN * SPH - AE[i] - AW[i]) * T0[i]
        AP[i] = (V/DT) * DEN * SPH + AE[i] + AW[i]
elif(IT==2):
    # FULLY IMPLICIT
    for i in range(1,N):
        B[i] = (V/DT) * DEN * SPH * T0[i]
        AP[i] = (V/DT) * DEN * SPH + AE[i] + AW[i]
elif(IT==3):
    # CRANK-NICOLSON
    for i in range(1,N):
        AE[i] = 0.5 * AE[i]
        AW[i] = 0.5 * AW[i]
        B[i] = AE[i] * T0[i+1] + AW[i] * T0[i-1] + ((V/DT) * DEN * SPH - AE[i] - AW[i]) * T0[i]
        AP[i] = (V/DT) * DEN * SPH + AE[i] + AW[i]

#   BOUNDARY CONDITIONS
#-------------------------------------------------------------

    B[1] = B[1] + AW[1] * T[0]
    B[N] = B[N] + AE[N] * T[N+1]

#-------------------------------------------------------------
#   ASSEMBLE COEFFICIENT MATRIX AND THE RIGHT HAND SIDE VECTOR
#-------------------------------------------------------------

# !!!!!!!!!!!!!!!!CALL TDMA !!!!!!!!!!!!!!!!!!!!

#-------------------------------------------------------------
#   PRINT RESULTS ON THE SCREEN AND ON THE OUTPUT FILE
#-------------------------------------------------------------

