from mpmath import mp, mpf

# Set precision
mp.dps = 50  # You can increase if needed

# Utility to initialize mpmath arrays
def mpzeros(length):
    return [mpf('0') for _ in range(length)]

# --- Initialization ---
NCV = 10
N = NCV

AE = mpzeros(NCV)
AW = mpzeros(NCV)
AP = mpzeros(NCV)
B = mpzeros(NCV)
T = mpzeros(NCV + 2)
TEX = mpzeros(NCV)
X = mpzeros(NCV + 2)
CON = mpzeros(NCV + 2)

ITGS = 0
GSMTOL = mpf('1e-20')  # Very tight tolerance for high precision

# --- Grid and Boundary Conditions ---
T[0] = mpf('10')
T[N + 1] = mpf('210')
X[0] = mpf('0')
X[N + 1] = mpf('10')
S = mpf('1')
DX = (X[N + 1] - X[0]) / N

X[1] = X[0] + DX / 2
for i in range(2, N + 1):
    X[i] = X[i - 1] + DX

# --- Solver Selection and Iteration ---
IS = 2  # 1 - TDMA, 2 - GSM
MAXIT = 1
CON0 = mpf('1')
for i in range(N + 2):
    CON[i] = CON0

# --- TDMA Solver ---
def TDMA():
    for i in range(1, N):
        AP[i] = AP[i] - AW[i] * AE[i - 1] / AP[i - 1]
        B[i] = B[i] + AW[i] * B[i - 1] / AP[i - 1]

    T[N] = B[N - 1] / AP[N - 1]
    for i in range(N - 1, 0, -1):
        T[i] = (B[i - 1] + AE[i - 1] * T[i + 1]) / AP[i - 1]

# --- GSM Solver ---
def GSM():
    global ITGS
    for i in range(1, N + 1):
        T[i] = mpf('0')

    AW[0] = mpf('0')
    AE[N - 1] = mpf('0')

    for ITGS in range(0, 1000):
        RES = mpf('0')
        for i in range(1, N + 1):
            RES += abs(AE[i - 1] * T[i + 1] + AW[i - 1] * T[i - 1] + B[i - 1] - AP[i - 1] * T[i])
        if ITGS == 0:
            RESN = RES
        if RES / RESN < GSMTOL:
            break
        for i in range(1, N + 1):
            T[i] = (AE[i - 1] * T[i + 1] + AW[i - 1] * T[i - 1] + B[i - 1]) / AP[i - 1]

# --- Main Solver Loop ---
NIT = 0
while NIT < MAXIT:
    NIT += 1

    for i in range(N):
        AE[i] = mpf('0.5') * (CON[i + 1] + CON[i + 2]) * S / DX
        AW[i] = mpf('0.5') * (CON[i + 1] + CON[i]) * S / DX

    AW[0] *= mpf('2')
    AE[N - 1] *= mpf('2')

    for i in range(N):
        AP[i] = AW[i] + AE[i]
        B[i] = mpf('0')

    B[0] += AW[0] * T[0]
    B[N - 1] += AE[N - 1] * T[N + 1]

    if IS == 1:
        TDMA()
    elif IS == 2:
        GSM()

    if MAXIT > 1:
        for i in range(1, N + 1):
            CON[i] = CON0 * T[i]

# --- Results ---
print("\nTDMA SOLVER\n" if IS == 1 else f"\nGAUSS-SEIDEL SOLVER: {ITGS} ITERATIONS\n")
print(f"{'I':<3} {'X':<10} {'T':<25} {'TEXACT':<25} {'ERROR':<25}")
ERROR = mpf('0')
for i in range(1, N + 1):
    if MAXIT == 1:
        TEX[i - 1] = ((T[N + 1] - T[0]) * (X[i] - X[0])) / (X[N + 1] - X[0]) + T[0]
    else:
        TEX[i - 1] = mp.sqrt(((T[i + 1] ** 2 - T[0] ** 2) * (X[i] - X[0])) / (X[i + 1] - X[0]) + T[0] ** 2)
    ERR = T[i] - TEX[i - 1]
    ERROR += abs(ERR)
    print(f"{i:<3} {mp.nstr(X[i], 10):<10} {mp.nstr(T[i], 20):<25} {mp.nstr(TEX[i - 1], 20):<25} {mp.nstr(ERR, 10):<25}")

ERROR=ERROR/N
print("\n AVERAGE ERROR : ", ERROR)