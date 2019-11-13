$SIZES LVR=50
$PROBLEM RUN friberg3.2 neutropenia PKPD example with Friberg-Karlsson PK-PD

$INPUT C ID TIME EVID AMT CMT DVO DV WT DVID
$DATA ../../../data/friberg.csv IGNORE=(C='C')


$ABBR DECLARE INTEGER FIRST_WRITE_PAR
$ABBR DECLARE INTEGER FIRST_WRITE_IPAR

$SUBROUTINES ADVAN13 TOL=6 OTHER = ../../extrasend.f90

$MODEL
NCOMPARTMENTS=8
COMP = (GUT, DEFDOSE)
COMP = (CENT)
COMP = (PERIPH)
COMP = (PROL)
COMP = (TRANS1)
COMP = (TRANS2)
COMP = (TRANS3)
COMP = (CIRC)

$PK
; Request extra information for Bayesian analysis.
; An extra call will then be made for accepted samples 
include '/opt/NONMEM/nm74gf/util/nonmem_reserved_general'
BAYES_EXTRA_REQUEST=1

nThin = THETA(12)

VWT = LOG(WT/70) ; normalized to 70 kg adult

MU_1 = THETA(1) + 0.75*VWT       ; CL
MU_2 = THETA(2) + 0.75*VWT       ; Q
MU_3 = THETA(3) + VWT            ; V1
MU_4 = THETA(4) + VWT            ; V2
MU_5 = THETA(5)                  ; KA
MU_6 = THETA(6)                  ; MTT
MU_7 = THETA(7)                  ; CIRC0
MU_8 = THETA(8)                  ; ALPHA 
MU_9 = THETA(9)                  ; GAMMA
MU_10 = THETA(10)                ; log(SIGMA) additive error PK
MU_11 = THETA(11)                ; log(SIGMA) additive error PD

" CALL EXTRASEND()

CL  = EXP(MU_1 +ETA(1))
Q   = EXP(MU_2 +ETA(2))
V1  = EXP(MU_3 +ETA(3))
V2  = EXP(MU_4 +ETA(4))
KA  = EXP(MU_5 +ETA(5))
MTT = EXP(MU_6 +ETA(6))
CIRC0=EXP(MU_7 +ETA(7))
ALPHA=EXP(MU_8 +ETA(8))

GAMMA= EXP(MU_9 + ETA(9))
SIG1 = EXP(MU_10 + ETA(10)) ; Additive Error PK
SIG2 = EXP(MU_11 + ETA(11)) ; Additive Error PD

; intermediate calculations
K10 = CL/V1
K12 = Q/V1
K21 = Q/V2
KTR = 4/MTT

KPROL = KTR
KCIRC = KTR

; initial PK conditions
A_0(1) = 0
A_0(2) = 0
A_0(3) = 0

; initial PD conditions
A_0(4) = CIRC0
A_0(5) = CIRC0
A_0(6) = CIRC0
A_0(7) = CIRC0
A_0(8) = CIRC0


S2 = V1
S3 = V2


IF(BAYES_EXTRA==1 .AND. NIREC==1 .AND. NDREC==1 .AND. &
     ITER_REPORT>0 .AND. &
     ITER_REPORT/nThin == INT(ITER_REPORT/nThin)) THEN 
IF(FIRST_WRITE_PAR==0) THEN
" OPEN(unit=52,FILE='./par.txt') 
    FIRST_WRITE_PAR=1
  ENDIF
  " WRITE(52,'(I12,1X,50(1X,1PG19.10E3))') ITER_REPORT, &
  " THETA(1), THETA(2), THETA(3), THETA(4), THETA(5), THETA(6), THETA(7), &
  " THETA(8), THETA(9), THETA(10), THETA(11), &  
  " OMEGA(1,1), OMEGA(2,1), OMEGA(2,2), &
  " OMEGA(3,1), OMEGA(3,2), OMEGA(3,3), &
  " OMEGA(4,1), OMEGA(4,2), OMEGA(4,3), OMEGA(4,4), & 
  " OMEGA(5,1), OMEGA(5,2), OMEGA(5,3), OMEGA(5,4), OMEGA(5,5), &
  " OMEGA(6,1), OMEGA(6,2), OMEGA(6,3), OMEGA(6,4), OMEGA(6,5), OMEGA(6,6), & 
  " OMEGA(7,1), OMEGA(7,2), OMEGA(7,3), OMEGA(7,4), OMEGA(7,5), OMEGA(7,6), OMEGA(7,7), & 
  " OMEGA(8,1), OMEGA(8,2), OMEGA(8,3), OMEGA(8,4), OMEGA(8,5), OMEGA(8,6), OMEGA(8,7), OMEGA(8,8)
ENDIF

IF(BAYES_EXTRA==1 .AND. NDREC==1 .AND. ITER_REPORT>0 .AND. &
  ITER_REPORT/nThin == INT(ITER_REPORT/nThin)) THEN 
  IF(FIRST_WRITE_IPAR==0) THEN
    " OPEN(unit=50,FILE='./ipar'//TFI(PNM_NODE_NUMBER)//'.txt') 
    FIRST_WRITE_IPAR=1
  ENDIF
  " WRITE(50,'(I12,1X,F14.0,13(1X,1PG12.5))') ITER_REPORT, ID, &
  " ETA(1), ETA(2), ETA(3), ETA(4), ETA(5), ETA(6), ETA(7), ETA(8)
ENDIF



$DES

DADT(1) = -KA*A(1) 
DADT(2) =  KA*A(1) - (K10 + K12)*A(2) + K21*A(3)
DADT(3) =  K12 *A(2) - K21*A(3) 

CONC = A(2)/V1 ; Drug concentration

EDRUG = ALPHA * CONC
IF(EDRUG > 1.0) EDRUG = 1.0

DADT(4) = KPROL *A(4) *(1 - EDRUG) * ((CIRC0 / A(8)) **GAMMA) - KTR *A(4)
DADT(5) = KTR * (A(4) - A(5))
DADT(6) = KTR * (A(5) - A(6))
DADT(7) = KTR * (A(6) - A(7))
DADT(8) = KTR*A(7) - KCIRC*A(8)


$ERROR (OBSERVATIONS ONLY)
IPRED = LOG(F)
IND = 0
IF(DVID.EQ.2) IND = 1 
YPK = IPRED+SIG1*EPS(1) ; Error for PK
YPD = IPRED+SIG2*EPS(2) ; Error for PD
Y = YPK*IND + YPD*(1-IND)

; Code below into chains

; Initial values of THETA
$THETA
2.3458475542171      ; Log CL = THETA(1)
2.78171156456303      ; Log Q = THETA(2)
3.75999205254463    ; Log V1 = THETA(3)
4.80182772590762    ; Log V2= THETA(4) 
1.1232933486655    ; Log Ka = THETA(5) 
4.81755731950135    ; Log MTT = THETA(6) 
1.51210724744491    ; Log CIRC0 = THETA(7) 
-8.85566751095978    ; Log ALPHA = THETA(8) 
-1.72297421286309    ; Log GAMMA = THETA(9) 
-1.69124241256553    ; Log PK Add Error = THETA(10) 
-1.35071411231885    ; Log PD Add Error = THETA(11) 
(1 FIX)
$OMEGA BLOCK(8) ;INITIAL values of OMEGAs
0.292207235255741
-0.0117007263304952
0.100151663503853
0.000864378508609732
0.000864378508609732
0.372178952161053
-0.0126926903411543
-0.0126926903411543
-0.0126926903411543
0.322822344860368
-0.00505749701238301
-0.00505749701238301
-0.00505749701238301
-0.00505749701238301
0.32017009273741
-0.0167849892262572
-0.0167849892262572
-0.0167849892262572
-0.0167849892262572
-0.0167849892262572
0.509879787251849
-0.00230743618966977
-0.00230743618966977
-0.00230743618966977
-0.00230743618966977
-0.00230743618966977
-0.00230743618966977
0.217594518728487
0.00736046369659033
0.00736046369659033
0.00736046369659033
0.00736046369659033
0.00736046369659033
0.00736046369659033
0.00736046369659033
0.0950651538009127
$OMEGA DIAG(3) ;INITIAL values of unused OMEGAs
0 FIX; GAMMA
0 FIX; Add PK
0 FIX; Add PD
$SIGMA  ;Initial value of SIGMA
1 FIX; Add PK
1 FIX; Add PD
$PRIOR NWPRI
$THETAP          ; Prior information of THETAS
(2.3 FIX)     ; CL Log(10)  
(2.7 FIX)     ; Q  Log(15)
(3.56 FIX)     ; V1 Log(35)
(4.7 FIX)    ; V2 Log(105)
(0.693 FIX)      ; KA log(2)
(4.83 FIX)    ; MTT Log(125)
(1.61 FIX)      ; CIRC0 Log(5) 
(-8.11 FIX)   ; ALPHA Log(3E-4)
(-1.78 FIX)   ; GAMMA Log(0.17)
(-1.60 FIX)   ; log Add error PK Log(0.20)
(-1.60 FIX)   ; log Add error PD Log(0.20)
$THETAPV BLOCK(11)     ;  variances for priors on THETAS (var-cov)
0.25   FIX          ; CL 
0.00   0.25            ; Q
0.00   0.00   0.25        ; V1
0.00   0.00   0.00   0.25      ; V2
0.00   0.00   0.00   0.00   0.25    ; KA 
0.00   0.00   0.00   0.00   0.00   0.04    ; MTT 
0.00   0.00   0.00   0.00   0.00   0.00   0.04    ; CIRC0 
0.00   0.00   0.00   0.00   0.00   0.00   0.00   1   ; ALPHA 
0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  0.04  ; GAMMA  
0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  0.00  1.00  ; Error PK  
0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  0.00  0.00  1.00  ; Error PD  
$OMEGAP BLOCK(8) ;prior information for OMEGA
0.045    FIX          ; CL 
0.00   0.045            ; Q
0.00   0.00   0.045        ; V1
0.00   0.00   0.00   0.045      ; V2
0.00   0.00   0.00   0.00   0.045    ; KA
0.00   0.00   0.00   0.00   0.00   0.045    ; MTT
0.00   0.00   0.00   0.00   0.00   0.00   0.045    ; CIRC0
0.00   0.00   0.00   0.00   0.00   0.00   0.00  0.045   ; ALPHA
$OMEGAPD (11, FIXED)     ; df for OMEGA prior
$ESTIMATION
METHOD=NUTS INTERACTION
AUTO= 2
CTYPE = 0
OLKJDF= 0
OVARF= 1
NUTS_DELTA= 0.8
NBURN = 250 NITER = 250 SEED = 4052
PRINT = 100 NOABORT FILE = /dev/null
$TABLE ID EVID TIME DV IPRED CWRES CWRESI NPDE WT
NOPRINT ONEHEADER FILE=./friberg3.2.tab
$TABLE ID WT CL Q V1 V2 KA MTT CIRC0 ALPHA ETAS(1:LAST)
NOPRINT ONEHEADER FILE=./friberg3.2par.tab
