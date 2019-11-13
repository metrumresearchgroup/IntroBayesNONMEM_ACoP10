$SIZES LVR=50
$PROBLEM RUN 1 neutropenia PKPD example with Friberg-Karlsson PK-PD

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

