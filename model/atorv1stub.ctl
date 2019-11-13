$PROBLEM RUN# atorv1
$INPUT C ID TIME AMT DV CMT DOSE EVID MDV ADDL II SS FORM WT PER BLQ
;$INPUT C ID SID=DROP STDY=DROP TIME DV DVR=DROP DVC=DROP AMT AMTO=DROP CMT TYPE=DROP
;       DOSE DOSN=DROP NTIM=DROP EVID MDV ADDL II SS FORM WT SEX=DROP AGE
;       HT=DROP RACE=DROP PER BLQ TRCD=DROP ANLT=DROP TXT=DROP 
;$DATA ../../../data/atorv1172atvb.csv IGNORE=(C='C', BLQ.EQ.1)
$DATA ../../../data/atorvWrkShop.csv IGNORE=(C='C', BLQ.EQ.1)

$ABBR DECLARE INTEGER FIRST_WRITE_PAR
$ABBR DECLARE INTEGER FIRST_WRITE_IPAR

$SUBROUTINES ADVAN4 TRANS4 OTHER = ../../extrasend.f90

$PK
; Request extra information for Bayesian analysis.
; An extra call will then be made for accepted samples 
include '/opt/NONMEM/nm74gf/util/nonmem_reserved_general'
BAYES_EXTRA_REQUEST=1

nThin = THETA(8)

VWT = LOG(WT/70) ; normalized to 70 kg adult
MU_1 = THETA(1) + VWT * 0.75 ; CL
MU_2 = THETA(2) + VWT        ; V2
MU_3 = THETA(3) + VWT * 0.75 ; Q
MU_4 = THETA(4) + VWT        ; V3
MU_5 = THETA(5)              ; ka
MU_6 = THETA(6)              ; F1chew
MU_7 = THETA(7)              ; SD

" CALL EXTRASEND()

;ATORVASTATIN

CL = EXP(MU_1 + ETA(1)) ; atv cl 
V2 = EXP(MU_2 + ETA(2))
Q = EXP(MU_3 + ETA(3))
V3 = EXP(MU_4 + ETA(4))
;deltaKA = EXP(MU_5 + ETA(5))
KA = EXP(MU_5 + ETA(5))
F1chew = EXP(MU_6 + ETA(6))
SD = EXP(MU_7 + ETA(7))

; Constrain KA > lambda1 (typical value) to enhance identifiability
;T1 = CL/V2
;T23 = Q/V2
;T32 = Q/V3
;L1 = ((T1+T23+T32)+SQRT((T1+T23+T32)**2-4*T1*T32))/2
;L2 = ((T1+T23+T32)-SQRT((T1+T23+T32)**2-4*T1*T32))/2

;KA = deltaKA + L2
S2 = V2/1000              ; DOSE IN uM & CONC IN nM/L
F1 = 1.0             ; TABLET AS REFERENCE
IF(FORM.EQ.2) F1 = F1chew  ; CHEWABLE F1

IF(BAYES_EXTRA==1 .AND. NIREC==1 .AND. NDREC==1 .AND. &
   ITER_REPORT>0 .AND. &
   ITER_REPORT/nThin == INT(ITER_REPORT/nThin)) THEN 
  IF(FIRST_WRITE_PAR==0) THEN
    " OPEN(unit=52,FILE='./par.txt') 
    FIRST_WRITE_PAR=1
  ENDIF
  " WRITE(52,'(I12,1X,20(1X,1PG19.10E3))') ITER_REPORT, &  
  " OMEGA(1,1), OMEGA(2,1), OMEGA(2,2), &
  " THETA(1), THETA(2), THETA(3), THETA(4), THETA(5), THETA(6), THETA(7)
ENDIF
IF(BAYES_EXTRA==1 .AND. NDREC==1 .AND. ITER_REPORT>0 .AND. &
  ITER_REPORT/nThin == INT(ITER_REPORT/nThin)) THEN 
  IF(FIRST_WRITE_IPAR==0) THEN
    " OPEN(unit=50,FILE='./ipar'//TFI(PNM_NODE_NUMBER)//'.txt') 
    FIRST_WRITE_IPAR=1
  ENDIF
  " WRITE(50,'(I12,1X,F14.0,10(1X,1PG19.10E3))') ITER_REPORT, ID, &
  " ETA(1), ETA(2)
ENDIF


$ERROR

LOQ = 0.45          ; nM/L
IPRED = F
DUM = (LOQ-IPRED)/(SD*IPRED)
CUMD=PHI(DUM)
IF(BLQ.EQ.0)THEN
  F_FLAG=0
  Y = F*(1+SD*ERR(1)) ; ATV
ELSE
  F_FLAG=1
  Y=CUMD
ENDIF

; Code below into chains

