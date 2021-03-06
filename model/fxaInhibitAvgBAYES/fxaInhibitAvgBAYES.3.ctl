$PROB fxaInhibitAvgBAYES.3 example, Emax model
$INPUT C ID DOSE DV CONC EVID CMT               
$DATA ../../../data/fxaInhAvg.csv     IGNORE=(C='C')

$SUBROUTINES OTHER = ../../extrasend.f90

$PRED
; Request extra information for Bayesian analysis.
; An extra call will then be made for accepted samples 
include '/opt/NONMEM/nm74gf/util/nonmem_reserved_general'
BAYES_EXTRA_REQUEST=1

nThin = THETA(5)

MU_1 = THETA(1)
MU_2 = THETA(2)
MU_3 = THETA(3)
;MU_4 = THETA(4)

" CALL EXTRASEND()

LGT_EMAX = MU_1 + ETA(1)

EMAX = 100 / (1 + EXP(-LGT_EMAX))
EC50 = EXP(MU_2 + ETA(2))
GAMMA = MU_3 + ETA(3)
SD = EXP(THETA(4))

RESP = EMAX * CONC**GAMMA / (EC50**GAMMA + CONC**GAMMA)
 
Y = RESP + SD * ERR(1) 

IF(BAYES_EXTRA==1 .AND. NIREC==1 .AND. NDREC==1 .AND. &
     ITER_REPORT > 0 .AND. &
     ITER_REPORT/nThin == INT(ITER_REPORT/nThin)) THEN 
IF(FIRST_WRITE_PAR==0) THEN
" OPEN(unit=52,FILE='./par.txt') 
    FIRST_WRITE_PAR=1
  ENDIF
  " WRITE(52,'(I12,1X,5(1X,1PG19.10E3))') ITER_REPORT, &
  " THETA(1), THETA(2), THETA(3), THETA(4) 
ENDIF

; Initial values of THETA
$THETA
2.26996133957687      ; EMAX = 100 * EXPIT(THETA(1))
4.76295777215047  ; EC50 = EXP(THETA(2))
(0,  1.86712858320982 )   ; GAMMA = THETA(3)
1.82466541224132   ; SD = EXP(THETA(4)) 
(1 FIX)
$OMEGA DIAG(3);INITIAL values of OMEGAs
(0.0 FIX)
(0.0 FIX)
(0.0 FIX)
$SIGMA ;INITIAL values of SIGMAs
(1 FIX)


$PRIOR NWPRI
$THETAP          ; Prior information of THETAS
(0 FIX)      ;  THETA(1) EMAX
(5.5 FIX)      ;  THETA(2) EC50
(1 FIX)      ;  THETA(3) gamma
(2.3 FIX)      ; THETA(4) SD
$THETAPV BLOCK(4)     ;  variances for priors on THETAS (var-cov)
2.0 FIX ; EMAX almost uniform
0.00 1  ; EC50
0.00 0.00 1  ; gamma
0.00 0.00 0.00 1   ; SD
$ESTIMATION
METHOD = BAYES INTERACTION LAPLACIAN
AUTO= 1
CTYPE = 0
OLKJDF= 0
OVARF= 1
NUTS_DELTA= 0.8
NBURN = 10000 NITER = 10000 SEED = 8663
PRINT = 100 NOABORT FILE = /dev/null
$TABLE ID EVID CONC DV PRED
NOPRINT ONEHEADER FILE=./fxaInhibitAvgBAYES.3.tab
$TABLE ID EMAX EC50 GAMMA SD
NOPRINT ONEHEADER FILE=./fxaInhibitAvgBAYES.3par.tab
