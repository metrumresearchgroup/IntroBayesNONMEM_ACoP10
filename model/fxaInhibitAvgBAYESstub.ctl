$PROB fxaInhIBITAvgBAYES example, Emax model
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

