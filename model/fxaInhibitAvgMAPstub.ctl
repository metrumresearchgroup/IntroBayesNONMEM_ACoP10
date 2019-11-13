$PROB fxaInhIBITAvgBAYES example, Emax model
$INPUT C ID DOSE DV CONC EVID CMT                

$DATA ../../../data/fxaInhAvg.csv     IGNORE=(C='C')

$PRED

MU_1 = THETA(1)
MU_2 = THETA(2)
MU_3 = THETA(3)
MU_4 = THETA(4)

LGT_EMAX = MU_1 + ETA(1)

EMAX = 100 / (1 + EXP(-LGT_EMAX))
EC50 = EXP(MU_2 + ETA(2))
GAMMA = EXP(MU_3 + ETA(3))
SD = EXP(MU_4 + ETA(4))

RESP = EMAX * CONC**GAMMA / (EC50**GAMMA + CONC**GAMMA)
 
Y = RESP + SD * ERR(1) 

