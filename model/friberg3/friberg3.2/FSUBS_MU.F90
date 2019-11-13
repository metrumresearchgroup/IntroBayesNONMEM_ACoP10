      SUBROUTINE MUMODEL2(THETA,MU_,ICALL,IDEF,NEWIND,&
      EVTREC,DATREC,IREV,NVNT,INDXS,F,G,H,IRGG,GG,NETAS)
      USE NMPRD4P
      USE DECLAREVARIABLES
      USE SIZES,     ONLY: DPSIZE,ISIZE
      USE PRDIMS,    ONLY: GPRD,HPRD,GERD,HERD,GPKD
      USE NMBAYES_REAL,    ONLY: PRIORINFO
      USE PRDATA, ONLY: MXSTEP=>MXSTP01
      USE NMPRD_INT, ONLY: IERPRD,IERPRDU,NETEXT,IQUIT
      USE NMPRD_CHAR,ONLY: ETEXT
      USE ROCM_REAL, ONLY: PR_Y
      USE ROCM_REAL, ONLY: PR_CT
      USE NMPRD_REAL,ONLY: ETA,EPS
      USE NMPR_INT,  ONLY: RPTI=>NRPT_IN,RPTO=>NRPT_OUT,RPTON=>NRPT_ON,PRDFL=>IUSEPRD
      USE ROCM_INT,  ONLY: NREP,IREP=>NCREP
      USE ROCM_INT,  ONLY: MIXNUM=>MIXCALL,MIXEST=>IMIXEST
      USE ROCM_REAL, ONLY: MIXP
      USE NMPRD_INT, ONLY: MSEC=>ISECDER,MFIRST=>IFRSTDER,COMACT,COMSAV,IFIRSTEM
      USE NMPRD_INT, ONLY: MDVRES,ETASXI,NPDE_MODE,NOFIRSTDERCODE
      USE NMPRD_REAL, ONLY: DV_LOQ,CDF_L,DV_LAQ,CDF_LA
      USE NMPR_REAL, ONLY: THSIMP=>THET_P,OMSIMP=>OMEG_P,SGSIMP=>SIGM_P,THSIMPR
      USE ROCM_REAL, ONLY: DEN_NP
      USE ROCM_REAL, ONLY: OMEGA=>VARNF
      USE ROCM_INT,  ONLY: NIREC=>NINDREC,NDREC=>NDATINDR
      USE ROCM_INT,  ONLY: NINDR=>NINDOBS,INDR1=>IDXOBSF,INDR2=>IDXOBSL
      USE ROCM_INT,  ONLY: LIREC=>NDATPASS
      USE NMPRD_INT, ONLY: NPROB,IPROB,IDXSUP,NITR_SUP,NITSUP
      USE NMPRD_INT, ONLY: NTHES_=>NWTHT,NETAS_=>NWETA,NEPSS_=>NWEPS
      USE NMPRD_INT, ONLY: COMACT,COMSAV
      USE ROCM_REAL, ONLY: IIDX=>DIDVALX
      USE ROCM_REAL, ONLY: CNTID=>OFV_IND
      USE CMNM1_INT, ONLY: NIND_7=>NIND
      USE PROCM_REAL,ONLY: TSTATE
      USE PROCM_INT, ONLY: A_0FLG
      USE PROCM_REAL,ONLY: A=>AMNT,DAETA,D2AETA
      USE PRMOD_REAL,ONLY: A_0,DA_0,D2A_0
      USE PROCM_REAL,ONLY: DOSTIM,DDOST,D2DOST
      USE PROCM_REAL,ONLY: DOSREC
      USE PRMOD_INT, ONLY: I_SS,ISSNOW,ISSMOD
      USE NMBAYES_REAL, ONLY: LDF
      USE PROCM_INT, ONLY: MNOW=>MTNOW,MPAST=>MTPAST,MNEXT=>MTNEXT
      USE PKERR_REAL,ONLY: MTIME
      USE PRMOD_INT, ONLY: MTDIFF
      USE NMPRD_INT, ONLY: COMACT,COMSAV
      USE NMPRD_INT, ONLY: U06,U29,UST
      USE ROCM_REAL, ONLY: MIXPT=>MIXP_RAW
      USE NM_INTERFACE, ONLY: TFI
      USE NM_INTERFACE, ONLY: TFD
      USE NMPRD_INT, ONLY: MDV100
      USE NMPRD_INT, ONLY: MDVI1
      USE NMPRD_INT, ONLY: MDVI2
      USE NMPRD_INT, ONLY: MDVI3
      USE NMPRD_INT, ONLY: NTHETA
      USE NMPRD_INT, ONLY: NREC
      USE NMPRD_REAL, ONLY: DB
      USE NMPRD_REAL, ONLY: LB
      USE NMPRD_REAL, ONLY: UB
      USE CMNM1_INT, ONLY: NIND
      USE CMNM2_INT, ONLY: ICOMPRS
      USE CMNM2_REAL, ONLY: SPEC
      USE CMNM2_INT, ONLY: NOBSIND
      USE CMNM2_INT, ONLY: NRECIND
      USE ROCM_INT, ONLY: FIRSTREC
      USE ROCM_INT, ONLY: LASTREC
      USE ROCM_INT, ONLY: FIRSTOBS
      USE ROCM_INT, ONLY: LASTOBS
      USE ROCM_INT, ONLY: FIRSTDOS
      USE ROCM_INT, ONLY: LASTDOS
      USE ROCM_INT, ONLY: EFIRSTREC
      USE ROCM_INT, ONLY: ELASTREC
      USE ROCM_INT, ONLY: EFIRSTOBS
      USE ROCM_INT, ONLY: ELASTOBS
      USE ROCM_INT, ONLY: EFIRSTDOS
      USE ROCM_INT, ONLY: ELASTDOS
      USE CMNM2_INT ,ONLY: IRECIDX
      USE NMPRD_REAL, ONLY: EPS73
      USE NMPRD_REAL, ONLY: REPS
      USE NMPRD_REAL, ONLY: TOL
      USE NMPRD_REAL, ONLY: INFNTY
      USE NMPRD_REAL, ONLY: SQINFN
      USE CMNM7_REAL, ONLY: SREPS
      USE CMNM7_REAL, ONLY: STOL
      USE CMNM4_INT, ONLY: MAXDREC
      USE CMNM4_REAL, ONLY: SMALL
      USE VERSION, ONLY: LEV
      USE NMPRD_INT, ONLY: NPROB
      USE NMPRD_INT, ONLY: IPROB
      USE ROCM_REAL, ONLY: SIGD
      USE NMPRD_INT, ONLY: IPRNV
      USE NMPRD_INT, ONLY: IPS
      USE ROCM_INT, ONLY: NINDREC
      USE ROCM_INT, ONLY: NDATINDR
      USE ROCM_INT, ONLY: NINDOBS
      USE NMPRD_INT, ONLY: IDPROB
      USE CMNM2_INT, ONLY: NOBSIND_MAX
      USE NM_BAYES_INT, ONLY: ISEED
      USE NM_BAYES_INT, ONLY: IRANM
      USE NM_BAYES_REAL, ONLY: MEPS
      USE NMDATA, ONLY: PI
      USE NM_BAYES_INT, ONLY: IKEY_PITER
      USE NM_BAYES_INT, ONLY: IKEY_TERM
      USE NMBAYES_INT, ONLY: ITERATION
      USE NMBAYES_INT, ONLY: ITER_REPORT
      USE NMBAYES_INT, ONLY: NPAT
      USE NMBAYES_INT, ONLY: SAEM_MODE
      USE NMBAYES_INT, ONLY: ITS_MODE
      USE NMBAYES_INT, ONLY: IMP_MODE
      USE NMBAYES_INT, ONLY: PITERPRINT
      USE NMBAYES_REAL, ONLY: OBJI
      USE NMBAYES_CHAR, ONLY: SBANNER
      USE NMBAYES_INT, ONLY: IBMETHOD
      USE NMBAYES_INT, ONLY: BAYES_METHOD
      USE EST_DEFS, ONLY: EST_STANDARD
      USE EST_DEFS, ONLY: EST_DEFAULT
      USE EST_DEFS, ONLY: EST_DIRECT
      USE EST_DEFS, ONLY: EST_BAYES
      USE EST_DEFS, ONLY: EST_ITS
      USE EST_DEFS, ONLY: EST_IMP
      USE EST_DEFS, ONLY: EST_IMPMAP
      USE EST_DEFS, ONLY: EST_SAEM
      USE EST_DEFS, ONLY: EST_CHAIN
      USE NMBAYES_INT, ONLY: EST_COUNTER
      USE NMBAYES_INT, ONLY: IEST_COUNTER
      USE NMBAYES_REAL, ONLY: DHPRIOR
      USE NMBAYES_REAL, ONLY: LDF
      USE NMBAYES_INT, ONLY: ITABLE
      USE NMBAYES_INT, ONLY: BAYES_EXTRA_REQUEST
      USE NMBAYES_INT, ONLY: BAYES_EXTRA
      USE NMBAYES_INT, ONLY: PATCOUNT
      USE NMBAYES_INT, ONLY: CONSTRAIN
      USE NM_BAYES_INT, ONLY: IKEY_SUBJ
      USE NM_BAYES_INT, ONLY: IKEY_BEOPRINT
      USE PNM_CONFIG, ONLY: PNM_RUN_MODE
      USE PNM_CONFIG, ONLY: PNM_NODE_NUMBER
      USE PNM_CONFIG, ONLY: PNM_NODES
      USE PNM_CONFIG, ONLY: PNM_COMPUTERS
      USE PNM_CONFIG, ONLY: PNM_PARAPRINT
      USE PNM_CONFIG, ONLY: PNM_SPLIT
      USE PNM_CONFIG, ONLY: PNM_PARSE_NUMBER
      USE PNM_CONFIG, ONLY: PNM_TRANSFER_TYPE
      USE PNM_CONFIG, ONLY: PNM_TIMEOUTI
      USE PNM_CONFIG, ONLY: PNM_TIMEOUT
      USE PNM_CONFIG, ONLY: PNM_COMMAND
      USE PNM_CONFIG, ONLY: PNM_WORKER_DIRS
      USE PNM_CONFIG, ONLY: PNM_MTOUCH
      USE PNM_CONFIG, ONLY: PNM_MSLEEP
      USE PNM_CONFIG, ONLY: PNM_WTOUCH
      USE PNM_CONFIG, ONLY: PNM_WSLEEP
      USE NMBAYES_CHAR, ONLY: RANMETHOD
      USE NMBAYES_CHAR, ONLY: SRANMETHOD
      USE NMBAYES_CHAR, ONLY: TRANMETHOD
      USE NMBAYES_CHAR, ONLY: CRANMETHOD
      USE NMBAYES_CHAR, ONLY: CCRANMETHOD
      USE NMBAYES_INT, ONLY: ETASTYPE
      USE NMPRD_INT, ONLY: NETAZ
      USE NMPRD_INT, ONLY: NEPSZ
      USE NMBAYES_REAL, ONLY: OMEGANNL
      USE NM_BAYES_INT, ONLY: NM_STEP,BASE_STEP,EST_STEP,COV_STEP,TABLE_STEP,SIML_STEP,INE_STEP,&
&NONP_STEP
      USE NMBAYES_INT, ONLY: MUFIRSTREC
      USE NMBAYES_INT, ONLY: OBJQUICK
      USE CVIDAROOT, ONLY: NRTFN,NRTCT,RTINFO,RTSIGNAL,TRTINFO,YRTINFO,YPRTINFO,GRTINFO
      IMPLICIT REAL(KIND=DPSIZE) (A-Z)
      REAL(KIND=DPSIZE)   :: MU_(*)
      INTEGER NEWIND
      REAL(KIND=DPSIZE) :: EVTREC
 991  FORMAT (35F14.4)
 992  FORMAT (35E15.7)
 998  FORMAT (1PE22.15,1X,1PE22.15/)
      SAVE
      INTEGER(KIND=ISIZE) :: FIRSTEM
      INTEGER(KIND=ISIZE) :: ICALL,IDEF,IREV,NVNT,INDXS,IRGG,NETAS
      DIMENSION :: IDEF(7,*),THETA(*),EVTREC(IREV,*),INDXS(*),GG(IRGG,GPKD+1,*)
      INTEGER(KIND=ISIZE), POINTER :: S1NUM, S2NUM, S1NIT, S2NIT, S1IT, S2IT
      REAL(KIND=DPSIZE), POINTER :: DEN_,CDEN_(:)
! MIXPT=>MIXP_RAW:  Mixture probabilities for data-average block
! TFI:  function converts integer argument i to text: TFI(i)
! TFD:  function converts double precision argument d to text: TFD(d)
! MDV100: =0 if no MDV data itme, =1 if original MDV in data set is 0 or 1, =2 if original MDV in data set is 100 or 101
! MDVI1: DO not ignore MDV>100 records in OBJ routine
! MDVI2: DO not ignore MDV>100 records in OBJ2 routine
! MDVI3: DO not ignore MDV>100 records in OBJ3 routine
! NTHETA: Length of theta; may be set to 1 when length is 0
! DB: Upper bound minus lower bound (UB(I)-LB(I))
! LB: LB(I)=Lower bound for THETA(I)
! UB: UB(I)=Upper bound for THETA(I)
! NIND: Total number of individuals; set in INPT
! ICOMPRS: 0 = usual format for var-covar output, 1 = compressed format
! modified such that 0 = usual format, 9 or fewer elements, 1 = usual format more than 9 elements and 2 = compressed format
! SPEC:  = 1E10, used as a flag in var-covar
! NOBSIND: Number of observation records in individual record
! NRECIND: No. of records / this individual;
! NREC: No. of records in entire dataset;
! number of data records in individual record
! the following 12 items can be compared against NDREC
! FIRSTREC: first record of subject
! LASTREC: last record of subject
! FIRSTOBS: first observation record of subject (for which MDV=0 or 100)
! LASTOBS: last observation record of subject (for which MDV=0 or 100)
! FIRSTDOS: First record of subject with EVID=1 or EVID=4.  FIRSTDOS=-1 when there are no dose records, or PREDPP is not used.
! LASTDOS: Last record of subject with EVID=1 or EVID=4.  LASTDOS=-1 when there are no dose records, or PREDPP is not used.
! EFIRSTREC: first record of subject during estimation (so, among records for which MDV=0 or MDV=1)
! ELASTREC: last record of subject during estimation (so, among records for which MDV=0 or MDV=1)
! EFIRSTOBS: first observation record of subject  during estimation (so, among records for which MDV=0)
! ELASTOBS: last observation record of subject during estimation (so, among records for which MDV=0)
! EFIRSTDOS: First record of subject with EVID=1 or EVID=4 during estimation (so, among records for which MDV=1).  EFIRSTDOS=-1 if no dose record PREDPP is not
! ELASTDOS: Last record of subject with EVID=1 or EVID=4  during estimation (so, among records for which MDV=1).  ELASTDOS=-1 if no dose record PREDPP is not us
! IRECIDX: IRECIDX+1 is starting absolute record number in the data set for the present subject (so while NDREC always sets to 1 for the first record of each su
! During ICALL=1 (INFN initialization), ICALL=3 (INFN finalization), or ICALL=4 (simulation)
! only FIRSTREC and LASTREC are available.  The other values will be set to -1.
! EPS73: Mantissa length.  This is the number 2**(-m), where m is the number of binary digits occurring in the double precision mantissa.
! REPS: Same as REPST; Machine precision. Smallest double precision number that when added to 1 yields a number that is not 1 (same as MEPS)
! TOL: same as TOLT
! INFNTY: SQRT (LARGET); LARGET is largest double precision number representable by the machine
! SQINFN: SQRT (INFNTY)
! SREPS: same as SREPST
! STOL: same as STOLT
! MAXDREC: Maximum number of data record for an individual up to & including this individual
! SMALL: SQRT(TOL*REPS)
! LEV: Level of release
! NPROB: Number of problems
! IPROB: The number of the current problem
! SIGD: minimum value of significant digits; min(DIFA)
! IPRNV: 0 = no printing of NONMEM input info. after 1st iteration of active 1st(2nd) level super problem
! 1 = normal printing for active 1st (2nd) level super problem.
! IPS: Indicator variable,
! 1-> Population data,
! 2-> Single-subject data
! NINDREC (alias NIREC): The number of the individual record at current call (this the subject number)
! NDATINDR (alias NDREC): The number of the data record within the individual record at the current call
! NINDOBS: Number of individual records containing an observation record
! IDPROB: Problem ID
! NOBSIND_MAX: Maximum number of observations for any subject
! ISEED: Present seed number of ranmethods 0-3
! IRANM: Random method for Monte Carlo EM methods
! MEPS: Machine double precision (typically about 1.0E-15)
! PI: Contains the number PI
! IKEY_PITER: =1 to printer iterations, in response to ctrl-I, or signal file iter.sig
! IKEY_TERM: =5 to end run, in response to ctrl-E, or stop.sig file
! =11 to end mode, in response to ctrl-K, or next.sig file
! ITERATION: Present iteration of an EM/Bayesian Method (always non-negative).
! ITER_REPORT: Iteration number that is reported to output (can be negative, if during a burn period).
! NPAT: Number of subjects that have data
! SAEM_MODE: =0 Bayes method
! =1 Stochastic period of SAEM, IMP, ITS, DIRECT, IMPMAP
! =2 reduced stochastic/accumulative period of SAEM
! ITS_MODE: =0 no MAP estimation
! =1 MAP estimation on first iteration
! =2 MAP estimation all iterations (IMPMAP)
! =3 ITS
! IMP_MODE: =0 no Monte Carlo  importance sampling
! =1 Monte Carlo importance sampling
! PITERPRINT: =1 print iteration information
! =0 don’t print iteration information
! OBJI: OBJI(NIREC,1)= objective function value for subject NIREC
! SBANNER: Text to stimation method
! IBMETHOD: Holds METHOD type
! EST_* are defined in module EST_DEFS
! EST_STANDARD (FO/FOCE/LaPLCE)
! EST_DIRECT
! EST_BAYES
! EST_ITS
! EST_SAEM
! EST_IMP
! EST_IMPMAP
! EST_CHAIN
! BAYES_METHOD: ==0 FOR REGULAR BAYES, =1 FOR NUTS
! EST_COUNTER: Total number of $EST statements.
! IEST_COUNTER: Present $EST statement being executed.
! DHPRIOR: Prior portion for EM methods
! LDF: LDF(i)=Loss of degrees of freedom for omega i (usually=1)..  Used in correction for EM analyses only.
! ITABLE: Table number (#TBLN in NONMEM report file, and table numbers in extra output files, such as .ext, .phi, .phm, etc.)
! BAYES_EXTRA_REQUEST: See example8
! BAYES_EXTRA: See example8
! PATCOUNT: Number of subjects with data
! CONSTRAIN: CONSTRAIN setting by which CONSTRAINT routine responds
! IKEY_SUBJ: IF IKEY_SUBJ=20, then ctrl-T or subject.sig signal
! IKEY_BEOPRINT: Toggle to ctrl-B or sending paraprint.sig file.
! PNM_RUN_MODE: May have values:
! PNM_MANAGER
! PNM_WORKER
! PNM_SINGLE
! PNM_NODE_NUMBER: 1=MANAGER,  or SINGLE
! 2-PNM_NODES=WORKER
! PNM_NODES: Number of NODES
! PNM_COMPUTERS: COMPUTERS setting in pnm file
! PNM_PARAPRINT: PARAPRINT setting in pnm file
! PNM_SPLIT: PARSE_TYPE setting in pnm file
! PNM_PARSE_NUMBER: PARSE_NUM
! PNM_TRANSFER_TYPE: TRANSFER_TYPE
! PNM_TIMEOUTI: TIMEOUTI
! PNM_TIMEOUT: TIMEOUT
! PNM_COMMAND: Array pnm commands
! PNM_WORKER_DIRS: Array with worker directories (1-nodes-1)
! PNM_MTOUCH: MTOUCH
! PNM_MSLEEP: MSLEEP
! PNM_WTOUCH: WTOUCH
! PNM_WSLEEP: WSLEEP
! RANMETHOD: $EST RANMETHOD
! SRANMETHOD: $SIML RANMETHOD
! TRANMETHOD: $TABLE RANMETHOD
! CRANMETHOD: $EST METHOD=CHAIN RANMETHOD
! CCRANMETHOD: $CHAIN RANMETHOD
! ETASTYPE: ETASTYPE
! NETAZ: Maximum eta index to have non-zero OMEGA diagonal
! NEPSZ: Maximum eps index to have non-zero SIGMA diagonal
! OMEGANNL: For each x,y in
! $ANNEAL x:y
! OMEGANNL(x)=y
! Those not defined in $ANNEAL are set to -1
! Used in CONSTRAINT.f90
! NM_STEP CAN HAVE THE FOLLOWING VALUES: BASE_STEP, EST_STEP,COV_STEP,TABLE_STEP,SIML_STEP,INE_STEP,NONP_STEP
! MUFIRSTREC: SET TO 1 IN $PK OR $PRED TO CAUSE COVARIATES OF ONLY FIRST RECORD OF SUBJECT TO BE USED IN ASSESSING MU
! OBJQUICK: SET TO 1 IN $PRED OR $PK TO USE QUICK CODE FOR NOUTURN SAMPLING.  ONLY SIMPLE MODELS, WITH NO $MIX, OR $LEVEL
! CVIDAROOT USED FOR ADVAN14 (SEE CVODEU.F90 COMMENTS), OR ADVAN15 (SEE IDAU.F90 COMMENTS) FOR ROOT FINDING
      S1NUM=>IDXSUP(1)
      S2NUM=>IDXSUP(2)
      S1NIT=>NITR_SUP(1)
      S2NIT=>NITR_SUP(2)
      S1IT=>NITSUP(1)
      S2IT=>NITSUP(2)
      DEN_=>DEN_NP(1)
      CDEN_=>DEN_NP(2:)
      FIRSTEM=IFIRSTEM
      IF (ICALL <= 1) THEN
      CALL ASSOCNMPRD4
      IDEF(   1,0001)=  -9
      IDEF(   1,0002)=  -1
      IDEF(   1,0003)=   1
      IDEF(   1,0004)=   0
      IDEF(   2,0003)=   0
      IDEF(   2,0004)=   0
      IDEF(   3,0002)=  12
      IDEF(   3,0003)=  13
      CALL GETETA(ETA)
      IF (IQUIT == 1) RETURN
      RETURN
      ENDIF
      IF (NEWIND /= 2) THEN
      IF (ICALL == 4) THEN
      CALL SIMETA(ETA)
      ELSE
      CALL GETETA(ETA)
      ENDIF
      IF (IQUIT == 1) RETURN
      ENDIF
 !  level            0
      WT=EVTREC(NVNT,009)
      BAYES_EXTRA_REQUEST=1.D0
      nThin=THETA(012)
      B000001=WT/70.D0
      B000002=DLOG(B000001)
      VWT=B000002
      MU_1=THETA(001)+0.75D0*VWT
      MU_(001)=MU_1
      MU_2=THETA(002)+0.75D0*VWT
      MU_(002)=MU_2
      MU_3=THETA(003)+VWT
      MU_(003)=MU_3
      MU_4=THETA(004)+VWT
      MU_(004)=MU_4
      MU_5=THETA(005)
      MU_(005)=MU_5
      MU_6=THETA(006)
      MU_(006)=MU_6
      MU_7=THETA(007)
      MU_(007)=MU_7
      MU_8=THETA(008)
      MU_(008)=MU_8
      MU_9=THETA(009)
      MU_(009)=MU_9
      MU_10=THETA(010)
      MU_(010)=MU_10
      MU_11=THETA(011)
      MU_(011)=MU_11
       RETURN
      CALL EXTRASEND()
      B000003=MU_1+ETA(001)
      B000004=DEXP(B000003)
      CL=B000004
      B000005=MU_2+ETA(002)
      B000006=DEXP(B000005)
      Q=B000006
      B000007=MU_3+ETA(003)
      B000008=DEXP(B000007)
      V1=B000008
      B000009=MU_4+ETA(004)
      B000010=DEXP(B000009)
      V2=B000010
      B000011=MU_5+ETA(005)
      B000012=DEXP(B000011)
      KA=B000012
      B000013=MU_6+ETA(006)
      B000014=DEXP(B000013)
      MTT=B000014
      B000015=MU_7+ETA(007)
      B000016=DEXP(B000015)
      CIRC0=B000016
      B000017=MU_8+ETA(008)
      B000018=DEXP(B000017)
      ALPHA=B000018
      B000019=MU_9+ETA(009)
      B000020=DEXP(B000019)
      GAMMA=B000020
      B000021=MU_10+ETA(010)
      B000022=DEXP(B000021)
      SIG1=B000022
      B000023=MU_11+ETA(011)
      B000024=DEXP(B000023)
      SIG2=B000024
      K10=CL/V1
      K12=Q/V1
      K21=Q/V2
      KTR=4.D0/MTT
      KPROL=KTR
      KCIRC=KTR
      IF(A_0FLG == 1)THEN
      A_0(1)=0.D0
      ENDIF
      IF(A_0FLG == 1)THEN
      A_0(2)=0.D0
      ENDIF
      IF(A_0FLG == 1)THEN
      A_0(3)=0.D0
      ENDIF
      IF(A_0FLG == 1)THEN
      A_0(4)=CIRC0
      ENDIF
      IF(A_0FLG == 1)THEN
      A_0(5)=CIRC0
      ENDIF
      IF(A_0FLG == 1)THEN
      A_0(6)=CIRC0
      ENDIF
      IF(A_0FLG == 1)THEN
      A_0(7)=CIRC0
      ENDIF
      IF(A_0FLG == 1)THEN
      A_0(8)=CIRC0
      ENDIF
      S2=V1
      S3=V2
      B000053=ITER_REPORT/nThin
      B000055=ITER_REPORT/nThin
      B000054=INT(B000055)
      IF(BAYES_EXTRA == 1.D0.AND.NIREC == 1.D0.AND.NDREC == 1.D0.AND. &
      ITER_REPORT >  0.D0.AND.B000053 == B000054)THEN
      IF(FIRST_WRITE_PAR == 0)THEN
      OPEN(unit=52,FILE='./par.txt')
      FIRST_WRITE_PAR=1
      ENDIF
      WRITE(52,'(I12,1X,50(1X,1PG19.10E3))') ITER_REPORT, &
      THETA(1), THETA(2), THETA(3), THETA(4), THETA(5), THETA(6), THETA(7), &
      THETA(8), THETA(9), THETA(10), THETA(11), &
      OMEGA(1,1), OMEGA(2,1), OMEGA(2,2), &
      OMEGA(3,1), OMEGA(3,2), OMEGA(3,3), &
      OMEGA(4,1), OMEGA(4,2), OMEGA(4,3), OMEGA(4,4), &
      OMEGA(5,1), OMEGA(5,2), OMEGA(5,3), OMEGA(5,4), OMEGA(5,5), &
      OMEGA(6,1), OMEGA(6,2), OMEGA(6,3), OMEGA(6,4), OMEGA(6,5), OMEGA(6,6), &
      OMEGA(7,1), OMEGA(7,2), OMEGA(7,3), OMEGA(7,4), OMEGA(7,5), OMEGA(7,6), OMEGA(7,7), &
      OMEGA(8,1), OMEGA(8,2), OMEGA(8,3), OMEGA(8,4), OMEGA(8,5), OMEGA(8,6), OMEGA(8,7), OMEGA(8,8&
&)
      ENDIF
      B000067=ITER_REPORT/nThin
      B000069=ITER_REPORT/nThin
      B000068=INT(B000069)
      IF(BAYES_EXTRA == 1.D0.AND.NDREC == 1.D0.AND.ITER_REPORT >  0.D0 &
      .AND.B000067 == B000068)THEN
      IF(FIRST_WRITE_IPAR == 0)THEN
      OPEN(unit=50,FILE='./ipar'//TFI(PNM_NODE_NUMBER)//'.txt')
      FIRST_WRITE_IPAR=1
      ENDIF
      WRITE(50,'(I12,1X,F14.0,13(1X,1PG12.5))') ITER_REPORT, EVTREC(NVNT,002), &
      ETA(1), ETA(2), ETA(3), ETA(4), ETA(5), ETA(6), ETA(7), ETA(8)
      ENDIF
      P000001=KA
      P000002=K10
      P000003=K12
      P000004=K21
      P000005=V1
      P000006=ALPHA
      P000007=CIRC0
      P000008=GAMMA
      P000009=KPROL
      P000010=KTR
      P000011=KCIRC
      IF (FIRSTEM == 1) THEN
!                      A000098 = DERIVATIVE OF SIG1 W.R.T. ETA(010)
      A000098=B000022
!                      A000103 = DERIVATIVE OF SIG2 W.R.T. ETA(011)
      A000103=B000024
      B000025=1.D0/V1
!                      A000106 = DERIVATIVE OF K10 W.R.T. ETA(001)
      A000106=B000025*B000004
      B000026=-CL/V1/V1
!                      A000107 = DERIVATIVE OF K10 W.R.T. ETA(003)
      A000107=B000026*B000008
      B000033=1.D0/V1
!                      A000122 = DERIVATIVE OF K12 W.R.T. ETA(002)
      A000122=B000033*B000006
      B000034=-Q/V1/V1
!                      A000123 = DERIVATIVE OF K12 W.R.T. ETA(003)
      A000123=B000034*B000008
      B000041=1.D0/V2
!                      A000138 = DERIVATIVE OF K21 W.R.T. ETA(002)
      A000138=B000041*B000006
      B000042=-Q/V2/V2
!                      A000139 = DERIVATIVE OF K21 W.R.T. ETA(004)
      A000139=B000042*B000010
      B000049=-4.D0/MTT/MTT
!                      A000154 = DERIVATIVE OF KTR W.R.T. ETA(006)
      A000154=B000049*B000014
      IF(A_0FLG == 1)THEN
      DA_0(004,007)=B000016
      ENDIF
      IF(A_0FLG == 1)THEN
      DA_0(005,007)=B000016
      ENDIF
      IF(A_0FLG == 1)THEN
      DA_0(006,007)=B000016
      ENDIF
      IF(A_0FLG == 1)THEN
      DA_0(007,007)=B000016
      ENDIF
      IF(A_0FLG == 1)THEN
      DA_0(008,007)=B000016
      ENDIF
!                      A000176 = DERIVATIVE OF S2 W.R.T. ETA(003)
      A000176=B000008
!                      A000178 = DERIVATIVE OF S3 W.R.T. ETA(004)
      A000178=B000010
!                      A000180 = DERIVATIVE OF P000001 W.R.T. ETA(005)
      A000180=B000012
!                      A000182 = DERIVATIVE OF P000002 W.R.T. ETA(003)
      A000182=A000107
!                      A000183 = DERIVATIVE OF P000002 W.R.T. ETA(001)
      A000183=A000106
!                      A000187 = DERIVATIVE OF P000003 W.R.T. ETA(003)
      A000187=A000123
!                      A000188 = DERIVATIVE OF P000003 W.R.T. ETA(002)
      A000188=A000122
!                      A000192 = DERIVATIVE OF P000004 W.R.T. ETA(004)
      A000192=A000139
!                      A000193 = DERIVATIVE OF P000004 W.R.T. ETA(002)
      A000193=A000138
!                      A000197 = DERIVATIVE OF P000005 W.R.T. ETA(003)
      A000197=B000008
!                      A000199 = DERIVATIVE OF P000006 W.R.T. ETA(008)
      A000199=B000018
!                      A000201 = DERIVATIVE OF P000007 W.R.T. ETA(007)
      A000201=B000016
!                      A000203 = DERIVATIVE OF P000008 W.R.T. ETA(009)
      A000203=B000020
!                      A000205 = DERIVATIVE OF P000009 W.R.T. ETA(006)
      A000205=A000154
!                      A000207 = DERIVATIVE OF P000010 W.R.T. ETA(006)
      A000207=A000154
!                      A000209 = DERIVATIVE OF P000011 W.R.T. ETA(006)
      A000209=A000154
      GG(0001,1,1)=P000001
      GG(0001,0006,1)=A000180
      GG(0002,1,1)=P000002
      GG(0002,0002,1)=A000183
      GG(0002,0004,1)=A000182
      GG(0003,1,1)=P000003
      GG(0003,0003,1)=A000188
      GG(0003,0004,1)=A000187
      GG(0004,1,1)=P000004
      GG(0004,0003,1)=A000193
      GG(0004,0005,1)=A000192
      GG(0005,1,1)=P000005
      GG(0005,0004,1)=A000197
      GG(0006,1,1)=P000006
      GG(0006,0009,1)=A000199
      GG(0007,1,1)=P000007
      GG(0007,0008,1)=A000201
      GG(0008,1,1)=P000008
      GG(0008,0010,1)=A000203
      GG(0009,1,1)=P000009
      GG(0009,0007,1)=A000205
      GG(0010,1,1)=P000010
      GG(0010,0007,1)=A000207
      GG(0011,1,1)=P000011
      GG(0011,0007,1)=A000209
      GG(0012,1,1)=S2
      GG(0012,0004,1)=A000176
      GG(0013,1,1)=S3
      GG(0013,0005,1)=A000178
      ELSE
      GG(0001,1,1)=P000001
      GG(0002,1,1)=P000002
      GG(0003,1,1)=P000003
      GG(0004,1,1)=P000004
      GG(0005,1,1)=P000005
      GG(0006,1,1)=P000006
      GG(0007,1,1)=P000007
      GG(0008,1,1)=P000008
      GG(0009,1,1)=P000009
      GG(0010,1,1)=P000010
      GG(0011,1,1)=P000011
      GG(0012,1,1)=S2
      GG(0013,1,1)=S3
      ENDIF
      IF (MSEC == 1) THEN
!                      A000100 = DERIVATIVE OF A000098 W.R.T. ETA(010)
      A000100=B000022
!                      A000105 = DERIVATIVE OF A000103 W.R.T. ETA(011)
      A000105=B000024
      B000029=-1.D0/V1/V1
!                      A000110 = DERIVATIVE OF B000025 W.R.T. ETA(003)
      A000110=B000029*B000008
!                      A000111 = DERIVATIVE OF A000106 W.R.T. ETA(003)
      A000111=B000004*A000110
!                      A000112 = DERIVATIVE OF A000106 W.R.T. ETA(001)
      A000112=B000025*B000004
      B000030=-1.D0/V1/V1
!                      A000113 = DERIVATIVE OF B000026 W.R.T. ETA(001)
      A000113=B000030*B000004
      B000031=CL/V1/V1/V1
!                      A000114 = DERIVATIVE OF B000026 W.R.T. ETA(003)
      A000114=B000031*B000008
      B000032=CL/V1/V1/V1
!                      A000115 = DERIVATIVE OF B000026 W.R.T. ETA(003)
      A000115=B000032*B000008+A000114
!                      A000116 = DERIVATIVE OF A000107 W.R.T. ETA(003)
      A000116=B000008*A000115
!                      A000117 = DERIVATIVE OF A000107 W.R.T. ETA(003)
      A000117=B000026*B000008+A000116
      B000037=-1.D0/V1/V1
!                      A000126 = DERIVATIVE OF B000033 W.R.T. ETA(003)
      A000126=B000037*B000008
!                      A000127 = DERIVATIVE OF A000122 W.R.T. ETA(003)
      A000127=B000006*A000126
!                      A000128 = DERIVATIVE OF A000122 W.R.T. ETA(002)
      A000128=B000033*B000006
      B000038=-1.D0/V1/V1
!                      A000129 = DERIVATIVE OF B000034 W.R.T. ETA(002)
      A000129=B000038*B000006
      B000039=Q/V1/V1/V1
!                      A000130 = DERIVATIVE OF B000034 W.R.T. ETA(003)
      A000130=B000039*B000008
      B000040=Q/V1/V1/V1
!                      A000131 = DERIVATIVE OF B000034 W.R.T. ETA(003)
      A000131=B000040*B000008+A000130
!                      A000132 = DERIVATIVE OF A000123 W.R.T. ETA(003)
      A000132=B000008*A000131
!                      A000133 = DERIVATIVE OF A000123 W.R.T. ETA(003)
      A000133=B000034*B000008+A000132
      B000045=-1.D0/V2/V2
!                      A000142 = DERIVATIVE OF B000041 W.R.T. ETA(004)
      A000142=B000045*B000010
!                      A000143 = DERIVATIVE OF A000138 W.R.T. ETA(004)
      A000143=B000006*A000142
!                      A000144 = DERIVATIVE OF A000138 W.R.T. ETA(002)
      A000144=B000041*B000006
      B000046=-1.D0/V2/V2
!                      A000145 = DERIVATIVE OF B000042 W.R.T. ETA(002)
      A000145=B000046*B000006
      B000047=Q/V2/V2/V2
!                      A000146 = DERIVATIVE OF B000042 W.R.T. ETA(004)
      A000146=B000047*B000010
      B000048=Q/V2/V2/V2
!                      A000147 = DERIVATIVE OF B000042 W.R.T. ETA(004)
      A000147=B000048*B000010+A000146
!                      A000148 = DERIVATIVE OF A000139 W.R.T. ETA(004)
      A000148=B000010*A000147
!                      A000149 = DERIVATIVE OF A000139 W.R.T. ETA(004)
      A000149=B000042*B000010+A000148
      B000051=4.D0/MTT/MTT/MTT
!                      A000156 = DERIVATIVE OF B000049 W.R.T. ETA(006)
      A000156=B000051*B000014
      B000052=4.D0/MTT/MTT/MTT
!                      A000157 = DERIVATIVE OF B000049 W.R.T. ETA(006)
      A000157=B000052*B000014+A000156
!                      A000158 = DERIVATIVE OF A000154 W.R.T. ETA(006)
      A000158=B000014*A000157
!                      A000159 = DERIVATIVE OF A000154 W.R.T. ETA(006)
      A000159=B000049*B000014+A000158
      IF(A_0FLG == 1)THEN
      D2A_0(004,007,007)=B000016
      ENDIF
      IF(A_0FLG == 1)THEN
      D2A_0(005,007,007)=B000016
      ENDIF
      IF(A_0FLG == 1)THEN
      D2A_0(006,007,007)=B000016
      ENDIF
      IF(A_0FLG == 1)THEN
      D2A_0(007,007,007)=B000016
      ENDIF
      IF(A_0FLG == 1)THEN
      D2A_0(008,007,007)=B000016
      ENDIF
!                      A000177 = DERIVATIVE OF A000176 W.R.T. ETA(003)
      A000177=B000008
!                      A000179 = DERIVATIVE OF A000178 W.R.T. ETA(004)
      A000179=B000010
!                      A000181 = DERIVATIVE OF A000180 W.R.T. ETA(005)
      A000181=B000012
!                      A000184 = DERIVATIVE OF A000182 W.R.T. ETA(003)
      A000184=A000117
!                      A000185 = DERIVATIVE OF A000183 W.R.T. ETA(001)
      A000185=A000112
!                      A000186 = DERIVATIVE OF A000183 W.R.T. ETA(003)
      A000186=A000111
!                      A000189 = DERIVATIVE OF A000187 W.R.T. ETA(003)
      A000189=A000133
!                      A000190 = DERIVATIVE OF A000188 W.R.T. ETA(002)
      A000190=A000128
!                      A000191 = DERIVATIVE OF A000188 W.R.T. ETA(003)
      A000191=A000127
!                      A000194 = DERIVATIVE OF A000192 W.R.T. ETA(004)
      A000194=A000149
!                      A000195 = DERIVATIVE OF A000193 W.R.T. ETA(002)
      A000195=A000144
!                      A000196 = DERIVATIVE OF A000193 W.R.T. ETA(004)
      A000196=A000143
!                      A000198 = DERIVATIVE OF A000197 W.R.T. ETA(003)
      A000198=B000008
!                      A000200 = DERIVATIVE OF A000199 W.R.T. ETA(008)
      A000200=B000018
!                      A000202 = DERIVATIVE OF A000201 W.R.T. ETA(007)
      A000202=B000016
!                      A000204 = DERIVATIVE OF A000203 W.R.T. ETA(009)
      A000204=B000020
!                      A000206 = DERIVATIVE OF A000205 W.R.T. ETA(006)
      A000206=A000159
!                      A000208 = DERIVATIVE OF A000207 W.R.T. ETA(006)
      A000208=A000159
!                      A000210 = DERIVATIVE OF A000209 W.R.T. ETA(006)
      A000210=A000159
      GG(0001,0006,0006)=A000181
      GG(0002,0002,0002)=A000185
      GG(0002,0004,0002)=A000186
      GG(0002,0004,0004)=A000184
      GG(0003,0003,0003)=A000190
      GG(0003,0004,0003)=A000191
      GG(0003,0004,0004)=A000189
      GG(0004,0003,0003)=A000195
      GG(0004,0005,0003)=A000196
      GG(0004,0005,0005)=A000194
      GG(0005,0004,0004)=A000198
      GG(0006,0009,0009)=A000200
      GG(0007,0008,0008)=A000202
      GG(0008,0010,0010)=A000204
      GG(0009,0007,0007)=A000206
      GG(0010,0007,0007)=A000208
      GG(0011,0007,0007)=A000210
      GG(0012,0004,0004)=A000177
      GG(0013,0005,0005)=A000179
      ENDIF
      RETURN
      END
