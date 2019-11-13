!*********************************COPYRIGHT******************************************
!                                                                                   !
!       THE NONMEM SYSTEM MAY BE DISTRIBUTED ONLY BY ICON DEVELOPMENT               !
!       SOLUTIONS.                                                                  !
!                                                                                   !
!       COPYRIGHT BY ICON DEVELOPMENT SOLUTIONS                                     !
!       2009-2017 ALL RIGHTS RESERVED.                                              !
!                                                                                   !
!       DO NOT ATTEMPT TO MODIFY CODE WITHOUT FIRST CONSULTING WITH                 !
!       ICON DEVELOPMENT SOLUTIONS.                                                 !
!                                                                                   !
!************************************************************************************
!
!-----------------------------HISTORY------------------------------------------------
! VERSION     : NONMEM VII
! AUTHOR      : ALISON J. BOECKMANN
! CREATED ON  : SEP/1986
! LANGUAGE    : FORTRAN 90/95
! LAST UPDATE : FEB/1991 - NUMERIC DERIVATIVES ELIMINATED 
!               JUL/2008 - COMMON BLOCKS REPLACED WITH MODULES
!               NOV/2008 - INTRODUCED HEADER INFORMATIONS AND RESTRUCTURED AS PER
!                          THE NONMEM STANDARDS
!               JAN/2009 - ASSUMED SIZE ARRAY'S REPLACED WITH ASSUMED SHAPE
!               FEB/2009 - MOVED DATA STATEMENTS TO PRDATA.F90 FILE
!               APR/2009 - INTRODUCED ERROR CHECK FOR FILE OPERATION
!               FEB/2010 - CHANGED SIZES to PRSIZES
!               FEB/2011 - INTEGRATED 7.2BETA5.8B MODIFICATIONS
!               NOV/2016 - INTEGRATED NONMEM7.4 CHANGES.  Using now ADVAN8 as template, rather than ADVAN9
!
!----------------------------- ADVAN13.f90 -------------------------------------------
!
! SUBROUTINE ADVAN(ICALL)
!
! DESCRIPTION : ADVAN6 and ADVAN13 are routines in PREDPP's library which implement the
!               general non-linear model. The general non-linear model is used for
!               systems in which a drug is distributed between compartments according
!               to first-order processes. ADVAN13 may be used in preference to ADVAN6
!               when the differential equations describing the processes are stiff.
!
! ARGUMENTS   : ICALL
!               IN     - ICALL
!                        ICALL - Values are 0,1,2 or 3. ICALL values 0,1,2,3 from NONMEM
!                                are unchanged. ICALL values >= 4 are converted to 2 for
!                                this routine.
!               OUT    - NONE
!               IN OUT - NONE
!
! CALLED BY   : FCN2   - FCN2 copies Y to the state vector and calls ADVAN to advance the
!                        system (SS6 already set up infusion array and DT). It inputs the
!                        SS dose over the advance interval. It then computes the 
!                        difference between the initial state vector and the post-advance
!                        state vector for multiple doses-interfaces with ADVAN.
!               SADVAN - Stands for Supervisor ADVAN.  This is an interface.
!                        PREDPP contains a library of routines, called ADVAN routines,
!                        which implement specific kinetic models. Exactly one ADVAN
!                        routine must be selected for each NONMEM / PREDPP run.
!                        Its function is to "ADVANCE" the kinetic system from one
!                        state time point to the next
!               SS6    - General routine for all models
!
! CALLS       : DES             - To solve diffrential equations
!               MODEL           - Model definition information
!               TOL             - The error which is controlled by way of the parameter TOL is
!                                 an estimate of the local truncation error, that is, the error
!                                 committed on taking a single step with the method, starting
!                                 with data regarded as exact.
!               lsoda
!               RUN_TIME_ERRORS - Set run-time error conditions
!               FCN1            - Evaluate YPRIME(1),...,YPRIME(N). FCN1 calls DES to obtain DADT. It
!                                 then adds in the rate of any infusions that are in progress. 
!               ERRORMSGS       - Writes error messages to a file JUNIT and sets IQUIT to 1 indicating
!                                 that NONMEM has to quit for non-super problems. And for super problems
!                                 calculation continues with next sub problem.
!
! ALGORITHM   : - IF (ICALL == 2) Normal Entry; else
!                 - IF (ICALL /= 0) allow initialization call to DES; RETURN; else
!                 - Initialization call (ICALL=0)
!                 - Set maximum values
!                 - Call for MODEL definition information - Call MODEL
!                 - Obtain the no. of significant digits - Call TOL
!                 - RETURN
!               - Normal entry
!                 - Compute for not a special call
!                 - Set I2FLAG=0. Set I2FLAG=1 if there will be partials wrt time due to
!                   term.infusion
!                 - HR(I,1)=Infusion rate into A(I) after terminating inf. ends
!                 - HR(I,2)=Infusion rate into A(I) before it ends
!                 - Save old values if needed for output compartment
!                 - SS6 sets MFLAG8=1 when the call is from ZSPOW
!                 - Set up now for first integration interval
!                 - Only FCN is called when MIT=2, so pass same name for FCN and FCNJ
!                 - If error in integration to IDA(ISPEC) or DT, call RUN_TIME_ERRORS.
!                 - All finished if no terminating infusion
!                 - Branch if numeric derivatives are not to be obtained
!                 - Adjustment for end point - FCN1 handles ON/OFF, adds both infusions
!                 - Compute for second set of derivatives
!                 - Set up for second integration interval
!                 - If error in integration to IDA(ISPEC) or DT, call RUN_TIME_ERRORS.
!                 - Finish integrating to TSTART+DT
!                 - Go back for another group of ETAS if necessary
!                 - Must also subtract IDEA from DTINT before re-doing the integ.
!                 - For call from SS, must return all the first derivatives
!                 - Adjust derivatives at end point - FCN1 handles ON/OFF, adds infusion
!                 - Any value greater than NCM1 will do
!                 - If error in integration to IDA(ISPEC) or DT, call RUN_TIME_ERRORS.
!                 - Calculate the output compartment
!               - Exit point
!
! MODULES USED: PRSIZES,PRDATA,NMPRD_INT,PRCOM_INT,PRCM_INT,PRCOM_REAL,PROCM_REAL,
!               PRCM_REAL,NMPRD_CHAR,PRCOM_LOG,PRCM_LOG,NM_INTERFACE,PR_INTERFACE
!
! CONTAINS    : RUN_TIME_ERRORS
!
! LOCAL'S     : CTEMP,CTEMPC,DADT,DADT2,H,I,I2FLAG,ICOMP,IER,IND,IOUT,IWK,J,K,KLOOP,L,
!               METH,METHOD,MIT,NPR,NY,OLDA,OLDAE2,OLDAET,PE,SAVSEC,TL,
!               WK,X,XEND,Y
!
!---------------------------- END OF HEADER -----------------------------------------
!
      SUBROUTINE ADVAN(ICALL)
!      
      USE PRSIZES,      ONLY: ISIZE,DPSIZE,PE,PC,MAXFCN,PG,PIR
      USE PRDATA,       ONLY: TEN,P00001
! INTEGER
      USE NMPRD_INT,    ONLY: IERPRD,IFIRSTEM,IFRSTDER
      USE PRCOM_INT,    ONLY: ADVID,BETA,CALCTR,IATT,IDNO,IH,IMAX,ISPEC,LINK,LOGUNT,    &
                              MAXCAL,MCOMP,MITER,NBP,NCM1,NPETAS2=>NPETAS,NRD,SV
      USE PRCM_INT,     ONLY: INRD,MFLAG8,LST,LEND,SLEND,SLST,SPW,SLOOP,IANRD,ITOTL,MAPINV, &
                              MFLAG1,MFLAG2,MFLAG3, &
                              MFLAG4,MFLAG5,MFLAG6,MFLAG7,MFLAG8,XLNCM1,XNCM1,IANRD,SPW, &
                              LST,LEND,SLEND,SLST,SLOOP
! REAL
      USE PRCOM_REAL,   ONLY: D2ADTE,D2TINT,D2TSTA,H2ETA,I2REA,I2DEA,R2E,ADTE,DA,DET,DP,&
                              DT,DTINT,DTSTAR,HETA,HR,IDA,IDEA,IRA,IREA,R,RE,TSTART,ZERO,G3,ONE
      USE PROCM_REAL,   ONLY: AMNT,DAETA,D2AETA
      USE PRCM_REAL,    ONLY: DAC,DPC,DTC
      USE PRMOD_INT,    ONLY: IPRINT
! CHARACTER
      USE NMPRD_CHAR,   ONLY: ETEXT
! LOGICAL
      USE PRCOM_LOG,    ONLY: NOETAS2=>NOETAS,SECOND,DTIME,LFLAG,DTDER
      USE PRCM_LOG,     ONLY: COMPAC,DIDCAA,DIDDES,DOFINL,MAPPED
! INTERFACE
      USE NM_INTERFACE, ONLY: ERRORMSGS
      USE PR_INTERFACE, ONLY: DES,JAC,MODEL,TOL,FCN1
      USE NM_BAYES_INT,    ONLY: TOLTOPRED,ADVID2
!      
      IMPLICIT NONE
!
      INTEGER(KIND=ISIZE), INTENT(IN) :: ICALL
! 
      SAVE
!
!------------------------------------------------------------------------------------
!     COMMON /NMPRD1/ IERPRD,NETEXT
!     COMMON /NMPRD2/ ETEXT(3)
!     INTEGER IERPRD,NETEXT
!     CHARACTER*132 ETEXT
!     COMMON /PRCOM0/ NP,NBP,YFORM
!     COMMON /PRCOM0/ MAXKF,IFORM
!     COMMON /PRCOM0/ IDC,IDO,MAXIC,ISV,IINST,ITURN
!     COMMON /PRCOM0/ JTIME,JCONT,JEVENT,JAMT,JRATE,JSS,JDELTA
!     COMMON /PRCOM0/ JCOMPT,JCOMPF,JERROR,SSC,KREC,JMORE,JDUM
!     COMMON /PRCM02/ ITOTL,INRD
!     INTEGER ITOTL(PC),INRD(PM)
!     COMMON /PRCOM1/ NOETAS,SECOND
!     COMMON /PRCOM3/ ITRANS,IRGG,IREV,NPETAS,NPEPS
!     COMMON /PRCOM4/ G3,HH,DELTA,DT,DTE
!     COMMON /PRCOM4/ YMULT,ZERO,ONE,XR,XD,TSTART,DTSTAR
!     COMMON /PRCOM4/ DDELTA,D2DELT,ADTE,D2ADTE
!     COMMON /PRCOM5/ ISPEC,DCTR,BETA,DD
!     COMMON /PRCOM5/ IP,IPOOL,IHEAD,INEXT,IBACK,SV
!     COMMON /PRCOM6/ IA,IAA,IAEA,IRA,IREA,IDA,IDEA,R,RE
!     COMMON /PRCOM6/ RHO,RHOE,SDEL,SDELE,SSA,SSAE,SSR,SSRE
!     COMMON /PRCOM6/ SAMT,SDEL1
!     COMMON /PRCOM6/ I2AEA,I2REA,I2DEA,R2E,D2DTE,D2TSTA
!     COMMON /PRCOM7/ ADVID,SSID
!     COMMON /PRCOMN/ LOGUNT,NC
!     DOUBLE PRECISION DELTA,G3,HH
!     DOUBLE PRECISION SDELE,RHOE,SSAE,SSRE,YMULT,ZERO,XR,XD
!     DOUBLE PRECISION ONE,TSTART,DTSTAR(PE)
!     DOUBLE PRECISION DDELTA(PE),D2DELT(PE,PE),ADTE(PE),D2ADTE(PE,PE)
!     DOUBLE PRECISION IA(90),IAA(90),IAEA(90,PE),IRA(90),IDA(90)
!     DOUBLE PRECISION IREA(90,PE),IDEA(90,PE),R(PC),RE(PC,PE)
!     DOUBLE PRECISION I2REA(90,PE,PE),I2DEA(90,PE,PE),I2AEA(90,PE,PE)
!     DOUBLE PRECISION R2E(PC,PE,PE),D2DTE(PE,PE)
!     DOUBLE PRECISION D2TSTA(PE,PE)
!     DOUBLE PRECISION DT,DTE(PE),RHO,SDEL,SSA,SSR
!     DOUBLE PRECISION SAMT,SDEL1
!     DIMENSION SDELE(PE),RHOE(PE),SSAE(PE),SSRE(PE)
!     DIMENSION G3(PG+1,PE+1,PE+1),HH(PE,PE)
!     INTEGER LOGUNT,IRGG,IREV,ITRANS,NPETAS,NPEPS
!     INTEGER JCONT,JTIME,JEVENT,JAMT,JRATE,JSS,JDELTA
!     INTEGER JCOMPT,JCOMPF,JERROR
!     INTEGER NC,IDC,IDO,NP,NBP,SSC,KREC,JMORE,JDUM
!     INTEGER ISV(PC),SV(PC)
!     INTEGER IINST(PC),ITURN(PC)
!     INTEGER ADVID,SSID,MAXKF,IFORM(PG+1),YFORM,MAXIC
!     INTEGER BETA(90),IPOOL(90),IP,IHEAD,INEXT(90),IBACK(90)
!     INTEGER ISPEC,DD(90),DCTR
!     LOGICAL NOETAS,SECOND
!     COMMON /PRCMX1/ GENMOD,MAPPED,COMPAC
!     LOGICAL GENMOD,MAPPED,COMPAC
!     COMMON /PROCM4/ A,DAETA,D2AETA
!     DOUBLE PRECISION A,DAETA,D2AETA
!     DIMENSION A(PC),DAETA(PC,PE),D2AETA(PC,PE,PE)
!     COMMON /PRCOMD/ DA,DP,DET,HR,HETA,H2ETA
!     COMMON /PRCOME/ NRD,MCOMP,NCM1,IH
!     COMMON /PRCOML/ IDNO,IATT,LINK,ILINK,INUM
!     DOUBLE PRECISION Y !,DA !,DP,HR,HETA,DET,H2ETA
!     DIMENSION Y(P99) ! ,DA(PM,PM),DP(PM,PG),HR(PM,2),DET(PM)
!     DIMENSION HETA(PM,PE,2),H2ETA(PM,PE,PE,2)
!     INTEGER MCOMP,NCM1,IH
!     INTEGER NY,IND,IER,NPR !,NRD,IDNO,,IATT,LINK,ILINK,INUM
!     DIMENSION IATT(PM,9),LINK(PM,PC),ILINK(PM,PC),INUM(PM)
!     COMMON /PRCOMT/ DTIME,LFLAG,DTDER
!     LOGICAL DTIME,LFLAG,DTDER
!     COMMON /PRCOMF/ MAXCAL,CALCTR
!     INTEGER MAXCAL,CALCTR
! FOR LOOPS ON ETAS IN ADVAN6, ADVAN13, SS6
!     COMMON /PRCMET/ LST,LEND,SLEND(PE),SLST(PE),SPW,SLOOP
!     INTEGER LST,LEND,SLEND,SLST,SPW,SLOOP
! THE COMPACT ARRAYS THEMSELVES
!     COMMON /PRCMDC/ DAC(PIR,1),DPC(PIR,1),DTC(PIR)
!     DOUBLE PRECISION DAC,DPC,DTC
!     COMMON /PRCOMG/ MITER,IDUM1,IMAX,ISTFLG,INTFLG
!     INTEGER MITER,IDUM1,IMAX,ISTFLG,INTFLG
!     COMMON /PRCOMX/ DTINT(PE),D2TINT(PE,PE)
!     DOUBLE PRECISION DTINT,D2TINT
!     COMMON /PRCMLX/ MFLAG1,MFLAG2,MFLAG3,MFLAG4,MFLAG5,MFLAG6
!     COMMON /PRCMLX/ MFLAG7,MFLAG8,SS3
!     INTEGER MFLAG1,MFLAG2,MFLAG3,MFLAG4,MFLAG5,MFLAG6
!     INTEGER MFLAG7,MFLAG8,SS3
!     COMMON /PRCM03/ DIDCAA,DIDDES,DIDAES,DOFINL
!     LOGICAL DIDCAA,DIDDES,DIDAES,DOFINL
!------------------------------------------------------------------------------------
!
! Local Variables
!
      INTEGER(KIND=ISIZE), PARAMETER :: PM=PC-1,PPW=2*PE*PM+PM,P99=PC*PE+PC-PE-1,P11=PE+1
      INTEGER(KIND=ISIZE) :: I,I2FLAG,ICOMP,IER,IND,IOUT,J,K,KLOOP,L,METH,METHOD,NPR,NY
!      
      REAL(KIND=DPSIZE)   :: DADT,DADT2,H,OLDA,OLDAE2,OLDAET,TL,X,XEND,Y

      REAL(KIND=DPSIZE)   :: RWORK(22+9*P99+P99*P99)

      INTEGER(KIND=ISIZE) :: ISTATE,ITASK,IWORK(20+P99),LIW,LRW,MF,MFLAGA,MFLAGX,&
                             MFLAGY,MFLAGZ,IREDO


!
      CHARACTER(LEN=4)    :: IST
      CHARACTER(LEN=3)    :: CTEMP
      CHARACTER(LEN=12)   :: CTEMPC
!
      LOGICAL :: SAVSEC
      INTEGER NPETAS
      LOGICAL NOETAS
      INTEGER JTOTL,ITOL,IOPT
      EXTERNAL JTOTL,RES2
!
      DIMENSION :: Y(P99),OLDA(PM),OLDAET(PM,PE),OLDAE2(PM,PE,PE),DADT(PPW),DADT2(PPW)

      REAL(KIND=DPSIZE)   :: ZNRD(0:PM),ZANRD(0:PM),ATOL(P99),RTOL(P99),AOTOL(P99),OLDDT
!      
! Variables used to advance by DT; For infusions into primary compartment
! which continue past DT:
!   R(IDC)        = Infusion rate into primary compartment
!   RE(IDC,K)     = Derivative of rate w.r.t ETA(K)
! If an infusion into primary compt. is terminating, ISPEC>0 and:
!   IRA(ISPEC)    = Its infusion rate
!   IREA(ISPEC,K) = Derivative of IRA w.r.t ETA(K)
!   IDA(ISPEC)    = Its duration
!   IDEA(ISPEC,K) = Derivative of IDEA w.r.t ETA(K)
!   BETA(ISPEC)   = Its compartment number.
!   TSTART        = Starting time of interval over which to advance
! Variables supplied by user's subroutine DES:
!   DADT(I)       = Derivative of A(I) w.r.t time T
!   DA(I,K)       = Derivative of DADT(I) w.r.t A(K)
!   DP(I,J)       = Derivative of DADT(I) w.r.t GG(J)
!
! DADT must be larger than PM for the adjustment when LFLAG is true
!
      IOPT=0
      IWORK(20+P99)=0.0d+00
      IF (ICALL /= 2) THEN
        IF (ICALL == 0) THEN    ! Initialization call
          ADVID=13
          ADVID2=ADVID

          INRD(0:MCOMP)=0
          IANRD(0:MCOMP)=0
          CALL TOLPROCESS(INRD,IANRD,PC-1)             ! Obtain the no. of significant digits
!          write(*,*) 'inrd ',inrd(0:mcomp)
!          write(*,*) 'ianrd ',ianrd(0:mcomp)
          IF (IERPRD >0) GO TO 999
          NRD=INRD(1)
          ZNRD(0)=ONE

          DO I=1,MCOMP
            IF (INRD(I) == 0) THEN
              ZNRD(I)=ZNRD(I-1)
            ELSE
              ZNRD(I)=TEN**(-INRD(I))
            ENDIF
          END DO
          ZANRD(0)=ZERO
          DO I=1,MCOMP
            IF (IANRD(I) == 0) THEN
              ZANRD(I)=ZANRD(I-1)
            ELSE
              ZANRD(I)=TEN**(-IANRD(I))
            ENDIF
          END DO

          TL=TEN**(-NRD)

          IF(TOLTOPRED==1) GO TO 999
          
          SPW=PPW
! Set maximum values     
          DADT(1:MCOMP)=ZERO
          FORALL(I=1:MCOMP,J=1:9) IATT(I,J)=0
          FORALL(I=1:MCOMP,J=1:PG) DP(I,J)=ZERO
          FORALL(I=1:MCOMP,J=1:MCOMP) DA(I,J)=ZERO
          FORALL(I=1:MCOMP,J=1:MCOMP+1) LINK(I,J)=0
! Call for MODEL definition information
          Y(1:P99)=ZERO
          CALL MODEL(IDNO,NCM1,NPR,MCOMP,IATT,LINK) 
          IF (IERPRD >0) GO TO 999
          WRITE (LOGUNT,125,ERR=502)
          NBP=NPR; IOUT=NCM1+1
!          METHOD=2
!          MITER=2

! No. sig. digits = 0 for compartment I; copy prior value for remainder   
! Set unchanging info. for LSODA
          ITOL=4
! 
          LRW=22+9*P99+P99*P99   ! LRW=11407; changed to give the true size of RWORK  3/2005 AJB
!       LIW = 119
          LIW=20+P99 ! Bauer changes
! Turn off debugging flag
          IPRINT=0    
!             
! Flags describing external events:
! MFLAG1: 1 when start of IR or reset has occurred
! MFLAG2: 1 when call to PK has occurred
! MFLAG3: 1 when a new bolus dose was added to state vector
! MFLAG4: 1 when an infusion was started
! MFLAG5: 1 when an infusion has ended (internal to ADVAN9)
! MFLAG6: 1 when status vector changed (compt ON)
! MFLAG7: 1 when status vector changed (compt OFF)
! MFLAG8: 1 when steady state is being computed        
     
! Set flags initially
          MFLAG1=0; MFLAG4=0; MFLAG7=0; OLDDT=0.
          MFLAG2=0; MFLAG5=0; MFLAG8=0
          MFLAG3=0; MFLAG6=0; MFLAGA=0

        ELSE    ! Added 8/94. allows initialization call to DES        
          CALL DES(AMNT,G3,TSTART,DADT,PIR,DAC,DPC,DTC)  
        END IF
        GO TO 999
      END IF  




! Normal Entry
      IF (DOFINL) THEN
        IF (COMPAC) THEN
          CALL DES(AMNT,G3,TSTART,DADT,PIR,DAC,DPC,DTC)
        ELSE
          CALL DES(AMNT,G3,TSTART,DADT,MCOMP,DA,DP,DET)
        END IF
        DIDDES=.FALSE.
        GO TO 999
      END IF
!
! Not a special call
      DIDCAA=.TRUE.
      DIDDES=.TRUE.
      IF(IFIRSTEM*IFRSTDER>0) THEN
        NPETAS=NPETAS2
        NOETAS=NOETAS2
        NY=NCM1*(NPETAS+1)
      ELSE
        NY=NCM1
        NPETAS=0
        NOETAS=.TRUE.
      ENDIF
!      
!      IF (MITER == -1) THEN
!        MIT=2
!      ELSE
!        MIT=MITER
!      END IF
!
      IF (IMAX == -1) THEN
        MAXCAL=MAXFCN
        IF (MAXCAL == 0) MAXCAL=1000000  ! Default added 11/2007 in case old SIZES is used
      ELSE
        MAXCAL=IMAX
      END IF
!
      IF (LFLAG .AND. (.NOT. NOETAS)) THEN
        DTINT(1:NPETAS)=DTSTAR(1:NPETAS)
      END IF
!
      I2FLAG=0   ! Set I2FLAG=1 if there will be partials w.r.t time due to term. infusion
      DTDER=LFLAG .OR. DTIME
!      
      IF (ISPEC /= 0 .AND. (.NOT. NOETAS)) THEN
        DO L=1,NPETAS
          IF (IDEA(ISPEC,L) /= ZERO) THEN 
            I2FLAG=1
            DTDER=.TRUE.
            EXIT
          END IF
        END DO
      END IF  
!
! HR(I,1)=infusion rate into A(I) after terminating inf. ends
! HR(I,2)=infusion rate into A(I) before it ends
!
      HR(1:NCM1,1)=R(1:NCM1)
      HR(1:NCM1,2)=R(1:NCM1)
!
      IF (.NOT. NOETAS) THEN
        FORALL(I=1:NCM1,L=1:NPETAS) HETA(I,L,1)=RE(I,L)
        FORALL(I=1:NCM1,L=1:NPETAS) HETA(I,L,2)=RE(I,L)
      END IF 
! 
      IF (ISPEC /= 0) THEN
        ICOMP=BETA(ISPEC)
        HR(ICOMP,2)=HR(ICOMP,2)+IRA(ISPEC)
        IF (.NOT. NOETAS) THEN
          DO L=1,NPETAS
             HETA(ICOMP,L,2)=HETA(ICOMP,L,2)+IREA(ISPEC,L)
          END DO   
        END IF  
      END IF  
!
! Save old values if needed for output compartment
      IF (SV(IOUT) /= 0) THEN
        DO I=1,NCM1
          OLDA(I)=AMNT(I)
          IF (NOETAS) CYCLE
          DO L=1,NPETAS
            OLDAET(I,L)=DAETA(I,L)
          END DO
        END DO
      END IF  
!
      IF (.NOT. NOETAS) THEN
        IF (SECOND) THEN
          IF (ISPEC /= 0 .AND. I2FLAG /= 1) THEN
            OUTER: DO L=1,NPETAS
              DO J=L,NPETAS
                IF (I2DEA(ISPEC,J,L) /= ZERO) THEN 
                  I2FLAG=1
                  DTDER=.TRUE.
                  EXIT OUTER
                END IF
              END DO
            END DO OUTER
          END IF  
!
          DO L=1,NPETAS
            DO J=L,NPETAS
              DO I=1,NCM1
                H2ETA(I,J,L,1)=R2E(I,J,L)
                H2ETA(I,J,L,2)=R2E(I,J,L)
              END DO
            END DO
          END DO
!
          IF (ISPEC /= 0) THEN
            DO L=1,NPETAS
              DO J=L,NPETAS
                H2ETA(ICOMP,J,L,2)=H2ETA(ICOMP,J,L,2)+I2REA(ISPEC,J,L)
              END DO
            END DO
          END IF
!
          IF (SV(IOUT) /= 0) THEN
            DO L=1,NPETAS
              DO J=L,NPETAS
                DO I=1,NCM1
                  OLDAE2(I,J,L)=D2AETA(I,J,L)
                END DO
              END DO
            END DO
          END IF
!
          IF (LFLAG) THEN
            DO L=1,NPETAS
              DO J=L,NPETAS
                D2TINT(J,L)=D2TSTA(J,L)
              END DO
            END DO
          END IF
        END IF  
!
! SS6 sets MFLAG8=1 when the call is from ZSPOW
        IF (MFLAG8 /= 1) THEN
          IF (.NOT. SECOND) THEN
            LST=1; LEND=NPETAS
          ELSE
            KLOOP=SLOOP
            LEND=SLEND(KLOOP)
            LST=SLST(KLOOP)
          END IF  
        END IF  
      END IF  


! Set up now for first integration interval
        MFLAGX=MFLAG1+MFLAG6+MFLAG7   ! MFLAGX is change of compartment status    
        MFLAGY=MFLAG2+MFLAG3          ! MFLAGY is call to PK or new bolus dose
        MFLAGZ=MFLAG4+MFLAG5          ! MFLAGZ is change in infusion

        IF (MFLAG2 /= 0) THEN
          MF=20
          IF (MITER == -1) THEN
            IF (NOETAS) THEN
              MF=MF+1
            ELSE
              MF=MF+2
            END IF
          ELSE
            MF=MF+MITER
          END IF
        END IF
! Get relative tolerances for each compartment and its partials
        IF (MFLAGX /= 0) THEN
          K=0
          IF (MAPPED) THEN
            DO J=0,NPETAS
              DO I=1,NCM1
                K=K+1
                RTOL(K)=ZNRD(MAPINV(I))
                AOTOL(K)=ZANRD(MAPINV(I))
              END DO
            END DO

            IF (SECOND) THEN
              DO L=1,NPETAS
                DO J=1,L
                  DO I=1,NCM1
                    K=K+1
                    RTOL(K)=ZNRD(MAPINV(I))
                    AOTOL(K)=ZANRD(MAPINV(I))
                  END DO
                END DO
              END DO
            END IF

          ELSE
            DO J=0,NPETAS
              DO I=1,NCM1
                K=K+1
                RTOL(K)=ZNRD(I)
                AOTOL(K)=ZANRD(I)
              END DO
            END DO

            IF (SECOND) THEN
              DO L=1,NPETAS
                DO J=1,L
                  DO I=1,NCM1
                    K=K+1
                    RTOL(K)=ZNRD(I)
                    AOTOL(K)=ZANRD(I)
                  END DO
                END DO
              END DO
            END IF

          END IF
        END IF
!
! Set up now for first integration interval
 1110 Y(1:NCM1)=AMNT(1:NCM1)
!
      IF (.NOT. NOETAS) THEN
        K=0
        DO L=1,LEND
          K=K+NCM1
          DO I=1,NCM1
            Y(K+I)=DAETA(I,L)
          END DO
        END DO
        IF (SECOND) THEN
          DO L=LST,LEND
            DO J=1,L
              K=K+NCM1
              DO I=1,NCM1
                Y(K+I)=D2AETA(I,L,J)
              END DO
            END DO
          END DO
          NY=K+NCM1
        END IF
      END IF  
!
      IH=2
      X=TSTART
      XEND=X+DT
      IF (ISPEC /= 0) XEND=X+IDA(ISPEC)
      IND=1
      H=(XEND-X)*P00001
!
! Only FCN is called when MIT=2, so pass same name for FCN and FCNJ
      IF (H > 0) THEN
        CALCTR=0
!        CALL DGEAR1(NY,FCN1,FCN3,X,H,Y,XEND,TL,METH,MIT,IND,IWK,WK,IER)
        CALL ATOL_SET(MF,ATOL,RTOL,AOTOL,NY)
        ISTATE=1
        ITASK=1
        IREDO=1
  340   CONTINUE        
        CALL LSODA(RES2,NY,Y,X,XEND,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,LIW,JAC,MF)
        IF(ISTATE==-1.AND.IREDO==0) THEN
          ISTATE=2
          ITASK=1
          IREDO=1
          GO TO 340
        ENDIF
        IF (IERPRD > 0) THEN
          IF (IERPRD == 999) CALL RUN_TIME_ERRORS(2)
          GO TO 999
        END IF
          IF (ISTATE < 0 .OR. ISTATE > 2) THEN 
            WRITE (IST,'(I4)') ISTATE
            ETEXT(2)='ERROR IN LSODA: CODE'//IST
            IERPRD=1; GO TO 999
          END IF
      END IF
!
! All finished if no terminating infusion
      IF (ISPEC /= 0)  THEN
! Branch if numeric derivatives are not to be obtained
        IF (I2FLAG /= 0) THEN
! Adjustment for end point - FCN1 handles on/off, adds both infusions
          SAVSEC=SECOND
          SECOND=.FALSE.
          CALCTR=0
!          CALL FCN1(Y,DADT,NY,XEND)
          CALL RES2(NY,XEND,Y,DADT)
          IF (IERPRD > 0) THEN
            IF (IERPRD == 999) CALL RUN_TIME_ERRORS(2)
            GO TO 999
          END IF  
!          
          FORALL(L=1:NPETAS) DTINT(L)=DTINT(L)+IDEA(ISPEC,L)
          K=0
          DO J=1,LEND
            K=K+NCM1
            DO I=1,NCM1
              Y(K+I)=Y(K+I)+DADT(I)*IDEA(ISPEC,J)
            END DO
          END DO
!
          IF (SAVSEC) THEN
! Second set of derivatives
            DO L=LST,LEND
              DO J=1,L
                D2TINT(L,J)=D2TINT(L,J)+I2DEA(ISPEC,L,J)
              END DO
            END DO
            CALCTR=0
!            CALL FCN1(Y,DADT2,NY,XEND)
            CALL RES2(NY,XEND,Y,DADT2)
            IF (IERPRD > 0) THEN
              IF (IERPRD == 999) CALL RUN_TIME_ERRORS(2)
              GO TO 999
            END IF  
!            
            DO L=LST,LEND
              DO J=1,L
                K=K+NCM1
                DO I=1,NCM1
                  Y(K+I)=Y(K+I)+DADT2(J*NCM1+I)*IDEA(ISPEC,L)+DADT(I)*I2DEA(ISPEC,L,J) &
                         +DADT(L*NCM1+I)*IDEA(ISPEC,J)
                END DO
              END DO
            END DO
          END IF  
          SECOND=SAVSEC
        END IF  
!
! Start of second integration interval
        X=TSTART+IDA(ISPEC)
        XEND=TSTART+DT
        IND=1
        H=(XEND-X)*P00001
        IH=1
        IF (H > 0) THEN
          CALCTR=0
!          CALL DGEAR1(NY,FCN1,FCN3,X,H,Y,XEND,TL,METH,MIT,IND,IWK,WK,IER)
        ISTATE=1
        ITASK=1
        CALL ATOL_SET(MF,ATOL,RTOL,AOTOL,NY)
        IREDO=1
 341    CONTINUE
        CALL LSODA(RES2,NY,Y,X,XEND,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,LIW,JAC,MF)
        IF(ISTATE==-1.AND.IREDO==0) THEN
          ISTATE=2
          ITASK=1
          IREDO=1
          GO TO 341
        ENDIF
          IF (IERPRD > 0) THEN          
            IF (IERPRD == 999) CALL RUN_TIME_ERRORS(2)
            GO TO 999
          END IF  
          IF (ISTATE < 0 .OR. ISTATE > 2) THEN 
            WRITE (IST,'(I4)') ISTATE
            ETEXT(2)='ERROR IN LSODA: CODE'//IST
            IERPRD=1; GO TO 999
          END IF
        END IF
      END IF  
!
! Finished integrating to TSTART+DT
      IF (.NOT. NOETAS) THEN
        K=NCM1*(LST-1)
        DO L=LST,LEND
          K=K+NCM1
          DO I=1,NCM1
            DAETA(I,L)=Y(K+I)
          END DO
        END DO
        IF (SECOND) THEN
          DO L=LST,LEND
            DO J=1,L
              K=K+NCM1
              DO I=1,NCM1
                D2AETA(I,L,J)=Y(K+I)
              END DO
            END DO
          END DO
!
! Go back for another group of ETAS if necessary
! Must also subtract IDEA from DTINT before re-doing the integ.
          IF (MFLAG8 == 0) THEN
            KLOOP=KLOOP-1
            IF (KLOOP >= 1) THEN
              IF (ISPEC /= 0) THEN
                FORALL(L=1:NPETAS) DTINT(L)=DTINT(L)-IDEA(ISPEC,L)
                DO L=LST,LEND
                  DO J=1,L
                    D2TINT(L,J)=D2TINT(L,J)-I2DEA(ISPEC,L,J)
                  END DO
                END DO
              END IF
              LEND=SLEND(KLOOP)
              LST=SLST(KLOOP)
              GO TO 1110
            END IF
          ELSE
! For call from SS, must return all the first derivatives
            IF (LST > 1) THEN
              K=0
              DO L=1,LST-1
                K=K+NCM1
                DO I=1,NCM1
                  DAETA(I,L)=Y(K+I)
                END DO
              END DO
            END IF
          END IF
        END IF  
      END IF  
!
      AMNT(1:NCM1)=Y(1:NCM1)
!
      IF (I2FLAG == 1 .OR. DTIME) THEN
! Adjust derivatives at end point - FCN1 handles on/off, adds infusion
        IF (MFLAG8 == 0) THEN
          LST=1; LEND=NPETAS
        END IF
        SAVSEC=SECOND
        SECOND=.FALSE.
!
! Any value greater than NCM1 will do
        NY=2*NCM1
        Y(1:NCM1)=AMNT(1:NCM1)
        K=0
        DO L=1,NPETAS
          K=K+NCM1
          DO I=1,NCM1
            Y(K+I)=DAETA(I,L)
          END DO
        END DO
!        
        CALCTR=0
!        CALL FCN1(Y,DADT,NY,XEND)
        CALL RES2(NY,XEND,Y,DADT)
        IF (IERPRD > 0) THEN
          IF (IERPRD == 999) CALL RUN_TIME_ERRORS(2)
          GO TO 999
        END IF  
!        
        DO L=1,NPETAS
          IF (ISPEC > 0) THEN
            DO I=1,NCM1
              DAETA(I,L)=DAETA(I,L)+DADT(I)*(ADTE(L)-IDEA(ISPEC,L))
            END DO
          ELSE
            DO I=1,NCM1
              DAETA(I,L)=DAETA(I,L)+DADT(I)*ADTE(L)
            END DO
          END IF
        END DO
!
        IF (SAVSEC) THEN
          K=0
          DO L=1,NPETAS
            K=K+NCM1
            DO I=1,NCM1
              Y(K+I)=DAETA(I,L)
            END DO
          END DO
          IF (ISPEC /= 0) THEN
            FORALL(L=1:NPETAS) DTINT(L)=DTINT(L)-IDEA(ISPEC,L)+ADTE(L)
          ELSE
            FORALL(L=1:NPETAS) DTINT(L)=DTINT(L)+ADTE(L)
          END IF
          CALCTR=0
!          CALL FCN1(Y,DADT2,NY,XEND)
          CALL RES2(NY,XEND,Y,DADT2)
          IF (IERPRD > 0) THEN
            IF (IERPRD == 999) CALL RUN_TIME_ERRORS(2)
            GO TO 999
          END IF  
!          
          DO L=LST,LEND
            DO J=1,L
              DO I=1,NCM1
                IF (ISPEC > 0) THEN
                  D2AETA(I,L,J)= D2AETA(I,L,J)+DADT2(J*NCM1+I)*(ADTE(L)-IDEA(ISPEC,L))  &
                                +DADT(I)*(D2ADTE(L,J)-I2DEA(ISPEC,L,J))                 &
                                +DADT(L*NCM1+I)*(ADTE(J)-IDEA(ISPEC,J))
                ELSE
                  D2AETA(I,L,J)= D2AETA(I,L,J)+DADT2(J*NCM1+I)*ADTE(L)                  &
                                +DADT(I)*D2ADTE(L,J)+DADT(L*NCM1+I)*ADTE(J)
                END IF
              END DO
            END DO
          END DO
        END IF  
        SECOND=SAVSEC
      END IF
!
! Calculate the output compartment
      IF (SV(IOUT) /= 0) THEN
        DO I=1,NCM1
          IF(JTOTL(I)/=0) CYCLE
          AMNT(IOUT)=AMNT(IOUT)+OLDA(I)-AMNT(I)+R(I)*DT
        END DO
!        
        IF (.NOT. NOETAS) THEN
          DO L=1,NPETAS
            DO I=1,NCM1
              IF(JTOTL(I)/=0) CYCLE
              DAETA(IOUT,L)=DAETA(IOUT,L)+OLDAET(I,L)-DAETA(I,L)+RE(I,L)*DT+R(I)*ADTE(L)
            END DO
          END DO
          IF (SECOND) THEN
            DO L=1,NPETAS
              DO J=L,NPETAS
                DO I=1,NCM1
                  IF(JTOTL(I)/=0) CYCLE
                  D2AETA(IOUT,J,L)=D2AETA(IOUT,J,L)+OLDAE2(I,J,L)-D2AETA(I,J,L)             &
                                   +R2E(I,J,L)*DT+RE(I,J)*ADTE(L)+RE(I,L)*ADTE(J)           &
                                   +R(I)*D2ADTE(J,L)
                END DO
              END DO
            END DO
          END IF
        END IF  
!        
        IF (ISPEC /= 0) THEN
          AMNT(IOUT)=AMNT(IOUT)+IRA(ISPEC)*IDA(ISPEC)
          IF (.NOT. NOETAS) THEN 
            DO L=1,NPETAS
              DAETA(IOUT,L)= DAETA(IOUT,L)+IRA(ISPEC)*IDEA(ISPEC,L)+IDA(ISPEC)*IREA(ISPEC,L)
            END DO  
            IF (SECOND) THEN 
              DO L=1,NPETAS
                DO J=L,NPETAS
                  D2AETA(IOUT,J,L)= D2AETA(IOUT,J,L)+I2DEA(ISPEC,J,L)*IRA(ISPEC)            &
                                   +IDEA(ISPEC,J)*IREA(ISPEC,L)+IDEA(ISPEC,L)*IREA(ISPEC,J) &
                                   +IDA(ISPEC)*I2REA(ISPEC,J,L)
                END DO
              END DO
            END IF  
          END IF  
        END IF  
      END IF  
!
! Exit point 
      IF (LFLAG .OR. I2FLAG /= 0 .OR. DTIME) THEN
        DTINT(1:NPETAS)=ZERO
        IF (SECOND) THEN
          DO L=1,NPETAS
            DO J=L,NPETAS
              D2TINT(J,L)=ZERO
            END DO
          END DO
        END IF
      END IF

! Reset the flags before return
! MFLAG5 is reset internally
! MFLAG8 is reset by SS9
      MFLAG1=0; MFLAG4=0; OLDDT=DT
      MFLAG2=0; MFLAG6=0
      MFLAG3=0; MFLAG7=0

!
  125 FORMAT (' GENERAL NONLINEAR KINETICS MODEL WITH STIFF/NONSTIFF EQUATIONS (LSODA, ADVAN13)')
!
  999 RETURN
!
  502 CALL ERRORMSGS(502,FILENAME='OUTPUT')
!
      CONTAINS
!
! RUN-TIME error conditions; Error in integration to IDA(ISPEC) or DT

        SUBROUTINE RUN_TIME_ERRORS(ITEMP)
!
          INTEGER, INTENT(IN) :: ITEMP
!
          SELECT CASE(ITEMP)
          CASE(1)
            ETEXT(1)='NUMERICAL DIFFICULTIES WITH INTEGRATION ROUTINE.'
            WRITE (CTEMP,'(I3)') NRD
            ETEXT(2)='NO. OF REQUIRED SIGNIFICANT DIGITS IN SOLUTION VECTOR'
            ETEXT(3)='TO DIFFERENTIAL EQUATIONS, '//CTEMP// ', MAY BE TOO LARGE.'
          CASE(2)  
            ETEXT(2)='NUMERICAL DIFFICULTIES WITH INTEGRATION ROUTINE.'
            WRITE (CTEMPC,'(I12)') MAXCAL
            ETEXT(3)='MAXIMUM NO. OF EVALUATIONS OF DIFFERENTIAL EQUATIONS, '            &
            //CTEMPC//', EXCEEDED.'
          END SELECT 
!
          IERPRD=1
!
  999     RETURN   
!
        END SUBROUTINE RUN_TIME_ERRORS
!
      END SUBROUTINE ADVAN

!-----------------------------HISTORY------------------------------------------------
! VERSION     : NONMEM VII
! AUTHOR      : ROBERT J. BAUER
! CREATED ON  : AUG/2016
! LANGUAGE    : FORTRAN 90/95
! LAST UPDATE : NOV/2016 - INTRODUCED HEADER INFORMATIONS AND RESTRUCTURED AS PER
!                          THE NONMEM STANDARDS
!
!------------------------------------ JTOTL.F90 -------------------------------------
!
! FUNCTION JTOTL(I)
!
! DESCRIPTION : special mapping routine
!
!
! ARGUMENTS   : I
!               IN     - I
!               OUT    - NONE
!               IN OUT - NONE
!
! CALLED BY   : ADVAN,ADVAN
!
! CALLS       : NONE
!
! ALGORITHM   : see description
!
! MODULES USED: SIZES,PRCM_INT,PRCM_LOG
!
! CONTAINS    : NONE
!
! LOCALS      : JTOTL
!
!---------------------------- END OF HEADER -----------------------------------------
!
      FUNCTION JTOTL(I)
!
      USE PRSIZES,        ONLY: ISIZE
!
! INTEGER
      USE PRCM_INT,     ONLY: ITOTL,MAPINV
! LOGICAL
      USE PRCM_LOG,     ONLY: MAPPED
!
      IMPLICIT NONE
!
      INTEGER(KIND=ISIZE), INTENT(IN) :: I
!
!
!------------------------------------------------------------------------------------
!
! Local Variables
!
      INTEGER(KIND=ISIZE) :: JTOTL
!
      IF(MAPPED) THEN
        JTOTL=ITOTL(MAPINV(I))
      ELSE
        JTOTL=ITOTL(I)
      ENDIF
 999  CONTINUE
      RETURN
      END FUNCTION JTOTL
   
!    
!-----------------------------HISTORY------------------------------------------------
! VERSION     : NONMEM VII
! AUTHOR      : ROBERT J. BAUER
! CREATED ON  : DEC/2008
! LANGUAGE    : FORTRAN 90/95
! LAST UPDATE : DEC/2008 - INTRODUCED HEADER INFORMATIONS AND RESTRUCTURED AS PER
!                          THE NONMEM STANDARDS
!               AUG/2012 - INTEGRATED NONMEM7.3ALPHA6.3 CHANGES
! LAST UPDATE : NOV/2016 - INTRODUCED HEADER INFORMATIONS AND RESTRUCTURED AS PER
!                          THE NONMEM STANDARDS.  Used ADVAN9 as template
!
!---------------------------- ATOL_SET.F90 ------------------------------------------
!
! SUBROUTINE ATOL_SET(MF,ATOL,RTOL,AOTOL,NY)
!
! DESCRIPTION : Sets the Absolute Tolerance parameter based on column dimension. 
!
! ARGUMENTS   : MF,ATOL,RTOL,AOTOL,NY
!               IN     - NY,RTOL,AOTOL
!                        NY   - Column dimension.
!                        RTOL - Relative Tolerance
!                        ATOL - original Absolute Tolerance
!               OUT    - NONE
!               IN OUT - MF,ATOL
!                        ATOL - Absolute Tolerance Parameter (scalar or vector)
!
! CALLED BY   : ADVAN13 - ADVAN13 is a routine in PREDPP's library which implements the general
!                         non-linear model with equilibrium compartments. This general non-
!                         linear model is used for systems in which a drug is distributed
!                         between compartments according to a system of first-order
!                         differential-algebraic processes.
!               ADVAN15  - ADVAN15 is a routine in PREDPP's library which implements  the  general
!                         non-linear  model  with  equilibrium  compartments.  This general non-
!                         linear model is used for  systems  in  which  a  drug  is  distributed
!                         between   compartments   according   to   a   system   of  first-order
!                         differential-algebraic processes.
!
! CALLS       : DES
!
! ALGORITHM   : - Loop over the column dimension and set the Absolute Tolerance parameter.
!
! MODULES USED: PRSIZES,NMBAYES_INT,NMPRD_INT
!
! CONTAINS    : NONE
!
! LOCAL'S     : I
!
!---------------------------- END OF HEADER -----------------------------------------
!
      SUBROUTINE ATOL_SET(MF,ATOL,RTOL,AOTOL,NY)     ! Bauer added this subroutine
!
      USE PRSIZES,     ONLY: ISIZE,DPSIZE,PIR
! INTEGER
      USE NMBAYES_INT, ONLY: ATOLU=>ATOL    ! 7.2b55a
      USE NMPRD_INT,    ONLY: IFIRSTEM,IFRSTDER,IFIRSTEMJAC
!
      IMPLICIT NONE
!
      INTEGER(KIND=ISIZE), INTENT(IN)     :: NY
!
      REAL(KIND=DPSIZE),   INTENT(IN)     :: RTOL(*),AOTOL(*)
      REAL(KIND=DPSIZE),   INTENT(IN OUT) :: ATOL(*)
      INTEGER(KIND=ISIZE), INTENT(IN OUT)     :: MF
!
!------------------------------------------------------------------------------------
!
! Local Variables
!
      INTEGER(KIND=ISIZE) :: I
      REAL(KIND=DPSIZE) :: DVAL1,DVAL2,DVAL3,DVAL4,DVAL5,DVAL6,DVAL7
!
! TEST THAT IF NO JACOBIAN AVAILABLE, THEN HAVE LSODI1 DETERMINE IT
      IFIRSTEMJAC=-2
! CALLING DES WITH IFIRSTEMJAC RESULTS IN IFIRSTEMJAC BEING SET, AND IMMEDIATE RETURN.
! SO DUMMY ARGUMENTS ARE SUFFICIENT.
      CALL DES(DVAL1,DVAL2,DVAL3,DVAL4,PIR,DVAL5,DVAL6,DVAL7)
      IF(IFIRSTEMJAC==0) MF=22      
      ATOL(1:NY)=AOTOL(1:NY)
!
  999 RETURN
!
      END SUBROUTINE ATOL_SET

!*********************************COPYRIGHT******************************************
!                                                                                   !
!       THE NONMEM SYSTEM MAY BE DISTRIBUTED ONLY BY ICON DEVELOPMENT               !
!       SOLUTIONS.                                                                  !
!                                                                                   !
!       COPYRIGHT BY ICON DEVELOPMENT SOLUTIONS                                     !
!       2009-2017 ALL RIGHTS RESERVED.                                              !
!                                                                                   !
!************************************************************************************
!
!-----------------------------HISTORY------------------------------------------------
! VERSION     : NONMEM VII
! AUTHOR      : ROBERT J. BAUER
! CREATED ON  : AUG/2016
! LANGUAGE    : FORTRAN 90/95
! LAST UPDATE : NOV/2016 - INTRODUCED HEADER INFORMATIONS AND RESTRUCTURED AS PER
!                          THE NONMEM STANDARDS
!
!------------------------------------ RES2.F90 --------------------------------------
!
! SUBROUTINE RES2(NEQ,T,X,XP)
!
! DESCRIPTION : connector between user ODE functions and LSODA
!
!
! ARGUMENTS   : NEQ,T,X,XP
!               IN     - NEQ
!               OUT    - NONE
!               IN OUT - T,X,XP
!
! CALLED BY   : ADVAN
!
! CALLS       : FCN1
!
! ALGORITHM   : see description
!
! MODULES USED: PRSIZES,SIZES
!
! CONTAINS    : NONE
!
! LOCALS      : NONE
!
!---------------------------- END OF HEADER -----------------------------------------
!
      SUBROUTINE RES2(NEQ,T,X,XP)
!
      USE PRSIZES,      ONLY: ISIZE,DPSIZE,PC,PE
!
! INTERFACE
      USE PR_INTERFACE, ONLY: RES,FCN1
!
      IMPLICIT NONE
!
      INTEGER(KIND=ISIZE), INTENT(IN)     :: NEQ
      REAL(KIND=DPSIZE),   INTENT(IN OUT)     :: T,X(*),XP(*)
!
!
!------------------------------------------------------------------------------------
!
! Local Variables
!
!
      CALL FCN1(X,XP,NEQ,T)
  999 RETURN
!
      END SUBROUTINE RES2
