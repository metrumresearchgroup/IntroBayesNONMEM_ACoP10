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
! CREATED ON  : AUG/2016
! LANGUAGE    : FORTRAN 90/95
! LAST UPDATE : NOV/2016 - INTRODUCED HEADER INFORMATIONS AND RESTRUCTURED AS PER
!                          THE NONMEM STANDARDS
!
!----------------------------- JAC.F90 ---------------------------------------------
! SUBROUTINE JAC(NY,PDD,NROWP)
!
! DESCRIPTION : Evaluate the Jacobian matrix of partial derivatives.  
!               Based on FCN3, except that NROWP was added
!
! ARGUMENTS   : NY,PDD,NROWP
!               IN     - NY,NROWP
!                        NY     - Size of Jabobian matrix
!                        NROWP  - Row dimension of PDD
!               OUT    - PDD
!                        PDD - Modified the Jacobian
!               IN OUT - NONE
!
! CALLED BY   : ADVAN13, ADVAN14
!
! CALLS       : NONE
!
! ALGORITHM   : - Copy DA as Jacobian (only valid for indiv. data) - for MITER=1
!                 - Check for compact arrays
!                 - Zero-out PDD because DGEAR passes a portion of some work array
!                 - Check if mapped
!
! MODULES USED: PRSIZES,PRCM_INT,PRCM_REAL,PRCOM_REAL,PRCM_LOG 
!
! LOCAL'S     : I,J
!
!---------------------------- END OF HEADER -----------------------------------------
!
      SUBROUTINE JAC(NY,PDD,NROWP)     
!      
      USE PRSIZES,    ONLY: ISIZE,DPSIZE
! INTEGER    
      USE PRCM_INT,   ONLY: AA,IAI,IAJ,MAPINV,NIAI1
! REAL
      USE PRCM_REAL,  ONLY: DAC
      USE PRCOM_REAL, ONLY: DA
! LOGICAL
      USE PRCM_LOG,   ONLY: MAPPED,COMPAC
!
      IMPLICIT NONE
!     
      INTEGER(KIND=ISIZE), INTENT(IN)  :: NY,NROWP
      REAL(KIND=DPSIZE),   INTENT(OUT) :: PDD(NROWP,*)
!      
      SAVE
!      
!------------------------------------------------------------------------
!      COMMON /PRCOMD/ DA,DP,DET,HR,HETA,H2ETA
!      DOUBLE PRECISION PDD !,DA,DP,HR,HETA,DET,H2ETA
!      DIMENSION PDD(NY,NY) !,DA(PM,PM),DP(PM,PG),HR(PM,2)
!      COMMON /PRCMX1/ GENMOD,MAPPED,COMPAC
!      LOGICAL GENMOD,MAPPED,COMPAC
!      COMMON /PRCMX2/ XNC,XNCM1,MAP,MAPINV,NBRON,XNBRON
!      INTEGER XNC,XNCM1,MAP(PC),MAPINV(PC),NBRON,XNBRON
!      DIMENSION HETA(PM,PE,2),H2ETA(PM,PE,PE,2),DET(PM)
!      COMMON /PRCMDE/ NIAI,NIPI,NIT,XNIAI,XNIPI,XNIT
!      INTEGER NIAI,NIPI,NIT,XNIAI,XNIPI,XNIT
!      COMMON /PRCMDE/ IAI(PIR),IAJ(PIR),IAK(PIR)
!      COMMON /PRCMDE/ IPI(PIR),IPJ(PIR),IPK(PIR)
!      COMMON /PRCMDE/ ITI(PIR),ITK(PIR)
!      COMMON /PRCMDE/ XIAI(PIR),XIAJ(PIR),XIAK(PIR)
!      COMMON /PRCMDE/ XIPI(PIR),XIPJ(PIR),XIPK(PIR)
!      COMMON /PRCMDE/ XITI(PIR),XITK(PIR)
!      COMMON /PRCMDE/ IAC(PIR),IPC(PIR),ITC(PIR)
!      COMMON /PRCMDE/ NIAI1,NIPI1,NIT1,XNIAI1,XNIPI1,XNIT1
!      COMMON /PRCMDE/ AA,PP,TT,XAA,XPP,XTT
!      INTEGER NIAI1,NIPI1,NIT1,XNIAI1,XNIPI1,XNIT1
!      INTEGER AA(PIR),PP(PIR),TT(PIR),XAA(PIR),XPP(PIR),XTT(PIR)
!      INTEGER IAI,IAJ,IAK,IPI,IPJ,IPK,ITI,ITK
!      INTEGER XIAI,XIAJ,XIAK,XIPI,XIPJ,XIPK,XITI,XITK
!      INTEGER IAC,IPC,ITC
!      COMMON /PRCMDC/ DAC(PIR,1),DPC(PIR,1),DTC(PIR)
!      DOUBLE PRECISION DAC,DPC,DTC
!------------------------------------------------------------------------
!
! Local Variables
!
      INTEGER(KIND=ISIZE) :: I,J
!
! For implementation of general ADVANS' compartment on/off feature
! GENMOD true for general ADVAN (Set by PREDI)
! XNC,XNCM1,XNBRON are original values of NC,NCM1,NBRON (Set by PREDI)
! NBRON tells how many compartments are now on (Not counting output)
! Mapped is true when a mapping in effect (Set by SADVAN, SSS)
! MAP,MAPINV: Maps 'REAL' compartment nos. to 'Reduced set' and V.V.
! Counters for indices for the COMPACT arrays (Versions are unmapped)
! Indices for the COMPACT arrays (Versions are unmapped)
! The COMPACT arrays themselves
!
! Copy DA as Jacobian (only valid for Indiv. data) - for MITER=1
!
    
      IF (COMPAC) THEN          ! Compact arrays.
        PDD(1:NY,1:NY)=0.0D0    ! Zero-out PDD because DGEAR passes a portion of some work array
        FORALL (I=1:NIAI1) PDD(IAI(I),IAJ(I))=DAC(AA(I),1) 
      ELSE 
        IF (MAPPED) THEN   
          FORALL(I=1:NY,J=1:NY) PDD(I,J)=DA(MAPINV(I),MAPINV(J))
        ELSE
          FORALL(I=1:NY,J=1:NY) PDD(I,J)=DA(I,J)
        END IF
      END IF
!      
  999 RETURN
!    
      END SUBROUTINE JAC
