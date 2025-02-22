! This file has been modified for the use in ICON

MODULE YOERRTA13

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

PUBLIC

SAVE

!     -----------------------------------------------------------------
!*    ** *YOERRTA13* - RRTM COEFFICIENTS FOR INTERVAL 13
!     BAND 13:  2080-2250 cm-1 (low - H2O,N2O; high - nothing)
!     ABozzo 201306 updated to rrtmg v4.85
!     -----------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: NG13 = 4

REAL(KIND=JPRB) :: FRACREFA(NG13,9)
REAL(KIND=JPRB) , DIMENSION(NG13) :: FRACREFB

REAL(KIND=JPRB) :: KA(9,5,13,NG13) ,ABSA(585,NG13)
REAL(KIND=JPRB) :: SELFREF(10,NG13)
REAL(KIND=JPRB) :: FORREF(4,NG13)
REAL(KIND=JPRB) :: KA_MCO2(9,19,NG13)
REAL(KIND=JPRB) :: KA_MCO(9,19,NG13)
REAL(KIND=JPRB) :: KB_MO3(19,NG13)

EQUIVALENCE (KA(1,1,1,1),ABSA(1,1))

!$ACC DECLARE CREATE(FRACREFA, FRACREFB, KA, ABSA, SELFREF, FORREF, KA_MCO2, &
!$ACC                KA_MCO, KB_MO3)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! ABSA    : REAL     absorption coefficient of major absorber for M reference tropospheric 
!                    pressures and N reference tropospheric temperatures 
! FRACREFA: REAL     distance from r and T reference tabulated points (troposphere)
! KA      : REAL     absorption coefficient of major absorber (equiv. to ABSA)   
! SELFREF : REAL     self broadening coefficient for water vapour
! STRRAT  : REAL     weighting factors for the transition between tropospheric 
!                    and stratospheric computations
!     -----------------------------------------------------------------
END MODULE YOERRTA13
