!*************************************************************************
!*                                                                       *
!*     ISIC - Stochastic Nonsmooth Optimization based                    *
!*     Incremental Clustering Software                                   *
!*                                                                       *
!*                                                                       *
!*     Computation of the value of the clustering problem and the        *
!*     corresponding subgradients. This file is specially designed       *
!*     to solve clustering problems with stochastic inexact LMBM         *
!*     used within ISIC software.                                        * 
!*                                                                       *
!*     This work was part of the PRIVIDEMA project funded by the         *
!*     European Union under Grant Agreement No. 101167964.               *
!*                                                                       *
!*     ISIC software is covered by the MIT license.                      *
!*                                                                       *
!*************************************************************************
!* 
!*     Modules included:
!*
!*     obj_fun             ! Computation of the clustering function
!*                         ! and its subgradient.
!*

MODULE obj_fun             ! Computation of the value and the subgradient of the 
                           ! objective function.

    USE r_precision, ONLY : prec  ! Precision for reals.
    IMPLICIT NONE

    PUBLIC :: &
        myfg,          &   ! S Computation of the value and subgradient of the objective (combined).
        clustfull2         ! S Computation of value of clustering function with full batch and whole data

CONTAINS
    
    !==============================================================
    ! Combined computation of objective value and subgradient
    !==============================================================

    SUBROUTINE myfg(n,x,f,g,iterm)

        USE param, ONLY : zero, one, two, large
        USE initclust, ONLY : &
            a => a_used, & ! Data matrix.
            nft, &         ! Number of features.
            nc, &          ! Current number of clusters.
            tnorm, &       ! Total number of norms computed.
            m              ! Number of variables in optimization.
        USE initislmb, ONLY : &
            nbatch, &      ! Size of the batch.
            ib, &          ! First index of the batch.
            batchind       ! Batch indices.

        IMPLICIT NONE

        REAL(KIND=prec), DIMENSION(n), INTENT(IN)  :: x      ! Vector of variables.
        REAL(KIND=prec), INTENT(OUT)               :: f      ! Value of the objective.
        REAL(KIND=prec), DIMENSION(n), INTENT(OUT) :: g      ! Subgradient.
        INTEGER, INTENT(IN)  :: n                            ! Number of variables.
        INTEGER, INTENT(OUT) :: iterm                        ! 0=ok, -3=failure.

        REAL(KIND=prec) :: dcent1(nc,nc), f10, f2, f3
        INTEGER :: i, i2, j, k, kmin
        INTRINSIC :: MIN

        iterm = 0
        f = zero
        DO i = 1, m
            g(i) = zero
        END DO

            DO i = 1, nc
                dcent1(i,i) = zero
            END DO
            DO i = 1, nc
                DO j = i+1, nc
                    dcent1(i,j) = zero
                    DO k = 1, nft
                        dcent1(i,j) = dcent1(i,j) + (x(k+(i-1)*nft) - x(k+(j-1)*nft))**2
                    END DO
                    dcent1(j,i) = dcent1(i,j)
                END DO
            END DO

            kmin = 1
            DO i2 = 1, nbatch
                i   = batchind(i2+ib)
                f2  = large
                f10 = large
                loop1: DO k = 1, nc
                    IF (k >= 2) THEN
                        IF (dcent1(k,kmin) >= f10) CYCLE loop1
                    END IF
                    tnorm = tnorm + one
                    f3 = zero
                    DO j = 1, nft
                        f3 = f3 + (a(j,i) - x(j+(k-1)*nft))**2
                    END DO
                    IF (f3 < f2) THEN
                        f2   = f3
                        f10  = 4.0_prec * f2
                        kmin = k
                    END IF
                END DO loop1
                f = f + f2
                DO j = 1, nft      ! gradient for nearest center kmin
                    g(j+(kmin-1)*nft) = g(j+(kmin-1)*nft) + two*(x(j+(kmin-1)*nft) - a(j,i))
                END DO
            END DO

        IF (f > large)  iterm = -3
        IF (f < zero) iterm = -3

    END SUBROUTINE myfg

    !=================================================================
    ! Computation of clustering problem with full batch and full data
    !=================================================================

    SUBROUTINE clustfull2(x,f)

        USE param, ONLY : zero, one, two, large
        USE initclust, ONLY : &
            maxdim, &      ! Maximum number of variables in optimization,
                           !   maxdim >= nclust * nft.
            a, &           ! Data matrix.
            nft, &         ! Number of features.
            nc, &          ! Current number of clusters, loops from 1 to nclust.
            nrecord, &     ! Number of records in data.
            tnorm          ! The total number of norms computed.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec), DIMENSION(maxdim), INTENT(IN) :: &
            x              ! Vector of variables.
        REAL(KIND=prec), INTENT(OUT) :: &
            f              ! Value of the clustering function.
        
        ! Local variables 
        REAL(KIND=prec) ::  &
            dcent1(nc,nc), &
            f10,f2,f3

        INTEGER :: &
            i,j,k,kmin

        DO i=1,nc
            dcent1(i,i)=zero
        END DO

        DO i=1,nc
            DO j=i+1,nc
                dcent1(i,j)=zero
                 DO k=1,nft
                    dcent1(i,j)=dcent1(i,j)+(x(k+(i-1)*nft)-x(k+(j-1)*nft))**2
                END DO
                dcent1(j,i)=dcent1(i,j)
            END DO
        END DO        

        f=zero
        kmin=1
        DO i=1,nrecord
    
            f2=large
            f10=large
            loop1: DO k=1,nc
                IF(k >= 2) THEN
                    IF(dcent1(k,kmin) >= f10) THEN
                        CYCLE loop1
                    END IF
                END IF                     
                tnorm=tnorm+one
                f3=zero
                DO j=1,nft
                    f3=f3+(a(j,i)-x(j+(k-1)*nft))**2
                END DO
                IF(f3 < f2) THEN
                    f2=f3
                    f10=4.0_prec*f2
                    kmin=k
                END IF
            END DO loop1
            f=f+f2
        END DO

        RETURN

    END SUBROUTINE clustfull2


END MODULE obj_fun
