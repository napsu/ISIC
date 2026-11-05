!*************************************************************************
!*                                                                       *
!*     ISIC - Stochastic Nonsmooth Optimization based                    *
!*     Incremental Clustering Software                                   *
!*                                                                       *
!*     Subroutines for clustering used within ISIC software.             *
!*                                                                       *
!*     Napsu Karmitsa 2026 (last modified 14.4.2026).                    *
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
!*     clusteringmod       ! 
!* 

MODULE clusteringmod       ! Subroutines for clustering software

    USE r_precision, ONLY : prec      ! Precision for reals.
    IMPLICIT NONE

    ! MODULE clusteringmod includes the following subroutines (S).

    PUBLIC :: &
        check, &           ! S  Checking the results and calculation of
                           !    validity indices.
        step1, &           ! S  Computation of the centroid and the value of
                           !    clustering function at centroid.
        step2              ! S  Select the best candidate solution.

CONTAINS

    !============================================================================
    !    Subroutine check summarizes the results, selects data points as initial
    !    cluster centers, and calculates cluster validity indices.
    !============================================================================
    SUBROUTINE check(x,f,db,db2,dn,dn3,sep)

        USE param, ONLY : zero, one, two, small, large
        USE initclust, ONLY : &
            a => a_used, &    ! Data matrix.
            nrecord => nrec_used, &   ! Number of data points.
            nft, &            ! Number of features.
            nc, &             ! Current number of clusters, loops from 1 to nclust.
            tnorm, &          ! The total number of norms computed.
            m1, &             ! Maximum number of initial solutions.
            nel, &            ! nel(nclust), nel(i)=number of records in cluster i
            nk, &             ! nk(nrecord,nclust), nk(j,i) (j=1,...,nel(i)) tells which data 
                              !   points are in cluster i
            ncand, &          ! Number of canditate points.
            lcand, &          ! lcand(nrecord), lcand(i) 
            list1, &          ! list1(nrecord), list1(i)=the cluster where point i belongs
            dminim, &         ! dminim(nrecord), the distance of a(i) and the nearest centroid
            dcent             ! dcent(nclust,nclust), Distance (affinity) matrix for cluster centers
 
        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) ::      &
            x                     ! Cluster centers
        REAL(KIND=prec), INTENT(OUT) ::      &
            f, &                  ! Value of the clustering function 
            db,db2,dn,dn3,sep     ! Davies-Boulding, Dunn, and separation indices

        ! Local variables
        REAL(KIND=prec) ::      &
            f10,                &
            f1,f2,              &
            fdb,                &
            fk2,fk3,            &
            fm,fm3,             &
            dc(nc,nc),          & ! Distances between cluster centers.
            dc1,dn1,dn2,radm,   &
            fk(nc),             &
            rad(nc),            & ! Average distance of data points in cluster i to its centroid.
            radmax(nc),         & ! Maximum distance of data points in cluster i to its centroid.
            trdist,             & ! 
            r1,r2,r3,r4,r5

        ! Local integer variables
        INTEGER ::  &
            i,j,k, &
            jmin, &
            k1,     &
            ncand1, &
            n1

        trdist = 4.0_prec
        !trdist = MAX(4.0_prec, LOG(REAL(nft,prec))) 
        
        nel = 0
        DO j=1,nc
            rad(j)=zero
            radmax(j)=zero
        END DO
        
        f = zero
        jmin = 1

        !=========================================================================
        ! Identifying data points that are cluster centers.
        !=========================================================================
        outerloop: DO k=1,nrecord
            f2=large
            f10=large
            innerloop: DO j=1,nc
                IF(j > 1) THEN
                    IF(dcent(j,jmin) >= f10) CYCLE innerloop
                END IF
                f1=zero
                DO k1=1,nft
                    f1=f1+(a(k1,k)-x(k1+(j-1)*nft))**2
                END DO
                tnorm=tnorm+one
                IF(f1 < f2) THEN
                    f2=f1
                    f10 = trdist * f1
                    jmin=j
                END IF
            END DO innerloop
            dminim(k)=f2
            f=f+f2
            nel(jmin)=nel(jmin)+1
            nk(nel(jmin),jmin)=k
            list1(k)=jmin
            rad(jmin)=rad(jmin)+f2
            radmax(jmin)=MAX(radmax(jmin),f2)

        END DO outerloop
        
        !=========================================================================
        ! Computing average distance of data points in cluster k to its centroid.
        !=========================================================================
        DO k=1,nc 
            IF(nel(k) > 0) THEN
                rad(k)=rad(k)/REAL(nel(k),prec)
            ELSE
                rad(k)=zero               
            END IF
        END DO
        
        !=========================================================================
        ! Selecting data points as initial centers.
        !=========================================================================
        ncand=0
        loop1: DO k=1,nc ! Step 1 in Finding set of starting points algorithm, returns set A_1.
            IF(nel(k) == 1) CYCLE loop1 ! If only one point in the cluster it is already the centroid
            IF(radmax(k) <= 1.0E-1_prec*f/nrecord) CYCLE loop1 ! If radmax is small the point is the centroid
            n1=MAX(1,nel(k)*m1/nrecord) ! Max number of data points selected from cluster k.
            ncand1=0
            r1=radmax(k)
            r2=rad(k)
            r5=1.0E-03_prec*(r1-r2)
            IF(ABS(r1-r2) < 1.0E-03_prec) r5=1.0E-06_prec
            r3=r1+r5
            DO k1=1,10000 ! Limit to prevent endless loop, this should never happen
                IF (k1 == 10000) PRINT*,'Warning: loop limit reached in starting point procedure.'
                r3=r3-r5
                r4=r3-r5
                IF(r3 < r2) CYCLE loop1 ! Only loop points whose distance is above the average
                DO j=1,nel(k)
                    i=nk(j,k)           ! Check data points that are in cluster k
                    IF((dminim(i) <= r3).and.(dminim(i) >= r4)) THEN ! Add the furthest points as canditate solutions
                        ncand=ncand+1
                        ncand1=ncand1+1
                        lcand(ncand)=i
                        IF(ncand1 >= n1) CYCLE loop1
                    END IF
                END DO
            END DO
        END DO loop1

        !=========================================================================
        ! Calculation of distances between cluster centers.
        !=========================================================================
        do i=1,nc
            dc(i,i) = zero
        end do

        do i=1,nc
            do j=i+1,nc
                dc1 = zero
                do k=1,nft
                    dc1=dc1+(x(k+(i-1)*nft)-x(k+(j-1)*nft))**2
                end do
                dc(i,j)=SQRT(dc1)
                dc(j,i)=dc(i,j)
            end do
        end do


        !=====================================================
        ! Calculation of Davies-Bouldin (DB) validity index
        ! (last changed 25 August 2017)
        ! db is the index using Euclidean distances while
        ! db2 uses squared Euclidean distances
        !=====================================================

        DO i=1,nc
            fk(i)=zero
        END DO

        fdb=zero
        db2=zero

        DO i=1,nc
            fk(i)=SQRT(rad(i))
        END DO


        DO k=1,nc
            fm=zero
            fm3=zero
            DO i=1,nc
                IF (i.ne.k) THEN
                    fk2=fk(i)+fk(k)
                    fk3=(rad(i)+rad(k))/(dc(i,k)*dc(i,k))
                    f2=fk2/dc(i,k)
                    fm=MAX(fm,f2)
                    fm3=MAX(fm3,fk3)
                END IF
            END DO
            fdb=fdb+fm
            db2=db2+fm3
        END DO
        db=fdb/REAL(nc,prec)
        db2=db2/REAL(nc,prec)


        !============================================================
        ! Calculation of Dunn validity index
        !============================================================
        radm= zero
        DO i=1,nc
            radm=MAX(radm,radmax(i)) ! maximum radius of all clusters
        END DO
        radm = SQRT(radm)


        dn3 = zero
        dn = large
        DO i=1,nc
            dn1 = large
            DO j=1,nc
                IF (j.NE.i) THEN
                    dn2=dc(i,j)/radm
                    dn1=MIN(dn2,dn1) ! distance of cluster i to the closest cluster j
                END IF
            END DO
            dn=MIN(dn,dn1)
            dn3=MAX(dn,dn1)
        END DO
        IF (nc == 1) dn=zero
        IF (nc == 1) dn3=zero

        !============================================================
        ! Calculation of the quality of separation
        !============================================================
        sep=zero
        !        do i=1,nc
        !            m2=nel(i)
        !            m3=0
        !            loop_sep: do j=1,m2
        !                k=nk(j,i)
        !                do k1=1,nc
        !                    if(k1.ne.i) then
        !                        d1=zero
        !                        do k2=1,nft
        !                            d1=d1+(a(k2,k)-x(k2+(k1-1)*nft))**2
        !                        end do
        !                        if(d1<radmax(k1)) then
        !                            m3=m3+1
        !                            CYCLE loop_sep
        !                        end if
        !                    end if
        !                end do
        !                end do loop_sep
        !                sep=sep+REAL(m3,prec)/REAL(nel(i),prec)
        !            end do
        !            sep=sep/REAL(nc,prec)
        !============================================================


        RETURN

    END SUBROUTINE check


    !===============================================================================
    !  Subroutine step1 computes the centroid and the value of clustering function at 
    !  centroid using sample of data (but not minibatch, Step 0 in ISIC algorithm)
    !===============================================================================
    SUBROUTINE step1(f,x)

        USE param, ONLY : zero, one
        USE initclust, ONLY : &
            a => a_used, &   ! Data matrix
            nft, &           ! Number of features
            nrecord => nrec_used, &     ! Number of data points
            tnorm            ! The total number of norms computed.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec), DIMENSION(nft), INTENT(OUT) ::      &
            x              ! Cluster centers
        REAL(KIND=prec), INTENT(OUT) ::  &
            f              ! Cluster function value

        ! Local variables
        REAL(KIND=prec) ::  &
            f1
        INTEGER ::      &
            i,j
        
        x = zero
        f = zero
        DO i=1,nrecord
            DO j=1,nft
                x(j)=x(j)+a(j,i) 
            END DO
        END DO
        x = x/REAL(nrecord,prec)

        DO i=1,nrecord
            f1=zero
            tnorm=tnorm+one
            DO j=1,nft
                f1=f1+(a(j,i)-x(j))*(a(j,i)-x(j))
            END DO
            f=f+f1
                   
        END DO
        RETURN

    END SUBROUTINE step1


    !=========================================================================
    !  Subroutine step2: Select the best candidate solution
    !=========================================================================
    SUBROUTINE step2(x2)

        USE param, ONLY : zero, one, large
        USE initclust, ONLY : &
            a => a_used, & ! Data matrix
            nft, &         ! Number of features
            nrecord => nrec_used, &     ! Number of records in data.
            dminim, &      ! dminim(nrecord), the distance of a(i) and the nearest centroid.
            ncand, &       ! Number of canditate points.
            lcand, &       ! Data points which are selected as initial cluster centers.
            l4, &          ! Auxiliary array
            tnorm, &       ! The total number of norms computed.
            iphi, &        ! 0 = full-sample Phi, 1 = subsampled Phi.
            mphi           ! Evaluation-sample size when iphi = 1.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(INOUT) ::  &
            x2             ! Initial cluster centers (x2(1:nft) = center of best candidate)

        ! Local variables
        REAL(KIND=prec) ::  &
            d1,d3,d21,  &
            fmin,       &
            f2,         &
            trdist,     &
            urand

        INTEGER ::  &
            ibest,  &
            nclose, &
            i1,     &
            i,j,    &
            k,l,    &
            mphi_use, &
            itmp,   &
            jswap

        INTEGER, ALLOCATABLE :: &
            iphi_ind(:)     ! Indices of the evaluation sample A_Phi.

        INTRINSIC :: MIN, MAX

        trdist = MAX(4.0_prec, LOG(REAL(nft,prec))) 

        IF (ncand == 0) THEN 
            i = 1
            f2 = dminim(1)
            DO j = 2, nrecord
                IF (dminim(j) > f2) THEN
                    f2 = dminim(j)
                    i = j
                END IF
            END DO
            ncand = 1
            lcand(1) = i
        END IF

        IF (ncand > 1) THEN
            fmin   = large
            ibest  = 1

            IF (iphi == 0) THEN
                ! Full-sample Phi: evaluate each candidate against all of A_S.
                DO i1=1,ncand
                    i   = lcand(i1)
                    f2  = dminim(i)
                    d21 = zero
                    loop1: DO l=1,nrecord
                        IF (f2 < (trdist * dminim(l))) THEN
                            d3 = zero
                            DO k=1,nft
                                d3 = d3 + (a(k,i)-a(k,l))**2
                            END DO
                            tnorm = tnorm + one
                            d21   = d21 + MIN(zero, d3-dminim(l))
                        END IF
                    END DO loop1
                    IF (d21 < fmin) THEN
                        fmin  = d21
                        ibest = i1
                    END IF
                END DO
            ELSE
                ! Subsampled Phi: evaluate each candidate on a uniform subsample A_Phi of A_S.
                mphi_use = MAX(1, MIN(mphi, nrecord))
                ALLOCATE(iphi_ind(nrecord))
                DO l=1,nrecord
                    iphi_ind(l) = l
                END DO
                IF (mphi_use < nrecord) THEN
                    ! Partial Fisher-Yates: first mphi_use entries are a uniform sample of 1..nrecord.
                    DO l=1,mphi_use
                        CALL random_number(urand)
                        jswap = l + INT(urand * REAL(nrecord - l + 1, prec))
                        IF (jswap < l) jswap = l
                        IF (jswap > nrecord) jswap = nrecord
                        itmp = iphi_ind(l)
                        iphi_ind(l) = iphi_ind(jswap)
                        iphi_ind(jswap) = itmp
                    END DO
                END IF

                DO i1=1,ncand
                    i   = lcand(i1)
                    f2  = dminim(i)
                    d21 = zero
                    loop_phi: DO j=1,mphi_use
                        l = iphi_ind(j)
                        IF (f2 < (trdist * dminim(l))) THEN
                            d3 = zero
                            DO k=1,nft
                                d3 = d3 + (a(k,i)-a(k,l))**2
                            END DO
                            tnorm = tnorm + one
                            d21   = d21 + MIN(zero, d3-dminim(l))
                        END IF
                    END DO loop_phi
                    IF (d21 < fmin) THEN
                        fmin  = d21
                        ibest = i1
                    END IF
                END DO
                DEALLOCATE(iphi_ind)
            END IF
        ELSE
            ibest = 1
        END IF

        ! Compute the centroid of the ball B(a) around the best candidate.
        i      = lcand(ibest)
        nclose = 0
        DO j=1,nrecord
            IF (dminim(i) < (trdist*dminim(j))) THEN
                d1 = zero
                DO k=1,nft
                    d1 = d1 + (a(k,i)-a(k,j))**2
                END DO
                tnorm = tnorm + one
                IF (d1 < 0.6_prec*dminim(j)) THEN
                !IF (d1 < dminim(j)) THEN ! testi 5.9.
                    nclose     = nclose + 1
                    l4(nclose) = j
                END IF
            END IF
        END DO
        IF (nclose > 0) THEN
            DO k=1,nft
                d3 = zero
                DO j=1,nclose
                    d3 = d3 + a(k, l4(j))
                END DO
                x2(k) = d3 / REAL(nclose, prec)
            END DO
        ELSE 
            DO k = 1, nft
                x2(k) = a(k, i)
            END DO
        END IF

        RETURN

    END SUBROUTINE step2

END MODULE clusteringmod
