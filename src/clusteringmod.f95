!*************************************************************************
!*                                                                       *
!*     ISIC - Stochastic Nonsmooth Optimization based                    *
!*     Incremental Clustering Software                                   *
!*                                                                       *
!*     Subroutines for clustering used within ISIC software.             *
!*                                                                       *
!*     Napsu Karmitsa 2026.                                              *
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
            dcent, &          ! dcent(nclust,nclust), Distance (affinity) matrix for cluster centers
            iimp, &           ! 0 = shell sweep, 1 = global sort, 2 = per-cluster sort.
            ord, &            ! Indices of data points in decreasing dminim order
            rtail             ! Suffix sums of dminim along ord
 
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
            trdist,             &
            r1,r2,r3,r4,r5

        ! Local integer variables
        INTEGER ::  &
            i,j,k, &
            jmin, &
            k1,     &
            s,      &
            n1,     &
            ncand1, &
            quota(nc), taken(nc), cdone(nc), head(nc), tail(nc)

        INTEGER, ALLOCATABLE :: nxt(:)

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

        IF (iimp == 1) THEN
            ! Sort working-sample points by dminim (descending) and suffix sums.
            DO i=1,nrecord
                ord(i) = i
            END DO
            CALL sort_ind_by_dminim_desc(ord, nrecord)
            rtail(nrecord+1) = zero
            DO s=nrecord,1,-1
                rtail(s) = rtail(s+1) + dminim(ord(s))
            END DO

            !=========================================================================
            ! Selecting data points as initial centers (furthest above-average per cluster).
            ! One pass of ord; lcand stays cluster-major, each cluster in decreasing r.
            !=========================================================================
            DO k=1,nc
                head(k)  = 0
                taken(k) = 0
                cdone(k) = 1
                quota(k) = 0
                IF (nel(k) > 1) THEN
                    IF (radmax(k) > 1.0E-1_prec*f/nrecord) THEN
                        cdone(k) = 0
                        quota(k) = MAX(1, nel(k)*m1/nrecord)
                    END IF
                END IF
            END DO

            ALLOCATE(nxt(nrecord))
            DO s=1,nrecord
                i = ord(s)
                k = list1(i)
                IF (cdone(k) == 1) CYCLE
                IF (dminim(i) < rad(k)) THEN
                    cdone(k) = 1          ! Later members of k are even closer
                    CYCLE
                END IF
                IF (head(k) == 0) THEN
                    head(k) = i
                ELSE
                    nxt(tail(k)) = i
                END IF
                tail(k) = i
                nxt(i)  = 0
                taken(k) = taken(k) + 1
                IF (taken(k) >= quota(k)) cdone(k) = 1
            END DO

            ncand = 0
            DO k=1,nc
                i = head(k)
                DO WHILE (i /= 0)
                    ncand = ncand + 1
                    lcand(ncand) = i
                    i = nxt(i)
                END DO
            END DO
            DEALLOCATE(nxt)
        ELSE IF (iimp == 2) THEN
            !=========================================================================
            ! Per-cluster decreasing dminim (pi_j). Furthest above-average candidates.
            ! No global ord / rtail. Ties may differ from the iimp == 1 heapsort.
            !=========================================================================
            DO k=1,nc
                IF (nel(k) > 1) CALL sort_cluster_by_dminim_desc(k)
            END DO
            ncand = 0
            DO k=1,nc
                IF (nel(k) <= 1) CYCLE
                IF (radmax(k) <= 1.0E-1_prec*f/nrecord) CYCLE
                n1 = MAX(1, nel(k)*m1/nrecord)
                ncand1 = 0
                DO j=1,nel(k)
                    i = nk(j,k)
                    IF (dminim(i) < rad(k)) EXIT
                    ncand = ncand + 1
                    lcand(ncand) = i
                    ncand1 = ncand1 + 1
                    IF (ncand1 >= n1) EXIT
                END DO
            END DO
        ELSE
            !=========================================================================
            ! Selecting data points as initial centers (shell sweep, no sort).
            !=========================================================================
            ncand=0
            loop1: DO k=1,nc
                IF(nel(k) == 1) CYCLE loop1
                IF(radmax(k) <= 1.0E-1_prec*f/nrecord) CYCLE loop1
                n1=MAX(1,nel(k)*m1/nrecord)
                ncand1=0
                r1=radmax(k)
                r2=rad(k)
                r5=1.0E-03_prec*(r1-r2)
                IF(ABS(r1-r2) < 1.0E-03_prec) r5=1.0E-06_prec
                r3=r1+r5
                DO k1=1,10000
                    IF (k1 == 10000) PRINT*,'Warning: loop limit reached in starting point procedure.'
                    r3=r3-r5
                    r4=r3-r5
                    IF(r3 < r2) CYCLE loop1
                    DO j=1,nel(k)
                        i=nk(j,k)
                        IF((dminim(i) <= r3).and.(dminim(i) >= r4)) THEN
                            ncand=ncand+1
                            ncand1=ncand1+1
                            lcand(ncand)=i
                            IF(ncand1 >= n1) CYCLE loop1
                        END IF
                    END DO
                END DO
            END DO loop1
        END IF

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
    SUBROUTINE step2(x, x2)

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
            mphi, &        ! Evaluation-sample size when iphi = 1.
            iimp, &        ! 0 = nested Phi loops, 1 = ImprovementPotential.
            ord            ! Indices of data points in decreasing dminim order.

        IMPLICIT NONE

        ! Subroutine arguments
        REAL(KIND=prec), DIMENSION(:), INTENT(IN) ::  &
            x              ! Current cluster centers, centers 1..nc-1. Same array as check.
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

!        trdist = MAX(4.0_prec, LOG(REAL(nft,prec))) 
        trdist = 4.0_prec ! 11.11.2026

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

        IF (iimp == 2) THEN
            CALL improvement_potential_acc(x, ibest, nclose)
        ELSE IF (ncand > 1) THEN
            fmin   = large
            ibest  = 1

            IF (iimp == 0) THEN
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
                CALL improvement_potential(trdist, ibest)
            END IF
        ELSE
            ibest = 1
        END IF

        ! Compute the centroid of the ball B(a) around the best candidate.
        ! iimp == 2 already filled l4 and nclose by one more pass on the winner.
        i       = lcand(ibest)
        IF (iimp /= 2) THEN
        nclose  = 0
        f2      = dminim(i)
        IF (iimp == 1) THEN
            ! Walk ord (decreasing r); stop when r* >= tau*r_j (cannot enter B).
            DO l=1,nrecord
                j = ord(l)
                IF (f2 >= trdist*dminim(j)) EXIT
                d1 = zero
                DO k=1,nft
                    d1 = d1 + (a(k,i)-a(k,j))**2
                END DO
                tnorm = tnorm + one
                IF (d1 < 0.6_prec*dminim(j)) THEN
                    nclose     = nclose + 1
                    l4(nclose) = j
                END IF
            END DO
            IF (nclose > 1) CALL sort_int_asc(l4, nclose)
        ELSE
            ! Sequential scan in index order. No sort of l4.
            DO j=1,nrecord
                IF (f2 < (trdist*dminim(j))) THEN
                    d1 = zero
                    DO k=1,nft
                        d1 = d1 + (a(k,i)-a(k,j))**2
                    END DO
                    tnorm = tnorm + one
                    IF (d1 < 0.6_prec*dminim(j)) THEN
                        nclose     = nclose + 1
                        l4(nclose) = j
                    END IF
                END IF
            END DO
        END IF
        END IF
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


    !=========================================================================
    !  ImprovementPotential: data-major ranking with hard cutoff (A),
    !  bound elimination (B), and global early stop (C).
    !  Evaluation set is A_S (iphi=0) or a uniform subsample A_Phi (iphi=1).
    !=========================================================================
    SUBROUTINE improvement_potential(trdist, ibest)

        USE param, ONLY : zero, one
        USE initclust, ONLY : &
            a => a_used, &
            nft, &
            nrecord => nrec_used, &
            dminim, &
            ncand, &
            lcand, &
            tnorm, &
            iphi, &
            mphi, &
            ord, &
            rtail

        IMPLICIT NONE

        REAL(KIND=prec), INTENT(IN) :: trdist
        INTEGER, INTENT(OUT) :: ibest

        REAL(KIND=prec) ::      &
            phi_best, d3, urand, r_pass, rj
        INTEGER ::              &
            meval, mphi_use,    &
            i, j, k, s, i1, ic, jp, nact, itmp, jswap, use_pack
        INTEGER, ALLOCATABLE :: &
            eval_ind(:),        &
            iact(:),            &
            s_freeze(:),        &
            cord(:)
        REAL(KIND=prec), ALLOCATABLE :: &
            phi_c(:), rtail_loc(:), aphi(:,:)

        INTRINSIC :: MIN, MAX

        use_pack = 0
        mphi_use = MAX(1, MIN(mphi, nrecord))
        IF (iphi == 0 .OR. mphi_use >= nrecord) THEN
            meval = nrecord
        ELSE
            meval = mphi_use
            ALLOCATE(eval_ind(nrecord))
            DO i=1,nrecord
                eval_ind(i) = i
            END DO
            IF (meval < nrecord) THEN
                DO i=1,meval
                    CALL random_number(urand)
                    jswap = i + INT(urand * REAL(nrecord - i + 1, prec))
                    IF (jswap < i) jswap = i
                    IF (jswap > nrecord) jswap = nrecord
                    itmp = eval_ind(i)
                    eval_ind(i) = eval_ind(jswap)
                    eval_ind(jswap) = itmp
                END DO
            END IF
            CALL sort_ind_by_dminim_desc(eval_ind, meval)
            ALLOCATE(rtail_loc(meval+1), aphi(nft, meval))
            rtail_loc(meval+1) = zero
            DO s=meval,1,-1
                rtail_loc(s) = rtail_loc(s+1) + dminim(eval_ind(s))
                DO j=1,nft
                    aphi(j,s) = a(j, eval_ind(s))
                END DO
            END DO
            use_pack = 1
        END IF

        ALLOCATE(iact(ncand), phi_c(ncand), s_freeze(ncand), cord(ncand))

        DO i1=1,ncand
            phi_c(i1)    = zero
            iact(i1)     = i1
            s_freeze(i1) = meval + 1
            cord(i1)     = i1
        END DO
        ! Larger r freezes sooner; one forward pointer along the sorted sample.
        CALL sort_cand_by_r_desc(cord, ncand)
        s = 1
        DO i1=1,ncand
            ic = cord(i1)
            DO WHILE (s <= meval)
                IF (use_pack == 0) THEN
                    rj = dminim(ord(s))
                ELSE
                    rj = dminim(eval_ind(s))
                END IF
                IF (dminim(lcand(ic)) >= trdist * rj) EXIT
                s = s + 1
            END DO
            s_freeze(ic) = s
        END DO

        nact     = ncand
        phi_best = zero
        ibest    = 1

        outer: DO s=1,meval
            IF (use_pack == 0) THEN
                jp = ord(s)
                rj = dminim(jp)
            ELSE
                jp = s
                rj = dminim(eval_ind(s))
            END IF

            ! Mechanism A: hard cutoff (process active list backwards).
            k = nact
            DO WHILE (k >= 1)
                ic = iact(k)
                IF (dminim(lcand(ic)) >= trdist * rj) THEN
                    iact(k) = iact(nact)
                    nact    = nact - 1
                    k       = k - 1
                ELSE
                    i  = lcand(ic)
                    d3 = zero
                    IF (use_pack == 0) THEN
                        DO j=1,nft
                            d3 = d3 + (a(j,i)-a(j,jp))**2
                        END DO
                    ELSE
                        DO j=1,nft
                            d3 = d3 + (a(j,i)-aphi(j,s))**2
                        END DO
                    END IF
                    tnorm     = tnorm + one
                    phi_c(ic) = phi_c(ic) + MIN(zero, d3-rj)
                    IF (phi_c(ic) < phi_best) phi_best = phi_c(ic)
                    k = k - 1
                END IF
            END DO

            ! Mechanism B: bound on remaining points that still pass tau.
            k = nact
            DO WHILE (k >= 1)
                ic = iact(k)
                IF (use_pack == 0) THEN
                    r_pass = rtail(s+1) - rtail(s_freeze(ic))
                ELSE
                    r_pass = rtail_loc(s+1) - rtail_loc(s_freeze(ic))
                END IF
                IF (r_pass < zero) r_pass = zero
                IF (phi_c(ic) - r_pass > phi_best) THEN
                    iact(k) = iact(nact)
                    nact    = nact - 1
                    k       = k - 1
                ELSE
                    k = k - 1
                END IF
            END DO

            ! Mechanism C: global stop.
            IF (nact == 0) EXIT outer
        END DO outer

        ibest = 1
        DO i1=2,ncand
            IF (phi_c(i1) < phi_c(ibest)) ibest = i1
        END DO

        DEALLOCATE(iact, phi_c, s_freeze, cord)
        IF (use_pack == 1) THEN
            DEALLOCATE(eval_ind, rtail_loc, aphi)
        END IF

    END SUBROUTINE improvement_potential


    ! In-place heapsort: dminim(ind(1)) >= ... >= dminim(ind(n)).
    SUBROUTINE sort_ind_by_dminim_desc(ind, n)

        USE initclust, ONLY : dminim

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(INOUT) :: ind(:)

        INTEGER :: i, start, bottom, tmp

        IF (n <= 1) RETURN

        DO start = n/2, 1, -1
            CALL siftdown(start, n)
        END DO
        DO bottom = n, 2, -1
            tmp         = ind(1)
            ind(1)      = ind(bottom)
            ind(bottom) = tmp
            CALL siftdown(1, bottom-1)
        END DO
        DO i=1,n/2
            tmp          = ind(i)
            ind(i)       = ind(n+1-i)
            ind(n+1-i)   = tmp
        END DO

    CONTAINS

        SUBROUTINE siftdown(start, bottom)
            INTEGER, INTENT(IN) :: start, bottom
            INTEGER :: root, child, tmp
            root = start
            DO
                child = 2*root
                IF (child > bottom) EXIT
                IF (child+1 <= bottom) THEN
                    IF (dminim(ind(child+1)) > dminim(ind(child))) THEN
                        child = child + 1
                    END IF
                END IF
                IF (dminim(ind(root)) >= dminim(ind(child))) EXIT
                tmp        = ind(root)
                ind(root)  = ind(child)
                ind(child) = tmp
                root       = child
            END DO
        END SUBROUTINE siftdown

    END SUBROUTINE sort_ind_by_dminim_desc


    ! Ascending heapsort of the first n integer entries (index order for ball B).
    SUBROUTINE sort_int_asc(ind, n)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(INOUT) :: ind(:)

        INTEGER :: start, bottom, tmp

        IF (n <= 1) RETURN

        DO start = n/2, 1, -1
            CALL siftdown_int(start, n)
        END DO
        DO bottom = n, 2, -1
            tmp         = ind(1)
            ind(1)      = ind(bottom)
            ind(bottom) = tmp
            CALL siftdown_int(1, bottom-1)
        END DO

    CONTAINS

        SUBROUTINE siftdown_int(start, bottom)
            INTEGER, INTENT(IN) :: start, bottom
            INTEGER :: root, child, tmp
            root = start
            DO
                child = 2*root
                IF (child > bottom) EXIT
                IF (child+1 <= bottom) THEN
                    IF (ind(child+1) > ind(child)) child = child + 1
                END IF
                IF (ind(root) >= ind(child)) EXIT
                tmp        = ind(root)
                ind(root)  = ind(child)
                ind(child) = tmp
                root       = child
            END DO
        END SUBROUTINE siftdown_int

    END SUBROUTINE sort_int_asc


    ! Sort candidate slots so dminim(lcand(ind(1))) >= ... >= dminim(lcand(ind(n))).
    SUBROUTINE sort_cand_by_r_desc(ind, n)

        USE initclust, ONLY : dminim, lcand

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n
        INTEGER, INTENT(INOUT) :: ind(:)

        INTEGER :: i, start, bottom, tmp

        IF (n <= 1) RETURN

        DO start = n/2, 1, -1
            CALL siftdown_cand(start, n)
        END DO
        DO bottom = n, 2, -1
            tmp         = ind(1)
            ind(1)      = ind(bottom)
            ind(bottom) = tmp
            CALL siftdown_cand(1, bottom-1)
        END DO
        DO i=1,n/2
            tmp         = ind(i)
            ind(i)      = ind(n+1-i)
            ind(n+1-i)  = tmp
        END DO

    CONTAINS

        SUBROUTINE siftdown_cand(start, bottom)
            INTEGER, INTENT(IN) :: start, bottom
            INTEGER :: root, child, tmp
            root = start
            DO
                child = 2*root
                IF (child > bottom) EXIT
                IF (child+1 <= bottom) THEN
                    IF (dminim(lcand(ind(child+1))) > dminim(lcand(ind(child)))) THEN
                        child = child + 1
                    END IF
                END IF
                IF (dminim(lcand(ind(root))) >= dminim(lcand(ind(child)))) EXIT
                tmp        = ind(root)
                ind(root)  = ind(child)
                ind(child) = tmp
                root       = child
            END DO
        END SUBROUTINE siftdown_cand

    END SUBROUTINE sort_cand_by_r_desc


    !=========================================================================
    !  iimp = 2. Per-cluster cutoff, U(a) branch-and-bound, early-abandon
    !  distances. wmax <= 0 scores every prefix point (theta = 1); wmax > 0
    !  draws a Bernoulli mask and scales the score. The search keeps only the
    !  scalar score. Ball B is rebuilt by one more pass on the winner.
    !  step2 runs with nc already incremented, so the assignment has nc-1 centers.
    !=========================================================================
    SUBROUTINE improvement_potential_acc(x, ibest, nclose)

        USE param, ONLY : zero, one, large
        USE initclust, ONLY : &
            a => a_used, &
            nft, &
            nrecord => nrec_used, &
            nc, &
            dminim, &
            nel, &
            nk, &
            ncand, &
            lcand, &
            l4, &
            tnorm, &
            wmax

        IMPLICIT NONE

        REAL(KIND=prec), DIMENSION(:), INTENT(IN) :: x   ! Centers 1..nc-1 from check.
        INTEGER, INTENT(OUT) :: ibest, nclose

        REAL(KIND=prec) ::      &
            phi, phi_best, g2, g, uk, delta, rj, acc, &
            sumv(nft), scale(nft), ac(nft), &
            w_tot, theta, urand
        INTEGER ::              &
            ncent,              &
            i1, ic, j, k, p, h, &
            ia, ii, nj, npos, ip, &
            kappa(nft)
        INTEGER, ALLOCATABLE :: &
            cord(:), njcut(:,:), pbase(:)
        LOGICAL, ALLOCATABLE :: &
            omega(:)
        REAL(KIND=prec), ALLOCATABLE :: &
            u(:), ps(:)

        INTRINSIC :: SQRT

        ncent    = nc - 1
        ibest    = 1
        nclose   = 0
        phi_best = large
        IF (ncand < 1 .OR. ncent < 1) RETURN

        ! Largest coordinate variance first. One pass per step2; both distance
        ! passes use kappa. Order does not change a finished sum of squares.
        sumv  = zero
        scale = zero
        DO ip=1,nrecord
            DO j=1,nft
                sumv(j)  = sumv(j) + a(j,ip)
                scale(j) = scale(j) + a(j,ip)*a(j,ip)
            END DO
        END DO
        DO j=1,nft
            scale(j) = scale(j) - sumv(j)*sumv(j)/REAL(nrecord,prec)
            kappa(j) = j
        END DO
        CALL sort_kappa_desc()

        ALLOCATE(cord(ncand), u(ncand), njcut(ncand, ncent), pbase(ncent), ps(nrecord))

        npos = 0
        DO j=1,ncent
            pbase(j) = npos + 1
            acc = zero
            DO p=1,nel(j)
                acc = acc + SQRT(dminim(nk(p,j)))
                npos = npos + 1
                ps(npos) = acc
            END DO
        END DO

        DO ic=1,ncand
            ia    = lcand(ic)
            u(ic) = zero
            DO j=1,ncent
                IF (nel(j) == 0) THEN
                    njcut(ic,j) = 0
                    CYCLE
                END IF
                g2 = zero
                DO k=1,nft
                    g2 = g2 + (a(k,ia) - x(k+(j-1)*nft))**2
                END DO
                nj = prefix_len(j, g2)
                njcut(ic,j) = nj
                IF (nj > 0) THEN
                    g = SQRT(g2)
                    u(ic) = u(ic) + (2.0_prec*g*ps(pbase(j)+nj-1) - REAL(nj,prec)*g2)
                END IF
            END DO
            cord(ic) = ic
        END DO

        DO i1=2,ncand
            ic = cord(i1)
            uk = u(ic)
            k  = i1
            DO WHILE (k > 1)
                IF (u(cord(k-1)) >= uk) EXIT
                cord(k) = cord(k-1)
                k = k - 1
            END DO
            cord(k) = ic
        END DO

        w_tot = zero
        DO ic=1,ncand
            DO j=1,ncent
                w_tot = w_tot + REAL(njcut(ic,j), prec)
            END DO
        END DO
        IF (REAL(wmax, prec) <= zero .OR. w_tot == zero) THEN
            theta = one
        ELSE
            theta = MIN(one, REAL(wmax, prec) / w_tot)
        END IF
        IF (theta < one) THEN
            ALLOCATE(omega(nrecord))
            DO ip=1,nrecord
                CALL random_number(urand)
                omega(ip) = urand < theta
            END DO
        END IF

        DO i1=1,ncand
            ic = cord(i1)
            IF (phi_best < large) THEN
                IF (theta < one) THEN
                    IF (u(ic) <= -theta*phi_best) EXIT
                ELSE
                    IF (u(ic) <= -phi_best) EXIT
                END IF
            END IF
            ia   = lcand(ic)
            DO h=1,nft
                ac(h) = a(kappa(h), ia)
            END DO
            phi  = zero
            DO j=1,ncent
                nj = njcut(ic,j)
                DO p=1,nj
                    ii = nk(p,j)
                    IF (theta < one) THEN
                        IF (.NOT. omega(ii)) CYCLE
                    END IF
                    rj    = dminim(ii)
                    delta = bounded_delta(ac, ii, rj)
                    IF (delta < rj) THEN
                        IF (theta < one) THEN
                            phi = phi + (delta - rj) / theta
                        ELSE
                            phi = phi + (delta - rj)
                        END IF
                    END IF
                END DO
            END DO
            IF (phi < phi_best) THEN
                phi_best = phi
                ibest    = ic
            END IF
        END DO

        ia     = lcand(ibest)
        DO h=1,nft
            ac(h) = a(kappa(h), ia)
        END DO
        nclose = 0
        DO j=1,ncent
            nj = njcut(ibest,j)
            DO p=1,nj
                ii    = nk(p,j)
                rj    = dminim(ii)
                delta = bounded_delta(ac, ii, rj)
                IF (delta < 0.6_prec*rj) THEN
                    nclose     = nclose + 1
                    l4(nclose) = ii
                END IF
            END DO
        END DO
        IF (nclose > 1) CALL sort_int_asc(l4, nclose)

        DEALLOCATE(cord, u, njcut, pbase, ps)
        IF (ALLOCATED(omega)) DEALLOCATE(omega)

    CONTAINS

        FUNCTION bounded_delta(ac, ii, rj) RESULT(delta)
            IMPLICIT NONE
            REAL(KIND=prec), INTENT(IN) :: ac(nft), rj
            INTEGER, INTENT(IN) :: ii
            REAL(KIND=prec) :: delta
            INTEGER :: h, kf
            delta = zero
            h = 1
            DO WHILE (h <= nft .AND. delta < rj)
                kf = kappa(h)
                delta = delta + (ac(h) - a(kf,ii))**2
                h = h + 1
            END DO
            tnorm = tnorm + one
        END FUNCTION bounded_delta

        FUNCTION prefix_len(j, g2) RESULT(nj)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: j
            REAL(KIND=prec), INTENT(IN) :: g2
            INTEGER :: nj, lo, hi, mid
            nj = 0
            lo = 1
            hi = nel(j)
            DO WHILE (lo <= hi)
                mid = (lo + hi) / 2
                IF (4.0_prec * dminim(nk(mid,j)) > g2) THEN
                    nj = mid
                    lo = mid + 1
                ELSE
                    hi = mid - 1
                END IF
            END DO
        END FUNCTION prefix_len

        ! kappa(1) has the largest scale. Max-heap extract leaves ascending
        ! order; the final swap restores decreasing scale.
        SUBROUTINE sort_kappa_desc()
            INTEGER :: i, start, bottom, tmp
            IF (nft <= 1) RETURN
            DO start = nft/2, 1, -1
                CALL siftdown_kappa(start, nft)
            END DO
            DO bottom = nft, 2, -1
                tmp            = kappa(1)
                kappa(1)       = kappa(bottom)
                kappa(bottom)  = tmp
                CALL siftdown_kappa(1, bottom-1)
            END DO
            DO i=1,nft/2
                tmp               = kappa(i)
                kappa(i)          = kappa(nft+1-i)
                kappa(nft+1-i)    = tmp
            END DO
        END SUBROUTINE sort_kappa_desc

        SUBROUTINE siftdown_kappa(start, bottom)
            INTEGER, INTENT(IN) :: start, bottom
            INTEGER :: root, child, tmp
            root = start
            DO
                child = 2*root
                IF (child > bottom) EXIT
                IF (child+1 <= bottom) THEN
                    IF (scale(kappa(child+1)) > scale(kappa(child))) THEN
                        child = child + 1
                    END IF
                END IF
                IF (scale(kappa(root)) >= scale(kappa(child))) EXIT
                tmp            = kappa(root)
                kappa(root)    = kappa(child)
                kappa(child)   = tmp
                root           = child
            END DO
        END SUBROUTINE siftdown_kappa

    END SUBROUTINE improvement_potential_acc


    ! In-place heapsort: dminim(nk(1,col)) >= ... >= dminim(nk(nel(col),col)).
    SUBROUTINE sort_cluster_by_dminim_desc(col)

        USE initclust, ONLY : dminim, nk, nel

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: col

        INTEGER :: n, i, start, bottom, tmp

        n = nel(col)
        IF (n <= 1) RETURN

        DO start = n/2, 1, -1
            CALL siftdown_cl(start, n)
        END DO
        DO bottom = n, 2, -1
            tmp           = nk(1,col)
            nk(1,col)     = nk(bottom,col)
            nk(bottom,col)= tmp
            CALL siftdown_cl(1, bottom-1)
        END DO
        DO i=1,n/2
            tmp              = nk(i,col)
            nk(i,col)        = nk(n+1-i,col)
            nk(n+1-i,col)    = tmp
        END DO

    CONTAINS

        SUBROUTINE siftdown_cl(start, bottom)
            INTEGER, INTENT(IN) :: start, bottom
            INTEGER :: root, child, tmp
            root = start
            DO
                child = 2*root
                IF (child > bottom) EXIT
                IF (child+1 <= bottom) THEN
                    IF (dminim(nk(child+1,col)) > dminim(nk(child,col))) THEN
                        child = child + 1
                    END IF
                END IF
                IF (dminim(nk(root,col)) >= dminim(nk(child,col))) EXIT
                tmp              = nk(root,col)
                nk(root,col)     = nk(child,col)
                nk(child,col)    = tmp
                root             = child
            END DO
        END SUBROUTINE siftdown_cl

    END SUBROUTINE sort_cluster_by_dminim_desc

END MODULE clusteringmod
