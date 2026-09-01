!*************************************************************************
!*                                                                       *
!*     ISIC - Stochastic Nonsmooth Optimization based                    *
!*     Incremental Clustering Software (bigClust version 2.0)            *
!*                                                                       *
!*     ISIC differs from original bigClust by using                      *
!*     stochastic inexactLMBM as an optimization engine and              *
!*     greedy starting point selection procedure.                         * 
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
!*
!*     Codes included:
!*
!*     isic.f95              - Mainprogram for clustering software
!*                             (this file).
!*     parameters.f95        - Parameters. Inludes modules:
!*                               - r_precision - Precision for reals,
!*                               - param - Parameters,
!*                               - exe_time - Execution time.
!*     initisic.f95          - initialization of clustering parameters and inexactLMBM.
!*                             Includes modules:
!*                               - initclust - Initialization of parameters for clustering.
!*                               - initislmb - Initialization of inexact SLMB.
!*     clusteringmod.f95     - Subroutines for clustering software.
!*     islmb.f95             - iSLMB - Stochastic inexact limited memory bundle method.
!*     objfun.f95            - Computation of the cluster function and subgradients values.
!*     subpro.f95            - subprograms for iSLMB.
!*
!*     Makefile              - makefile.
!*
!*
!*     To use the software modify initisic.f95 as needed.
!*
!*
!*     References:
!*
!*     for ISIC:
!*       N. Karmitsa, J. Lampainen, K. Joki, M.M. Mäkelä, "ISIC — Incremental Clustering Software 
!*       Based on Stochastic Inexact Nonsmooth Optimization", 2026
!*
!*     for BigClust:
!*       N. Karmitsa, V.-P. Eronen, M.M. Mäkelä, A. Airola, T. Pahikkala, "Stochastic limited 
!*       memory bundle method for clustering in big data", Pattern Recognition, Vol. 165, 111654, 2025.
!*
!*     for inexactLMBM:
!*       J. Lampainen, K. Joki, N. Karmitsa, M.M. Mäkelä, "Inexact limited memory bundle method",
!*       submitted 2026, https://arxiv.org/abs/2604.08067
!*
!*     for original LMBM:
!*       N. Haarala, K. Miettinen, M.M. Mäkelä, "Globally Convergent Limited Memory Bundle Method  
!*       for Large-Scale Nonsmooth Optimization", Mathematical Programming, Vol. 109, No. 1,
!*       pp. 181-205, 2007. DOI 10.1007/s10107-006-0728-2.
!*
!*       M. Haarala, K. Miettinen, M.M. Mäkelä, "New Limited Memory Bundle Method for Large-Scale 
!*       Nonsmooth Optimization", Optimization Methods and Software, Vol. 19, No. 6, pp. 673-692, 2004. 
!*       DOI 10.1080/10556780410001689225.
!*
!*     for clustering:
!*       A. Bagirov, N. Karmitsa, S Taheri, "Partitional Clustering via Nonsmooth Optimization",
!*       Springer, 2020 (2025 2nd edition).
!*
!*     for NSO:
!*       A. Bagirov, N. Karmitsa, M.M. Mäkelä, "Introduction to nonsmooth optimization: theory, 
!*       practice and software", Springer, 2014.
!*
!*
!******************************************************************************
!*
!*     * PROGRAM isic *
!*
!*     Main program for stochastic nonsmooth clustering software with iSLMB.
!*
!******************************************************************************

PROGRAM isic 

    USE r_precision, ONLY : prec        ! Precision for reals.
    USE param, ONLY : zero, one, large  ! Parameters.

    USE initclust, ONLY : &             ! Initialization of clustering parameters.
        infile, &                       ! Dataset file.
        outfile0, &                     ! Result file with cluster centers.
        outfile1, &                     ! Result file with function values,
                                        !   validity indices, and cpu-time.
        ibig, &                         ! Switch for the size of data, if ibig==1, no computation with whole data are made.
        ipc, &                          ! Printout parameter for cluster centers.
        maxsize, &                      ! Maximum number of candidate points of data set.
        a, &                            ! Data matrix a(nft,nrecord), from input file.
        a_used, &                       ! Partial data matrix a(nft,nrec_used), from input file.
        xbest, &                        ! Best solution obtained
        f_prev_nc, &                    ! Objective value from the previous nc (nc-1 clusters) 
        dcent, &                        ! Distance (affinity) matrix for cluster centers
        tlimit, &                       ! Time limit, from user.
        tnorm, &                        ! The total number of norms computed.
        nclust, &                       ! Maximum number of clusters, from user.
        nft, &                          ! Number of features in data, from user.
        nrecord, &                      ! Number of instances in data, from user.
        nrec_used, &                    ! Number of data points used at the time.
        nc, &                           ! Current number of clusters, loops from 1 to nclust.
        ns, &                           ! Switch for auxiliary and real clustering problem.
        m, &                            ! Number of variables in optimization:
                                        !   m = nft    if ns = 1,
                                        !   m = nft*nc if ns = 2.
        m1, &                           ! Maximum number of initial solutions.
        lcand, &
        def_clustpar                    ! Default values of clustering parameters.
    USE initislmb, ONLY : &             ! 
        nbatch, &                       ! Size of the batch 1 <= nbatch <= nrecords
        maxbatch, &                     ! Maximum number of different batches
        nth, &                          ! Number of iterations after which a new batch is selected
        batchtype, &                    ! Type for selecting the batches:
                                        !     0  - Totally random batches
                                        !     1  - Random batches such that all indices are used at least once      
        batchind
    USE obj_fun, ONLY : clustfull,clustfull2      ! Computation of full cluster function value
    USE clusteringmod                   ! Subprograms for clustering.
    USE islmb_mod, ONLY : islmb_init      ! SLMB method for optimization.
    USE exe_time, ONLY : getime         ! Execution time.

    IMPLICIT NONE


    REAL(KIND=prec), DIMENSION(:), allocatable :: &
        x, x2
    REAL(KIND=prec) :: &
        barf, &
        db,db2,dn,dn3,sep, & !validity indices
        f, &
        fold, &
        fclust1, &
        fcurrent 
    REAL :: &
        time1, &
        time3, &
        time4, &
        inittime, inittime2, &
        ns1time, &
        ns2time
    INTEGER :: &
        i,j,k,k1,j1, &
        nseed
    REAL :: urand    ! Random number
    !character(len=10):: junk  ! Needed e.g. with cord_19_embeddings_2022-06-02.csv
    REAL :: cputotal,usertimetotal

    INTEGER, allocatable :: seed(:)
    CHARACTER(LEN=20) :: seed_arg
    INTEGER :: seed_val

    !==================================================================
    !   Initialization
    !==================================================================
    call random_seed(size = nseed)
    allocate(seed(nseed))

    seed_val = 12345                              ! default seed
    IF (command_argument_count() >= 1) THEN
        CALL get_command_argument(1, seed_arg)
        READ(seed_arg, *) seed_val
    END IF
    seed = seed_val
    call random_seed(put=seed)

    CALL def_clustpar()
    do i=1,nrec_used
        batchind(i)=i
    end do

    allocate(x(nclust*nft),lcand(nrec_used))
    allocate(x2(maxsize))


    OPEN(40,file=outfile0) ! Cluster centers
    OPEN(42,file=outfile1) ! Summary of results
    OPEN(78,file=infile,status='old',form='formatted')


    !==================================================================
    !   Initializing output
    !==================================================================
    WRITE(40, *) 'iSIC: Cluster centers'
    WRITE(40, *)
    WRITE(40, *) 'This file contains the final cluster centers with k =',nclust,' clusters and'
    WRITE(40, *) 'the solutions to all intermediate l-clustering problems with l = 1,...,k-1. '
    if (ipc == 0) then
        WRITE(40, *) 
        WRITE(40, *) 'Writing down the cluster centers is canceled. Select ipc = 1 in initisic to enable them.'
    end if
    WRITE(40, *)
    WRITE(42, *) 'Summary of results of iSIC:'
    WRITE(42, *)
    WRITE(42, *) 'Data: ',infile
    WRITE(42, *) 'No. data points:         ',nrecord
    WRITE(42, *) 'No. used data points:    ',nrec_used
    WRITE(42, *) 'No. features:            ',nft
    WRITE(42, *) 'Size of the batch:       ',nbatch
    WRITE(42, *) 'Max no. batches:         ',maxbatch
    WRITE(42, *) 'No. iterations after which a new batch is selected: ', nth
    WRITE(42, *) 'Type of the batch:       ',batchtype
                !     0  - Totally random batches
                !     1  - Random batches such that all indices are used at least once        
    WRITE(42, *) 'Type of data:            ',ibig
                !     0  - Large data
                !     1  - Big data (no computations with full data)
    WRITE(42, *) 'Max. canditate starting points: ',m1 
    WRITE(42, *) 'Seed for random numbers: ',seed(1)
    
    WRITE(42, *) 'Optimization with inexact SLMB.'
    WRITE(42, *)
    WRITE(42, *) 'nclust |', ' f                      |',' D-B       |',' Dunn      |',&
    &' Norms       |', ' CPU-time |', ' Init. time |', ' Valid. time |',&
    &' Aux.Clust. time |', ' Clust. time |'


    !==================================================================
    !   Reading the data
    !==================================================================
    
    tnorm = zero

    !if (ibig == 0) then
        CALL getime(cputotal) 
        !READ(78,*) ! read the first line (needed e.g. with Range-Queries-Aggregates.csv)
        DO i=1,nrecord 
            !READ(78,'(A)', ADVANCE='NO') junk ! read the first char (needed e.g. with cord_19_embeddings_2022-06-02.csv)
            READ(78,*) (a(k,i),k=1,nft)
        END DO
        CALL getime(usertimetotal) 
        print*,'Data read.',usertimetotal-cputotal
    !else !(ibig == 1)   
        ! We always read the whole data
    !end if 
    
    CALL getime(time1) 

    ! Read random partial data matrix
    if (nrecord /= nrec_used) then
        do i = 1, nrec_used
            do k = 1, nft
                a_used(k,i) = a(k,i)
            end do
        end do
        do i = nrec_used + 1, nrecord
            call random_number(urand)
            j = 1 + INT(REAL(i, prec) * urand)
            if (j <= nrec_used) then
                do k = 1, nft
                    a_used(k,j) = a(k,i)
                end do
            end if
        end do
    else if (ibig == 0) then
        do i=1,nrec_used 
            do k=1,nft
                a_used(k,i) = a(k,i)
            end do
        end do
    else ! ibig == 1
        STOP 'Error: Big data with nrecord == nrec_used.'
    end if    

    
    !================================================================
    !   Computing cluster centers incrementally
    !================================================================
    outerloop: DO nc=1,nclust  ! number of clusters


        PRINT 42,nc
42      FORMAT('Cluster No.:',i10)

        IF(nc > 1) THEN
            if (nrecord*nft > 10000000) then
                m1  = max(m1-50,100)
            else if (nrecord*nft > 1000000) then
                m1  = max(m1-50,200)
            end if
            f_prev_nc = f !
            CALL getime(inittime) 
            CALL step2(x2) !  Step2 computes clusters for each data point and gives inital points
            CALL getime(time4) 
            inittime = time4 - inittime
            ns=1 ! Auxiliary clustering problem

            CALL islmb_init(x2(1:nft),barf,time4) ! Call for inexact SLMB
            ns1time=time4
            
            m=nft*nc
            ns2time=0.0

            DO i=1,nft
                x(i+(nc-1)*nft)=x2(i)
            END DO
            ns=2
            CALL islmb_init(x,fcurrent,time4) ! Call for inexact SLMB
            ! Returns fcurrent in batch.

            ns2time=time4
            CALL clustfull(x, f) !

            DO j1=1,m
                xbest(j1)=x(j1)
            END DO

            
            !================================================================
            !    Write cluster centers to outfile
            !================================================================
            WRITE(40, *) '____________________________________________________'
            WRITE(40, *)
            WRITE(40, *) '    Number of clusters: ',nc
            WRITE(40, *) '____________________________________________________'
            WRITE(40, *)
            if (ipc == 1) then        
                WRITE(40, *)
    
                DO j=1,nc
                    WRITE(40, *) 'Center of cluster No.',j
                    WRITE(40, *)
                    WRITE(40,49) (x(i+(j-1)*nft),i=1,nft)
                END DO
49              FORMAT(5f16.8)
                WRITE(40, *)
                
            end if

            !================================================================

        ELSE  ! nc=1, Step 0 in SLMB-Clust algorithm
            ns2time=0.0
            CALL step1(f,x) ! computes the centroid and the value of clustering function at the centroid
            fclust1=f 
            DO j1=1,nft ! xbest for aux-problem
                xbest(j1)=x(j1)
            END DO
            

            !================================================================
            !    Write the centroid to outfile
            !================================================================
            WRITE(40, *)
            WRITE(40, *) '____________________________________________________'
            WRITE(40, *)
            WRITE(40, *) '    Number of clusters: ',nc
            WRITE(40, *) '____________________________________________________'
            WRITE(40, *)

            if (ipc == 1) then
                WRITE(40, *) 'Center of cluster No.',nc
                WRITE(40, *)
                WRITE(40,48) (x(i),i=1,nft)
48              FORMAT(5f16.8)
                WRITE(40, *)
            end if
            !================================================================

        END IF

        ! Read random partial data matrix
        if (nrecord /= nrec_used) then
            do i = 1, nrec_used
                do k = 1, nft
                    a_used(k,i) = a(k,i)
                end do
            end do
            do i = nrec_used + 1, nrecord
                call random_number(urand)
                j = 1 + INT(REAL(i, prec) * urand)
                if (j <= nrec_used) then
                    do k = 1, nft
                        a_used(k,j) = a(k,i)
                    end do
                end if
            end do
        end if
    
        DO k=1,nc
            dcent(k,k)=zero
        END DO
        DO k=1,nc
            DO k1=k+1,nc
                dcent(k,k1)=zero
                DO j=1,nft
                    dcent(k,k1)=dcent(k,k1)+(x(j+(k-1)*nft)-x(j+(k1-1)*nft))**2
                END DO
                dcent(k1,k)=dcent(k,k1)
            END DO
        END DO

        ! Subsample resampled between optimization and evaluation
        CALL getime(inittime2) 
        CALL check(x,f,db,db2,dn,dn3,sep) ! checks the results and computes validity indices
        CALL getime(time4) 
        
        inittime2 = time4 - inittime2

        !if (ibig == 0) then ! We do not compute this with big data
        fold = f    
        if (nrecord /= nrec_used) then ! or if entire data set is used
            call clustfull2(x,fold)
        end if
        

        CALL getime(time3)

        WRITE(40,*)
        WRITE(40,543) f
        WRITE(40,*)
        WRITE(40,142) INT(tnorm)
        time4=time3-time1
        if (time4 > tlimit) then
            print*,' Warning: time limit exceed. Returned final solution has only ',nc,' clusters.'
            EXIT outerloop
        end if
        WRITE(40,*)
        WRITE(40,141) time4
        WRITE(40,*)
        WRITE(40,*)
        WRITE(42,603) nc,fold,db,dn,tnorm,time4,inittime,inittime2,ns1time,ns2time
    END DO outerloop
    WRITE(42, *)
    
    WRITE(42, *) 'Data read in ',usertimetotal-cputotal,'seconds.'
    WRITE(42, *)
    WRITE(42, *)

    
    
141 FORMAT('CPU time:                                         ',f12.3)
142 FORMAT('The total number of norms:                  ',i18)
543 FORMAT('The value of the cluster function:',f28.6)
603 FORMAT(i4,f38.4,f12.4,f12.4,f15.0,f10.2,f10.2,f10.2,f10.2,f10.2)
    
    CLOSE(40)
    CLOSE(42)
    CLOSE(78)
        
    deallocate(x,lcand)
    deallocate(x2)
    deallocate(seed)

    STOP

END PROGRAM isic
