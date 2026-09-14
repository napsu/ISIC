!*************************************************************************
!*                                                                       *
!*     ISIC - Stochastic Nonsmooth Optimization based                    *
!*     Incremental Clustering Software                                   *
!*                                                                       *
!*     ISIC is designed for solving minimum sum-of-squares clustering    *
!*     problems in very large and big datasets. It uses a greedy         *
!*     starting-point algorithm and stochastic inexactLMBM as its        *
!*     optimization engine.                                              *
!*                                                                       *
!*     Napsu Karmitsa 2026 (last modified 14.9.2026).                    *
!*                                                                       *
!*     This work was part of the PRIVIDEMA project funded by the         *
!*     European Union under Grant Agreement No. 101167964.               *
!*                                                                       *
!*     ISIC software is covered by the MIT license.                      *
!*                                                                       *
!*************************************************************************
!*
!*     Initialization of parameters for ISIC                              
!*
!*     Modules included: 
!* 
!*     initclust        ! Initialization of parameters for clustering.
!*     initislmb        ! Initialization of inexact SLMB -solver.
!*

MODULE initclust  ! Initialization of parameters for clustering codes.

    USE r_precision, ONLY : prec  ! Precision for reals.
    IMPLICIT NONE

    ! Names of input and output files:
    CHARACTER(LEN=80), SAVE :: &
        infile = '../data/Skin_NonSkin.txt',&     ! Input data 
        outfile0 = 'isic_centers.txt', &          ! Result file with cluster centers.
        outfile1 = 'isic_results.txt'             ! Result file with function values,
                                                  ! Davies-Bouldin (DB) validity index,
                                                  ! and cpu-time.
        
    ! Input integer parameters. Give values here.
    INTEGER, PARAMETER :: &
        nclust = 25, &          ! Maximum number of clusters, from user.
        nft = 3, &                  ! Number of features in data, from user.
        nrecord = 245057, &         ! Number of records in data, from user.
        !
        nrec_used = nrecord, &      ! Number of data points used at once in computations, nrec_used <= nrecord.
                                    !   E.g. nrec_used = nrecord if nrecord < 100000, nrec_used = MAX(nrecord/50,10000), otherwise.        
        !nrec_used = int(MAX(nrecord/10.0,10000.0)), &   ! 
        ibig     = 0, &             ! parameter for the size of data, ibig=1 => no computations with whole data will be made
                                    ! note that this version 0.1 still reads the full data.
        ipc      = 0, &             ! Printout parameter for cluster centers (0 = no printout).
        iaux     = 0, &             ! Auxiliary problem:
                                    !   0  - skip aux; use the new center directly in ns=2 (default),
                                    !   1  - solve aux clustering on the new center.
        maxdim   = nclust*nft, &    ! Maximum number of variables in optimization (do not change).
        maxsize  = 200*maxdim       ! Space for the candidate points,
                                    !   maxsize > m1*maxdim (usually equality is enough).


    ! Input real parameters. Give values here.
    REAL, SAVE :: &
        tlimit    = 7200.0          ! Maximum CPU-time in seconds, from user.



    ! Other Real parameters.
    REAL(KIND=prec), SAVE :: &      !
        a(nft,nrecord), &           ! Data matrix, from input file, give the name of the file above.
        a_used(nft,nrec_used), &    ! Partial data matrix.
        tnorm, &                    ! The total number of norms computed.
        xbest(maxdim), &            ! Best solution obtained.
        dcent(nclust,nclust), &     ! Distance (affinity) matrix for cluster centers
        dminim(nrec_used), &        ! Distance of a(i) and the nearest centroid.
        f_prev_nc                   ! Objective value from the previous nc (nc-1 clusters)

    ! Other integer parameters.
    INTEGER, SAVE :: &              !
        nc, &                       ! Current number of clusters, loops from 1 to nclust
        ns, &                       ! Switch for auxiliary and real clustering problem.
        m, &                        ! Number of variables in optimization:
                                    !   m = nft    if ns = 1,
                                    !   m = nft*nc if ns = 2.
        m1, &                       ! Maximum number of datapoints selected as initial 
                                    !   solutions, m1 <= nrecord. Initialized in def_clusterpar below.
        ncand, &                    ! Number of canditate points.
        list1(nrec_used), &         ! list1(i) gives the cluster where point i belongs.
        nel(nclust), &              ! nel(i) is the number of data points in cluster i.
        nk(nrec_used,nclust), &     ! nk(j,i) (j=1,...,nel(i)) tells which data 
                                    !   points are in cluster i.
        l4(nrec_used)               ! Auxiliary array

    ! Allocatable tables
    INTEGER, SAVE, DIMENSION(:), allocatable :: &
        lcand                       ! lcand(ncand), Data points which are selected as initial cluster centers.


CONTAINS

    SUBROUTINE def_clustpar()       ! Default values for parameters.

        IMPLICIT NONE
        INTEGER :: size

        size = nrecord*nft
        IF (size < 1000000) THEN
            m1=min(nrecord,500)  ! Maximum number of data points selected as initial solutions.                      
                                 ! m1 <= nrecord, Note that m1 is adaptive -> 200.
        ELSE IF (size <= 10000000) THEN
            m1=500
        ELSE 
            m1=200 
        END IF

    END SUBROUTINE def_clustpar

END MODULE initclust


MODULE initislmb  ! Initialization of parameters for inexact SLMB.

    USE r_precision, ONLY : prec   ! Precision for reals.
    USE param, ONLY : zero, one, small   ! Parameters.
    USE initclust, ONLY : &
        maxdim, &                  ! Maximum dimension of x.
        nrec_used                  ! Number of used records in data.
    IMPLICIT NONE

    ! Parameters
    INTEGER, PARAMETER :: &
        na      = 10, &            ! Maximum bundle dimension, na >= 2. HIHU, onko tämä järkevä stokastisessa?
        mcu     = 7, &             ! Maximum number of stored corrections, mcu >=1.
        mcinit  = 7, &             ! Initial maximum number of stored corrections, mcu >= mcinit >= 3.
                                   ! If mcinit <= 0, the default value mcinit = 3 will be used.
                                   ! However, the value mcinit = 7 is recommented.
        inma    = 1, &             ! Selection of optimization method:
                                   !   inma = 0, InexactLMBM with simple step size selection,
                                   !   inma = 1, InexactLMBM with nonmonotone step size selection,
        mnma    = 5                ! Maximum number of function values used in nonmonotone line search.
        

    ! Real parameters (if parameter value <= 0.0 the default value of the parameter will be used).
    REAL(KIND=prec), SAVE :: &
        tolb    = -small, &        ! Tolerance for the function value, 
                                   !   if tolb == 0 the default value -large will be used.
        tolf    = 1.0E-08_prec, &  ! Tolerance for change of function values (default = 1.0E-8).
        tolf2   = 1.0E-06_prec, &  ! Tolerance for change of gradient norms in null steps.
                                   !   - If tolf2 < 0 the the parameter and the corresponding termination 
                                   !   criterion will be ignored.
                                   !   - If tolf2 = 0 the default value 1.0E-6 will be used 
                                   !     (with inma=0, 1.0E-08 is recommented). 
        tolg    = 1.0_prec, &      ! Tolerance for the termination criterion: 
                                   !   default = 1.0E-5, for clustering larger values are recommented.
        tolg2   = 1.0E+3_prec, &   ! Tolerance for the second termination criterion:
                                   !   default = 1.0E-3, for clustering larger values are recommended.
        eta     = 0.50_prec, &     ! Distance measure parameter, eta > 0.
                                   !   - If eta < 0  the default value 0.5 will be used.
        epsl    = 0.24_prec        ! Decrease parameter, 0 < epsl < 0.25 (default = 0.24);

    ! Integer parameters (if value <= 0 the default value of the parameter will be used).
    INTEGER, SAVE :: & 
        n       = maxdim,  &      ! Number of variables (do not change)
        nbatch  = int(max(SQRT(real(nrec_used))*10.0,1000.0)), & ! Size of the batch 1 <= nbatch <= nrecords
                                  ! The √nrecord scaling is motivated by the fact that gradient variance scales as 1/nbatch, 
                                  ! so to keep the relative noise level constant as the dataset grows you need nbatch ∝ nrecord. 
                                  ! The square root is a compromise that avoids making batches prohibitively large for very large datasets.
       !nbatch  = int(max(nrec_used/20.0,1000.0)), & ! Size of the batch 1 <= nbatch <= nrecords
        maxbatch = 100, &         ! Number of different batches (change in subroutine below)
        nth     = 10, &           ! Number of iterations after which a new batch is selected
                                  !    (change in subroutine below)
        batchtype = 0, &          ! Type for selecting the batches:
                                  !    0  - Totally random batches
                                  !    1  - Random batches such that all indices are used at least once      
        mit     = 1000, &         ! Maximum number of iterations: mit = maxbatch*nth
        mfe     = 1000, &         ! Maximum number of function evaluations: mfe = maxbatch*nth
        mtesf   = MAX(mnma+5,10), & ! Maximum number of iterations with changes of
                                  !   function values smaller than tolf (default = 10)
        iprint  =    1, &         ! Printout specification (only for real clustering problem, 
                                  !   for aux clustering problem the printout is cancelled in islmb.f95,   
                                  !   select either 0 or 1):
                                  !    -1  - No printout
                                  !     0  - Only the error messages
                                  !     1  - The final values of the objective function
                                  !          (default used if iprint < -1)
                                  !     2  - The final values of the objective function and the
                                  !          most serious warning messages
                                  !     3  - The whole final solution
                                  !     4  - At each iteration values of the objective function
                                  !     5  - At each iteration the whole solution
        iscale  =    0            ! Selection of the scaling with SLMB:
                                  !     0  - Scaling at every iteration with STU/UTU (default)
                                  !     1  - Scaling at every iteration with STS/STU
                                  !     2  - Interval scaling with STU/UTU
                                  !     3  - Interval scaling with STS/STU
                                  !     4  - Preliminary scaling with STU/UTU
                                  !     5  - Preliminary scaling with STS/STU
                                  !     6  - No scaling


    INTEGER, SAVE :: ib                             ! Index for batch 
    REAL(KIND=prec), DIMENSION(maxdim), SAVE :: x   ! Vector of variables
    INTEGER :: i
    INTEGER, DIMENSION(nrec_used), SAVE :: batchind ! Batch indices 


CONTAINS

    SUBROUTINE defaults()  ! Default values for parameters.

        USE param, ONLY: small, large, zero, one, half
        IMPLICIT NONE

        IF (nbatch <= 0) nbatch  = int(max(nrec_used/50.0,1000.0)) ! Size of the batch.
        IF (iprint < -1) iprint  = 1               ! Printout specification.
        IF (mit   <= 0) mit      = 5000            ! Maximum number of iterations.
        IF (mfe   <= 0) mfe      = n*mit           ! Maximum number of function evaluations.
        IF (tolf  <= zero) tolf  = 1.0E-08_prec    ! Tolerance for change of function values.
        IF (tolf2 == zero) tolf2 = 1.0E-06_prec    ! Second tolerance for change of function values.
        IF (tolb  == zero) tolb  = -large + small  ! Tolerance for the function value.
        IF (tolg  <= zero) tolg  = 1.0E-05_prec    ! Tolerance for the termination criterion.
        IF (tolg2 <= zero) tolg2 = 1.0E-03_prec    ! Tolerance for the second termination criterion.
        IF (eta   <  zero) eta   = 0.5_prec     ! Distance measure parameter
        IF (epsl  <= zero) epsl  = 0.24_prec       ! Line search parameter,
        IF (mtesf <= 0) mtesf    = 10              ! Maximum number of iterations with changes
                                                   ! of function values smaller than tolf.
        IF (iscale > 6 .OR. iscale < 0) iscale = 0 ! Selection of the scaling.

    END SUBROUTINE defaults

    SUBROUTINE init_slmbpar()  ! Subroutine for further initialization of parameters.

        USE initclust, ONLY : ns, nc, nrec_used
        IMPLICIT NONE

        IF (ns == 1) THEN
            maxbatch = 1        ! Number of different batches (aux. problem)
            nth     = 20        ! Number of iterations after which a new batch is selected
        ELSE 
            maxbatch = MAX(20, INT(0.15 * SQRT(REAL(nrec_used)))) ! Number of different batches (clustering problem)
            nth     = 10        ! Number of iterations after which a new batch is selected
        END IF
                
    END SUBROUTINE init_slmbpar

END MODULE initislmb