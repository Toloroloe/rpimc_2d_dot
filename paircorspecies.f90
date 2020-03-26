module paircorspecies

  ! Pair correlation functions between each species as defined in the INPUT file.
  !
  ! This module also calculates the 2SEM error estimates, which does
  ! increase the simulation time... especially with large grids and many species.
  ! However, the inspection of results speeds up.

  use math
  use init
  implicit none
  integer, private, save  :: maarita_rajat    
  integer, private, save  :: histkoko
  integer(kind = isoluku), private, save :: nblocks
  integer(kind = isoluku), private, save :: nblock
  integer(kind = isoluku), private, save :: modpc
  integer(kind = isoluku), private, save :: dpc_ind
  integer(kind = isoluku), private, save :: dpc_indmax

  
 
  REAL(kind =rk), DIMENSION(:,:,:), allocatable, SAVE  :: pk_rajat
  REAL(kind = rk), DIMENSION(:,:), allocatable, SAVE  :: dr_pk
  LOGICAL, SAVE :: pcD

  INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE, private  :: whichpc

  integer, save, private :: NumPairCor
  type table_paircor
     character(len=10) :: str
     integer(kind = isoluku), DIMENSION(:), ALLOCATABLE    :: pc
     integer(kind = isoluku), DIMENSION(:), ALLOCATABLE    :: pcTOT
     real(kind=rk), DIMENSION(:), ALLOCATABLE    :: r
     real(kind=rk) :: dr
     integer(kind=isoluku) :: nCountsBlock
     integer(kind=isoluku) :: nCountsTOT

     real(kind = rk) :: Rave        ! <R>
     real(kind = rk) :: R2ave       ! <R**2>
     real(kind = rk) :: Rave1per    ! <1/R>
     
  end type table_paircor

  real(kind = rk) :: RaveMPI        ! <R>
  real(kind = rk) :: R2aveMPI       ! <R**2>
  real(kind = rk) :: Rave1perMPI    ! <1/R>   
  
  type(table_paircor), dimension(:), allocatable, private :: pcdata

  real(kind=rk), DIMENSION(:,:,:), ALLOCATABLE  :: dpc

  real(kind = rk), dimension(:), allocatable :: meanMPI, twoSEMMPI
  

contains
  
  subroutine initPairCorrelation()
    implicit none
    integer :: nblocks,i,NumPC
   
    maarita_rajat  = -1

    histkoko = PCGridSize
    nblocks = NumOfBlocks

    NumPC=getNumPairCorrelationFunctions(NumSpecies,NumSpeciesEQ1)
    NumPairCor=NumPC

    ! "frequency" in which the error estimate will be calculated
    dpc_indmax = 100
    modpc=nblocks/dpc_indmax
    if (modpc==0) then
       modpc=1
       allocate(dpc(nblocks,NumPC,histkoko))
    elseif (modpc==1 .and. nblocks>dpc_indmax) then
       write(*,*) 'WARNING: modulo(NumberOfBlocks,dpc_indmax) should be 0'
       allocate(dpc(dpc_indmax,NumPC,histkoko))
    else
       allocate(dpc(dpc_indmax,NumPC,histkoko))
    end if
    dpc_ind = 1
    
    !allocate(dpc(nblocks,NumPC,histkoko))
    dpc = 0.0_rk

    allocate(pcdata(NumPC))

    do i=1,NumPC
       ALLOCATE(pcdata(i)%pc(histkoko))
       ALLOCATE(pcdata(i)%pcTOT(histkoko))
       ALLOCATE(pcdata(i)%r(histkoko))

       pcdata(i)%pc = 0
       pcdata(i)%pcTOT = 0
       pcdata(i)%r = 0.0_rk
       pcdata(i)%str = ' '
       pcdata(i)%nCountsBlock = 0
       pcdata(i)%nCountsTOT = 0
       pcdata(i)%dr = 0.0_rk
       pcdata(i)%Rave     = 0.0_rk
       pcdata(i)%R2ave    = 0.0_rk
       pcdata(i)%Rave1per = 0.0_rk
    end do

    allocate(meanMPI(histkoko))
    allocate(twoSEMMPI(histkoko))

    allocate(whichpc(hi_lkm(3),hi_lkm(3)))
    call makePCtable(NumSpecies)

    ALLOCATE(pk_rajat(hi_lkm(3),hi_lkm(3),2))
    ALLOCATE(dr_pk(hi_lkm(3),hi_lkm(3)))

    nblock = 0
    pk_rajat   = 0.0_rk
    !pk_rajat(:,:,1) = 10000.0_rk    
    pk_rajat(:,:,1) = PCLowLimit
    !pk_rajat(:,:,2) = PCUppLimit
    pk_rajat(:,:,2) = 2.0_rk*sqrt(-2.0_rk/hiuk(1)%meff/sqrt(omega2)*log(5.0e-5_rk))
    dr_pk = 0.0_rk
    
    pcD = .FALSE.

  end subroutine initPairCorrelation

  subroutine BlockInitPairCorrelation
    implicit none
    integer :: i

    do i=1,NumPairCor
       pcdata(i)%pc = 0
       pcdata(i)%nCountsBlock = 0
       pcdata(i)%Rave     = 0.0_rk
       pcdata(i)%R2ave    = 0.0_rk
       pcdata(i)%Rave1per = 0.0_rk
    end do
    nblock   = nblock + 1

    if (nblock==block_cut) then
       do i=1,NumPairCor
          pcdata(i)%pcTOT = 0
          pcdata(i)%nCountsTOT = 0
       end do
    end if

  end subroutine BlockInitPairCorrelation

  subroutine endPairCorrelation()
    implicit none 

    deallocate(dpc)
    deallocate(whichpc)
    deallocate(pcdata)
    DEALLOCATE(pk_rajat)   
    DEALLOCATE(dr_pk)     

  end subroutine endPairCorrelation

    
  SUBROUTINE calcPairCorrelation(sgn_1,raja)
    use rtable
    use refbeadmodule
    IMPLICIT NONE
    
    integer(kind = isoluku), intent(in) :: sgn_1
    INTEGER(kind=isoluku) :: h, i, j, k, alku, loppu,raja,ModTrot
    integer               :: hkoko, indeksi
    REAL(kind =rk)        :: dr,dr2,oneperdr,up,down
    logical, save         :: NotBeenThere = .true.
    REAL(kind=rk)         :: rksgn
    integer               :: whichpc2
    integer               :: n_x,n_y,n_z
    integer(kind = isoluku), dimension(2) :: refbeads

    alku = 1
    loppu = hi_lkm(3)

    rksgn = real(sgn_1,kind=rk)

    if (exchange_on) then
       call getRefBead(refbeads(1),refbeads(2))
    end if
    
    if (raja>maarita_rajat) then
              
       hkoko = histkoko
       up = 1.0_rk
       down = 1.0_rk

       DO h = alku,loppu-1
          DO i = h+1,loppu
             if (NotBeenThere) then
                if (pk_rajat(h,i,1)<0.5) then
                   dr_pk(h,i) = up*pk_rajat(h,i,2)/hkoko
                   pk_rajat(h,i,1) = 0.0_rk
                else                   
                   dr_pk(h,i) = (up*pk_rajat(h,i,2)-down*pk_rajat(h,i,1))/hkoko
                end if
                if (h==loppu-1 .and. i==loppu) then
                   NotBeenThere = .false.
                end if
             end if

             whichpc2 = whichpc(h,i)
             

             !if (exchange_on) then
             !   ModTrot=2
             !else
             ModTrot=hiuk(h)%Trotter
             !end if
             DO k = 1,ModTrot
                
                !if (exchange_on) then
                !   dr=evaldr(refbeads(k),i,h)
                !   oneperdr=evaloneperdr(refbeads(k),i,h)
                !   dr2=dr**2
                !else                   
                dr=evaldr(k,i,h)
                oneperdr=evaloneperdr(k,i,h)
                dr2=dr**2
                !end if
                
                pcdata(whichpc2)%R2ave = pcdata(whichpc2)%R2ave + rksgn*dr2
                pcdata(whichpc2)%Rave = pcdata(whichpc2)%Rave + rksgn*dr
                pcdata(whichpc2)%Rave1per = pcdata(whichpc2)%Rave1per + rksgn*oneperdr
                
                indeksi = FLOOR((dr-down*pk_rajat(h,i,1))/dr_pk(h,i)) + 1
                
                !IF (indeksi <= histkoko .and. indeksi>0) THEN
                if (index_check(indeksi,1,histkoko)) then
                   pcdata(whichpc2)%pc(indeksi) = pcdata(whichpc2)%pc(indeksi) + 1
                   pcdata(whichpc2)%pcTOT(indeksi) = pcdata(whichpc2)%pcTOT(indeksi) + 1
                end if
                !Ntot(h,i) = Ntot(h,i) + sgn_1
                pcdata(whichpc2)%nCountsBlock = pcdata(whichpc2)%nCountsBlock + 1
                pcdata(whichpc2)%nCountsTOT = pcdata(whichpc2)%nCountsTOT + 1
                
             END DO
          
             IF (.NOT. pcD) THEN
                DO j = 1,histkoko
                   pcdata(whichpc2)%r(j) = down*pk_rajat(h,i,1) &
                        + REAL(j-1,kind=rk)*dr_pk(h,i)
                END DO
                pcdata(whichpc2)%dr=dr_pk(h,i)
             END IF
             
          END DO
       END DO
       
       pcD = .TRUE.

    else
       
       DO h = alku,loppu-1
          DO i = h+1,loppu
             
             whichpc2 = whichpc(h,i)

             !do k = 1,2
             DO k = 1, hiuk(h)%Trotter                
                dr=evaldr(k,i,h)
                oneperdr=evaloneperdr(k,i,h)
                dr2=dr**2                

                pcdata(whichpc2)%R2ave = pcdata(whichpc2)%R2ave + rksgn*dr2
                pcdata(whichpc2)%Rave = pcdata(whichpc2)%Rave + rksgn*dr
                pcdata(whichpc2)%Rave1per = pcdata(whichpc2)%Rave1per + rksgn*oneperdr

                pcdata(whichpc2)%nCountsBlock = pcdata(whichpc2)%nCountsBlock + 1

                if (dr>pk_rajat(h,i,2)) then
                   pk_rajat(h,i,2) = dr
                end if

                if (dr<pk_rajat(h,i,1)) then
                   pk_rajat(h,i,1) = dr
                end if

             END DO

          END DO
       END DO

    end if
    
  END SUBROUTINE calcPairCorrelation

  subroutine printPairCorrelation(sgn_2)
    implicit none
    integer(kind = isoluku), intent(in) :: sgn_2
    integer(kind = isoluku) :: tt
    integer :: i,j
    CHARACTER( LEN = 5 )  :: blockChar
    real(kind = rk) :: Volume,mean_pc

    tt = sgn_2
    
    WRITE( blockChar, '(I5)' ) nblock

    if (dpc_ind<=dpc_indmax) then
       call blockvaluesPairCorrelation
       if (modulo(nblock,modpc)==0) then
          dpc(dpc_ind,:,:) = dpc(dpc_ind,:,:)/modpc
          dpc_ind = dpc_ind + 1
       end if
    end if


    do i=1,NumPairCor
       OPEN( 10, FILE=TRIM(tiedAlku_e) // '.paircor.' &
            //  trim(pcdata(i)%str) & 
            // '.dat',STATUS = 'REPLACE')
!       OPEN( 101, FILE= 'pc_data/' // TRIM(tiedAlku_e) // '.paircor.' // trim(pcdata(i)%str) & 
!                // '.b' // TRIM(ADJUSTL(blockChar)) // '.dat',STATUS = 'REPLACE')
       OPEN( 110, FILE=TRIM(tiedAlku_e) // &
            '.DistAve.' // trim(pcdata(i)%str) // '.dat', position = 'append')
       WRITE(110,'((E20.10),(E20.10),(E20.10))') &
            pcdata(i)%Rave/pcdata(i)%nCountsBlock, &
            pcdata(i)%R2ave/pcdata(i)%nCountsBlock, &
            pcdata(i)%Rave1per/pcdata(i)%nCountsBlock
       
       do j=1,histkoko
          !Volume = 4.0_rk/3*pi*((pcdata(i)%r(j)+pcdata(i)%dr)**3-pcdata(i)%r(j)**3)
          Volume = pi*((pcdata(i)%r(j)+pcdata(i)%dr)**2-pcdata(i)%r(j)**2)
          mean_pc = real(pcdata(i)%pcTOT(j),kind=rk)/pcdata(i)%nCountsTOT/Volume
          WRITE( 10, '((F15.6),(E15.6),(E15.6))') pcdata(i)%r(j)+pcdata(i)%dr/2, mean_pc, 2.0_rk*SEM_dpc(mean_pc,i,j)
!          WRITE( 101, '((F15.6),(E15.6))') pcdata(i)%r(j), real(pcdata(i)%pc(j),kind=rk)/pcdata(i)%nCountsBlock/Volume
       end do
       
       close(10)
!       close(101)
       close(110)
    end do

  end subroutine printPairCorrelation

  subroutine printPairCorrelationMPI(sgn_2)
    use mpi
    implicit none
    integer(kind = isoluku), intent(in) :: sgn_2
    integer(kind = isoluku) :: tt
    integer :: i,j
    CHARACTER( LEN = 5 )  :: blockChar
    real(kind = rk) :: Volume,mean_pc,twoSEM
    integer :: myid, numprocs, ierr
    
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)


    tt = sgn_2
    
    WRITE( blockChar, '(I5)' ) nblock

    if (dpc_ind<=dpc_indmax) then
       call blockvaluesPairCorrelation
       if (modulo(nblock,modpc)==0) then
          dpc(dpc_ind,:,:) = dpc(dpc_ind,:,:)/modpc
          dpc_ind = dpc_ind + 1
       end if
    end if


    do i=1,NumPairCor
       call MPI_Reduce(pcdata(i)%Rave,RaveMPI,1,MPI_DOUBLE_PRECISION &
            ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       call MPI_Reduce(pcdata(i)%R2ave,R2aveMPI,1,MPI_DOUBLE_PRECISION &
            ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       call MPI_Reduce(pcdata(i)%Rave1per,Rave1perMPI,1,MPI_DOUBLE_PRECISION &
            ,MPI_SUM,0,MPI_COMM_WORLD,ierr)

       if (myid .eq. 0) then
          OPEN( 110, FILE=TRIM(tiedAlku_e) // &
               '.DistAve.' // trim(pcdata(i)%str) // '.dat', position = 'append')
          WRITE(110,'((E20.10),(E20.10),(E20.10))') &
               RaveMPI/pcdata(i)%nCountsBlock/numprocs, &
               R2aveMPI/pcdata(i)%nCountsBlock/numprocs, &
               Rave1perMPI/pcdata(i)%nCountsBlock/numprocs
          close(110)
       end if

       do j=1,histkoko
          !Volume = 4.0_rk/3*pi*((pcdata(i)%r(j)+pcdata(i)%dr)**3-pcdata(i)%r(j)**3)
          Volume = pi*((pcdata(i)%r(j)+pcdata(i)%dr)**2-pcdata(i)%r(j)**2)
          mean_pc = real(pcdata(i)%pcTOT(j),kind=rk)/pcdata(i)%nCountsTOT/Volume
          call MPI_Reduce(mean_pc,meanMPI(j),1,MPI_DOUBLE_PRECISION &
               ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          twoSEM=2.0_rk*SEM_dpc(mean_pc,i,j)
          call MPI_Reduce(twoSEM**2,twoSEMMPI(j),1,MPI_DOUBLE_PRECISION &
               ,MPI_SUM,0,MPI_COMM_WORLD,ierr)
       end do
       
       if (myid .eq. 0) then
          OPEN( 10, FILE=TRIM(tiedAlku_e) // '.paircor.' &
               //  trim(pcdata(i)%str) & 
               // '.dat',STATUS = 'REPLACE')
          do j=1,histkoko
             WRITE( 10, '((F15.6),(E15.6),(E15.6))') pcdata(i)%r(j)+pcdata(i)%dr/2, &
                  meanMPI(j)/numprocs, sqrt(twoSEMMPI(j))
          end do
          close(10)
       end if
       
     
    end do

  end subroutine printPairCorrelationMPI


  subroutine blockvaluesPairCorrelation
    implicit none
    integer :: i,j
    real(kind = rk) :: Volume

    do i=1,NumPairCor
       do j=1,histkoko
          
          !Volume = 4.0_rk/3*pi*((pcdata(i)%r(j)+pcdata(i)%dr)**3-pcdata(i)%r(j)**3)
          Volume = pi*((pcdata(i)%r(j)+pcdata(i)%dr)**2-pcdata(i)%r(j)**2)
          dpc(dpc_ind,i,j) = dpc(dpc_ind,i,j) + real(pcdata(i)%pc(j),kind=rk)/pcdata(i)%nCountsBlock/Volume
          
       end do
    end do
    
  end subroutine blockvaluesPairCorrelation



!   function mean_dpaircorr(ind,hi1,hi2) result(ka)
!     implicit none
!     integer, intent(in) :: ind, hi1,hi2
!     real(kind = rk) :: ka
!     integer :: i

!     ka = 0.0_rk
!     do i=1,nblock
!        ka = ka + dpaircorr(i)%pc(ind,hi1,hi2)
!     end do
!     ka = ka / nblock

!   end function mean_dpaircorr


  function SEM_dpc(ka_pc,ispecies,igrid) result(SEM)
    implicit none
    integer, intent(in) :: ispecies,igrid
    real(kind = rk), intent(in) :: ka_pc
    real(kind = rk) :: SEM,varE,kappa,coef
    integer :: i
    integer(kind = isoluku) :: nblock2

    if (modulo(nblock,modpc)==0) then        
       nblock2=minval((/nblock/modpc,dpc_indmax/))
    else
       nblock2=dpc_ind-1
    end if

    if (nblock2<=1) then
       SEM = 0.0_rk
    else
       
       varE = 0.0_rk
       do i=1,nblock2
          varE = varE + (dpc(i,ispecies,igrid) - ka_pc)**2
       end do
       varE = varE/(nblock2-1)
       
       if (varE<1.0e-40_rk) then
          SEM = 0.0_rk
       else
          kappa = 0.0_rk
          do i = 1,nblock2-1
             coef = Cfunc(i,nblock2,ka_pc,varE,ispecies,igrid)
             if (coef<0) then
                kappa = 1.0_rk + 2.0_rk*kappa
                SEM = sqrt(varE*kappa/nblock2)
                return
             else
                kappa = kappa + coef
             end if
          end do
          kappa = 1.0_rk + 2.0_rk*kappa
          
          SEM = sqrt(varE*kappa/nblock2)
       end if
    end if

  end function SEM_dpc

  function Cfunc(i,dim,meanE,varE,ispecies,igrid) result(coeff)
    implicit none
    integer, intent(in) :: ispecies,igrid,i
    integer(kind = isoluku), intent(in) :: dim
    !real(kind = rk), dimension(dim) :: E
    real(kind = rk), intent(in) :: meanE, varE
    real(kind = rk) :: coeff
    integer(kind = isoluku) :: j
    
    coeff = 0.0_rk
    do j = 1,dim-i
       coeff = coeff + (dpc(i,ispecies,igrid)-meanE)*&
            (dpc(j+1,ispecies,igrid)-meanE)
    end do
    coeff = coeff/(dim-i)/varE
    
  end function Cfunc


  function getNumPairCorrelationFunctions(nSpecies,nSpeciesEQ1) result(Nn)
    implicit none
    integer :: i, nSpecies,nSpeciesEQ1,Nn

    Nn=1
    if (nSpecies>1) then
       do i=2,nSpecies
          Nn = Nn + i
       end do
    end if
    Nn = Nn-nSpeciesEQ1
    
  end function getNumPairCorrelationFunctions

  subroutine makePCtable(NumPC)
    implicit none
    integer :: i,j,pcNum,NumPC
    integer, dimension(NumPC,NumPC) :: table
    character(len=10) :: str1,str2

    table=0
    pcNum = 1
    do i=1,hi_lkm(3)-1
       do j=i+1,hi_lkm(3)
          
          if (table(hiuk(i)%species,hiuk(j)%species)==0) then
             table(hiuk(i)%species,hiuk(j)%species)=pcNum
             pcNum = pcNum + 1
          end if
          
          whichpc(i,j) = table(hiuk(i)%species,hiuk(j)%species)
          whichpc(j,i) = table(hiuk(i)%species,hiuk(j)%species)
          
          WRITE( str1 , '(I5)'  ) hiuk(i)%species
          WRITE( str2 , '(I5)'  ) hiuk(j)%species

          pcdata(table(hiuk(i)%species,hiuk(j)%species))%str = &
               TRIM(ADJUSTL(str1)) // '.' // &
               TRIM(ADJUSTL(str2))
          
       end do
    end do

    
  end subroutine makePCtable


end module paircorspecies
