module observables
  
  use math
  use init
  use paircorspecies
  use ag_density
  implicit none
  
  
  INTEGER, PARAMETER                   :: NumObserv = 4
  REAL(KIND=RK), DIMENSION(NumObserv), save :: Observ
  CHARACTER(LEN=10), DIMENSION(NumObserv),SAVE :: obs_name
  INTEGER, PARAMETER           :: i_pot = 1    !  V
  INTEGER, PARAMETER           :: i_tot = 2    !  Etot
  INTEGER, PARAMETER           :: i_kin = 3    !  Ekin 
  INTEGER, PARAMETER           :: i_ee = 4    !  EextAG
  REAL(KIND=RK), DIMENSION(NumObserv), save :: ObservMPI

  real(kind=rk), save :: sgnMPI
  real(kind=rk), private, save :: VextAG

  integer, private, save :: nblock
 
  INTEGER(kind = isoluku), DIMENSION(2), private, save  :: mittausten_lkm  
  INTEGER(kind = isoluku), DIMENSION(3), SAVE   :: sgn

  integer(kind = isoluku),dimension(:),allocatable,save :: BlockAcceptedCount
  integer(kind = isoluku),dimension(:),allocatable,save :: ParticleMovedCount
    
  real(kind = rk), dimension(3), private, save :: sysCM

contains


  subroutine initObservables()
    implicit none    
    integer :: i,xi,yi    

    allocate(BlockAcceptedCount(hi_lkm(3)+1))
    allocate(ParticleMovedCount(hi_lkm(3)+1))
    
    BlockAcceptedCount = 0
    ParticleMovedCount = 0
    
    obs_name(i_pot) = 'E_pot'
    obs_name(i_tot) = 'E_tot'
    obs_name(i_kin) = 'E_kin'
    obs_name(i_ee) = 'E_ee'
    
    mittausten_lkm = 0
    sgn = 0

    if (PCGridSize>0) then
       call initPairCorrelation()
    end if

    call init_ag_densityMPI()

    nblock=0
    
  end subroutine initObservables  

  subroutine BlockInitObservables()
    implicit none

    BlockAcceptedCount = 0
    ParticleMovedCount = 0
    mittausten_lkm(1)  = 0
    sgn(2)             = 0
    if (PCGridSize>0) then
       call BlockInitPairCorrelation()
       call block_init_ag_density()
    end if
    
    nblock = nblock + 1
    Observ             = 0.0_rk

  end subroutine BlockInitObservables

  subroutine endObservables()
    implicit none

    deallocate(BlockAcceptedCount)
    deallocate(ParticleMovedCount)
    
    if (PCGridSize>0) then
       call endPairCorrelation()
       call end_ag_density()
    end if
    
  end subroutine endObservables
  
  subroutine calcObservables()
    use rtable
    implicit none
    real(kind = rk) :: Upot, Ekin,Eee,Eharm

    if (PCGridSize>0 .or. .not. ThermalEstimator) then
       call initrtable_PBC()
    end if    
    if (PCGridSize>0) then
       CALL calcPairCorrelation(sgn(1),mittausten_lkm(2))
       call calc_ag_density(sgn(1))
    end if


    Ekin=0.0_rk
    
    if (ThermalEstimator) then
       Upot=ThermalPotEnergyEstimator()
       Ekin = ThermalKineticEnergyEstimator()
       Observ(i_tot) = Observ(i_tot) + (Ekin + Upot)/sgn(1)       
       Observ(i_ee) = Observ(i_ee) + eePotentialEnergy()/sgn(1)       
    else
       !Upot = PotentialEnergy()
       !Observ(i_tot) = Observ(i_tot) + Upot/2/sgn(1)
       Eee = eePotentialEnergy()
       Eharm = harmPotentialEnergy()
       Ekin = -Eee/2+Eharm
       Upot = Eee + Eharm
       Observ(i_tot) = Observ(i_tot) + (Eee+4.0_rk*Eharm)/2/sgn(1)
       Observ(i_ee) = Observ(i_ee) + Eee/sgn(1)       
    end if
    Observ(i_pot) = Observ(i_pot) + Upot/sgn(1)
    Observ(i_kin) = Observ(i_kin) + Ekin/sgn(1)
    
    ! Number of measurements in a block
    mittausten_lkm(1) = mittausten_lkm(1) + 1
    ! Number of "measurements" in a block taking
    ! into account the sign
    sgn(2) = sgn(2) + sgn(1)

  end subroutine calcObservables


  function eePotentialEnergy() result(V)
    use rtable
    implicit none
    integer(kind = isoluku) :: i,k,kk
    real(kind=rk) :: oneperdr,ZZ,V1
    real(kind=rk) :: V
    
    V = 0.0_rk
    V1 = 0.0_rk
    
    do k = 1,hi_lkm(3)-1
       do kk = k+1,hi_lkm(3)
          ZZ = hiuk(k)%Z * hiuk(kk)%Z
          do  i = 1,hiuk(k)%Trotter
             oneperdr=evaloneperdr(i,kk,k)
             V1 = V1 + ZZ*oneperdr
          end do
          V = V + V1/hiuk(k)%Trotter
          V1 = 0.0_rk
       end do
    end do
    V = V/dielectric_constant
    
  end function eePotentialEnergy

  function harmPotentialEnergy() result(V)
    use rtable
    implicit none
    integer(kind = isoluku) :: i,k,kk
    real(kind=rk) :: oneperdr,ZZ,V1
    real(kind=rk) :: V
    
    V = 0.0_rk
    V1 = 0.0_rk
    
    do k = 1,hi_lkm(3)
       V1 = 0.0_rk
       do  i = 1,hiuk(k)%Trotter
          V1 = V1 + hiuk(k)%meff*omega2*sum(hiuk(k)%upaikka(:,i)**2)/2
       end do
       V = V + V1/hiuk(k)%Trotter
    end do

  end function harmPotentialEnergy


  subroutine printObservablesMPI(NMax)
    use mpi
    use refbeadmodule
    implicit none
    integer :: i,j
    integer(kind=isoluku) :: NMax
    CHARACTER( LEN = 4 )  :: hichar
    CHARACTER( LEN = 6 )  :: blockchar,eChar
    integer(kind = isoluku) :: ref_bead,beta2_bead
    integer :: myid, numprocs, ierr
    
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
  
    if (myid .eq. 0) then
       if (exchange_on) then
          call getRefBead(ref_bead,beta2_bead)
          do i=1,size(exchange,1)
             WRITE( eChar, '(I6)' ) i
             OPEN( 10, FILE='NumPositives.' // trim(adjustl(eChar)) // '.dat', &
                  position = 'append')
             write(10,*) ref_bead, sum(exchange(i)%PathSign)
             close(10)    
          end do
          do i=1,hi_lkm(1)
             if (hiuk(i)%ParticleKind .ne. 0) then
                WRITE( eChar, '(I6)' ) i
                OPEN( 10, FILE='Perm.' // trim(adjustl(eChar)) // '.dat', &
                     status = 'replace')
                do j=1,hiuk(1)%Trotter
                   write(10,*) hiuk(i)%perm(j)
                end do
                close(10)
             end if
          end do
          open(123,file='ref.bead', status='replace')
          write(123,*) ref_bead
          close(123)
          open(125,file='total_sign.dat', status='replace')
          write(125,*) sgn(1)
          close(125)

       end if
    end if
       
       
    ! total number of measurements
    mittausten_lkm(2) = mittausten_lkm(2) + mittausten_lkm(1)
    ! total "number" of measurements taking the sign into account
    sgn(3) = sgn(3) + sgn(2)

    Observ = Observ/mittausten_lkm(1)

    ObservMPI=0.0_rk
    call MPI_Reduce(Observ,ObservMPI,NumObserv,MPI_DOUBLE_PRECISION &
         ,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    call MPI_Reduce(real(sgn(2),kind=rk),sgnMPI,1,MPI_DOUBLE_PRECISION &
         ,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    
    ! print energy observables
    if (myid .eq. 0) then
       do i=1,NumObserv
          OPEN( 10, FILE=TRIM(adjustl(obs_name(i))) // '.dat', &
               position = 'append')
          write(10,*) ObservMPI(i)/numprocs
          close(10)
       end do
       

       OPEN( 10, FILE='average_sign.dat', &
            position = 'append')
       write(10,*) sgnMPI/mittausten_lkm(1)/numprocs
       close(10)
       
       
       
       ! print the acceptance prob. for the block
       OPEN( 10, FILE='BlockAcceptProb.dat', &
            position = 'append')
       write(10,*) real(BlockAcceptedCount(hi_lkm(3)+1),kind=rk)/&
            ParticleMovedCount(hi_lkm(3)+1)
       close(10)    
       
       do i=1,hi_lkm(3)
          ! print the acceptance prob. for a particle per block
          WRITE( hichar, '(I4)' ) i
          OPEN( 10, FILE='PartBlockAcceptProb.' // trim(adjustl(hichar)) &
               // '.dat', &
               position = 'append')
          if (ParticleMovedCount(i)==0) then
             write(10,*) 0.0_rk
          else
             write(10,*) real(BlockAcceptedCount(i),kind=rk)/ParticleMovedCount(i)
          end if
          close(10)
       end do
       
       if (GetBlockConf) then
          WRITE( blockchar, '(I6)' ) nblock
          
          ! Electron configuration
          OPEN( 10, FILE=TRIM(tiedAlku_e) // '.electrons', STATUS = 'REPLACE')
          OPEN( 13, FILE='Block.' // TRIM(adjustl(blockchar)) // '.electrons', STATUS = 'REPLACE')
          do i = 1,hi_lkm(1)
             do j = 1,hiuk(i)%Trotter
                WRITE( 10, *) hiuk(i)%vpaikka(1,j),hiuk(i)%vpaikka(2,j)
                WRITE( 13, *) hiuk(i)%vpaikka(1,j),hiuk(i)%vpaikka(2,j)
             end do
          end do
          close(13)
          CLOSE( 10 )
          
          ! Nuclei configuration    
          OPEN( 10, FILE=TRIM(tiedAlku_y) // '.nuclei', STATUS = 'REPLACE')
          OPEN( 13, FILE='Block.' // TRIM(adjustl(blockchar)) // '.nuclei', STATUS = 'REPLACE')
          do i = hi_lkm(1)+1,hi_lkm(3)
             do j = 1,hiuk(i)%Trotter
                WRITE( 10, *) hiuk(i)%vpaikka(1,j),hiuk(i)%vpaikka(2,j)
                WRITE( 13, *) hiuk(i)%vpaikka(1,j),hiuk(i)%vpaikka(2,j)
             end do
          end do
          close(13)
          CLOSE( 10 ) 
          
       else
          ! Electron configuration
          OPEN( 10, FILE=TRIM(tiedAlku_e) // '.electrons', STATUS = 'REPLACE')
          do i = 1,hi_lkm(1)
             do j = 1,hiuk(i)%Trotter
                WRITE( 10, *) hiuk(i)%vpaikka(1,j),hiuk(i)%vpaikka(2,j)
             end do
          end do
          CLOSE( 10 )
          
          ! Nuclei configuration    
          OPEN( 10, FILE=TRIM(tiedAlku_y) // '.nuclei', STATUS = 'REPLACE')
          do i = hi_lkm(1)+1,hi_lkm(3)
             do j = 1,hiuk(i)%Trotter
                WRITE( 10, *) hiuk(i)%vpaikka(1,j),hiuk(i)%vpaikka(2,j)
             end do
          end do
          CLOSE( 10 )    
          
       end if
       
       
    end if

    if (PCGridSize>0) then
       call printPairCorrelationMPI(sgn(2))
       call print_ag_densityMPI()
    end if
    
       

  end subroutine printObservablesMPI


  subroutine printConfBackUp
    use refbeadmodule
    implicit none
    integer :: i,j
    CHARACTER( LEN = 6 )  :: eChar
    integer(kind = isoluku) :: ref_bead,beta2_bead
        
    ! Electron configuration
    OPEN( 10,file='e.backup', STATUS = 'REPLACE')
    do i = 1,hi_lkm(1)
       do j = 1,hiuk(i)%Trotter
          WRITE( 10, *) hiuk(i)%vpaikka(1,j),hiuk(i)%vpaikka(2,j)
       end do
    end do
    CLOSE( 10 )

    ! Nuclei configuration    
    OPEN( 10, file='n.backup', STATUS = 'REPLACE')
    do i = hi_lkm(1)+1,hi_lkm(3)
       do j = 1,hiuk(i)%Trotter
          WRITE( 10, *) hiuk(i)%vpaikka(1,j),hiuk(i)%vpaikka(2,j)
       end do
    end do
    CLOSE( 10 )    

    if (exchange_on) then
       call getRefBead(ref_bead,beta2_bead)
       do i=1,hi_lkm(3)
          if (hiuk(i)%ParticleKind==1) then
             WRITE( eChar, '(I6)' ) i
             OPEN( 10, FILE='Perm.' // trim(adjustl(eChar)) // '.backup', &
                  status = 'replace')
             do j=1,hiuk(1)%Trotter
                write(10,*) hiuk(i)%perm(j)
             end do
             close(10)
          end if
       end do
       open(123,file='ref.bead.backup', status='replace')
       write(123,*) ref_bead
       close(123)
    end if



  end subroutine printConfBackUp
  


  function PotentialEnergy() result(V)
    use rtable
    implicit none
    integer(kind = isoluku) :: i,k,kk
    real(kind=rk) :: oneperdr,ZZ,V1
    real(kind=rk) :: V
    
    V = 0.0_rk
    V1 = 0.0_rk
    
    do k = 1,hi_lkm(3)-1
       do kk = k+1,hi_lkm(3)
          ZZ = hiuk(k)%Z * hiuk(kk)%Z
          do  i = 1,hiuk(k)%Trotter
             oneperdr=evaloneperdr(i,kk,k)
             V1 = V1 + ZZ*oneperdr
          end do
          V = V + V1/hiuk(k)%Trotter
          V1 = 0.0_rk
       end do
    end do
    V = V/dielectric_constant
    
    do k = 1,hi_lkm(3)
       V1 = 0.0_rk
       do  i = 1,hiuk(k)%Trotter
          V1 = V1 + 2.0_rk*hiuk(k)%meff*omega2*sum(hiuk(k)%upaikka(:,i)**2)
       end do
       V = V + V1/hiuk(k)%Trotter
    end do


  end function PotentialEnergy

  function ThermalPotEnergyEstimator() result(E)
    implicit none
    integer(kind=isoluku) :: hi,kk,j,jj,he1,he2,i,hi1p,hi2p
    real(kind = rk) :: V,V1,r1,r2,x,y,s,E,ZZ,Vag1,Vag2
    real(kind = rk), dimension(dimensio) :: r1v,r2v
    
    V = 0.0_rk
    V1 = 0.0_rk
    E = 0.0_rk
    
    !write(*,*) 'Thermal Estimator'
    
    do hi=1,hi_lkm(3)-1       
       do kk=hi+1,hi_lkm(3)
          
          do he1=1,hiuk(hi)%Trotter
          
             if (he1==hiuk(hi)%Trotter) then
                he2 = 1
             else
                he2 = he1 + 1
             end if

             ! Trotter numbers must be the same
             ! for all moving particles
             if (hiuk(kk)%Trotter == 1) then
                j = 1
                jj = 1
             else
                j = he1
                jj = he2
             end if
           
             ZZ = hiuk(hi)%Z * hiuk(kk)%Z
             
             
             r1v = hiuk(hi)%upaikka(:,he1)-&
                  hiuk(kk)%upaikka(:,j)
             r2v = hiuk(hi)%upaikka(:,he2)-&
                  hiuk(kk)%upaikka(:,jj)
             
             r1 = sqrt(sum(r1v**2))
             r2 = sqrt(sum(r2v**2))
             
             if ( ZZ < 0.0) then          
                ! ....
             else                            
                V1 = primAction2(r1,r2,ZZ)
             end if
             V = V + V1
             V1 = 0.0_rk             
          end do
          E = E + V/hiuk(hi)%Trotter
          V = 0.0_rk
          
       end do
    end do
    
    do i=1,hi_lkm(3)
       V=0.0_rk
       do j=1,hiuk(i)%Trotter          
          V = V+hiuk(i)%meff*omega2*&
               sum(hiuk(i)%upaikka(:,j)**2)/2
          if (.not. NoMagneticField) then
             V = V+hiuk(i)%meff*omegaL2*&
                  sum(hiuk(i)%upaikka(:,j)**2)/2
          end if
       end do
       E=E+V/hiuk(i)%Trotter
    end do
    

  end function ThermalPotEnergyEstimator




  function ThermalKineticEnergyEstimator() result(E)
    implicit none
    real(kind=rk) :: IntE,dr2,E
    integer :: hi,nn
    integer(kind=isoluku) :: j,hi1p,hi2p
    
    E = 0.0_rk
    do hi=1,s_lkm
       IntE = 0.0_rk
       if (hiuk(hi)%lambda>1.0e-10_rk) then          
          hi1p=hi
          hi2p=hi
          do j = 1,hiuk(hi)%Trotter
             hi1p=hi2p
             hi2p=hiuk(hi1p)%perm(j)
             if (j < hiuk(hi)%Trotter) then
                dr2 = sum((hiuk(hi2p)%upaikka(:,j+1)-&
                     hiuk(hi1p)%upaikka(:,j))**2)
             else
                dr2 = sum((hiuk(hi1p)%upaikka(:,j)-&
                     hiuk(hi2p)%upaikka(:,1))**2)
             end if
             IntE = IntE + &
                  dr2/hiuk(hi)%lambda/4/tau_in**2
          end do
          IntE = real(dimensio,kind=rk)/tau_in/2 - IntE/hiuk(hi)%Trotter
       end if
       E = E + IntE
    end do

  end function ThermalKineticEnergyEstimator

  function primAction2(r,rp,ZZ) result(V)
    implicit none
    real(kind = rk) :: r,rp,ZZ,V,x,y
    
    if (r<1.0e-7) then
       x = 1.0e-7_rk
    else
       x = r
    end if
    if (rp<1.0e-7) then
       y = 1.0e-7_rk
    else
       y = rp
    end if
    
    V = 0.5_rk*ZZ*(1.0_rk/x + 1.0_rk/y)/dielectric_constant
    
  end function primAction2


  
  subroutine UpdateMultilevelL()    
    implicit none
    integer :: i
    real(kind =rk) :: prob

    do i=1,s_lkm
       if (ParticleMovedCount(i)>3) then
          prob = real(BlockAcceptedCount(i),kind=rk)/ParticleMovedCount(i)
          if (prob>AccHigh .and. 2**hiuk(i)%multi_L<hiuk(i)%Trotter/2) then
             hiuk(i)%multi_L = hiuk(i)%multi_L + 1
          elseif (prob<AccLow .and. hiuk(i)%multi_L>1) then
             hiuk(i)%multi_L = hiuk(i)%multi_L - 1
          end if
       end if
    end do

    

  end subroutine UpdateMultilevelL
  
 

  

end module observables
