

MODULE init

  use math
  
  IMPLICIT NONE
  
  TYPE:: hi_tiedot
     
     INTEGER(kind = isoluku) :: Trotter
     INTEGER :: lkm ! number of similar particles
     INTEGER :: species ! number of the species
     integer :: first
     integer :: last
     INTEGER                 :: multi_L
     INTEGER                 :: ParticleKind !0->Boltzmannon, 1->Fermion
     CHARACTER( len=3 )      :: tyyppi
     CHARACTER( len=11 )     :: statistics !Boltzmannon or Fermion
     REAL(KIND=RK)               :: Z
     REAL(KIND=RK)               :: meff
     REAL(KIND=RK)               :: lambda
     REAL(KIND=RK)               :: spin
     REAL(KIND=RK), DIMENSION(:,:), ALLOCATABLE  :: vpaikka
     REAL(KIND=RK), DIMENSION(:,:), ALLOCATABLE  :: upaikka
     INTEGER(kind = isoluku), DIMENSION(:), ALLOCATABLE :: perm
     INTEGER(kind = isoluku), DIMENSION(:), ALLOCATABLE :: perm_old
     INTEGER(kind = isoluku), DIMENSION(:), ALLOCATABLE :: pn
     INTEGER(kind = isoluku), DIMENSION(:), ALLOCATABLE :: po
     
     !SVA@CSC siirretytpaikat contains a list of beads 
     !which have moved during one MC trial step. Used to
     !selectively copy bead positions between upaikka nad vpaikka
     INTEGER(kind = isoluku), DIMENSION(:), ALLOCATABLE ::siirretytpaikat 
     !SVA@CSC  siirretytpaikatlkm contains the number bead indexes in siirretytpaikat 
     INTEGER:: siirretytpaikatlkm     
     
  END TYPE hi_tiedot
   
  TYPE( hi_tiedot ), DIMENSION(:), ALLOCATABLE,SAVE  :: hiuk


  character(len=11), dimension(4), private :: ParticleKind = &
       (/'Boltzmannon','Fermion    ','Boson      ','ExactFermi '/)

  
  TYPE:: exchangeMatrix
    
    INTEGER(kind = isoluku)            :: Trotter
    INTEGER            :: lkm ! number of ex particles
    INTEGER            :: nauha
    INTEGER            :: first
    INTEGER            :: last
    REAL(KIND=RK)               :: Z
    REAL(KIND=RK)               :: meff
    REAL(KIND=RK)               :: lambda
    REAL(KIND=RK)               :: spin
    integer, DIMENSION(:), ALLOCATABLE  :: PathSign
    integer, DIMENSION(:), ALLOCATABLE  :: PathSignOld
    REAL(KIND=RK), DIMENSION(:,:), ALLOCATABLE  :: Matrix
    
  END TYPE exchangeMatrix
  
  TYPE( exchangeMatrix ), DIMENSION(:), ALLOCATABLE,SAVE  :: exchange

  
  TYPE:: effpot_tiedot
     
     REAL(KIND=RK)    :: ZZ
     REAL(KIND=RK)    :: m_red
     REAL(KIND=RK)    :: lambda_red
     integer          :: modcusp 
     integer          :: potnro = 0     
     
  END TYPE effpot_tiedot
  
  TYPE( effpot_tiedot ), DIMENSION(:,:), ALLOCATABLE,SAVE  :: effpot_data

  integer, save :: NumPairPots
  integer, save :: pairpot_dim
  integer(kind = isoluku), save :: NumSquares

  logical, save :: GetBlockConf

  ! if DisplaceMoves=0.1 then about 10% 
  ! of the moves are displace moves
  ! and 90% are multilevel bisection moves
  real(kind = rk), save :: DisplaceMoves



  integer(kind=isoluku), dimension(:), allocatable :: MovedBeads
  integer(kind=isoluku) :: NumMovedBeads     


  integer, save       :: desired_tau_on
  real(kind=rk), save :: desired_tau
  integer(kind=isoluku), private, save :: ModTrotter

  logical, save :: NoMagneticField
  real(kind=rk), save :: omegaL ! Landau frequency
  real(kind=rk), save :: omegaL2 ! Landau frequency squared

  integer(kind=isoluku) :: block_cut

  logical, save :: isfixednode

  integer(kind=isoluku), save :: sign_multiply

  integer, save :: NumSpecies
  integer, save :: NumSpeciesEQ1

  integer(kind = isoluku),save   :: NumOfBlocks
  integer, save                  :: PCGridSize
  real(kind=rk), save            :: PCLowLimit
  real(kind=rk), save            :: PCUppLimit

  logical, save :: LevelUpdateOn
  real(kind=rk), save :: AccLow
  real(kind=rk), save :: AccHigh       
  

  real(kind=rk), save            :: beta,tau_in
  REAL(KIND=RK), SAVE            :: Temperature
  INTEGER, SAVE                  :: otosvali
  
  REAL(KIND=RK), PARAMETER       :: a0 = 0.529177208e-10_rk
  REAL(KIND=RK), PARAMETER       :: ev = 3.6749326014e-2_rk
  REAL(KIND=RK), PARAMETER       :: hviiva =  1.0_rk

  INTEGER, DIMENSION(3), SAVE  :: hi_lkm
  INTEGER, SAVE                :: s_lkm

  LOGICAL, save                :: alusta_satunnais
  LOGICAL, SAVE                :: bisection_hyv
  LOGICAL, SAVE                :: exchange_on
  LOGICAL, SAVE                :: ThermalEstimator
  CHARACTER(LEN=30), SAVE      :: otsikko
  character(len=40), save      :: tiedAlku_e, tiedAlku_y





CONTAINS

  SUBROUTINE alustamuuttujat
    use math
    IMPLICIT NONE
    
    INTEGER               :: i, j
    LOGICAL               :: olemassa
  
    CHARACTER(LEN=20)     :: TrotterChar,NParticlesChar
    CHARACTER(LEN=20)     :: nimi,int_str
    real(kind=rk), dimension(2) :: g_rand
    real(kind=rk), dimension(hi_lkm(1),2) :: r_init
    real(kind=rk) :: dr,dr2

    
    ! --------------------------
    ! Elektronit (or negative particles)
    ! --------------------------

    !elektroneilla sama Trotterin luku
    WRITE( TrotterChar, '(I10)'  ) hiuk(1)%Trotter 
    WRITE( NParticlesChar, '(I5)'  ) hi_lkm(1)

       
    nimi = 'M' // TRIM(ADJUSTL(TrotterChar)) // 'N' //&
         TRIM(ADJUSTL(NParticlesChar))
    
    INQUIRE ( FILE = TRIM(nimi) // '.electrons', EXIST=olemassa )

    tiedAlku_e = nimi

    IF ( olemassa ) THEN
       
       OPEN( 10, FILE=TRIM(nimi) // '.electrons' )
       
       DO i = 1, hi_lkm(1) 
          DO j = 1, hiuk(i)%Trotter             
             READ( 10, * ) hiuk(i)%vpaikka(:,j)             
          END DO
       END DO
       
       CLOSE( 10 )
       
    ELSE       
       
       WRITE(*,*) "------------------------------------------------"
       WRITE(*,*) "       Generating initial configuration         "
       WRITE(*,*) "------------------------------------------------"

       DO i = 1, hi_lkm(1) 
          g_rand=gaussinen()          
          do j = 1, hiuk(i)%Trotter
             hiuk(i)%vpaikka(:,j) = g_rand*&
                  sqrt(2.0_rk*hiuk(i)%lambda/sqrt(omega2))
          end do
          hiuk(i)%upaikka = hiuk(i)%vpaikka
          hiuk(i)%siirretytpaikatlkm=0
       END DO

       !DO i = 1, hi_lkm(1) 
       !   hiuk(i)%vpaikka = 0.0_rk
       !   hiuk(i)%upaikka = 0.0_rk
       !   hiuk(i)%siirretytpaikatlkm=0
       !END DO
       
    END IF

   
    do i = 1, hi_lkm(1) !only electrons in this code
       WRITE( int_str, '(I5)'  ) i
       INQUIRE ( FILE = 'Perm.' // trim(adjustl(int_str)) // '.dat', EXIST=olemassa )
       if (olemassa) then
          open( 10, file='Perm.' // trim(adjustl(int_str)) // '.dat')
          do j = 1, hiuk(i)%Trotter
             read(10,*) hiuk(i)%perm(j)
          end do
          close(10)
          hiuk(i)%perm_old=hiuk(i)%perm
       end if
    end do


    ! -------------------
    ! Ytimet (positive particles)
    ! -------------------

    ! ytimilla sama Trotterin luku
    if (hi_lkm(2)==0) then
       WRITE( TrotterChar, '(I10)'  ) 0 
    else
       WRITE( TrotterChar, '(I10)'  ) hiuk(hi_lkm(1)+1)%Trotter 
    end if
    WRITE( NParticlesChar, '(I5)'  ) hi_lkm(2)
        
   !WRITE( TrotterChar, '(I10)'  ) hiuk(hi_lkm(1)+1)%Trotter 
   WRITE( NParticlesChar, '(I5)'  ) hi_lkm(2)

       
   nimi = 'M' // TRIM(ADJUSTL(TrotterChar)) // 'N' //&
        & TRIM(ADJUSTL(NParticlesChar))
    
   INQUIRE ( FILE = TRIM(nimi) // '.nuclei', EXIST=olemassa )

   tiedAlku_y = nimi
    
   IF ( olemassa ) THEN
       
      OPEN( 10, FILE=TRIM(nimi) // '.nuclei' )
       
      DO i = hi_lkm(1)+1,hi_lkm(3) 
         DO j = 1, hiuk(i)%Trotter            
            READ( 10, * ) hiuk(i)%vpaikka(:,j)            
         END DO
      END DO

      CLOSE( 10 )
           
   ELSE

      if (hi_lkm(2) .ne. 0) then
         WRITE(*,*) "--------------------------------------------"
         WRITE(*,*) "Olis syytä olla ytimien paikat jo tiedossa!!"
         WRITE(*,*) "--------------------------------------------"
         
         DO i = hi_lkm(1)+1, hi_lkm(3) 
            hiuk(i)%vpaikka = 0.0_rk
            hiuk(i)%upaikka = 0.0_rk
            hiuk(i)%siirretytpaikatlkm=0
         END DO
      end if
      
   END IF
        
  END SUBROUTINE alustamuuttujat
  

  !-----------------------------------------------
  ! Alustetaan satunnaislukugeneraattori
  !

  SUBROUTINE alustasatunnaisgeneraattori( siemen1 ) 
    
    IMPLICIT NONE
    
    INTEGER, INTENT(OUT)               :: siemen1
    INTEGER, DIMENSION(:), ALLOCATABLE :: siemen
    INTEGER, DIMENSION(8)              :: t
    INTEGER                            :: koko, tila


    !Tehdään vähän satunnaisuutta
    CALL RANDOM_SEED(size=koko)

    ALLOCATE(siemen(koko), stat=tila)
    IF ( tila > 0 ) STOP 'Tilanvaraus epäonnistui'

    IF ( alusta_satunnais ) THEN
       CALL DATE_AND_TIME(values = t)
       !siemen = 100*t(7) + t(8)/10 + 1
       siemen = 1000*t(7) + t(8) + 1

       IF (siemen(1) == 0) THEN
          siemen = 1768
       END IF
    ELSE
       siemen = 1768
    END IF

    CALL RANDOM_SEED( put=siemen )

    siemen1 = siemen(1) 

  END SUBROUTINE alustasatunnaisgeneraattori

  SUBROUTINE alustasatunnaisgeneraattoriMPI( siemen1, ind ) 
    
    IMPLICIT NONE
    
    INTEGER, INTENT(OUT)               :: siemen1
    INTEGER, INTENT(in)                :: ind
    INTEGER, DIMENSION(:), ALLOCATABLE :: siemen
    INTEGER, DIMENSION(8)              :: t
    INTEGER                            :: koko, tila


    !Tehdään vähän satunnaisuutta
    CALL RANDOM_SEED(size=koko)

    ALLOCATE(siemen(koko), stat=tila)
    IF ( tila > 0 ) STOP 'Tilanvaraus epäonnistui'

    IF ( alusta_satunnais ) THEN
       CALL DATE_AND_TIME(values = t)
       !siemen = 100*t(7) + t(8)/10 + 1
       siemen = 1000*t(7) + t(8)*ind

       IF (siemen(1) == 0) THEN
          siemen = 1768 + ind
       END IF
    ELSE
       siemen = 1769
    END IF

    CALL RANDOM_SEED( put=siemen )

    siemen1 = siemen(1) 

  END SUBROUTINE alustasatunnaisgeneraattoriMPI





  SUBROUTINE read_INPUT(blokkeja,NMax)
    use mpi
    use math
    IMPLICIT NONE
    
    INTEGER(kind = isoluku) :: blokkeja,NMax,Trotter
    integer                 :: i,nP1,nP_tot
    INTEGER               :: lkm, lkm2,spin, multi_L,ex_lkm
    REAL(KIND=RK)         :: Z, m, coef,Tfermi,unit_conv
    real(kind = rk)       :: Lx,Ly,Lz, temp,omega,Tdesired,omega_in_meV
    integer               :: tauORtemp,intParticleKind
    integer               :: intDisplaceMoveOn, ReadFile, desired_tau_units
    INTEGER               :: allocstat, ex_im_in
    !CHARACTER( len=5 )    :: exchange_on_text
    LOGICAL               :: olemassa
    integer               :: ThermEstim,LRon,GetConf,Lupdateon,eh_ex
    CHARACTER( len=20 )   :: estimator,paircors,helptext,LongRange,Conf,Lupdate
    integer :: myid, numprocs, ierr
    
    
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
    
    ThermalEstimator = .false.
    
    INQUIRE ( FILE = 'INPUT', EXIST=olemassa )
    
    IF ( olemassa ) THEN
       
       lkm2 = 0
       ex_lkm = 0
       hi_lkm = 0

       OPEN( 10, FILE='INPUT')
       
       READ( 10, *) otsikko
       READ( 10, *) nP1, nP_tot

       NumSpecies = nP1

       hi_lkm(3) = nP_tot

       ! varataan muistia hiukkasille
       ALLOCATE( hiuk( nP_tot ), STAT = allocstat )
       IF ( allocstat /= 0 ) THEN
          WRITE(*,*) 'Virhe muistinvarauksessa'
          STOP
       END IF

       !luetaan hiukkasten datat
       NumSpeciesEQ1=0
       DO i = 1,nP1
          READ( 10, *) lkm, Z, m, spin, Trotter, multi_L,intParticleKind
          if (lkm==1) then
             NumSpeciesEQ1 = NumSpeciesEQ1 + 1
          end if
             
          CALL alustaHiukkaset ( lkm2, lkm, Z, m, &
               spin, Trotter, multi_L,intParticleKind,i)
          lkm2 = lkm2 + lkm
          IF (lkm>1 .and. intParticleKind .ne. 0) THEN
             ex_lkm = ex_lkm + 1
          END IF
       END DO
       
       
       READ( 10, *) s_lkm
       READ( 10, *) otosvali
       READ( 10, *) temp, tauORtemp
       READ( 10, *) blokkeja, NMax
       READ( 10, *) estimator, ThermEstim
       READ( 10, *) paircors, PCLowLimit, PCUppLimit, PCGridSize, block_cut
       READ( 10, *) Conf, GetConf
       READ( 10, *) Lupdate,Lupdateon,AccLow,AccHigh
       read( 10, *) helptext, dielectric_constant, omega_in_meV, omegaL
       read( 10, *) helptext, desired_tau_on, desired_tau, desired_tau_units
              
       CLOSE( 10 )

       unit_conv = dielectric_constant**2/hiuk(1)%meff
       omega=omega_in_meV*meV2Hartree
       !write(*,*) omega
       omega2=omega**2
       if (omegaL>0.0_rk) then
          NoMagneticField=.false.
          omegaL2=omegaL**2
       else
          NoMagneticField=.true.
       end if
       
       if (hiuk(1)%ParticleKind==1) then
          isfixednode=.true.
       else
          isfixednode=.false.
       end if

       if (Lupdateon==1) then
          LevelUpdateOn = .true.
       else
          LevelUpdateOn = .false.
       end if
       
       if (GetConf==1) then
          GetBlockConf = .true.
       else
          GetBlockConf = .false.
       end if
       
       if (ThermEstim==1) then
          ThermalEstimator = .true.
       else
          ThermalEstimator = .false.
       end if

       if (hi_lkm(1)==2 .and. ex_lkm==1) then
          ! spin polarized for two electrons
          coef = 2.0_rk
       elseif (hi_lkm(1)==hiuk(1)%lkm .and. ex_lkm==1) then
          !spin polarized
          coef = 2.0_rk
       else
          !unpolarized spin (or something else)
          coef = 1.0_rk
       end if
       
       if (tauORtemp==0) then
          tau_in = temp
       elseif (tauORtemp==1) then
          tau_in = 1.0_rk/kB/temp/hiuk(1)%Trotter
       else
          temp=temp*omega/kB*sqrt(coef*hi_lkm(1))
          tau_in = 1.0_rk/kB/temp/hiuk(1)%Trotter
       end if
              
       Tfermi=omega/kB*sqrt(coef*hi_lkm(1))
           
    ELSE
       WRITE(*,*) 'You need to have an INPUT file'
    END IF

       
    Tdesired = 1.0_rk/(tau_in*hiuk(1)%Trotter)
    if (desired_tau_on==1) then
       if (desired_tau_units==0) then          
          ModTrotter=floor(1.0_rk/desired_tau/Tdesired*0.999999999999_rk)+1
       else
          ModTrotter=floor(1.0_rk/desired_tau/Tdesired*0.999999999999_rk*omega*sqrt(coef*hi_lkm(1)))+1
       end if
       if (modulo(ModTrotter,2) .ne. 0) then
          ModTrotter=ModTrotter+1
       end if
       tau_in=1.0_rk/(ModTrotter*Tdesired)
       call alustaHiukkaset2()
       !if (ex_lkm>0) then
       !   call alustaExchange2()
       !end if
    end if

    if (myid .eq. 0) then
       open(123,file='general_info.txt', status='replace')
       write(123,*) '2d Harmonic dot'
       write(123,*) ' Num of electrons ', hi_lkm(3)
       if (coef<1.5_rk) then
          write(123,*) ' unpolarized '
       else
          write(123,*) ' polarized '
       end if
       write(123,*) ' omega (Ha)       ', omega
       write(123,*) ' omega (eff Ha)   ', omega*unit_conv
       write(123,*) ' omega (meV)      ', omega_in_meV
       write(123,*) ' T(K)             ', 1.0_rk/(tau_in*kB*hiuk(1)%Trotter)
       write(123,*) ' T_F(K)           ', Tfermi
       write(123,*) ' T/T_F            ', 1.0_rk/(tau_in*hiuk(1)%Trotter)/omega/sqrt(coef*hi_lkm(1))
       write(123,*) ' '
       write(123,*) ' time-step      ', tau_in
       write(123,*) ' Trotter number ', hiuk(1)%Trotter
       write(123,*) ' Notice: time-step in (E_F)**(-1) ', tau_in*omega*sqrt(coef*hi_lkm(1))
       if (desired_tau_on==1) then
          write(123,*) ' '
          write(123,*) ' Trotter Number was changed to match the'
          if (desired_tau_units==0) then
             write(123,*) ' desired time-step in Ha**(-1)   ', desired_tau
             write(123,*) '             abs(desired-actual) ', abs(desired_tau-tau_in)
             write(123,*) ' desired time-step in (E_F)**(-1)   ', desired_tau*omega*sqrt(coef*hi_lkm(1))
             write(123,*) '             abs(desired-actual) ', abs(desired_tau-tau_in)*omega
          else
             write(123,*) ' desired time-step in Ha**(-1)   ', desired_tau/omega/sqrt(coef*hi_lkm(1))
             write(123,*) '             abs(desired-actual) ', abs(desired_tau/omega/sqrt(coef*hi_lkm(1))-tau_in)
             write(123,*) ' desired time-step in (E_F)**(-1)   ', desired_tau
             write(123,*) '             abs(desired-actual) ', abs(desired_tau-tau_in*omega*sqrt(coef*hi_lkm(1)))
          end if
       end if
       close(123)
    end if

    dielectric_constant=1.0_rk
    omega = omega*unit_conv
    omega2 = omega**2
    tau_in = tau_in/unit_conv
    do i=1,hi_lkm(1)
       hiuk(i)%meff = 1.0_rk
       hiuk(i)%lambda = 0.5_rk
    end do

    IF (ex_lkm > 0) THEN
       exchange_on = .TRUE.
       ALLOCATE( exchange( ex_lkm ), STAT = allocstat )
       IF ( allocstat /= 0 ) THEN
          WRITE(*,*) 'Virhe muistinvarauksessa'
          STOP
       END IF
       CALL alustaExchange()
    ELSE
       exchange_on = .FALSE.
    END IF
    

    allocate(MovedBeads(hiuk(1)%Trotter))
    MovedBeads=0
    NumMovedBeads=0   

  END SUBROUTINE read_INPUT

  SUBROUTINE alustaHiukkaset2()
    use math
    IMPLICIT NONE
    INTEGER :: hi,SpeciesNum
    INTEGER :: lkm, lkm2, spin, multi_L,intParticleKind
    integer(kind = isoluku) :: Trotter
    REAL(KIND=RK)    :: Z, m
    INTEGER :: allocstat

    DO hi = 1, hi_lkm(3)
       
       hiuk(hi)%Trotter = ModTrotter
       
       DEALLOCATE(hiuk(hi)%vpaikka)
       DEALLOCATE(hiuk(hi)%upaikka)
       ALLOCATE(hiuk(hi)%vpaikka(dimensio,hiuk(hi)%Trotter),&
            STAT=allocstat)
       ALLOCATE(hiuk(hi)%upaikka(dimensio,hiuk(hi)%Trotter),&
            STAT=allocstat)
       
       hiuk(hi)%vpaikka = 0.0_rk
       hiuk(hi)%upaikka = 0.0_rk
       
       DEALLOCATE(hiuk(hi)%siirretytpaikat)
       ALLOCATE(hiuk(hi)%siirretytpaikat(hiuk(hi)%Trotter),&
            STAT=allocstat)
       hiuk(hi)%siirretytpaikat = 0
       hiuk(hi)%siirretytpaikatlkm = 0

       DEALLOCATE(hiuk(hi)%perm)
       DEALLOCATE(hiuk(hi)%perm_old)
       DEALLOCATE(hiuk(hi)%pn)
       DEALLOCATE(hiuk(hi)%po)
       ALLOCATE(hiuk(hi)%perm(hiuk(hi)%Trotter),&
            STAT=allocstat)
       ALLOCATE(hiuk(hi)%perm_old(hiuk(hi)%Trotter),&
            STAT=allocstat)
       ALLOCATE(hiuk(hi)%pn(hiuk(hi)%Trotter),&
            STAT=allocstat)
       ALLOCATE(hiuk(hi)%po(hiuk(hi)%Trotter),&
            STAT=allocstat)
       hiuk(hi)%perm=hi
       hiuk(hi)%perm_old=hi
       hiuk(hi)%pn=hi
       hiuk(hi)%po=hi       

     END DO
     
    IF ( allocstat /= 0 ) THEN
       WRITE(*,*) 'Virhe muistinvarauksessa'
       STOP
    END IF
    
  END SUBROUTINE alustaHiukkaset2

  SUBROUTINE alustaExchange2
    use math
    IMPLICIT NONE
    INTEGER :: i, j
    INTEGER :: allocstat

    j = 1
    i = 1
    do j=1,size(exchange)
       exchange(j)%Trotter = ModTrotter

       DEALLOCATE(exchange(j)%PathSign)
       DEALLOCATE(exchange(j)%PathSignOld)       
       ALLOCATE(exchange(j)%PathSign(hiuk(i)%Trotter),&
            STAT=allocstat)
       ALLOCATE(exchange(j)%PathSignOld(hiuk(i)%Trotter),&
            STAT=allocstat)
       
       exchange(j)%PathSign = 0
       exchange(j)%PathSignOld = 0
       exchange(j)%Matrix = 0.0_rk
    end do

    IF ( allocstat /= 0 ) THEN
       WRITE(*,*) 'Virhe muistinvarauksessa'
       STOP
    END IF
    
  END SUBROUTINE alustaExchange2



  SUBROUTINE alustaHiukkaset( lkm2, lkm, Z, m, spin, Trotter, multi_L,intParticleKind,SpeciesNum)
    use math
    IMPLICIT NONE
    INTEGER :: hi,SpeciesNum
    INTEGER :: lkm, lkm2, spin, multi_L,intParticleKind
    integer(kind = isoluku) :: Trotter
    REAL(KIND=RK)    :: Z, m
    INTEGER :: allocstat

    DO hi = lkm2 + 1, lkm2 + lkm

       hiuk(hi)%Trotter = Trotter
       hiuk(hi)%lkm     = lkm
       hiuk(hi)%first   = lkm2+1
       hiuk(hi)%last    = lkm2+lkm
       hiuk(hi)%species = SpeciesNum
       hiuk(hi)%meff    = m
       hiuk(hi)%multi_L = multi_L
       hiuk(hi)%spin    = spin
       hiuk(hi)%Z       = Z
       hiuk(hi)%ParticleKind = intParticleKind
       hiuk(hi)%statistics = ParticleKind(intParticleKind+1)
       IF ( hiuk(hi)%Z < 0.) THEN
          hiuk(hi)%tyyppi  = 'neg'
          hi_lkm(1) = hi_lkm(1) + 1
       ELSE
          hiuk(hi)%tyyppi  = 'pos'
          hi_lkm(2) = hi_lkm(2) + 1
       END IF       
       
       if (hiuk(hi)%meff>1.0e10_rk .or. hiuk(hi)%meff<1.0e-10_rk) then
          hiuk(hi)%lambda = 0.0_rk
       else
          hiuk(hi)%lambda = 1.0_rk/hiuk(hi)%meff/2
       end if
       
       ALLOCATE(hiuk(hi)%vpaikka(dimensio,hiuk(hi)%Trotter),&
            STAT=allocstat)
       ALLOCATE(hiuk(hi)%upaikka(dimensio,hiuk(hi)%Trotter),&
            STAT=allocstat)

       hiuk(hi)%vpaikka = 0.0_rk
       hiuk(hi)%upaikka = 0.0_rk

       
       !ALLOCATE(hiuk(hi)%siirretytpaikat(2**hiuk(hi)%multi_L),&
       !    STAT=allocstat)
       ALLOCATE(hiuk(hi)%siirretytpaikat(hiuk(hi)%Trotter),&
            STAT=allocstat)
       hiuk(hi)%siirretytpaikat = 0
       hiuk(hi)%siirretytpaikatlkm = 0

       ALLOCATE(hiuk(hi)%perm(hiuk(hi)%Trotter),&
            STAT=allocstat)
       ALLOCATE(hiuk(hi)%perm_old(hiuk(hi)%Trotter),&
            STAT=allocstat)
       ALLOCATE(hiuk(hi)%pn(hiuk(hi)%Trotter),&
            STAT=allocstat)
       ALLOCATE(hiuk(hi)%po(hiuk(hi)%Trotter),&
            STAT=allocstat)
       hiuk(hi)%perm=hi
       hiuk(hi)%perm_old=hi
       hiuk(hi)%pn=hi
       hiuk(hi)%po=hi       

     END DO
     
    IF ( allocstat /= 0 ) THEN
       WRITE(*,*) 'Virhe muistinvarauksessa'
       STOP
    END IF
    
  END SUBROUTINE alustaHiukkaset

  SUBROUTINE alustaExchange
    use math
    IMPLICIT NONE
    INTEGER :: i, j
    INTEGER :: allocstat

    j = 1
    i = 1
    DO WHILE (i<hi_lkm(3))
       IF (hiuk(i)%lkm > 1 .and. hiuk(i)%ParticleKind .ne. 0) THEN
          exchange(j)%Trotter = hiuk(i)%Trotter
          exchange(j)%meff    = hiuk(i)%meff
          exchange(j)%lambda  = hiuk(i)%lambda
          exchange(j)%lkm     = hiuk(i)%lkm
          exchange(j)%nauha   = i
          exchange(j)%first   = hiuk(i)%first
          exchange(j)%last   = hiuk(i)%last
          exchange(j)%Z       = hiuk(i)%Z
          
          ALLOCATE(exchange(j)%Matrix(hiuk(i)%lkm,hiuk(i)%lkm),&
               STAT=allocstat)
          ALLOCATE(exchange(j)%PathSign(hiuk(i)%Trotter),&
               STAT=allocstat)
          ALLOCATE(exchange(j)%PathSignOld(hiuk(i)%Trotter),&
               STAT=allocstat)

          exchange(j)%PathSign = 0
          exchange(j)%PathSignOld = 0
          exchange(j)%Matrix = 0.0_rk
          j = j + 1
          i = i + hiuk(i)%lkm
       ELSE
          i = i+1
       END IF
    END DO

    IF ( allocstat /= 0 ) THEN
       WRITE(*,*) 'Virhe muistinvarauksessa'
       STOP
    END IF
    
  END SUBROUTINE alustaExchange

  

END MODULE init
