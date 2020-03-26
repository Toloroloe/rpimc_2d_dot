PROGRAM PathIntegral
  
!Physics is like sex: sure, it may give some practical results,
!but that's not why we do it. --Richard Feynman.

  use mpi
  use math
  USE init
  USE montecarlo
  USE observables
  use rtable
  use refbeadmodule

  IMPLICIT NONE
  
  INTEGER(kind = isoluku)              :: blokkeja
  INTEGER(kind = isoluku)              :: NMax
  INTEGER, PARAMETER                   :: kanava = 6
  INTEGER                              :: allocstat, i, j, blk
  INTEGER                              :: siemen1
  CHARACTER(LEN=8)                     :: pvm
  CHARACTER(LEN=10)                    :: aloitusaika
  CHARACTER(LEN=30)                    :: teksti
  CHARACTER(LEN=10)                    :: n1,n2
  real(kind = rk), dimension(dimensio) :: ag_box
  logical :: olemassa
  integer(kind=isoluku) :: refbead_in
  integer :: myid, numprocs, ierr
  

  call MPI_INIT(ierr) 
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

  CALL read_INPUT(blokkeja,NMax)
      
  NumOfBlocks = blokkeja
  
  beta = tau_in*real(hiuk(1)%Trotter,kind=rk)
  Temperature = 1.0_rk/beta/kB 

  if (numprocs>1) then
     
     alusta_satunnais  = .false.
     siemen1 = 0
     
     CALL alustasatunnaisgeneraattoriMPI( siemen1, myid)  
     
     CALL alustamuuttujat()
     
     alusta_satunnais  = .TRUE.
     siemen1 = 0
     
     CALL alustasatunnaisgeneraattoriMPI( siemen1, myid)
     
  else
     alusta_satunnais  = .TRUE.
     siemen1 = 0
     
     CALL alustasatunnaisgeneraattori(siemen1)  
     
     CALL alustamuuttujat()

  end if
  
  

  if (exchange_on) then
     if (isfixednode) then
        inquire(file='ref.bead', exist=olemassa)
        if (olemassa) then
           open(123,file='ref.bead')
           read(123,*) refbead_in
           close(123)
           call initRefBead(refbead_in,&
                modulo(refbead_in+hiuk(1)%Trotter/2-2,hiuk(1)%Trotter)+1)
        else
           call initRefBead(int(1,kind=isoluku),hiuk(1)%Trotter/2)
        end if
     end if
  end if

  inquire(file='total_sign.dat', exist=olemassa)
  if (olemassa) then
     open(123,file='total_sign.dat')
     read(123,*) sign_multiply
     close(123)
  else
     sign_multiply=1
  end if
  

  call initObservables()
    
  sgn     = 0
  sgn(1) = sign_multiply

  DO i = 1, hi_lkm(3)
     do j = 1,hiuk(i)%Trotter
        hiuk(i)%upaikka(:,j) = hiuk(i)%vpaikka(:,j)
     end do
  END DO

  if ( myid .eq. 0 ) then
     call printConfBackUp()
  end if
  

  if ( myid .eq. 0 ) then
     CALL DATE_AND_TIME(date = pvm)
     CALL DATE_AND_TIME(time = aloitusaika)
     WRITE( n2, '(I10)' ) blokkeja  
     WRITE( *, * ) 'Date: ', pvm(7:8),'.',pvm(5:6), '.',&
          & pvm(1:4)
     WRITE( *, * ) 'Time: ', aloitusaika(1:2), ':',&
          & aloitusaika(3:4), '.', aloitusaika(5:6)     
  end if



  DO blk = 1,blokkeja
     
     if (modulo(blk,10)==0 .and. myid .eq. 0) then
        call printConfBackUp()
     end if
     
     call BlockInitObservables()     

     ! initiate a table of distances
     call initrtable_PBC()     
     
     ! Do a Monte Carlo block
     CALL mc_PBC( NMax )
     
     ! Print block averages
     call printObservablesMPI( NMax )
     
     if (myid .eq. 0) then
        WRITE( n1, '(I10)' ) blk
        
        teksti = 'Block    ' // trim(adjustl(n1)) &
             // ' / ' // trim(adjustl(n2))
        
        write(*,*) teksti
        CALL DATE_AND_TIME(date = pvm)
        CALL DATE_AND_TIME(time = aloitusaika)
        WRITE( *, * ) '   Date: ', pvm(7:8),'.',pvm(5:6), '.',&
             & pvm(1:4)
        WRITE( *, * ) '   Time: ', aloitusaika(1:2), ':',&
             & aloitusaika(3:4), '.', aloitusaika(5:6)     
     end if

     if (LevelUpdateOn) then
        call UpdateMultilevelL()
     end if
     
          
  END DO
  
  

  
  !----------------------------------------------------
  ! Tyhjataan muisti
  
  call endObservables()
  call deallocatertable()
  !call lopetamuuttujat()
  !call lopetaCSpline2d2()  
  deallocate(hiuk)
  !if (ThermalEstimator) then
  !   call lopeta_dU_dbeta()
  !end if
  IF (exchange_on) THEN
     DEALLOCATE(exchange, STAT = allocstat)
  END IF
  
  call MPI_FINALIZE(ierr)
  stop


  
END PROGRAM PathIntegral

