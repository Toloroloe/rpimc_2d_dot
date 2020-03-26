MODULE rtable
  use math
  IMPLICIT NONE
  REAL(KIND=RK),DIMENSION(:,:,:), ALLOCATABLE,PRIVATE::drtable
  REAL(KIND=RK),DIMENSION(:,:,:), ALLOCATABLE,PRIVATE::oneperdrtable
  REAL(KIND=RK),PARAMETER :: oneperlimit=1.0e-10_rk

  INTEGER,DIMENSION(:,:,:), ALLOCATABLE,PRIVATE::jtable
CONTAINS

  subroutine deallocatertable()
    implicit none

    deallocate(drtable)
    deallocate(oneperdrtable)
    deallocate(jtable)

  end subroutine deallocatertable


  !this initializes the dr tale based on the old positions
  SUBROUTINE initrtable() 
    USE math
    use init
    LOGICAL,SAVE::firsttime=.TRUE.
    INTEGER(kind = isoluku) :: i,j,k,kk,Trot,maxTrot
    REAL(KIND=RK) :: dr2
    INTEGER               :: allocstat
    
    IF (firsttime) THEN
       WRITE(*,*) "INITIALIZING rtable"
       
       maxTrot=1
       DO k = 1,hi_lkm(3)
          IF(maxTrot<hiuk(k)%Trotter) THEN
             maxTrot=hiuk(k)%Trotter 
          END IF
       END DO
       
       ALLOCATE(drtable(maxTrot,hi_lkm(3),hi_lkm(3)), STAT = allocstat)
       IF ( allocstat /= 0 ) THEN
          WRITE(*,*) 'Virhe muistinvarauksessa initdr()'
          STOP
       END IF
       ALLOCATE(oneperdrtable(maxTrot,hi_lkm(3),hi_lkm(3)), STAT = allocstat)
       IF ( allocstat /= 0 ) THEN
          WRITE(*,*) 'Virhe muistinvarauksessa initdr()'
          STOP
       END IF
       ALLOCATE(jtable(maxTrot,hi_lkm(3),hi_lkm(3)), STAT = allocstat)
       IF ( allocstat /= 0 ) THEN
          WRITE(*,*) 'Virhe muistinvarauksessa initdr()'
          STOP
       END IF
       firsttime=.FALSE.

    END IF

    drtable = 0.0_rk
    oneperdrtable = 1.0_rk
    jtable = 1

    DO k = 1,hi_lkm(3)-1
       DO kk = k+1,hi_lkm(3)
          IF (k .NE. kk) THEN
             j = 1 
             Trot = hiuk(k)%Trotter / hiuk(kk)%Trotter
             !          write(*,*) 'Trot ',Trot
             !          pause
             DO  i = 1, hiuk(k)%Trotter
                IF (i > j*Trot) THEN
                   j = j + 1
                END IF
                
                dr2 = (hiuk(k)%vpaikka(1,i) - hiuk(kk)%vpaikka(1,j))**2 + &
                     (hiuk(k)%vpaikka(2,i) - hiuk(kk)%vpaikka(2,j))**2 + &
                     (hiuk(k)%vpaikka(3,i) - hiuk(kk)%vpaikka(3,j))**2
                drtable(i,k,kk) = SQRT(dr2)
                !            write(*,*) i,j,dr2
                IF(drtable(i,k,kk)<oneperlimit) THEN
                   oneperdrtable(i,k,kk) = 1.0_rk/oneperlimit
                ELSE
                   oneperdrtable(i,k,kk) = (1.0_rk/drtable(i,k,kk))
                END IF
                
                jtable(i,k,kk)=j

                drtable(i,kk,k) = drtable(i,k,kk)
                oneperdrtable(i,kk,k) = oneperdrtable(i,k,kk)

                jtable(i,kk,k)=jtable(i,k,kk)
                
             END DO
             !         pause
          END IF
       END DO
    END DO

  END SUBROUTINE initrtable

    !if bead he for particle hi has moved then update only thos dr values which need 
  !updating
  SUBROUTINE updatedr(hi,he)
    USE math
    use init
    INTEGER(kind = isoluku),INTENT(in):: hi
    INTEGER(kind = isoluku), intent(in) ::he
    INTEGER(kind = isoluku) :: i,j,k,kk,ii,jj,ind,ind2,ind3
    REAL(KIND=RK) :: dr2,Trot2

    k=hi
    i=he
    
    DO kk = 1,hi_lkm(3)
       IF (k .NE. kk) THEN
          Trot2 = real(hiuk(k)%Trotter,kind=rk) &
               / real(hiuk(kk)%Trotter,kind=rk)
          
          if (Trot2<1.0) then
             ii = ceiling(real(i,kind=rk)/Trot2)-1
             jj = ceiling(real(i,kind=rk)/Trot2) 
             ind2 = kk
             ind3 = k
          else
             ii = ceiling(real(i,kind=rk)/Trot2)
             jj = ii
             ind2 = k
             ind3 = kk
          end if
          
                
          ind = i
          do j = ii,jj
             if (Trot2<1.0) then
                ind = j
             end if
             dr2 = (hiuk(k)%upaikka(1,i) - hiuk(kk)%upaikka(1,j))**2 + &
                  (hiuk(k)%upaikka(2,i) - hiuk(kk)%upaikka(2,j))**2 + &
                  (hiuk(k)%upaikka(3,i) - hiuk(kk)%upaikka(3,j))**2
             
             drtable(ind,ind2,ind3) = SQRT(dr2)
!             write(*,*) ind,kk,k,drtable(ind,kk,k)
!             pause
             IF(drtable(ind,ind2,ind3)<oneperlimit) THEN
                oneperdrtable(ind,ind2,ind3) = 1.0_rk/oneperlimit
             ELSE
                oneperdrtable(ind,ind2,ind3) = (1.0_rk/drtable(ind,ind2,ind3))
             END IF
             
             drtable(ind,ind3,ind2) = drtable(ind,ind2,ind3)
             oneperdrtable(ind,ind3,ind2) = oneperdrtable(ind,ind2,ind3)
             
          end do
          
       END IF
    END DO


  END SUBROUTINE updatedr



  !this initializes the dr tale based on the old positions
  SUBROUTINE initrtable_PBC() 
    USE math
    use init
    use mpi
    LOGICAL,SAVE::firsttime=.TRUE.
    INTEGER(kind = isoluku) :: i,j,k,kk,Trot,maxTrot
    INTEGER               :: allocstat
    integer :: myid, numprocs, ierr
    
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
    
    IF (firsttime) THEN
       if (myid .eq. 0) then
          WRITE(*,*) "INITIALIZING rtable"
       end if
       
       maxTrot=1
       DO k = 1,hi_lkm(3)
          IF(maxTrot<hiuk(k)%Trotter) THEN
             maxTrot=hiuk(k)%Trotter 
          END IF
       END DO
       
       ALLOCATE(drtable(maxTrot,hi_lkm(3),hi_lkm(3)), STAT = allocstat)
       IF ( allocstat /= 0 ) THEN
          WRITE(*,*) 'Virhe muistinvarauksessa initdr()'
          STOP
       END IF
       ALLOCATE(oneperdrtable(maxTrot,hi_lkm(3),hi_lkm(3)), STAT = allocstat)
       IF ( allocstat /= 0 ) THEN
          WRITE(*,*) 'Virhe muistinvarauksessa initdr()'
          STOP
       END IF
       ALLOCATE(jtable(maxTrot,hi_lkm(3),hi_lkm(3)), STAT = allocstat)
       IF ( allocstat /= 0 ) THEN
          WRITE(*,*) 'Virhe muistinvarauksessa initdr()'
          STOP
       END IF
       firsttime=.FALSE.

    END IF

    drtable = 0.0_rk
    oneperdrtable = 1.0_rk
    jtable = 1

    DO k = 1,hi_lkm(3)-1
       DO kk = k+1,hi_lkm(3)
          IF (k .NE. kk) THEN
             j = 1 
             Trot = hiuk(k)%Trotter / hiuk(kk)%Trotter
             !          write(*,*) 'Trot ',Trot
             !          pause
             DO  i = 1, hiuk(k)%Trotter
                IF (i > j*Trot) THEN
                   j = j + 1
                END IF
                
                drtable(i,k,kk) = sqrt(sum((hiuk(k)%vpaikka(:,i)- &
                     hiuk(kk)%vpaikka(:,j))**2))
                
                !            write(*,*) i,j,dr2
                IF(drtable(i,k,kk)<oneperlimit) THEN
                   oneperdrtable(i,k,kk) = 1.0_rk/oneperlimit
                ELSE
                   oneperdrtable(i,k,kk) = (1.0_rk/drtable(i,k,kk))
                END IF
                
                jtable(i,k,kk)=j

                drtable(i,kk,k) = drtable(i,k,kk)
                oneperdrtable(i,kk,k) = oneperdrtable(i,k,kk)

                jtable(i,kk,k)=jtable(i,k,kk)
                
             END DO
             !         pause
          END IF
       END DO
    END DO


  END SUBROUTINE initrtable_PBC

    !if bead he for particle hi has moved then update only thos dr values which need 
  !updating
  SUBROUTINE updatedr_PBC(hi,he)
    USE math
    use init
    INTEGER(kind = isoluku),INTENT(in):: hi
    INTEGER(kind = isoluku), intent(in) ::he
    INTEGER(kind = isoluku) :: i,j,k,kk,ii,jj,ind,ind2,ind3
    REAL(KIND=RK) :: Trot2

    k=hi
    i=he
    
    DO kk = 1,hi_lkm(3)
       IF (k .NE. kk) THEN
          Trot2 = real(hiuk(k)%Trotter,kind=rk) &
               / real(hiuk(kk)%Trotter,kind=rk)
          
          if (Trot2<1.0) then
             ii = ceiling(real(i,kind=rk)/Trot2)-1
             jj = ceiling(real(i,kind=rk)/Trot2) 
             ind2 = kk
             ind3 = k
          else
             ii = ceiling(real(i,kind=rk)/Trot2)
             jj = ii
             ind2 = k
             ind3 = kk
          end if
          
                
          ind = i
          do j = ii,jj
             if (Trot2<1.0) then
                ind = j
             end if
 
             drtable(ind,ind2,ind3) = sqrt(sum((hiuk(k)%upaikka(:,i)-&
                  hiuk(kk)%upaikka(:,j))**2))
!             write(*,*) ind,kk,k,drtable(ind,kk,k)
!             pause
             IF(drtable(ind,ind2,ind3)<oneperlimit) THEN
                oneperdrtable(ind,ind2,ind3) = 1.0_rk/oneperlimit
             ELSE
                oneperdrtable(ind,ind2,ind3) = (1.0_rk/drtable(ind,ind2,ind3))
             END IF
             
             drtable(ind,ind3,ind2) = drtable(ind,ind2,ind3)
             oneperdrtable(ind,ind3,ind2) = oneperdrtable(ind,ind2,ind3)
             
          end do
          
       END IF
    END DO


  END SUBROUTINE updatedr_PBC


  PURE FUNCTION evaloneperdr(i,kk,k) RESULT(dr)
    use math
    REAL(KIND=RK)::dr
    INTEGER(kind = isoluku),INTENT(in):: i,kk,k
    dr=oneperdrtable(i,kk,k)
  END FUNCTION evaloneperdr
    
  PURE FUNCTION evaldr(i,kk,k) RESULT(dr)
    use math
    REAL(KIND=RK)::dr
    INTEGER(kind = isoluku),INTENT(in):: i,kk,k
    dr=drtable(i,kk,k)
  END FUNCTION evaldr
  
  PURE FUNCTION getj(i,kk,k) RESULT(jt)
    use math
    integer(kind = isoluku) :: jt
    INTEGER(kind = isoluku),INTENT(in):: i,kk,k
    jt=jtable(i,kk,k)
  END FUNCTION getj

  FUNCTION evaldrdebug(j,i,kk,k) RESULT(dr)
    USE math
    use init
    REAL(KIND=RK)::dr,dr2true,drtrue
    INTEGER:: i,j,kk,k
    dr=drtable(i,kk,k)

    dr2true = (hiuk(k)%upaikka(1,i) - hiuk(kk)%upaikka(1,j))**2 + &
        & (hiuk(k)%upaikka(2,i) - hiuk(kk)%upaikka(2,j))**2 + &
        & (hiuk(k)%upaikka(3,i) - hiuk(kk)%upaikka(3,j))**2
    drtrue = SQRT(dr2true)
    IF(ABS(dr-drtrue)>1e-10) THEN
      WRITE(*,*) "WARNING evaldrdebug: ERROR in dr table"
      WRITE(*,*) "WARNING evaldrdebug: dr=",dr,"drtrue=",drtrue
      WRITE(*,*) "WARNING evaldrdebug: jtable=",jtable(i,kk,k),"j,i,kk,k=",j,i,kk,k
      dr=drtrue
    END IF
    
    
  END FUNCTION evaldrdebug
    
END MODULE rtable
