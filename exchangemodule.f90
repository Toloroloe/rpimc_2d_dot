module exchangemodule

contains
  
  SUBROUTINE make_ex_matrix_PBC(A,he1,he2,jakaja,nauha,lkm,tau)

    USE math
    IMPLICIT NONE
    REAL(KIND=RK), DIMENSION(:,:), INTENT(out) :: A
    REAL(KIND=RK), INTENT(out) :: jakaja
    INTEGER(kind = isoluku), INTENT(in) :: he1,he2
    integer, intent(in) :: nauha,lkm
    INTEGER :: i,j,m,n
    REAL(KIND=RK) :: E,tau

    E = 0.0_rk
    A = 0.0_rk
    jakaja = 1.0_rk
    m = 1
    n = 1

    DO i = nauha,nauha+lkm-1
      IF (m+nauha-1>nauha+lkm-1) THEN
        m = 1
      END IF
      DO j = nauha,nauha+lkm-1
        IF (n+nauha-1>nauha+lkm-1) THEN
          n = 1
        END IF
        CALL T_action_PBC( E,tau,i,j,he1,he2 )
        A(n,m) = EXP( -E )

        IF (i == j) THEN
          jakaja = jakaja*A(n,m)
        END IF
        n = n+1
      END DO
      m = m+1
    END DO

  END SUBROUTINE make_ex_matrix_PBC

  SUBROUTINE make_ex_matrix2(A,he1,he2,jakaja,first,last,lkm,tau)
    USE math
    IMPLICIT NONE
    REAL(KIND=RK), DIMENSION(lkm,lkm), INTENT(out) :: A
    REAL(KIND=RK), INTENT(out) :: jakaja
    INTEGER(kind = isoluku), INTENT(in) :: he1,he2
    integer, intent(in) :: first,last,lkm
    INTEGER :: i,j,m,n
    REAL(KIND=RK) :: E,tau

    E = 0.0_rk
    A = 0.0_rk
    jakaja = 1.0_rk
    m = 1
    n = 1

    DO i = first,last
       n = i-first+1
       DO j = first,last
          m = j-first+1
          CALL T_action_PBC( E,tau,i,j,he1,he2 )
          A(n,m) = EXP( -E )
          
          IF (m == n) THEN
             jakaja = jakaja*A(n,m)
          END IF
       END DO
    END DO
    
  END SUBROUTINE make_ex_matrix2

  SUBROUTINE make_ex_matrix3(A,he1,he2,first,last,lkm,tau)

    USE math
    IMPLICIT NONE
    REAL(KIND=RK), DIMENSION(lkm,lkm), INTENT(out) :: A
    INTEGER(kind = isoluku), INTENT(in) :: he1,he2
    integer, intent(in) :: first,last,lkm
    INTEGER :: i,j,m,n
    REAL(KIND=RK) :: E,tau

    E = 0.0_rk
    A = 0.0_rk
    m = 1
    n = 1

    DO i = first,last
       n = i-first+1
       DO j = first,last
          m = j-first+1
          CALL T_action2( E,tau,i,j,he1,he2 )
          A(n,m) = EXP( -E )          
       END DO
    END DO
    
  END SUBROUTINE make_ex_matrix3
  


  FUNCTION exchangeAction_PBC(he1,he2,tau) RESULT( EX )

    USE math
    use init

    IMPLICIT NONE
    INTEGER(kind = isoluku), INTENT(in) :: he1,he2
    INTEGER :: i
    REAL(KIND=RK) :: EX,jakaja,tau,jakaja1
    
    EX = 1.0_rk
    jakaja1 = 1.0_rk
    jakaja = 1.0_rk

    if (he1==he2) then
       EX=1.0_rk
    else
       
       DO i = 1,SIZE(exchange,1)
          !Oikeastaan tarvitsee muodostaa vain se matriisi,
          !jossa on liikutettu hiukkanen. Muut determinantit
          !kumoutuvat vertailussa (muuttumattomina) pois.
          !CALL make_ex_matrix_PBC(exchange(i)%Matrix,&
          !    he1,he2,jakaja1,&
          !    exchange(i)%nauha,exchange(i)%lkm,tau)
          CALL make_ex_matrix2(exchange(i)%Matrix,&
               he1,he2,jakaja1,exchange(i)%first,exchange(i)%last,&
               exchange(i)%lkm,tau)
          
          EX = det(exchange(i)%Matrix,exchange(i)%lkm)
          if (isnan(EX)) then
             EX=-1.0_rk
          end if
          if (EX<0.0_rk) then
             return
          end if
          jakaja = jakaja*jakaja1
       END DO
       
       EX = EX/jakaja
    end if

!     i=excycle
!     CALL make_ex_matrix_PBC(exchange(i)%Matrix,&
!          he1,he2,exchange(i)%Trotter,jakaja,&
!          exchange(i)%nauha,exchange(i)%lkm,tau)
    
!     EX = det(exchange(i)%Matrix,exchange(i)%lkm)/jakaja


  END FUNCTION exchangeAction_PBC


  
  subroutine PermutationCycle(species,nptcls,indices,indicesCount,weight,bead1,bead2,tau)
    ! http://arxiv.org/pdf/physics/0506020v4.pdf
    USE math
    use init
    IMPLICIT NONE
    integer, intent(in) :: nptcls
    integer, intent(out) :: indicesCount    
    integer, dimension(nptcls),intent(out) :: indices 
    INTEGER(kind = isoluku), INTENT(in) :: species,bead1,bead2
    real(kind=rk), intent(out) :: weight
    real(kind=rk), intent(in) :: tau
    INTEGER :: ptcl, i,ii,first,j
    real(kind=rk) :: randnum,QvalTot,Qval,Prob
    real(kind=rk), dimension(nptcls,nptcls) :: exM,exM2
    logical :: open_cycle

    exM=0.0_rk
    exM2=0.0_rk
    weight=0.0_rk
    indicesCount=0
    open_cycle=.true.
    
    ii=species
    first=exchange(ii)%first
  
    CALL make_ex_matrix3(exM,&
         bead1,bead2,exchange(ii)%first,exchange(ii)%last,&
         nptcls,tau)
    exM2=exM
    
        
    indices=0
    call random_number(randnum)
    ptcl=first+floor(randnum*0.99_rk*real(nptcls,kind=rk))

    indicesCount=indicesCount+1  
    indices(indicesCount)=ptcl
    
    do while (open_cycle) 
       
       ! take away returns to previous particles
       do i=1,indicesCount
          exM2(indices(indicesCount)-first+1,indices(i)-first+1) = 0.0_rk
       end do
       
       ! allow cycle to close
       if (indicesCount>0) then
          exM2(indices(indicesCount)-first+1,indices(1)-first+1)=&
               exM(indices(indicesCount)-first+1,indices(1)-first+1)
       end if
       
       !if (isfixednode) then
       !   ! with restricted pimc even permutations are not allowed to close
       !   if (modulo(indicesCount,2)==0) then
       !      exM2(indices(indicesCount)-first+1,indices(1)-first+1)=0.0_rk
       !   end if
       !end if

       QvalTot=sum(exM(indices(indicesCount)-first+1,1:nptcls))
       Qval=sum(exM2(indices(indicesCount)-first+1,1:nptcls))

       call random_number(randnum)       
       if (Qval/QvalTot<randnum) then
          weight=0.0_rk
          indices=0
          indicesCount=0
          return
       end if

       call random_number(randnum)       
       Prob=0.0_rk
       j=0
       do while (Prob<=randnum .and. j<nptcls)
          j=j+1
          Prob=Prob+exM2(indices(indicesCount)-first+1,j)/Qval
       end do
       ptcl = first+j-1
       if (ptcl==indices(1)) then
          ! cycle is closed
          open_cycle=.false.
       !elseif (indicesCount==nptcls) then
       !   indices=0
       !   indicesCount=0
       !   weight=0.0_rk
       !   return
       else
          indicesCount=indicesCount+1  
          indices(indicesCount)=ptcl
       end if
       
       
    end do
    
    if (isfixednode) then
       ! with restricted pimc even permutations are not allowed to close
       if (modulo(indicesCount,2)==0) then
          indices=0
          indicesCount=0
          weight=0.0_rk
          return
       end if
    end if
    
    weight=1.0_rk
    if (indicesCount>1) then
       do i=1,indicesCount-1
          weight=weight*exM(indices(i)-first+1,indices(i+1)-first+1)/&
               exM(indices(i)-first+1,indices(i)-first+1)
       end do
       weight=weight*exM(indices(indicesCount)-first+1,indices(1)-first+1)/&
            exM(indices(indicesCount)-first+1,indices(indicesCount)-first+1)
    end if

  end subroutine PermutationCycle



  FUNCTION exchangeAction2_PBC(excycle,he1,he2,tau) RESULT( EX )

    USE math
    use init

    IMPLICIT NONE
    INTEGER(kind = isoluku), INTENT(in) :: he1,he2,excycle
    INTEGER(kind=isoluku) :: i
    REAL(KIND=RK) :: EX,jakaja,tau,jakaja1
    
    EX = 1.0_rk
    jakaja1 = 1.0_rk
    jakaja = 1.0_rk   

    i=excycle
    
    CALL make_ex_matrix2(exchange(i)%Matrix,&
         he1,he2,jakaja,exchange(i)%first,exchange(i)%last,&
         exchange(i)%lkm,tau)
    
    EX = det(exchange(i)%Matrix,exchange(i)%lkm)/jakaja

  END FUNCTION exchangeAction2_PBC

  FUNCTION exchangeAction2(excycle,he1,he2,tau) RESULT( EX )

    USE math
    use init

    IMPLICIT NONE
    INTEGER(kind = isoluku), INTENT(in) :: he1,he2,excycle
    INTEGER(kind=isoluku) :: i
    REAL(KIND=RK) :: EX,jakaja,tau,jakaja1
    
    EX = 1.0_rk
    jakaja1 = 1.0_rk
    jakaja = 1.0_rk   

    i=excycle
    
    CALL make_ex_matrix3(exchange(i)%Matrix,&
         he1,he2,exchange(i)%first,exchange(i)%last,&
         exchange(i)%lkm,tau)
    
    EX = det(exchange(i)%Matrix,exchange(i)%lkm)

  END FUNCTION exchangeAction2


  SUBROUTINE T_action2( E, tau, hi1, hi2,he1,he2)
    use math
    use init
    implicit none
    real(kind=rk), intent(out) :: E
    real(kind=rk) :: tau,dr2
    integer :: hi1p,hi2p
    integer, intent(in) :: hi1,hi2
    integer(kind = isoluku) :: he1,he2,i,j


    hi1p=hiuk(hi1)%pn(he1) ! should be equal to hi1 for start bead
    hi2p=hiuk(hi2)%pn(he2)

    E = 0.0_rk
    
    dr2 = sum((hiuk(hi1p)%upaikka(:,he1)-&
         hiuk(hi2p)%upaikka(:,he2))**2)
    E = dr2/(hiuk(hi1p)%lambda*tau)/4

  end SUBROUTINE T_action2

  SUBROUTINE T_action_PBC( E, tau, hi1, hi2,helmi1,helmi2)
    USE math
    use init
    IMPLICIT NONE
    REAL(KIND=RK) :: tau, dr2,ff,image_corr
    REAL(KIND=RK), INTENT(out) :: E
    INTEGER(kind = isoluku) :: helmi1,helmi2,Trotter
    real(kind=rk), dimension(2) :: ddr,Box,r1v,r2v
    integer :: hi1,hi2,hi1p,hi2p,hip,i
    real(kind=rk) :: oo, r1,r2
    
    E = 0.0_rk
    
    hi1p=hi1
    hi2p=hi2
    
    dr2 = sum((hiuk(hi1p)%upaikka(:,helmi1)-&
         hiuk(hi2p)%upaikka(:,helmi2))**2)
    E = dr2/(hiuk(hi1p)%lambda*tau)/4
    

  END SUBROUTINE T_action_PBC

  function FindExchangeCycle(particle) result(excycle)
    use math
    use init
    implicit none
    integer(kind = isoluku) :: particle, excycle,i,j
    !integer :: i,j
    
    excycle=1
    do i=1,size(exchange,1)
       if (particle>=exchange(i)%first .and. &
            particle<=exchange(i)%last) then
          excycle=i
          return
       end if
    end do
    

  end function FindExchangeCycle

  
end module exchangemodule
