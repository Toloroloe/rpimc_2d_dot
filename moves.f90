module moves

  use math
  use init
  implicit none  

  logical, private, save :: positivePathFirstTime = .false.

contains 

  SUBROUTINE makeAmove(particle,bead)    
    IMPLICIT NONE
    INTEGER(kind = isoluku), INTENT(out)   :: particle,bead
    REAL(kind=rk)          :: sat
    
    
    ! Which particle??
    CALL RANDOM_NUMBER( sat )
    DO WHILE (sat == 1.)
       CALL RANDOM_NUMBER( sat )
    END DO
    
    particle = CEILING(sat*s_lkm)   
    

    ! Which bead??
    ! -- if multilevel, then this is the
    !    first bead of the "chain" (not moved)
    CALL RANDOM_NUMBER( sat )
    DO WHILE (sat == 1.)
       CALL RANDOM_NUMBER( sat )
    END DO
    
    bead = CEILING(hiuk(particle)%Trotter*sat)
    
    call Move(particle,bead)
    

    
  END SUBROUTINE makeAmove

  subroutine UpdateMovingLabels(start_bead,bead_end)
    use rtable 
    implicit none
    integer(kind = isoluku) :: particle
    integer(kind = isoluku), intent(in) :: start_bead,bead_end
    integer(kind =isoluku) :: hi1p,i,j,hip,hipo

    do particle=1,hi_lkm(1)
       hip=particle
       hipo=particle
       
       i=start_bead
       hiuk(particle)%pn(i)=hip
       hiuk(particle)%po(i)=hipo       
       do while (i .ne. bead_end)
          !if (hiuk(hip)%perm(i) .ne. hip) then
          hip=hiuk(hip)%perm(i)
          !end if
          !if (hiuk(hipo)%perm_old(i) .ne. hipo) then
          hipo=hiuk(hipo)%perm_old(i)
          !end if
          i=modulo(i,hiuk(hip)%Trotter)+1
          hiuk(particle)%pn(i)=hip
          hiuk(particle)%po(i)=hipo
       end do
    end do

  end subroutine UpdateMovingLabels


  subroutine Move(particle,bead)
    use exchangemodule
    use refbeadmodule
    implicit none
    integer(kind = isoluku) :: particle, bead, ref_bead
    real(kind = rk) :: sat
    logical :: positivePath
    
    
    if (hiuk(particle)%ParticleKind .ne. 0) then          
       
       if (isfixednode) then
          positivePath=.false.
          if (.not. positivePathFirstTime) then
             call CheckPath(positivePath)
             if (positivePath) then
                positivePathFirstTime = .true.
             end if
          end if
          if (positivePathFirstTime) then
             CALL RANDOM_NUMBER( sat )
             if (sat<0.1_rk) then                 
                CALL RANDOM_NUMBER( sat )
                ref_bead = floor(sat*real(hiuk(particle)%Trotter,kind=rk))+1
                call setRefBead(ref_bead,&
                     modulo(ref_bead+hiuk(particle)%Trotter/2-1,hiuk(particle)%Trotter)+1)
                call CheckPath2(positivePath)
                if (positivePath) then
                   call setNewRefBead2Old()                   
                else
                   call setOldRefBead2New()
                end if
             end if
             
          end if
       end if
       
       
       CALL multimoveEX_PBC(particle,bead,hiuk(particle)%multi_L)
          
    else
       CALL multimove2_PBC(particle,bead,hiuk(particle)%multi_L)
    end if
    
    
  end subroutine Move


  subroutine swapMoveStart(indices,indCount,start_bead)
    use init
    use rtable 
    implicit none
    integer, intent(in) :: indCount
    integer, dimension(indCount), intent(in) :: indices
    integer(kind = isoluku), intent(in) :: start_bead
    integer(kind=isoluku), dimension(indCount) :: ii
    integer :: i

    do i=1,indCount-1
       ii(i)=hiuk(indices(i+1))%perm(start_bead)
    end do
    ii(indCount)=hiuk(indices(1))%perm(start_bead)
    
    do i=1,indCount
       hiuk(indices(i))%perm(start_bead)=ii(i)
    end do

  end subroutine swapMoveStart


  subroutine BisectionMove_PBC(particle,bead_start,bead2move,bead_end,lambda,LevelTau)
    use rtable 
    implicit none
    integer(kind = isoluku), intent(in) :: particle
    integer(kind = isoluku), intent(in) :: bead_start,bead2move,bead_end
    !integer, intent(in) :: bead_start,bead2move,bead_end
    real(kind = rk), intent(in) :: lambda, LevelTau
    integer(kind =isoluku) :: hi1p,hi2p,i,j,hip

    hi1p=hiuk(particle)%pn(bead_start)
    hip=hiuk(particle)%pn(bead2move)
    hi2p=hiuk(particle)%pn(bead_end)

    hiuk(hip)%upaikka(:,bead2move) = &
         (hiuk(hi1p)%upaikka(:,bead_start)+&
         hiuk(hi2p)%upaikka(:,bead_end))/2 + &
         gaussinen()*sqrt(lambda*LevelTau)
  
    hiuk(hip)%siirretytpaikatlkm=hiuk(hip)%siirretytpaikatlkm+1
    hiuk(hip)%siirretytpaikat(hiuk(hip)%siirretytpaikatlkm)=bead2move
    
  end subroutine BisectionMove_PBC
  




  subroutine BisectionMove2_PBC(particle,bead_start,bead2move,bead_end,lambda,LevelTau,sampleProb)
    use rtable 
    implicit none
    integer(kind = isoluku), intent(in) :: particle
    integer(kind = isoluku), intent(in) :: bead_start,bead2move,bead_end
    !integer, intent(in) :: bead_start,bead2move,bead_end
    real(kind = rk), intent(in) :: lambda, LevelTau
    real(kind = rk) :: sigma2,logT_old,logT_new
    real(kind = rk), intent(inout) :: sampleProb
    real(kind = rk), dimension(dimensio) :: mean_new,mean_old,deltaR
    integer(kind =isoluku) :: hi1p,hi2p,i,j,hip,hi1po,hi2po,hipo
    
    hi1p=hiuk(particle)%pn(bead_start)
    hip=hiuk(particle)%pn(bead2move)
    hi2p=hiuk(particle)%pn(bead_end)
    !hi1po=hiuk(particle)%po(bead_start)
    !hipo=hiuk(particle)%po(bead2move)
    !hi2po=hiuk(particle)%po(bead_end)


    sigma2 = lambda*LevelTau

    mean_new = (hiuk(hi1p)%upaikka(:,bead_start)+&
         hiuk(hi2p)%upaikka(:,bead_end))/2
    
    !mean_old = (hiuk(hi1po)%vpaikka(:,bead_start)+&
    !     hiuk(hi2po)%vpaikka(:,bead_end))/2
    
    deltaR = gaussinen()*sqrt(sigma2)
    
    hiuk(hip)%upaikka(:,bead2move) = mean_new + deltaR    
    
    hiuk(hip)%siirretytpaikatlkm=&
         hiuk(hip)%siirretytpaikatlkm+1
    hiuk(hip)%siirretytpaikat(&
         hiuk(hip)%siirretytpaikatlkm)=bead2move
    
    !logT_old = -sum((mean_old-hiuk(&
    !     hipo)%vpaikka(:,bead2move))**2)/2/sigma2
    !logT_new = -sum((mean_new-hiuk(&
    !     hip)%upaikka(:,bead2move))**2)/2/sigma2
    
    !sampleProb = sampleProb-logT_new+logT_old
    
  end subroutine BisectionMove2_PBC

  
  subroutine multimoveEX_PBC(particle,start_bead,LevelMax)
    use exchangemodule
    use refbeadmodule
    implicit none
    integer, intent(in) :: LevelMax
    integer(kind = isoluku), intent(in) :: particle
    integer(kind = isoluku), intent(in) :: start_bead
    real(kind = rk) :: lambda, LevelTau
    integer(kind = isoluku) :: bead_start,bead2move,bead_end,Trotter
    integer(kind = isoluku) :: ref_bead,diff,beta2_bead,diff2,excycle,ptcl
    integer :: i,j,i_max
    integer, dimension(hi_lkm(3)) :: indices
    integer :: indicesCount
    real(kind = rk) :: dU,dU_old_stage,randnum,ex_new,ex_old,ex_new2,sampleProb
    logical :: checkMove, positivePath,swapped
    real(kind=rk) :: UP, UI,FW,BW,weight,w2

    positivePath=.false.
    checkMove = .false.
    swapped = .false.
    sampleProb=0.0_rk
    indices=0
    indicesCount=0
    
    Trotter = hiuk(particle)%Trotter
    lambda = hiuk(particle)%lambda

    dU_old_stage = 0.0_rk
    dU = 0.0_rk
    
    ex_new = 1.0_rk
    ex_old = 1.0_rk
    
    call getRefBead(ref_bead,beta2_bead)

    excycle=FindExchangeCycle(particle)

    if (hiuk(particle)%ParticleKind .ne. 0) then
       LevelTau = real(2**LevelMax,kind = rk)*tau_in
       bead_end=modulo(start_bead+2**LevelMax-1,Trotter) + 1
       
       call UpdateMovingLabels(start_bead,bead_end)         
       weight = 0.0_rk
       
       do while (indicesCount==0) 
          call PermutationCycle(excycle,hiuk(particle)%lkm,&
               indices(1:hiuk(particle)%lkm),indicesCount,&
               weight,start_bead,bead_end,LevelTau)
       end do
       
       if (indicesCount>1) then          
          call swapMoveStart(indices(1:indicesCount),indicesCount,start_bead)
          swapped = .true.
          !dU=-log(weight) 
          dU=0.0_rk  ! already taken care of
       elseif (indicesCount==1) then
          ! nothing to be swapped, since identity permutation
          swapped = .true.
       end if
    end if
    
    if (.not. swapped) then
       bisection_hyv = .false.
       return
    end if
    
    bead_end=modulo(start_bead+2**LevelMax-1,Trotter) + 1
    call UpdateMovingLabels(start_bead,bead_end)
    
    do i = LevelMax,1,-1
       
       LevelTau = real(2**(i-1),kind = rk)*tau_in
       
       dU_old_stage = dU
       dU = 0.0_rk
       sampleProb = 0.0_rk
       do j = start_bead,start_bead+2**LevelMax-2,2**i
          bead_start = modulo(j-1,Trotter) + 1
          bead2move = modulo(bead_start+2**(i-1)-1,Trotter) + 1
          bead_end = modulo(bead_start+2**i-1,Trotter) + 1
          
          do ptcl=1,indicesCount
             call BisectionMove2_PBC(int(indices(ptcl),kind=isoluku),&
                  bead_start,bead2move,bead_end,lambda,LevelTau,&
                  sampleProb)
          end do
             
          if (isfixednode) then
             diff = abs(bead2move-ref_bead)
             diff = minval((/diff,Trotter-diff/))
             
             if (diff .ne. 0) then   
                ex_new = exchangeAction2_PBC(excycle,ref_bead,bead2move,real(diff,kind=rk)*tau_in)                
                if (ex_new>0.0_rk) then
                   !
                else
                   bisection_hyv = .false.
                   return
                end if
                ex_new = exchangeAction2_PBC(excycle,ref_bead,bead2move,real(diff,kind=rk)*tau_in+tau_in)                
                if (ex_new>0.0_rk) then
                   !
                else
                   bisection_hyv = .false.
                   return
                end if
                ex_new = exchangeAction2_PBC(excycle,ref_bead,bead2move,maxval((/real(diff,kind=rk)*tau_in-tau_in,0.02_rk*tau_in/)))
                if (ex_new>0.0_rk) then
                   !
                else
                   bisection_hyv = .false.
                   return
                end if
                
             else
                
                ! ref slice moved
                checkMove = .true.
             end if
             !if (exchangeAction2(excycle,bead_start,bead2move,LevelTau)<0.0_rk&
             !     .or. exchangeAction2(excycle,bead2move,bead_end,LevelTau)&
             !     <0.0_rk) then
             !   bisection_hyv = .false.
             !   return
             !end if
          end if

          if (indicesCount>1) then
             dU = dU &
                  + calc_stage_dU_All(&
                  bead_start,bead2move,bead_end,LevelTau,i)
          else
             dU = dU &
                  + calc_stage_dU_PBC(int(indices(1),kind=isoluku),&
                  bead_start,bead2move,bead_end,LevelTau,i)
          end if
          NumMovedBeads=NumMovedBeads+1
          MovedBeads(NumMovedBeads)=bead2move
       end do
       
       
       call random_number(randnum)  

       if (exp(sampleProb-dU+dU_old_stage)<randnum) then
          bisection_hyv = .false.
          return
       end if
       
    end do

    
    if (checkMove .and. isfixednode) then       
       call CheckPath22(positivePath,excycle)
       if (.not. positivePath) then          
          bisection_hyv = .false.
          return
       end if
    end if
    
    
    bisection_hyv = .true.
    if (modulo(indicesCount,2)==0 .and. hiuk(particle)%ParticleKind==3) then
       sign_multiply = -1
    else
       sign_multiply = 1
    end if
    

    
  end subroutine multimoveEX_PBC


  subroutine multimove_PBC(particle,start_bead,LevelMax)
    use refbeadmodule
    implicit none
    integer, intent(in) :: LevelMax
    integer(kind = isoluku), intent(in) :: particle
    integer(kind = isoluku), intent(in) :: start_bead
    real(kind = rk) :: lambda, LevelTau
    integer(kind = isoluku) :: bead_start,bead2move,bead_end,Trotter,diff,ref_bead,beta2_bead
    integer :: i,j
    real(kind = rk) :: dU,dU_old_stage,randnum,ex_new,sat,eh_ex

    Trotter = hiuk(particle)%Trotter
    lambda = hiuk(particle)%lambda

    dU_old_stage = 0.0_rk
    dU = 0.0_rk
        
    do i = LevelMax,1,-1
       
       LevelTau = real(2**(i-1),kind = rk)*tau_in
 
       dU_old_stage = dU
       dU = 0.0_rk      
       do j = start_bead,start_bead+2**LevelMax-2,2**i
          bead_start = modulo(j-1,Trotter) + 1
          bead2move = modulo(bead_start+2**(i-1)-1,Trotter) + 1
          bead_end = modulo(bead_start+2**i-1,Trotter) + 1
          
          call BisectionMove_PBC(particle,&
               bead_start,bead2move,bead_end,lambda,LevelTau)
          
          dU = dU &
               + calc_stage_dU_PBC(particle,&
               bead_start,bead2move,bead_end,LevelTau,i)
          NumMovedBeads=NumMovedBeads+1
          MovedBeads(NumMovedBeads)=bead2move
       end do
       
       call random_number(randnum)
       if (exp(-dU+dU_old_stage)<randnum) then
          bisection_hyv = .false.
          return
       end if
       
    end do
    
    bisection_hyv = .true.
    
  end subroutine multimove_PBC

  subroutine multimove2_PBC(particle,start_bead,LevelMax)
    use refbeadmodule
    implicit none
    integer, intent(in) :: LevelMax
    integer(kind = isoluku), intent(in) :: particle
    integer(kind = isoluku), intent(in) :: start_bead
    real(kind = rk) :: lambda, LevelTau
    integer(kind = isoluku) :: bead_start,bead2move,bead_end,Trotter,diff,ref_bead,beta2_bead
    integer :: i,j
    real(kind = rk) :: dU,dU_old_stage,randnum,sat,SampleProb

    Trotter = hiuk(particle)%Trotter
    lambda = hiuk(particle)%lambda

    dU_old_stage = 0.0_rk
    dU = 0.0_rk
        
    do i = LevelMax,1,-1
       
       LevelTau = real(2**(i-1),kind = rk)*tau_in
 
       dU_old_stage = dU       
       dU = 0.0_rk      
       SampleProb = 0.0_rk
       do j = start_bead,start_bead+2**LevelMax-2,2**i
          bead_start = modulo(j-1,Trotter) + 1
          bead2move = modulo(bead_start+2**(i-1)-1,Trotter) + 1
          bead_end = modulo(bead_start+2**i-1,Trotter) + 1
          
          call BisectionMove2_PBC(particle,&
               bead_start,bead2move,bead_end,lambda,LevelTau,SampleProb)
          
          dU = dU &
               + calc_stage_dU_PBC(particle,&
               bead_start,bead2move,bead_end,LevelTau,i)

          NumMovedBeads=NumMovedBeads+1
          MovedBeads(NumMovedBeads)=bead2move
       end do
       

       call random_number(randnum)
       if (exp(SampleProb-dU+dU_old_stage)<randnum) then
          bisection_hyv = .false.
          return
       end if
       
    end do
    
    bisection_hyv = .true.
    
  end subroutine multimove2_PBC

  function calc_stage_dU_All(bead_start,bead2move,bead_end,LevelTau,Level) result(dU)
    use action
    use energies
    implicit none
    integer, intent(in) :: Level
    integer(kind = isoluku), intent(in) :: bead_start,bead2move,bead_end
    real(kind = rk), intent(in) :: LevelTau
    real(kind = rk) :: dU
    integer(kind=isoluku) :: i
    
    if (Level==1) then
       dU = PairEnergyAll(bead_start,bead2move,LevelTau,.true.) &
            - PairEnergyAll(bead_start,bead2move,LevelTau,.false.) &
            + PairEnergyAll(bead2move,bead_end,LevelTau,.true.) &
            - PairenergyAll(bead2move,bead_end,LevelTau,.false.)
       !do i=1,s_lkm
       !   dU = dU+KineticActionNew_PBC(i,bead_start,bead2move,LevelTau)&
       !        -KineticActionOld_PBC(i,bead_start,bead2move,LevelTau)&
       !        +KineticActionNew_PBC(i,bead2move,bead_end,LevelTau)&
       !        -KineticActionOld_PBC(i,bead2move,bead_end,LevelTau)
       !end do
    else
       dU = 0.0_rk
       !do i=1,s_lkm
       !   dU = dU+KineticActionNew_PBC(i,bead_start,bead2move,LevelTau)&
       !        -KineticActionOld_PBC(i,bead_start,bead2move,LevelTau)&
       !        +KineticActionNew_PBC(i,bead2move,bead_end,LevelTau)&
       !        -KineticActionOld_PBC(i,bead2move,bead_end,LevelTau)
       !end do
    end if
    
  end function calc_stage_dU_All


  function calc_stage_dU_PBC(particle,bead_start,bead2move,bead_end,LevelTau,Level) result(dU)
    use action
    implicit none
    integer, intent(in) :: Level
    integer(kind = isoluku), intent(in) :: particle
    integer(kind = isoluku), intent(in) :: bead_start,bead2move,bead_end
    real(kind = rk), intent(in) :: LevelTau
    real(kind = rk) :: dU
    
    if (Level==1) then
       dU = PairActionNew_PBC(particle,bead_start,bead2move,LevelTau) &
            - PairActionOld_PBC(particle,bead_start,bead2move,LevelTau) &
            + PairActionNew_PBC(particle,bead2move,bead_end,LevelTau) &
            - PairActionOld_PBC(particle,bead2move,bead_end,LevelTau)
       !dU = dU+KineticActionNew_PBC(particle,bead_start,bead2move,LevelTau)&
       !     -KineticActionOld_PBC(particle,bead_start,bead2move,LevelTau)&
       !     +KineticActionNew_PBC(particle,bead2move,bead_end,LevelTau)&
       !     -KineticActionOld_PBC(particle,bead2move,bead_end,LevelTau)
    else
       dU = 0.0_rk
       !dU = dU+KineticActionNew_PBC(particle,bead_start,bead2move,LevelTau)&
       !     -KineticActionOld_PBC(particle,bead_start,bead2move,LevelTau)&
       !     +KineticActionNew_PBC(particle,bead2move,bead_end,LevelTau)&
       !     -KineticActionOld_PBC(particle,bead2move,bead_end,LevelTau)
    end if
    
  end function calc_stage_dU_PBC


  function calc_stage_dK_PBC(particle,bead_start,bead2move,bead_end,LevelTau,Level) result(dS)
    use action
    implicit none
    integer, intent(in) :: Level
    integer(kind = isoluku), intent(in) :: particle
    integer(kind = isoluku), intent(in) :: bead_start,bead2move,bead_end
    !integer, intent(in) :: bead_start,bead2move,bead_end
    real(kind = rk), intent(in) :: LevelTau
    real(kind = rk) :: dS
    
    dS = KineticActionNew_PBC(particle,bead_start,bead2move,LevelTau) &
         - KineticActionOld_PBC(particle,bead_start,bead2move,LevelTau) &
         + KineticActionNew_PBC(particle,bead2move,bead_end,LevelTau) &
         - KineticActionOld_PBC(particle,bead2move,bead_end,LevelTau)
    
  end function calc_stage_dK_PBC

  subroutine CheckPath(positivePath)
    use mpi
    use exchangemodule
    use refbeadmodule
    implicit none
    logical, intent(out) :: positivePath
    integer(kind=isoluku) :: ref1,ref2,i,Trotter,diff,diff2,bead,k
    logical, save :: inFirstTime = .true.
    real(kind = rk) :: ex_new,ex_old,ex_new2,ex_new3
    integer :: ind
    integer :: myid, numprocs, ierr
    
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

    

    call getRefBead(ref1,ref2)    
    positivePath = .true. 
    do k=1,size(exchange,1)
       Trotter = exchange(k)%Trotter
       do i=1,Trotter        
          if (inFirstTime) then
             !open(123,file='exTest.dat',position='append')
             diff = abs(i-ref1)          
             diff = minval((/diff,Trotter-diff/))
             
             if (diff .ne. 0) then
                bead = i
                ex_new = exchangeAction2_PBC(k,ref1,&
                     bead,real(diff,kind=rk)*tau_in)  
                ex_new2 = exchangeAction2_PBC(k,ref1,&
                     bead,real(diff,kind=rk)*tau_in+tau_in)  
                ex_new3 = exchangeAction2_PBC(k,ref1,&
                     bead,maxval((/real(diff,kind=rk)*tau_in-tau_in,0.02_rk*tau_in/)))  
               
                if (ex_new>0.0_rk .and. ex_new2>0.0_rk&
                     .and. ex_new3>0.0_rk) then                             
                   exchange(k)%PathSign(i)=1
                   exchange(k)%PathSignOld(i)=1
                else
                   exchange(k)%PathSign(i)=0
                   exchange(k)%PathSignOld(i)=0
                end if
             else
                exchange(k)%PathSign(i)=1
                exchange(k)%PathSignOld(i)=1
             end if
             if (exchange(k)%PathSign(i)==0) then
                positivePath = .false.
             end if
          else
             if (exchange(k)%PathSign(i)==0) then
                positivePath = .false.
                return
             end if
          end if
       end do
       
       if (myid .eq. 0) then
          write(*,*) 'cycle            = ', k
          write(*,*) 'sum(PathSign)    = ', sum(exchange(k)%PathSign)
          write(*,*) 'sum(PathSignOld) = ', sum(exchange(k)%PathSignOld)
          write(*,*) ' '
       end if
    end do
    
    inFirstTime=.false.
    
  end subroutine CheckPath
  
  
  subroutine CheckPath2(positivePath)
    use exchangemodule
    use refbeadmodule
    implicit none
    logical, intent(out) :: positivePath
    integer(kind=isoluku) :: ref1,ref2,i,Trotter,diff,diff2,bead,sumN,sumO,k
    real(kind = rk) :: ex_new,ex_old,ex_new2,ex_new3
    integer :: ind
    
    
    call getRefBead(ref1,ref2)
    
    positivePath = .true.
    
    do k=1,size(exchange,1)
       Trotter = exchange(k)%Trotter
       do i=1,Trotter       
          diff = abs(i-ref1)
          diff = minval((/diff,Trotter-diff/))                          
          
          if (diff .ne. 0) then
             bead = i
             ex_new = exchangeAction2_PBC(k,ref1,&
                  bead,real(diff,kind=rk)*tau_in)
             ex_new2 = exchangeAction2_PBC(k,ref1,&
                  bead,real(diff,kind=rk)*tau_in+tau_in)
             ex_new3 = exchangeAction2_PBC(k,ref1,&
                  bead,maxval((/real(diff,kind=rk)*tau_in-tau_in,0.02_rk*tau_in/)))

             if (ex_new>0.0_rk .and. ex_new2>0.0_rk&
                  .and. ex_new3>0.0_rk) then                
                !exchange(k)%PathSign(i)=0
             else
                positivePath = .false.
                return
                !exchange(k)%PathSign(i)=1
             end if
          else
             !exchange(k)%PathSign(i)=1
          end if
       end do
       
       ! sumN=sum(exchange(k)%PathSign)
       ! sumO=sum(exchange(k)%PathSignOld)
       
       ! if (sumN>=sumO) then 
       !    exchange(k)%PathSignOld = exchange(k)%PathSign
       !    positivePath = .true.
       ! else
       !    exchange(k)%PathSign = exchange(k)%PathSignOld
       !    positivePath = .false.
       !    return
       ! end if
    end do
    
  end subroutine CheckPath2


  subroutine CheckPath22(positivePath,excycle)
    use exchangemodule
    use refbeadmodule
    implicit none
    logical, intent(out) :: positivePath
    integer(kind=isoluku) :: ref1,ref2,i,Trotter,diff,diff2,bead,sumN,sumO,k
    integer(kind=isoluku), intent(in) :: excycle
    real(kind = rk) :: ex_new,ex_old,ex_new2,ex_new3
    integer :: ind
    
    
    call getRefBead(ref1,ref2)
    
    positivePath = .true.
    
    !do k=1,size(exchange,1)
    k=excycle
    Trotter = exchange(k)%Trotter
    do i=1,Trotter       
       diff = abs(i-ref1)
       diff = minval((/diff,Trotter-diff/))                          
       
       if (diff .ne. 0) then
          bead = i
          ex_new = exchangeAction2_PBC(k,ref1,&
               bead,real(diff,kind=rk)*tau_in)
          ex_new2 = exchangeAction2_PBC(k,ref1,&
               bead,real(diff,kind=rk)*tau_in+tau_in)
          ex_new3 = exchangeAction2_PBC(k,ref1,&
               bead,maxval((/real(diff,kind=rk)*tau_in-tau_in,0.02_rk*tau_in/)))
          
          if (ex_new>0.0_rk .and. ex_new2>0.0_rk&
               .and. ex_new3>0.0_rk) then             
             !exchange(k)%PathSign(i)=1             
          else
             positivePath = .false.
             return
             !exchange(k)%PathSign(i)=0
          end if
       else
          !exchange(k)%PathSign(i)=1
       end if
    end do
    
    !sumN=sum(exchange(k)%PathSign)
    !sumO=sum(exchange(k)%PathSignOld)
    
    !if (sumN>=sumO) then 
    !   positivePath = .true.
    !else
    !   positivePath = .false.
    !   return
    !end if
    !end do
    
  end subroutine CheckPath22

end module moves
