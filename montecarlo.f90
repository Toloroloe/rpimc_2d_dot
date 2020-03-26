MODULE montecarlo


CONTAINS
  
  SUBROUTINE mc_PBC( NMax )
    
    USE init
    USE observables
    USE rtable
    USE moves

    IMPLICIT NONE


    INTEGER(kind = isoluku), INTENT(IN) :: NMax
    INTEGER  :: s_indx, otosind, i,ii
    integer(kind = isoluku) :: nauha, helmi
    INTEGER :: sgn_old

    !---------------------------------------------------
    !
    
    !----------------------------------------------------
    ! Alustetaan muuttujia
    
    otosind  = 0
    nauha    = 1

    sgn_old = sgn(1)
    !sgn(1) = 1

            
    !---------------------------------------------
    ! Monte Carlo -silmukka alkaa
    !
    
    !s_indx=0
    !DO while (s_indx <= NMax)
    DO s_indx=1,NMax
       bisection_hyv = .FALSE.
       
       !--------------------------------------------
       ! Valitaan siirrettava hiukkanen ja arvotaan uudet paikat
       
       CALL makeAmove(nauha,helmi)
              
       !-------------------------------------------
       ! Metropoliksen vertailu

       !ii=0
       IF ( bisection_hyv ) THEN
          !-------------------------------------------------------------
          !Accepted

          !! This should be used if more than one particle moved at the same time
          do i=1,s_lkm 
             if (hiuk(i)%siirretytpaikatlkm>0) then
                hiuk(i)%vpaikka(:,hiuk(i)%siirretytpaikat(1:hiuk(i)%siirretytpaikatlkm))=&
                     hiuk(i)%upaikka(:,hiuk(i)%siirretytpaikat(1:hiuk(i)%siirretytpaikatlkm))
                
                !Particle acceptance
                BlockAcceptedCount(i) =  BlockAcceptedCount(i) + 1
                ParticleMovedCount(i) =  ParticleMovedCount(i) + 1
                ParticleMovedCount(hi_lkm(3)+1) =  ParticleMovedCount(hi_lkm(3)+1) + 1
                BlockAcceptedCount(hi_lkm(3)+1) =  BlockAcceptedCount(hi_lkm(3)+1) + 1                    
                hiuk(i)%siirretytpaikatlkm=0
                !ii=1
             end if
             hiuk(i)%perm_old = hiuk(i)%perm
          end do
          sgn(1)=sgn(1)*sign_multiply
          sgn_old = sgn(1)

          NumMovedBeads=0
          MovedBeads=0

          
       ELSE
          !Rejected

          !! This should be used if more than one particle moved at the same time
          do i=1,s_lkm
             if (hiuk(i)%siirretytpaikatlkm>0) then
                hiuk(i)%upaikka(:,hiuk(i)%siirretytpaikat(1:hiuk(i)%siirretytpaikatlkm))=&
                     hiuk(i)%vpaikka(:,hiuk(i)%siirretytpaikat(1:hiuk(i)%siirretytpaikatlkm))
                
                hiuk(i)%siirretytpaikatlkm=0
                ParticleMovedCount(i) =  ParticleMovedCount(i) + 1
                ParticleMovedCount(hi_lkm(3)+1) =  ParticleMovedCount(hi_lkm(3)+1) + 1
                !ii=1
             end if
             hiuk(i)%perm = hiuk(i)%perm_old
          end do

          NumMovedBeads=0
          MovedBeads=0
          
          sgn(1) = sgn_old
       END IF              


       !if (ii==1) then
       !   s_indx = s_indx+1
       otosind = otosind + 1
       !end if
       
       
       ! Make a "measurement" every once in a while
       IF ( otosind == otosvali ) THEN
          
          CALL calcObservables()
          otosind = 0          
          
       END  IF
       
       
       
    END DO
    
  END SUBROUTINE mc_PBC
  
END MODULE montecarlo
