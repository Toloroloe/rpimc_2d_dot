module Energies

contains

  function PairEnergy(hi,he1,he2,tau,newORold) result(V)
    use math
    use init
    implicit none
    integer(kind = isoluku) :: j,kk,hi,he1,he2,jj
    real(kind=rk) :: r1,r2,ZZ,V1,tau
    real(kind=rk) :: V,s,x,y,Vnodal
    real(kind=rk), dimension(dimensio) :: r1v,r2v,re
    logical :: newORold
    
    s = 0.0_rk
    r1v = 0.0_rk
    r2v = 0.0_rk
    r1 = 0.0_rk
    r2 = 0.0_rk
    V = 0.0_rk
    V1 = 0.0_rk
    j = 1 
    jj = 1
    
    
    do kk = 1,hi_lkm(3)
       if (hi .ne. kk) then
          
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

          if (newORold) then
             r1 = sqrt(sum((hiuk(hiuk(hi)%pn(he1))%upaikka(:,he1)-&
                  hiuk(hiuk(kk)%pn(j))%upaikka(:,j))**2))
             r2 = sqrt(sum((hiuk(hiuk(hi)%pn(he2))%upaikka(:,he2)-&
                  hiuk(hiuk(kk)%pn(jj))%upaikka(:,jj))**2))
             
             
          else
             r1 = sqrt(sum((hiuk(hiuk(hi)%po(he1))%vpaikka(:,he1)-&
                  hiuk(hiuk(kk)%po(j))%vpaikka(:,j))**2))
             r2 = sqrt(sum((hiuk(hiuk(hi)%po(he2))%vpaikka(:,he2)-&
                  hiuk(hiuk(kk)%po(jj))%vpaikka(:,jj))**2))
          
          end if
          

          if ( ZZ < 0.0) then          
             ! should not come here with this code
          else
             ! Just Coulombic 1/r-type repulsion
             V1 = primActionOld(r1,r2,tau,ZZ)
          end if
          
          V = V + V1
          V1 = 0.0_rk
       end if
    end do
    
    if (newORold) then
       re=hiuk(hiuk(hi)%pn(he1))%upaikka(:,he1)
    else
       re=hiuk(hiuk(hi)%po(he1))%vpaikka(:,he1)
    end if
    
    V=V+hiuk(hi)%meff*omega2*sum(re**2)/2*tau
    
    if (.not. NoMagneticField) then
       V=V+hiuk(hi)%meff*omegaL2*sum(re**2)/2*tau
    end if

    
  end function PairEnergy

  function PairEnergyAll(he1,he2,tau,newORold) result(V)
    use math
    use init
    implicit none
    integer(kind = isoluku) :: j,kk,hi,he1,he2,jj,i
    real(kind=rk) :: r1,r2,ZZ,V1,tau
    real(kind=rk) :: V,s,x,y,Vnodal
    real(kind=rk), dimension(dimensio) :: r1v,r2v,re
    logical :: newORold
    
    s = 0.0_rk
    r1v = 0.0_rk
    r2v = 0.0_rk
    r1 = 0.0_rk
    r2 = 0.0_rk
    V = 0.0_rk
    V1 = 0.0_rk
    j = 1 
    jj = 1
    
    do hi=1,hi_lkm(3)-1
       do kk = hi+1,hi_lkm(3)
      
          
          ! Trotter numbers must be the same
          ! for all moving particles
          if (hiuk(kk)%Trotter == 1) then
             j = 1
             jj = 1
          else
             j = he1
             jj = he2
          end if
          
          if (hiuk(hiuk(hi)%pn(he1))%siirretytpaikatlkm>0 .or. &
               hiuk(hiuk(kk)%pn(j))%siirretytpaikatlkm>0 .or. &
               hiuk(hiuk(hi)%pn(he2))%siirretytpaikatlkm>0 .or. &
               hiuk(hiuk(kk)%pn(jj))%siirretytpaikatlkm>0) then
             
             
             ZZ = hiuk(hi)%Z * hiuk(kk)%Z
             
             
             if (newORold) then
                r1 = sqrt(sum((hiuk(hiuk(hi)%pn(he1))%upaikka(:,he1)-&
                     hiuk(hiuk(kk)%pn(j))%upaikka(:,j))**2))
                r2 = sqrt(sum((hiuk(hiuk(hi)%pn(he2))%upaikka(:,he2)-&
                     hiuk(hiuk(kk)%pn(jj))%upaikka(:,jj))**2))
                
                
             else
                r1 = sqrt(sum((hiuk(hiuk(hi)%po(he1))%vpaikka(:,he1)-&
                     hiuk(hiuk(kk)%po(j))%vpaikka(:,j))**2))
                r2 = sqrt(sum((hiuk(hiuk(hi)%po(he2))%vpaikka(:,he2)-&
                     hiuk(hiuk(kk)%po(jj))%vpaikka(:,jj))**2))
                
             end if

             
             if ( ZZ < 0.0) then          
                ! should not come here
             else            
                V1 = primActionOld(r1,r2,tau,ZZ)
             end if
             
             V = V + V1
             V1 = 0.0_rk
          end if
       end do
    end do
       
       
    do i=1,hi_lkm(3)    
       if (newORold) then
          re=hiuk(i)%upaikka(:,he1)
       else
          re=hiuk(i)%vpaikka(:,he1)
       end if
       V=V+hiuk(i)%meff*omega2*sum(re**2)/2*tau
       if (.not. NoMagneticField) then
          V=V+hiuk(hi)%meff*omegaL2*sum(re**2)/2*tau
       end if
    end do
    
  end function PairEnergyAll

  
  function primActionOld(r,rp,tau,ZZ) result(V)
    use math
    use init
    implicit none
    real(kind = rk) :: r,rp,tau,ZZ,V,x,y

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
    
    V = 0.5_rk*tau*ZZ*(1.0_rk/x + 1.0_rk/y)/dielectric_constant
    
  end function primActionOld


end module Energies
