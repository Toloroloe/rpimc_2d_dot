module Action

contains

  function KineticActionNew_PBC(hi,he1,he2,tau) result(T)
    use math
    use init
    implicit none
    real(kind=rk) :: T,tau,dr2
    !integer :: hi
    integer(kind = isoluku) :: hi,he1,he2,hi1p,hi2p

    hi1p=hiuk(hi)%pn(he1)
    hi2p=hiuk(hi)%pn(he2)

    T = 0.0_rk
    dr2 = sum((hiuk(hi1p)%upaikka(:,he1)-hiuk(hi2p)%upaikka(:,he2))**2)
    T = dr2/(hiuk(hi)%lambda*tau)/4
    
  end function KineticActionNew_PBC

  function KineticActionOld_PBC(hi,he1,he2,tau) result(T)
    use math
    use init
    implicit none
    real(kind=rk) :: T,tau,dr2
    !integer :: hi
    integer(kind = isoluku) :: hi,he1,he2,hi1p,hi2p

    hi1p=hiuk(hi)%po(he1)
    hi2p=hiuk(hi)%po(he2)

    T = 0.0_rk
    dr2 = sum((hiuk(hi1p)%vpaikka(:,he1)-hiuk(hi2p)%vpaikka(:,he2))**2)
    T = dr2/(hiuk(hi)%lambda*tau)/4
    
  end function KineticActionOld_PBC

  function PairActionNew_PBC(hi,he1,he2,tau) result(U)
    use math
    use init
    use energies
    implicit none
    real(kind=rk) :: U,tau
    !integer :: hi
    integer(kind = isoluku) :: hi,he1,he2

    U = PairEnergy(hi,he1,he2,tau,.true.)    
    
  end function PairActionNew_PBC

  function PairActionOld_PBC(hi,he1,he2,tau) result(U)
    use math
    use init
    use energies
    implicit none
    real(kind=rk) :: U,tau
    !integer :: hi
    integer(kind = isoluku) :: hi,he1,he2

    U = PairEnergy(hi,he1,he2,tau,.false.)    
    
  end function PairActionOld_PBC

end module Action

