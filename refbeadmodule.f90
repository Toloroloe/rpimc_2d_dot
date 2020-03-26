module refbeadmodule
  
  use math
  implicit none
  integer(kind = isoluku),private :: refbead,beta2bead
  integer(kind = isoluku),private :: refbeadOld,beta2beadOld
  


contains

  subroutine initRefBead(ref_bead,beta2_bead)
    implicit none
    integer(kind = isoluku), intent(in) :: ref_bead, beta2_bead
    
    refbead = ref_bead
    beta2bead = beta2_bead
    refbeadOld = ref_bead
    beta2beadOld = beta2_bead
    
  end subroutine initRefBead

  subroutine setOldRefBead2New()
    implicit none
    
    refbead = refbeadOld
    beta2bead = beta2beadOld
    
  end subroutine setOldRefBead2New

  subroutine setNewRefBead2Old()
    implicit none
    
    refbeadOld = refbead
    beta2beadOld = beta2bead
    
  end subroutine setNewRefBead2Old

  subroutine setRefBead(ref_bead,beta2_bead)
    implicit none
    integer(kind = isoluku), intent(in) :: ref_bead, beta2_bead
    
    refbead = ref_bead
    beta2bead = beta2_bead
    
  end subroutine setRefBead

  subroutine getRefBead(ref_bead,beta2_bead)
    implicit none
    integer(kind = isoluku), intent(out) :: ref_bead, beta2_bead
    
    ref_bead = refbead
    beta2_bead = beta2bead

  end subroutine getRefBead

end module refbeadmodule
