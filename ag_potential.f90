module ag_potential
  use math
  implicit none
  real(kind=rk), dimension(:), allocatable, private :: trapCoef
  real(kind=rk), dimension(:), allocatable, private :: depth
  real(kind=rk), private :: latticeConstant
  real(kind=rk), private :: a_in_nm
  real(kind=rk), private :: depth_in_meV
  real(kind=rk), private :: trapCoef_in
  real(kind=rk), dimension(dimensio), private :: simulationCell
  real(kind=rk), dimension(:,:), allocatable, private :: TrapCenter
  integer, private :: NumTraps
  integer, private :: cell_multiply_x,cell_multiply_y
  
  logical, save :: gridORanal=.false.
  
  
  

contains

  subroutine getCellMultiplyers(xi,yi)
    implicit none
    integer, intent(out) :: xi,yi
    
    xi = cell_multiply_x
    yi = cell_multiply_y

  end subroutine getCellMultiplyers

  function AGpotential(particle,bead,newORold) result(Vpot)
    implicit none
    integer(kind=isoluku),intent(in) :: particle,bead
    logical, intent(in) :: newORold
    real(kind = rk) :: Vpot

    if (gridORanal) then
       Vpot=AGgrid(particle,bead,newORold)
    else
       Vpot=AGanalytical(particle,bead,newORold)
    end if
    
  end function AGpotential
  
  function AGanalytical(particle,bead,newORold) result(Vpot)
    use init
    use periodicbc
    implicit none
    integer(kind=isoluku),intent(in) :: particle,bead
    logical, intent(in) :: newORold
    real(kind=rk) :: Vpot
    real(kind=rk), dimension(dimensio) :: re
    integer :: i

    if (newORold) then
       re=hiuk(hiuk(particle)%pn(bead))%upaikka(:,bead)
    else
       re=hiuk(hiuk(particle)%po(bead))%vpaikka(:,bead)
    end if

    Vpot=0.0_rk
    do i=1,NumTraps
       Vpot=Vpot+depth(i)*&
            exp(-(periodic_distance2(re,TrapCenter(i,:))/trapCoef(i)**2)**4);
    end do
    

  end function AGanalytical

  function AGgrid(particle,bead,newORold) result(Vpot)
    implicit none
    integer(kind=isoluku),intent(in) :: particle,bead
    logical, intent(in) :: newORold
    real(kind=rk) :: Vpot
    
  end function AGgrid
  
  subroutine setAGparameters(crystal_mom)
    ! for analytical one
    implicit none
    real(kind=rk), dimension(dimensio), intent(inout) :: crystal_mom
    real(kind=rk), dimension(dimensio) :: dr
    integer :: mult,i,j,k
    logical :: olemassa
    character(len=15) :: helptext

    olemassa = .false.
    inquire(file='ag_params.txt', exist=olemassa)
    if (olemassa) then
       open(123,file='ag_params.txt')
       read(123,*) helptext, a_in_nm
       read(123,*) helptext, depth_in_meV
       read(123,*) helptext, trapCoef_in
       read(123,*) helptext, cell_multiply_x
       read(123,*) helptext, cell_multiply_y
       close(123)
    else
       write(*,*) 'AG potential: using defaults'
       write(*,*) '       a = 150nm            '
       write(*,*) '       depth = 0.6meV       '
       write(*,*) '       R = 0.4a               '
       write(*,*) '       mult = 1               '
       open(123,file='ag_params.txt',status='replace')
       write(123,*) 'a_in_nm ', 150.0_rk
       write(123,*) 'depth_in_meV ', 0.6_rk
       write(123,*) 'trap_scale ', 0.4_rk
       write(123,*) 'cell_multiply_x ', 1
       write(123,*) 'cell_multiply_y ', 1
       close(123)
       a_in_nm = 150.0_rk
       depth_in_meV = 0.6_rk
       trapCoef_in=0.4_rk
       cell_multiply_x=1
       cell_multiply_y=1
    end if


    latticeConstant = a_in_nm*nm2bohr

    simulationCell = (/3.0_rk*latticeConstant, sqrt(3.0_rk)*latticeConstant/)

    crystal_mom(1)=2.0_rk*pi*crystal_mom(1)/simulationCell(1)
    crystal_mom(2)=2.0_rk*pi*crystal_mom(2)/simulationCell(2)

    NumTraps=4*cell_multiply_x*cell_multiply_y
    allocate(TrapCenter(NumTraps,dimensio))
    TrapCenter=0.0_rk
    TrapCenter(1,1:2)=(/latticeConstant/2, latticeConstant/4*sqrt(3.0_rk)/)
    TrapCenter(2,1:2)=(/latticeConstant, latticeConstant/4*3.0_rk*sqrt(3.0_rk)/)
    TrapCenter(3,1:2)=(/latticeConstant*2.0_rk,latticeConstant/4*3.0_rk*sqrt(3.0_rk)/)
    TrapCenter(4,1:2)=(/2.5_rk*latticeConstant,latticeConstant/4*sqrt(3.0_rk)/)

    k=5
    ! x-direction
    do i=0,cell_multiply_x-1
       ! y-direction
       do j=0,cell_multiply_y-1
          if (.not. (i==0 .and. j==0)) then
             dr=(/real(i,kind=rk)*simulationCell(1),&
                  real(j,kind=rk)*simulationCell(2)/)
             TrapCenter(k,1:2)=TrapCenter(1,1:2)+dr
             TrapCenter(k+1,1:2)=TrapCenter(2,1:2)+dr
             TrapCenter(k+2,1:2)=TrapCenter(3,1:2)+dr
             TrapCenter(k+3,1:2)=TrapCenter(4,1:2)+dr
             k=k+4;
          end if
       end do
    end do

    simulationCell(1)=simulationCell(1)*real(cell_multiply_x,kind=rk)
    simulationCell(2)=simulationCell(2)*real(cell_multiply_y,kind=rk)

    open(123,file='CellInfo.txt')
    write(123,*) 'Box ', simulationCell
    write(123,*) 'NumTraps ', NumTraps
    write(123,*) 'Traps at'
    do i=1,NumTraps
       write(123,*) TrapCenter(i,:)
    end do
    close(123)
    
    allocate(depth(NumTraps))
    depth = -depth_in_meV*meV2Hartree    
    
    allocate(trapCoef(NumTraps))
    trapCoef = trapCoef_in*latticeConstant
    !trapCoef = trapCoef_in*150.0_rk*nm2bohr
    ! keep the size of the trap constant, but move the traps   
    
    
  end subroutine setAGparameters

  function getAGcell() result(Cell)
    implicit none
    real(kind=rk), dimension(dimensio) :: Cell
    
    Cell = simulationCell
    
  end function getAGcell



  subroutine writeAGparameters()
    implicit none
    
  end subroutine writeAGparameters
  
  subroutine readAGparameters()
    ! for analytical one
    implicit none
    
    
  end subroutine readAGparameters

  subroutine readGridPotential()
    implicit none

  end subroutine readGridPotential


end module ag_potential
