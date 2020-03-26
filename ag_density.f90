module ag_density
  ! 2d electron density
  use math
  use init
  implicit none
  
  integer(kind=isoluku), dimension(:,:), allocatable,private :: density
  integer(kind=isoluku), private :: densityCount

  real(kind=rk), dimension(:,:), allocatable,private :: densityMPI
  !integer(kind=isoluku), private :: densityCountMPI


  real(kind = rk), dimension(:), allocatable,private :: xgrid
  real(kind = rk), dimension(:), allocatable,private :: ygrid
  integer, private :: xdim
  integer, private :: ydim
  real(kind=rk),private :: dx
  real(kind=rk),private :: dy

  integer(kind = isoluku), private, save :: nblock
  
contains
  
  
  subroutine init_ag_densityMPI()
    use mpi
    implicit none
    !integer, intent(in) :: mult1,mult2
    real(kind = rk) :: limits
    integer :: myid, numprocs, ierr
    
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

    xdim = 100
    ydim = 100

    allocate(xgrid(xdim+1))
    allocate(ygrid(ydim+1))

    limits=sqrt(-2.0_rk/hiuk(1)%meff/sqrt(omega2)*log(5.0e-5_rk))

    xgrid=linspace(-limits,limits,xdim+1)
    ygrid=linspace(-limits,limits,ydim+1)

    dx = xgrid(2)-xgrid(1)
    dy = ygrid(2)-ygrid(1)
    
    allocate(density(xdim,ydim))
    density=0
    densityCount=0
    
    !if (myid .eq. 0) then
    allocate(densityMPI(xdim,ydim))
    densityMPI=0.0_rk
       !densityCountMPI=0
    !end if
    
    nblock = 0

  end subroutine init_ag_densityMPI
  
  subroutine block_init_ag_density()
    implicit none
    
    nblock = nblock + 1

    if (nblock==block_cut) then
       density=0
       densityCount=0
    end if
    
  end subroutine block_init_ag_density


  subroutine end_ag_density()
    implicit none
    
    deallocate(xgrid)
    deallocate(ygrid)
    deallocate(density)

  end subroutine end_ag_density

  subroutine calc_ag_density(sgn)
    implicit none
    integer :: indx,indy,i,j
    integer(kind = isoluku), intent(in) :: sgn
    
    do i=1,hi_lkm(1) ! over electrons (or negative particles) only
       do j=1,hiuk(i)%Trotter
          indx=floor(hiuk(i)%upaikka(1,j)/dx) + xdim/2 + 1
          indy=floor(hiuk(i)%upaikka(2,j)/dy) + ydim/2 + 1
          if (indx>0 .and. indx<=xdim .and. indy>0 .and. indy<=ydim) then
             density(indx,indy)=density(indx,indy)+sgn
          end if
          densityCount=densityCount + sgn
       end do
    end do
    
  end subroutine calc_ag_density
  
  subroutine print_ag_density()
    implicit none
    integer :: i,j

    ! gnuplot format
    open(123,file='ag_density.dat')
    do i=1,xdim
       do j=1,ydim
          write(123,'(3(ES18.7))') xgrid(i)+dx/2,ygrid(j)+dy/2,&
               real(density(i,j),kind=rk)/densityCount
       end do
       write(123,'(A)') ' '
    end do
    close(123)
    
  end subroutine print_ag_density
  

  subroutine print_ag_densityMPI()
    use mpi
    implicit none
    integer :: i,j
    integer :: myid, numprocs, ierr
    
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

    densityMPI=0.0_rk
    call MPI_Reduce(real(density,kind=rk)/densityCount,&
         densityMPI,xdim*ydim,MPI_DOUBLE_PRECISION &
         ,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    ! gnuplot format
    if (myid .eq. 0) then
       open(123,file='ag_density.dat')
       do i=1,xdim
          do j=1,ydim
             write(123,'(3(ES18.7))') xgrid(i)+dx/2,ygrid(j)+dy/2,&
                  densityMPI(i,j)/numprocs
          end do
          write(123,'(A)') ' '
       end do
       close(123)
    end if
    
  end subroutine print_ag_densityMPI
  


end module ag_density
