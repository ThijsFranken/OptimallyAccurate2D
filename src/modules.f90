module parameters

  implicit none

  
  character(120) :: filename
  integer :: tmpint
  character(80) :: modelname
  
  integer, parameter :: times = 1 ! this can make the dx,dz,dt finer  

  ! switch OPT / CONV
  logical :: optimise

  ! switch video
  logical :: videoornot 

  ! writing strains
  logical :: writingStrain
  integer :: iReceiverStart,iReceiverInterval,nReceiver
  integer :: izReceiverStart
  integer, allocatable :: nrx(:),nrz(:) ! Receiver positions


  integer :: iSourceStart,iSourceInterval,nSource
  integer :: izSourceStart
  integer, allocatable, dimension(:) :: iisx, iisz

  integer :: iSource, iReceiver
  integer :: maxnx,maxnz,maxnt

  double precision, parameter :: pi=3.1415926535897932d0 
  double precision, parameter :: ZERO = 0.d0
    
  
  ! parameters for the gridding
  double precision dt,dx,dz
  ! parameters for the wavefield
  integer nt,nx,nz,it,ist,isx,isz,ix,iz,recl_size
  ! Attention ! nx and nz are modified with absorbing boundaries
  double precision, allocatable, dimension(:,:) :: ux,uz,ux1,ux2,uz1,uz2
  double precision, allocatable, dimension(:,:) :: e1,e2,e3,e4,e5,e6,e7,e8
  double precision, allocatable, dimension(:,:) :: e13,e14,e15,e16,e17,e18,e19,e20
  double precision, allocatable, dimension(:,:) :: f1,f2,f3,f4,f5,f6,f7,f8
  double precision, allocatable, dimension(:,:) :: f13,f14,f15,f16,f17,f18,f19,f20
  double precision, allocatable, dimension(:,:) :: work

  ! for discontinuities
  
  double precision, allocatable, dimension(:,:) :: ee12,ee34,ee56,ee65,ee78,ee87
  double precision, allocatable, dimension(:,:) :: ff12,ff34,ff56,ff65,ff78,ff87

  
  ! parameter for the structure
  
  character(80) :: vpfile, vsfile, rhofile   ! modelname
  double precision, allocatable, dimension(:,:) :: rho,lam,mu,fx,fz,vs,vp

  ! Courant number
  double precision :: cp ! maxvalue of vp
  double precision :: Courant_number

  
  ! parameter for the receiver
  integer :: ir,j
  real, allocatable, dimension(:,:) :: synx,synz
  real, allocatable, dimension(:) :: time
  character(200) :: outfile
 
  
  ! parameter for the waveform
  double precision t
  !parameter for video  
  real, allocatable, dimension(:,:) :: video
  real, allocatable, dimension(:,:) :: snapux,snapuz
  integer :: IT_DISPLAY
 
  integer(2) head(1:120)
  character(80) :: routine

  ! switch C-PML 
  logical, parameter :: USE_PML_XMIN = .true.
  logical, parameter :: USE_PML_XMAX = .true.
  logical, parameter :: USE_PML_YMIN = .true.
  logical, parameter :: USE_PML_YMAX = .true.
  ! thickness of the PML layer in grid points
  integer, parameter :: NPOINTS_PML = 100*times
  double precision, parameter :: CerjanRate = 0.0015
  double precision, allocatable, dimension(:,:) :: weightBC
  ! Cerjan boundary condition
  integer :: lmargin(1:2),rmargin(1:2)
  

  ! Ricker wavelets source
  double precision f0,t0
  !double precision tp,ts


  ! for evolution of total energy in the medium
  double precision epsilon_xx,epsilon_yy,epsilon_xy
  double precision, allocatable, dimension(:) :: total_energy_kinetic,total_energy_potential
  
  ! power to compute d0 profile
  double precision, parameter :: NPOWER = 2.d0

  double precision, parameter :: K_MAX_PML = 1.d0 ! from Gedney page 8.11
  double precision :: ALPHA_MAX_PML

  ! for water
  
  integer, allocatable, dimension(:,:) :: liquidmarkers

  ! for discontinuities

  integer, allocatable, dimension(:,:) :: markers
  integer :: nDiscon ! number of discontinuities
  integer :: lengthDiscon ! with x,z coordinates
  double precision, allocatable :: dscr(:,:,:) ! discontinuity coordinates
  double precision :: tmpvaluex,tmpvaluez

  ! for free surface
  
  integer,allocatable, dimension(:,:) :: zerodisplacement
  integer :: lengthFreeSurface ! with x,z coordinates
  double precision, allocatable :: free(:,:)
 
  

  ! for waveform inversion
  
  
  real(kind(0e0)), allocatable, dimension(:,:):: singleStrainDiagonal,singleStrainShear,tmpsingleStrain

  
  character(340) :: commandline
  

  

end module parameters


module paramFrechet
  
  implicit none
  integer :: i1Source, i2Source ! 2 sources for cross correlations
  integer :: isx1,isx2,isz1,isz2,it1,it2
  real, allocatable, dimension(:,:) :: singleStrainForward,singleStrainBack
  real, allocatable, dimension(:,:) :: singleKernelP,singleKernelS
  double precision, allocatable, dimension (:,:) :: strainForward,strainBack
  double precision, allocatable, dimension (:,:) :: kernelP,kernelS
  
end module paramFrechet
