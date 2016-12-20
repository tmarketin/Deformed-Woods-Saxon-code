module params

  type state ! basis state
    integer :: n,nz,nperp,nrho,l,lambda,sigma,omega
    character(1) :: parity
  end type state
  
  type estate ! eigenstate
    integer :: omega,bdim ! angular momentum projection and basis dimension
    real(kind=8) :: energy,uu,vv,vvtilde
    real(kind=8), dimension(50) :: vector
    character(1) :: parity
  end type estate

! program parameters
  integer, parameter :: lout = 11 ! general output file
  integer, parameter :: lin = 12 ! general input file
  integer, parameter :: ltmp = 13 ! temporary i/o file, used only inside subroutines
  integer, parameter :: lplot = 14
  integer, parameter :: lplot2 = 15
  integer, parameter :: lenergy = 16 ! binding energy

! numerical constants
  real(kind=8), parameter :: epsint = 1e-4 ! integration precision
  real(kind=8), parameter :: nreps = 1e-10 ! numerical recipes precision

! mathematical constants
  real(kind=8), parameter :: half = 0.5
  real(kind=8), parameter :: third = 1.0/3.0
  real(kind=8), parameter :: pi = 3.14159262
  real(kind=8), parameter :: sq54pi = sqrt(5/(4*pi))
  real(kind=8), parameter :: pi4 = 0.7511255444649425
  
! physical constants
  real(kind=8), parameter :: hqc = 197.328284 ! hbar*c in MeV * fm
  real(kind=8), parameter :: alphai = 137.03602 

! nuclear physics constants  
  real(kind=8), parameter :: nucleon_mass = 939.0 ! in MeV
  integer, parameter :: itn = 1 ! neutron isospin index
  integer, parameter :: itp = 2 ! proton isospin index
  
! nucleus parameters
  integer, parameter :: n_proton = 38
  integer, parameter :: n_neutron = 42
  integer, parameter :: n_mass = n_proton + n_neutron

! parameters for Hermite polynomials  
  integer, parameter :: nzmax = 18 ! 12
  integer, parameter :: nrzmax = 1000
  real(kind=8), parameter :: rzmin = -20.0 ! in fm
  real(kind=8), parameter :: rzmax = 20.0 ! in fm
  real(kind=8), parameter :: rzstep = (rzmax-rzmin)/nrzmax ! in fm
  
! parameters for Laguerre polynomials
  integer, parameter :: nrmax = 18 ! 12
  integer, parameter :: nlmax = 18 ! 12
  integer, parameter :: nrrmax = 1000
  integer, parameter :: npol = (nrmax+1)*(nlmax+1) ! number of polynomials
  real(kind=8), parameter :: rrmin = 0.0 ! in fm
  real(kind=8), parameter ::  rrmax = 25.0 ! in fm
  real(kind=8), parameter :: rrstep = rrmax/nrrmax ! in fm
  
! oscillator basis parameters
  integer, parameter :: nmax = 12 ! highest oscillator shell 0,1,2,...,nmax
  integer, parameter :: omegamax = 2*nmax + 1 ! maximum value of omega
  integer, parameter :: nomega = (omegamax+1)/2 ! number of omega states
  integer, parameter :: nstatemax = (nmax + 1)*(nmax + 2) ! maximum number of states of a particular omega
  
! Woods-Saxon potential parameters
! universal parameters  
!  real(kind=8), dimension(2), parameter :: r_c = (/1.347,1.275/) ! in fm
!  real(kind=8), dimension(2), parameter :: a_c = (/0.7,0.7/)  ! in fm
!  real(kind=8), parameter :: v_c0 = -49.6 ! in MeV
!  real(kind=8), parameter :: kappa = 0.86  
!  real(kind=8), parameter :: lambda_so = 36.0 ! Woods-Saxon spin-orbit strength
!  real(kind=8), dimension(2), parameter :: r_so = (/1.32,1.32/) ! in fm
!  real(kind=8), dimension(2), parameter :: a_so = (/0.7,0.7/)  ! in fm
!  real(kind=8), parameter :: r0 = 1.20 ! in fm
! parameters from Isakov
  real(kind=8), dimension(2), parameter :: r_c = (/1.260,1.260/) ! in fm
  real(kind=8), dimension(2), parameter :: a_c = (/0.662,0.662/)  ! in fm
  real(kind=8), parameter :: v_c0 = -52.06 ! in MeV
  real(kind=8), parameter :: kappa = 0.639 
  real(kind=8), parameter :: lambda_so = 24.1 ! Woods-Saxon spin-orbit strength
  real(kind=8), dimension(2), parameter :: r_so = (/1.160,1.160/) ! in fm
  real(kind=8), dimension(2), parameter :: a_so = (/0.662,0.662/)  ! in fm
  real(kind=8), parameter :: r0 = 1.2 ! in fm
  
! Nilsson model spin-orbit parameters
  real(kind=8), parameter :: kappa_so = 0.08 ! nilsson spin-orbit strength
  
! Gauss integration parameters
  integer, parameter :: ngauss_z = 32 ! number of Gauss points in z coordinate integration
  integer, parameter :: ngauss_r = 32 ! number of Gauss points in radial integration
  
! Liquid drop energy parameters
  real(kind=8), parameter :: av = -15.68
  real(kind=8), parameter :: as = 18.56
  real(kind=8), parameter :: ac = 0.717  
  real(kind=8), parameter :: ai = 28.1
  real(kind=8), parameter :: di = 34.0
  
! *** variables ***

! oscillator basis parameters
  real(kind=8) :: hbar_omega0,hbar_omegaz,hbar_omegar
  real(kind=8) :: b0,bz,br
  type(state), dimension(nstatemax,omegamax) :: basis
  integer, dimension(nomega) :: nbasis ! number of states for each omega

! wavefunctions
  real(kind=8), dimension(0:nzmax,ngauss_z) :: wfzg,dwfzg,ddwfzg ! in gauss-hermite points of rzbz
  real(kind=8), dimension(0:nrmax,0:nlmax,ngauss_r) :: wfrg,dwfrg,ddwfrg ! in gauss-laguerre points of rrbr
  
! coordinates
  real(kind=8), dimension(0:nrzmax) :: rz,rzbz ! rzbz = rz/bz  
  real(kind=8), dimension(0:nrrmax) :: rr,rrbr ! rrbr = (rr/br)**2
  
! energies and states
  real(kind=8), dimension(2,nomega,nstatemax) :: energy ! eigenvalues(isospin,omega block,state)  
  real(kind=8), dimension(2,nomega,nstatemax,nstatemax) :: evectors ! eigenvectors(isospin,omega block,state,vector)
  real(kind=8), dimension(2,nomega,nstatemax,nstatemax) :: ham ! hamiltonian for diagonalization
  type(estate), dimension(2,nomega*nstatemax) :: energy_full ! all eigenstates

! Woods-Saxon potential
  real(kind=8), dimension(2) :: v_c ! depth of proton and neutron potentials
  real(kind=8), dimension(2) :: v_so ! strength of neutron and proton spin-orbit potentials
  real(kind=8) :: sfactor ! scaling factor c(beta)
  real(kind=8), dimension(ngauss_r,ngauss_z) :: wspotential ! Woods-Saxon potential in Gauss meshpoints
  
! Coulomb potential  
  real(kind=8) :: cdens ! charge density constant
  real(kind=8), dimension(ngauss_r,ngauss_z) :: coulomb ! Coulomb potential
  
! deformation parameter
  real(kind=8) :: beta
  
! binding energy and pairing
  real(kind=8) :: benergy  
  real(kind=8), dimension(2) :: lambda,lambdat ! chemical potentials
  real(kind=8), dimension(2) :: delta ! pairing gap

! Gauss integration variables
  real(kind=8), dimension(ngauss_z) :: gxz
  real(kind=8), dimension(ngauss_z,2) :: gwz
  real(kind=8), dimension(ngauss_r) :: gxr
  real(kind=8), dimension(ngauss_r,2) :: gwr

end module params
