module phys

  use params
  use func
  use eigen

  contains
  
  subroutine initialize
  
    integer :: i
  
    write(lout,*) 'Deformation parameter beta = ',beta
    
    hbar_omega0 = 41.0/(n_mass**third)
    hbar_omegaz = hbar_omega0*exp(-sq54pi*beta)
    hbar_omegar = hbar_omega0*exp(half*sq54pi*beta)
    
    b0 = hqc/sqrt(nucleon_mass*hbar_omega0) ! in fm
    bz = hqc/sqrt(nucleon_mass*hbar_omegaz) ! in fm
    br = hqc/sqrt(nucleon_mass*hbar_omegar) ! in fm
    
    write(lout,*) 'hbar*omega_0 = ',hbar_omega0,'  b0 = ',b0
    write(lout,*) 'hbar*omega_z = ',hbar_omegaz,'  bz = ',bz
    write(lout,*) 'hbar*omega_r = ',hbar_omegar,'  br = ',br
    
    write(lout,*) 'Highest oscillator shell = ',nmax
    write(lout,'(a17,i3,a2)') 'Highest Omega = ',omegamax,'/2'
    write(lout,*) 'Maximum number of states in a omega block = ',nstatemax
    
    v_c(1) = v_c0*(1.0 - kappa*(n_neutron-n_proton)/(n_neutron+n_proton))
    v_c(2) = v_c0*(1.0 + kappa*(n_neutron-n_proton)/(n_neutron+n_proton))
    v_so = lambda_so*v_c
    
    write(lout,'(a47,2f14.6)') 'Depth of neutron and proton potential (MeV) = ',v_c
    write(lout,'(a58,2f14.6)') 'Depth of neutron and proton spin-orbit potential (MeV) = ',v_so
    
    sfactor = scaling_factor(beta)
    write(lout,*) 'Scaling factor c(beta) = ',sfactor
    
    open(unit=lin,file="gaussh_32.dat",status="old")
    do i=1,ngauss_z
      read(lin,*) gxz(i),gwz(i,1),gwz(i,2)
    end do ! i
    close(lin)
    
    write(lout,*)
    write(lout,*) 'Gauss-Hermite points and weights for z coordinate integration'
    do i=1,size(gxz)
      write(lout,'(3e14.6)') gxz(i),gwz(i,1),gwz(i,2)
    end do
    
    open(unit=lin,file="gaussl_32.dat",status="old")
    do i=1,ngauss_r
      read(lin,*) gxr(i),gwr(i,1),gwr(i,2)
    end do ! i
    close(lin)
    
    write(lout,*)
    write(lout,*) 'Gauss-Laguerre points and weights for z coordinate integration'
    do i=1,size(gxr)
      write(lout,'(3e14.6)') gxr(i),gwr(i,1),gwr(i,2)
    end do
    
    cdens = cdensity(itp)
    write(lout,*)
    write(lout,*) 'Charge density = ',cdens
    call coulpot(coulomb)
    
    delta(1) = 1.0_8
    delta(2) = 1.0_8
    
    call wspot_output
  
  end subroutine initialize
  
  subroutine hamiltonian(base,nbase,ham)
  
    type(state), dimension(:,:), intent(in) :: base
    integer, dimension(:), intent(in) :: nbase
    real(kind=8), dimension(2,nomega,nstatemax,nstatemax), intent(out) :: ham
    
    logical :: isprint
    integer :: i1,i2,j,k,l,isospin
    integer :: nr1,l1,nz1,s1,nr2,l2,nz2,s2,n1,n2
    character(1) :: p1,p2
    real(kind=8) :: cc
    
    isprint = .false.
    
    do isospin=1,2
      do k=1,nomega ! loop over omega blocks
        do i1=1,nbase(k) ! loop over states in a block
          n1 = base(i1,k)%n
          nr1 = base(i1,k)%nrho
          l1 = base(i1,k)%lambda
          nz1 = base(i1,k)%nz
          s1 = base(i1,k)%sigma
          p1 = base(i1,k)%parity
          do i2=1,i1 ! loop over states in a block
            n2 = base(i2,k)%n
            nr2 = base(i2,k)%nrho
            l2 = base(i2,k)%lambda
            nz2 = base(i2,k)%nz
            s2 = base(i2,k)%sigma
            p2 = base(i2,k)%parity
            ham(isospin,k,i1,i2) = 0.0_8
            if(p1.eq.p2) then
              ham(isospin,k,i1,i2) = ham(isospin,k,i1,i2) + kinetic(isospin,nr1,l1,nz1,s1,nr2,l2,nz2,s2)

!              ham(isospin,k,i1,i2) = ham(isospin,k,i1,i2) + nilsson(isospin,nr1,l1,nz1,s1,nr2,l2,nz2,s2)
!              ham(isospin,k,i1,i2) = ham(isospin,k,i1,i2) + nilsson_so(isospin,nr1,l1,nz1,s1,nr2,l2,nz2,s2)

              ham(isospin,k,i1,i2) = ham(isospin,k,i1,i2) + woods_saxon(isospin,nr1,l1,nz1,s1,nr2,l2,nz2,s2)
              ham(isospin,k,i1,i2) = ham(isospin,k,i1,i2) + ws_so(isospin,n1,nr1,l1,nz1,s1,n2,nr2,l2,nz2,s2,i1,i2)
              if(isospin.eq.2) then
                ham(isospin,k,i1,i2) = ham(isospin,k,i1,i2) + coulmat(isospin,nr1,l1,nz1,s1,nr2,l2,nz2,s2)
              end if
            end if ! parity check
            ham(isospin,k,i2,i1) = ham(isospin,k,i1,i2)
          end do ! i2
        end do ! i1
        if(isprint) then
          do i1=1,nbase(k)
            print '(40f12.6)',ham(isospin,k,i1,1:nbase(k))
          end do ! i1
          print *
        end if
      end do ! k
    end do ! isospin
    
  end subroutine hamiltonian
  
! kinetic part of the matrix element
  function kinetic(iso,nr1,l1,nz1,s1,nr2,l2,nz2,s2)
  
    integer, intent(in) :: iso
    integer, intent(in) :: nr1,l1,nz1,s1,nr2,l2,nz2,s2
    
    real(kind=8) :: kinetic
    
    real(kind=8) :: tsum,fact
    
    fact = -(hqc**2)/(2.0*nucleon_mass)
    tsum = 0.0
    if(l1.eq.l2.and.s1.eq.s2) then ! radial terms, abs temporary
      tsum = tsum + matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),ddwfrg(nr2,l2,:),wfzg(nz2,:),unity)
      tsum = tsum + matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),dwfrg(nr2,l2,:),wfzg(nz2,:),ri)
    end if
    if(l1.eq.l2.and.s1.eq.s2) then ! angular term
      tsum = tsum - (l2**2)*matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),wfrg(nr2,l2,:),wfzg(nz2,:),rri)
    end if
    if(l1.eq.l2.and.s1.eq.s2) then ! z term
      tsum = tsum + matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),wfrg(nr2,l2,:),ddwfzg(nz2,:),unity)
    end if
    
    kinetic = fact*tsum
  
  end function kinetic
  
  function nilsson(iso,nr1,l1,nz1,s1,nr2,l2,nz2,s2)
 
    integer, intent(in) :: iso
    integer, intent(in) :: nr1,l1,nz1,s1,nr2,l2,nz2,s2
    
    real(kind=8) :: nilsson
    
    real(kind=8) :: vsum,fz,fr
   
    fz = hbar_omegaz/(2.0*bz*bz)
    fr = hbar_omegar/(2.0*br*br)
    vsum = 0.0
    if(l1.eq.l2.and.s1.eq.s2) then
      vsum = vsum + fz*matel(iso,gxr,gxz,wfrg(nr1,abs(l1),:),wfzg(nz1,:),wfrg(nr2,abs(l2),:),wfzg(nz2,:),zz2)
      vsum = vsum + fr*matel(iso,gxr,gxz,wfrg(nr1,abs(l1),:),wfzg(nz1,:),wfrg(nr2,abs(l2),:),wfzg(nz2,:),rr2)
    end if
    
    nilsson = vsum
    
  end function nilsson
  
! spin-orbit interaction in the Nilsson model
  function nilsson_so(iso,nr1,l1,nz1,s1,nr2,l2,nz2,s2)
 
    integer, intent(in) :: iso
    integer, intent(in) :: nr1,l1,nz1,s1,nr2,l2,nz2,s2
    
    real(kind=8) :: nilsson_so
    
    real(kind=8) :: vsum,f0,t1,t2,t3
   
    f0 = 2.0*hbar_omega0*kappa_so 
    vsum = 0.0
    if(l1-1.eq.l2.and.s1.eq.s2-2) then ! L+ S-
      t1 = matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),wfrg(nr2,l2,:),dwfzg(nz2,:),r1)
      t2 = matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),dwfrg(nr2,l2,:),wfzg(nz2,:),z1)
      t3 = l2*matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),wfrg(nr2,l2,:),wfzg(nz2,:),zr1)
      vsum = vsum + 0.5*(-t1 + t2 - t3)
    end if
    if(l1+1.eq.l2.and.s1.eq.s2+2) then ! L- S+
      t1 = matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),wfrg(nr2,l2,:),dwfzg(nz2,:),r1)
      t2 = matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),dwfrg(nr2,l2,:),wfzg(nz2,:),z1)
      t3 = l2*matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),wfrg(nr2,l2,:),wfzg(nz2,:),zr1)
      vsum = vsum + 0.5*(t1 - t2 - t3)
    end if
    if(l1.eq.l2.and.s1.eq.s2) then ! Lz Sz
      vsum = vsum - l2*(0.5*s2)*matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),wfrg(nr2,l2,:),wfzg(nz2,:),unity)
    end if
    
    nilsson_so = f0*vsum
    
  end function nilsson_so

! matrix elements of the Woods-Saxon potential
  function woods_saxon(iso,nr1,l1,nz1,s1,nr2,l2,nz2,s2)
  
    integer, intent(in) :: iso
    integer, intent(in) :: nr1,l1,nz1,s1,nr2,l2,nz2,s2
    
    real(kind=8) :: woods_saxon
    
    woods_saxon = 0.0
    if(l1.eq.l2.and.s1.eq.s2) then ! abs temporary
      woods_saxon = matel(iso,gxr,gxz,wfrg(nr1,abs(l1),:),wfzg(nz1,:),wfrg(nr2,abs(l2),:),wfzg(nz2,:),wspot)
    end if
  
  end function woods_saxon
    
  function ws_so(iso,n1,nr1,l1,nz1,s1,n2,nr2,l2,nz2,s2,i1,i2)
 
    integer, intent(in) :: iso
    integer, intent(in) :: n1,nr1,l1,nz1,s1,n2,nr2,l2,nz2,s2
    integer, intent(in) :: i1,i2
    
    real(kind=8) :: ws_so
    
    real(kind=8) :: vsum,t1,t2,t3,fact
    
    fact = 2.0*(hqc/(2.0*nucleon_mass))**2 ! factor 2 appears due to spin operators eigenvalues
    vsum = 0.0
    if(l1+1.eq.l2.and.s1.eq.s2+2) then ! sigma+ exp(-i\phi)
      t1 = matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),wfrg(nr2,l2,:),dwfzg(nz2,:),vsodr) 
      t2 = matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),dwfrg(nr2,l2,:),wfzg(nz2,:),vsodz) 
      t3 = matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),wfrg(nr2,l2,:),wfzg(nz2,:),vsodzrho) 
      vsum = vsum + (t1 - t2 - l2*t3)
    end if
    if(l1-1.eq.l2.and.s1.eq.s2-2) then ! sigma- exp(i\phi)
      t1 = matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),wfrg(nr2,l2,:),dwfzg(nz2,:),vsodr) 
      t2 = matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),dwfrg(nr2,l2,:),wfzg(nz2,:),vsodz) 
      t3 = matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),wfrg(nr2,l2,:),wfzg(nz2,:),vsodzrho)
      vsum = vsum - (t1 - t2 + l2*t3)
    end if
    if(l1.eq.l2.and.s1.eq.s2) then ! Lz sigma_z
      vsum = vsum + s2*l2*matel(iso,gxr,gxz,wfrg(nr1,l1,:),wfzg(nz1,:),wfrg(nr2,l2,:),wfzg(nz2,:),vsodrrho)
    end if
    ws_so = -0.5*fact*vsum*v_so(iso)
    
  end function ws_so

! matrix element of the Coulomb interaction
  function coulmat(iso,nr1,l1,nz1,s1,nr2,l2,nz2,s2)
  
    integer, intent(in) :: iso
    integer, intent(in) :: nr1,l1,nz1,s1,nr2,l2,nz2,s2
    
    real(kind=8) :: coulmat
    
    integer :: i
    real(kind=8) :: sumr,sumz,fact
    
    sumr = 0.0
    if(l1.eq.l2.and.s1.eq.s2) then
      do i=1,ngauss_r
        sumz = bz*sum(wfzg(nz1,:)*coulomb(i,:)*wfzg(nz2,:)*gwz(:,2))
        sumr = sumr + wfrg(nr1,l1,i)*sumz*wfrg(nr2,l2,i)*gwr(i,2)
      end do ! i
      sumr = 0.5*br*br*sumr
    end if
    
    coulmat = sumr
  
  end function coulmat

! matrix elements of the Coulomb interaction
! for a spherical homogenous charge distribution
  function coulmatsph(iso,nr1,l1,nz1,s1,nr2,l2,nz2,s2)
  
    integer, intent(in) :: iso
    integer, intent(in) :: nr1,l1,nz1,s1,nr2,l2,nz2,s2
    
    real(kind=8) :: coulmatsph
    
    coulmatsph = 0.0
    if(l1.eq.l2.and.s1.eq.s2) then ! abs temporary
      coulmatsph = matel(iso,gxr,gxz,wfrg(nr1,abs(l1),:),wfzg(nz1,:),wfrg(nr2,abs(l2),:),wfzg(nz2,:),coulpotsph)
    end if
  
  end function coulmatsph

! calculates the scaling factor c(beta) of the radius 
! in order to keep the volume of the nucleus constant
  function scaling_factor(beta)
  
    real(kind=8), intent(in) :: beta
    
    real(kind=8) :: scaling_factor
    
    integer :: i,j
    real(kind=8) :: fsum,fstep,y20,fact,theta
    
    fstep = pi/300.0
    fsum = 0.0
    do i=0,300
      theta = i*fstep
      y20 = 0.25*sqrt(5.0/pi)*(3*(cos(theta)**2)-1.0)
      fact = (1.0 + beta*y20)**3
      if(i.eq.0.or.i.eq.nrzmax) then
        fsum = fsum + fact*sin(theta)
      else
        if(odd(i)) then
          fsum = fsum + 4*fact*sin(theta)
        else
          fsum = fsum + 2*fact*sin(theta)
        end if
      end if
    end do
    fsum = fstep*third*fsum
    
    scaling_factor = (2.0/fsum)**third

  end function scaling_factor


! distance, now obsolete
! should use dist_simple  
  function dist(rho,z,rshape,isospin)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: rho,z
    
    interface
      function rshape(ct,isospin)
        integer, intent(in) :: isospin
        real(kind=8), intent(in) :: ct
        real(kind=8) :: rshape
      end function rshape
    end interface
    
    real(kind=8) :: dist
    
    integer :: i
    real(kind=8) :: r,ctheta,rsigma,stheta
    real(kind=8) :: ct1,ct2,ctmid,r1,r2,rmid,d1,d2,dmid,st1,st2,stmid,dmin
    
    r = sqrt(rho**2 + z**2)
    ctheta = z/r
    stheta = sqrt(1.0 - ctheta**2)
    if(r.lt.1e-5) then
      dist = -rshape(0.0_8,isospin)
    else
! bisection
      i = 0
      ct1 = 1.0
      st1 = 0.0
      r1 = rshape(ct1,isospin)
      d1 = sqrt(r**2 + r1**2 - 2.0*r*r1*(ctheta*ct1+stheta*st1))
      ct2 = -1.0
      st2 = 0.0
      r2 = rshape(ct2,isospin)    
      d2 = sqrt(r**2 + r2**2 - 2.0*r*r2*(ctheta*ct2+stheta*st2))
      ctmid = 0.0
      stmid = 1.0
      rmid = rshape(ctmid,isospin)
      dmid = sqrt(r**2 + rmid**2 - 2.0*r*rmid*(ctheta*ctmid+stheta*stmid))
      do ! infinite loop
        i = i + 1
        dmin = min(d1,d2,dmid)
        if(d1.lt.d2) then
          ct2 = ctmid
          st2 = stmid
          r2 = rmid
          d2 = dmid
        else
          ct1 = ctmid
          st1 = stmid
          r1 = rmid
          d1 = dmid
        end if
        ctmid = 0.5*(ct1 + ct2)
        stmid = sqrt(1.0 - ctmid**2)
        rmid = rshape(ctmid,isospin)
        dmid = sqrt(r**2 + rmid**2 - 2.0*r*rmid*(ctheta*ctmid+stheta*stmid))
        if(abs(dmid-dmin).lt.0.001.or.i.gt.20) then
          exit
        end if
      end do 
      
      if(r.gt.rmid) then
        dist = dmid
      else
        dist = -dmid
      end if
    end if
  
  end function dist
  
! Woods-Saxon potential  
  function wspot(isospin,rho,z)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: rho,z
    
    real(kind=8) :: wspot
    
    wspot = v_c(isospin)/(1.0+exp(dist_simple(rho,z,shape_quadrupole,isospin)/a_c(isospin)))
  
  end function wspot
  
  subroutine coulpot(coulomb)
  
    real(kind=8), dimension(ngauss_r,ngauss_z) :: coulomb ! Coulomb potential
    
    integer :: i,j
    real(kind=8) :: rr,rz
    
    do i=1,ngauss_r
      rr = br*sqrt(gxr(i)) ! transformation from eta to r
      do j=1,ngauss_z/2
        rz = bz*gxz(j) ! transformation from xi to z
        coulomb(i,j) = coulpotpoint(rr,rz)
        coulomb(i,ngauss_z-j+1) = coulomb(i,j)
      end do ! j
    end do ! i
  
  end subroutine coulpot

  
! precalculates Coulomb potential for a deformed homogenous charge distribution
! for a single (rho,z) point
  function coulpotpoint(rho,z)
  
    real(kind=8), intent(in) :: rho,z
    
    real(kind=8) :: coulpotpoint
  
    integer :: i,j
    real(kind=8) :: sumr,sumz,fact,x,eli,den,rprime,zprime
    
    sumr = 0.0
    do i=1,ngauss_r
      rprime = br*sqrt(gxr(i))
      sumz = 0.0
      do j=1,ngauss_z
        zprime = bz*gxz(j)
        den = sqrt((z-zprime)**2 + (rho + rprime)**2)
        x = sqrt(4.0*rho*rprime)/den
        eli = elliptic2(x)
        fact = den*eli*ddws(itp,rprime,zprime)
        sumz = sumz + fact*gwz(j,2)
      end do ! j
      sumz = bz*sumz
      sumr = sumr + sumz*gwr(i,2)
    end do ! i
    sumr = 0.5*br*br*sumr
   
    coulpotpoint = 2.0*sumr*cdens*hqc/alphai
  
  end function coulpotpoint

! Coulomb potential for a homogenous charged sphere
  function coulpotsph(isospin,rho,z)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: rho,z
    
    real(kind=8) :: coulpotsph
    
    real(kind=8) :: r,ctheta,rsurf,pfact,ctmp
    
    pfact = n_proton/alphai
    
    r = sqrt(rho**2 + z**2)
    if(abs(r).lt.0.00001.and.abs(z).lt.0.00001) then 
      ctheta = 1.0
    else
      ctheta = z/r
    endif
    rsurf = shape_quadrupole(ctheta,isospin)
    if(r.gt.rsurf) then
      ctmp = 1.0/r
    else
      ctmp = (3.0*rsurf**2-r**2)/(2.0*rsurf**3)
    end if
    coulpotsph = pfact*ctmp*hqc
  
  end function coulpotsph

! outputs proton and neutron potentials
  subroutine wspot_output
  
    integer :: iso,i,j,n
    character(14), dimension(2) :: fname = (/"wspot_prot.dat","wspot_neut.dat"/)
    real(kind=8) :: x,y,r
    
    n = 200
    r = 8.0
    do iso=1,2
      open(unit=ltmp,file=fname(iso),status="replace")
      do i=0,n
        do j=0,n
          x = i*r/n
          y = -r + j*2.0*r/n
          write(ltmp,'(3f12.6)') x,y,wspot(iso,x,y)
        end do ! j
      end do ! i
      close(ltmp)
    end do ! iso
  
  end subroutine wspot_output
  
! driver routine for diagonalization of each omega 
! block hamiltonian  
  subroutine solve(nbase,ham,energy,evectors)
  
    integer, dimension(:), intent(in) :: nbase
    real(kind=8), dimension(2,nomega,nstatemax,nstatemax), intent(in) :: ham
    real(kind=8), dimension(2,nomega,nstatemax),  intent(out) :: energy
    real(kind=8), dimension(2,nomega,nstatemax,nstatemax), intent(out) :: evectors 
    
    logical :: isene, isvec
    integer :: iso,i,n,j
    real(kind=8), dimension(:), allocatable :: hene
    real(kind=8), dimension(:,:), allocatable :: hdiag,hvec
    
    isene = .false.
    isvec = .false.
    
    do iso=1,2
      do i=1,nomega
        n = nbase(i)
        allocate(hdiag(n,n))
        allocate(hene(n))
        allocate(hvec(n,n))
        
        hdiag = ham(iso,i,1:n,1:n)
        call diagonalize(hdiag,hene,hvec)
        energy(iso,i,1:n) = hene
        evectors(iso,i,1:n,1:n) = hvec
        
        if(isene) then
          print '(i3,a2)',(2*i-1),'/2'
          print '(10f12.6)',hene/hbar_omega0
          print *
        end if
        
        if(isvec) then
          do j=1,n
            print '(40f12.6)',hvec(j,:)
          end do ! j
          print *
        end if
        
        deallocate(hdiag)
        deallocate(hene)
        deallocate(hvec)
      end do ! i
    end do ! iso
  
!    allocate(hene(sum(nbase(1:nomega))))
    
!    j = 0
!    do i=1,nomega
!      hene(j+1:j+1+nbase(i)) = energy(1,i,1:nbase(i))
!      j = j + nbase(i)
!    end do ! i
!    write(lplot,'(30f14.6)') hene
    
!    deallocate(hene)
  
  end subroutine solve
  
  subroutine fill_orbits(nbase,energy,evectors,etot)
  
    integer, dimension(:), intent(in) :: nbase
    real(kind=8), dimension(2,nomega,nstatemax),  intent(in) :: energy
    real(kind=8), dimension(2,nomega,nstatemax,nstatemax), intent(in) :: evectors 
    type(estate), dimension(2,nomega*nstatemax), intent(out) :: etot 
    
    integer :: iso,i,j,n,cnt,k,ir
    integer, dimension(2) :: pcnt
    real(kind=8) :: be
    type(estate) :: rra
    
! assign energy and quantum numbers to eigenstates    
    n = sum(nbase)
    do iso=1,2
      cnt = 1
      do i=1,nomega
        do j=1,nbase(i)
          etot(iso,cnt)%omega = 2*i-1
          etot(iso,cnt)%energy = energy(iso,i,j)
          do k=1,nbase(i)
            if(abs(evectors(iso,i,j,k)).gt.0.2) then
              etot(iso,cnt)%parity = basis(k,i)%parity
              exit
            end if
          end do ! k
          etot(iso,cnt)%bdim = nbase(i)
          etot(iso,cnt)%vector = evectors(iso,i,1:nbase(i),j)
          cnt = cnt + 1
        end do ! j
      end do ! i
    end do ! iso
    
! sort eigenstates by energy using heap sort
    do iso=1,2
      l = n/2 + 1
      ir = n
      do 
        if(l.gt.1)then
          l = l - 1
          rra = etot(iso,l)
        else          
          rra = etot(iso,ir)
          etot(iso,ir) = etot(iso,1)
          ir = ir - 1     
          if(ir.eq.1)then
            etot(iso,1) = rra
            exit
          endif
        endif
        i = l
        j = 2*l
        do while(j.le.ir)
          if(j.lt.ir.and.etot(iso,j)%energy.gt.etot(iso,j+1)%energy) j = j + 1
          if(rra%energy.gt.etot(iso,j)%energy) then
            etot(iso,i) = etot(iso,j)
            i = j
            j = 2*j
          else
            j = ir + 1
          end if
          etot(iso,i) = rra
        end do ! do while loop
      end do ! infinite loop
    end do ! iso
    
! occupation probabilities    
    pcnt(1) = n_neutron
    pcnt(2) = n_proton
    do iso=1,2
      do i=n,1,-1
        if(pcnt(iso).eq.0) then
          etot(iso,i)%vv = 0.0
          etot(iso,i)%uu = 1.0
        elseif(pcnt(iso).eq.1) then
          etot(iso,i)%vv = 0.5
          etot(iso,i)%uu = 0.5
          pcnt(iso) = 0
        else
          etot(iso,i)%vv = 1.0
          etot(iso,i)%uu = 0.0
          pcnt(iso) = pcnt(iso) - 2
        end if
      end do
    end do ! iso
    
  end subroutine fill_orbits
  
! finds chemical potentials via particle number
  subroutine fchem(lambda,nbase,etot,gamma,npart)
  
    integer, dimension(:), intent(in) :: nbase
    type(estate), dimension(2,nomega*nstatemax), intent(out) :: etot 
    real(kind=8), intent(in) :: gamma
    real(kind=8), dimension(2), intent(out) :: lambda
    real(kind=8), dimension(2), intent(out) :: npart
    
    integer :: n,iso,i
    real(kind=8) :: llow,lhigh,lmid,l4,lold
    real(kind=8) :: flow,fhigh,fmid,f4,sqn
  
    n = sum(nbase)
    npart(1) = n_neutron
    npart(2) = n_proton
    do iso=1,2 ! Ridders' method for finding chemical potentials
      do i=n,1,-1 ! bracketing first
        if(etot(iso,i)%vv.lt.0.4) then
          llow = etot(iso,i)%energy - 20.0
          lhigh = etot(iso,i)%energy + 20.0
          exit
        end if
      end do ! i
      lold = 0.5*(llow - lhigh)
      flow = npart(iso) - pnum(iso,llow,nbase,etot,gamma)
      fhigh = npart(iso) - pnum(iso,lhigh,nbase,etot,gamma)
      do 
        lmid = 0.5*(llow + lhigh)
        fmid = npart(iso) - pnum(iso,lmid,nbase,etot,gamma)
        sqn = sqrt(fmid**2 - flow*fhigh)
        l4 = lmid + (lmid - llow)*sign(1.0_8,flow - fhigh)*fmid/sqn
        if(abs(l4-lold).lt.0.0001) then
          lambda(iso) = l4
          exit
        end if
        f4 = npart(iso) - pnum(iso,l4,nbase,etot,gamma)
        if(abs(f4).lt.1.0e-8) then
          lambda(iso) = l4
          exit
        end if
        if(sign(fmid,f4).ne.fmid) then
          llow = lmid
          flow = fmid
          lhigh = l4
          fhigh = f4
        elseif(sign(flow,f4).ne.fl) then
          lhigh = l4
          fhigh = f4
        else
          llow = l4
          flow = f4
        end if
        lold = l4
      end do ! infinite loop
      f4 = pnum(iso,lambda(iso),nbase,etot,gamma)
!      print *,'f4 = ',iso,l4,f4
    end do ! isospin
  
  end subroutine fchem
  
! calculates particle number for a given
! chemical potential
  function pnum(isospin,lam,nbase,etot,gamma)
  
    integer, intent(in) :: isospin
    integer, dimension(:), intent(in) :: nbase
    type(estate), dimension(2,nomega*nstatemax), intent(out) :: etot 
    real(kind=8), intent(in) :: gamma
    real(kind=8), intent(in) :: lam
        
    real(kind=8) :: pnum
    
    integer :: n
    real(kind=8) :: ti,stmp
    
    n = sum(nbase)
    stmp = 0.0
    do i=1,n
      ti = (lam - etot(isospin,i)%energy)/gamma
      etot(isospin,i)%vvtilde = gonint(ti)
      stmp = stmp + 2.0*etot(isospin,i)%vvtilde
    end do ! i
    pnum = stmp
    
  end function pnum
  
! generalized occupation numbers integral  
  function gonint(ti)
  
    real(kind=8), intent(in) :: ti
    
    real(kind=8) :: gonint
    
    integer :: i
    real(kind=8) :: stmp,step,x
    
    stmp = 0.0
    if(ti.gt.5.0) then
      stmp = 1.0
    elseif(ti.lt.4.and.ti.gt.-5.0) then
      step = (ti+5.0)/1500
      do i=0,1500
        x = -5.0 + i*step
        if(i.eq.0.or.i.eq.1500) then
          stmp = stmp + sdf(x)
        elseif(odd(i)) then
          stmp = stmp + 4.0*sdf(x)
        else
          stmp = stmp + 2.0*sdf(x)
        end if
      end do ! i
      stmp = stmp*step*third
    end if
    gonint = stmp
  
  end function gonint
  
end module phys
