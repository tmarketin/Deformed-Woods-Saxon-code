module wavefunction

  use params
  use func
  use phys

  contains

! driver routine for initialization
  subroutine initw
  
    logical :: test
  
! generate coordinates  
    test = coordinates(rz,rzbz,rr,rrbr)
    
! generate z-coordinate wavefunctions  
    test = phi_z(gxz,wfzg,dwfzg,ddwfzg)
  
! generate radial wavefunctions
    test = phi_radial(gxr,wfrg,dwfrg,ddwfrg)
    
! generate basis of the axially deformed harmonic oscillator
    call create_basis(basis,nbasis)
  
! test norm and orthogonality of wavefunctions    
!    call wf_test(basis,nbasis,gxr,gxz,wfzg,wfrg)
  
  end subroutine initw
  
! makes an array holding values of z coordinate
! used to calculate hermite polynomials  
  function coordinates(rz,rzbz,rr,rrbr)
  
    real(kind=8), dimension(0:nrzmax), intent(out) :: rz
    real(kind=8), dimension(0:nrzmax), intent(out) :: rzbz
    real(kind=8), dimension(0:nrrmax), intent(out) :: rr
    real(kind=8), dimension(0:nrrmax), intent(out) :: rrbr
    
    logical :: coordinates
    
    integer :: i
    
    do i=0,nrzmax
      rz(i) = rzmin + i*rzstep
      rzbz(i) = rz(i)/bz
    end do
    
    rr(0) = rrstep/10.0
    rrbr(0) = (rr(0)/br)**2
    do i=1,nrrmax
      rr(i) = i*rrstep
      rrbr(i) = (rr(i)/br)**2
    end do ! i
    
    open(unit=12,file="coords_z.dat",status="replace")
    do i=0,nrzmax
      write(12,'(2f12.6)'),rz(i),rzbz(i)
    enddo
    close(12)
    
    open(unit=12,file="coords_rho.dat",status="replace")
    do i=0,nrrmax
      write(12,'(2f12.6)'),rr(i),rrbr(i)
    enddo
    close(12)
    
    print *,'Coordinates generated.'
    write(lout,*) 'Coordinates generated.'
   
    coordinates = .true.
  
  end function coordinates
  
! calculates hermite polynomials and puts them in an 
! array, first index is the order of the polynomials, 
! second index denotes the coordinate
! recursion relation is used
! H_{n+1}(x) = 2*x*H_{n}(x) - 2*n*H_{n-1}(x)
  function hermite(n,x)
  
    integer, intent(in) :: n
    real(kind=8), intent(in) :: x
    
    real(kind=8) :: hermite
    
    integer :: i,j
    real(kind=8) :: h0,h1,h2
    
    h1 = 1.0
    h2 = 2.0*x
    if(n.eq.0) then
      hermite = h1
    elseif(n.eq.1) then
      hermite = h2
    else
      do i=2,n
        h0 = h1
        h1 = h2
        h2 = 2*x*h1-2*(i-1)*h0
      end do ! i
      hermite = h2
    endif
    
  end function hermite
  
! normalization factor of the z-coordinate wave function
  function hnorm(n)
  
  	integer, intent(in) :: n
  	
  	real(kind=8) :: hnorm
  	
  	integer :: i
  	real(kind=8) :: fact,nrm
  	
  	fact = 1.0/sqrt(sqrt(pi))
  	nrm = 1.0
  	if(n.eq.0) then
  		hnorm = fact*nrm
  	else
			do i=1,n
				nrm = nrm/sqrt(2.0*i)
			end do ! i
			hnorm = fact*nrm
  	end if
  
  end function hnorm
  
! calculates wavefunctions in the z coordinate  
! norm is 1/\sqrt{bz * 2**nz * n! * sqrt{pi}}
! wf = norm * H_nz * exp(-0.5*rzbz**2)
! derivatives are with respect to z, not z/bz
  function phi_z(zeta,wfg,dwfg,ddwfg)
  
    real(kind=8), dimension(ngauss_z), intent(in) :: zeta
    real(kind=8), dimension(0:nzmax,ngauss_z), intent(out) :: wfg,dwfg,ddwfg
      
    logical :: phi_z
    
    logical :: int_test
    integer :: i,j
    real(kind=8) :: norm,sbz,sspi,ezhalf,wfsum
    
    phi_z = .true.
    int_test = .false.
    
    sbz = 1.0/sqrt(bz) 
    sspi = 1.0/sqrt(sqrt(pi))
    
! wavefunctions
    norm = sspi*sbz
    do j=0,nzmax
      if(j.gt.0) then
        norm = norm/(sqrt(2.0*j))
      end if
      do i=1,ngauss_z
        wfg(j,i) = norm*hermite(j,zeta(i))*exp(-(zeta(i)**2)/2.0)
      end do ! i
    end do ! j
    
! normalization test    
    if(int_test) then
      do i=0,nzmax
        print '(i3,f18.12)',i,sum(wfg(i,:)*wfg(i,:)*gwz(:,2))*bz
      end do
      print *
    end if
    
! first derivative
! phi'(n) = \sqrt{2*n}*phi(n-1)/bz - rzbz*phi(n)/bz
!    dwfg(0,:) = -zeta*wfg(0,:)/bz
!    do i=1,nzmax
!      dwfg(i,:) = -zeta*wfg(i,:)/bz + sqrt(2.0*i)*wfg(i-1,:)/bz
!    end do ! i
    
! first derivative, alternate
    norm = sspi/(bz**1.5)
    do j=0,nzmax
      if(j.gt.0) then
        norm = norm/(sqrt(2.0*j))
      end if
      do i=1,ngauss_z
        dwfg(j,i) = norm*(zeta(i)*hermite(j,zeta(i))-hermite(j+1,zeta(i)))*exp(-0.5*(zeta(i)**2))
      end do ! i
    end do ! j
    
! first derivative integration test
    if(int_test) then
      do i=0,nzmax
        print '(i3,f18.12)',i,sum(dwfg(i,:)*dwfg(i,:)*gwz(:,2))*bz
      end do
      print *
    end if    
    
! second derivative
! phi''(n) = \sqrt{2*n}*phi'(n-1) - phi(n) - rzbz*phi'(n)
    ddwfg(0,:) = -wfg(0,:)/(bz**2)-zeta*dwfg(0,:)/bz
    do i=1,nzmax
      ddwfg(i,:) = -wfg(i,:)/(bz**2) - zeta*dwfg(i,:)/bz + sqrt(2.0*i)*dwfg(i-1,:)/bz
    end do ! i

! second derivative integration test
    if(int_test) then
      do i=0,nzmax
        print '(i3,f18.12)',i,sum(ddwfg(i,:)*ddwfg(i,:)*gwz(:,2))*bz
      end do
      print *
    end if    
    
    print *,'Z coordinate wave functions and derivatives generated.'
    write(lout,*) 'Z coordinate wave functions and derivatives generated.'
  
  end function phi_z
  
! generalized Laguerre polynomials
! f(n,l,k) is L^{l-1}_{n-1}(r_k)
! first L(0,0) = 1.0, L(1,0) = 1 - r
! recursion relation
! L_{k+1} = [(2k+1-x) L_{k} - k L_{k-1}]/(k+1)
! and
! L^{l+1}_{n} = \sum_{i=0}^{n} L^{l}_{i}
  function laguerre(r,f)
  
    real(kind=8), dimension(0:nrrmax), intent(in) :: r
    real(kind=8), dimension(0:nrmax,0:nlmax+1,0:nrrmax), intent(out) :: f
    
    logical :: laguerre
    
    integer :: n,l,k
    real(kind=8), dimension(0:npol) :: hprint
    
    f(0,0,:) = 1.0_8
    f(1,0,:) = 1.0_8 - r
    do n=2,nrmax
      f(n,0,:) = ((2*(n-1)+1-r)*f(n-1,0,:)-(n-1)*f(n-2,0,:))/n
    end do ! n
    do l=1,nlmax+1
      do n=0,nrmax
        f(n,l,:) = 0.0_8
        do k=0,n
          f(n,l,:) = f(n,l,:) + f(k,l-1,:)
        end do ! k
      end do ! n
    end do ! l
    
    open(unit=12,file="laguerre.dat",status="replace")
    do j=0,nrrmax
      hprint(0) = r(j)
      do n=0,nrmax
        do l=0,nlmax
          hprint(n*(nlmax+1)+l+1) = f(n,l,j)
        end do ! l
      end do ! n
      write (12,'(f10.6,20f18.10)') hprint
    end do
    close(12)
    
    print *,'Laguerre polynomials generated.'
    write(lout,*) 'Laguerre polynomials generated.'
    laguerre = .true.
    
  end function laguerre
  
! another method of calculating Laguerre polynomials  
  function lag(n,alpha,x)
  
    integer, intent(in) :: n
    integer, intent(in) :: alpha
    real(kind=8), intent(in) :: x
    
    real(kind=8) :: lag
    
    integer :: i
    real(kind=8) :: l0,l1,lfin
    
    if(n.eq.0) then
      lfin = 1.0
    else
      l1 = 0.0
      lfin = 1.0
      do i=1,n
        l0 = l1
        l1 = lfin
        lfin = ((2.0*i - 1 + alpha - x)*l1 - (i - 1 + alpha)*l0)/ i
      end do
    end if
    
    lag = lfin
  
  end function lag
  
! norm of the radial wavefunction  
  function lagnorm(n,l)
  
    integer, intent(in) :: n
    integer, intent(in) :: l
    
    real(kind=8) :: lagnorm
    
    integer :: i
    real(kind=8) :: tmp
    
    tmp = 1.0
    if(l.eq.0) then
      lagnorm = tmp
    else
      do i=n+1,n+l
        tmp = tmp/sqrt(real(i))
      end do ! i
        lagnorm = tmp
    end if
  
  end function lagnorm
  
  function phi_radial(eta,wfg,dwfg,ddwfg)
  
    real(kind=8), dimension(ngauss_r), intent(in) :: eta
    real(kind=8), dimension(0:nrmax,0:nlmax,ngauss_r), intent(out) :: wfg,dwfg,ddwfg
    
    logical :: phi_radial
    
    logical :: int_test
    integer :: i,j,k
    real(kind=8) :: eeta,poly,fact,fd1,fd2,fdr,norm
    
    phi_radial = .true.
    int_test = .false.
    
! wavefunctions    
    do k=1,ngauss_r ! loop over gauss points
      eeta = exp(-0.5*eta(k))
      poly = 1.0
      do j=0,nlmax ! loop over angular momentum projection
        if(j.gt.0) then
          poly = poly*sqrt(eta(k))
        endif
        do i=0,nrmax ! loop over n_r
          fact = lagnorm(i,j)*sqrt(2.0)/br
          wfg(i,j,k) = fact*poly*eeta*lag(i,j,eta(k))
        end do ! i
      end do ! j
    end do ! k
    
! normalization test
    if(int_test) then
      print *,'Normalization test :'
      do i=0,nrmax
        do j=0,nlmax
          norm = sum(wfg(i,j,:)*gwr(:,2)*wfg(i,j,:))*(br**2)/2
          print '(2i3,f18.12)',i,j,norm
          if(abs(norm-1.0).gt.1e-5) then
            phi_radial = .false.
          endif
        end do ! j
      end do ! i    
      print *
    end if
    
! first derivatives
    do i=0,nrmax
      do j=0,nlmax
        do k=1,ngauss_r
          fd1 = 0.5*(j+2.0*i-eta(k))
          fdr = 2.0/(br*sqrt(eta(k)))
          if(i.eq.0) then
            dwfg(i,j,k) = fd1*wfg(i,j,k)*fdr
          else
            fd2 = sqrt(1.0*i*(i+j))
            dwfg(i,j,k) = (fd1*wfg(i,j,k) - fd2*wfg(i-1,j,k))*fdr
          endif
        end do ! k
      end do ! j
    end do ! i
    
! integration test
    if(int_test) then
      print *,'First derivative integration test :'
      do i=0,nrmax
        do j=0,nlmax
          norm = sum(dwfg(i,j,:)*gwr(:,2)*dwfg(i,j,:))*(br**2)/2
          print '(2i3,f18.12)',i,j,norm
        end do ! j
      end do ! i    
      print *
    end if

! second derivatives
    do i=0,nrmax
      do j=0,nlmax
        do k=1,ngauss_r
          fd1 = (j+2.0*i-1.0-eta(k))/(br*sqrt(eta(k)))
          fdr = 2.0/(br**2)
          if(i.eq.0) then
            ddwfg(i,j,k) = fd1*dwfg(i,j,k)-fdr*wfg(i,j,k)
          else
            fd2 = 2.0*sqrt(1.0*i*(i+j))/(br*sqrt(eta(k)))
            ddwfg(i,j,k) = fd1*dwfg(i,j,k) - fd2*dwfg(i-1,j,k) - fdr*wfg(i,j,k)
          endif
        end do ! k
      end do ! j
    end do ! i    
    
! integration test
    if(int_test) then
      print *,'Second derivative integration test :'
      do i=0,nrmax
        do j=0,nlmax
          norm = sum(ddwfg(i,j,:)*gwr(:,2)*ddwfg(i,j,:))*(br**2)/2
          print '(2i3,f18.12)',i,j,norm
        end do ! j
      end do ! i    
      print *
    end if
  
  end function phi_radial

! generates radial wavefunctions
! for l=0 wf = norm*L^{0}_{nr}(rrbr)*exp(-0.5*rrbr)  
! for l>0 use recursion relation
! x L^{l+1}_{n}(x) = (n+l)L^{l}_{n-1}(x) - (n-x)L^{l}_{n}
! derivatives with respect to r, not eta
  function phi_r(r,eta,l_poly,wf,dwf,ddwf)
  
    real(kind=8), dimension(0:nrrmax) :: r
    real(kind=8), dimension(0:nrrmax) :: eta
    real(kind=8), dimension(0:nrmax,0:nlmax+1,0:nrrmax), intent(in) :: l_poly  
    real(kind=8), dimension(0:nrmax,0:nlmax,0:nrrmax), intent(out) :: wf
    real(kind=8), dimension(0:nrmax,0:nlmax,0:nrrmax), intent(out) :: dwf
    real(kind=8), dimension(0:nrmax,0:nlmax,0:nrrmax), intent(out) :: ddwf
    
    logical phi_r
    
    integer :: n,l,i,j
    real(kind=8) :: norm,fact,wfsum
    real(kind=8), dimension(0:nrrmax) :: ehalf,sqeta,sqetai,term1,term2,term3
    real(kind=8), dimension(0:npol) :: hprint
    real(kind=8), dimension(0:nrmax,0:nlmax,0:nrrmax):: dwf_test
    
    phi_r = .true.
    
    do i=0,nrrmax
      ehalf(i) = exp(-0.5*eta(i))
      sqeta(i) = sqrt(eta(i))
      sqetai(i) = 1.0
    end do ! i
    
    norm = sqrt(2.0)/br 
    do n=0,nrmax
      wf(n,0,:) = norm*l_poly(n,0,:)*ehalf
    end do ! n
    
    do l=1,nlmax
      sqetai = sqetai*sqeta
      do n=0,nrmax
        fact = norm*lagnorm(n,l)
        wf(n,l,:) = fact*sqetai*l_poly(n,l,:)*ehalf
      end do ! n
    end do ! l
    
! first derivative, completely numerically
    do n=0,nrmax
      do l=0,nlmax
        dwf(n,l,0) = (wf(n,l,1)-wf(n,l,0))/rrstep
        do i=1,nrrmax-1
          dwf(n,l,i) = (wf(n,l,i+1)-wf(n,l,i-1))/(2.0*rrstep)
        end do ! i
        dwf(n,l,0) = 0.0
      end do ! l
    end do ! n

! first derivative
!    dwf(0,0,:) = sqeta*(-1)*wf(0,0,:)/br
!    sqetai = 1.0
!    do l=1,nlmax
!      dwf(0,l,:) = sqrt(2.0)*l*sqetai*l_poly(n,l,:)*ehalf*lagnorm(0,l)/(br**2) - sqeta*wf(0,l,:)/br
!      sqetai = sqetai * sqeta
!    end do ! l
!    do n=1,nrmax
!      dwf(n,0,:) = -2*sqeta*(0.5*wf(n,0,:) + lagnorm(n,0)*norm*l_poly(n-1,1,:)*ehalf)/br
!      sqetai = 1.0
!      do l=1,nlmax
!        term1 = lagnorm(n,l)*l*sqrt(2.0)*sqetai*l_poly(n,l,:)*ehalf/(br**2)
!        term2 = sqeta*wf(n,l,:)/br
!        term3 = lagnorm(n,l)*2.0*sqrt(2.0)*sqetai*sqeta*sqeta*l_poly(n-1,l+1,:)*ehalf/(br**2)
!        dwf(n,l,:) = term1 - term2 - term3
!        sqetai = sqetai * sqeta
!      end do ! l
!    end do ! n
    
! second derivative, completely numerically
    do n=0,nrmax
      do l=0,nlmax
        ddwf(n,l,0) = (dwf(n,l,1)-dwf(n,l,0))/rrstep
        do i=1,nrrmax-1
          ddwf(n,l,i) = (dwf(n,l,i+1)-dwf(n,l,i-1))/(2.0*rrstep)
        end do ! i
        ddwf(n,l,0) = 0.0
      end do ! l
    end do ! n
    
    print *,'Radial wave functions and derivatives generated.'
    write(lout,*) 'Radial wave functions and derivatives generated.'

    open(unit=12,file="wf_r.dat",status="replace")
    do j=0,nrrmax
      hprint(0) = eta(j)
      do n=0,nrmax
        do l=0,nlmax
          hprint(n*(nlmax+1)+l+1) = wf(n,l,j)
        end do ! l
      end do ! n
      write (12,'(f10.6,100e18.10)') hprint
    end do
    close(12)

    open(unit=12,file="dwf_r.dat",status="replace")
    do j=0,nrrmax
      hprint(0) = eta(j)
      do n=0,nrmax
        do l=0,nlmax
          hprint(n*(nlmax+1)+l+1) = dwf(n,l,j)
        end do ! l
      end do ! n
      write (12,'(f10.6,100f18.10)') hprint
    end do
    close(12)

    open(unit=12,file="ddwf_r.dat",status="replace")
    do j=0,nrrmax
      hprint(0) = eta(j)
      do n=0,nrmax
        do l=0,nlmax
          hprint(n*(nlmax+1)+l+1) = ddwf(n,l,j)
        end do ! l
      end do ! n
      write (12,'(f10.6,100f18.10)') hprint
    end do
    close(12)

! normalization test    
!    do n=0,nrmax
!      do l=0,nlmax
!        wfsum = r(0)*(wf(n,l,0)**2) + r(nrrmax)*(wf(n,l,nrrmax)**2)
!        do i=2,nrrmax-1,2
!          wfsum = wfsum + r(i)*2*(wf(n,l,i)**2)
!        end do ! i
!        do i=1,nrrmax-1,2
!          wfsum = wfsum + r(i)*4*(wf(n,l,i)**2)
!        end do ! i
!      wfsum = rrstep*third*wfsum
!      if(abs(wfsum-1.0).gt.epsint) phi_r = .false.
!      write(lout,*) 'Integral of n ',n,' l ',l,' = ',wfsum
!      end do ! l
!    end do ! n
    
!    if(phi_r) then
!      write(lout,*) 'Radial wavefunctions normalized.'
!    else
!      write(lout,*) 'Radial wavefunctions not normalized.'
!    end if
  
  end function phi_r
  
  subroutine create_basis(base,nbase)
  
    type(state), dimension(:,:), intent(out) :: base
    integer, dimension(nomega), intent(out) :: nbase
    
    integer :: nz,nrho,nperp,lambda,sigma,l,parity,n,ome,omega,counter
    character(16) :: sstate
    character(10) :: somega,sspin
    
    open(ltmp,file="basis.dat",status="replace")
    
    write(lout,*) 'Number of states per omega value: '
    do omega=1,omegamax,2
      counter = 0
      write(ltmp,*) 'Omega = ',omega,'/2'
      write(ltmp,'(9a10)') 'N |','nz |','Lambda |','Omega |','Sigma |','l |','Parity |','nperp |','nrho |'
      write(ltmp,'(a90)') repeat('-',90)
      do n=0,nmax
        do nz = n,0,-1
          nperp = n - nz
          do lambda = -nperp,nperp,2
            nrho = (nperp - abs(lambda))/2
            l = nz + abs(lambda)
            if(odd(l)) then ! parity pi = (-1)**(nz + |Lambda|)
              parity = -1
            else
              parity = 1
            end if 
            do sigma=-1,1,2
              ome = 2*lambda+sigma
              if(ome.eq.omega) then ! abs(ome)
                write(sspin,'(i8,a1,i1)') sigma,'/',2
                write(somega,'(i8,a1,i1)') ome,'/',2
                if(nmax.lt.10) then
                  write(sstate,'(a1,3i1,a1,i3,a1,i1)') '[',n,nz,abs(lambda),']',abs(ome),'/',2
                  write(ltmp,'(i7,2i10,4a10,2i10,a19)') n,nz,lambda,somega,sspin,orbit(l),cparity(parity),nperp,nrho,sstate
                else
                  write(sstate,'(a1,3i3,a1,i3,a1,i1)') '[',n,nz,abs(lambda),']',abs(ome),'/',2
                  write(ltmp,'(i7,2i10,4a10,2i10,a19)') n,nz,lambda,somega,sspin,orbit(l),cparity(parity),nperp,nrho,sstate
                end if
                counter = counter + 1
                base(counter,(omega+1)/2)%n = n
                base(counter,(omega+1)/2)%nz = nz
                base(counter,(omega+1)/2)%nperp = nperp
                base(counter,(omega+1)/2)%nrho = nrho
                base(counter,(omega+1)/2)%lambda = lambda
                base(counter,(omega+1)/2)%l = l
                base(counter,(omega+1)/2)%sigma = sigma ! 1 or -1
                base(counter,(omega+1)/2)%omega = ome ! only numerator
                base(counter,(omega+1)/2)%parity = cparity(parity) ! character
              end if
            end do ! sigma
          end do ! lambda
        end do ! nz
      enddo ! n
      nbase((omega+1)/2) = counter
      write(lout,'(a9,i3,a3,a20,i4)') 'Omega = ',omega,'/2,','  number of states = ',nbase((omega+1)/2)
      write(ltmp,*) 
    enddo ! omega
    print *,'Axially deformed harmonic oscillator basis generated.'
    write(lout,*) 'Axially deformed harmonic oscillator basis generated.'
    
    close(ltmp)
  
  end subroutine create_basis
  
  subroutine wf_test(base,nbase,rr,rz,wfzg,wfrg)
  
    type(state), dimension(:,:), intent(out) :: base
    integer, dimension(nomega), intent(out) :: nbase
    real(kind=8), dimension(ngauss_r), intent(in) :: rr
    real(kind=8), dimension(ngauss_z), intent(in) :: rz
    real(kind=8), dimension(0:nzmax,ngauss_z), intent(in) :: wfzg
    real(kind=8), dimension(0:nrmax,0:nlmax,ngauss_r), intent(in) :: wfrg

    logical :: isnorm,isort
    integer i1,j1,nr1,l1,nz1,s1
    integer i2,j2,nr2,l2,nz2,s2
    real(kind=8) :: sproduct

    isnorm = .true.
    isort = .true.
    do i1=1,omegamax,2
      do j1=1,nbase((i1+1)/2)
        s1 = base(j1,(i1+1)/2)%sigma
        nz1 = base(j1,(i1+1)/2)%nz
        nr1 = base(j1,(i1+1)/2)%nrho
        l1 = base(j1,(i1+1)/2)%lambda
        do i2=1,omegamax,2
          do j2=1,nbase((i2+1)/2)
            s2 = base(j2,(i2+1)/2)%sigma
            nz2 = base(j2,(i2+1)/2)%nz
            nr2 = base(j2,(i2+1)/2)%nrho
            l2 = base(j2,(i2+1)/2)%lambda
            if(l1.ne.l2.or.s1.ne.s2) then ! angular and spin part of the product
              sproduct = 0.0 
            else
              sproduct = matel(1,rr,rz,wfrg(nr1,abs(l1),:),wfzg(nz1,:),wfrg(nr2,abs(l2),:),wfzg(nz2,:),unity)
            end if
            if(i1.eq.i2.and.j1.eq.j2) then
              if(abs(sproduct-1.0).gt.epsint) then
                isnorm = .false.
                write(lout,*) 'State (nr,lambda,nz) ',nr1,l1,nz1,' not normalized, norm = .',sproduct
              end if
            else
              if(abs(sproduct-0.0).gt.epsint) then
                isort = .false.
                print *,i1,j1,i2,j2,sproduct,l1,l2
              end if
            end if
          end do ! j2
        end do ! i2
      end do ! j1
    end do ! i1
    
    if(isnorm) then
      print *,'All states normalized.'
      write(lout,*) 'All states normalized.'
    end if
    if(isort) then
      print *,'All states orthogonal.'
      write(lout,*) 'All states orthogonal.'
    end if

  end subroutine wf_test
  
end module wavefunction
