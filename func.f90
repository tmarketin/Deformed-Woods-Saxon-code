module func

  use params

  contains

  function odd(i)
  
    integer, intent(in) :: i
    
    logical :: odd
    
    odd = .true.
    if(mod(i,2).eq.0) then
      odd = .false.
    end if
  
  end function odd
  
  function orbit(l)
  
    integer, intent(in) :: l
    
    character(1) :: orbit
    
    select case (l)
    case(0)
      orbit = 's'
    case(1)
      orbit = 'p'
    case(2)
      orbit = 'd'
    case(3)
      orbit = 'f'
    case(4)
      orbit = 'g'
    case(5)
      orbit = 'h'
    case(6)
      orbit = 'i'
    case(7)
      orbit = 'j'
    case(8)
      orbit = 'k'
    case(9)
      orbit = 'l'
    case(10)
      orbit = 'm'
    case(11)
      orbit = 'n'
    end select
  
  end function orbit
  
  function cparity(p)
  
    integer, intent(in) :: p
    
    character(1) :: cparity
    
    if(p.eq.-1) then
      cparity = '-'
    else
      cparity = '+'
    end if
  
  end function cparity

  function unity(isospin,r,z)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: r,z
    real(kind=8) :: unity
    
    unity = 1.0
  
  end function unity
  
! returns 1/rho  
  function ri(isospin,r,z)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: r,z
    real(kind=8) :: ri
    
    ri = 1.0/r
  
  end function ri

! returns (1/rho)^{2}
  function rri(isospin,r,z)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: r,z
    real(kind=8) :: rri
    
    rri = 1.0/(r**2)
  
  end function rri
  
! returns (rho)^{2}
  function rr2(isospin,r,z)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: r,z
    real(kind=8) :: rr2
    
    rr2 = r**2
  
  end function rr2
  
! returns (z)^{2}
  function zz2(isospin,r,z)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: r,z
    real(kind=8) :: zz2
    
    zz2 = z**2
  
  end function zz2
  
! returns z
  function z1(isospin,r,z)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: r,z
    real(kind=8) :: z1
    
    z1 = z
  
  end function z1
  
! returns rho
  function r1(isospin,r,z)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: r,z
    real(kind=8) :: r1
    
    r1 = r
  
  end function r1
  
! returns z/rho
  function zr1(isospin,r,z)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: r,z
    real(kind=8) :: zr1
    
    zr1 = z/r
  
  end function zr1
  
! returns z/rho
  function testfunction1(isospin,r,z)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: r,z
    real(kind=8) :: testfunction1
    
    real(kind=8) :: dst
    
    
    dst = dist_simple(r,z,shape_quadrupole_so,isospin)
    
    testfunction1 = dst
  
  end function testfunction1

! returns derivative of the spin-orbit potential 
! with respect to rho
  function vsodr(isospin,r,z)
    
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: r,z
    real(kind=8) :: vsodr
    
    real(kind=8) :: dst,edst,fact,vder,rz2,cst
    
    dst = dist_simple(r,z,shape_quadrupole_so,isospin)
    edst = exp(dst/a_so(isospin))
    fact = -edst/(a_so(isospin)*((1.0 + edst)**2))
    
    cst = r_so(isospin)*(n_mass**third)*sfactor*beta*1.5*sqrt(5.0/pi)
    rz2 = r**2 + z**2
    vder = r/sqrt(rz2) + cst*r*(z**2)/(rz2**2)
    
    if(r.lt.1e-6.and.abs(z).lt.1e-6) then
    	vsodr = 0.0
   	else
   		vsodr = fact*vder
    end if
    
  end function vsodr
  
! returns derivative of the spin-orbit potential 
! with respect to rho and divided by rho
  function vsodrrho(isospin,r,z)
    
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: r,z
    real(kind=8) :: vsodrrho
    
    real(kind=8) :: dst,edst,fact,vder,rz2,cst
    
    dst = dist_simple(r,z,shape_quadrupole_so,isospin)
    edst = exp(dst/a_so(isospin))
    fact = -edst/(a_so(isospin)*((1.0 + edst)**2))
    
    cst = r_so(isospin)*(n_mass**third)*sfactor*beta*1.5*sqrt(5.0/pi)
    rz2 = r**2 + z**2
    vder = r/sqrt(rz2) + cst*r*(z**2)/(rz2**2)
    
    vsodrrho = fact*vder/r
  
  end function vsodrrho

! returns derivative of the spin-orbit potential 
! with respect to z
  function vsodz(isospin,r,z)
    
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: r,z
    real(kind=8) :: vsodz
    
    real(kind=8) :: dst,edst,fact,vder,rz2,cst1,cst2
    
    dst = dist_simple(r,z,shape_quadrupole_so,isospin)
    edst = exp(dst/a_so(isospin))
    fact = -edst/(a_so(isospin)*((1.0 + edst)**2))
    
    cst = r_so(isospin)*(n_mass**third)*sfactor*beta*1.5*sqrt(5.0/pi)
    rz2 = r**2 + z**2
    vder = z/sqrt(rz2) - cst*z*(r**2)/(rz2**2)
    
    vsodz = fact*vder
  
  end function vsodz
  
! returns derivative of the spin-orbit potential 
! with respect to z and divided by rho
  function vsodzrho(isospin,r,z)
    
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: r,z
    real(kind=8) :: vsodzrho
    
    real(kind=8) :: dst,edst,fact,vder,rz2,cst
    
    dst = dist_simple(r,z,shape_quadrupole_so,isospin)
    edst = exp(dst/a_so(isospin))
    fact = -edst/(a_so(isospin)*((1.0 + edst)**2))
    
    cst = r_so(isospin)*(n_mass**third)*sfactor*beta*1.5*sqrt(5.0/pi)
    rz2 = r**2 + z**2
    vder = z/sqrt(rz2) - cst*z*(r**2)/(rz2**2)
    
    vsodzrho = fact*vder/r
  
  end function vsodzrho
  
! calculates the laplacian of the Woods-Saxon potential
! it is used in calculating the Coulomb potential
! see Vautherin, PRC 7, 296 (1973), p. 303
  function ddws(isospin,r,z)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: r,z
    
    real(kind=8) :: ddws
    
    real(kind=8) :: dst,edst,edsti,cst,rz2
    real(kind=8) :: ddz,d2dz,ddr,d2dr,lapr,lapz,lmid
    
    dst = dist_simple(r,z,shape_quadrupole,isospin)
    edst = exp(dst/a_c(isospin))
    edsti = 1.0/(1.0 + edst)
    cst = r_c(isospin)*(n_mass**third)*sfactor*beta*1.5*sqrt(5.0/pi)
    rz2 = r**2 + z**2
    
    ddr = r/sqrt(rz2) + cst*r*(z**2)/(rz2**2)
    d2dr = (z**2)/(rz2**1.5) + cst*(z**4-3*(r**2)*(z**2))/(rz2**3)
    ddz = z/sqrt(rz2) - cst*z*(r**2)/(rz2**2)
    d2dz = r**2/(rz2**1.5) - cst*(r**4-3*(r**2)*(z**2))/(rz2**3)
    
    lapr = (edsti**2)*edst*(edsti*(edst-1.0)*(ddr**2)-a_c(isospin)*d2dr)/(a_c(isospin)**2)
    lapz = (edsti**2)*edst*(edsti*(edst-1.0)*(ddz**2)-a_c(isospin)*d2dz)/(a_c(isospin)**2)
    lmid = -(edsti**2)*edst*ddr/(a_c(isospin)*r)
    
    if(abs(r).lt.0.0001.and.abs(z).lt.0.0001) then
      ddws = 0.0
    else
      ddws = lapr + lmid + lapz
    endif
  
  end function ddws
  
! calculates charge density
  function cdensity(isospin)
  
    integer, intent(in) :: isospin
    
    real(kind=8) :: cdensity
    
    integer :: i
    real(kind=8) :: sumr,rr,fact
    real(kind=8), dimension(ngauss_z) :: ftmp
    
    sumr = 0.0
    do i=1,ngauss_r
      rr = br*sqrt(gxr(i))
      fact = rr/(1.0+exp((rr-1.2*(n_mass**third))/a_c(itp)))
      sumr = sumr + fact*gwr(i,2)
    end do ! i
    sumr = 4.0*pi*0.5*br*br*sumr
   
    cdensity = n_proton/sumr
    
  end function cdensity

! for a given scaling factor and theta
! calculates the radius of the nuclear surface
  function shape_quadrupole(ct,isospin)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: ct
    
    real(kind=8) :: shape_quadrupole
    
    real :: y20
    
    y20 = 0.25*sqrt(5.0/pi)*(3*(ct**2)-1.0)
    shape_quadrupole = r_c(isospin)*(n_mass**third)*sfactor*(1.0+beta*y20)

  end function shape_quadrupole
  
! for a given scaling factor and theta
! calculates the radius of the nuclear surface
! for the spin-orbit potential
  function shape_quadrupole_so(ct,isospin)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: ct
    
    real(kind=8) :: shape_quadrupole_so
    
    real :: y20
    
    y20 = 0.25*sqrt(5.0/pi)*(3*(ct**2)-1.0)
    shape_quadrupole_so = r_so(isospin)*(n_mass**third)*sfactor*(1.0+beta*y20)

  end function shape_quadrupole_so
  
! calculates the distance from a point (rho,z) to 
! the equivalent point on the 
! nuclear surface defined by rshape  
  function dist_simple(rho,z,rshape,isospin)
  
    integer, intent(in) :: isospin
    real(kind=8), intent(in) :: rho,z
    
    interface
      function rshape(ct,isospin)
        integer, intent(in) :: isospin
        real(kind=8), intent(in) :: ct
        real(kind=8) :: rshape
      end function rshape
    end interface
    
    real(kind=8) :: dist_simple
    
    real(kind=8) :: r,ctheta,rsigma
    
    r = sqrt(rho**2 + z**2)
    ctheta = z/r
    rsigma = rshape(ctheta,isospin)
    
    dist_simple = r - rsigma
  
  end function dist_simple
  
! calculates matrix elements  
  function matel(iso,reta,rxi,wfr1,wfz1,wfr2,wfz2,fnc)
  
    integer, intent(in) :: iso
    real(kind=8), dimension(ngauss_r), intent(in) :: reta
    real(kind=8), dimension(ngauss_z), intent(in) :: rxi
    real(kind=8), dimension(ngauss_r), intent(in) :: wfr1,wfr2
    real(kind=8), dimension(ngauss_z), intent(in) :: wfz1,wfz2

    interface
      function fnc(isospin,r,z)
        integer, intent(in) :: isospin
        real(kind=8), intent(in) :: r,z
        real(kind=8) :: fnc
      end function fnc
    end interface

    real(kind=8) :: matel
    
    integer :: i,j,k
    real(kind=8) :: rsum,zsum,rr,rz
    real(kind=8), dimension(ngauss_z) :: ftmp
    
    rsum = 0.0
    do j=1,ngauss_r ! integration over eta coordinate
      rr = br*sqrt(reta(j)) ! transformation from eta to r
      do i=1,ngauss_z ! integration over xi coordinate
        rz = bz*rxi(i) ! transformation from xi to z
        ftmp(i) = fnc(iso,rr,rz)
      end do ! i
      zsum = bz*sum(wfz1*gwz(:,2)*ftmp*wfz2)
      rsum = rsum + wfr1(j)*gwr(j,2)*wfr2(j)*zsum
!      print '(i3,4f12.6)',j,rr,gwr(j,2),wfr1(j),wfr2(j)
!      print *,zsum,rsum
    end do ! j
    rsum = 0.5*br*br*rsum
    
    matel = rsum

  end function matel
  
! Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
! arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
! Legendre n-point quadrature formula. 
  subroutine gauleg(x1,x2,x,w,n)
  
    integer, intent(in) :: n
    real(kind=8), intent(in) :: x1,x2
    real(kind=8), dimension(n), intent(out) :: x,w
    
    integer :: m,j,i
    real(kind=8) :: z1,z,xm,xl,pp,p3,p2,p1;
    
    m = (n+1)/2
    xm = 0.5*(x2+x1)
    xl = 0.5*(x2-x1)
    do i=1,m
      z = cos(pi*(i-0.25)/(n+0.5))
      do ! infinite loop
        p1 = 1.0
        p2 = 0.0
        do j=1,n
          p3 = p2
          p2 = p1
          p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
        end do ! j
        pp = n*(z*p1-p2)/(z**2 - 1.0)
        z1 = z
        z = z1 - p1/pp
        if(abs(z1-z).lt.epsint) then
          exit
        end if
      end do
      x(i) = xm-xl*z
      x(n+1-i) = xm+xl*z
      w(i) = 2.0*xl/((1.0-z**2)*(pp**2))
      w(n+1-i) = w(i)
    end do ! m
  
  end subroutine gauleg
  
  subroutine gauher(x,w,n)
  
    integer, intent(in) :: n
    real(kind=8), dimension(n), intent(out) :: x,w
    
    integer, parameter :: nrmaxit = 10
    
    integer :: i,its,j,m
    real(kind=8) :: p1,p2,p3,pp,z,z1
    
    m = (n+1)/2
    do i=1,m
      select case(i) ! pocetne vrijednosti nultocaka
        case(1)
          z = sqrt(2.0*n+1.0)-1.85575*((2.0*n+1.0)**(-0.16667))
        case(2)
          z = z - 1.14*n**0.426/z ! provjeriti zagradu
        case(3)
          z = 1.86*z - 0.86*x(1)
        case(4)
          z = 1.91*8 - 0.91*x(2)
        case default
          z = 2.0*z - x(i-2)
      end select
      do its=1,nrmaxit
        p1 = pi4
        p2 = 0.0_8
        do j=1,n
          p3 = p2
          p2 = p1
          p1 = z*sqrt(2.0_8/j)*p2 - sqrt((j-1.0_8)/j)*p3
        end do ! j
        pp = sqrt(2.0_8*n)*p2
        z1 = z
        z = z1 - p1/pp
        if(abs(z-z1).le.nreps) then
          exit ! precision achieved
        end if
      end do ! its
      x(i) = z
      x(n+1-i) = -z
      w(i) = 2.0_8/(pp*pp)
      w(n+1-i) = w(i)
    end do ! i
  
  end subroutine gauher
  
! polynomial interpolation
  function polint(xa,ya,x)
  
    real(kind=8), intent(in) :: x
    real(kind=8), dimension(:), intent(in) :: xa,ya
    
    real(kind=8) :: polint
    
    integer :: i,j,n
    real(kind=8) :: x1,x2
    real(kind=8), dimension(size(xa)) :: p
    
    n = size(p)
    p = ya
    do i=1,n-1
      do j=1,n-i
        x1 = xa(j)
        x2 = xa(j+i)
        p(j) = ((x-x2)*p(j)+(x1-x)*p(j+1))/(x1-x2)
      end do ! j
    end do ! i
    
    polint = p(1)
  
  end function polint
  
! function used in calculating shell corrections to binding energy  
! polynomial x Gaussian
  function sdf(x)
    
    real(kind=8), intent(in) :: x
    
    real(kind=8) :: sdf
    
    integer :: i
    real(kind=8) :: w,p
    real(kind=8), dimension(4) :: a,xx
    
    a(1) = 35.0_8/16.0_8
    xx(1) = 1.0
    a(2) = -35.0_8/8.0_8
    xx(2) = x**2
    a(3) = 7.0_8/4.0_8
    xx(3) = x**4
    a(4) = -1.0_8/6.0_8
    xx(4) = x**6
    
    w = exp(-x**2)/sqrt(pi)
    p = sum(a*xx)
    
    sdf =  p*w
  
  end function sdf
  
! complete elliptic integral of the first kind
! using polynomial approximation  
  function elliptic1(m)
  
    real(kind=8), intent(in) :: m
    
    real(kind=8) :: elliptic1
    
    integer :: i
    real(kind=8) :: m1
    real(kind=8), dimension(5) :: a = (/1.38629436112_8,0.09666344259_8,0.03590092383_8, &
                                                      0.03742563713_8,0.01451196212_8/)
    real(kind=8), dimension(5) :: b = (/0.5_8,0.12498593597_8,0.06880248576_8, &
                                              0.03328355346_8,0.00441787012_8/)
    real(kind=8), dimension(5) :: m1a
    
    m1 = 1.0 - m
    m1a(1) = 1.0_8
    do i=2,5
      m1a(i) = m1a(i-1)*m1
    end do ! i
    
    elliptic1 = sum(a*m1a) + sum(b*m1a)*log(1.0/m1)
  
  end function elliptic1
  
  function elliptic2(m)
    
    real(kind=8), intent(in) :: m
    
    real(kind=8) :: elliptic2
    
    integer :: i
    real(kind=8) :: m1
    real(kind=8), dimension(4) :: a = (/0.44325141463_8,0.06260601220_8, &
                                        0.04757383546_8,0.01736506451_8/)
    real(kind=8), dimension(4) :: b = (/0.24990368310_8,0.09200180037_8, &
                                        0.04069697526_8,0.00526449639_8/)
    real(kind=8), dimension(4) :: m1a
    
    m1 = 1.0 - m
    m1a(1) = m1
    do i=2,4
      m1a(i) = m1a(i-1)*m1
    end do ! i
    
    if(abs(m1).lt.0.00001) then
      elliptic2 = 1.0
    else
      elliptic2 = 1.0 + sum(a*m1a) + sum(b*m1a)*log(1.0/m1)
    endif
    
  end function elliptic2
  
end module func
