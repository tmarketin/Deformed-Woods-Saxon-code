module eigen

  use params
  
  contains
  
  subroutine diagonalize(matrix,evalues,evectors)
  
    real(kind=8), dimension(:,:), intent(in) :: matrix
    real(kind=8), dimension(:), intent(out) :: evalues
    real(kind=8), dimension(:,:), intent(out) :: evectors
    
    logical :: t,t1
    
    integer :: ierror
    real(kind=8), dimension(size(evalues)) :: subdiag
    
    t = issym(matrix)
    if(t) then
      t1 = tred2(matrix,evalues,subdiag,evectors)
!      print *,evalues
      t1 = tql2(evalues,subdiag,evectors,ierror)
!      print *, evalues
    end if
  
  end subroutine diagonalize
  
! examines the symmetry of the matrix to be diagonalized 
  function issym(a)
  
    real(kind=8), dimension(:,:), intent(in) :: a
    
    logical :: issym
    
    integer :: n1,n2,i,j
    
    issym = .true.
    n1 = size(a(:,1))
    n2 = size(a(1,:))
    if(n1.ne.n2) then
      print *,'Dimensions of matrix a are not equal'
      issym = .false.
      return
    end if
    do i=1,n1
      do j=i,n1
        if(abs(a(i,j)-a(j,i)).gt.0.01) then
          print *,i,j
          issym = .false.
        end if
      end do ! j
    end do ! i
    
    if(issym) then
!      print *,'Matrix is symmetric.'
    else
      print *,'Matrix is not symmetric.'
    end if
  
  end function issym
  
! reduction of a real symmetric matrix to tridiagonal form
! a(input) is the symmetric matrix
! d(output) - diagonal elements of the tridiagonal matrix
! e(output) - elements on the subdiagonal
! v(output) - orthogonal transformation matrix
  function tred2(a,d,e,v)
  
    real(kind=8), dimension(:,:), intent(in) :: a
    real(kind=8), dimension(:), intent(out) :: d
    real(kind=8), dimension(:), intent(out) :: e
    real(kind=8), dimension(:,:), intent(out) :: v
    
    logical :: tred2
  
    integer :: i,j,k,l,n,jp1
    real(kind=8) :: scl,h,f,g,hh
    
    
    n = size(a(1,:))
    
    v = a
    d = v(n,:)
    
    do i=n,2,-1
      l = i - 1
    
! Scale to avoid under/overflow.
      h = 0.0
      scl = 0.0
      if(l.gt.1) then
        scl = sum(abs(d(1:l)))
      end if
      
      if(scl.eq.0.0) then
        e(i) = d(l)
        d(1:l) = v(l,1:l)
        v(i,1:l) = 0.0
        v(1:l,i) = 0.0
      else
        d(1:l) = d(1:l)/scl
        h = sum(d(1:l)*d(1:l))

        f = d(l)
        g = -sign(sqrt(h),f)
        e(i) = scl*g
        h = h - f*g
        d(l) = f - g
        
! form a*u
        e(1:l) = 0.0
        do j=1,l
          f = d(j)
          v(j,i) = f
          g = e(j) + v(j,j)*f
          jp1 = j + 1
          if(l.ge.jp1) then
            g = g + sum(v(jp1:l,j)*d(jp1:l))
            e(jp1:l) = e(jp1:l) + v(jp1:l,j)*f
          end if
          e(j) = g
        end do ! j
      
! form p
        f = 0.0
        e(1:l) = e(1:l)/h
        f = sum(e(1:l)*d(1:l))
        hh = f/(h + h)
        
! form q
        e(1:l) = e(1:l) - hh*d(1:l)
        
! form reduced a
        do j=1,l
          f = d(j)
          g = e(j)
          v(j:l,j) = v(j:l,j) - f*e(j:l) - g*d(j:l)
          d(j) = v(l,j)
          v(i,j) = 0.0
        end do ! j
        
      end if
    
      d(i) = h
    
    end do ! i
      
      do i=2,n
        l = i - 1
        v(n,l) = v(l,l)
        v(l,l) = 1.0
        h = d(i)
        if(h.ne.0.0) then
          d(1:l) = v(1:l,i)/h
          do j=1,l
            g = 0.0
            g = sum(v(1:l,i)*v(1:l,j))
            v(1:l,j) = v(1:l,j) - g*d(1:l)
          end do ! j
        end if
        
        v(1:l,i) = 0.0
      end do ! i
      
      d(1:n) = v(n,1:n)
      v(n,1:n) = 0.0
      v(n,n) = 1.0
      e(1) = 0.0
      
      tred2 = .true.
      
  end function tred2
  
!     finds dsqrt(a**2+b**2) without overflow or destructive underflow
  function pythag(a,b)
  
    real(kind=8), intent(in) :: a,b
    
    real(kind=8) :: pythag
    
    real(kind=8) :: p,r,s,t,u
    
    p = max(abs(a),abs(b))
    if(p.ne.0.0) then
      r = (min(abs(a),abs(b))/p)**2
      do
        t = 4.0 + r
        if(t.eq.4.0) then
          exit
        end if
        s = r/t
        u = 1.0 + 2.0*s
        p = u*p
        r = ((s/u)**2)*r
      end do ! infinite loop
    end if
    pythag = p
  
  end function pythag
  
! calculates eigenvalues and eigenvectors of a tridiagonal matrix
! d(inout) - diagonal elements of the matrix, on output holds the eigenvalues
! e(input) - elements on the subdiagonal
! v(inout) - transformation matrix from real, symmetric to tridiagonal, on
! output holds the eigenvectors as columns  
  function tql2(d,e,v,ierr)
  
    integer, intent(out) :: ierr
    real(kind=8), dimension(:), intent(inout) :: d
    real(kind=8),dimension(:), intent(inout) :: e
    real(kind=8), dimension(:,:), intent(inout) :: v
    
    logical :: tql2
    
    integer :: i,j,k,l,m,n,ii,l1,l2,nm,mml
    real(kind=8) :: c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2
    
    n = size(d)
    
    tql2 = .true.
    ierr = 0
    
    if(n.eq.1) then
      ierr = 1
      write(lout,*) 'Order of matrix = 1'
      return
    end if
    
    e(1:n-1) = e(2:n)
    
    f = 0.0
    tst1 = 0.0
    e(n) = 0.0
    
    do l=1,n
      j = 0
      tst1 = max(tst1,abs(d(l)) + abs(e(l)))
      
      do m=l,n
        tst2 = tst1 + abs(e(m))
        if(abs(e(m)).lt.tst1*1.0e-8) then
          exit
        end if
      end do ! m
      
      if(m.gt.l) then
        do
        
! check convergence
          if(j.eq.30) then
            tql2 = .false.
            print *,'Number of iteration greater than 30.'
            return
          end if
          j = j + 1
          
! form shift
          l1 = l + 1
          l2 = l1 + 1
          g = d(l)
          p = (d(l1) - g)/(2.0*e(l))
          r = pythag(p,1.0_8)
          if(p.lt.0.0) then ! replaces sign function in expressions
            r = -r
          end if
          d(l) = e(l)/(p + r)
          d(l1) = e(l)*(p + r)
          dl1 = d(l1)
          h = g - d(l)
          if(l2.le.n) then
            d(l2:n) = d(l2:n) - h
          end if
          f = f + h
          
! ql transformation
          p = d(m)
          c = 1.0
          c2 = c
          el1 = e(l1)
          s = 0.0
          s2 = 0.0
          mml = m - l
          
          do ii=1,mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c*e(i)
            h = c*p
            r = pythag(p,e(i))
            e(i+1) = s*r
            s = e(i)/r
            c = p/r
            p = c*d(i) - s*g
            d(i+1) = h + s*(c*g + s*d(i))
! form vector
            do k=1,n
              h = v(k,i+1)
              v(k,i+1) = s*v(k,i) + c*h
              v(k,i) = c*v(k,i) - s*h
            end do ! k
          end do ! ii
          
          p = -s*s2*c3*el1*e(l)/dl1
          e(l) = s*p
          d(l) = c*p
          tst2 = tst1 + abs(e(l))
          if(tst2.le.tst1) then
            exit ! infinite loop exit
          end if
        end do ! infinite loop
      
      end if
      
      d(l) = d(l) + f
      e(l) = 0.0
    end do ! l
    
! sorting eigenvalues
    do i=1,n-1
      k = i
      p = d(i)
      do j=i+1,n
        if(d(j).le.p) then
          k = j
          p = d(j)
        end if
      end do ! j
      
      if(k.ne.i) then
        d(k) = d(i)
        d(i) = p
        do j=1,n
          p = v(j,i)
          v(j,i) = v(j,k)
          v(j,k) = p
        end do ! j
      end if
    end do ! i
  
  end function tql2
  
end module eigen