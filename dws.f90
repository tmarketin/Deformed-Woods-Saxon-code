program dfs

  use params
  use wavefunction
  use eigen
  use phys
  use func

  integer :: i,j,k,n,cnt
  real(kind=8), dimension(:), allocatable :: hene

  open(unit=lplot,file="energy_plot_neutron.dat",status="unknown")
  open(unit=lplot2,file="energy_plot_proton.dat",status="unknown")
  
  do i=1,41
    beta = -0.40_8 + (i-1)*0.02_8
    print '(a24,f5.2)','Deformation parameter = ',beta
    
    open(unit=lout,file="output.dat",status="replace")
    
  ! initializes parameters  
    call initialize
    
  ! initialize wavefunctions and basis
    call initw  
    
  ! calculate Hamiltonian matrix per omega block
    call hamiltonian(basis,nbasis,ham)
  
  ! calculate eigenenergies and eigenvectors
    call solve(nbasis,ham,energy,evectors)  
    
  ! fill eigenstates with nucleons and sort by energy
    call fill_orbits(nbasis,energy,evectors,energy_full)
    
    close(lout)    
    n = sum(nbasis)
    allocate(hene(n))
    cnt = 1
    do k=1,nomega
      hene(cnt:cnt + nbasis(k)) = energy(1,k,1:nbasis(k))
      cnt = cnt + nbasis(k)
    end do ! k
    write(lplot,'(f5.2,1000f12.6)') beta,hene
    
    cnt = 1
    do k=1,nomega
      hene(cnt:cnt + nbasis(k)) = energy(2,k,1:nbasis(k))
      cnt = cnt + nbasis(k)
    end do ! k
    write(lplot2,'(f5.2,1000f12.6)') beta,hene
    
    deallocate(hene)
  
  end do ! i
  
  close(lplot)
  close(lplot2)
  
end program dfs
