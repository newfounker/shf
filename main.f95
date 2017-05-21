module mathsub
implicit none
contains
function dot(a, b, d) result(dotproduct)
    integer :: i, d
    double precision, dimension(d), intent(in) :: a, b
    double precision :: dotproduct

    dotproduct = 0
    do i = 1, d
        dotproduct = dotproduct + a(i)*b(i)
    enddo
end function dot
function inv(matrix,d) result(inverse)
    integer :: d
    double precision, dimension(d,d), intent(in) :: matrix
    double precision, dimension(d,d) :: inverse

    double precision, dimension(d) :: work 
    integer, dimension(d) :: ipiv  
    integer :: info

    inverse = matrix

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(d, d, inverse, d, ipiv, info)

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(d, inverse, d, ipiv, work, d, info)
    inverse = inverse
end function inv
end module mathsub

!!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!!

program main
use, intrinsic :: iso_fortran_env
implicit none
    character(len=50) :: dirname,line
    integer :: nargs,i,j,k,l,m,nbase,stat,lcount
    double precision  :: enuc,val,rmsD,deltaE
    double precision,dimension(:),allocatable       :: diagS,en
    double precision,dimension(:),allocatable       :: energy
    double precision,dimension(:,:),allocatable     :: S,T,V,H,F,C,tranS,itraS,work
    double precision,dimension(:,:,:),allocatable   :: D
    double precision,dimension(:,:,:,:),allocatable :: eri
    double precision,parameter :: ONE = 1.0, ZERO = 0.0,&
        delta1 = 0.000000000001,delta2=0.00000000001
    integer,parameter :: maxiter=100
    logical :: converged

    write(*,*) "* Fortran test program *"
    write(*,*)

! Check if command line argument is present
    nargs = command_argument_count()
    if (nargs == 0) then
        write(*,*) "Usage: hf-scf dirname"
        ! TODO: exit program
    endif

! Get dirname from commandline argument
    call get_command_argument(1,dirname)
    write(*,*) "* Reading nuclear repulsion energy from: ", &
        trim(dirname)//'enuc.dat'
    open(1,file=trim(dirname)//'enuc.dat',status='old')
    read(1,*) enuc
    close(1)

!!!------------------------------------------------------------------------!!!
    write(*,*) "* Reading overlap integrals from: ", trim(dirname)//'s.dat'
    open(2,file=trim(dirname)//'s.dat',status='old')
    lcount = 0
    nbase = 0
    countloop: do
    read (2,*,iostat=stat) i,j,line
    if ( stat /= 0 ) then
        if ( stat == iostat_end ) then
            exit countloop
        else
            write(*,*) stat
            stop
        endif
    endif
    lcount= lcount+1
    nbase = max(nbase,i,j)
    enddo countloop
    allocate(S(nbase,nbase))
    rewind(2)
    do i = 1, lcount
        read(2,*) j,k,val
        S(j,k) = val
        if(j/=k) S(k,j) = S(j,k)
    enddo
    close(2)

!!!------------------------------------------------------------------------!!!
    write(*,*) "* Reading kinetic energy from: ", trim(dirname)//'t.dat'
    open(3,file=trim(dirname)//'t.dat',status='old')
    allocate(T(nbase,nbase))
    do i = 1, lcount
        read(3,*) j,k,val
        T(j,k) = val
        if(k/=j) T(k,j) = T(j,k)
    enddo
    close(3)

!!!------------------------------------------------------------------------!!!
    write(*,*) "* Reading nuclear attraction integrals from: ", &
        trim(dirname)//'v.dat'
    open(4,file=trim(dirname)//'v.dat',status='old')
    allocate(V(nbase,nbase))
    do i = 1, lcount
        read(4,*) j,k,val
        V(j,k) = val
        if(k/=j) V(k,j) = V(j,k)
    enddo
    close(4)

!!!------------------------------------------------------------------------!!!
    write(*,*) "* Constructing approximated Hamiltonian..."
    allocate(H(nbase,nbase))
    do i = 1, nbase
    do j = i, nbase
        H(i,j) = T(i,j) + V(i,j)
        if(i/=j) H(j,i) = H(i,j)
    enddo
  !     write(*,*) H(i,:)
    enddo
    deallocate(T)
    deallocate(V)

!!!------------------------------------------------------------------------!!!
    write(*,*) "* Reading two electron repulsion integrals from: ", &
        trim(dirname)//'eri.dat'
    allocate(eri(nbase,nbase,nbase,nbase))
    open(5,file=trim(dirname)//'eri.dat',status='old')
    lcount = 0
    countloop2: do
    read (5,*,iostat=stat) i,j,line
    if ( stat /= 0 ) then
        if ( stat == iostat_end ) then
            exit countloop2
        else
            write(*,*) stat
            stop
        endif
    endif
    lcount= lcount+1
    enddo countloop2
    rewind(5)
    do i = 1, lcount
        read(5,*) j,k,l,m,val
        eri(j,k,l,m) = val
        eri(k,j,l,m) = val
        eri(j,k,m,l) = val
        eri(k,j,m,l) = val
        eri(l,m,j,k) = val
        eri(m,l,j,k) = val
        eri(l,m,k,j) = val
        eri(m,l,k,j) = val
    enddo
    close(5)

!!!------------------------------------------------------------------------!!!
    write(*,*) "* Diagonalizing of overlap matrix..."
    allocate(tranS(nbase,nbase))
    allocate(itraS(nbase,nbase))
    allocate(diagS(nbase))
    call eigenvectors(S,tranS,nbase,diagS)
    write(*,*) "* Build symmetric orthoganalization matrix..."

    itraS = transpose(tranS)
    do i = 1, nbase
        diagS(i) = 1/sqrt(diagS(i))
        do j = 1, nbase
            itraS(i,j) = itraS(i,j) * diagS(i)
        enddo
    enddo
    call dgemm('N','N',nbase,nbase,nbase,ONE,tranS,nbase,itraS,nbase,ZERO,S,nbase)
  ! do i = 1, nbase
  !     write(*,*) S(i,:)
  ! enddo
    deallocate(tranS)
    deallocate(itraS)
    deallocate(diagS)

!!!------------------------------------------------------------------------!!!
    write(*,*) "* Build inital Fock matrix..."
    allocate(F(nbase,nbase))
    allocate(work(nbase,nbase))
    call dgemm('T','N',nbase,nbase,nbase,ONE,S,nbase,H,nbase,ZERO,work,nbase)
    call dgemm('N','N',nbase,nbase,nbase,ONE,work,nbase,S,nbase,ZERO,F,nbase)
    deallocate(work)
  ! do i = 1, nbase
  !     write(*,*) F(i,:)
  ! enddo

!!!------------------------------------------------------------------------!!!
    write(*,*) "* Diagonalize Fock matrix..."
    allocate(C(nbase,nbase))
    allocate(en(nbase))
    call eigenvectors(F,C,nbase,en)
  ! write(*,*) en

!!!------------------------------------------------------------------------!!!
    write(*,*) "* Transform eigenvectors into AO-basis..."
    allocate(work(nbase,nbase))
    call dgemm('N','N',nbase,nbase,nbase,ONE,S,nbase,C,nbase,ZERO,work,nbase)
    C = transpose(work)
    deallocate(work)
  ! do i = 1, nbase
  !     write(*,*) C(i,:)
  ! enddo
    
!!!------------------------------------------------------------------------!!!
    write(*,*) "* Build initial density..."
    allocate(D(0:maxiter,nbase,nbase))
    do i = 1, nbase
    do j = i, nbase
        D(0,i,j) = 0.0
        do k = 1, nbase
            D(0,i,j) = D(0,i,j) + ( C(i,k) * C(j,k) )
        enddo
        if(i/=j) D(0,j,i) = D(0,i,j)
    enddo
    enddo
  ! do i = 1, nbase
  !     write(*,*) D(0,i,:)
  ! enddo

    write(*,*)
!!!------------------------------------------------------------------------!!!
    write(*,*) "* Calculate the inital SCF energy"
    allocate(energy(0:maxiter))

    energy(0) = enuc
    do i = 1, nbase
    do j = 1, nbase
        energy(0) = energy(0) + D(0,i,j) * ( H(i,j) + F(i,j) )
    enddo
    enddo
    write(*,*) 0, energy(0)

!!!------------------------------------------------------------------------!!!
    write(*,*) "* Start SCF iteration"
    scfloop: do m = 1, maxiter
    ! Calculate new Fock matrix
        do i = 1, nbase
        do j = i, nbase
            F(i,j) = H(i,j)
            do k = j, nbase
            do l = k, nbase
                F(i,j) = F(i,j) + D(m-1,i,j) * ( 2*eri(i,j,k,l) - eri(i,k,j,l) )
            enddo
            enddo
            if(i/=j) F(j,i) = F(i,j)
        enddo
        enddo
    ! Orthogonalize Fock Matrix
        allocate(work(nbase,nbase))
        call dgemm('T','N',nbase,nbase,nbase,ONE,S,nbase,F,nbase,ZERO,work,nbase)
        call dgemm('N','N',nbase,nbase,nbase,ONE,work,nbase,S,nbase,ZERO,F,nbase)
        deallocate(work)
    ! Diagonalize Fock Matrix
        call eigenvectors(F,C,nbase,en)
    ! Back-transform AO-Basis
        allocate(work(nbase,nbase))
        call dgemm('N','N',nbase,nbase,nbase,ONE,S,nbase,C,nbase,ZERO,work,nbase)
        C = transpose(work)
        deallocate(work)
    ! Build new density matrix
        rmsD = ZERO
        do i = 1, nbase
        do j = i, nbase
            D(m,i,j) = ZERO
            do k = 1, nbase
                D(m,i,j) = D(m,i,j) + ( C(i,k) * C(j,k) )
            enddo
            if(i/=j) D(m,j,i) = D(m,i,j)
            rmsD = rmsD + ( D(m,i,j) - D(m-1,i,j) )**2
        enddo
        enddo
        rmsD = sqrt(rmsD)
    ! Compute new SCF energy
        energy(m) = enuc
        do i = 1, nbase
        do j = 1, nbase
            energy(m) = energy(m) + D(m,i,j) * ( H(i,j) + F(i,j) )
        enddo
        enddo
        deltaE = energy(m)-energy(m-1)
        write(*,*) m, energy(m), deltaE, rmsD
    ! Check for convergence
        if(( deltaE < delta1 ).and.( rmsD < delta2 ))then
            write(*,*) "Hurray! SCF is converged!"
            converged = .true.
            exit scfloop
        endif
        if(m==maxiter) then
            write(*,*) "Warning: SCF is not converged!"
            converged = .false.
        endif
    enddo scfloop

!!!------------------------------------------------------------------------!!!
    deallocate(S)
    deallocate(H)
    deallocate(F)
    deallocate(C)
    deallocate(D)
    deallocate(en)
    deallocate(eri)
    deallocate(energy)
end program main

!!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!!

subroutine eigenvalues(matrix,d,eigv)
implicit none
    integer,intent(in)  :: d
    double precision, intent(inout) :: matrix(d,d), eigv(d)
    double precision, allocatable   :: work(:)
    integer             :: lwork, info

    lwork = max(1,3*d-1)
    allocate(work(lwork))

    call dsyev('N','U',d,matrix,d,eigv,work,lwork,info)
    deallocate(work)
end subroutine eigenvalues
subroutine eigenvectors(matrix,dmatrix,d,eigv)
implicit none
    integer,intent(in)  :: d
    double precision, intent(inout) :: matrix(d,d),dmatrix(d,d), eigv(d)
    double precision, allocatable   :: work(:)
    integer             :: lwork, info
    dmatrix = matrix

    lwork = max(1,3*d-1)
    allocate(work(lwork))

    call dsyev('V','U',d,dmatrix,d,eigv,work,lwork,info)
    deallocate(work)
end subroutine eigenvectors
