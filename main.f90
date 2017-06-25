subroutine get_atomnumber(filename,nat)
implicit none
    character(len=50) :: filename
    integer :: nat
    write(*,*) "* Reading geometry from: ", trim(filename)
    open(0,file=trim(filename),status='old')
    read(0,*) nat
end subroutine get_atomnumber
subroutine get_geometry(nat,nel,atom,x,y,z)
implicit none
    integer :: nat,i
    integer :: nel,atom(nat)
    real*8  :: x(nat),y(nat),z(nat)
  ! open(0,file=trim(filename),status='old')
  ! read(0,*) natoms
    nel = 0
    do i = 1, nat
        read(0,*) atom(i), x(i), y(i), z(i)
        nel = nel + atom(i)
    enddo
    close(0)
end subroutine get_geometry

program soulhf
use, intrinsic :: iso_fortran_env
use HF
use MP
use CC
implicit none
    character(len=50) :: dirname,line
    integer :: nargs,i,j,k,l,m,n,nbase,stat,lcount,nelec,natoms
    real*8  :: enuc,val,rmsP,deltaE,enn,eno,toten,scfen,mp2en,encc
    integer,dimension(:),allocatable      :: atom
    real*8,dimension(:),allocatable       :: x,y,z,diagS,en
    real*8,dimension(:,:),allocatable     :: S,T,V,H,F,C,P0,Pn,&
        tranS,itraS,work,tmp1,tmp2
    real*8,dimension(:,:,:,:),allocatable :: eri
    real*8,parameter :: ONE = 1.0, ZERO = 0.0,&
        delta1 = 0.000000000001,delta2=0.00000000001
    integer,parameter :: maxiter=100
    logical :: converged,debug
    debug = .true.

    write(*,*) "* Fortran test program *"

! Check if command line argument is present
    nargs = command_argument_count()
    if (nargs == 0) then
        write(*,*) "Usage: hf-scf dirname"
        write(*,*)
        stop
    endif
    write(*,*)

! Get dirname from commandline argument
    call get_command_argument(1,dirname)

!!!------------------------------------------------------------------------!!!
    call get_atomnumber(trim(dirname)//'geom.dat',natoms)
    allocate(atom(natoms))
    allocate(x(natoms))
    allocate(y(natoms))
    allocate(z(natoms))
    call get_geometry(natoms,nelec,atom,x,y,z)

!!!------------------------------------------------------------------------!!!
    write(*,*) "* Reading nuclear repulsion energy from: ", &
        trim(dirname)//'enuc.dat'
    open(1,file=trim(dirname)//'enuc.dat',status='old')
    read(1,*) enuc
    if(debug) write(*,*) enuc
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
    if(debug) then
        do i = 1, nbase
            write(*,*) S(i,:)
        enddo
    endif

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
    if(debug) then
        do i = 1, nbase
            write(*,*) T(i,:)
        enddo
    endif

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
    if(debug) then
        do i = 1, nbase
            write(*,*) V(i,:)
        enddo
    endif

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
    do i= 1, nbase
        do j= 1, nbase
            do k= 1, nbase
                do l= 1, nbase
                    eri(i,j,k,l) = ZERO
                enddo
            enddo
        enddo
    enddo
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
    allocate(C(nbase,nbase))
    allocate(H(nbase,nbase))
    allocate(en(nbase))
    call rhf(scfen,H,V,T,S,eri,enuc,C,en,nbase,nelec)
    toten = scfen

    write(*,*)
    write(*,*) "SCF energy: ", toten

!!!------------------------------------------------------------------------!!!
    call mp_two(mp2en,eri,C,en,nbase,nelec/2)
    toten = toten + mp2en
    write(*,*) "MP2 energy correction: ", mp2en
    write(*,*) "MP2 energy: ", toten
    
    call ccsd(encc,mp2en,H,C,eri,nelec,nbase)
    toten = toten + encc
    write(*,*)
    write(*,*) "CCSD energy:", toten
!!!------------------------------------------------------------------------!!!
    deallocate(H,V,T,S,eri,C,en)
    deallocate(atom,x,y,z)
end program soulhf
