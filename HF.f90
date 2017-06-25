module HF
implicit none
contains

!!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!!
  subroutine rhf(en,H,V,T,S,eri,enuc,C_ao,en_mo,nao,nel)
    integer :: nao,nel,i
    real*8  :: enuc,en,val
    real*8,dimension(nao)             :: en_mo
    real*8,dimension(nao,nao)         :: S,O,H,V,T,F,C_ao,C_mo,P0,Pn
    real*8,dimension(nao,nao,nao,nao) :: eri
    real*8,parameter :: delta1=0.000000000001, delta2=0.00000000001
    integer :: maxiter
    logical :: converged,debug
    debug = .false.
    en = 0.
    ! costruct approximated Hamiltonian
    H = V + T
    ! diagonalize overlap matrix
    call constr_orth(S,O,nao)
    ! build coefficent matrix and transform it to AO basis
    call constr_coef(C_ao,C_mo,H,O,nao,en_mo)
    ! build inital density
    call constr_dens(P0,P0,C_ao,val,nao,nel/2,.false.)
    ! calculate inital SCF energy
    call scfenergy(en,H,H,P0,nao)
    en = enuc + en
    write(*,*)
    write(*,*) 0,en
    ! start SCF iteration
    call scf(en,enuc,H,F,O,C_ao,C_mo,en_mo,P0,Pn,eri,nao,nel/2,100,converged)
    if(converged) then
      write(*,*) "Hurray! SCF is converged!"
    else
      write(*,*) "Warning: SCF is not converged!"
    endif
  end subroutine rhf

!!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!!
  subroutine scf(enn,enuc,H,F,O,C_ao,C_mo,en_mo,P0,Pn,eri,nao,nocc,maxiter,converged)
    integer :: m,nao,nocc
    real*8  :: enuc,rmsP,deltaE,enn,eno
    real*8,dimension(nao)             :: en_mo
    real*8,dimension(nao,nao)         :: O,H,F,C_ao,C_mo,P0,Pn
    real*8,dimension(nao,nao,nao,nao) :: eri
    real*8,parameter :: delta1 = 0.000000000001,delta2=0.00000000001
    integer :: maxiter
    logical :: converged,debug
    debug = .false.

    write(*,*) "* Start SCF iteration"
    converged = .false.
    scfloop: do m = 1, maxiter
      ! Calculate new Fock matrix
      call constr_fock(F,H,P0,eri,nao)
      ! Diagonalize Fock Matrix
      call constr_coef(C_ao,C_mo,F,O,nao,en_mo)
      ! Build new density matrix
      call constr_dens(Pn,P0,C_ao,rmsP,nao,nocc,.true.)
      P0 = Pn
      ! Compute new SCF energy
      call scfenergy(enn,H,F,Pn,nao)
      enn = enn + enuc
      deltaE = abs(enn - eno)
      write(*,*) m, enn, deltaE, rmsP
      eno = enn
      ! Check for convergence
      if(( deltaE < delta1 ).and.( rmsP < delta2 ))then
        converged = .true.
        exit scfloop
      endif
    enddo scfloop
  end subroutine scf

!!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!!
  subroutine constr_dens(P,P0,C,rmsP,nao,nocc,lscf)
    integer :: nao,nocc,i,j,k
    real*8  :: rmsP
    real*8,dimension(nao,nao) :: P,P0,C
    logical :: lscf
    rmsP = 0.
    do i = 1, nao
      do j = i, nao
        P(i,j) = 0.
        do k = 1, nocc
          P(i,j) = P(i,j) + ( C(i,k) * C(j,k) )
        enddo
        if (i.ne.j) P(j,i) = P(i,j)
        if (lscf) rmsP = rmsP + ( P(i,j) - P0(i,j) )**2
      enddo
    enddo
    if (lscf) rmsP = sqrt(rmsP)
  end subroutine constr_dens

!!!------------------------------------------------------------------------!!!
  subroutine constr_coef(C_ao,C_mo,F,O,nao,moen)
    integer :: nao
    real*8, dimension(nao,nao) :: C_ao,C_mo,F,O,tmp1,tmp2
    real*8 :: moen(nao)
    tmp1 = matmul(transpose(O),F)
    tmp2 = matmul(tmp1,O)
    call eigv(tmp2,C_mo,nao,moen)
    C_ao = matmul(O,C_mo)
  end subroutine constr_coef

!!!------------------------------------------------------------------------!!!
  subroutine constr_orth(S,O,nao)
    integer :: nao,i,j
    real*8  :: S(nao,nao),O(nao,nao)
    real*8  :: invS(nao,nao),diaS(nao,nao),eigs(nao)
    ! diagonalize overlap matrix
    call eigv(S,diaS,nao,eigS)
    ! build symmetric orthogonalization matrix
    invS = transpose(diaS)
    do i = 1, nao
    eigS(i) = 1/sqrt(eigS(i))
      do j = 1, nao
        invS(i,j) = invS(i,j) * eigS(i)
      enddo
    enddo
    O = matmul(diaS,invS)
  end subroutine constr_orth

!!!------------------------------------------------------------------------!!!
  subroutine constr_fock(F,H,P,eri,d)
    integer, intent(in) :: d
    real*8, intent(in)  :: H(d,d),P(d,d),eri(d,d,d,d)
    real*8, intent(out) :: F(d,d)
    integer             :: i,j,k,l
    do i = 1, d
      do j = i, d
        F(i,j) = H(i,j)
        do k = 1, d
          do l = 1, d
            F(i,j) = F(i,j) + P(k,l) * ( 2*eri(i,j,k,l) - eri(i,k,j,l) )
          enddo
        enddo
        if(i/=j) F(j,i) = F(i,j)
      enddo
    enddo
  end subroutine constr_fock

!!!------------------------------------------------------------------------!!!
  subroutine scfenergy(en,H,F,P,d)
  implicit none
    integer, intent(in) :: d
    real*8, intent(in)  :: H(d,d),F(d,d),P(d,d)
    real*8, intent(out) :: en
    integer             :: i,j
    en = 0.0
    do i = 1, d
      do j = 1, d
        en = en + P(i,j) * ( H(i,j) + F(i,j) )
      enddo
    enddo
  end subroutine scfenergy

!!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!!!
  subroutine eigv(matrix,dmatrix,d,lambda)
  implicit none
    integer,intent(in)    :: d
    real*8, intent(inout) :: matrix(d,d), dmatrix(d,d), lambda(d)
    integer,parameter     :: lwmax = 10000
    real*8                :: work(lwmax)
    integer               :: lwork, info
    dmatrix = matrix
    lwork = -1
    call dsyev('V','U',d,dmatrix,d,lambda,work,lwork,info)
    lwork = min(lwmax,int(work(1))) * 2*d

    call dsyev('V','U',d,dmatrix,d,lambda,work,lwork,info)

    if (info.ne.0) write(*,*) "Error: eigenvalue problem not solvable"
  end subroutine eigv
end module HF
