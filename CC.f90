module cc
implicit none
contains
  subroutine ccsd(encc,enmp2,H,C_mo,eri_mo,nel,nao)
    integer :: nao,nel,i
    real*8  :: H(nao,nao),eri_mo(nao,nao,nao,nao),C_mo(nao,nao)
    real*8  :: F_so(nao*2,nao*2),eri_so(nao*2,nao*2,nao*2,nao*2)
    real*8  :: encc,enmp2
    logical :: converged 

    call spinmo(eri_so,eri_mo,nao,nao*2)
    call spinF(F_so,H,C_mo,eri_so,nao,nao*2,nel)
    do i = 1, nao*2
      write(*,*) F_so(i,:)
    enddo

    call ccsd_scf(encc,enmp2,F_so,eri_so,nel,nao*2,converged)
    write(*,*) "leaving ccsd"
    if (converged) then
      write(*,*) "Hurray! CC calculation converged!"
    else
      write(*,*) "Warning! CC is not converged!"
    endif
  end subroutine ccsd

  subroutine ccsd_scf(enncc,enmp2,F,eri,nocc,nso,converged)
    integer :: i,j,k,l,m
    integer :: nocc,nso
    integer,parameter :: maxiter=100
    real*8  :: enmp2,enocc,enncc,delta_E
    real*8,parameter  :: delta3 = 0.000000000001
    real*8,dimension(nocc+1:nso,nocc) :: ts,tsn,D2
    real*8,dimension(nocc+1:nso,nocc+1:nso,nocc,nocc) :: td,tdn,D4
    real*8,dimension(nso,nso) :: F
    real*8,dimension(nso,nso,nso,nso) :: eri
    real*8  :: F_ae(nocc+1:nso,nocc+1),F_mi(nocc,nocc),F_me(nocc,nocc+1:nso)
    real*8  :: W_mnij(nocc,nocc,nocc,nocc),W_abef(nocc+1:nso,nocc+1:nso,nocc+1:nso,nocc+1:nso),&
    &       W_mbej(nocc,nocc+1:nso,nocc,nocc+1:nso)
    logical :: converged
    D2 = 0. ! denominater array from eq. 12 and inital guess T1; Stanton et al. J. Chem. Phys. 94 (1990)
    do i = 1, nocc
      do k = nocc+1, nso
        D2(k,i) = F(i,i) - F(k,k)
        ts(k,i) = F(k,i) / D2(k,i)
        write(5,*) k,i,D2(k,i),ts(k,i)
      enddo
    enddo
!   ts = 0.

    D4 = 0. ! denominater array from eq. 13 and inital guess T2; Stanton et al. J. Chem. Phys. 94 (1990)
    do i = 1, nocc
      do j = 1, nocc
        do k = nocc+1, nso
          do l = nocc+1, nso
            D4(k,l,i,j) = F(i,i) + F(j,j) - F(k,k) - F(l,l)
            td(k,l,i,j) = td(k,l,i,j) + eri(i,j,k,l) / D4(k,l,i,j)
            write(7,*) k,l,i,j,D4(k,l,i,j),td(k,l,i,j)
          enddo
        enddo
      enddo
    enddo

    ! set inital energy to MP2 energy
    enocc = enmp2
    write(*,*) 0, enocc

    converged=.false.
    ccloop: do m = 1, maxiter
      call constr_intm(F,eri,ts,td,F_ae,F_mi,F_me,W_mnij,W_abef,W_mbej,nocc,nso)
      call constr_T1(tsn,ts,td,D2,F,eri,F_ae,F_mi,F_me,nocc,nso)
      do i = nocc+1,nso
        write(*,*) ts(i,:)
      enddo
      write(*,*)
      ts = tsn
      do i = nocc+1,nso
        write(*,*) ts(i,:)
      enddo
      call constr_T2(tdn,ts,td,D4,F,eri,F_ae,F_mi,F_me,W_mnij,W_abef,W_mbej,nocc,nso)
      td = tdn
      call correnergy(enncc,ts,td,F,eri,nocc,nso)
      delta_E = abs(enncc-enocc)
      enocc = enncc
      write(*,*) m, enocc, delta_E
      exit ccloop
      if (delta_E < delta3) then
        converged=.true.
        exit ccloop
      endif
    enddo ccloop
    write(*,*) "leaving loop"
  end subroutine ccsd_scf

  subroutine correnergy(en,ts,td,F,eri,nocc,nso)
    integer :: a,b,i,j
    integer :: nocc,nso
    real*8  :: eri(nso,nso,nso,nso),ts(nocc+1:nso,nocc),td(nocc+1:nso,nocc+1:nso,nocc,nocc),F(nso,nso)
    real*8  :: en,en_
    real*8,parameter :: one_half = 0.
    en_ = 0.
    do i = 1, nocc
      do a = nocc+1, nso
        en_ = en_ + F(i,a) * ts(a,i)
        do j = 1, nocc
          do b = nocc+1, nso
            en_ = en_ + one_half * eri(i,j,a,b) * ( one_half * td(a,b,i,j) + ts(a,i) * ts(b,j) )
          enddo
        enddo
      enddo
    enddo
    en = en_
  end subroutine correnergy

  subroutine constr_T1(tsn,ts,td,D,Fock,eri,F_ae,F_mi,F_me,nocc,nso)
    integer :: a,b,e,f
    integer :: m,n,i,j
    integer :: nocc,nso
    real*8,dimension(nocc+1:nso,nocc) :: tsn,ts,D
    real*8,dimension(nocc+1:nso,nocc+1:nso,nocc,nocc) :: td
    real*8 :: Fock(nso,nso), eri(nso,nso,nso,nso),ts_
    real*8 :: F_ae(nocc+1:nso,nocc+1),F_mi(nocc,nocc),F_me(nocc,nocc+1:nso)
    real*8,parameter :: one=1.,one_half=0.5, one_forth=0.25
    tsn = 0.
    do a = nocc+1, nso
      do i = 1, nocc
        ts_ = Fock(i,a)
        do e = nocc+1, nso
          ts_ = ts_ + ts(e,i) * F_ae(a,e)
        enddo
        do m = 1, nocc
          ts_ = ts_ - ts(a,m) * F_mi(m,i)
          do e = nocc+1, nso
            ts_ = ts_ + td(a,e,i,m) * F_me(m,e)
            do f = nocc+1, nso
              ts_ = ts_ - one_half * td(e,f,i,m) * eri(m,a,e,f)
            enddo
            do n = 1, nocc
              ts_ = ts_ - one_half * td(a,e,m,n) * eri(n,m,e,i)
            enddo
          enddo
        enddo
        do n = 1, nocc
          do f = nocc+1, nso
            ts_ = ts_ - ts(f,n) * eri(n,a,i,f)
          enddo
        enddo
        tsn(a,i) = ts_ / D(a,i)
      enddo
    enddo
  end subroutine constr_T1

  subroutine constr_T2(tdn,ts,td,D,Fock,eri,F_ae,F_mi,F_me,W_mnij,W_abef,W_mbej,nocc,nso)
    integer :: a,b,e,f
    integer :: m,n,i,j
    integer :: nocc,nso
    real*8,dimension(nocc+1:nso,nocc) :: ts
    real*8,dimension(nocc+1:nso,nocc+1:nso,nocc,nocc) :: tdn,td,D
    real*8 :: Fock(nso,nso), eri(nso,nso,nso,nso),td_
    real*8 :: F_ae(nocc+1:nso,nocc+1),F_mi(nocc,nocc),F_me(nocc,nocc+1:nso)
    real*8 :: W_mnij(nocc,nocc,nocc,nocc),W_abef(nocc+1:nso,nocc+1:nso,nocc+1:nso,nocc+1:nso),&
    &      W_mbej(nocc,nocc+1:nso,nocc,nocc+1:nso)
    real*8,parameter :: one=1.,one_half=0.5, one_forth=0.25
    tdn = 0.
    do a = nocc+1, nso
      do b = nocc+1, nso
        do i = 1, nocc
          do j = 1, nocc
            td_ = eri(i,j,a,b)
            do e = nocc+1, nso
              td_ = td_ + td(a,e,i,j) * F_ae(b,e) - td(b,e,i,j) * F_ae(a,e) &
                        + ts(e,i) * eri(a,b,e,j)  - ts(e,j) * eri(a,b,e,i)
              do m = 1, nocc
                td_ = td_ - one_half * td(a,e,i,j) * F_me(m,e) * ( ts(b,m) - ts(a,m) )
              enddo
              do f = nocc+1, nso
                td_ = td_ + one_half * W_abef(a,b,e,f) &
                &           * ( td(e,f,i,j) + ts(e,i)*ts(f,j) - ts(e,j)*ts(f,i) )
              enddo
            enddo
            do m = 1, nocc
              td_ = td_ - td(a,b,i,m) * F_mi(m,j) + td(a,b,j,m) * F_mi(m,i) &
              &         - ts(a,m) * eri(m,b,i,j)  + ts(b,m) * eri(m,a,i,j)
              do e = nocc+1, nso
                td_ = td_ - one_half * td(a,b,i,m) * F_me(m,e) * ( ts(e,j) - ts(e,i) )
                td_ = td_ + td(a,e,i,m) * W_mbej(m,b,e,j) - td(a,e,j,m) * W_mbej(m,b,e,i) &
                &         + td(b,e,i,m) * W_mbej(m,a,e,j) - td(b,e,j,m) * W_mbej(m,a,e,i)
                td_ = td_ - ts(e,i)*ts(a,m) * eri(m,b,e,j) + ts(e,j)*ts(a,m) * eri(m,b,e,i) &
                &         - ts(e,i)*ts(b,m) * eri(m,a,e,j) + ts(e,j)*ts(b,m) * eri(m,a,e,i)
              enddo
              do n = 1, nocc
                td_ = td_ + one_half * W_mnij(m,n,i,j) &
                &           * ( td(a,b,m,n) + ts(a,m)*ts(b,n) - ts(b,m)*ts(a,n) )
              enddo
            enddo
            tdn(a,b,i,j) = td_ / D(a,b,i,j)
          enddo
        enddo
      enddo
    enddo
  end subroutine constr_T2

  subroutine constr_intm(Fock,eri,ts,td,F_ae,F_mi,F_me,W_mnij,W_abef,W_mbej,nocc,nso)
    integer :: m,n,i,j
    integer :: a,b,e,f
    integer :: nocc,nso
    real*8  :: delta
    real*8  :: Fock(nso,nso),eri(nso,nso,nso,nso),&
    &       ts(nocc+1:nso,nocc), td(nocc+1:nso,nocc+1:nso,nso,nso)
    real*8  :: F_,F_ae(nocc+1:nso,nocc+1:nso),F_mi(nocc,nocc),F_me(nocc,nocc+1:nso)
    real*8  :: W_,W_mnij(nocc,nocc,nocc,nocc),W_mbej(nocc,nocc+1:nso,nocc,nocc+1:nso),&
    &       W_abef(nocc+1:nso,nocc+1:nso,nocc+1:nso,nocc+1:nso)
    real*8,parameter :: one=1.,one_half=0.5, one_forth=0.25
    ! may use one F(nso,nso) and one W(nso,nso,nso,nso) to store all intermediates
    F_ae = 0. ! intermediate F_ae from eq. 3 in Stanton et al. J. Chem. Phys. 94 (1990)
    do a = nocc+1, nso
      do e = nocc+1, nso
        F_ = (one - delta(a,e)) * Fock(a,e)
        do m = 1, nocc
          F_ = F_ - one_half * Fock(m,e) * ts(a,m)
          do f  = nocc+1, nso
            F_ = F_ + ts(f,m) * eri(m,a,f,e)
            do n = 1, nocc
              F_ = F_ - one_half * eri(m,n,e,f) &
              &    * ( td(a,f,m,n) + one_half * ( ts(a,m) * ts(f,n) - ts(a,n) * ts(f,m) ) )
            enddo
          enddo
        enddo
        F_ae(a,e) = F_
      enddo
    enddo
    F_mi = 0. ! intermediate F_mi from eq. 4 in Stanton et al. J. Chem. Phys. 94 (1990)
    do m = 1, nocc
      do i = 1, nocc
        F_ = (one - delta(m,i)) * Fock(m,i)
        do e = nocc+1, nso
          F_ = F_ + one_half * Fock(m,e) * ts(e,i)
          do n  = 1, nocc
            F_ = F_ + ts(e,n) * eri(m,n,i,e)
            do f = nocc+1, nso
              F_ = F_ + one_half * eri(m,n,e,f) &
              &    * ( td(e,f,i,n) + one_half * ( ts(e,i) * ts(f,n) - ts(e,n) * ts(f,i) ) )
            enddo
          enddo
        enddo
        F_mi(m,i) = F_
      enddo
    enddo
    F_me = 0. ! intermediate F_mi from eq. 5 in Stanton et al. J. Chem. Phys. 94 (1990)
    do m = 1, nocc
      do e = nocc+1, nso
        F_ = Fock(m,e)
        do n = 1, nocc
          do f  = nocc+1, nso
            F_ = F_ + ts(f,n) * eri(m,n,e,f)
          enddo
        enddo
        F_me(m,e) = F_
      enddo
    enddo
    W_mnij = 0. ! intermediate W_mnij from eq. 6 in Stanton et al. J. Chem. Phys. 94 (1990)
    do m = 1, nocc
      do n = 1, nocc
        do i = 1, nocc
          do j = 1, nocc
            W_ = eri(m,n,i,j)
            do e = nocc+1, nso
              W_ = W_ + ts(e,j) * eri(m,n,i,e) - ts(e,i) * eri(m,n,j,e)
              do f = nocc+1, nso
                W_ = W_ + one_forth * eri(m,n,e,f) &
                &    * ( td(e,f,i,j) + ts(e,i) * ts(f,j) - ts(e,j) * ts(f,i) )
              enddo
            enddo
            W_mnij(m,n,i,j) = W_
          enddo
        enddo
      enddo
    enddo
    W_abef = 0. ! intermediate W_abef from eq. 7 in Stanton et al. J. Chem. Phys. 94 (1990)
    do a = nocc+1, nso
      do b = nocc+1, nso
        do e = nocc+1, nso
          do f = nocc+1, nso
            W_ = eri(a,b,e,f)
            do m = 1, nocc
              W_ = W_ + ts(a,m) * eri(b,m,e,f) - ts(b,m) * eri(a,m,e,f)
              do n = 1, nocc
                W_ = W_ + one_forth * eri(m,n,e,f) &
                &         * ( td(a,b,m,n) + ts(a,n) * ts(b,m) - ts(a,m) * ts(b,n) )
              enddo
            enddo
            W_abef(a,b,e,f) = W_
          enddo
        enddo
      enddo
    enddo
    W_mbej = 0. ! intermediate W_mbej from eq. 8 in Stanton et al. J. Chem. Phys. 94 (1990)
    do m = 1, nocc
      do b = nocc+1, nso
        do e = nocc+1, nso
          do j = 1, nocc
            W_ = eri(m,b,e,j)
            do f = nocc+1, nso
              W_ = W_ + ts(f,j) * eri(m,b,e,f)
            enddo
            do n = 1, nocc
              W_ = W_ - ts(b,n) * eri(m,n,e,j)
              do f = nocc+1, nso
                W_ = W_ - one_half * td(f,b,j,n) + ts(f,j) * ts(b,n) *  eri(m,n,e,f) 
              enddo
            enddo
            W_mbej(m,b,e,j) = W_
          enddo
        enddo
      enddo
    enddo
  end subroutine constr_intm

  function kron(a,b)
    integer :: a,b
    real*8  :: kron
    if (mod(a,2)==mod(b,2)) then
      kron = 1.
    else
      kron = 0.
    endif
    return
  end function kron

  subroutine spinF(F,H,C,eri,nao,nso,nel)
    integer :: i,j,m,p,q,nao,nso,nel
    real*8,dimension(nao,nao) :: H,H_mo,C,tmp1,tmp2
    real*8  :: eri(nso,nso,nso,nso)
    real*8  :: F(nso,nso),kron
    H_mo=0.
!   do i = 1, nao
!     write(*,*) H(1,:)
!   enddo
    tmp2=transpose(C)
!   do i = 1, nao
!     write(*,*) tmp2(1,:)
!   enddo
    tmp1=matmul(tmp2,H)
    H_mo=matmul(tmp1,C)
    do i=1,nso
      p = ceiling(i/2.)
      do j=1,nso
        q = ceiling(j/2.)
        F(i,j) = H_mo(p,q) * kron(i,j)
        do m=1,nel
          F(i,j) = F(i,j) + eri(i,m,j,m)
        enddo
        write(3,*) i,p,j,q,kron(i,j),F(i,j),H_mo(p,q)
      enddo
    enddo
  end subroutine spinF

  subroutine spinmo(eri_so,eri_mo,nao,nso)
    integer :: nao,nso
    integer :: i,j,k,l,p,q,r,s
    real*8  :: eri_mo(nao,nao,nao,nao),eri_so(nso,nso,nso,nso),kron
    do i=1,nso
      p = ceiling(i/2.)
      do j=1,nso
        q = ceiling(j/2.)
        do k=1,nso
          r = ceiling(k/2.)
          do l=1,nso
            s = ceiling(l/2.)
            eri_so(i,j,k,l) = eri_mo(p,r,q,s)   * kron(p,r) * kron(q,s) &
            &                 - eri_mo(p,q,r,s) * kron(p,s) * kron(q,r)
          enddo
        enddo
      enddo
    enddo
  end subroutine spinmo

end module cc

function delta(i,j)
  integer :: i,j
  real*8  :: delta
  if (i.eq.j) then
    delta = 1.
  else
    delta = 0.
  end if
  return
end function delta
function kron(a,b)
  integer :: a,b
  real*8  :: kron
  if (mod(a,2)==mod(b,2)) then
    kron = 1.
  else
    kron = 0.
  endif
  return
end function kron


