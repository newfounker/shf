module MP
implicit none
contains
  subroutine mp_two(mp2en,eri_ao,C_mo,en_mo,nao,nocc)
    integer :: nao,nocc
    real*8,dimension(nao,nao,nao,nao) :: eri_ao,eri_mo
    real*8  :: C_mo(nao,nao),en_mo(nao),mp2en
    mp2en=0.

    call transf_eri(eri_ao,eri_mo,C_mo,nao)
    call mp_energy(mp2en,eri_mo,en_mo,nao,nocc)
  end subroutine mp_two

  subroutine mp_energy(mp2en,eri,en,nao,nocc)
    integer :: nao,nocc
    real*8  :: eri(nao,nao,nao,nao),en(nao),mp2en
    integer :: i,j,k,l
    mp2en = 0.
    do i = 1, nocc
      do k = nocc+1, nao
        do j = 1, nocc
          do l = nocc+1, nao
            mp2en = mp2en + eri(i,k,j,l) * ( 2*eri(i,k,j,l) - eri(i,l,j,k) ) &
            &               / ( en(i) + en(j) - en(k) - en(l) )
          enddo
        enddo
      enddo
    enddo
  end subroutine mp_energy

  subroutine transf_eri(eri_ao,eri_mo,C,nao)
    integer :: nao
    real*8,dimension(nao,nao,nao,nao) :: eri_ao,eri_mo,temp
    real*8,dimension(nao,nao)         :: C,tmp1,tmp2
    integer :: i,j,k,l
    eri_mo=0.
    tmp1=0.
    tmp2=0.
    temp=0.
    do i = 1, nao
      do j = 1, i!-1
        do k = 1, nao
          do l = 1, k!-1
            tmp1(k,l) = eri_ao(i,j,k,l)
            tmp1(l,k) = tmp1(k,l)
          enddo
        enddo
        tmp2=0.
        tmp2 = matmul(transpose(C),tmp1)
        tmp1=0.
        tmp1 = matmul(tmp2,C)
        do k = 1, nao
          do l = 1, k!-1
            temp(k,l,i,j) = tmp1(k,l)
          enddo
        enddo
      enddo
    enddo
    do k = 1, nao
      do l = 1, k!-1
        tmp1=0.
        tmp2=0.
        do i = 1, nao
          do j = 1, i!-1
            tmp1(i,j) = temp(k,l,i,j)
            tmp1(j,i) = tmp1(i,j)
          enddo
        enddo
        tmp2=0.
        tmp2 = matmul(transpose(C),tmp1)
        tmp1=0.
        tmp1 = matmul(tmp2,C)
        do i = 1, nao
          do j = 1, i!-1
            eri_mo(k,l,i,j) = tmp1(i,j)
            eri_mo(k,l,j,i) = tmp1(i,j)
            eri_mo(l,k,i,j) = tmp1(i,j)
            eri_mo(l,k,j,i) = tmp1(i,j)
            eri_mo(i,j,k,l) = tmp1(i,j)
            eri_mo(i,j,l,k) = tmp1(i,j)
            eri_mo(j,i,k,l) = tmp1(i,j)
            eri_mo(j,i,l,k) = tmp1(i,j)
          enddo
        enddo
      enddo
    enddo
  end subroutine transf_eri
end module MP
