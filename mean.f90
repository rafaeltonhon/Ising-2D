program mean
    implicit none

    ! defining our types
    integer, parameter :: r8=selected_real_kind(15,8)
    integer, parameter :: i4=selected_int_kind(8)

    ! defining ou arrays
    real(kind=r8), allocatable, dimension(:) :: a, sigmab, ab
    integer(kind=i4), allocatable, dimension(:) :: y, yk

    ! defining ou variables
    real(kind=r8) :: abar, aux
    integer(kind=i4) :: nf, i, ib, jb, nb, b, bmax
    integer(kind=i4) :: ttrans, tmax
    character(len=20) :: filename

    ! reading the input values
    open(unit=1,file='autocorrelation-in.dat')
    read(1,*) filename ! reading the input file name
    read(1,*) ttrans                   ! tranzient time
    read(1,*) tmax                       ! maximum time to be used in the correlation function
    read(1,*) nf                       ! size of the file
    read(1,*) bmax                         ! maximum number of bins
    close(1)
    
    allocate(a(1:nf),ab(1:nf),sigmab(1:nf))
    ! read the input values
    open(unit=2,file=filename)
    do i=1,nf
        read(2,*) aux, a(i)
    enddo
    open(unit=2,file='mean.dat',position='append')
    write(2,*) meana(a,nf), sigmaa(a,meana(a,nf),nf)
    close(2)
contains
function meana(a,np)
    real(kind=r8) :: a(np), mia, meana
    integer(kind=i4) :: np, i
    mia=0.0_r8
    do i=1,np
        mia=mia+a(i)
    enddo
    meana=mia/np
endfunction meana

function sigmaa(a,mia,np)
    real(kind=r8) :: a(np), mia, siga, sigmaa
    integer(kind=i4) :: np, i
    sigmaa=0.0_r8
    do i=1,np
        siga=siga+(mia-a(i))**2.0_r8
    enddo
    sigmaa=dsqrt(siga/(np*(np-1)))
endfunction sigmaa
endprogram mean