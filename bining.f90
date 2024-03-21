program bining
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
    abar=0.0_r8
    do i=1,nf
        read(2,*) aux, a(i)
        abar=abar+a(i)
    enddo
    abar=abar/nf

    close(2)
    open(unit=3,file='bining.dat')
    ! make the block method
    do b=1,bmax
        ! we are making nb bins, then each bin has b points
        ab=0.0_r8
        nb=nf/b
        ! construct the nb blocks
        do ib=0,nb-1
            do jb=1,b
                ab(ib+1)=ab(ib+1)+a(jb+ib*b)
            enddo
        enddo
        ab=ab/b

        ! compute the standar deviation of this set of size
        sigmab(nb)=0.0_r8
        do ib=1,nb
        
            sigmab(nb)=sigmab(nb)+(ab(ib)-abar)**2.0_r8
        enddo
        sigmab(nb)=sigmab(nb)/(nb-1)
        write(3,*) b,sigmab(nb)
    enddo
    close(3)
end program bining