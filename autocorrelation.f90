program autocorrelation
    implicit none

    ! defining our types
    integer, parameter :: r8=selected_real_kind(15,8)
    integer, parameter :: i4=selected_int_kind(8)

    ! defining ou arrays
    real(kind=r8), allocatable, dimension(:) :: x, chi
    integer(kind=i4), allocatable, dimension(:) :: y

    ! defining ou variables
    real(kind=r8) :: chi0, xbar, mu, sum, mu1, mu2, tau, a, b
    integer(kind=i4) :: nfile, ttrans, i, j, n, t, tmax
    character(len=20) :: filename

    ! reading the input values
    open(unit=1,file='autocorrelation-in.dat')
    read(1,*) filename ! reading the input file name
    read(1,*) ttrans                   ! tranzient time
    read(1,*) tmax                       ! maximum time to be used in the correlation function
    close(1)
    
    ! now we read the shit that it contains
    open(unit=1,file=filename)
    nfile=filesize(1)
    close(1)
    
    open(unit=2,file=filename)
    ! allocatting the needed memory
    n=nfile-ttrans
    
    allocate(x(1:n),y(1:n),chi(1:n))
    ! reading the data and calculating the mean value
    mu=0.0_r8
    j=1
    !print*,n, nfile, ttrans
    do i=1,nfile
        read(2,*) a, b
        if(i.gt.ttrans)then
            y(j)=a
            x(j)=b
            mu=mu+x(j)
            j=j+1
        endif
    enddo
    close(2)
    mu=mu/(j-1)

    ! now we calculate the autocorrelation function
    chi=0.0_r8
    do t=1,n!tmax
    !print*,t
        chi0=0.0_r8
        mu1=0.0_r8
        mu2=0.0_r8
        do i=1,n-t!tmax-t
            chi0=chi0+x(i)*x(i+t)
            mu1=mu1+x(i)
            mu2=mu2+x(i+t)
        enddo
        chi(t)=chi0/(n-t)-mu1*mu2/(n-t)**2.0_r8
    enddo
    chi=chi/chi(1)
    !print*,chi  

    ! get out with the data
    open(unit=3,file='chi.dat')
    tau=0.50_r8
    do t=1,tmax
        write(3,*) t-1, chi(t), log(chi(t))
        tau=tau+chi(t)
    enddo
    open(unit=4,file='tauint.dat')
    write(4,*) tau  
    close(3)
    contains
    ! function that return us the size of the input file
    function filesize(fileunit)
        real(kind=r8) :: input
        integer(kind=i4) :: filesize, fileunit, i, io
        i=-1
        io=0
        do while (io.eq.0)
            read(fileunit,*,iostat=io) input
            i=i+1
        enddo
        filesize=i
    endfunction filesize
end program autocorrelation