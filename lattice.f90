! program that solves the 2-dimentional ising model
! we will calculate only the behavior with temperature for now
program ising2d
    implicit none

    ! defining our types
    integer, parameter :: r8=selected_real_kind(15,8)
    integer, parameter :: i4=selected_int_kind(8)
    
    ! defining our vectors
    real(kind=r8), allocatable, dimension(:,:) :: spins

    ! defining our variables
    real(kind=r8) :: epsilon, temp, r
    integer(kind=i4) :: nspins, nmcsteps, nterm, ncorr, i, j

    ! reading our input values from the input file
    open(unit=1,file='lattice-in.dat')
    read(1,*) nspins, nmcsteps
    read(1,*) temp
    read(1,*) nterm
    read(1,*) ncorr
    close(1)

    ! allocating the needed memory and initiallyzing the vectors
    allocate(spins(0:nspins-1,0:nspins-1))
    spins=1.0_r8
    !do i=0,nspins-1
    !do j=0,nspins-1
    !    call random_number(r)
    !    if (r.lt.0.50_r8) then
    !        spins(i,j)=1
    !    else
    !        spins(i,j)=-1
    !    endif
    !enddo
    !enddo
    
    call mclattice(spins,temp,nmcsteps,nspins)
    ! deallocating the memory
    deallocate(spins)
    contains

    ! function that return us the variantion of the energy when we flip the j spin
    function deltaene(i,j,spins,delta)
        real(kind=r8) :: spins(0:nspins-1,0:nspins-1), deltaene
        integer(kind=i4) :: i, j, ip1, im1, jp1, jm1, delta
        ip1=i+1
        im1=i-1
        jp1=j+1
        jm1=j-1
        ! appling the periodic boundary conditions
        if(i.eq.0) im1=nspins-1
        if(i.eq.nspins-1) ip1=0
        if(j.eq.0) jm1=nspins-1
        if(j.eq.nspins-1) jp1=0
        deltaene=(spins(ip1,j)+spins(im1,j)+spins(i,jp1)+spins(i,jm1))*(delta*spins(i,j))
    endfunction deltaene

    ! subroutine that do the monte carlo method in one time step on the lattice
    subroutine mclattice(spins,temp,nmcsteps,nspins)
        real(kind=r8), intent(inout) :: spins(0:nspins-1,0:nspins-1)
        real(kind=r8), intent(in) :: temp
        real(kind=r8) :: denes, energy, magnetization
        integer(kind=i4), intent(in) :: nmcsteps, nspins
        integer(kind=i4) :: i, j, k, l
        
        ! initial values
        open(unit=2,file='lattice-mag.dat')
        open(unit=3,file='lattice-ene.dat')
        
        ! we do k passages though the lattice
        do k=1,nmcsteps+nterm
            do l=1,ncorr ! discart correlated configurations
                energy=0.0_r8
                magnetization=0.0_r8
                do i=0,nspins-1
                do j=0,nspins-1
                    ! for each spin we calculate the variantion on the erngy for flipping it
                    denes=0.0_r8
                    denes=deltaene(i,j,spins,2)

                    ! if this variation is negative, we accept the change
                    ! if not, we apply the metropolis criteria
                    call random_number(epsilon)
                    if((denes.lt.0._r8)) then
                        spins(i,j)=-spins(i,j)
                    endif
                    if((denes.gt..0_r8).and.((exp(-denes/temp).gt.epsilon)))then
                        spins(i,j)=-spins(i,j)
                    endif

                    magnetization=magnetization+spins(i,j)
                    energy=energy+0.5_r8*deltaene(i,j,spins,-1)
                enddo
                enddo
            enddo
            if(k.gt.nterm)then
                write(2,*) k, (magnetization/nspins**2.0_r8)
                write(3,*) k, energy/nspins**2.0_r8
            endif
        enddo
        close(2)
        close(3)
    endsubroutine mclattice
endprogram ising2d
