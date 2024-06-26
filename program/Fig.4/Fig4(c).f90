module Lorenz
    implicit none
    real,parameter :: h=0.001
    integer,parameter :: MaxT=2500000,T_trans=2000000,N=3,M=100,cycle_number=200
    real :: x(M,N),epsilon
    integer,allocatable :: neighbour_matrix(:,:)
contains

    subroutine x0()
        implicit none
        integer :: i
        real :: x1,x2,x3
        call random_seed()
        do i=1,M,1
            call random_number(x1)
            call random_number(x2)
            call random_number(x3)
            x(i,1)=31.0*x1-14.0
            x(i,2)=44.0*x2-20.0
            x(i,3)=23.0*x3-1
            !write(5,*) x(i,1),x(i,2),x(i,3)
        end do
    end subroutine x0

    subroutine fnf(xx,fx,i)
        implicit none
        real :: xx(N),fx(N),coupling
        integer :: t,i,j
        real :: a,b,c
        a=10.0
        b=100.0
        c=-5.0
        coupling=0.0
        do j=1,M,1
            coupling=coupling+neighbour_matrix(i,j)*(x(j,1)-x(i,1))
        end do
        fx(1)=a*(xx(2)-xx(1))+epsilon*coupling
        fx(2)=-xx(1)*xx(3)-c*xx(2)
        fx(3)=xx(1)*xx(2)-b
        return
    end subroutine fnf

    !Runge Kuta of the fourth order
    subroutine rk4(x)
        implicit none
        integer :: i
        real :: x(M,N),xx(N),fx(N)
        real :: kx1(N),kx2(N),kx3(N),kx4(N)
        do i=1,M,1
            xx=x(i,:)
            call fnf(xx,fx,i)
            kx1=h*fx
            xx=x(i,:)+0.5*kx1
            call fnf(xx,fx,i)
            kx2=h*fx
            xx=x(i,:)+0.5*kx2
            call fnf(xx,fx,i)
            kx3=h*fx
            xx=x(i,:)+kx3
            call fnf(xx,fx,i)
            kx4=h*fx
            x(i,:)=x(i,:)+(kx1+2.0*kx2+2.0*kx3+kx4)/6.0
        end do
        return
    end subroutine rk4

    !Umax and Umin
    subroutine U_max_min(xxx,U_max,U_min)
        implicit none
        real :: U_max(M),U_min(M),xxx(M)
        integer :: i
        do i=1,M,1
            if(U_max(i)<xxx(i)) then
                U_max(i)=xxx(i)
            end if
            if(U_min(i)>xxx(i)) then
                U_min(i)=xxx(i)
            end if
        end do
        return
    end subroutine U_max_min

    !Locally coupled one-dimensional ring
    subroutine neighbour(K)
        implicit none
        integer :: K
        integer :: i,j,number
        allocate(neighbour_matrix(M,M))
        neighbour_matrix=0
        do i=1,M-1,1
            do j=i+1,M,1
                if(abs(i-j)<=K/2.or.abs(i-j)>=M-K/2) then
                    neighbour_matrix(i,j)=1
                    neighbour_matrix(j,i)=1
                end if
            end do
        end do
        !Output network matrix
        do i=1,M,1
            do j=1,M,1
                write(10,"(I2)",advance='no') neighbour_matrix(i,j)
            end do
            write(10,*)
        end do
        return
    end subroutine neighbour


end module Lorenz



program main
    use Lorenz
    implicit none
    integer :: t,i,ii,j,counter,p=1,num_large_counter,num_small_counter,num_large(cycle_number),num_small(cycle_number),num_blend(cycle_number)
    real :: U_max(M),U_min(M),num_large_average,num_small_average,num_blend_average
    open(5,file="x_y_z.txt")
    open(10,file="neighbour_matrix.txt")
    open(20,file="epsilon_large_small_blend_p.txt")
    open(30,file="xmax_xmin.txt")
    call neighbour(2*p)
    do epsilon=0.5, 30, 0.5
        num_large=0
        num_small=0
        num_blend=0
        num_large_average=0.0
        num_small_average=0.0
        num_blend_average=0.0    
        do ii=1,cycle_number,1
            print *, epsilon, ii
            U_max=-10000.0
            U_min=10000.0
            num_large_counter=0
            num_small_counter=0
            call x0()
            do t=1,MaxT,1
                call rk4(x)
                if(t>T_trans) then
                    call U_max_min(x(:,1),U_max,U_min)
                end if
            end do
            do i=1,M,1
                write(30,*) epsilon,ii,i,U_max(i),U_min(i)
                if((abs(U_max(i))-0.0>5.0).and.(abs(U_min(i))-0.0>5.0)) then
                    num_large_counter=num_large_counter+1
                elseif((abs(U_max(i))-0.0<5.0).or.(abs(U_min(i))-0.0<5.0)) then
                    num_small_counter=num_small_counter+1
                end if
            end do
            if(num_large_counter==M) then
                num_large(ii)=num_large(ii)+1
            else if(num_small_counter==M) then
                num_small(ii)=num_small(ii)+1
            else if((num_large_counter>0.and.num_large_counter<M).and.(num_small_counter>0.and.num_small_counter<M)) then
                num_blend(ii)=num_blend(ii)+1
            end if
        end do
        do j=1,cycle_number,1
            num_large_average=num_large_average+num_large(j)
            num_small_average=num_small_average+num_small(j)
            num_blend_average=num_blend_average+num_blend(j)
        end do
        write(20,"(5(F15.6))") epsilon,num_large_average/real(cycle_number),&
        &num_small_average/real(cycle_number),num_blend_average/real(cycle_number)
    end do
    close(5)
    close(10)
    close(20)
    close(30)
    deallocate(neighbour_matrix)
end program main