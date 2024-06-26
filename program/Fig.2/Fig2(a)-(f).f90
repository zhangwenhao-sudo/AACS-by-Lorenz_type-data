module Lorenz
    implicit none
    real,parameter :: h=0.001
    integer,parameter :: MaxT=2500000,T_trans=2400000,N=3,M=100
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
            x(i,3)=23.0*x3-1.0
            write(10,*) x(i,1),x(i,2),x(i,3)
        end do
    end subroutine x0

    subroutine fnf(xx,fx,i)
        implicit none
        real :: xx(N),fx(N),coupling
        integer :: t,i,j
        real :: a,b,c
        a=10.0
        b=100.0
        c=5.0
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
        real :: x(M,N),xx(N),fx(N),epsilon
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
        integer :: K !k:Neighbor number
        integer :: i,j,number
        allocate(neighbour_matrix(M,M))
        !initialize
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
                write(50,"(I2)",advance='no') neighbour_matrix(i,j)
            end do
            write(50,*)
        end do
        return
    end subroutine neighbour


end module Lorenz



program main
    use Lorenz
    implicit none
    integer :: t,i,ii,j,p=1
    real :: U_max(M)=-1000.0,U_min(M)=1000.0
    open(10,file="x0.txt")
    open(20,file="x_50.txt")
    open(50,file="neighbour_matrix.txt")
    open(60,file="i_t_x.txt")
    open(70,file="x_max_min.txt")
    call neighbour(2*p)
    do ii=100,100
        call x0()
        write(*,*) "Please enter the coupling strength: epsilon="
        read(*,*) epsilon
        do t=1,MaxT,1
            call rk4(x)
            if(t>=T_trans) then
                write(20,*) (t-T_trans)*h,x(50,1),x(50,2),x(50,3)
                call U_max_min(x(:,1),U_max,U_min)
                do j=1,M,1
                    if(mod(t,100)==0) then
                        write(60,*) j,(t-T_trans)*h,x(j,1),x(j,2),x(j,3)
                    end if
                end do
            end if
        end do
        do i=1,M,1
            write(70,*) i,U_max(i),U_min(i)
        end do
    end do
    close(10)
    close(20)
    close(50)
    close(60)
    close(70)
    deallocate(neighbour_matrix)
end program main