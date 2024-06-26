!calculate phase diagram
module Lorenz
    implicit none
    real,parameter :: h=0.001
    integer,parameter :: MaxT=2500000,T_trans=2400000,N=3,M=1
    real :: x(M,N),c
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
            x(i,1)=24.0*x1-12.0
            x(i,2)=30.0*x2-15.0
            x(i,3)=20.0*x3-1.0
            write(10,*) x(i,1),x(i,2),x(i,3)
        end do
    end subroutine x0

    subroutine fnf(xx,fx,i)
        implicit none
        real :: xx(N),fx(N)
        integer :: t,i,j
        real :: a,b
        a=10.0
        b=100.0
        fx(1)=a*(xx(2)-xx(1))
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


end module Lorenz



program main
    use Lorenz
    implicit none
    integer :: t,i,j
    real :: U_max(M),U_min(M)
    open(10,file="x0_y0_z0.txt")
    open(20,file="c_x_max_min.txt")
    open(30,file="i_t_x.txt")
    open(40,file="t_x_1.txt")
    call x0()
    U_max=-1000.0
    U_min=1000.0
    write(*,*) "Please enter the value of c="
    read(*,*) c
    do t=1,MaxT,1
        call rk4(x)
        if(t>T_trans) then
            call U_max_min(x(:,1),U_max,U_min)
            !phase diagram
            write(40,*) (t-T_trans)*h,x(1,1),x(1,2),x(1,3)
        end if
    end do
    do i=1,M,1
        write(20,*) c,i,U_max(i),U_min(i)
    end do
    close(10)
    close(20)
end program main