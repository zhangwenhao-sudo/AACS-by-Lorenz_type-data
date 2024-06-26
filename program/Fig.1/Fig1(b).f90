!------------Calculate Lyaponov Exponent Spectrum of Lorenze_like System----------------------
!-----------------------Lorenze Equation:-------------------
!------------------------dx/dt=aa*(y-x)----------------------
!------------------------dy/dt=cc*x-x*z-y--------------------
!------------------------dz/dt=x*y-bb*z----------------------
module head
    implicit none
    real(kind=8), parameter :: aa=10.0D0, bb=100.0D0, dt=0.001D0
    real(kind=8) :: cc
    integer(kind=4), parameter :: st=1e5,n=1e6
    integer(kind=4), parameter :: m=3
contains

    subroutine initial(in)
        implicit none
        real(kind=8), dimension(m) :: in
        call random_seed()
        call random_number(in)
    end subroutine

    subroutine fnf1(in,out)
        implicit none
        real(kind=8), dimension(m) :: in, out
        out(1)=dt*(aa*(in(2)-in(1)))
        out(2)=dt*(-cc*in(2)-in(1)*in(3))
        out(3)=dt*(in(1)*in(2)-bb)
    end subroutine

    subroutine fnf2(in,out,jj)
        implicit none
        real(kind=8), dimension(m,m) :: in, out
        real(kind=8), dimension(m,m) :: jj
        integer(kind=4) :: i
        out=dt*matmul(jj,in) 
    end subroutine

end module

program main
    use head
    implicit none
    integer(kind=4) :: step,i,j,ii
    real(kind=8), dimension(m) :: fnorm, flya, fnorm1,fnorm2
    real(kind=8), dimension(m) :: y,k1,k2,k3,k4
    real(kind=8), dimension(m,m) :: jj,w,v,k11,k22,k33,k44
    real(kind=8), external ::funcdj
    real(kind=8), external ::funcfs

    open(unit=10,file='flya.dat')
    open(unit=11,file='x-y.dat')

    fnorm=0.0D0
    flya=0.0D0

    do i=1,m,1
        do j=1,m,1
            if(i==j) then
                w(i,j)=1.0D0
            else
                w(i,j)=0.0D0
            end if
        end do
    end do

    v=0.0D0

    do ii=0,4000,1
        call initial(y)
        cc=real(ii*0.01-20.0)
        do step=1,n,1

            !-----------Evolution of Single Lorenz_like system-------------------------
            call fnf1(y,k1)
            call fnf1(y+K1/2,K2)
            call fnf1(y+K2/2,K3)
            call fnf1(y+K3,K4)
            y=y+(k1+2.0D0*k2+2.0D0*k3+k4)/6.0d0

            !---------Jacobian Matrix of Single Lorenz_like System-------------------
            jj(1,1)=-aa
            jj(1,2)=aa
            jj(1,3)=0.0D0
            jj(2,1)=-y(3)
            jj(2,2)=-cc
            jj(2,3)=-y(1)
            jj(3,1)=y(2)
            jj(3,2)=y(1)
            jj(3,3)=0.0D0
            !-------------------------------------------------------------------

            !-----------Evolution of Jacobian Matrix----------------------------
            call fnf2(w,k11,jj)
            call fnf2(w+K11/2,K22,jj)
            call fnf2(w+K22/2,K33,jj)
            call fnf2(w+K33,K44,jj)
            w=w+(k11+2.0D0*k22+2.0D0*k33+k44)/6.0d0
            !---------------------------------------------------------------------

            !-----------Gram-Schmit Normalization-----------------------------------
            do j=2,m,1
                do i=1,j-1,1
                    w(:,j)=w(:,j)-funcdj(w(:,j),w(:,i))/funcfs(w(:,i))*w(:,i)
                end do
            end do
            !-------------------------------------------------------------------------

            !-----------Calculate the module of every vector direction----------------
            do j=1,m,1
                fnorm(j)=funcfs(w(:,j))
                fnorm(j)=sqrt(fnorm(j))
                w(:,j)=w(:,j)/fnorm(j) !---------Renormalization every vector
                fnorm(j)=log(fnorm(j))
            end do
            !--------------------------------------------------------------------------

            !----------Calculate Lyaponov Exponent of every vector direction------------
            if(step>st) then
                do i=1,m,1
                    flya(i)=flya(i)+fnorm(i)
                end do
            end if
            !----------------------------------------------------------------------------

            if(step>n-st) then
                write(11,'(f15.9,2x,f15.9,2x,f15.9)') y(1),y(2),y(3)
            end if
        end do


        flya=flya/(n-st)/dt

        write(10,'(f15.9,2x,f15.9,2x,f15.9,2x,f15.9)') cc,flya(1),flya(2),flya(3)
    end do
    stop
end program

!--------------------calculate dot metrix of two vector: a.b-------------------
real(kind=8) function funcdj(u,v)
    use head
    implicit none
    real(kind=8), dimension(m) :: u, v
    integer(kind=4) :: i

    funcdj=0.0d0
    do i=1,m,1
        funcdj=funcdj+u(i)*v(i)
    end do

    return
end function
!--------------------------------------------------------------------------------

!------------------calculate module of one vector direction---------------------
real(kind=8) function funcfs(u)
    use head
    implicit none
    real(kind=8), dimension(m) :: u
    integer(kind=4) :: i

    funcfs=0.0d0
    do i=1,m,1
        funcfs=funcfs+u(i)*u(i)
    end do

    return
end function
!-------------------------------------------------------------------------------