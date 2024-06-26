module Lorenz
    implicit none
    real,parameter :: h=0.001,PI=3.1415926
    integer,parameter :: MaxT=2500000,T_trans=2000000,N=3,M=100,number_cycle=200
    real :: x(M,N),R,omega(M),epsilon,c=5.0
    integer :: number_cycle_state(9),gradient_number_omega
    integer,allocatable :: neighbour_matrix(:,:)
contains

    subroutine x0()
        implicit none
        integer :: i
        real :: x1,x2,x3
        call random_seed()
        open(5,file="x_y_z.txt")
        do i=1,M,1
            call random_number(x1)
            call random_number(x2)
            call random_number(x3)
            x(i,1)=31.0*x1-14.0
            x(i,2)=44.0*x2-20.0
            x(i,3)=23.0*x3-1.0
            !            write(5,*) x(i,1),x(i,2),x(i,3)
        end do
        close(5)
    end subroutine x0

    subroutine fnf(xx,fx,i)
        implicit none
        real :: xx(N),fx(N),coupling
        integer :: t,i,j
        real :: a,b
        a=10.0
        b=100.0
        coupling=0.0
        do j=1,M,1
            coupling=coupling+neighbour_matrix(i,j)*(x(j,1)-x(i,1))
        end do
        fx(1)=a*(xx(2)-xx(1))+epsilon*coupling
        fx(2)=-xx(1)*xx(3)-c*xx(2)
        fx(3)=xx(1)*xx(2)-b
        return
    end subroutine fnf


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
    !Umax1 andUmin1(10s)
    subroutine U_max_min_1(xxx1,U_max1,U_min1)
        implicit none
        real :: U_max1(M),U_min1(M),xxx1(M)
        integer :: i
        do i=1,M,1
            if(U_max1(i)<xxx1(i)) then
                U_max1(i)=xxx1(i)
            end if
            if(U_min1(i)>xxx1(i)) then
                U_min1(i)=xxx1(i)
            end if
        end do
        return
    end subroutine U_max_min_1

    !omega
    subroutine Order_Parameter_Omega(sequence)
        implicit none
        integer,parameter :: cycle_Number=10   !Number of records
        real :: sequence(M,MaxT-T_trans)
        integer :: i,j,jj,periodT(M),number(M),t_start(M),t_end(M),counter
        number=0
        t_start=0
        t_end=0
        periodT=0
        do jj=1,M,1
            counter=1
            do j=1,MaxT-T_trans-1,1
                if(sequence(jj,j)>0.0.and.sequence(jj,j+1)<0.0) then
                    if(counter==1) then
                        t_start(jj)=j
                    else if(counter==10) then
                        t_end(jj)=j
                        exit
                    end if
                    counter=counter+1
                end if
            end do
            if(t_start(jj)==0.and.t_end(jj)==0) then
                periodT(jj)=0
            else
                periodT(jj)=t_end(jj)-t_start(jj)
            end if
        end do
        do i=1,M,1
            if(periodT(i)==0) then
                omega(i)=0
            else
                omega(i)=(2.0*PI*cycle_Number)/(periodT(i)*h)
            end if
            write(50,*) epsilon,i,omega(i)
        end do
        return
    end subroutine Order_Parameter_Omega

    !Synchronization coefficient
    subroutine Order_Parameter_R(sequence)
        implicit none
        integer :: i,j,t,ii
        real :: ansRsub,x_average(MaxT-T_trans),sequence(M,MaxT-T_trans)
        ansRsub=0.0
        x_average=0.0
        do t=1,MaxT-T_trans,1
            x_average(t)=sum(sequence(:,t))/real(M)
        end do
        do i=1,M,1
            ansRsub=ansRsub+sum(sequence(i,:)**2)/real(MaxT-T_trans)-(sum(sequence(i,:))/&
                    real(MaxT-T_trans))**2
        end do
        if(ansRsub<1.0e-5) then
            R=1.0
        else
            R=(sum(x_average(:)**2)/real(MaxT-T_trans)-(sum(x_average(:))/real(MaxT-T_trans))**2)/&
                    (ansRsub/real(M))
        end if
        return
    end subroutine Order_Parameter_R

    !ladder
    subroutine gradient(array1)
        implicit none
        real :: gradient_array(M),array1(M)
        integer :: count,i
        count=0
        do i=1,M,1
            if(all(abs(array1(i) - gradient_array(:count))>0.002)) then
                count = count +1
                gradient_array(count)=array1(i)
            end if
        end do
        gradient_number_omega=count
        return
    end subroutine gradient

    subroutine Order_Parameter_Multiple_states(U_max,U_min,U_max1,U_min1)
        implicit none
        integer :: i,j,total_Rest,Traveling_number,state(9),num_small_counter
        real :: U_max(M),U_min(M),U_max1(M),U_min1(M),state_epsilon(3)
        total_Rest=0
        Traveling_number=0
        num_small_counter=0
        state=0
        state_epsilon=0
        do i=1,M,1
            if(abs(omega(i))<=1.0e-5) then
                total_Rest=total_Rest+1
            else if(abs(omega(i)-omega(1))<=1.0e-3.and.abs(omega(i))>1.0e-3.and.R<0.8) then
                Traveling_number=Traveling_number+1
            end if
        end do
        do i=1,M,1
          if((abs(U_max(i))-0.0<5.0).or.(abs(U_min(i))-0.0<5.0)) then
                num_small_counter=num_small_counter+1
          end if
        end do
        call gradient(omega)
        !        write(60,*) epsilon,total_Rest/real(M)
        if(total_Rest==M.and.all(abs(U_max(:)-10.0)<0.1).and.all(abs(U_min(:)-10.0)<0.1)) then
            !Resting state 1
            state(1)=1
            number_cycle_state(1)=number_cycle_state(1)+1
        else if(total_Rest==M.and.all(abs(U_max(:)+10.0)<0.1).and.all(abs(U_min(:)+10.0)<0.1)) then
            !Resting state 2
            state(2)=1
            number_cycle_state(2)=number_cycle_state(2)+1
        else if(Traveling_number==M) then
            !Traveling
            state(3)=1
            number_cycle_state(3)=number_cycle_state(3)+1
            state_epsilon(2)=1
            !write(100,*) c,epsilon,state_epsilon(1),state_epsilon(2),state_epsilon(3)
            !stop
        else if(R>0.95) then
            !synchronization
            state(4)=1
            number_cycle_state(4)=number_cycle_state(4)+1
        else if(num_small_counter==M) then
            !Small local oscillation
            state(8)=1
            number_cycle_state(8)=number_cycle_state(8)+1
            state_epsilon(1)=1
            write(100,*) c,epsilon,state_epsilon(1),state_epsilon(2),state_epsilon(3)
            stop
        else if(total_Rest/=M.and.all(abs(U_max1(:)-U_min1(:))<0.1).and.all(abs(U_max(:)-U_min(:))>1.0)) then
            state(6)=1
            number_cycle_state(6)=number_cycle_state(6)+1
        else if(gradient_number_omega==2) then
            !Double frequency
            state(7)=1
            number_cycle_state(7)=number_cycle_state(7)+1
        else if(0<total_Rest.and.total_Rest<M) then
            !blend
            state(5)=1
            number_cycle_state(5)=number_cycle_state(5)+1
            state_epsilon(3)=1
            write(110,*) c,epsilon,state_epsilon(1),state_epsilon(2),state_epsilon(3)
        else
            !other
            state(9)=1
            number_cycle_state(9)=number_cycle_state(9)+1
        end if
        write(70,*) c,epsilon,state(1),state(2),state(3),state(4),state(5),state(6),state(7),state(8),state(9)
        return
    end subroutine Order_Parameter_Multiple_states

    !Locally coupled one-dimensional ring
    subroutine neighbour(K)
        implicit none
        integer :: K   !k:Neighbor number
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
    integer :: t,i,ii,j,jj,p=1
    real :: U_max(M),U_min(M),sequence(M,MaxT-T_trans),U_max1(M),U_min1(M)
    open(10,file="neighbour_matrix.txt")
    open(20,file="i_t_x_1.txt")
    open(30,file="xmax_xmin.txt")
    open(40,file="i_t_x.txt")
    open(50,file="i_omega.txt")
    open(60,file="epsilon_total_Rest_P.txt")
    open(70,file="epsilon_state.txt")
    open(80,file="epsilon_R.txt")
    open(90,file="epsilon_number_cycle_state_P.txt")
    open(100,file="c_epsilon_stateEpsilon.txt")
    open(110,file="c_epsilon_stateEpsilon_blend.txt")
    call neighbour(2*p)
    do epsilon=1, 50, 1

        !Cycle number
        number_cycle_state=0
        do j=1,number_cycle,1
            write(*,*) "epsilon=",epsilon,"Cycle number=",j
            U_max=-10000.0
            U_min=10000.0
            U_max1=-10000.0
            U_min1=10000.0
            sequence=0.0
            omega=0.0
            R=0.0
            gradient_number_omega=0
            call x0()
            do t=1,MaxT,1
                call rk4(x)
                if(t>T_trans) then
                    sequence(:,t-T_trans)=x(:,1)
                    !write(20,*) (t-T_trans)*h,x(1,1),x(1,2),x(1,3)
                    call U_max_min(x(:,1),U_max,U_min)
                    if(t>MaxT-10000) then
                        call U_max_min_1(x(:,1),U_max1,U_min1)
                    end if
                end if
            end do
            call Order_Parameter_Omega(sequence)
            call Order_Parameter_R(sequence)
            write(80,*) c,R
            call Order_Parameter_Multiple_states(U_max,U_min,U_max1,U_min1)
            do ii=1,M,1
                write(30,*)  ii,U_max(ii),U_min(ii)
            end do
        end do
        write(90,"(11(F15.6))") c,epsilon,number_cycle_state(1)/real(number_cycle),number_cycle_state(2)/real(number_cycle),&
                number_cycle_state(3)/real(number_cycle),number_cycle_state(4)/real(number_cycle),&
                number_cycle_state(5)/real(number_cycle),number_cycle_state(6)/real(number_cycle),&
                number_cycle_state(7)/real(number_cycle),number_cycle_state(8)/real(number_cycle),&
                number_cycle_state(9)/real(number_cycle)
    end do
    close(10)
    close(20)
    close(30)
    close(40)
    close(50)
    close(60)
    close(70)
    close(80)
    close(90)
    deallocate(neighbour_matrix)
end program main