program background_Tools
    implicit none
    integer :: i,j
    real :: data(4001,2)
    OPEN(UNIT=10, FILE='c_state.txt', STATUS='OLD')
    OPEN(UNIT=20, FILE='c_range_state.txt', STATUS='new')
    data=0.0
    do i=1,4001,1
        read(10,*) data(i,1),data(i,2)
        do j=-30,15,1
            write(20,*) data(i,1),j,data(i,2)
        end do
    end do
end program
