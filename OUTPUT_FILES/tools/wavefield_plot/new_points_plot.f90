program new_point

implicit none

integer :: i,j
!integer(8) :: temp

!temp = 3600000000_8
!!!here I create bunch of points which have constant interval and is bounded insided a circle
open(11,file='new_point.txt',action='write')
do i = -10000,10000,100
   do j = -10000,10000,100
      write(11,*) i,j
   enddo
enddo
close(11)

end program new_point
