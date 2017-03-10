program new_point

implicit none

integer :: i,j
integer(8) :: temp

temp = 3600000000_8
!!!here I create bunch of points which have constant interval and is bounded insided a circle
open(11,file='new_point.txt',action='write')
do i = -60000,60000,200
   do j = -60000,60000,200
      if(i**2 + j**2 <= temp) write(11,*) i,j
   enddo
enddo
close(11)

end program new_point

