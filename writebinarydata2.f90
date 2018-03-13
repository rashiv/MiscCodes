program binarydata
   integer :: i, var1
   real :: var2
   open(unit = 10, status='replace',file='test.bin',form='unformatted')
   do i = 1, 10
      write(10) i*i
   end do
   do i = 1, 5
      write(10) i*33.0
   end do
   close(10)
   open(unit = 10, status='old',file='test.bin',form='unformatted')
   do i = 1, 10
      read(10) var1
      print *, var1
   end do
   do i = 1, 5
      read(10) var2
      print *, var2
   end do
end program binarydata