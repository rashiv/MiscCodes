program testsort

integer:: n,dist(9)
integer, allocatable :: sdist(:)
n=9
data dist /10,20,60,40,50,30,30,10,20/

allocate(sdist(9))

print *,dist
call sort(dist,n,sdist)

print *, dist
print *, sdist

contains

subroutine sort(a,length,b)
!subroutine to sort a 1-D array in descending order
integer :: length,i,j,a(length),b(length),temp

do i=1,length
  b(i) = i
end do

do i = 2, length
  j = i-1
  temp = a(i)
  do while (j>=1.and.a(j)<temp)
    a(j+1) = a(j)
    b(j+1) = b(j)
    j = j-1 
  end do
  a(j+1) = temp
  b(j+1) = i
end do

end subroutine sort
end program
