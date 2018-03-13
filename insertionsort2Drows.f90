program insertionsort

integer:: n,dist(3,3)
integer, allocatable :: sdist(:,:),sindex(:,:)
n=3
data dist /10,20,60,40,50,30,30,10,20/

allocate(sdist(3,3))
allocate(sindex(3,3))

call sortrows(dist,n)

contains

subroutine sortrows(a,length)
!subroutine to sort a 2-D array row-wise
integer :: length,row,i,j,a(length,length),temp
integer, allocatable :: sdist(:,:),sindex(:,:)

allocate(sdist(length,length))
allocate(sindex(length,length))

sdist = a

do i=1,length
  do j=1,length
    sindex(i,j) = j
  end do
end do

print *, "original matrix:"
print *, sdist
print *, "original indices:"
print *, sindex

do row = 1, length
  do i = 2, length
    j = i-1
    temp = sdist(row,i)
    do while (j>=1.and.sdist(row,j)>temp)
        sdist(row,j+1) = sdist(row,j)
        sindex(row,j+1) = sindex(row,j)
        j = j-1
    end do
    sdist(row,j+1) = temp
    sindex(row,j+1) = i
  end do
end do

print *, "sorted matrix:"
print *, sdist
print *, "sorted indices:"
print *, sindex

end subroutine sortrows

end program testsort
