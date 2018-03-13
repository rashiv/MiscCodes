program SpherHarm
implicit none

real, parameter :: pi = 3.1415927
complex, parameter :: iota = (0,1)

contains
!-----------------------------------------------------------!

function fac(q,steps)
integer :: q, steps, n
real :: fac
fac = 1.0
if(q.ge.2) then
  do n = q,2,-steps
    fac = fac*n
  end do
end if
end function fac
!--------------------------------------------------------------!
function Ylm(theta,phi,l,m)

real :: prefact, theta, phi
complex :: Ylm
integer :: l,m

prefact = -1**((m+abs(m))/2.0)*sqrt((2*l+1)*fac(l-abs(m),1)/(4*pi*fac(l+abs(m),1)))

Ylm = prefact*Plm(theta,l,abs(m))*exp(iota*m*phi)

end function Ylm 
!---------------------------------------------------------------!

function Plm(theta,l1,m1)
!======================================
! calculates Legendre polynomials Plm(x)
! using the recurrence relation
! if l > 100 the function returns 0.0
!======================================
real :: theta, x1, Plm
real, allocatable :: P(:,:)
integer :: l,m,l1,m1

allocate(P(0:l1,-l1:l1))

x1 = cos(theta)
P(0,0) = 1.0
P(1,0) = x1
do m = 1,l1
  P(m,m) = -1**m*fac((2*m-1),2)*(1-x1**2)**(m/2.0)
  if(m<l1) P(m+1,m) = x1*(2*m+1)*P(m,m)
end do
do l = 2,l1
  do m = 0,l-2
    P(l,m) = (x1*(2*l-1)*P(l-1,m)-(l+m-1)*P(l-2,m))/(l-m)    
  end do
end do
do l = 1,l1
  do m = 1,l
    P(l,-m) = -1**m*fac(l-m,1)/fac(l+m,1)*P(l,m)
  end do
end do

Plm = P(l1,m1)

end function Plm
!--------------------------------------------------------------------!
end program SpherHarm
