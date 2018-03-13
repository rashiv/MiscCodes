module clusterstuff
! Wang, Gould, Klein, Physical Review E 76, 031604 (2007)
! 25 August 2014. program based on Rashi's program, but added subroutines
! 8 September 2014. Add subroutines to analyze hot particles
implicit none
public :: initialize, computeSeparations, computeQlm, determineSolidlikeParticles
public :: findSolidClusters, findHotClusters
public :: analyzeParticles
public :: factorial, Ylm, Plm, findMinLabel, assignLabels, reduce
integer, parameter :: r8 = selected_real_kind(8)
real, parameter :: pi = 3.141592653589793238462643383279502884197169399
real (kind=r8), public :: L
integer, public :: N
!complex, parameter :: iota = (0.0,1.0)
complex, allocatable, public :: qlmtilde(:,:), qlm(:,:)
logical, allocatable, public :: solid(:), hot(:)

contains

subroutine initialize(x, y, z, nndist, cijThreshold, coherentNeighbors)
   real (kind=r8), allocatable, intent(out) :: x(:), y(:), z(:)
   real (kind=r8), intent(out) :: nndist, cijThreshold
   integer, intent(out) :: coherentNeighbors
   real (kind=r8) :: density, a, temperature, multiplier
   integer :: i, dummy, stat
   character(len=70) :: dum
   character(len = 80) :: line
   open(10,file='colloidConf4.5.dat',status='old')
   read(10,*) N, L, L, L, temperature      ! for  colloidConf1.dat
   print *, "N=", N
   print *, "L=", L
   density = N/L**3
   print *, "number density=", density
   a = (2.0/density)**(1.0/3.0)    ! lattice spacing
   print *, "lattice spacing=", a
   print *, "temperature=", temperature
   allocate(x(N))
   allocate(y(N))
   allocate(z(N))
   allocate(solid(N))       ! true if particle is solid-like
   allocate(hot(N))         ! true according to Jan's criterion
   solid = .false.
   hot = .false.
   multiplier = 1.13
   nndist = multiplier*a*sqrt(3.0)/2.0      ! nearest neighbor distance
   print *, "nearest neighbor distance multiplier=", multiplier
   cijThreshold = 0.5
   print *, "cijThreshold=", cijThreshold
   coherentNeighbors = 8
   print *, "coherentNeighbors=", coherentNeighbors
   i = 1
   do while(.true.)
      read (unit=10,fmt = "(a80)",iostat=stat) line
      if (stat /= 0) exit
      read(line, *) x(i), y(i), z(i), dummy
      if (dummy == 1) hot(i) = .true.
      i = i + 1
      end do
   close(unit=10)
end subroutine initialize

subroutine computeSeparations(x, y, z, nndist, nni, nnparticles)
   ! compute separation between all atoms
   real (kind=r8), dimension(:), intent(in) :: x, y, z
   real (kind=r8), intent(in) :: nndist
   integer, intent(out) :: nni(:), nnparticles(:,:)
   real (kind=r8) :: dx, dy, dz, r2, r
   real (kind=r8) :: phi, theta
   real (kind=r8) :: ql
   integer :: i, j, m
   nni = 0      ! number of nearest neighbors of particle i
   nnparticles = 0     ! stores nearest neighbors of particle i
   qlm =  (0.0,0.0)
   do i = 1, N
      ql = 0.0
      do j = 2, N
         if (i /= j) then
            dx = x(i) - x(j)
            dx = dx - L*anint(dx/L)
            dy = y(i) - y(j)
            dy = dy - L*anint(dy/L)
            dz = z(i) - z(j)
            dz = dz - L*anint(dz/L)
            r2 = dx**2 + dy**2 + dz**2
            r = sqrt(r2)
            ! unit vectors
            dx = dx/r
            dy = dy/r
            dz = dz/r
            if ( (r .le. nndist) .and. (i .ne. j) ) then
               nni(i) = nni(i) + 1       ! number of nearest neighbors of particle i
               nnparticles(i, nni(i)) = j     ! nearest neighbor indices of particle i
               phi = acos(dy)
               theta = atan(dx/dz)
               do m = -6, 6
                  qlm(i,m) = qlm(i,m) + Ylm(theta,phi,6,m)    ! q6 analysis for l = 6. Eq. (5)
               end do
            end if
         end if
      end do
   end do
end subroutine computeSeparations

subroutine computeQlm(nni, nnparticles)
   integer, intent(in) :: nni(:), nnparticles(:,:)
   real (kind=r8) :: phi, theta
   real (kind=r8) :: ql
   integer :: i, j, m
   qlmtilde = (0.0,0.0)
   do i = 1, N
      ql = 0.0
      do m = -6, 6
         if (nni(i) > 0) qlm(i,m) = qlm(i,m)/nni(i)      ! Eq. (5) \overline{q}_{lm}
      end do
      do m = -6, 6
         ql = ql + qlm(i,m)*conjg(qlm(i,m))      ! Eq. (6)
      end do
      do m = -6, 6
         qlmtilde(i,m) = qlm(i,m)/sqrt(ql)          ! Eq. (9)
      end do
      ql = sqrt(ql*4.0*pi/(2.0*l + 1.0))         ! Eq. (6)
   end do
end subroutine computeQlm

subroutine determineSolidlikeParticles(nnparticles, nni, cijThreshold, coherentNeighbors, nsolid)
   integer, intent(in) :: nnparticles(:,:), nni(:)
   real (kind=r8), intent(in) :: cijThreshold
   real (kind=r8) :: meanhotcij, meansolidcij, meanliquidcij
   integer, intent(in) :: coherentNeighbors
   integer, intent(out) :: nsolid
   integer :: i, j, m, ncoherent, nhot
   complex :: cij
   nsolid = 0
   nhot = 0
   meanhotcij = 0.0
   meansolidcij = 0.0
   do i = 1, N
      ncoherent = 0   ! number of coherent neighbors of particle i
      do j = 1, nni(i)     ! compare ith atom with its neighbors
         cij = (0.0,0.0)
         do m = -6, 6
            cij = cij + qlmtilde(i,m)*conjg(qlmtilde(nnparticles(i,j),m))   ! Eq. (10)
        end do
        ! criterion for j being a coherent neighbor of i
        if ( (real(cij) > cijThreshold) .and. (i .ne. j) ) ncoherent = ncoherent + 1
      end do     ! finish all neighbors
      if (ncoherent .ge. coherentNeighbors) then
         solid(i) = .true.      ! particle i is solid-like
         nsolid = nsolid + 1
      end if
      ! compute statistics for Jan's hot particles
      if (hot(i)) then
         meanhotcij = meanhotcij + real(cij)
         nhot = nhot + 1
      end if
      if (solid(i)) then
         meansolidcij = meansolidcij + real(cij)
      else
         meanliquidcij = meanliquidcij + real(cij)
      end if
   end do
   print *, "number of solid-like particles=", nsolid
   print *, "number of liquid-like particles=", N - nsolid
   print *, "number of Jan's hot particles=", nhot
   if (nhot > 0) print *, "mean cij hot=", meanhotcij/nhot
   print *, "mean cij solid=", meansolidcij/nsolid
   print *, "mean cij liquid-like=", meanliquidcij/(N-nsolid)
end subroutine determineSolidlikeParticles
   
subroutine findSolidClusters(nni, nnparticles)
   integer, intent(in) :: nni(:), nnparticles(:,:)
   integer, allocatable :: clustersize(:), label(:), pointer(:)
   integer :: next, i, j, nn, min_label, labelMax, iLabel
   allocate(clustersize(N))
   allocate(label(N))
   allocate(pointer(N))
   label = 0
   pointer = 0
   clustersize = 0
   next = 1   ! next cluster label
   do i = 1, N
      if (solid(i)) then
         nn = nni(i)    ! number of nearest neighbors of particle i
         ! find minimum label of particle i and its neighbors
         call findMinLabel(i, nn, nnparticles, label, pointer, min_label)
         if (min_label == N + 1) then
            ! new cluster
            label(i) = next
            pointer(next) = next
            min_label = label(i)
            next = next + 1
         end if
         ! assign labels to neighbors
         call assignLabels(i, nn, nnparticles, label, pointer, min_label)
      end if
   end do
   labelMax = 0
   clustersize = 0
   do i = 1, N
      if (solid(i)) then
         j = reduce(label(i), pointer)    ! just to make sure
         if (j > labelMax) labelMax = j
         clustersize(j) = clustersize(j) + 1
      end if
   end do
   do i = 1, N
      if (clustersize(i) > 0) then
         print *, "clusters", i, clustersize(i)
         ! uncomment the following to print labels of particles in each cluster
         !do j = 1, N
            !if (label(j) == label(i)) print *, j  ! labels of solid-like particle in each cluster
         !end do
      end if
   end do
end subroutine findSolidClusters

subroutine findHotClusters(nni, nnparticles)
   integer, intent(in) :: nni(:), nnparticles(:,:)
   integer, allocatable :: clustersize(:), label(:), pointer(:)
   integer :: next, i, j, nn, min_label, labelMax, iLabel, liquidLike
   allocate(clustersize(N))
   allocate(label(N))
   allocate(pointer(N))
   label = 0
   pointer = 0
   clustersize = 0
   liquidLike = 0
   next = 1   ! next cluster label
   do i = 1, N
      if (.not. solid(i)) then
      liquidLike = liquidLike + 1
         !print *, i
         nn = nni(i)    ! number of nearest neighbors of particle i
         ! find minimum label of particle i and its neighbors
         call findMinLabel(i, nn, nnparticles, label, pointer, min_label)
         if (min_label == N + 1) then
            ! new cluster
            label(i) = next
            pointer(next) = next
            min_label = label(i)
            next = next + 1
         end if
         ! assign labels to neighbors
         call assignLabels(i, nn, nnparticles, label, pointer, min_label)
      end if
   end do
   labelMax = 0
   clustersize = 0
   do i = 1, N
      if (.not. solid(i)) then
         j = reduce(label(i), pointer)    ! just to make sure
         if (j > labelMax) labelMax = j
         clustersize(j) = clustersize(j) + 1
      end if
   end do
   do i = 1, N
      if (clustersize(i) > 0) then
         !print *, "cluster", i, clustersize(i)
         !do j = 1, N
            !if (label(j) == label(i)) print *, j
         !end do
      end if
   end do
end subroutine findHotClusters

subroutine analyzeParticles(nni)
   integer, intent(in) :: nni(:)
   real (kind=r8) :: meanall, meansolid, meanliquid, meanhot
   integer, dimension(16) :: histogram, histogramhot, histogramsolid, histogramliquid
   integer :: i, nn, nsolid, nliquid, nhot, nnmax
   histogram = 0
   histogramsolid = 0
   histogramliquid = 0
   histogramhot = 0
   nsolid = 0
   nliquid = 0
   nhot = 0
   nnmax = 16
   do i = 1, N
      nn = nni(i)    ! number of nearest neighbors of particle i
      histogram(nn) = histogram(nn) + 1     ! all particles
      if (hot(i)) then
         histogramhot(nn) = histogramhot(nn) + 1
         nhot = nhot + 1
      end if
       if (solid(i)) then
          histogramsolid(nn) = histogramsolid(nn) + 1
          nsolid = nsolid + 1
       else
          histogramliquid(nn) = histogramliquid(nn) + 1
          nliquid = nliquid + 1
       end if
   end do
   meanhot = 0.0
   print *, "nearest neighbor histogram of hot particles"
   do i = 1, nnmax
      if (histogramhot(i) > 0) print *, i, histogramhot(i)
      meanhot = meanhot + i*histogramhot(i)
   end do
   if (sum(histogramhot) > 0) then
      print *, "mean number of nearest neighbors of hot particles=", meanhot/nhot
   end if
   meanall = 0.0
   meansolid = 0.0
   meanliquid = 0.0
   print *, "neighbors       all particles    solid      liquid"
   do i = 1, nnmax
      if (histogram(i) > 0) print *, i, histogram(i), histogramsolid(i), histogramliquid(i)
      meanall = meanall + i*histogram(i)
      meansolid = meansolid + i*histogramsolid(i)
      meanliquid = meanliquid + i*histogramliquid(i)
   end do
   print *, "mean number of neighbors of all particles=", meanall/N
   print *, "mean number of neighbors of solid particles=", meansolid/nsolid
   print *, "mean number of neighbors of liquid particles=", meanliquid/nliquid
end subroutine analyzeParticles

subroutine findMinLabel(i, nn, nnparticles, label, pointer, min_label)
   ! find minimum nonzero cluster label of particle i and its neighbors
   integer, intent(in) :: i, nn
   integer, intent(in) :: nnparticles(:,:)
   integer, intent(inout) :: label(:), pointer(:)
   integer, intent(out) :: min_label
   integer :: k, neighbor
   if (label(i) > 0) then
      label(i) = reduce(label(i), pointer)
      min_label = label(i)
   else
      min_label= N + 1     ! 1 if no nearest neighbors
   end if
   do k = 1, nn      ! nn is number of neighbors of particle i
      neighbor = nnparticles(i, k)
      if (solid(neighbor)) then
         if (label(neighbor) > 0) then
            label(neighbor) = reduce(label(neighbor), pointer)
            if (label(neighbor) < min_label) then
               min_label = label(neighbor)
            end if
         end if
      end if
   end do
end subroutine findMinLabel

subroutine assignLabels(i, nn, nnparticles, label, pointer, min_label)
   integer, intent(in) :: i, nn
   integer, intent(in) :: nnparticles(:,:)
   integer, intent(inout) :: label(:), pointer(:)
   integer, intent(in) :: min_label
   integer :: k, neighbor
   if (label(i) > 0) pointer(label(i)) = min_label
   label(i) = min_label
   do k = 1, nn
      neighbor = nnparticles(i, k)
      if (solid(neighbor)) then
         if (label(neighbor) > 0) pointer(label(neighbor)) = min_label
         label(neighbor) = min_label
      end if
   end do
end subroutine assignLabels

real (kind=8) function factorial(q,steps)
   integer, intent(in) :: q, steps
   integer :: k
   factorial = 1.0
   if (q .ge. 2) then
      do k = q, 2,-steps
      factorial = factorial*k
     end do
   end if
end function factorial

complex function Ylm(theta, phi, l, m)
   real (kind=8), intent(in) :: theta, phi
   integer, intent(in) :: l, m
   complex :: iota = (0.0,1.0)
   real (kind=8) :: prefact
   prefact = -1**((m + abs(m))/2.0)*sqrt((2*l+1)*factorial(l-abs(m),1)/(4*pi*factorial(l+abs(m),1)))
   Ylm = prefact*Plm(theta,l,abs(m))*exp(iota*m*phi)
end function Ylm

function Plm(theta, l1, m1)
   ! compute Legendre polynomials P_{lm}(x) using recurrence relation
   ! if l > 100, the function returns 0.0
   real (kind=8), intent(in) :: theta
   integer, intent(in) :: l1, m1
   real (kind=8), allocatable :: P(:,:)
real (kind=8) :: x1, Plm
   integer :: l, m
   allocate(P(0:l1,-l1:l1))
   x1 = cos(theta)
   P(0,0) = 1.0
   P(1,0) = x1
   do m = 1, l1
      P(m,m) = -1**m*factorial((2*m-1),2)*(1-x1**2)**(m/2.0)
     if (m < l1) P(m+1,m) = x1*(2*m+1)*P(m,m)
   end do
   do l = 2, l1
      do m = 0, l-2
          P(l,m) = (x1*(2*l-1)*P(l-1,m)-(l+m-1)*P(l-2,m))/(l-m)
     end do
   end do
   do l = 1, l1
      do m = 1, l
         P(l,-m) = -1**m*factorial(l-m,1)/factorial(l+m,1)*P(l,m)
     end do
   end do
   Plm = P(l1,m1)
end function Plm

integer function reduce(k, pointer)
   integer, intent(inout) :: k
   integer, dimension(:), intent(in) :: pointer
   if (k .ne. 0) then
      do while (pointer(k) /= k)
         k = pointer(k)
      end do
   end if
   reduce = k
end function reduce

end module clusterstuff

program solidLike
   use clusterstuff
   implicit none
   real (kind=r8), allocatable :: x(:), y(:), z(:)
   integer, allocatable :: nni(:), nnparticles(:,:)
   real (kind=r8) :: nndist, cijThreshold
   integer :: coherentNeighbors, nsolid
   integer :: i, noverlap
   call initialize(x, y, z, nndist, cijThreshold, coherentNeighbors)
   allocate(nni(N))         ! number of nearest neighbors of particle i
   allocate(nnparticles(N,100))     ! list of neighbors of particle i
   allocate(qlm(N,-6:6))
   allocate(qlmtilde(N,-6:6))
   call computeSeparations(x, y, z, nndist, nni, nnparticles)
   call computeQlm(nni, nnparticles)
   deallocate(qlm)
   deallocate(x)
   deallocate(y)
   deallocate(z)
   call determineSolidlikeParticles(nnparticles, nni, cijThreshold, coherentNeighbors, nsolid)
   deallocate(qlmtilde)
   call findHotClusters(nni, nnparticles)
   !call findSolidClusters(nni, nnparticles)
   noverlap = 0     ! # particles that are both hot and liquid-like   ! 
   do i = 1, N
      if ( (.not. solid(i)) .and. (hot(i)) ) then
            noverlap = noverlap + 1
            !print *, i      ! labels of particles that are hot and liquid-like
      end if
   end do
   print *, "number of particles that are both liquid-like and hot=", noverlap
   print *, "percentaage overlap of liquid-like and hot particles=", real(noverlap)/real(N - nsolid)
   call analyzeParticles(nni)
end program solidLike