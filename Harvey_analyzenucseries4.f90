program analyzenucseries4
! similar to analyzenucseries2, with subroutine to find the transition probability W(s,s')
   ! similar to analyzenucseries3 but written for 3d data
   ! find histogram of changes in cluser size after every MC step per spin
   ! read time series of s, M
   implicit none
   integer, parameter :: r8 = selected_real_kind(8)
   ! P, probability of largest cluster during run
   integer, dimension(:), allocatable :: P
   ! probability that order parameter goes from s to s'
   real (kind=8), dimension(:,:), allocatable :: W
   real (kind=8) :: normW, left, right, ps, psprime
   real (kind=8) :: M
   integer :: stat, s, smax, tmax, L, N, sold
   integer :: t, i, nskip, sp
   real (kind=8) :: Pdecay, Psame, Pgrowth, normP
   character(len = 80) :: line
   tmax = 158000
   nskip = 0
   smax = 1939   ! largest cluster
   allocate(P(smax))
   allocate(W(smax,smax))
   P = 0
   normP = 0       ! counts # times that largest cluster is counted
   W = 0
   open(10,file='Rashidata/time_clust.txt',status='old')
   do i = 1, nskip
      read (unit=10,fmt = '(a80)',iostat=stat) line
   end do
   t = 0
   do while(.true.)
      read (unit=10,fmt = '(a80)',iostat=stat) line
      if (stat /= 0) exit
      !read(line, *) s, M
      read(line, *) M, s
      P(s) = P(s) + 1
      normP = normP + 1
      if (t >= 1) then
         W(sold, s) = W(sold, s) + 1
      end if
      sold = s
      t = t + 1
   end do
   close(unit=10)
   print *, "transition probabilities"
   print *, "s      s'      W(s,s')"
   ! W(s, s') is number of times s \to s'
   do s = 1, smax
      do sp = 1, smax
         if (W(s, sp) > 0) print *, s, sp, W(s,sp)
      end do
   end do
   print *, "s,   probability of growth, decay, no change"
   do s = 1, smax
      Pdecay = 0.
      Pgrowth = 0.
      Psame = 0.
      normP = 0.
      do sp = 1, smax
         if (W(s, sp) > 0) then
            if (sp < s) then
               Pdecay = Pdecay + W(s,sp)
            else if (sp == s) then
               Psame = Psame + W(s,sp)
            else
               Pgrowth = Pgrowth + W(s,sp)
            end if
            normP = Pgrowth + Pdecay + Psame
            Pgrowth = Pgrowth/normP
            Pdecay = Pdecay/normP
            Psame = Psame/normP
         end if
      end do     ! end sum over all s'
   if (normP > 0) print *, s, Pgrowth, Pdecay, Psame
   end do
   print *, "s     P(s)"
   normP = sum(P)
   do s = 1, smax
      if (P(s) > 0) print *, s, P(s)/normP
   end do
   print *, "check on detailed balance"
   do s = 1, smax
      Ps = P(s)/dble(normP)
      do sp = 2, smax
         if ( W(s, sp) > 0 .and. (W(sp,s) > 0) .and. (s .ne. sp) ) then
            Psprime = P(sp)/dble(normP)
            normW = sum(W(s,:))   ! number of times from s to any s'
            left = Ps*W(s,sp)/dble(normW)
            normW = sum(W(:,sp))   ! mumber of times from sp to any s
            right = Psprime*W(sp,s)/dble(normW)
            normW =  sum(W(:,sp))
         print *, s, sp, left, right
         end if
      end do
   end do
end program analyzenucseries4