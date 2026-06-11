
 subroutine fv3wam_dry_adjust(levs, t, kappa, prsl, prsi, delp)
!     
! vay oct 2023 for conv-instability in the MLT 
! make dry adjustment of tkin(levs) +> convectively stable configuration
! "warm the conv-unstable layers": relax unstable T-re profile => T-stab
!  Apply GW-scheme only for in the conv-stable atmosphere
!  

      implicit none
      integer             , parameter ::    k_trop = 40       ! SMT-domain
      integer             , parameter ::    niter =  15                
      integer, intent(in)  :: levs
      real, intent(in)    ::  kappa(levs), prsl(levs), delp(levs)
      real, intent(in)    ::  prsi(levs+1)      
      real, intent(inout) ::  t(levs)  
! local
!
!
      real                      ::    gammad, dtdp , zeps, zepsdp, zgamma, rdenom    
      logical                   ::     dodad, ilconv 
      integer                   ::   k, jiter
      real, dimension(levs)     ::   c1dad, c2dad, c3dad, c4dad
      real                      ::   dtm, dtins
! from exobase
! to k_trop define the "instability" in the column
!        

       zeps = 2.0e-5     ! set convergence criteria     
       dtdp = (t(levs-1)-t(levs))/(prsl(levs-1)-prsl(levs))
       gammad = .25*(kappa(levs-1)+kappa(levs))/(t(levs-1)+t(levs))/prsi(levs)  
       dodad = (dtdp + zeps) .gt. gammad
!  do i = 1, ncol
!      cappa = 0.5*(cappav(i,2) + cappav(i,1))
!      gammad = cappa*0.5*(t(i,2) + t(i,1))/pint(i,2)
!      dtdp = (t(i,2) - t(i,1))/(pmid(i,2) - pmid(i,1))
!      dodad(i) = (dtdp + zeps) .gt. gammad
!   end do
       
   do k=2, k_trop   
         gammad = .25*(kappa(k-1)+kappa(k))/(t(k-1)+t(k))/prsi(k)
         dtdp = (t(k-1) - t(k))/(prsl(k-1) - prsl(k))
         dodad = dodad .or. (dtdp + zeps).gt.gammad
   end do
   
!         do k = 1, nlvdry
!            c1dad(k) = cappa*0.5*(pmid(i,k+1)-pmid(i,k))/pint(i,k+1)
!            c2dad(k) = (1. - c1dad(k))/(1. + c1dad(k))
!            rdenom = 1./(pdel(i,k)*c2dad(k) + pdel(i,k+1))
!            c3dad(k) = rdenom*pdel(i,k)
!            c4dad(k) = rdenom*pdel(i,k+1)
!         end do
   
      if( .not. dodad) return
!
! dodat=.T. => make dry adjustment can be "multiple" layers in FV3       
         zeps = 2.0e-5
         do k=1, k_trop
            c1dad(k) = 0.25*(kappa(k+1)+kappa(k))*(prsl(k+1)-prsl(k))/prsi(k)
            c2dad(k) = (1. - c1dad(k))/(1. + c1dad(k))
            rdenom = 1./(delp(k)*c2dad(k) + delp(k-1))
            c3dad(k) = rdenom*delp(k)
            c4dad(k) = rdenom*delp(k+1)
         end do
	 
50       do jiter=1,niter
            ilconv = .true.
            do k=2,k_trop
               zepsdp = zeps*(prsl(k+1) - prsl(k))
               zgamma = c1dad(k)*(t(k+1) + t(k))
	       dtins = zgamma + zepsdp
	       dtm = t(k-1)-t(k)
               if ( dtm >= dtins) then
                  ilconv = .false.
                  t(k) = t(k)*c3dad(k) + t(k+1)*c4dad(k)
                  t(k+1) = c2dad(k)*t(k)
               endif
            end do
            if (ilconv) go to 80           ! convergence => next longitude
         end do
!
! Crude (double) convergence criterion if no convergence in niter iterations
!
          zeps = zeps + zeps
	 
         if (zeps > 1.e-4) then
!            write(6,*)'DADADJ: No convergence in dry adiabatic adjustment'
           else
!            write(6,810) zeps
            go to 50
         endif
      
80    continue
!
!
810   format(//,'FV3WAM DADADJ: Convergence doubled to EPS=',E9.4)
!   COL: do i = 1, ncol
!      if (dodad(i)) then
!         zeps = 2.0e-5_r8
!         do k = 1, nlvdry
!            c1dad(k) = cappa*0.5_r8*(pmid(i,k+1)-pmid(i,k))/pint(i,k+1)
!            c2dad(k) = (1._r8 - c1dad(k))/(1._r8 + c1dad(k))
!            rdenom = 1._r8/(pdel(i,k)*c2dad(k) + pdel(i,k+1))
!            c3dad(k) = rdenom*pdel(i,k)
!            c4dad(k) = rdenom*pdel(i,k+1)
!         end do
!50       continue
!
!         do jiter = 1, niter
!            ilconv = .true.

!!            do k = 1, nlvdry
!               zepsdp = zeps*(pmid(i,k+1) - pmid(i,k))
!              zgamma = c1dad(k)*(t(i,k) + t(i,k+1))
!
!               if ((t(i,k+1)-t(i,k)) >= (zgamma+zepsdp)) then
!                  ilconv = .false.
!                  t(i,k+1) = t(i,k)*c3dad(k) + t(i,k+1)*c4dad(k)
!                  t(i,k) = c2dad(k)*t(i,k+1)
!               end if
!
!            end do
!            if (ilconv) cycle COL ! convergence => next longitude
!         end do
	     
     return
     end subroutine fv3wam_dry_adjust
