!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!>@brief The module 'fv_mapz' contains the vertical mapping routines \cite lin2004vertically
!>@note April 12, 2012 -SJL: This revision may actually produce rounding level differences
!! due to the elimination of KS to compute pressure level for remapping.

module fv_mapz_mod

! Modules Included:
! <table>
! <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
! <table>
!   <tr>
!     <td>constants_mod</td>
!     <td>radius, pi=>pi_8, rvgas, rdgas, grav, hlv, hlf, cp_air, cp_vapor</td>
!   </tr>
!     <td>field_manager_mod</td>
!     <td>MODEL_ATMOS</td>
!   </tr>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td>fv_grid_type</td>
!   </tr>
!   <tr>
!     <td>fv_fill_mod</td>
!     <td>fillz</td>
!   </tr>
!   <tr>
!     <td>fv_grid_utils_mod</td>
!     <td>g_sum, ptop_min</td>
!   </tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>is_master</td>
!   </tr>
!   <tr>
!     <td>ccpp_static_api</td>
!     <td>ccpp_physics_run</td>
!   </tr>
!   <tr>
!     <td>CCPP_data</td>
!     <td>ccpp_suite, cdata_tile, GFDL_interstitial</td>
!   </tr>
!   <tr>
!     <td>fv_timing_mod</td>
!     <td>timing_on, timing_off</td>
!   </tr>
!   <tr>
!     <td>fv_tracer2d_mod</td>
!     <td>tracer_2d, tracer_2d_1L, tracer_2d_nested</td>
!   </tr>
!   <tr>
!     <td>mpp_mod/td>
!     <td>NOTE, mpp_error, get_unit, mpp_root_pe, mpp_pe</td>
!   </tr>
!   <tr>
!     <td>mpp_domains_mod/td>
!     <td> mpp_update_domains, domain2d</td>
!   </tr>
!   <tr>
!     <td>tracer_manager_mod</td>
!     <td>get_tracer_index</td>
!   </tr>
! </table>

  use constants_mod,     only: radius, pi=>pi_8, rvgas, rdgas, grav, hlv, hlf, cp_air, cp_vapor
  use tracer_manager_mod,only: get_tracer_index
  use field_manager_mod, only: MODEL_ATMOS
  use fv_grid_utils_mod, only: g_sum, ptop_min, cubed_to_latlon, update_dwinds_phys
  use fv_fill_mod,       only: fillz
  use mpp_domains_mod,   only: mpp_update_domains, domain2d
  use mpp_mod,           only: NOTE, FATAL, mpp_error, get_unit, mpp_root_pe, mpp_pe
  use fv_arrays_mod,     only: fv_grid_type, fv_grid_bounds_type, R_GRID, inline_mp_type
  use fv_timing_mod,     only: timing_on, timing_off
  use fv_mp_mod,         only: is_master, mp_reduce_min, mp_reduce_max
  ! CCPP fast physics
  use ccpp_static_api,   only: ccpp_physics_run
  use CCPP_data,         only: ccpp_suite
  use CCPP_data,         only: cdata => cdata_tile
  use CCPP_data,         only: GFDL_interstitial
#ifdef MULTI_GASES
  use multi_gases_mod,  only:  virq, virqd, vicpqd, vicvqd, num_gas
#endif
  use molecular_diffusion_mod, only : mdz_uv, mdz_tk,  mdz_w,   mdz_q4
  
  implicit none
  real, parameter:: consv_min= 0.001         !< below which no correction applies
  real, parameter:: t_min= 104.              !< below which applies stricter constraint 184K FV3..... WAM => 104K
  real, parameter:: r3 = 1./3., r23 = 2./3., r12 = 1./12.
  real, parameter:: cv_vap = 3.*rvgas        !< 1384.5
  real, parameter:: cv_air =  cp_air - rdgas !< = rdgas * (7/2-1) = 2.5*rdgas=717.68
! real, parameter:: c_ice = 2106.            !< heat capacity of ice at 0.C
  real, parameter:: c_ice = 1972.            !< heat capacity of ice at -15.C
  real, parameter:: c_liq = 4.1855e+3        !< GFS: heat capacity of water at 0C
! real, parameter:: c_liq = 4218.            !< ECMWF-IFS
  real, parameter:: cp_vap = cp_vapor        !< 1846.
  real, parameter:: tice = 273.16

  real, parameter :: w_max = 90.
  real, parameter :: w_min = -90.
  logical, parameter :: w_limiter = .false. ! doesn't work so well??

  real(kind=4) :: E_Flux = 0.
  private

  public compute_total_energy, Lagrangian_to_Eulerian, moist_cv, moist_cp,   &
         rst_remap, mappm, E_Flux, remap_2d, map_scalar

contains

!>@brief The subroutine 'Lagrangian_to_Eulerian' remaps deformed Lagrangian layers back to the reference Eulerian coordinate.
!>@details It also includes the entry point for calling fast microphysical processes. This is typically calle on the k_split loop.
 subroutine Lagrangian_to_Eulerian(last_step, consv, ps, pe, delp, pkz, pk,   &
                                   mdt, pdt, npx, npy, km, is,ie,js,je, isd,ied,jsd,jed,       &
                      nq, nwat, sphum, q_con, u, v, w, delz, pt, q, hs, r_vir, cp,  &
                      akap, cappa, kord_mt, kord_wz, kord_tr, kord_tm,  peln, te0_2d,        &
                      ng, ua, va, omga, te, ws, fill, reproduce_sum, out_dt, dtdt,      &
                      ptop, ak, bk, pfull, gridstruct, domain, do_sat_adj, &
                      hydrostatic, phys_hydrostatic, hybrid_z, do_omega, adiabatic, do_adiabatic_init, &
                      do_inline_mp, inline_mp, c2l_ord, bd, fv_debug, &
                      moist_phys)
		      
  logical, intent(in):: last_step
  logical, intent(in):: fv_debug
  real,    intent(in):: mdt                   !< remap time step
  real,    intent(in):: pdt                   !< phys time step
  integer, intent(in):: npx, npy
  integer, intent(in):: km
  integer, intent(in):: nq                     !< number of tracers (including h2o)
  integer, intent(in):: nwat
  integer, intent(in):: sphum                  !< index for water vapor (specific humidity)
  integer, intent(in):: ng
  integer, intent(in):: is,ie,isd,ied          !< starting & ending X-Dir index
  integer, intent(in):: js,je,jsd,jed          !< starting & ending Y-Dir index
  integer, intent(in):: kord_mt                !< Mapping order for the vector winds
  integer, intent(in):: kord_wz                !< Mapping order/option for w
  integer, intent(in):: kord_tr(nq)            !< Mapping order for tracers
  integer, intent(in):: kord_tm                !< Mapping order for thermodynamics
  integer, intent(in):: c2l_ord
  real, intent(in):: consv                     !< factor for TE conservation
  real, intent(in):: r_vir
  real, intent(in):: cp
  real, intent(in):: akap
  real, intent(in):: hs(isd:ied,jsd:jed)       !< surface geopotential
  real, intent(inout):: te0_2d(is:ie,js:je)
  real, intent(in):: ws(is:ie,js:je)

  logical, intent(in):: do_sat_adj
  logical, intent(in):: do_inline_mp
  logical, intent(in):: fill                    !< fill negative tracers
  logical, intent(in):: reproduce_sum
  
  logical, intent(in):: do_omega, adiabatic, do_adiabatic_init
  
  real, intent(in) :: ptop
  real, intent(in) :: ak(km+1)
  real, intent(in) :: bk(km+1)
  real, intent(in):: pfull(km)
  type(fv_grid_type), intent(IN), target :: gridstruct
  type(domain2d), intent(INOUT) :: domain
  type(fv_grid_bounds_type), intent(IN) :: bd

! INPUT/OUTPUT
  real, intent(inout):: pk(is:ie,js:je,km+1)          !< pe to the kappa
  real, intent(inout):: q(isd:ied,jsd:jed,km,*)
  real, intent(inout):: delp(isd:ied,jsd:jed,km)      !< pressure thickness
  real, intent(inout)::  pe(is-1:ie+1,km+1,js-1:je+1) !< pressure at layer edges
  real, intent(inout):: ps(isd:ied,jsd:jed)           !< surface pressure

! u-wind will be ghosted one latitude to the north upon exit
  real, intent(inout)::  u(isd:ied  ,jsd:jed+1,km)   !< u-wind (m/s)
  real, intent(inout)::  v(isd:ied+1,jsd:jed  ,km)   !< v-wind (m/s)
  real, intent(inout)::  w(isd:     ,jsd:     ,1:)   !< vertical velocity (m/s)
  real, intent(inout):: pt(isd:ied  ,jsd:jed  ,km)   !< cp*virtual potential temperature
                                                     !< as input; output: temperature
  real, intent(inout), dimension(isd:,jsd:,1:):: q_con, cappa
  real, intent(inout), dimension(is:,js:,1:)::delz
  logical, intent(in):: hydrostatic, phys_hydrostatic
  logical, intent(in):: hybrid_z
  logical, intent(in):: out_dt
  logical, intent(in):: moist_phys

  real, intent(inout)::   ua(isd:ied,jsd:jed,km)   !< u-wind (m/s) on physics grid
  real, intent(inout)::   va(isd:ied,jsd:jed,km)   !< v-wind (m/s) on physics grid
  real, intent(inout):: omga(isd:ied,jsd:jed,km)   !< vertical press. velocity (pascal/sec)
  real, intent(inout)::   peln(is:ie,km+1,js:je)   !< log(pe)
  real, intent(inout)::   dtdt(is:ie,js:je,km)
  real, intent(out)::    pkz(is:ie,js:je,km)       !< layer-mean pk for converting t to pt
  real, intent(out)::     te(isd:ied,jsd:jed,km)

  type(inline_mp_type), intent(inout):: inline_mp
  logical     :: remap_te=.false.
! !DESCRIPTION:
!
! !REVISION HISTORY:
! SJL 03.11.04: Initial version for partial remapping
!
!-----------------------------------------------------------------------
  real, allocatable, dimension(:,:,:) :: dp0, u0, v0
  real, allocatable, dimension(:,:,:) :: u_dt, v_dt
  real, dimension(is:ie,js:je):: te_2d, zsum0, zsum1
  real, dimension(is:ie,km)  :: q2, dp2, t0, w2
  real, dimension(is:ie,km+1):: pe1, pe2, pk1, pk2, pn2, phis
  real, dimension(isd:ied,jsd:jed,km):: pe4
  real, dimension(is:ie+1,km+1):: pe0, pe3, pe3u
  
  real, dimension(is:ie):: gsize, gz, cvm, qv
  
  real, dimension(is:ie+1, km) :: dpv, dpu, v2dis    
  real, dimension(is:ie,km+1)  :: vumol, ktmol, dfmol, rhomol, wgrav
  real, dimension(is:ie+1, km) :: u2dis
  real, dimension(is:ie,   km) :: ptdry, pkzdry
  real rcp, rg, rrg, bkh, dtmp, k1k, tpe, dlnp
  real :: delp_max, ps_max, ps_min, psij, ps_scal

  integer:: i,j,k
  integer:: kdelz
  integer:: ind_h2o, ind_o2, ind_o3p, ind_o3 
  integer:: nt, liq_wat, ice_wat, rainwat, snowwat, cld_amt, graupel, hailwat, ccn_cm3, iq, n, kmp, kp, k_next
  integer :: ierr

      ccpp_associate: associate( fast_mp_consv => GFDL_interstitial%fast_mp_consv, &
                                 kmp           => GFDL_interstitial%kmp            )

       k1k = rdgas/cv_air   ! akap / (1.-akap) = rg/Cv=0.4
        rg = rdgas
       rcp = 1./ cp
       rrg = -rdgas/grav
       ps_max   = 1080.e2
       ps_min   = 459.e2       
       delp_max = 40.5e2
       
       
       remap_te = .false.
       
!ind_h2o, ind_o2, ind_o3p, ind_o3 

       ind_o3 = get_tracer_index (MODEL_ATMOS, 'spo3')
       ind_o2 = get_tracer_index (MODEL_ATMOS, 'spo2')
       ind_o3p = get_tracer_index (MODEL_ATMOS, 'spo')
       ind_h2o = get_tracer_index (MODEL_ATMOS, 'sphum') 
       
       liq_wat = get_tracer_index (MODEL_ATMOS, 'liq_wat')
       ice_wat = get_tracer_index (MODEL_ATMOS, 'ice_wat')
       rainwat = get_tracer_index (MODEL_ATMOS, 'rainwat')
       snowwat = get_tracer_index (MODEL_ATMOS, 'snowwat')
       graupel = get_tracer_index (MODEL_ATMOS, 'graupel')
       hailwat = get_tracer_index (MODEL_ATMOS, 'hailwat')
       cld_amt = get_tracer_index (MODEL_ATMOS, 'cld_amt')
       ccn_cm3 = get_tracer_index (MODEL_ATMOS, 'ccn_cm3')

       if ( do_adiabatic_init .or. do_sat_adj ) then
            fast_mp_consv = (.not.do_adiabatic_init) .and. consv>consv_min
       endif

!$OMP parallel do default(none) shared(is,ie,js,je,km,pe,ptop,kord_tm,hydrostatic, &
!$OMP                                cp, pt,pk,rg,peln,q,nwat,liq_wat,rainwat,ice_wat,snowwat, &
!$OMP                                  graupel,hailwat,q_con,sphum,cappa,r_vir,rcp,k1k,delp, &
!$OMP                                  delz,akap,pkz,te,u,v,ps, gridstruct, last_step, remap_te, &
!$OMP                                  ak,bk,nq,isd,ied,jsd,jed,kord_tr,fill, adiabatic, &
#ifdef MULTI_GASES
!$OMP                                  num_gas,                                          &
#endif
!$OMP                   mdt, rhomol, wgrav ,vumol, ktmol, dfmol,ind_h2o, ind_o2, ind_o3p, ind_o3,   & 
!$OMP                   mdz_uv, mdz_tk,  mdz_w, mdz_q4, u2dis, v2dis, dpu, dpv, &
!$OMP                          delp_max, ps_max, ps_min, psij, ps_scal, &
!$OMP                          dlnp, tpe, hs,w,ws,kord_wz,do_omega,omga,rrg,kord_mt,pe4)    &
!$OMP                          private(qv,gz,cvm,kp,k_next,bkh,dp2, ptdry, pkzdry,  &
!$OMP                                      pe0,pe1,pe2,pe3,pe3u, pk1,pk2,pn2,phis,q2,w2)
  do 1000 j=js,je+1

     do k=1,km+1
        do i=is,ie
           pe1(i,k) = pe(i,k,j)
        enddo
     enddo
     
     do i=is,ie
        pe2(i,   1) = ptop
        pe2(i,km+1) =pe1(i,km+1)   ! pe(i,km+1,j)
     enddo

  if ( j /= (je+1) ) then
  
    if (  .not. remap_te ) then
       if ( kord_tm < 0 ) then
! Note: pt at this stage is Theta_v
! Transform virtual pt to virtual Temp

             if ( hydrostatic ) then
                 do k=1,km
                   do i=is,ie
#ifdef MULTI_GASES
                      pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
                      pkz(i,j,k) = exp(virqd(q(i,j,k,1:num_gas))/vicpqd(q(i,j,k,1:num_gas))*log(pkz(i,j,k)))
                      pt(i,j,k) = pt(i,j,k)*pkz(i,j,k)
#else
                      pt(i,j,k) = pt(i,j,k)*(pk(i,j,k+1)-pk(i,j,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
#endif
                   enddo
                 enddo
             else
!                              !NH Transform "density pt" to "density temp"
               do k=1,km
#ifdef MOIST_CAPPA
                 call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                                ice_wat, snowwat, graupel, hailwat, q, gz, cvm)
                  do i=is,ie
                     q_con(i,j,k) = gz(i)
#ifdef MULTI_GASES
                     cappa(i,j,k) = rdgas / ( rdgas + cvm(i)/virq(q(i,j,k,1:num_gas)) )
#else
                     cappa(i,j,k) = rdgas / ( rdgas + cvm(i)/(1.+r_vir*q(i,j,k,sphum)) )
#endif
                     pt(i,j,k) = pt(i,j,k)*exp(cappa(i,j,k)/(1.-cappa(i,j,k))*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
                 enddo
#else

                  do i=is,ie
#ifdef MULTI_GASES
                     pt(i,j,k) = pt(i,j,k)*exp(k1k*(virqd(q(i,j,k,1:num_gas))/vicvqd(q(i,j,k,1:num_gas))*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
#else
                     pt(i,j,k) = pt(i,j,k)*exp(k1k*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
#endif
                  enddo
! Using dry pressure for the definition of the virtual potential temperature
!                    pt(i,j,k) = pt(i,j,k)*exp(k1k*log(rrg*(1.-q(i,j,k,sphum))*delp(i,j,k)/delz(i,j,k)*    &
!                                              pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum))))
#endif
               enddo       ! k-loop
             endif         ! hydro test
           endif           ! kord_tm
	   
 else                      ! remap_te	   	     
     if ( hydrostatic ) then
     
           call pkez(km, is, ie, js, je, j, pe, pk, akap, peln, pkz, ptop)
! Compute cp_air*Tm + KE
           do k=1,km
                 do i=is,ie
#ifdef MULTI_GASES
                    pkz(i,j,k) = exp(virqd(q(i,j,k,1:num_gas))/vicpqd(q(i,j,k,1:num_gas))*log(pkz(i,j,k)))
#endif
                    te(i,j,k) = 0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                                                 v(i,j,k)**2+v(i+1,j,k)**2 -  &
                               (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j))  &
                              + cp_air*pt(i,j,k)*pkz(i,j,k)
                 enddo
           enddo
	   
     else                       !NH-energy  
              do k=km,1,-1
                 do i=is,ie
                    phis(i,k) = phis(i,k+1) - grav*delz(i,j,k)
                 enddo
#ifdef MOIST_CAPPA
                 call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                      ice_wat, snowwat, graupel,hailwat, q, gz(is:ie), cvm(is:ie))	
                 do i=is,ie		      	 
                    q_con(i,j,k) = gz(i)
                    cappa(i,j,k) = rdgas / ( rdgas + cvm(i)/(1.+r_vir*q(i,j,k,sphum)) )
                    pkz(i,j,k) = exp(cappa(i,j,k)/(1.-cappa(i,j,k))*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))		    
#ifdef MULTI_GASES
                    pkz(i,j,k) = exp(virqd(q(i,j,k,1:num_gas))/vicpqd(q(i,j,k,1:num_gas))*log(pkz(i,j,k)))
#endif		    
		    
                    te(i,j,k) = cvm(i)*pt(i,j,k)*pkz(i,j,k)/((1.+r_vir*q(i,j,k,sphum))*(1.-gz(i))) +     &
                         0.5 * w(i,j,k)**2 + 0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                         v(i,j,k)**2+v(i+1,j,k)**2 -  &
                         (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j)) +         &
                         0.5*(phis(i,k+1)+phis(i,k))
                 enddo
#else
                 do i=is,ie
                    pkz(i,j,k) = exp(k1k*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
#ifdef MULTI_GASES
                    pkz(i,j,k) = exp(virqd(q(i,j,k,1:num_gas))/vicpqd(q(i,j,k,1:num_gas))*log(pkz(i,j,k)))
#endif			    
                    te(i,j,k) = cv_air*pt(i,j,k)*pkz(i,j,k)/(1.+r_vir*q(i,j,k,sphum)) +     &
                         0.5 * w(i,j,k)**2 + 0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                         v(i,j,k)**2+v(i+1,j,k)**2 -  &
                         (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j)) +         &
                         0.5*(phis(i,k+1)+phis(i,k))
                 enddo
#endif
            enddo      !k-loop
           endif       !end hyd/NH choices	   
	   
       endif           ! .not.remap_te-IF
       
       
     if ( .not. hydrostatic ) then
           do k=1,km
              do i=is,ie
                 delz(i,j,k) = -delz(i,j,k) / delp(i,j,k) ! ="specific volume"/grav
              enddo
           enddo
      endif

! update ps
      do i=is,ie
         ps(i,j) = pe1(i,km+1)
!	 if (ps(i,j) > ps_max) then 
!	    ps(i,j) = ps_max
!	 endif 
	   pe2(i,km+1) = ps(i,j)	 
      enddo
!
! Hybrid sigma-P coordinate:
!
        do k=2,km
           do i=is,ie
              pe2(i,k) = ak(k) + bk(k)*ps(i,j)
           enddo
        enddo
        do k=1,km
           do i=is,ie
              dp2(i,k) = pe2(i,k+1) - pe2(i,k)
           enddo
        enddo

!------------
! update delp
!------------
      do k=1,km
         do i=is,ie
            delp(i,j,k) = dp2(i,k)
         enddo
      enddo

!------------------
! Compute p**Kappa
!------------------
   
      do i=is,ie
         pk1(i,1) = pk(i,j,1)
	 pk1(i,km+1) = exp(akap*alog(ps(i,j))) 
	 pn2(i,km+1) = log(ps(i,j))
      enddo
  
!============================================= BC-points 1 & km+1
   do i=is,ie
      pn2(i,   1) = peln(i,   1,j)
!      pn2(i,km+1) = peln(i,km+1,j)
      pk2(i,   1) = pk1(i,   1)
      pk2(i,km+1) = pk1(i,km+1)
   enddo

   do k=2,km
      do i=is,ie
         pn2(i,k) = log(pe2(i,k))
         pk2(i,k) = exp(akap*pn2(i,k))
      enddo
   enddo
!===================================== GFDL-way Dec 2023   
!      1) Remap Tv, thetav, or TE
!===================================== GFDL-way Dec 2023  
 if ( remap_te ) then
!----------------------------------
! Map TE   kord_tm == 0 "map1_cubic" >0 map_scalar
!----------------------------------  
         if ( kord_tm==0 ) then
!----------------------------------
! map Total Energy using GMAO cubic
!----------------------------------
            call map1_cubic (km,   pe1,  te,       &
                 km,   pe2,  te,       &
                 is, ie, j, isd, ied, jsd, jed, akap, T_VAR=1, conserv=.true.)
         else
            call map_scalar(km,  peln(is,1,j),  te, gz(is:ie),   &
                 km,  pn2,           te,              &
                 is, ie, j, isd, ied, jsd, jed, 1, abs(kord_tm), cp_air*t_min)
         endif
	 
 else	
!----------------------------------
! Map PT/T   kord_tm<0 
!---------------------------------- 
   if ( kord_tm<0 ) then
!----------------------------------
! Map t using logp
!----------------------------------
         call map_scalar(km,  peln(is,1,j),  pt, gz,   &
                         km,  pn2,           pt,              &
                         is, ie, j, isd, ied, jsd, jed, 1, abs(kord_tm), t_min)
   else
! Map pt using pe
         call map1_ppm (km,  pe1,  pt,  gz,       &
                        km,  pe2,  pt,                  &
                        is, ie, j, isd, ied, jsd, jed, 1, abs(kord_tm))
   endif
endif

!----------------
! Map constituents:  mapn_tracer
!----------------
      if( nq > 5 ) then
           call mapn_tracer(nq, km, pe1, pe2, q, dp2, kord_tr, j,     &
                            is, ie, isd, ied, jsd, jed, 0., fill)
      elseif ( nq > 0 ) then
! Remap one tracer at a time
         do iq=1,nq
             call map1_q2(km, pe1, q(isd,jsd,1,iq),     &
                          km, pe2, q2, dp2,             &
                          is, ie, 0, kord_tr(iq), j, isd, ied, jsd, jed, 0.)
            if (fill) call fillz(ie-is+1, km, 1, q2, dp2)
            do k=1,km
               do i=is,ie
                  q(i,j,k,iq) = q2(i,k)
               enddo
            enddo
         enddo
      endif
!----------------
! Map vert winds:  map1_ppm
!----------------
   if ( .not. hydrostatic ) then
! Remap vertical wind:
        call map1_ppm (km,   pe1,  w,  ws(is,j),   &
                       km,   pe2,  w,              &
                       is, ie, j, isd, ied, jsd, jed, -2, kord_wz)
		       
		       
! Remap delz for hybrid sigma-p coordinate
!----------------
! Map delz =1/dens:  map1_ppm
!----------------
        call map1_ppm (km,   pe1, delz,  gz,   & ! works
                       km,   pe2, delz,              &
                       is, ie, j, is,  ie,  js,  je,  1, abs(kord_tm))
        do k=1,km
           do i=is,ie
              delz(i,j,k) = -delz(i,j,k)*dp2(i,k)
           enddo
        enddo
!-------------------
!Fix excessive w - momentum conserving --- sjl
! gz(:) used here as a temporary array
!-------------------
        if ( w_limiter ) then
           do k=1,km
              do i=is,ie
                 w2(i,k) = w(i,j,k)
              enddo
           enddo
           do k=1, km-1
              do i=is,ie
                 if ( w2(i,k) > w_max ) then
                    gz(i) = (w2(i,k)-w_max) * dp2(i,k)
                    w2(i,k  ) = w_max
                    w2(i,k+1) = w2(i,k+1) + gz(i)/dp2(i,k+1)
                    print*, ' W_LIMITER down: ', i,j,k, w2(i,k:k+1), w(i,j,k:k+1)
                 elseif ( w2(i,k) < w_min ) then
                    gz(i) = (w2(i,k)-w_min) * dp2(i,k)
                    w2(i,k  ) = w_min
                    w2(i,k+1) = w2(i,k+1) + gz(i)/dp2(i,k+1)
                    print*, ' W_LIMITER down: ', i,j,k, w2(i,k:k+1), w(i,j,k:k+1)
                 endif
              enddo
           enddo
           do k=km, 2, -1
              do i=is,ie
                 if ( w2(i,k) > w_max ) then
                    gz(i) = (w2(i,k)-w_max) * dp2(i,k)
                    w2(i,k  ) = w_max
                    w2(i,k-1) = w2(i,k-1) + gz(i)/dp2(i,k-1)
                    print*, ' W_LIMITER up: ', i,j,k, w2(i,k-1:k), w(i,j,k-1:k)
                 elseif ( w2(i,k) < w_min ) then
                    gz(i) = (w2(i,k)-w_min) * dp2(i,k)
                    w2(i,k  ) = w_min
                    w2(i,k-1) = w2(i,k-1) + gz(i)/dp2(i,k-1)
                    print*, ' W_LIMITER up: ', i,j,k, w2(i,k-1:k), w(i,j,k-1:k)
                 endif
              enddo
           enddo
           do i=is,ie
              if (w2(i,1) > w_max*2. ) then
                 w2(i,1) = w_max*2 ! sink out of the top of the domain
                 print*, ' W_LIMITER top limited: ', i,j,1, w2(i,1), w(i,j,1)
              elseif (w2(i,1) < w_min*2. ) then
                 w2(i,1) = w_min*2.
                 print*, ' W_LIMITER top limited: ', i,j,1, w2(i,1), w(i,j,1)
              endif
           enddo
           do k=1,km
              do i=is,ie
                 w(i,j,k) = w2(i,k)
              enddo
           enddo
        endif
   endif

!----------
! Update pk
!----------
   do k=1,km+1
      do i=is,ie
         pk(i,j,k) = pk2(i,k)
      enddo
   enddo

!----------------
   if ( do_omega ) then
!                                Start do_omega
!                                Copy omega field to pe3
      do i=is,ie
         pe3(i,1) = 0.
      enddo
      do k=2,km+1
         do i=is,ie
            pe3(i,k) = omga(i,j,k-1)
         enddo
      enddo
   endif

   do k=1,km+1
      do i=is,ie
          pe0(i,k)   = peln(i,k,j)
         peln(i,k,j) =  pn2(i,k)
      enddo
   enddo

!------------
!
! Compute pkz
!
!< pk is pe**kappa(=rgas/cp_air), but pkz=plyr**kappa(=r/cp)
!------------
if ( .not. remap_te ) then

 if ( hydrostatic ) then
      do k=1,km
         do i=is,ie
            pkz(i,j,k) = (pk2(i,k+1)-pk2(i,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
#ifdef MULTI_GASES
            pkz(i,j,k) = exp(virqd(q(i,j,k,1:num_gas))/vicpqd(q(i,j,k,1:num_gas))*log(pkz(i,j,k)))
#endif
         enddo
      enddo
      
 else                !NH Note: pt at this stage is T_v or T_m
                     
         do k=1,km
#ifdef MOIST_CAPPA

            call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                          ice_wat, snowwat, graupel, hailwat, q, gz, cvm)
            do i=is,ie
               q_con(i,j,k) = gz(i)
#ifdef MULTI_GASES
               cappa(i,j,k) = rdgas / ( rdgas + cvm(i)/virq(q(i,j,k,1:num_gas)) )
#else
               cappa(i,j,k) = rdgas / ( rdgas + cvm(i)/(1.+r_vir*q(i,j,k,sphum)) )
#endif
               pkz(i,j,k) = exp(cappa(i,j,k)*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
	       pkzdry(i,k) =pkz(i,j,k)
	       ptdry(i,k) = pt(i,j,k)/pkzdry(i,k) 
            enddo
	    
#else
!                       dry cases

         if ( kord_tm < 0 ) then
           do i=is,ie
#ifdef MULTI_GASES
              pkz(i,j,k) = exp(akap*virqd(q(i,j,k,1:num_gas))/vicpqd(q(i,j,k,1:num_gas))*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
#else
              pkz(i,j,k) = exp(akap*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
#endif
           enddo
	   
         else               !kord_tm > 0  PT-?
	 
           do i=is,ie
#ifdef MULTI_GASES
              pkz(i,j,k) = exp(k1k*virqd(q(i,j,k,1:num_gas))/vicvqd(q(i,j,k,1:num_gas))*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
#else
              pkz(i,j,k) = exp(k1k*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
#endif
           enddo
         endif	
!================================
!    Need Tv for energy calculations 
!================================	    
         if ( kord_tm > 0 ) then                       !Need Tv for energy calculations 
            do k=1,km
               do i=is,ie
                  pt(i,j,k) = pt(i,j,k)*pkz(i,j,k)     !Need Tv for energy calculations
               enddo
            enddo
         endif	
	    
!        Bug  ???         if ( last_step .and. (.not.adiabatic) ) then
!                         do i=is,ie
!                          pt(i,j,k) = pt(i,j,k)*pkz(i,j,k)
!                       enddo
!       Bug              endif

	 
#endif
      enddo       !k-index
   endif          ! HYD or NHYD
   
 endif            ! endif not remap_te	 

 if ( last_step ) then
! Interpolate omega/pe3 (defined at pe0) to remapped cell center (dp2)
   if ( do_omega ) then
   do k=1,km
      do i=is,ie
         dp2(i,k) = 0.5*(peln(i,k,j) + peln(i,k+1,j))
      enddo
   enddo
   do i=is,ie
       k_next = 1
       do n=1,km
          kp = k_next
          do k=kp,km
             if( dp2(i,n) <= pe0(i,k+1) .and. dp2(i,n) >= pe0(i,k) ) then
                 omga(i,j,n) = pe3(i,k)  +  (pe3(i,k+1) - pe3(i,k)) *    &
                       (dp2(i,n)-pe0(i,k)) / (pe0(i,k+1)-pe0(i,k) )
                 k_next = k
                 exit
             endif
          enddo
       enddo
   enddo
    endif     ! end do_omega
   endif     ! end last_step
   
  endif !(j < je+1)

      do i=is,ie+1
         pe0(i,1) = pe(i,1,j)
      enddo
!------
! map u
!------
      do k=2,km+1
         do i=is,ie
            pe0(i,k) = 0.5*(pe(i,k,j-1)+pe1(i,k))
         enddo
      enddo

      do k=1,km+1
         bkh = 0.5*bk(k)
         do i=is,ie
            pe3(i,k) = ak(k) + bkh*(pe(i,km+1,j-1)+pe1(i,km+1))
	    pe3u(i,k) = pe3(i,k)
         enddo
      enddo
     if (.not. mdz_uv) then 
      call map1_ppm( km, pe0(is:ie,:),   u,   gz,   &
                     km, pe3(is:ie,:),   u,               &
                     is, ie, j, isd, ied, jsd, jed+1, -1, kord_mt)
     else
       call map1_ppm_dpwind( km, pe0(is:ie,:),   u,   gz,       &
                     km, pe3u(is:ie,:),   u, dpu(is:ie,:),               &
                     is, ie, j, isd, ied, jsd, jed+1, -1, kord_mt)		     
     endif		     

   if (j < je+1) then
!------
! map v
!------
       do i=is,ie+1
          pe3(i,1) = ak(1)
       enddo

       do k=2,km+1
          bkh = 0.5*bk(k)
          do i=is,ie+1
             pe0(i,k) =         0.5*(pe(i-1,k,   j)+pe(i,k,   j))
             pe3(i,k) = ak(k) + bkh*(pe(i-1,km+1,j)+pe(i,km+1,j))
          enddo
       enddo
     if (.not. mdz_uv) then 
       call map1_ppm (km, pe0,  v, gz,    &
                      km, pe3,  v, is, ie+1,    &
                      j, isd, ied+1, jsd, jed, -1, kord_mt)
     else
        call map1_ppm_dpwind (km, pe0,  v, gz,    &
                      km, pe3,  v, dpv, is, ie+1,    &
                      j, isd, ied+1, jsd, jed, -1, kord_mt)  		     
     endif		      
		      
!========================================		      
! 4a) update Tv and pkz from total energy 
!      (if remapping total energy)
!========================================
    if ( remap_te ) then
         do i=is,ie
            phis(i,km+1) = hs(i,j)
         enddo
         ! calculate Tv from TE
         if ( hydrostatic ) then
            do k=km,1,-1
               do i=is,ie
                  tpe = te(i,j,k) - phis(i,k+1) - 0.25*gridstruct%rsin2(i,j)*(    &
                       u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                       (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j) )
                  dlnp = rg*(peln(i,k+1,j) - peln(i,k,j))
                  pt(i,j,k)= tpe / (cp - pe2(i,k)*dlnp/delp(i,j,k))
		  
                  pkz(i,j,k) = (pk2(i,k+1)-pk2(i,k))/(akap*(peln(i,k+1,j)-peln(i,k,j)))
                  phis(i,k) = phis(i,k+1) + dlnp*pt(i,j,k)
		  
!    To do list add multi_gases -option	to TPE: pkz-akap
	  
               enddo
            enddo           ! end k-loop
        else                ! NH-case
            do k=km,1,-1
#ifdef MOIST_CAPPA
               call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                    ice_wat, snowwat, graupel, hailwat, q, gz(is:ie), cvm(is:ie))
               do i=is,ie
                  q_con(i,j,k) = gz(i)
#ifdef MULTI_GASES
                  cappa(i,j,k) = rdgas / ( rdgas + cvm(i)/virq(q(i,j,k,1:num_gas)) )
#else
                  cappa(i,j,k) = rdgas / ( rdgas + cvm(i)/(1.+r_vir*q(i,j,k,sphum)) )
#endif	       
               enddo
#endif
               do i=is,ie
                  phis(i,k) = phis(i,k+1) - delz(i,j,k)*grav
                  tpe = te(i,j,k) - 0.5*(phis(i,k)+phis(i,k+1)) - 0.5*w(i,j,k)**2 - 0.25*gridstruct%rsin2(i,j)*(    &
                       u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                       (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j) )
#ifdef MOIST_CAPPA
                  pt(i,j,k)= tpe / cvm(i)*(1.+r_vir*q(i,j,k,sphum))*(1.-gz(i))
                  pkz(i,j,k) = exp(cappa(i,j,k)*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
#else
                  pt(i,j,k)= tpe / cv_air *(1.+r_vir*q(i,j,k,sphum))
                  pkz(i,j,k) = exp(akap*log(rrg*delp(i,j,k)/delz(i,j,k)*pt(i,j,k)))
#endif
               enddo

            enddo           ! end k-loop
         endif		    ! NH or NYD  
      endif	            ! remap_te			      
		      
       if ( mdz_tk) then
         call  get_moldiff(mdt, vumol, ktmol, dfmol, rhomol, wgrav,   &
	            q, pt, w, u, v, dp2, pe2, dpu(is:ie,:), dpv, pe3u(is:ie,:), pe3, cappa,   &
                              j, je, is, ie, isd, ied, jsd, jed, km, nq,  &
	                      ind_h2o, ind_o2, ind_o3p, ind_o3 ) 
       endif
       			      
   endif ! (j < je+1)
!-----------------------------------------------------------------------------------------   
! j=je+1: we need update U(:j+1,:) by molecular diffusion   
!-----------------------------------------------------------------------------------------
     do k=1,km
        do i=is,ie
           pe4(i,j,k) = pe2(i,k+1)
        enddo
     enddo
     
       if ( mdz_tk .and. (j == je+1)) then
!
!  j=je+1: apply v-molecular viscosity only for u(:je+1,:) ONLY-!!!
!
        do i=is,ie+1
          pe3(i,1) = ak(1)
	  dpv(i,1) = pe3(i,1) 
       enddo

       do k=2,km+1
          bkh = 0.5*bk(k)
          do i=is,ie+1
            pe3(i,k) = ak(k) + bkh*(pe(i-1,km+1,j)+pe(i,km+1,j))
          enddo
       enddo  
      do k=2,km
       do i=is,ie+1      
        dpv(i,k)  = pe3(i,k+1)-pe3(i,k)
       enddo	
      enddo           
         call  get_moldiff(mdt, vumol, ktmol, dfmol, rhomol, wgrav,   &
	  q, pt, w, u, v, dpU(is:ie,:), pe3U(is:ie,:), dpu(is:ie,:), dpv, pe3u(is:ie,:), pe3, cappa,   &
                              j, je, is, ie, isd, ied, jsd, jed, km, nq,  &
	                      ind_h2o, ind_o2, ind_o3p, ind_o3 ) 
       endif
1000  continue  !j-loop

!===================
!6) Energy fixer
!===================

!$OMP parallel default(none) shared(is,ie,js,je,km,kmp,ptop,u,v,pe,ua,va,isd,ied,jsd,jed,kord_mt, &
!$OMP                               te_2d,te,delp,hydrostatic,hs,rg,pt,peln, adiabatic,        &
!$OMP                               cp,delz,nwat,rainwat,liq_wat,ice_wat,snowwat,              &
!$OMP                               graupel,hailwat,q_con,r_vir,sphum,w,pk,pkz,last_step,consv,        &
!$OMP                               do_adiabatic_init,zsum1,zsum0,te0_2d,domain,               &
!$OMP                               ng,gridstruct,E_Flux,pdt,dtmp,reproduce_sum,q,             &
!$OMP                               mdt,cld_amt,cappa,dtdt,out_dt,rrg,akap,do_sat_adj,         &
!$OMP                               fast_mp_consv,kord_tm, pe4,npx,npy, ccn_cm3,               &
!$OMP                               u_dt,v_dt,c2l_ord,bd,dp0,ps,cdata,GFDL_interstitial)        &
!$OMP                        shared(ccpp_suite)                                                &
#ifdef MULTI_GASES
!$OMP                        shared(num_gas)                                                   &
#endif
!$OMP                       private(q2,pe0,pe1,pe2,pe3,qv,cvm,gz,gsize,phis,kdelz,dp2,t0, ierr)

!$OMP do
  do k=2,km
     do j=js,je
        do i=is,ie
           pe(i,k,j) = pe4(i,j,k-1)
        enddo
     enddo
  enddo

  dtmp = 0.
if( last_step .and. (.not.do_adiabatic_init)  ) then

if ( consv > consv_min ) then

!$OMP do
    do j=js,je
    
       if ( hydrostatic ) then
            do i=is,ie
               gz(i) = hs(i,j)
               do k=1,km
                  gz(i) = gz(i) + rg*pt(i,j,k)*(peln(i,k+1,j)-peln(i,k,j))
               enddo
            enddo
            do i=is,ie
               te_2d(i,j) = pe(i,km+1,j)*hs(i,j) - pe(i,1,j)*gz(i)
            enddo

            do k=1,km
            do i=is,ie
               te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cp*pt(i,j,k) +   &
                            0.25*gridstruct%rsin2(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                                                        v(i,j,k)**2+v(i+1,j,k)**2 -  &
                           (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j)))
            enddo
            enddo
       else
           do i=is,ie
              te_2d(i,j) = 0.
              phis(i,km+1) = hs(i,j)
           enddo
           do k=km,1,-1
              do i=is,ie
                 phis(i,k) = phis(i,k+1) - grav*delz(i,j,k)
              enddo
           enddo

           do k=1,km
#ifdef MOIST_CAPPA
              call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                            ice_wat, snowwat, graupel, hailwat, q, gz, cvm)
              do i=is,ie
! KE using 3D winds:
              q_con(i,j,k) = gz(i)
#ifdef MULTI_GASES
              te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cvm(i)*pt(i,j,k)/virq(q(i,j,k,1:num_gas)) + &
#else
              te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cvm(i)*pt(i,j,k)/((1.+r_vir*q(i,j,k,sphum))*(1.-gz(i))) + &
#endif
                              0.5*(phis(i,k)+phis(i,k+1) + w(i,j,k)**2 + 0.5*gridstruct%rsin2(i,j)*( &
                              u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                             (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j))))
              enddo
#else
              do i=is,ie
#ifdef MULTI_GASES
                 te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cv_air*pt(i,j,k)/virq(q(i,j,k,1:num_gas)) + &
#else
                 te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cv_air*pt(i,j,k)/(1.+r_vir*q(i,j,k,sphum)) + &
#endif
                                 0.5*(phis(i,k)+phis(i,k+1) + w(i,j,k)**2 + 0.5*gridstruct%rsin2(i,j)*( &
                                 u(i,j,k)**2+u(i,j+1,k)**2 + v(i,j,k)**2+v(i+1,j,k)**2 -  &
                                (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*gridstruct%cosa_s(i,j))))
              enddo
#endif
           enddo   ! k-loop
       endif       ! end non-hydro

         do i=is,ie
            te_2d(i,j) = te0_2d(i,j) - te_2d(i,j)
            zsum1(i,j) = pkz(i,j,1)*delp(i,j,1)
         enddo
         do k=2,km
            do i=is,ie
               zsum1(i,j) = zsum1(i,j) + pkz(i,j,k)*delp(i,j,k)
            enddo
         enddo
         if ( hydrostatic ) then
            do i=is,ie
               zsum0(i,j) = ptop*(pk(i,j,1)-pk(i,j,km+1)) + zsum1(i,j)
            enddo
         endif

    enddo   ! j-loop

!$OMP single
      dtmp = consv*g_sum(domain, te_2d, is, ie, js, je, ng, gridstruct%area_64, 0, reproduce=.true.)
      E_Flux = dtmp / (grav*pdt*4.*pi*radius**2)    ! unit: W/m**2
                                                   ! Note pdt is "phys" time step
    if ( hydrostatic ) then
      dtmp = dtmp / (cp*    g_sum(domain, zsum0, is, ie, js, je, ng, gridstruct%area_64, 0, reproduce=.true.))
     else
       dtmp = dtmp / (cv_air*g_sum(domain, zsum1, is, ie, js, je, ng, gridstruct%area_64, 0, reproduce=.true.))
    endif
!$OMP end single

  elseif ( consv < -consv_min ) then

!$OMP do
      do j=js,je
         do i=is,ie
            zsum1(i,j) = pkz(i,j,1)*delp(i,j,1)
         enddo
         do k=2,km
            do i=is,ie
               zsum1(i,j) = zsum1(i,j) + pkz(i,j,k)*delp(i,j,k)
            enddo
         enddo
         if ( hydrostatic ) then
            do i=is,ie
               zsum0(i,j) = ptop*(pk(i,j,1)-pk(i,j,km+1)) + zsum1(i,j)
            enddo
         endif
      enddo

      E_Flux = consv
!$OMP single
      if ( hydrostatic ) then
           dtmp = E_flux*(grav*pdt*4.*pi*radius**2) /    &
                 (cp*g_sum(domain, zsum0,  is, ie, js, je, ng, gridstruct%area_64, 0, reproduce=.true.))
      else
           dtmp = E_flux*(grav*pdt*4.*pi*radius**2) /    &
                 (cv_air*g_sum(domain, zsum1,  is, ie, js, je, ng, gridstruct%area_64, 0, reproduce=.true.))
      endif
!$OMP end single
  endif        ! end consv check
endif          ! end last_step check

! Note: pt at this stage is T_v
! if ( (.not.do_adiabatic_init) .and. do_sat_adj ) then

  if ( do_sat_adj ) then
                                           call timing_on('sat_adj2')
    ! Call to CCPP fast_physics group
    if (cdata%initialized()) then
      call ccpp_physics_run(cdata, suite_name=trim(ccpp_suite), group_name='fast_physics', ierr=ierr)
      if (ierr/=0) then
        call mpp_error(NOTE, trim(cdata%errmsg))
        call mpp_error(FATAL, "Call to ccpp_physics_run for group 'fast_physics' failed")
      endif
    else
      call mpp_error (FATAL, 'Lagrangian_to_Eulerian: can not call CCPP fast physics because CCPP not initialized')
    endif
                                           call timing_off('sat_adj2')
  endif   ! do_sat_adj

       dtmp = 0.
  if ( last_step ) then
                               ! Output temperature Tk =TV/(1.+r_vir*q) if last_step
!$OMP do
        do k=1,km
           do j=js,je
	   
           if (hydrostatic) then !This is re-factored from AM4 so answers may be different
              do i=is,ie
                 pt(i,j,k) = (pt(i,j,k)+dtmp/cp*pkz(i,j,k)) / (1.+r_vir*q(i,j,k,sphum))
              enddo
           else	   
#ifdef USE_COND	      
                 call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                               ice_wat, snowwat, graupel, hailwat, q, gz, cvm)
                 do i=is,ie
#ifdef MULTI_GASES
                    pt(i,j,k) = (pt(i,j,k)+dtmp*pkz(i,j,k)) / virq(q(i,j,k,1:num_gas))
#else
                    pt(i,j,k) = (pt(i,j,k)+dtmp*pkz(i,j,k)) / ((1.+r_vir*q(i,j,k,sphum))*(1.-gz(i)))
#endif
                 enddo
           endif
	      
#else
              if ( .not. adiabatic ) then
                do i=is,ie
#ifdef MULTI_GASES
                    pt(i,j,k) = (pt(i,j,k)+dtmp*pkz(i,j,k)) / virq(q(i,j,k,1:num_gas))
#else
                    pt(i,j,k) = (pt(i,j,k)+dtmp*pkz(i,j,k)) / (1.+r_vir*q(i,j,k,sphum))
#endif
                 enddo
              endif
#endif
           enddo   ! j-loop
        enddo  ! k-loop
	
  else  ! last_step
!======================================= not last step======  
!$OMP do
       do k=1,km
          do j=js,je
             do i=is,ie
                pt(i,j,k) = pt(i,j,k)/pkz(i,j,k)
             enddo
          enddo
       enddo

  endif   ! "for" not last_step
!$OMP end parallel

  end associate ccpp_associate

 end subroutine Lagrangian_to_Eulerian


!>@brief The subroutine 'compute_total_energy' performs the FV3-consistent computation of the global total energy.
!>@details It includes the potential, internal (latent and sensible heat), kinetic terms.
 subroutine compute_total_energy(is, ie, js, je, isd, ied, jsd, jed, km,       &
                                 u, v, w, delz, pt, delp, q, qc, pe, peln, hs, &
                                 rsin2_l, cosa_s_l, &
                                 r_vir,  cp, rg, hlv, te_2d, ua, va, teq, &
                                 moist_phys, nwat, sphum, liq_wat, rainwat, ice_wat, snowwat, graupel, hailwat, hydrostatic, id_te)
!------------------------------------------------------
! Compute vertically integrated total energy per column
!------------------------------------------------------
! !INPUT PARAMETERS:
   integer,  intent(in):: km, is, ie, js, je, isd, ied, jsd, jed, id_te
   integer,  intent(in):: sphum, liq_wat, ice_wat, rainwat, snowwat, graupel, hailwat, nwat
   real, intent(inout), dimension(isd:ied,jsd:jed,km):: ua, va
   real, intent(in), dimension(isd:ied,jsd:jed,km):: pt, delp
   real, intent(in), dimension(isd:ied,jsd:jed,km,*):: q
   real, intent(in), dimension(isd:ied,jsd:jed,km):: qc
   real, intent(inout)::  u(isd:ied,  jsd:jed+1,km)
   real, intent(inout)::  v(isd:ied+1,jsd:jed,  km)
   real, intent(in)::  w(isd:,jsd:,1:)   !< vertical velocity (m/s)
   real, intent(in):: delz(is:,js:,1:)
   real, intent(in):: hs(isd:ied,jsd:jed)  !< surface geopotential
   real, intent(in)::   pe(is-1:ie+1,km+1,js-1:je+1) !< pressure at layer edges
   real, intent(in):: peln(is:ie,km+1,js:je)  !< log(pe)
   real, intent(in):: cp, rg, r_vir, hlv
   real, intent(in) :: rsin2_l(isd:ied, jsd:jed)
   real, intent(in) :: cosa_s_l(isd:ied, jsd:jed)
   logical, intent(in):: moist_phys, hydrostatic
!! Output:
   real, intent(out):: te_2d(is:ie,js:je)   !< vertically integrated TE
   real, intent(out)::   teq(is:ie,js:je)   !< Moist TE
!! Local
   real, dimension(is:ie,km):: tv
   real  phiz(is:ie,km+1)
   real  cvm(is:ie), qd(is:ie)
   integer i, j, k

!----------------------
! Output lat-lon winds:
!----------------------
!  call cubed_to_latlon(u, v, ua, va, dx, dy, rdxa, rdya, km, flagstruct%c2l_ord)

!$OMP parallel do default(none) shared(is,ie,js,je,isd,ied,jsd,jed,km,hydrostatic,hs,pt,qc,rg,peln,te_2d, &
!$OMP                                  pe,delp,cp,rsin2_l,u,v,cosa_s_l,delz,moist_phys,w, &
#ifdef MULTI_GASES
!$OMP                                  num_gas,                                           &
#endif
!$OMP                                  q,nwat,liq_wat,rainwat,ice_wat,snowwat,graupel,hailwat,sphum)   &
!$OMP                          private(phiz, tv, cvm, qd)
  do j=js,je

     if ( hydrostatic ) then

        do i=is,ie
           phiz(i,km+1) = hs(i,j)
        enddo
        do k=km,1,-1
           do i=is,ie
                tv(i,k) = pt(i,j,k)*(1.+qc(i,j,k))
              phiz(i,k) = phiz(i,k+1) + rg*tv(i,k)*(peln(i,k+1,j)-peln(i,k,j))
           enddo
        enddo

        do i=is,ie
           te_2d(i,j) = pe(i,km+1,j)*phiz(i,km+1) - pe(i,1,j)*phiz(i,1)
        enddo

        do k=1,km
           do i=is,ie
              te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*(cp*tv(i,k) +            &
                           0.25*rsin2_l(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +      &
                                            v(i,j,k)**2+v(i+1,j,k)**2 -      &
                       (u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s_l(i,j)))
           enddo
        enddo

     else
!-----------------
! Non-hydrostatic:
!-----------------
     do i=is,ie
        phiz(i,km+1) = hs(i,j)
        do k=km,1,-1
           phiz(i,k) = phiz(i,k+1) - grav*delz(i,j,k)
        enddo
     enddo
     do i=is,ie
        te_2d(i,j) = 0.
     enddo
     if ( moist_phys ) then
     do k=1,km
#ifdef MOIST_CAPPA
        call moist_cv(is,ie,isd,ied,jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                      ice_wat, snowwat, graupel, hailwat, q, qd, cvm)
#endif
        do i=is,ie
#ifdef MOIST_CAPPA
           te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cvm(i)*pt(i,j,k) +  &
#else
#ifdef MULTI_GASES
           te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cv_air*vicvqd(q(i,j,k,1:num_gas) )*pt(i,j,k) +  &
#else
           te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cv_air*pt(i,j,k) +  &
#endif
#endif
                        0.5*(phiz(i,k)+phiz(i,k+1)+w(i,j,k)**2+0.5*rsin2_l(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                        v(i,j,k)**2+v(i+1,j,k)**2-(u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s_l(i,j))))
        enddo
     enddo
     else
       do k=1,km
          do i=is,ie
#ifdef MULTI_GASES
             te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cv_air*vicvqd(q(i,j,k,1:num_gas))*pt(i,j,k) +  &
#else
             te_2d(i,j) = te_2d(i,j) + delp(i,j,k)*( cv_air*pt(i,j,k) +  &
#endif
                          0.5*(phiz(i,k)+phiz(i,k+1)+w(i,j,k)**2+0.5*rsin2_l(i,j)*(u(i,j,k)**2+u(i,j+1,k)**2 +  &
                          v(i,j,k)**2+v(i+1,j,k)**2-(u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))*cosa_s_l(i,j))))
          enddo
       enddo
     endif
     endif
  enddo

!-------------------------------------
! Diganostics computation for moist TE
!-------------------------------------
  if( id_te>0 ) then
!$OMP parallel do default(none) shared(is,ie,js,je,teq,te_2d,moist_phys,km,hlv,sphum,q,delp)
      do j=js,je
         do i=is,ie
            teq(i,j) = te_2d(i,j)
         enddo
         if ( moist_phys ) then
           do k=1,km
              do i=is,ie
                 teq(i,j) = teq(i,j) + hlv*q(i,j,k,sphum)*delp(i,j,k)
              enddo
           enddo
         endif
      enddo
   endif

  end subroutine compute_total_energy

  subroutine pkez(km, ifirst, ilast, jfirst, jlast, j, &
                  pe, pk, akap, peln, pkz, ptop)

! INPUT PARAMETERS:
   integer, intent(in):: km, j
   integer, intent(in):: ifirst, ilast        !< Latitude strip
   integer, intent(in):: jfirst, jlast        !< Latitude strip
   real, intent(in):: akap
   real, intent(in):: pe(ifirst-1:ilast+1,km+1,jfirst-1:jlast+1)
   real, intent(in):: pk(ifirst:ilast,jfirst:jlast,km+1)
   real, intent(IN):: ptop
! OUTPUT
   real, intent(out):: pkz(ifirst:ilast,jfirst:jlast,km)
   real, intent(inout):: peln(ifirst:ilast, km+1, jfirst:jlast)   !< log (pe)
! Local
   real pk2(ifirst:ilast, km+1)
   real pek
   real lnp
   real ak1
   integer i, k

   ak1 = (akap + 1.) / akap

        pek = pk(ifirst,j,1)
        do i=ifirst, ilast
           pk2(i,1) = pek
        enddo

        do k=2,km+1
           do i=ifirst, ilast
!             peln(i,k,j) =  log(pe(i,k,j))
              pk2(i,k) =  pk(i,j,k)
           enddo
        enddo

!---- GFDL modification
       if( ptop < ptop_min ) then
           do i=ifirst, ilast
               peln(i,1,j) = peln(i,2,j) - ak1
           enddo
       else
           lnp = log( ptop )
           do i=ifirst, ilast
              peln(i,1,j) = lnp
           enddo
       endif
!---- GFDL modification

       do k=1,km
          do i=ifirst, ilast
             pkz(i,j,k) = (pk2(i,k+1) - pk2(i,k) )  /  &
                          (akap*(peln(i,k+1,j) - peln(i,k,j)) )
          enddo
       enddo

 end subroutine pkez



 subroutine remap_z(km, pe1, q1, kn, pe2, q2, i1, i2, iv, kord)

! INPUT PARAMETERS:
      integer, intent(in) :: i1                !< Starting longitude
      integer, intent(in) :: i2                !< Finishing longitude
      integer, intent(in) :: kord              !< Method order
      integer, intent(in) :: km                !< Original vertical dimension
      integer, intent(in) :: kn                !< Target vertical dimension
      integer, intent(in) :: iv

      real, intent(in) ::  pe1(i1:i2,km+1)     !< height at layer edges
                                               !! (from model top to bottom surface)
      real, intent(in) ::  pe2(i1:i2,kn+1)     !< hieght at layer edges
                                               !! (from model top to bottom surface)
      real, intent(in) ::  q1(i1:i2,km)        !< Field input

! INPUT/OUTPUT PARAMETERS:
      real, intent(inout)::  q2(i1:i2,kn)      !< Field output

! LOCAL VARIABLES:
      real   qs(i1:i2)
      real  dp1(  i1:i2,km)
      real   q4(4,i1:i2,km)
      real   pl, pr, qsum, delp, esl
      integer i, k, l, m, k0

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)      ! negative
            q4(1,i,k) = q1(i,k)
         enddo
      enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) <= pe1(i,l) .and. pe2(i,k) >= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) >= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) < pe1(i,m+1) ) then
! Whole layer..
                    qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                    delp = pe2(i,k+1)-pe1(i,m)
                    esl = delp / dp1(i,m)
                    qsum = qsum + delp*(q4(2,i,m)+0.5*esl*               &
                         (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                    k0 = m
                 goto 123
                 endif
              enddo
              goto 123
           endif
      endif
100   continue
123   q2(i,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
1000  continue

 end subroutine remap_z

 subroutine map_scalar( km,   pe1,    q1,   qs,           &
                        kn,   pe2,    q2,   i1, i2,       &
                         j,  ibeg, iend, jbeg, jend, iv,  kord, q_min)
! iv=1
 integer, intent(in) :: i1                !< Starting longitude
 integer, intent(in) :: i2                !< Finishing longitude
 integer, intent(in) :: iv                !< Mode: 0 == constituents 1 == temp 2 == remap temp with cs scheme
 integer, intent(in) :: kord              !< Method order
 integer, intent(in) :: j                 !< Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                !< Original vertical dimension
 integer, intent(in) :: kn                !< Target vertical dimension
 real, intent(in) ::   qs(i1:i2)       !< bottom BC
 real, intent(in) ::  pe1(i1:i2,km+1)  !< pressure at layer edges
                                       !! (from model top to bottom surface)
                                       !! in the original vertical coordinate
 real, intent(in) ::  pe2(i1:i2,kn+1)  !< pressure at layer edges
                                       !! (from model top to bottom surface)
                                       !! in the new vertical coordinate
 real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) !< Field input
! !INPUT/OUTPUT PARAMETERS:
 real, intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) !< Field output
 real, intent(in):: q_min

! DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
! LOCAL VARIABLES:
   real    dp1(i1:i2,km)
   real   q4(4,i1:i2,km)
   real    pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q1(i,j,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call scalar_profile( qs, q4, dp1, km, i1, i2, iv, kord, q_min )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

  do i=i1,i2
     k0 = 1
     do 555 k=1,kn
      do l=k0,km
! locate the top edge: pe2(i,k)
      if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if( pe2(i,k+1) <= pe1(i,l+1) ) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
         else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if( pe2(i,k+1) > pe1(i,m+1) ) then
! Whole layer
                     qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                     esl = dp / dp1(i,m)
                     qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                           (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                     k0 = m
                     goto 123
                 endif
              enddo
              goto 123
         endif
      endif
      enddo
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
  enddo

 end subroutine map_scalar


 subroutine map1_ppm( km,   pe1,    q1,   qs,           &
                      kn,   pe2,    q2,   i1, i2,       &
                      j,    ibeg, iend, jbeg, jend, iv,  kord)
 integer, intent(in) :: i1                !< Starting longitude
 integer, intent(in) :: i2                !< Finishing longitude
 integer, intent(in) :: iv                !< Mode: 0 == constituents 1 == ??? 2 == remap temp with cs scheme
 integer, intent(in) :: kord              !< Method order
 integer, intent(in) :: j                 !< Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                !< Original vertical dimension
 integer, intent(in) :: kn                !< Target vertical dimension
 real, intent(in) ::   qs(i1:i2)       !< bottom BC
 real, intent(in) ::  pe1(i1:i2,km+1)  !< pressure at layer edges
                                       !! (from model top to bottom surface)
                                       !! in the original vertical coordinate
 real, intent(in) ::  pe2(i1:i2,kn+1)  !< pressure at layer edges
                                       !! (from model top to bottom surface)
                                       !! in the new vertical coordinate
 real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) !< Field input
! !INPUT/OUTPUT PARAMETERS:
 real, intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) !< Field output

! DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate

! LOCAL VARIABLES:
   real    dp1(i1:i2,km)
   real   q4(4,i1:i2,km)
   real    pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q1(i,j,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

  do i=i1,i2
     k0 = 1
     do 555 k=1,kn
      do l=k0,km
! locate the top edge: pe2(i,k)
      if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if( pe2(i,k+1) <= pe1(i,l+1) ) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
         else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if( pe2(i,k+1) > pe1(i,m+1) ) then
! Whole layer
                     qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                     esl = dp / dp1(i,m)
                     qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                           (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                     k0 = m
                     goto 123
                 endif
              enddo
              goto 123
         endif
      endif
      enddo
123   q2(i,j,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
555   continue
  enddo

 end subroutine map1_ppm


 subroutine mapn_tracer(nq, km, pe1, pe2, q1, dp2, kord, j,     &
                        i1, i2, isd, ied, jsd, jed, q_min, fill)
! INPUT PARAMETERS:
      integer, intent(in):: km                !< vertical dimension
      integer, intent(in):: j, nq, i1, i2
      integer, intent(in):: isd, ied, jsd, jed
      integer, intent(in):: kord(nq)
      real, intent(in)::  pe1(i1:i2,km+1)     !< pressure at layer edges
                                              !! (from model top to bottom surface)
                                              !! in the original vertical coordinate
      real, intent(in)::  pe2(i1:i2,km+1)     !< pressure at layer edges
                                              !! (from model top to bottom surface)
                                              !! in the new vertical coordinate
      real, intent(in)::  dp2(i1:i2,km)
      real, intent(in)::  q_min
      logical, intent(in):: fill
      real, intent(inout):: q1(isd:ied,jsd:jed,km,nq) ! Field input
! LOCAL VARIABLES:
      real:: q4(4,i1:i2,km,nq)
      real:: q2(i1:i2,km,nq) !< Field output
      real:: qsum(nq)
      real:: dp1(i1:i2,km)
      real:: qs(i1:i2)
      real:: pl, pr, dp, esl, fac1, fac2
      integer:: i, k, l, m, k0, iq

      do k=1,km
         do i=i1,i2
            dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         enddo
      enddo

      do iq=1,nq
         do k=1,km
            do i=i1,i2
               q4(1,i,k,iq) = q1(i,j,k,iq)
            enddo
         enddo
         call scalar_profile( qs, q4(1,i1,1,iq), dp1, km, i1, i2, 0, kord(iq), q_min )
      enddo

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,km
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            fac1 = pr + pl
            fac2 = r3*(pr*fac1 + pl*pl)
            fac1 = 0.5*fac1
            do iq=1,nq
               q2(i,k,iq) = q4(2,i,l,iq) + (q4(4,i,l,iq)+q4(3,i,l,iq)-q4(2,i,l,iq))*fac1  &
                                         -  q4(4,i,l,iq)*fac2
            enddo
            k0 = l
            goto 555
          else
! Fractional area...
            dp = pe1(i,l+1) - pe2(i,k)
            fac1 = 1. + pl
            fac2 = r3*(1.+pl*fac1)
            fac1 = 0.5*fac1
            do iq=1,nq
               qsum(iq) = dp*(q4(2,i,l,iq) + (q4(4,i,l,iq)+   &
                              q4(3,i,l,iq) - q4(2,i,l,iq))*fac1 - q4(4,i,l,iq)*fac2)
            enddo
            do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
               if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                  do iq=1,nq
                     qsum(iq) = qsum(iq) + dp1(i,m)*q4(1,i,m,iq)
                  enddo
               else
                  dp = pe2(i,k+1)-pe1(i,m)
                  esl = dp / dp1(i,m)
                  fac1 = 0.5*esl
                  fac2 = 1.-r23*esl
                  do iq=1,nq
                     qsum(iq) = qsum(iq) + dp*( q4(2,i,m,iq) + fac1*(         &
                                q4(3,i,m,iq)-q4(2,i,m,iq)+q4(4,i,m,iq)*fac2 ) )
                  enddo
                  k0 = m
                  goto 123
               endif
            enddo
            goto 123
          endif
      endif
100   continue
123   continue
      do iq=1,nq
         q2(i,k,iq) = qsum(iq) / dp2(i,k)
      enddo
555   continue
1000  continue

  if (fill) call fillz(i2-i1+1, km, nq, q2, dp2)

  do iq=1,nq
!    if (fill) call fillz(i2-i1+1, km, 1, q2(i1,1,iq), dp2)
     do k=1,km
        do i=i1,i2
           q1(i,j,k,iq) = q2(i,k,iq)
        enddo
     enddo
  enddo

 end subroutine mapn_tracer


 subroutine map1_q2(km,   pe1,   q1,            &
                    kn,   pe2,   q2,   dp2,     &
                    i1,   i2,    iv,   kord, j, &
                    ibeg, iend, jbeg, jend, q_min )


! INPUT PARAMETERS:
      integer, intent(in) :: j
      integer, intent(in) :: i1, i2
      integer, intent(in) :: ibeg, iend, jbeg, jend
      integer, intent(in) :: iv                !< Mode: 0 ==  constituents 1 == ???
      integer, intent(in) :: kord
      integer, intent(in) :: km                !< Original vertical dimension
      integer, intent(in) :: kn                !< Target vertical dimension

      real, intent(in) ::  pe1(i1:i2,km+1)     !< pressure at layer edges
                                               !! (from model top to bottom surface)
                                               !! in the original vertical coordinate
      real, intent(in) ::  pe2(i1:i2,kn+1)     !< pressure at layer edges
                                               !! (from model top to bottom surface)
                                               !! in the new vertical coordinate
      real, intent(in) ::  q1(ibeg:iend,jbeg:jend,km) ! Field input
      real, intent(in) ::  dp2(i1:i2,kn)
      real, intent(in) ::  q_min
! INPUT/OUTPUT PARAMETERS:
      real, intent(inout):: q2(i1:i2,kn) !< Field output
! LOCAL VARIABLES:
      real   qs(i1:i2)
      real   dp1(i1:i2,km)
      real   q4(4,i1:i2,km)
      real   pl, pr, qsum, dp, esl

      integer i, k, l, m, k0

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            q4(1,i,k) = q1(i,j,k)
         enddo
      enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  scalar_profile( qs, q4, dp1, km, i1, i2, iv, kord, q_min )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

! Mapping
      do 1000 i=i1,i2
         k0 = 1
      do 555 k=1,kn
      do 100 l=k0,km
! locate the top edge: pe2(i,k)
      if(pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1)) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
          else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                    qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                    esl = dp / dp1(i,m)
                   qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                       (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                   k0 = m
                   goto 123
                 endif
              enddo
              goto 123
          endif
      endif
100   continue
123   q2(i,k) = qsum / dp2(i,k)
555   continue
1000  continue

 end subroutine map1_q2



 subroutine remap_2d(km,   pe1,   q1,        &
                     kn,   pe2,   q2,        &
                     i1,   i2,    iv,   kord)
   integer, intent(in):: i1, i2
   integer, intent(in):: iv               !< Mode: 0 ==  constituents 1 ==others
   integer, intent(in):: kord
   integer, intent(in):: km               !< Original vertical dimension
   integer, intent(in):: kn               !< Target vertical dimension
   real, intent(in):: pe1(i1:i2,km+1)     !< pressure at layer edges
                                          !! (from model top to bottom surface)
                                          !! in the original vertical coordinate
   real, intent(in):: pe2(i1:i2,kn+1)     !< pressure at layer edges
                                          !! (from model top to bottom surface)
                                          !! in the new vertical coordinate
   real, intent(in) :: q1(i1:i2,km) !< Field input
   real, intent(out):: q2(i1:i2,kn) !< Field output
! !LOCAL VARIABLES:
   real   qs(i1:i2)
   real   dp1(i1:i2,km)
   real   q4(4,i1:i2,km)
   real   pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
          dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q1(i,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

   do i=i1,i2
      k0 = 1
      do 555 k=1,kn
#ifdef OLD_TOP_EDGE
         if( pe2(i,k+1) <= pe1(i,1) ) then
! Entire grid above old ptop
             q2(i,k) = q4(2,i,1)
         elseif( pe2(i,k) < pe1(i,1) .and. pe2(i,k+1)>pe1(i,1) ) then
! Partially above old ptop:
             q2(i,k) = q1(i,1)
#else
         if( pe2(i,k) <= pe1(i,1) ) then
! above old ptop:
             q2(i,k) = q1(i,1)
#endif
         else
           do l=k0,km
! locate the top edge: pe2(i,k)
           if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
               pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
               if(pe2(i,k+1) <= pe1(i,l+1)) then
! entire new grid is within the original grid
                  pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
                  q2(i,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                          *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
                  k0 = l
                  goto 555
               else
! Fractional area...
                 qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                         q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                        (r3*(1.+pl*(1.+pl))))
                 do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                    if(pe2(i,k+1) > pe1(i,m+1) ) then
                                                   ! Whole layer..
                       qsum = qsum + dp1(i,m)*q4(1,i,m)
                    else
                       dp = pe2(i,k+1)-pe1(i,m)
                      esl = dp / dp1(i,m)
                      qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                            (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                      k0 = m
                      goto 123
                    endif
                 enddo
                 goto 123
               endif
           endif
           enddo
123        q2(i,k) = qsum / ( pe2(i,k+1) - pe2(i,k) )
         endif
555   continue
   enddo

 end subroutine remap_2d


 subroutine scalar_profile(qs, a4, delp, km, i1, i2, iv, kord, qmin)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      !< vertical dimension
 integer, intent(in):: iv      !< iv =-1: winds iv = 0: positive definite scalars iv = 1: others
 integer, intent(in):: kord
 real, intent(in)   ::   qs(i1:i2)
 real, intent(in)   :: delp(i1:i2,km)     !< Layer pressure thickness
 real, intent(inout):: a4(4,i1:i2,km)     !< Interpolated values
 real, intent(in):: qmin
!-----------------------------------------------------------------------
 logical, dimension(i1:i2,km):: extm, ext5, ext6
 real  gam(i1:i2,km)
 real    q(i1:i2,km+1)
 real   d4(i1:i2)
 real   bet, a_bot, grat
 real   pmp_1, lac_1, pmp_2, lac_2, x0, x1
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km)
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo
 else
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo

  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif

!----- Perfectly linear scheme --------------------------------
 if ( abs(kord) > 16 ) then
  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo
  return
 endif
!----- Perfectly linear scheme --------------------------------
!------------------
! Apply constraints
!------------------
  im = i2 - i1 + 1

! Apply *large-scale* constraints
  do i=i1,i2
     q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
     q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1)
     enddo
  enddo

! Interior:
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
          else
! There exists a local min
               q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k))
          endif
        endif
     enddo
  enddo

! Bottom:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
     q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
  enddo

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
     enddo
  enddo

  do k=1,km
     if ( k==1 .or. k==km ) then
       do i=i1,i2
          extm(i,k) = (a4(2,i,k)-a4(1,i,k)) * (a4(3,i,k)-a4(1,i,k)) > 0.
       enddo
     else
       do i=i1,i2
          extm(i,k) = gam(i,k)*gam(i,k+1) < 0.
       enddo
     endif
     if ( abs(kord) > 9 ) then
       do i=i1,i2
          x0 = 2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k))
          x1 = abs(a4(2,i,k)-a4(3,i,k))
          a4(4,i,k) = 3.*x0
          ext5(i,k) = abs(x0) > x1
          ext6(i,k) = abs(a4(4,i,k)) > x1
       enddo
     endif
  enddo

!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  if ( iv==0 ) then
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  elseif ( iv==-1 ) then
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
  elseif ( iv==2 ) then
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  endif

  if ( iv/=2 ) then
     do i=i1,i2
        a4(4,i,1) = 3.*(2.*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
     enddo
     call cs_limiters(im, extm(i1,1), a4(1,i1,1), 1)
  endif

! k=2
   do i=i1,i2
      a4(4,i,2) = 3.*(2.*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
   enddo
   call cs_limiters(im, extm(i1,2), a4(1,i1,2), 2)

!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
  do k=3,km-2
     if ( abs(kord)<9 ) then
       do i=i1,i2
! Left  edges
          pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
          lac_1 = pmp_1 + 1.5*gam(i,k+2)
          a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                         max(a4(1,i,k), pmp_1, lac_1) )
! Right edges
          pmp_2 = a4(1,i,k) + 2.*gam(i,k)
          lac_2 = pmp_2 - 1.5*gam(i,k-1)
          a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                         max(a4(1,i,k), pmp_2, lac_2) )

          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==9 ) then
       do i=i1,i2
          if ( extm(i,k) .and. extm(i,k-1) ) then
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. extm(i,k+1) ) then
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. a4(1,i,k)<qmin ) then
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
            endif
          endif
       enddo
     elseif ( abs(kord)==10 ) then
       do i=i1,i2
          if( ext5(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
              elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                   pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                   lac_1 = pmp_1 + 1.5*gam(i,k+2)
                   a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                  max(a4(1,i,k), pmp_1, lac_1) )
                   pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                   lac_2 = pmp_2 - 1.5*gam(i,k-1)
                   a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                  max(a4(1,i,k), pmp_2, lac_2) )
              endif
          elseif( ext6(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
                  a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                 max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
                  a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                 max(a4(1,i,k), pmp_2, lac_2) )
              endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==12 ) then
       do i=i1,i2
          if( extm(i,k) ) then
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==13 ) then
       do i=i1,i2
          if( ext6(i,k) ) then
             if ( ext6(i,k-1) .and. ext6(i,k+1) ) then
! grid-scale 2-delta-z wave detected
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
             endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==14 ) then

       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==15 ) then   ! Revised abs(kord)=9 scheme
       do i=i1,i2
          if ( ext5(i,k) .and. ext5(i,k-1) ) then
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
          else if ( ext5(i,k) .and. ext5(i,k+1) ) then
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
          else if ( ext5(i,k) .and. a4(1,i,k)<qmin ) then
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
          elseif( ext6(i,k) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==16 ) then
       do i=i1,i2
          if( ext5(i,k) ) then
             if ( ext5(i,k-1) .or. ext5(i,k+1) ) then
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
             elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                 ! Left  edges
                 pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                 lac_1 = pmp_1 + 1.5*gam(i,k+2)
                 a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                     max(a4(1,i,k), pmp_1, lac_1) )
                 ! Right edges
                 pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                 lac_2 = pmp_2 - 1.5*gam(i,k-1)
                 a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                     max(a4(1,i,k), pmp_2, lac_2) )
             endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     else      ! kord = 11, 13
       do i=i1,i2
         if ( ext5(i,k) .and. (ext5(i,k-1).or.ext5(i,k+1).or.a4(1,i,k)<qmin) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
     endif

! Additional constraint to ensure positivity
     if ( iv==0 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  if ( iv==0 ) then
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  elseif ( iv .eq. -1 ) then
      do i=i1,i2
         if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
      enddo
  endif

  do k=km-1,km
     do i=i1,i2
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
     if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
     if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
  enddo

 end subroutine scalar_profile

 subroutine cs_profile(qs, a4, delp, km, i1, i2, iv, kord)
! Optimized vertical profile reconstruction:
! Latest: Apr 2008 S.-J. Lin, NOAA/GFDL
 integer, intent(in):: i1, i2
 integer, intent(in):: km      !< vertical dimension
 integer, intent(in):: iv      !< iv =-1: winds
                               !< iv = 0: positive definite scalars
                               !< iv = 1: others
 integer, intent(in):: kord
 real, intent(in)   ::   qs(i1:i2)
 real, intent(in)   :: delp(i1:i2,km)     !< layer pressure thickness
 real, intent(inout):: a4(4,i1:i2,km)     !< Interpolated values
!-----------------------------------------------------------------------
 logical, dimension(i1:i2,km):: extm, ext5, ext6
 real  gam(i1:i2,km)
 real    q(i1:i2,km+1)
 real   d4(i1:i2)
 real   bet, a_bot, grat
 real   pmp_1, lac_1, pmp_2, lac_2, x0, x1
 integer i, k, im

 if ( iv .eq. -2 ) then
      do i=i1,i2
         gam(i,2) = 0.5
           q(i,1) = 1.5*a4(1,i,1)
      enddo
      do k=2,km-1
         do i=i1, i2
                  grat = delp(i,k-1) / delp(i,k)
                   bet =  2. + grat + grat - gam(i,k)
                q(i,k) = (3.*(a4(1,i,k-1)+a4(1,i,k)) - q(i,k-1))/bet
            gam(i,k+1) = grat / bet
         enddo
      enddo
      do i=i1,i2
            grat = delp(i,km-1) / delp(i,km)
         q(i,km) = (3.*(a4(1,i,km-1)+a4(1,i,km)) - grat*qs(i) - q(i,km-1)) /  &
                   (2. + grat + grat - gam(i,km))
         q(i,km+1) = qs(i)
      enddo
      do k=km-1,1,-1
        do i=i1,i2
           q(i,k) = q(i,k) - gam(i,k+1)*q(i,k+1)
        enddo
      enddo
 else
  do i=i1,i2
         grat = delp(i,2) / delp(i,1)   ! grid ratio
          bet = grat*(grat+0.5)
       q(i,1) = ( (grat+grat)*(grat+1.)*a4(1,i,1) + a4(1,i,2) ) / bet
     gam(i,1) = ( 1. + grat*(grat+1.5) ) / bet
  enddo

  do k=2,km
     do i=i1,i2
           d4(i) = delp(i,k-1) / delp(i,k)
             bet =  2. + d4(i) + d4(i) - gam(i,k-1)
          q(i,k) = ( 3.*(a4(1,i,k-1)+d4(i)*a4(1,i,k)) - q(i,k-1) )/bet
        gam(i,k) = d4(i) / bet
     enddo
  enddo

  do i=i1,i2
         a_bot = 1. + d4(i)*(d4(i)+1.5)
     q(i,km+1) = (2.*d4(i)*(d4(i)+1.)*a4(1,i,km)+a4(1,i,km-1)-a_bot*q(i,km))  &
               / ( d4(i)*(d4(i)+0.5) - a_bot*gam(i,km) )
  enddo

  do k=km,1,-1
     do i=i1,i2
        q(i,k) = q(i,k) - gam(i,k)*q(i,k+1)
     enddo
  enddo
 endif
!----- Perfectly linear scheme --------------------------------
 if ( abs(kord) > 16 ) then
  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
  enddo
  return
 endif
!----- Perfectly linear scheme --------------------------------

!------------------
! Apply constraints
!------------------
  im = i2 - i1 + 1

! Apply *large-scale* constraints
  do i=i1,i2
     q(i,2) = min( q(i,2), max(a4(1,i,1), a4(1,i,2)) )
     q(i,2) = max( q(i,2), min(a4(1,i,1), a4(1,i,2)) )
  enddo

  do k=2,km
     do i=i1,i2
        gam(i,k) = a4(1,i,k) - a4(1,i,k-1)
     enddo
  enddo

! Interior:
  do k=3,km-1
     do i=i1,i2
        if ( gam(i,k-1)*gam(i,k+1)>0. ) then
! Apply large-scale constraint to ALL fields if not local max/min
             q(i,k) = min( q(i,k), max(a4(1,i,k-1),a4(1,i,k)) )
             q(i,k) = max( q(i,k), min(a4(1,i,k-1),a4(1,i,k)) )
        else
          if ( gam(i,k-1) > 0. ) then
! There exists a local max
               q(i,k) = max(q(i,k), min(a4(1,i,k-1),a4(1,i,k)))
          else
! There exists a local min
                 q(i,k) = min(q(i,k), max(a4(1,i,k-1),a4(1,i,k)))
               if ( iv==0 ) q(i,k) = max(0., q(i,k))
          endif
        endif
     enddo
  enddo

! Bottom:
  do i=i1,i2
     q(i,km) = min( q(i,km), max(a4(1,i,km-1), a4(1,i,km)) )
     q(i,km) = max( q(i,km), min(a4(1,i,km-1), a4(1,i,km)) )
  enddo

  do k=1,km
     do i=i1,i2
        a4(2,i,k) = q(i,k  )
        a4(3,i,k) = q(i,k+1)
     enddo
  enddo

  do k=1,km
     if ( k==1 .or. k==km ) then
       do i=i1,i2
          extm(i,k) = (a4(2,i,k)-a4(1,i,k)) * (a4(3,i,k)-a4(1,i,k)) > 0.
       enddo
     else
       do i=i1,i2
          extm(i,k) = gam(i,k)*gam(i,k+1) < 0.
       enddo
     endif
     if ( abs(kord) > 9 ) then
       do i=i1,i2
          x0 = 2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k))
          x1 = abs(a4(2,i,k)-a4(3,i,k))
          a4(4,i,k) = 3.*x0
          ext5(i,k) = abs(x0) > x1
          ext6(i,k) = abs(a4(4,i,k)) > x1
       enddo
     endif
  enddo

!---------------------------
! Apply subgrid constraints:
!---------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
! Top 2 and bottom 2 layers always use monotonic mapping

  if ( iv==0 ) then
     do i=i1,i2
        a4(2,i,1) = max(0., a4(2,i,1))
     enddo
  elseif ( iv==-1 ) then
      do i=i1,i2
         if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
      enddo
  elseif ( iv==2 ) then
     do i=i1,i2
        a4(2,i,1) = a4(1,i,1)
        a4(3,i,1) = a4(1,i,1)
        a4(4,i,1) = 0.
     enddo
  endif

  if ( iv/=2 ) then
     do i=i1,i2
        a4(4,i,1) = 3.*(2.*a4(1,i,1) - (a4(2,i,1)+a4(3,i,1)))
     enddo
     call cs_limiters(im, extm(i1,1), a4(1,i1,1), 1)
  endif

! k=2
   do i=i1,i2
      a4(4,i,2) = 3.*(2.*a4(1,i,2) - (a4(2,i,2)+a4(3,i,2)))
   enddo
   call cs_limiters(im, extm(i1,2), a4(1,i1,2), 2)

!-------------------------------------
! Huynh's 2nd constraint for interior:
!-------------------------------------
  do k=3,km-2
     if ( abs(kord)<9 ) then
       do i=i1,i2
! Left  edges
          pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
          lac_1 = pmp_1 + 1.5*gam(i,k+2)
          a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                         max(a4(1,i,k), pmp_1, lac_1) )
! Right edges
          pmp_2 = a4(1,i,k) + 2.*gam(i,k)
          lac_2 = pmp_2 - 1.5*gam(i,k-1)
          a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                         max(a4(1,i,k), pmp_2, lac_2) )

          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==9 ) then
       do i=i1,i2
          if ( extm(i,k) .and. extm(i,k-1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else if ( extm(i,k) .and. extm(i,k+1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
               a4(2,i,k) = a4(1,i,k)
               a4(3,i,k) = a4(1,i,k)
               a4(4,i,k) = 0.
          else
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==10 ) then
       do i=i1,i2
          if( ext5(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                   a4(2,i,k) = a4(1,i,k)
                   a4(3,i,k) = a4(1,i,k)
              elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                   pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                   lac_1 = pmp_1 + 1.5*gam(i,k+2)
                   a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                  max(a4(1,i,k), pmp_1, lac_1) )
                   pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                   lac_2 = pmp_2 - 1.5*gam(i,k-1)
                   a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                  max(a4(1,i,k), pmp_2, lac_2) )
              endif
          elseif( ext6(i,k) ) then
              if( ext5(i,k-1) .or. ext5(i,k+1) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
                  a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                                 max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
                  a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                                 max(a4(1,i,k), pmp_2, lac_2) )
              endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==12 ) then
       do i=i1,i2
          if( extm(i,k) ) then
! grid-scale 2-delta-z wave detected
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
          else        ! not a local extremum
            a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
! Check within the smooth region if subgrid profile is non-monotonic
            if( abs(a4(4,i,k)) > abs(a4(2,i,k)-a4(3,i,k)) ) then
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
              a4(4,i,k) = 6.*a4(1,i,k) - 3.*(a4(2,i,k)+a4(3,i,k))
            endif
          endif
       enddo
     elseif ( abs(kord)==13 ) then
       do i=i1,i2
          if( ext6(i,k) ) then
             if ( ext6(i,k-1) .and. ext6(i,k+1) ) then
! grid-scale 2-delta-z wave detected
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
             endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==14 ) then

       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo

     elseif ( abs(kord)==15 ) then   ! revised kord=9 scehem
       do i=i1,i2
          if ( ext5(i,k) ) then  ! c90_mp122
             if ( ext5(i,k-1) .or. ext5(i,k+1) ) then  ! c90_mp122
! grid-scale 2-delta-z wave detected
                  a4(2,i,k) = a4(1,i,k)
                  a4(3,i,k) = a4(1,i,k)
             endif
          elseif( ext6(i,k) ) then
! Check within the smooth region if subgrid profile is non-monotonic
                  pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                  lac_1 = pmp_1 + 1.5*gam(i,k+2)
              a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),  &
                                             max(a4(1,i,k), pmp_1, lac_1) )
                  pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                  lac_2 = pmp_2 - 1.5*gam(i,k-1)
              a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),  &
                                             max(a4(1,i,k), pmp_2, lac_2) )
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     elseif ( abs(kord)==16 ) then
       do i=i1,i2
          if( ext5(i,k) ) then
             if ( ext5(i,k-1) .or. ext5(i,k+1) ) then
                 a4(2,i,k) = a4(1,i,k)
                 a4(3,i,k) = a4(1,i,k)
             elseif ( ext6(i,k-1) .or. ext6(i,k+1) ) then
                 ! Left  edges
                 pmp_1 = a4(1,i,k) - 2.*gam(i,k+1)
                 lac_1 = pmp_1 + 1.5*gam(i,k+2)
                 a4(2,i,k) = min(max(a4(2,i,k), min(a4(1,i,k), pmp_1, lac_1)),   &
                                     max(a4(1,i,k), pmp_1, lac_1) )
                 ! Right edges
                 pmp_2 = a4(1,i,k) + 2.*gam(i,k)
                 lac_2 = pmp_2 - 1.5*gam(i,k-1)
                 a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), pmp_2, lac_2)),    &
                                     max(a4(1,i,k), pmp_2, lac_2) )
             endif
          endif
       enddo
       do i=i1,i2
          a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
       enddo
     else      ! kord = 11
       do i=i1,i2
         if ( ext5(i,k) .and. (ext5(i,k-1) .or. ext5(i,k+1)) ) then
! Noisy region:
              a4(2,i,k) = a4(1,i,k)
              a4(3,i,k) = a4(1,i,k)
              a4(4,i,k) = 0.
         else
              a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         endif
       enddo
     endif

! Additional constraint to ensure positivity
     if ( iv==0 ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 0)

  enddo      ! k-loop

!----------------------------------
! Bottom layer subgrid constraints:
!----------------------------------
  if ( iv==0 ) then
     do i=i1,i2
        a4(3,i,km) = max(0., a4(3,i,km))
     enddo
  elseif ( iv .eq. -1 ) then
      do i=i1,i2
         if ( a4(3,i,km)*a4(1,i,km) <= 0. )  a4(3,i,km) = 0.
      enddo
  endif

  do k=km-1,km
     do i=i1,i2
        a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
     enddo
     if(k==(km-1)) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 2)
     if(k== km   ) call cs_limiters(im, extm(i1,k), a4(1,i1,k), 1)
  enddo

 end subroutine cs_profile


 subroutine cs_limiters(im, extm, a4, iv)
 integer, intent(in) :: im
 integer, intent(in) :: iv
 logical, intent(in) :: extm(im)
 real , intent(inout) :: a4(4,im)   !< PPM array
! LOCAL VARIABLES:
 real  da1, da2, a6da
 integer i

 if ( iv==0 ) then
! Positive definite constraint
    do i=1,im
    if( a4(1,i)<=0.) then
        a4(2,i) = a4(1,i)
        a4(3,i) = a4(1,i)
        a4(4,i) = 0.
    else
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
         if( (a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12) < 0. ) then
! local minimum is negative
             if( a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i) ) then
                 a4(3,i) = a4(1,i)
                 a4(2,i) = a4(1,i)
                 a4(4,i) = 0.
             elseif( a4(3,i) > a4(2,i) ) then
                 a4(4,i) = 3.*(a4(2,i)-a4(1,i))
                 a4(3,i) = a4(2,i) - a4(4,i)
             else
                 a4(4,i) = 3.*(a4(3,i)-a4(1,i))
                 a4(2,i) = a4(3,i) - a4(4,i)
             endif
         endif
      endif
    endif
    enddo
 elseif ( iv==1 ) then
    do i=1,im
      if( (a4(1,i)-a4(2,i))*(a4(1,i)-a4(3,i))>=0. ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 else
! Standard PPM constraint
    do i=1,im
      if( extm(i) ) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
    enddo
 endif
 end subroutine cs_limiters



 subroutine ppm_profile(a4, delp, km, i1, i2, iv, kord)

! !INPUT PARAMETERS:
 integer, intent(in):: iv      !< iv =-1: winds
                               !! iv = 0: positive definite scalars
                               !! iv = 1: others
                               !! iv = 2: temp (if remap_t) and w (iv=-2)
 integer, intent(in):: i1      !< Starting longitude
 integer, intent(in):: i2      !< Finishing longitude
 integer, intent(in):: km      !< vertical dimension
 integer, intent(in):: kord    !< Order (or more accurately method no.):
                               !!
 real , intent(in):: delp(i1:i2,km)     !< layer pressure thickness

! !INPUT/OUTPUT PARAMETERS:
 real , intent(inout):: a4(4,i1:i2,km)  !< Interpolated values

! DESCRIPTION:
!
!   Perform the piecewise parabolic reconstruction
!
! !REVISION HISTORY:
! S.-J. Lin   revised at GFDL 2007
!-----------------------------------------------------------------------
! local arrays:
      real    dc(i1:i2,km)
      real    h2(i1:i2,km)
      real  delq(i1:i2,km)
      real   df2(i1:i2,km)
      real    d4(i1:i2,km)

! local scalars:
      integer i, k, km1, lmt, it
      real  fac
      real  a1, a2, c1, c2, c3, d1, d2
      real  qm, dq, lac, qmp, pmp

      km1 = km - 1
       it = i2 - i1 + 1

      do k=2,km
         do i=i1,i2
            delq(i,k-1) =   a4(1,i,k) - a4(1,i,k-1)
              d4(i,k  ) = delp(i,k-1) + delp(i,k)
         enddo
      enddo

      do k=2,km1
         do i=i1,i2
                 c1  = (delp(i,k-1)+0.5*delp(i,k))/d4(i,k+1)
                 c2  = (delp(i,k+1)+0.5*delp(i,k))/d4(i,k)
            df2(i,k) = delp(i,k)*(c1*delq(i,k) + c2*delq(i,k-1)) /      &
                                    (d4(i,k)+delp(i,k+1))
            dc(i,k) = sign( min(abs(df2(i,k)),              &
                            max(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))-a4(1,i,k),  &
                  a4(1,i,k)-min(a4(1,i,k-1),a4(1,i,k),a4(1,i,k+1))), df2(i,k) )
         enddo
      enddo

!-----------------------------------------------------------
! 4th order interpolation of the provisional cell edge value
!-----------------------------------------------------------

      do k=3,km1
         do i=i1,i2
            c1 = delq(i,k-1)*delp(i,k-1) / d4(i,k)
            a1 = d4(i,k-1) / (d4(i,k) + delp(i,k-1))
            a2 = d4(i,k+1) / (d4(i,k) + delp(i,k))
            a4(2,i,k) = a4(1,i,k-1) + c1 + 2./(d4(i,k-1)+d4(i,k+1)) *    &
                      ( delp(i,k)*(c1*(a1 - a2)+a2*dc(i,k-1)) -          &
                        delp(i,k-1)*a1*dc(i,k  ) )
         enddo
      enddo

!     if(km>8 .and. kord>4) call steepz(i1, i2, km, a4, df2, dc, delq, delp, d4)

! Area preserving cubic with 2nd deriv. = 0 at the boundaries
! Top
      do i=i1,i2
         d1 = delp(i,1)
         d2 = delp(i,2)
         qm = (d2*a4(1,i,1)+d1*a4(1,i,2)) / (d1+d2)
         dq = 2.*(a4(1,i,2)-a4(1,i,1)) / (d1+d2)
         c1 = 4.*(a4(2,i,3)-qm-d2*dq) / ( d2*(2.*d2*d2+d1*(d2+3.*d1)) )
         c3 = dq - 0.5*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,2) = qm - 0.25*c1*d1*d2*(d2+3.*d1)
! Top edge:
!-------------------------------------------------------
         a4(2,i,1) = d1*(2.*c1*d1**2-c3) + a4(2,i,2)
!-------------------------------------------------------
!        a4(2,i,1) = (12./7.)*a4(1,i,1)-(13./14.)*a4(1,i,2)+(3./14.)*a4(1,i,3)
!-------------------------------------------------------
! No over- and undershoot condition
         a4(2,i,2) = max( a4(2,i,2), min(a4(1,i,1), a4(1,i,2)) )
         a4(2,i,2) = min( a4(2,i,2), max(a4(1,i,1), a4(1,i,2)) )
         dc(i,1) =  0.5*(a4(2,i,2) - a4(1,i,1))
      enddo

! Enforce monotonicity  within the top layer

      if( iv==0 ) then
         do i=i1,i2
            a4(2,i,1) = max(0., a4(2,i,1))
            a4(2,i,2) = max(0., a4(2,i,2))
         enddo
      elseif( iv==-1 ) then
         do i=i1,i2
            if ( a4(2,i,1)*a4(1,i,1) <= 0. ) a4(2,i,1) = 0.
         enddo
      elseif( abs(iv)==2 ) then
         do i=i1,i2
            a4(2,i,1) = a4(1,i,1)
            a4(3,i,1) = a4(1,i,1)
         enddo
      endif

! Bottom
! Area preserving cubic with 2nd deriv. = 0 at the surface
      do i=i1,i2
         d1 = delp(i,km)
         d2 = delp(i,km1)
         qm = (d2*a4(1,i,km)+d1*a4(1,i,km1)) / (d1+d2)
         dq = 2.*(a4(1,i,km1)-a4(1,i,km)) / (d1+d2)
         c1 = (a4(2,i,km1)-qm-d2*dq) / (d2*(2.*d2*d2+d1*(d2+3.*d1)))
         c3 = dq - 2.0*c1*(d2*(5.*d1+d2)-3.*d1*d1)
         a4(2,i,km) = qm - c1*d1*d2*(d2+3.*d1)
! Bottom edge:
!-----------------------------------------------------
         a4(3,i,km) = d1*(8.*c1*d1**2-c3) + a4(2,i,km)
!        dc(i,km) = 0.5*(a4(3,i,km) - a4(1,i,km))
!-----------------------------------------------------
!        a4(3,i,km) = (12./7.)*a4(1,i,km)-(13./14.)*a4(1,i,km-1)+(3./14.)*a4(1,i,km-2)
! No over- and under-shoot condition
         a4(2,i,km) = max( a4(2,i,km), min(a4(1,i,km), a4(1,i,km1)) )
         a4(2,i,km) = min( a4(2,i,km), max(a4(1,i,km), a4(1,i,km1)) )
         dc(i,km) = 0.5*(a4(1,i,km) - a4(2,i,km))
      enddo


! Enforce constraint on the "slope" at the surface

#ifdef BOT_MONO
      do i=i1,i2
         a4(4,i,km) = 0
         if( a4(3,i,km) * a4(1,i,km) <= 0. ) a4(3,i,km) = 0.
         d1 = a4(1,i,km) - a4(2,i,km)
         d2 = a4(3,i,km) - a4(1,i,km)
         if ( d1*d2 < 0. ) then
              a4(2,i,km) = a4(1,i,km)
              a4(3,i,km) = a4(1,i,km)
         else
              dq = sign(min(abs(d1),abs(d2),0.5*abs(delq(i,km-1))), d1)
              a4(2,i,km) = a4(1,i,km) - dq
              a4(3,i,km) = a4(1,i,km) + dq
         endif
      enddo
#else
      if( iv==0 ) then
          do i=i1,i2
             a4(2,i,km) = max(0.,a4(2,i,km))
             a4(3,i,km) = max(0.,a4(3,i,km))
          enddo
      elseif( iv<0 ) then
          do i=i1,i2
             if( a4(1,i,km)*a4(3,i,km) <= 0. )  a4(3,i,km) = 0.
          enddo
      endif
#endif

   do k=1,km1
      do i=i1,i2
         a4(3,i,k) = a4(2,i,k+1)
      enddo
   enddo

!-----------------------------------------------------------
! f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
!-----------------------------------------------------------
! Top 2 and bottom 2 layers always use monotonic mapping
      do k=1,2
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

      if(kord >= 7) then
!-----------------------
! Huynh's 2nd constraint
!-----------------------
      do k=2,km1
         do i=i1,i2
! Method#1
!           h2(i,k) = delq(i,k) - delq(i,k-1)
! Method#2 - better
            h2(i,k) = 2.*(dc(i,k+1)/delp(i,k+1) - dc(i,k-1)/delp(i,k-1))  &
                     / ( delp(i,k)+0.5*(delp(i,k-1)+delp(i,k+1)) )        &
                     * delp(i,k)**2
! Method#3
!!!            h2(i,k) = dc(i,k+1) - dc(i,k-1)
         enddo
      enddo

      fac = 1.5           ! original quasi-monotone

      do k=3,km-2
        do i=i1,i2
! Right edges
!        qmp   = a4(1,i,k) + 2.0*delq(i,k-1)
!        lac   = a4(1,i,k) + fac*h2(i,k-1) + 0.5*delq(i,k-1)
!
         pmp   = 2.*dc(i,k)
         qmp   = a4(1,i,k) + pmp
         lac   = a4(1,i,k) + fac*h2(i,k-1) + dc(i,k)
         a4(3,i,k) = min(max(a4(3,i,k), min(a4(1,i,k), qmp, lac)),    &
                                        max(a4(1,i,k), qmp, lac) )
! Left  edges
!        qmp   = a4(1,i,k) - 2.0*delq(i,k)
!        lac   = a4(1,i,k) + fac*h2(i,k+1) - 0.5*delq(i,k)
!
         qmp   = a4(1,i,k) - pmp
         lac   = a4(1,i,k) + fac*h2(i,k+1) - dc(i,k)
         a4(2,i,k) = min(max(a4(2,i,k),  min(a4(1,i,k), qmp, lac)),   &
                     max(a4(1,i,k), qmp, lac))
!-------------
! Recompute A6
!-------------
         a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
        enddo
! Additional constraint to ensure positivity when kord=7
         if (iv == 0 .and. kord >= 6 )                      &
             call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 2)
      enddo

      else

         lmt = kord - 3
         lmt = max(0, lmt)
         if (iv == 0) lmt = min(2, lmt)

         do k=3,km-2
            if( kord /= 4) then
              do i=i1,i2
                 a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
              enddo
            endif
            if(kord/=6) call ppm_limiters(dc(i1,k), a4(1,i1,k), it, lmt)
         enddo
      endif

      do k=km1,km
         do i=i1,i2
            a4(4,i,k) = 3.*(2.*a4(1,i,k) - (a4(2,i,k)+a4(3,i,k)))
         enddo
         call ppm_limiters(dc(i1,k), a4(1,i1,k), it, 0)
      enddo

 end subroutine ppm_profile


 subroutine ppm_limiters(dm, a4, itot, lmt)

! INPUT PARAMETERS:
      real , intent(in):: dm(*)     !< Linear slope
      integer, intent(in) :: itot      !< Total Longitudes
      integer, intent(in) :: lmt       !< 0: Standard PPM constraint 1: Improved full monotonicity constraint
                                       !< (Lin) 2: Positive definite constraint
                                       !< 3: do nothing (return immediately)
! INPUT/OUTPUT PARAMETERS:
      real , intent(inout) :: a4(4,*)   !< PPM array AA <-- a4(1,i) AL <-- a4(2,i) AR <-- a4(3,i) A6 <-- a4(4,i)
! LOCAL VARIABLES:
      real  qmp
      real  da1, da2, a6da
      real  fmin
      integer i

! Developer: S.-J. Lin

      if ( lmt == 3 ) return

      if(lmt == 0) then
! Standard PPM constraint
      do i=1,itot
      if(dm(i) == 0.) then
         a4(2,i) = a4(1,i)
         a4(3,i) = a4(1,i)
         a4(4,i) = 0.
      else
         da1  = a4(3,i) - a4(2,i)
         da2  = da1**2
         a6da = a4(4,i)*da1
         if(a6da < -da2) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         elseif(a6da > da2) then
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
      endif
      enddo

      elseif (lmt == 1) then

! Improved full monotonicity constraint (Lin 2004)
! Note: no need to provide first guess of A6 <-- a4(4,i)
      do i=1, itot
           qmp = 2.*dm(i)
         a4(2,i) = a4(1,i)-sign(min(abs(qmp),abs(a4(2,i)-a4(1,i))), qmp)
         a4(3,i) = a4(1,i)+sign(min(abs(qmp),abs(a4(3,i)-a4(1,i))), qmp)
         a4(4,i) = 3.*( 2.*a4(1,i) - (a4(2,i)+a4(3,i)) )
      enddo

      elseif (lmt == 2) then

! Positive definite constraint
      do i=1,itot
      if( abs(a4(3,i)-a4(2,i)) < -a4(4,i) ) then
      fmin = a4(1,i)+0.25*(a4(3,i)-a4(2,i))**2/a4(4,i)+a4(4,i)*r12
         if( fmin < 0. ) then
         if(a4(1,i)<a4(3,i) .and. a4(1,i)<a4(2,i)) then
            a4(3,i) = a4(1,i)
            a4(2,i) = a4(1,i)
            a4(4,i) = 0.
         elseif(a4(3,i) > a4(2,i)) then
            a4(4,i) = 3.*(a4(2,i)-a4(1,i))
            a4(3,i) = a4(2,i) - a4(4,i)
         else
            a4(4,i) = 3.*(a4(3,i)-a4(1,i))
            a4(2,i) = a4(3,i) - a4(4,i)
         endif
         endif
      endif
      enddo

      endif

 end subroutine ppm_limiters



 subroutine steepz(i1, i2, km, a4, df2, dm, dq, dp, d4)
 integer, intent(in) :: km, i1, i2
   real , intent(in) ::  dp(i1:i2,km)       !< Grid size
   real , intent(in) ::  dq(i1:i2,km)       !< Backward diff of q
   real , intent(in) ::  d4(i1:i2,km)       !< Backward sum:  dp(k)+ dp(k-1)
   real , intent(in) :: df2(i1:i2,km)       !< First guess mismatch
   real , intent(in) ::  dm(i1:i2,km)       !< Monotonic mismatch
! INPUT/OUTPUT PARAMETERS:
      real , intent(inout) ::  a4(4,i1:i2,km)  !<First guess/steepened
! LOCAL VARIABLES:
      integer i, k
      real  alfa(i1:i2,km)
      real     f(i1:i2,km)
      real   rat(i1:i2,km)
      real   dg2

! Compute ratio of dq/dp
      do k=2,km
         do i=i1,i2
            rat(i,k) = dq(i,k-1) / d4(i,k)
         enddo
      enddo

! Compute F
      do k=2,km-1
         do i=i1,i2
            f(i,k) =   (rat(i,k+1) - rat(i,k))                          &
                     / ( dp(i,k-1)+dp(i,k)+dp(i,k+1) )
         enddo
      enddo

      do k=3,km-2
         do i=i1,i2
         if(f(i,k+1)*f(i,k-1)<0. .and. df2(i,k)/=0.) then
            dg2 = (f(i,k+1)-f(i,k-1))*((dp(i,k+1)-dp(i,k-1))**2          &
                   + d4(i,k)*d4(i,k+1) )
            alfa(i,k) = max(0., min(0.5, -0.1875*dg2/df2(i,k)))
         else
            alfa(i,k) = 0.
         endif
         enddo
      enddo

      do k=4,km-2
         do i=i1,i2
            a4(2,i,k) = (1.-alfa(i,k-1)-alfa(i,k)) * a4(2,i,k) +         &
                        alfa(i,k-1)*(a4(1,i,k)-dm(i,k))    +             &
                        alfa(i,k)*(a4(1,i,k-1)+dm(i,k-1))
         enddo
      enddo

 end subroutine steepz

!>@brief The subroutine 'rst_remap' remaps all variables required for a restart.
!>@details npz_restart /= npz (i.e., when the number of vertical levels is
!! changed at restart)
 subroutine rst_remap(km, kn, is,ie,js,je, isd,ied,jsd,jed, nq, ntp, &
                      delp_r, u_r, v_r, w_r, delz_r, pt_r, q_r, qdiag_r, &
                      delp,   u,   v,   w,   delz,   pt,   q,   qdiag,   &
                      ak_r, bk_r, ptop, ak, bk, hydrostatic, make_nh, &
                      domain, square_domain)
!------------------------------------
! Assuming hybrid sigma-P coordinate:
!------------------------------------
! INPUT PARAMETERS:
  integer, intent(in):: km                    !< Restart z-dimension
  integer, intent(in):: kn                    !< Run time dimension
  integer, intent(in):: nq, ntp               !< Number of tracers (including H2O)
  integer, intent(in):: is,ie,isd,ied         !< Starting & ending X-Dir index
  integer, intent(in):: js,je,jsd,jed         !< Starting & ending Y-Dir index
  logical, intent(in):: hydrostatic, make_nh, square_domain
  real, intent(IN) :: ptop
  real, intent(in) :: ak_r(km+1)
  real, intent(in) :: bk_r(km+1)
  real, intent(in) :: ak(kn+1)
  real, intent(in) :: bk(kn+1)
  real, intent(in):: delp_r(is:ie,js:je,km) !< Pressure thickness
  real, intent(in)::   u_r(is:ie,  js:je+1,km)   !< u-wind (m/s)
  real, intent(in)::   v_r(is:ie+1,js:je  ,km)   !< v-wind (m/s)
  real, intent(inout)::  pt_r(is:ie,js:je,km)
  real, intent(in)::   w_r(is:ie,js:je,km)
  real, intent(in)::   q_r(is:ie,js:je,km,1:ntp)
  real, intent(in)::   qdiag_r(is:ie,js:je,km,ntp+1:nq)
  real, intent(inout)::delz_r(is:ie,js:je,km)
  type(domain2d), intent(INOUT) :: domain
! Output:
  real, intent(out):: delp(isd:ied,jsd:jed,kn) !< Pressure thickness
  real, intent(out)::  u(isd:ied  ,jsd:jed+1,kn)   !< u-wind (m/s)
  real, intent(out)::  v(isd:ied+1,jsd:jed  ,kn)   !< v-wind (m/s)
  real, intent(out)::  w(isd:     ,jsd:     ,1:)   !< Vertical velocity (m/s)
  real, intent(out):: pt(isd:ied  ,jsd:jed  ,kn)   !< Temperature
  real, intent(out):: q(isd:ied,jsd:jed,kn,1:ntp)
  real, intent(out):: qdiag(isd:ied,jsd:jed,kn,ntp+1:nq)
  real, intent(out):: delz(is:,js:,1:)   !< Delta-height (m)
!-----------------------------------------------------------------------
  real r_vir, rgrav
  real ps(isd:ied,jsd:jed)  !< Surface pressure
  real  pe1(is:ie,km+1)
  real  pe2(is:ie,kn+1)
  real  pv1(is:ie+1,km+1)
  real  pv2(is:ie+1,kn+1)

  integer i,j,k , iq
  integer, parameter:: kord=4

#ifdef HYDRO_DELZ_REMAP
  if (is_master() .and. .not. hydrostatic) then
     print*, ''
     print*, ' REMAPPING IC: INITIALIZING DELZ WITH HYDROSTATIC STATE  '
     print*, ''
  endif
#endif

#ifdef HYDRO_DELZ_EXTRAP
  if (is_master() .and. .not. hydrostatic) then
     print*, ''
     print*, ' REMAPPING IC: INITIALIZING DELZ WITH HYDROSTATIC STATE ABOVE INPUT MODEL TOP  '
     print*, ''
  endif
#endif

#ifdef ZERO_W_EXTRAP
  if (is_master() .and. .not. hydrostatic) then
     print*, ''
     print*, ' REMAPPING IC: INITIALIZING W TO ZERO ABOVE INPUT MODEL TOP  '
     print*, ''
  endif
#endif

  r_vir = rvgas/rdgas - 1.
  rgrav = 1./grav

!$OMP parallel do default(none) shared(is,ie,js,je,ps,ak_r)
  do j=js,je
     do i=is,ie
        ps(i,j) = ak_r(1)
     enddo
  enddo

! this OpenMP do-loop setup cannot work in it's current form....
!$OMP parallel do default(none) shared(is,ie,js,je,km,ps,delp_r)
  do j=js,je
     do k=1,km
        do i=is,ie
           ps(i,j) = ps(i,j) + delp_r(i,j,k)
        enddo
     enddo
  enddo

! only one cell is needed
  if ( square_domain ) then
      call mpp_update_domains(ps, domain,  whalo=1, ehalo=1, shalo=1, nhalo=1, complete=.true.)
  else
      call mpp_update_domains(ps, domain, complete=.true.)
  endif

! Compute virtual Temp
!$OMP parallel do default(none) shared(is,ie,js,je,km,pt_r,r_vir,q_r)
  do k=1,km
     do j=js,je
        do i=is,ie
#ifdef MULTI_GASES
           pt_r(i,j,k) = pt_r(i,j,k) * virq(q_r(i,j,k,:))
#else
           pt_r(i,j,k) = pt_r(i,j,k) * (1.+r_vir*q_r(i,j,k,1))
#endif
        enddo
     enddo
  enddo

!$OMP parallel do default(none) shared(is,ie,js,je,km,ak_r,bk_r,ps,kn,ak,bk,u_r,u,delp, &
!$OMP                                  ntp,nq,hydrostatic,make_nh,w_r,w,delz_r,delp_r,delz, &
!$OMP                                  pt_r,pt,v_r,v,q,q_r,qdiag,qdiag_r) &
!$OMP                          private(pe1,  pe2, pv1, pv2)
  do 1000 j=js,je+1
!------
! map u
!------
     do k=1,km+1
        do i=is,ie
           pe1(i,k) = ak_r(k) + 0.5*bk_r(k)*(ps(i,j-1)+ps(i,j))
        enddo
     enddo

     do k=1,kn+1
        do i=is,ie
           pe2(i,k) = ak(k) + 0.5*bk(k)*(ps(i,j-1)+ps(i,j))
        enddo
     enddo

     call remap_2d(km, pe1, u_r(is:ie,j:j,1:km),       &
                   kn, pe2,   u(is:ie,j:j,1:kn),       &
                   is, ie, -1, kord)

  if ( j /= (je+1) )  then

!---------------
! Hybrid sigma-p
!---------------
     do k=1,km+1
        do i=is,ie
           pe1(i,k) = ak_r(k) + bk_r(k)*ps(i,j)
        enddo
     enddo

     do k=1,kn+1
        do i=is,ie
           pe2(i,k) =   ak(k) + bk(k)*ps(i,j)
        enddo
     enddo

!-------------
! Compute delp
!-------------
      do k=1,kn
         do i=is,ie
            delp(i,j,k) = pe2(i,k+1) - pe2(i,k)
         enddo
      enddo

!----------------
! Map constituents
!----------------
      if( nq /= 0 ) then
          do iq=1,ntp
             call remap_2d(km, pe1, q_r(is:ie,j:j,1:km,iq:iq),  &
                           kn, pe2,   q(is:ie,j:j,1:kn,iq:iq),  &
                           is, ie, 0, kord)
          enddo
          do iq=ntp+1,nq
             call remap_2d(km, pe1, qdiag_r(is:ie,j:j,1:km,iq:iq),  &
                           kn, pe2,   qdiag(is:ie,j:j,1:kn,iq:iq),  &
                           is, ie, 0, kord)
          enddo
      endif

      if ( .not. hydrostatic .and. .not. make_nh) then
! Remap vertical wind:
         call remap_2d(km, pe1, w_r(is:ie,j:j,1:km),       &
                       kn, pe2,   w(is:ie,j:j,1:kn),       &
                       is, ie, -1, kord)

#ifdef ZERO_W_EXTRAP
       do k=1,kn
       do i=is,ie
          if (pe2(i,k) < pe1(i,1)) then
             w(i,j,k) = 0.
          endif
       enddo
       enddo
#endif

#ifndef HYDRO_DELZ_REMAP
! Remap delz for hybrid sigma-p coordinate
         do k=1,km
            do i=is,ie
               delz_r(i,j,k) = -delz_r(i,j,k)/delp_r(i,j,k) ! ="specific volume"/grav
            enddo
         enddo
         call remap_2d(km, pe1, delz_r(is:ie,j:j,1:km),       &
                       kn, pe2,   delz(is:ie,j:j,1:kn),       &
                       is, ie, 1, kord)
         do k=1,kn
            do i=is,ie
               delz(i,j,k) = -delz(i,j,k)*delp(i,j,k)
            enddo
         enddo
#endif
      endif

! Geopotential conserving remap of virtual temperature:
       do k=1,km+1
          do i=is,ie
             pe1(i,k) = log(pe1(i,k))
          enddo
       enddo
       do k=1,kn+1
          do i=is,ie
             pe2(i,k) = log(pe2(i,k))
          enddo
       enddo

       call remap_2d(km, pe1, pt_r(is:ie,j:j,1:km),       &
                     kn, pe2,   pt(is:ie,j:j,1:kn),       &
                     is, ie, 1, kord)

#ifdef HYDRO_DELZ_REMAP
       !initialize delz from the hydrostatic state
       do k=1,kn
       do i=is,ie
          delz(i,j,k) = (rdgas*rgrav)*pt(i,j,k)*(pe2(i,k)-pe2(i,k+1))
       enddo
       enddo
#endif
#ifdef HYDRO_DELZ_EXTRAP
       !initialize delz from the hydrostatic state
       do k=1,kn
       do i=is,ie
          if (pe2(i,k) < pe1(i,1)) then
             delz(i,j,k) = (rdgas*rgrav)*pt(i,j,k)*(pe2(i,k)-pe2(i,k+1))
          endif
       enddo
       enddo
#endif
!------
! map v
!------
       do k=1,km+1
          do i=is,ie+1
             pv1(i,k) = ak_r(k) + 0.5*bk_r(k)*(ps(i-1,j)+ps(i,j))
          enddo
       enddo
       do k=1,kn+1
          do i=is,ie+1
             pv2(i,k) = ak(k) + 0.5*bk(k)*(ps(i-1,j)+ps(i,j))
          enddo
       enddo

       call remap_2d(km, pv1, v_r(is:ie+1,j:j,1:km),       &
                     kn, pv2,   v(is:ie+1,j:j,1:kn),       &
                     is, ie+1, -1, kord)

  endif !(j < je+1)
1000  continue

!$OMP parallel do default(none) shared(is,ie,js,je,kn,pt,r_vir,q)
  do k=1,kn
     do j=js,je
        do i=is,ie
#ifdef MULTI_GASES
           pt(i,j,k) = pt(i,j,k) / virq(q(i,j,k,:))
#else
           pt(i,j,k) = pt(i,j,k) / (1.+r_vir*q(i,j,k,1))
#endif
        enddo
     enddo
  enddo

 end subroutine rst_remap

!>@brief The subroutine 'mappm' is a general-purpose routine for remapping
!! one set of vertical levels to another.
 subroutine mappm(km, pe1, q1, kn, pe2, q2, i1, i2, iv, kord, ptop)

! IV = 0: constituents
! IV = 1: potential temp
! IV =-1: winds

! Mass flux preserving mapping: q1(im,km) -> q2(im,kn)

! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate

 integer, intent(in):: i1, i2, km, kn, kord, iv
 real, intent(in ):: pe1(i1:i2,km+1), pe2(i1:i2,kn+1) !< pe1: pressure at layer edges from model top to bottom
                                                      !!      surface in the ORIGINAL vertical coordinate
                                                      !< pe2: pressure at layer edges from model top to bottom
                                                      !!      surface in the NEW vertical coordinate
! Mass flux preserving mapping: q1(im,km) -> q2(im,kn)
 real, intent(in )::  q1(i1:i2,km)
 real, intent(out)::  q2(i1:i2,kn)
 real, intent(IN) :: ptop
! local
      real  qs(i1:i2)
      real dp1(i1:i2,km)
      real a4(4,i1:i2,km)
      integer i, k, l
      integer k0, k1
      real pl, pr, tt, delp, qsum, dpsum, esl

      do k=1,km
         do i=i1,i2
             dp1(i,k) = pe1(i,k+1) - pe1(i,k)
            a4(1,i,k) = q1(i,k)
         enddo
      enddo

      if ( kord >7 ) then
           call  cs_profile( qs, a4, dp1, km, i1, i2, iv, kord )
      else
           call ppm_profile( a4, dp1, km, i1, i2, iv, kord )
      endif

!------------------------------------
! Lowest layer: constant distribution
!------------------------------------
#ifdef NGGPS_SUBMITTED
      do i=i1,i2
         a4(2,i,km) = q1(i,km)
         a4(3,i,km) = q1(i,km)
         a4(4,i,km) = 0.
      enddo
#endif

      do 5555 i=i1,i2
         k0 = 1
      do 555 k=1,kn

         if(pe2(i,k) .le. pe1(i,1)) then
! above old ptop
            q2(i,k) = q1(i,1)
         elseif(pe2(i,k) .ge. pe1(i,km+1)) then
! Entire grid below old ps
#ifdef NGGPS_SUBMITTED
            q2(i,k) = a4(3,i,km)   ! this is not good.
#else
            q2(i,k) = q1(i,km)
#endif
         else

         do 45 L=k0,km
! locate the top edge at pe2(i,k)
         if( pe2(i,k) .ge. pe1(i,L) .and.        &
             pe2(i,k) .le. pe1(i,L+1)    ) then
             k0 = L
             PL = (pe2(i,k)-pe1(i,L)) / dp1(i,L)
             if(pe2(i,k+1) .le. pe1(i,L+1)) then

! entire new grid is within the original grid
               PR = (pe2(i,k+1)-pe1(i,L)) / dp1(i,L)
               TT = r3*(PR*(PR+PL)+PL**2)
               q2(i,k) = a4(2,i,L) + 0.5*(a4(4,i,L)+a4(3,i,L)  &
                       - a4(2,i,L))*(PR+PL) - a4(4,i,L)*TT
              goto 555
             else
! Fractional area...
              delp = pe1(i,L+1) - pe2(i,k)
              TT   = r3*(1.+PL*(1.+PL))
              qsum = delp*(a4(2,i,L)+0.5*(a4(4,i,L)+            &
                     a4(3,i,L)-a4(2,i,L))*(1.+PL)-a4(4,i,L)*TT)
              dpsum = delp
              k1 = L + 1
             goto 111
             endif
         endif
45       continue

111      continue
         do 55 L=k1,km
         if( pe2(i,k+1) .gt. pe1(i,L+1) ) then

! Whole layer..

            qsum  =  qsum + dp1(i,L)*q1(i,L)
            dpsum = dpsum + dp1(i,L)
         else
           delp = pe2(i,k+1)-pe1(i,L)
           esl  = delp / dp1(i,L)
           qsum = qsum + delp * (a4(2,i,L)+0.5*esl*            &
                 (a4(3,i,L)-a4(2,i,L)+a4(4,i,L)*(1.-r23*esl)) )
          dpsum = dpsum + delp
           k0 = L
           goto 123
         endif
55       continue
        delp = pe2(i,k+1) - pe1(i,km+1)
        if(delp > 0.) then
! Extended below old ps
#ifdef NGGPS_SUBMITTED
           qsum = qsum + delp * a4(3,i,km)    ! not good.
#else
           qsum = qsum + delp * q1(i,km)
#endif
          dpsum = dpsum + delp
        endif
123     q2(i,k) = qsum / dpsum
      endif
555   continue
5555  continue

 end subroutine mappm


!>@brief The subroutine 'moist_cv' computes the FV3-consistent moist heat capacity under constant volume,
!! including the heating capacity of water vapor and condensates.
!>@details See \cite emanuel1994atmospheric for information on variable heat capacities.
 subroutine moist_cv(is,ie, isd,ied, jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                     ice_wat, snowwat, graupel, hailwat, q, qd, cvm, t1)
  integer, intent(in):: is, ie, isd,ied, jsd,jed, km, nwat, j, k
  integer, intent(in):: sphum, liq_wat, rainwat, ice_wat, snowwat, graupel, hailwat
#ifdef MULTI_GASES
  real, intent(in), dimension(isd:ied,jsd:jed,km,num_gas):: q
#else
  real, intent(in), dimension(isd:ied,jsd:jed,km,nwat):: q
#endif
  real, intent(out), dimension(is:ie):: cvm, qd  !< qd is q_con
  real, intent(in), optional:: t1(is:ie)
!
  real, parameter:: t_i0 = 15.
  real, dimension(is:ie):: qv, ql, qs
  integer:: i

  select case (nwat)

   case(2)
     if ( present(t1) ) then  ! Special case for GFS physics
        do i=is,ie
           qd(i) = max(0., q(i,j,k,liq_wat))
           if ( t1(i) > tice ) then
                qs(i) = 0.
           elseif ( t1(i) < tice-t_i0 ) then
                qs(i) = qd(i)
           else
                qs(i) = qd(i)*(tice-t1(i))/t_i0
           endif
           ql(i) = qd(i) - qs(i)
           qv(i) = max(0.,q(i,j,k,sphum))
#ifdef MULTI_GASES
           cvm(i) = (1.-(qv(i)+qd(i)))*cv_air*vicvqd(q(i,j,k,1:num_gas)) + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
#else
           cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
#endif
        enddo
     else
        do i=is,ie
           qv(i) = max(0.,q(i,j,k,sphum))
           qs(i) = max(0.,q(i,j,k,liq_wat))
           qd(i) = qs(i)
#ifdef MULTI_GASES
           cvm(i) = (1.-qv(i))*cv_air*vicvqd(q(i,j,k,1:num_gas)) + qv(i)*cv_vap
#else
           cvm(i) = (1.-qv(i))*cv_air + qv(i)*cv_vap
#endif
        enddo
     endif
  case (3)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat)
        qs(i) = q(i,j,k,ice_wat)
        qd(i) = ql(i) + qs(i)
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
     enddo
  case(4)              ! K_warm_rain with fake ice
     do i=is,ie
       qv(i) = q(i,j,k,sphum)
       ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
       qs(i) = q(i,j,k,ice_wat)
       qd(i) = ql(i) + qs(i)
#ifdef MULTI_GASES
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air*vicvqd(q(i,j,k,1:num_gas)) + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
#else
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
#endif
     enddo
  case(5)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        qs(i) = q(i,j,k,ice_wat) + q(i,j,k,snowwat)
        qd(i) = ql(i) + qs(i)
#ifdef MULTI_GASES
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air*vicvqd(q(i,j,k,1:num_gas)) + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
#else
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
#endif
     enddo

  case(6)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        qs(i) = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel)
        qd(i) = ql(i) + qs(i)
#ifdef MULTI_GASES
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air*vicvqd(q(i,j,k,1:num_gas)) + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
#else
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
#endif
     enddo
  case(7)
     do i=is,ie 
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat) 
        qs(i) = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel) + q(i,j,k,hailwat)
        qd(i) = ql(i) + qs(i)
#ifdef MULTI_GASES
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air*vicvqd(q(i,j,k,1:num_gas)) + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
#else
        cvm(i) = (1.-(qv(i)+qd(i)))*cv_air + qv(i)*cv_vap + ql(i)*c_liq + qs(i)*c_ice
#endif
     enddo
  case default
     !call mpp_error (NOTE, 'fv_mapz::moist_cv - using default cv_air')
     do i=is,ie
         qd(i) = 0.
#ifdef MULTI_GASES
        cvm(i) = cv_air*vicvqd(q(i,j,k,1:num_gas))
#else
        cvm(i) = cv_air
#endif
     enddo
 end select

 end subroutine moist_cv

!>@brief The subroutine 'moist_cp' computes the FV3-consistent moist heat capacity under constant pressure,
!! including the heating capacity of water vapor and condensates.
 subroutine moist_cp(is,ie, isd,ied, jsd,jed, km, j, k, nwat, sphum, liq_wat, rainwat,    &
                     ice_wat, snowwat, graupel, hailwat, q, qd, cpm, t1)

  integer, intent(in):: is, ie, isd,ied, jsd,jed, km, nwat, j, k
  integer, intent(in):: sphum, liq_wat, rainwat, ice_wat, snowwat, graupel, hailwat
#ifdef MULTI_GASES
  real, intent(in), dimension(isd:ied,jsd:jed,km,num_gas):: q
#else
  real, intent(in), dimension(isd:ied,jsd:jed,km,nwat):: q
#endif
  real, intent(out), dimension(is:ie):: cpm, qd
  real, intent(in), optional:: t1(is:ie)
!
  real, parameter:: t_i0 = 15.
  real, dimension(is:ie):: qv, ql, qs
  integer:: i

  select case (nwat)

  case(2)
     if ( present(t1) ) then  ! Special case for GFS physics
        do i=is,ie
           qd(i) = max(0., q(i,j,k,liq_wat))
           if ( t1(i) > tice ) then
                qs(i) = 0.
           elseif ( t1(i) < tice-t_i0 ) then
                qs(i) = qd(i)
           else
                qs(i) = qd(i)*(tice-t1(i))/t_i0
           endif
           ql(i) = qd(i) - qs(i)
           qv(i) = max(0.,q(i,j,k,sphum))
#ifdef MULTI_GASES
           cpm(i) = (1.-(qv(i)+qd(i)))*cp_air * vicpqd(q(i,j,k,:)) + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
#else
           cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
#endif
        enddo
     else
     do i=is,ie
        qv(i) = max(0.,q(i,j,k,sphum))
        qs(i) = max(0.,q(i,j,k,liq_wat))
        qd(i) = qs(i)
#ifdef MULTI_GASES
        cpm(i) = (1.-qv(i))*cp_air*vicpqd(q(i,j,k,:)) + qv(i)*cp_vapor
#else
        cpm(i) = (1.-qv(i))*cp_air + qv(i)*cp_vapor
#endif
     enddo
     endif

  case(3)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat)
        qs(i) = q(i,j,k,ice_wat)
        qd(i) = ql(i) + qs(i)
#ifdef MULTI_GASES
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air*vicpqd(q(i,j,k,:)) + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
#else
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
#endif
     enddo
  case(4)    ! K_warm_rain scheme with fake ice
     do i=is,ie
       qv(i) = q(i,j,k,sphum)
       ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
       qs(i) = q(i,j,k,ice_wat)
       qd(i) = ql(i) + qs(i)
#ifdef MULTI_GASES
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air*vicpqd(q(i,j,k,:)) + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
#else
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
#endif
     enddo
  case(5)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        qs(i) = q(i,j,k,ice_wat) + q(i,j,k,snowwat)
        qd(i) = ql(i) + qs(i)
#ifdef MULTI_GASES
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air*vicpqd(q(i,j,k,:)) + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
#else
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
#endif
     enddo

  case(6)
     do i=is,ie
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat)
        qs(i) = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel)
        qd(i) = ql(i) + qs(i)
#ifdef MULTI_GASES
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air*vicpqd(q(i,j,k,:)) + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
#else
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
#endif
     enddo

  case(7)
     do i=is,ie 
        qv(i) = q(i,j,k,sphum)
        ql(i) = q(i,j,k,liq_wat) + q(i,j,k,rainwat) 
        qs(i) = q(i,j,k,ice_wat) + q(i,j,k,snowwat) + q(i,j,k,graupel) + q(i,j,k,hailwat)
        qd(i) = ql(i) + qs(i)
#ifdef MULTI_GASES
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air*vicpqd(q(i,j,k,:)) + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
#else
        cpm(i) = (1.-(qv(i)+qd(i)))*cp_air + qv(i)*cp_vapor + ql(i)*c_liq + qs(i)*c_ice
#endif
     enddo

  case default
     !call mpp_error (NOTE, 'fv_mapz::moist_cp - using default cp_air')
     do i=is,ie
        qd(i) = 0.
#ifdef MULTI_GASES
        cpm(i) = cp_air*vicpqd(q(i,j,k,:))
#else
        cpm(i) = cp_air
#endif
     enddo
  end select

 end subroutine moist_cp
!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  map1_cubic --- Cubic Interpolation for vertical re-mapping
!
! !INTERFACE:
  subroutine map1_cubic( km,   pe1,    q1,                 &
                         kn,   pe2,    q2,   i1, i2,       &
                         j,    ibeg, iend, jbeg, jend, akap, T_VAR, conserv)
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in) :: i1                ! Starting longitude
      integer, intent(in) :: i2                ! Finishing longitude
      real, intent(in) :: akap
      integer, intent(in) :: T_VAR             ! Thermodynamic variable to remap
                                               !     1:TE  2:T  3:PT
      logical, intent(in) :: conserv
      integer, intent(in) :: j                 ! Current latitude
      integer, intent(in) :: ibeg, iend, jbeg, jend
      integer, intent(in) :: km                ! Original vertical dimension
      integer, intent(in) :: kn                ! Target vertical dimension

      real, intent(in) ::  pe1(i1:i2,km+1)  ! pressure at layer edges
                                               ! (from model top to bottom surface)
                                               ! in the original vertical coordinate
      real, intent(in) ::  pe2(i1:i2,kn+1)  ! pressure at layer edges
                                               ! (from model top to bottom surface)
                                               ! in the new vertical coordinate

      real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) ! Field input
! !INPUT/OUTPUT PARAMETERS:
      real, intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) ! Field output

! !DESCRIPTION:
!
!     Perform Cubic Interpolation a given latitude
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate
!
! !REVISION HISTORY:
!    2005.11.14   Takacs    Initial Code
!    2016.07.20   Putman    Modified to make genaric for any thermodynamic variable
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      real       qx(i1:i2,km)
      real   logpl1(i1:i2,km)
      real   logpl2(i1:i2,kn)
      real   dlogp1(i1:i2,km)
      real    vsum1(i1:i2)
      real    vsum2(i1:i2)
      real   am2,am1,ap0,ap1,P,PLP1,PLP0,PLM1,PLM2,DLP0,DLM1,DLM2

      integer i, k, LM2,LM1,LP0,LP1

! Initialization
! --------------

      select case (T_VAR)
      case(1)
       ! Total Energy Remapping in Log(P)
        do k=1,km
            qx(:,k) = q1(i1:i2,j,k)
        logpl1(:,k) = log( 0.5*(pe1(:,k)+pe1(:,k+1)) )
        enddo
        do k=1,kn
        logpl2(:,k) = log( 0.5*(pe2(:,k)+pe2(:,k+1)) )
        enddo

        do k=1,km-1
        dlogp1(:,k) = logpl1(:,k+1)-logpl1(:,k)
        enddo

      case(2)
       ! Temperature Remapping in Log(P)
        do k=1,km
            qx(:,k) = q1(i1:i2,j,k)
        logpl1(:,k) = log( 0.5*(pe1(:,k)+pe1(:,k+1)) )
        enddo
        do k=1,kn
        logpl2(:,k) = log( 0.5*(pe2(:,k)+pe2(:,k+1)) )
        enddo

        do k=1,km-1
        dlogp1(:,k) = logpl1(:,k+1)-logpl1(:,k)
        enddo

      case(3)
       ! Potential Temperature Remapping in P^KAPPA
        do k=1,km
            qx(:,k) = q1(i1:i2,j,k)
        logpl1(:,k) = exp( akap*log( 0.5*(pe1(:,k)+pe1(:,k+1))) )
        enddo
        do k=1,kn
        logpl2(:,k) = exp( akap*log( 0.5*(pe2(:,k)+pe2(:,k+1))) )
        enddo

        do k=1,km-1
        dlogp1(:,k) = logpl1(:,k+1)-logpl1(:,k)
        enddo

      end select

      if (conserv) then
! Compute vertical integral of Input TE
! -------------------------------------
        vsum1(:) = 0.0
        do i=i1,i2
        do k=1,km
        vsum1(i) = vsum1(i) + qx(i,k)*( pe1(i,k+1)-pe1(i,k) )
        enddo
        vsum1(i) = vsum1(i) / ( pe1(i,km+1)-pe1(i,1) )
        enddo

      endif

! Interpolate TE onto target Pressures
! ------------------------------------
      do i=i1,i2
      do k=1,kn
         LM1 = 1
         LP0 = 1
         do while( LP0.le.km )
            if (logpl1(i,LP0).lt.logpl2(i,k)) then
               LP0 = LP0+1
            else
               exit
            endif
         enddo
         LM1 = max(LP0-1,1)
         LP0 = min(LP0, km)

! Extrapolate Linearly in LogP above first model level
! ----------------------------------------------------
         if( LM1.eq.1 .and. LP0.eq.1 ) then
             q2(i,j,k) = qx(i,1) + ( qx(i,2)-qx(i,1) )*( logpl2(i,k)-logpl1(i,1) ) &
                                                      /( logpl1(i,2)-logpl1(i,1) )

! Extrapolate Linearly in LogP below last model level
! ---------------------------------------------------
         else if( LM1.eq.km .and. LP0.eq.km ) then
             q2(i,j,k) = qx(i,km) + ( qx(i,km)-qx(i,km-1) )*( logpl2(i,k )-logpl1(i,km  ) ) &
                                                           /( logpl1(i,km)-logpl1(i,km-1) )

! Interpolate Linearly in LogP between levels 1 => 2 and km-1 => km
! -----------------------------------------------------------------
         else if( LM1.eq.1 .or. LP0.eq.km ) then
             q2(i,j,k) = qx(i,LP0) + ( qx(i,LM1)-qx(i,LP0) )*( logpl2(i,k  )-logpl1(i,LP0) ) &
                                                            /( logpl1(i,LM1)-logpl1(i,LP0) )
! Interpolate Cubicly in LogP between other model levels
! ------------------------------------------------------
         else
              LP1 = LP0+1
              LM2 = LM1-1
             P    = logpl2(i,k)
             PLP1 = logpl1(i,LP1)
             PLP0 = logpl1(i,LP0)
             PLM1 = logpl1(i,LM1)
             PLM2 = logpl1(i,LM2)
             DLP0 = dlogp1(i,LP0)
             DLM1 = dlogp1(i,LM1)
             DLM2 = dlogp1(i,LM2)

              ap1 = (P-PLP0)*(P-PLM1)*(P-PLM2)/( DLP0*(DLP0+DLM1)*(DLP0+DLM1+DLM2) )
              ap0 = (PLP1-P)*(P-PLM1)*(P-PLM2)/( DLP0*      DLM1 *(     DLM1+DLM2) )
              am1 = (PLP1-P)*(PLP0-P)*(P-PLM2)/( DLM1*      DLM2 *(DLP0+DLM1     ) )
              am2 = (PLP1-P)*(PLP0-P)*(PLM1-P)/( DLM2*(DLM1+DLM2)*(DLP0+DLM1+DLM2) )

             q2(i,j,k) = ap1*qx(i,LP1) + ap0*qx(i,LP0) + am1*qx(i,LM1) + am2*qx(i,LM2)

         endif

      enddo
      enddo
      if (conserv) then

! Compute vertical integral of Output TE
! --------------------------------------
        vsum2(:) = 0.0
        do i=i1,i2
        do k=1,kn
        vsum2(i) = vsum2(i) + q2(i,j,k)*( pe2(i,k+1)-pe2(i,k) )
        enddo
        vsum2(i) = vsum2(i) / ( pe2(i,kn+1)-pe2(i,1) )
        enddo

! Adjust Final TE to conserve
! ---------------------------
        do i=i1,i2
        do k=1,kn
           q2(i,j,k) = q2(i,j,k) + vsum1(i)-vsum2(i)
!          q2(i,j,k) = q2(i,j,k) * vsum1(i)/vsum2(i)
        enddo
        enddo

      endif

      return
!EOC
 end subroutine map1_cubic
!============================
!
subroutine map1_ppm_dpwind( km,   pe1,    q1,   qs,           &
                      kn,   pe2,    q2,   dpu,  i1, i2,       &
                      j,ibeg, iend, jbeg, jend,    iv,  kord)
 integer, intent(in) :: i1                !< Starting longitude
 integer, intent(in) :: i2                !< Finishing longitude
 integer, intent(in) :: iv                !< Mode: 0 == constituents 1 == ??? 2 == remap temp with cs scheme
 integer, intent(in) :: kord              !< Method order
 integer, intent(in) :: j                 !< Current latitude
 integer, intent(in) :: ibeg, iend, jbeg, jend
 integer, intent(in) :: km                !< Original vertical dimension
 integer, intent(in) :: kn                !< Target vertical dimension
 real, intent(in) ::   qs(i1:i2)       !< bottom BC
 real, intent(in) ::  pe1(i1:i2,km+1)  !< pressure at layer edges
                                       !! (from model top to bottom surface)
                                       !! in the original vertical coordinate
 real, intent(in) ::  pe2(i1:i2,kn+1)  !< pressure at layer edges
                                       !! (from model top to bottom surface)
                                       !! in the new vertical coordinate
 real, intent(in) ::    q1(ibeg:iend,jbeg:jend,km) !< Field input
! !INPUT/OUTPUT PARAMETERS:
 real, intent(inout)::  q2(ibeg:iend,jbeg:jend,kn) !< Field output
 real, intent(out) ::   dpu(i1:i2, kn) 

! DESCRIPTION:
! IV = 0: constituents
! pe1: pressure at layer edges (from model top to bottom surface)
!      in the original vertical coordinate
! pe2: pressure at layer edges (from model top to bottom surface)
!      in the new vertical coordinate

! LOCAL VARIABLES:
   real    dp1(i1:i2,km)
   
   real   q4(4,i1:i2,km)
   real    pl, pr, qsum, dp, esl
   integer i, k, l, m, k0

   do k=1,km
      do i=i1,i2
         dp1(i,k) = pe1(i,k+1) - pe1(i,k)
         q4(1,i,k) = q1(i,j,k)
      enddo
   enddo

! Compute vertical subgrid distribution
   if ( kord >7 ) then
        call  cs_profile( qs, q4, dp1, km, i1, i2, iv, kord )
   else
        call ppm_profile( q4, dp1, km, i1, i2, iv, kord )
   endif

  do i=i1,i2
     k0 = 1
     do 555 k=1,kn
        dpu(i,k) = pe2(i,k+1) - pe2(i,k)  
      do l=k0,km
! locate the top edge: pe2(i,k)
      if( pe2(i,k) >= pe1(i,l) .and. pe2(i,k) <= pe1(i,l+1) ) then
         pl = (pe2(i,k)-pe1(i,l)) / dp1(i,l)
         if( pe2(i,k+1) <= pe1(i,l+1) ) then
! entire new grid is within the original grid
            pr = (pe2(i,k+1)-pe1(i,l)) / dp1(i,l)
            q2(i,j,k) = q4(2,i,l) + 0.5*(q4(4,i,l)+q4(3,i,l)-q4(2,i,l))  &
                       *(pr+pl)-q4(4,i,l)*r3*(pr*(pr+pl)+pl**2)
               k0 = l
               goto 555
         else
! Fractional area...
            qsum = (pe1(i,l+1)-pe2(i,k))*(q4(2,i,l)+0.5*(q4(4,i,l)+   &
                    q4(3,i,l)-q4(2,i,l))*(1.+pl)-q4(4,i,l)*           &
                     (r3*(1.+pl*(1.+pl))))
              do m=l+1,km
! locate the bottom edge: pe2(i,k+1)
                 if( pe2(i,k+1) > pe1(i,m+1) ) then
! Whole layer
                     qsum = qsum + dp1(i,m)*q4(1,i,m)
                 else
                     dp = pe2(i,k+1)-pe1(i,m)
                     esl = dp / dp1(i,m)
                     qsum = qsum + dp*(q4(2,i,m)+0.5*esl*               &
                           (q4(3,i,m)-q4(2,i,m)+q4(4,i,m)*(1.-r23*esl)))
                     k0 = m
                     goto 123
                 endif
              enddo
              goto 123
         endif
      endif
      enddo
123   continue    
      q2(i,j,k) = qsum / dpu(i,k)            !( pe2(i,k+1) - pe2(i,k) )
555   continue
  enddo

 end subroutine map1_ppm_dpwind
!=============================== 
!
 subroutine get_coef_mdif(is, ie, km, t2, qo, qo2, qo3, qh2o, dp2, pe, grav,  &
                          vumol, ktmol, dfmol, rhomol,cp_mu, am_mol, zgeo)	
 use constants_mod,       only:  rdgas
 implicit none
 
 integer, intent(in)   ::  is, ie , km
 real, dimension(is:ie,km) :: qo, qo2, qh2o, qo3  
 real, dimension(is:ie,km) :: dp2, t2
 real, dimension(is:ie,km+1) :: pe, grav,  vumol, ktmol, dfmol, rhomol
 real, dimension(is:ie,km+1) :: cp_mu, zgeo
 real, dimension(is:ie,km  ) :: am_mol
 real, parameter::  amo=15.9994, amo2=2.*amo, amo3= 3.*amo      
 real, parameter::  amn2=28.013,  amh2o=18.0154    !g/mol
 
!< muo3 and muh2o are not precise, correct later
 real, parameter:: muo=3.9e-7, muo2=4.03e-7,  muo3=4.03e-7     !kg/m/s
 real, parameter::             mun2=3.43e-7,  muh2o=3.43e-7    !kg/m/s
!< lao3 is not precise values, but o3_n is very small
  real, parameter:: lao=75.9e-5, lao2=56.e-5,  lao3=36.e-5     !kg/m/s
  real, parameter::              lan2=56.e-5,  lah2o=55.e-5    !kg/m/s
  real, parameter:: cpo=1299.185, cpo2=918.0969, cpo3=820.2391
  real, parameter:: cpn2=1031.108, cph2o=1846.00
  real, parameter:: avgd=6.0221415e23  ! Avogadro constant
  real, parameter:: bz=1.3806505e-23   ! Boltzmann constant J/K
  real, parameter:: a12=9.69e18 ! O-O2 diffusion params
  real, parameter:: s12=0.774, s121 = s12+1.
  real, parameter:: ktdep = 0.69+1.  
  real :: a12bz,   avgdbz, t69, pmid, mumol, ravgdbz
  real :: vu, kt, df, rhom, tdep
  
  real :: n_o, n_o2, n_o3, n_h2, n_n2, qn2
  real :: cpmult, cvmult,rdmult
  real :: ramo, ramo2, ramo3, ramh2, ramn2, ratio, dpc, dpm 
  real :: vuz(km+1), ktz(km+1), dfz(km+1), rhoz(km+1), wcof, cpx(km+1), dpe(km+1)
  integer  :: i, k
  real :: Runiv, RE, Re2, g981, re2g
      runiv = 8314.
      re = 6370.e3; re2 = re*re
      g981 =9.8065; re2g = re2 * g981
      a12bz = a12 * bz
      avgdbz= avgd * bz *1000.
      ravgdbz= 1./avgdbz
      ramo =1./amo
      ramo2 =.5*ramo 
      ramo3 =ramo/3. 
      ramh2 = 1./amh2o      
      ramn2 =1./amn2 
      
      do i = is, ie
          dpe(1) = pe(i,1)
	  zgeo(i,km+1) = 0.
       do k=km, 1, -1
          dpe(k+1) = dp2(i,k)        
          qn2 = 1. -qo(i,k)-qo2(i,k)-qo3(i,k)-qh2o(i,k)
          mumol = 1./(qo(i,k)*ramo +qo2(i,k)*ramo2+qo3(i,k)*ramo3+qh2o(i,k)*ramh2+qn2*ramn2)
	  am_mol(i,k) = mumol
	  cpx(k+1) = qo(i,k)*cpo +qo2(i,k)*cpo2+qo3(i,k)*cpo3+qh2o(i,k)*cph2o+qn2*cpn2
	  n_o = qo(i,k)* mumol*ramo
	  n_o2= qo2(i,k)*mumol*ramo2
	  n_o3= qo3(i,k)*mumol*ramo3
	  n_h2= qh2o(i,k)*mumol*ramh2
	  n_n2= qn2*mumol*ramn2
          pmid = .5*(pe(i,k)+pe(i,k+1))
          tdep = t2(i,k)  ** ktdep
	  df = a12bz*t2(i,k)**s121/pmid
	  
          vu =  n_o*muo + n_o2*muo2 + n_o3*muo3 + n_n2*mun2 + n_h2*muh2o
          kt =  n_o*lao + n_o2*lao2 + n_o3*lao3 + n_n2*lan2 + n_h2*lah2o	  
          rhom = tdep*runiv/(mumol *pmid) 
	  rhoz(k+1) = 	pmid *mumol/t2(i,k)/runiv      !1e-3 * am * plyr(n)/temp(n) / avgdbz 
	  zgeo(i,k) = zgeo(i,k+1) + dp2(i,k)/rhoz(k+1)/g981	
	  grav(i,k) = re2g/(zgeo(i,k)+re)/(zgeo(i,k)+re)
          vuz(k+1) = vu * rhom 
          ktz(k+1) = kt * rhom 
          dfz(k+1) = df
       enddo
       
       ratio = 1.              !vuz(2)/vuz(3)
       vuz(1) =vuz(2) *ratio 
       ktz(1) =ktz(2) *ratio 
       dfz(1) =dfz(2) *ratio  
          k=1          
          qn2 = 1. -qo(i,k)-qo2(i,k)-qo3(i,k)-qh2o(i,k)
          mumol = 1./(qo(i,k)*ramo +qo2(i,k)*ramo2+qo3(i,k)*ramo3+qh2o(i,k)*ramh2+qn2*ramn2)
          rhoz(1) = pe(i,1) *mumol/t2(i,1)/runiv 
	  cpx(1) = qo(i,k)*cpo +qo2(i,k)*cpo2+qo3(i,k)*cpo3+qh2o(i,k)*cph2o+qn2*cpn2
	  vumol(i,1) = vuz(1)
	  ktmol(i,1) = ktz(1)/cpx(1)
	  dfmol(i,1) = dfz(1)
	  cp_mu(i,k) =cpx(1)
	  do k =1, km         
	    cpmult  =  .5*(cpx(k)+cpx(k+1))
	    vumol(i,k+1) = .5*(vuz(k)+vuz(k+1))
	    ktmol(i,k+1) = .5*(ktz(k)+ktz(k+1))/cpmult
	    dfmol(i,k+1) = .5*(dfz(k)+dfz(k+1))
	    wcof = .5*(rhoz(k)+rhoz(k+1))
	    dpc= .5*(dpe(k)+dpe(k+1))
	    rhomol(i,k+1) =wcof*wcof*grav(i,k+1)/dpc
	    cp_mu(i,k+1) =cpmult
	  enddo   	  
	  k= 1
	  rhomol(i,k) = 0.25*rhoz(k)*rhoz(k)*grav(i,k)/dp2(i,k)   
      enddo

! print *, 'get_coef_mdif'
 end  subroutine get_coef_mdif 
 
 subroutine get_coef_turb(is, ie, km, t2, u2, v2, dp2, pm2, pe2, gvz, cp_mu, am, zm, vuedd, ktedd, dfedd, &
             con_adj, nfixed, nrest)
 integer , intent(in) :: is, ie, km
 integer ::  nfixed, nrest
 
 logical, intent(in)  ::  con_adj
 real, intent(in), dimension(is:ie,km) :: t2, u2, v2 
 real, intent(in)     ::  dp2(is:ie, km), pm2(is:ie, km), pe2(is:ie, km+1)
 real, intent(in)     ::  am(is:ie, km)  
 real, intent(in), dimension(is:ie,km+1)      ::  gvz,  cp_mu, zm
 real, intent(out), dimension(is:ie,km+1)     :: vuedd, ktedd, dfedd
 integer :: i, k, ii, kk ,j
 
 real :: t(km), u(km), v(km), pe(km), dp(km), grav(km),  kap(km)
 
 real, parameter :: lturb = 35., vumin = 0.03, bn2min = 4.e-8, ric = 0.25
 real, parameter :: dw2min = 0.01, dked_min = 1., dked_max= 1.e6
 real, parameter :: lz_smag = 35., cs_smag = 0.53
 real, parameter :: ruvac = 0.25/8314.
 
 real  :: uz, vz, tz, tc, bn2, dzg, gcp, rik, shr2, ritur
 real  :: bn_uns, zmetk, zlturb, w1, kamp, zgrow, dtemp, ht 
 real  :: t1(km), pt(km), tadj(km), kzz_edd(km+1)
 real  :: lturbmax
 integer  :: ins_con(km), num_ins
 vuedd(:,:) = vumin
 ktedd(:,:) = vumin
 dfedd(:,:) = vumin
 
 nfixed =0
 nrest=0
 lturbmax =300.
 do i=is, ie 
     ins_con(1:km) = 0
     num_ins =0 
     kzz_edd(1:km+1) = vumin
  do k=1, km-1
     dzg =  zm(i,k) - zm(i,k+1) 
     uz =  (u2(i,k)     - u2(i,k+1)) 
     vz =  (v2(i,k)     - v2(i,k+1)) 
     tz =  (t2(i,k)     - t2(i,k+1)) /dzg 
     tc = .5*(t2(i,k)   + t2(i,k+1))
     gcp = gvz(i,k)/cp_mu(i,k)
     bn2 = gvz(i,k)*(tz +gcp)/tc
     shr2 = max(uz*uz+vz*vz, dw2min)
     shr2 = shr2/dzg/dzg
     rik = bn2/shr2
     ht = ruvac/tc*gvz(i,k)*am(i,k)
  if (rik < ric)  then
     if(bn2 < bn2min) then
       bn_uns = bn2
       bn2 = bn2min
	if (bn_uns < 0) then 
	 ins_con(k) = 1
	 num_ins = num_ins + 1
	endif 
     endif  
  endif   
    zgrow = exp(zm(i,k)*Ht)   
    lturbmax = dzg*0.25     
    zlturb = min(lturb*zgrow, lturbmax)
    
           kamp = sqrt(shr2)* zlturb*zlturb
	   ritur = bn2/shr2 *dzg *dzg
           w1 = 1./(1. + 5.*ritur)
	   w1 = w1*w1  
	   	   
!          w1 = 1./(1. + 10.*ritur)/(1.+8.*ritur)  !HB-93 
	            
       kzz_edd(k)= min(max(kamp * w1, dked_min), dked_max)  
  enddo
  
      nfixed = nfixed + num_ins  
!
! check con_adjust_t1   nrest
!      
!    do k=2, km
!       dzg =  zm(i,k) - zm(i,k+1) 
!       dtemp =t2(i,k)     - t2(i,k+1)
!       tz =  dtemp /dzg 
!       tc = .5*(t2(i,k)   + t2(i,k+1))  
!       gcp = gvz(i,k)/cp_mu(i,k)   
!       if ((tz +gcp).lt.0.) nrest= nrest+1       
!    enddo
     
    do k=2,km-1
        w1 = .25*( kzz_edd(k-1)+2.*kzz_edd(k)+kzz_edd(k+1) )     
        ktedd(i,k) =w1
        vuedd(i,k)= w1
        dfedd(i,k)= w1
    enddo	 
 enddo         !i-loop
 
 end subroutine get_coef_turb
 
 subroutine con_adjust_t2( is, ie, km, t2, dp2, pe2, cap2)
 
 
 integer            :: km, is, ie 
 real, intent(in)   ::  dp2(is:ie, km), pe2(is:ie, km+1), cap2(is:ie, km)
 real, intent(inout)::  t2(is:ie, km)
  
 integer :: i, k, kk, ktrop, jiter, niter
 
 real, dimension(km)   ::  c1dad,  c2dad, c3dad,  c4dad
 real :: zeps, rdenom, kappa, dpmid, zepsdp, zgamma, dtnext, tc1dad 
 real :: tins, atins
 logical :: stable
     ktrop = 72
     zeps = 2.0e-4
     niter= 1
     
  do 80 i=is, ie   
       stable = .true.
       
 77    do   jiter=1, niter     
   do k= 2, ktrop 
   
         kappa =.25*(cap2(i, k)+ cap2(i,k+1))
	 dpmid = .5*(pe2(i,k+2)-pe2(i,k))
         c1dad(k) = kappa*dpmid/pe2(i,k+1)
	 c2dad(k) = (1. - c1dad(k))/(1. + c1dad(k))
	 rdenom = 1./(dp2(i,k)*c2dad(k) + dp2(i,k+1))
	 c3dad(k) = rdenom*dp2(i,k)
	 c4dad(k) = rdenom*dp2(i,k+1)
	 zepsdp = zeps*dpmid
	 
!c1dad = kappa*dpmid/pe2(k+1)
	 	 
	 zgamma = c1dad(k)*(t2(i,k) + t2(i,k+1)) 
	   tins = t2(i,k+1)-t2(i,k)
	   if (tins .ge.  (zgamma+zepsdp)) then	
	     stable = .false. 
!	     print *, ' dri_adj1 ', k, tins, t2(i,k+1), t2(i,k), zgamma, zepsdp	     
	     t2(i,k+1) = t2(i,k)*c3dad(k) + t2(i,k+1)*c4dad(k)
	     t2(i,k) =   c2dad(k)*t2(i,k+1)  
!	     atins = t2(i,k+1)-t2(i,k)
!	     print *, ' dri_adj2 ', k, atins, t2(i,k+1), t2(i,k)
	  endif  
	    
      enddo
         if(stable) goto 80  ! next i-loop
!	 zeps = zeps +zeps
!              ! next iteration to check stability
      enddo    
      
  80 continue  
  
 end subroutine con_adjust_t2 
 
!====================================================
!q, pt, w, u, v, dp2, pe2, dpu, dpv,   
!
! implicit driver for the vertical M-E dissipation
! + convective adjustment "dry"
!====================================================
 subroutine  get_moldiff(dtin, vumol, ktmol, dfmol, rhomol, wgrav,   &
                         q, t, w, u, v, dp2, pe2, dpu, dpv, peu, pev, cap3, &
                         j, je, is, ie, isd, ied, jsd, jed, km, nq,  &
	                 ind_h2o, ind_o2, ind_o3p, ind_o3 )
 implicit none
 logical  :: con_adj			 
 integer, intent(in)   :: j, je		  
 integer, intent(in)   :: is, ie, isd,ied, jsd,jed, km, nq
 integer, intent(in)   ::  ind_h2o, ind_o2, ind_o3p, ind_o3 
 
 real, intent(inout)::  q(isd :ied,jsd:jed, km, nq) 
 real, intent(inout)::  t(isd :ied,jsd:jed, km)
 real, intent(in)   ::  cap3(isd :ied,jsd:jed, km)
 real, intent(inout)::  w(isd :ied,jsd:jed, km) 
 real, intent(inout)::  u(isd :ied,jsd:jed+1, km)  
 real, intent(inout)::  v(isd :ied+1,jsd:jed, km)   
 real, intent(in)   ::  dp2(is:ie, km), pe2(is:ie, km+1)
 real, intent(in)   ::  dpu(is:ie, km),dpv(is:ie+1, km)
 real, intent(in)   ::  peu(is:ie, km+1),pev(is:ie+1, km+1) 
 real, intent(in)   ::  dtin 

!   real, intent(inout)::  u(isd:ied  ,jsd:jed+1,km)   !< u-wind (m/s)-DDEBUG=ON
!   real, intent(inout)::  v(isd:ied+1,jsd:jed  ,km)   !< v-wind (m/s)
 real, intent(inout), dimension(is :ie, km+1) ::  vumol, ktmol, dfmol, rhomol, wgrav
 character(len=128) :: file_dis 
 logical                       :: mvis_debug
 real, dimension(is :ie, km+1) :: vuedd, ktedd, dfedd
 real, dimension(is :ie, km+1) :: cp_mu, zgeo
 real, dimension(is :ie, km)   :: am_mol, kion 
 real :: dt_ic, dt 
 real :: Ne_prof(km)
 real :: nef2, hef2, hch_f2, nede, hkde, hpde, bdde, bede,zchap
!
!
! 
 real :: q2(is : ie, km), t2(is : ie, km), ut(is:ie, km), vt(is:ie, km)
 real :: pm2(is : ie, km), cap2(is : ie, km) 
 real, dimension(is:ie,km) :: qo, qo2, qh2o, qo3 
 real, dimension(is:ie+1,km) :: v2
 real, dimension(is:ie,km)   :: v2c  
 integer :: itr, i, k, ii, kk, ju, nfixed, nrest
!
! km  1/km and cm-3
! 
!  dt_ic = 180.
  dt = dtin
  
  con_adj   =.false.
  mvis_debug=.false.
  
 nef2 = 7.10264e5
 hef2 = 300.e3
 hch_f2 = 39.e3
 nede = 1.43e13
 hkde = 95.e3
 hpde = 70.e3
 bdde = 0.65e-3
 bede =0.14e-3
 if (j == je+1) then 
    ju = je
 else
    ju = j
 endif    
 
 qo(is : ie, 1:km)  = q(is:ie, ju, 1:km, ind_o3p) 
 qo2(is : ie, 1:km) = q(is:ie, ju, 1:km, ind_o2) 
 qh2o(is : ie, 1:km)= q(is:ie, ju, 1:km, ind_h2o) 
 qo3(is : ie, 1:km) = q(is:ie, ju, 1:km, ind_o3) 
 

 wgrav(is:ie, 1:km+1) = 9.8065
 
 do ii = is, ie
 do kk = 1,km
   t2(ii, kk) = t(ii, ju, kk) 
   ut(ii, kk) = u(ii, j, kk)
   vt(ii, kk) = .5*(v(ii, ju, kk) +  v(ii+1, ju, kk))
   cap2(ii,kk) =cap3(ii, ju, kk)
   pm2(ii, kk) = .5*(pe2(ii,kk+1)+pe2(ii,kk))
!   (pe2(ii,kk+1)-pe2(ii,kk))/log(pe2(ii,kk+1)/pe2(ii,kk))

 enddo
  
 enddo  
  
 if (con_adj) call con_adjust_t2( is, ie, km, t2, dp2, pe2, cap2) 
 
 
 call get_coef_mdif(is, ie, km, t2, qo, qo2, qo3, qh2o, dp2, pe2, wgrav,  vumol, ktmol, dfmol, rhomol, &
        cp_mu, am_mol, zgeo)
	
 do ii = is, ie
  do kk = km, 1,-1  
   zchap = (zgeo(ii,kk)-hef2)/hch_f2 
   Ne_prof(kk) = nef2*exp(0.5*(1- zchap -exp(-zchap)))
   kion(ii,kk) = 7.22e-11* (t2(ii, kk) )** 0.37 *Ne_prof(kk) 
!   kion(ii,kk) = 0.
 enddo 
 enddo	
	
 call get_coef_turb(is, ie, km, t2, ut, vt,  dp2, pm2, pe2, wgrav, cp_mu, am_mol, zgeo, vuedd, ktedd, dfedd,&
                   con_adj, nfixed, nrest) 
   if((mvis_debug) .and. (is_master())) then
!      print *, 'isd, ied, jsd, jed ', isd, ied, jsd, jed
!      print *, ' is, ie, km ',    is, ie, km 
     file_dis = trim('zgeo'//'_dedug.form') 
     open(unit=77, file=file_dis, status='unknown', form='formatted')
     write(77,*) is
     write(77,*) ie
     write(77,*) isd, ied, jsd, jed
     write(77,*) km
     write(77,*) vumol
     write(77,*) vuedd  
     write(77,*) ktmol
     write(77,*) kion     
     write(77,*) zgeo
     write(77,*) cp_mu
     write(77,*) am_mol
     write(77,*) wgrav 
     write(77,*) ut
     write(77,*) vt 
     write(77,*) t2  
     write(77,*)  dp2
     write(77,*)  pm2
     write(77,*)  pe2  
     close(77)                              
   endif  		   
  vumol =vumol +  vuedd
  ktmol =ktmol +  ktedd*0.0
  dfmol =dfmol +  dfedd*0.0
  
 if (is_master() .and. j == je) then
!      print *, ' mdif-kion ', maxval(kion), 1./maxval(kion)*86400.
!      print *, ' mdif-zgeo/grav ', maxval(zgeo)*1.e-3, minval(wgrav)
      print *, ' mdif-ktedd ', maxval(ktedd), minval(ktedd)   
      print *, ' mdif-nfixed ', nfixed, nrest                 
 endif

    if (j == je+1) then
    
    ut(is:ie, 1:km) = u(is:ie, j, 1:km) 
     call uvget_molviscosity(dt, kion, vumol, rhomol, dpu(is:ie,:), peU(is:ie,:), wgrav, is, ie, isd, ied, jsd, jed, km, ut, 'Udis') 
    u(is:ie, j, 1:km) = ut(is:ie, 1:km)     
    
    else 
!
! Update all variables j <=je
!    
    call get_molviscosity(dt, ktmol, rhomol, dp2, pe2, wgrav, is, ie, isd, ied, jsd, jed, km, t2, 'Tdis') 
 do ii = is, ie
   do kk = 1,km
!
! consider to apply the dry-adjustment for t2(ii, kk)
! 
   t(ii, j, kk) = t2(ii, kk) 
!  if (t2(ii,kk) < 90. )  print *, 'fv_molvisZ Tkin', t2(ii,kk), kk
 enddo
 enddo 
!w    
 t2(is:ie, 1:km) = w(is:ie, j, 1:km)
 call uvget_molviscosity(dt, kion, vumol, rhomol, dp2, pe2, wgrav, is, ie, isd, ied, jsd, jed, km, t2, 'Wdis') 
 w(is:ie, j, 1:km) = t2(is:ie, 1:km)
 
! u
 t2(is:ie, 1:km) = u(is:ie, j, 1:km) 
 call uvget_molviscosity(dt, kion, vumol, rhomol, dpu(is:ie,:), peU(is:ie,:), wgrav, is, ie, isd, ied, jsd, jed, km, t2, 'Udis') 
 u(is:ie, j, 1:km) = t2(is:ie, 1:km) 
! 
! v
 v2c(is:ie, 1:km) = v(is:ie, j, 1:km) 
 call uvget_molviscosity(dt, kion, vumol, rhomol, dpv(is:ie,:), peV(is:ie,:), wgrav, &
       is, ie, isd, ied, jsd, jed, km, v2c, 'Vdis') 
 v(is:ie, j, 1:km) = v2c(is:ie, 1:km) 
 v(ie+1, j, 1:km) = .5*(v2c(ie, 1:km) +v(ie+1, j, 1:km))
 
 call get_molviscosity(dt, dfmol,rhomol, dp2, pe2, wgrav, is, ie, isd, ied, jsd, jed, km, qo, 'Opdis')
 q(is:ie, j, 1:km, ind_o3p)   = qo(is : ie, 1:km)
 
 call get_molviscosity(dt, dfmol, rhomol,dp2, pe2, wgrav, is, ie, isd, ied, jsd, jed, km, qo2,'O2dis')
 q(is:ie, j, 1:km, ind_o2)   = qo2(is : ie, 1:km) 
 
 call get_molviscosity(dt, dfmol,rhomol, dp2, pe2, wgrav, is, ie, isd, ied, jsd, jed, km, qo3, 'O3dis')
 q(is:ie, j, 1:km, ind_o3)   = qo3(is : ie, 1:km) 
  
 call get_molviscosity(dt, dfmol, rhomol,dp2, pe2, wgrav, is, ie, isd, ied, jsd, jed, km, qh2o, 'H2dis')
 q(is:ie, j, 1:km, ind_h2o)   = qh2o(is : ie, 1:km)
! print *, 'get_moldiff' 
 
 endif    ! j=je+1
  
 end subroutine get_moldiff
!=============================
 subroutine  uvget_molviscosity(dt, kion, vum, rhomol, dp2, pe, wgrav, is, ie, isd, ied, jsd, jed, km, u,strdis)
 	
  implicit none
  integer, intent(in)   ::  is, ie, isd,ied, jsd,jed, km
  real, intent(in)      ::  dt
  real, intent(inout)   ::  u(is :ie, km)
  real, intent(in)      ::   kion(is:ie, km)
  real, intent(in)      ::  vum(is:ie, km+1), dp2(is:ie, km), wgrav(is:ie, km+1), rhomol(is:ie, km+1)
  real, intent(in)      ::  pe(is:ie, km+1)
  character(len=*),  intent(in)      ::   strdis  
  integer, parameter :: dir_fac = 1
  
  integer  :: i, k
  real :: ak(km), bk(km),  ck(km), dk(km), ed(km), fd(km)
  real :: dta, rdt, dpc, dpf, gpc, gpf, rh2c, rh2a, tint, rinv  
  real :: uu(km)
  logical :: mvis_debug
  character(len=128) :: file_dis 
!  
!  
   mvis_debug = .false.
   dta = abs(dt)
   rdt = 1./dta
   
   if((mvis_debug) .and. (is_master())) then
!      print *, 'isd, ied, jsd, jed ', isd, ied, jsd, jed
!      print *, ' is, ie, km ',    is, ie, km 
     file_dis = trim(strdis//'_dedug.form') 
     open(unit=77, file=file_dis, status='unknown', form='formatted')
     write(77,*) is
     write(77,*) ie
     write(77,*) isd, ied, jsd, jed
     write(77,*) km
     write(77,*) vum
     write(77,*) rhomol
     write(77,*) dp2
     write(77,*) pe
     write(77,*) wgrav
     write(77,*) U    
   endif
  do i= is, ie
    ak(1) = 0.
    ck(km) =0.
    dpf = dp2(i,1)
    dpc = pe(i,1)
    gpf = .5*(wgrav(i,1)+wgrav(i,2))
    
    ck(1) = dta*vum(i,2)*gpf*wgrav(i,2)/dpf*rhomol(i,2)
    bk(1) = ck(1) + 1. +kion(i,1)*dta
    
    dk(1) = u(i,1)
    do k =2, km-1
       ak(k) = ck(k-1)
       dk(k) = u(i,k)
       gpf = .5*(wgrav(i,k)+wgrav(i,k+1))*vum(i,k+1)      
       ck(k) = dta*gpf/dp2(i,k)*rhomol(i,k) 
       bk(k) = ak(k) + ck(k) + 1.+kion(i,k)*dta
    enddo
    
       k= km
       ak(k) = ck(k-1)
       dk(k) = u(i,k)
       bk(k) = 1. +ak(k) +kion(i,k)*dta
    if (dir_fac == 1 ) then  
       ed(k) =  ak(k)/bk(k)
       fd(k)=   dk(k)/bk(k)
!
!   Y[k] = ed*Y[k-1] +fd[k]  ed = a(km)/b(km) fd = D(km)/b(km)
!       
       do k = km-1, 1, -1
          rinv = 1./(bk(k)-ck(k)*ed(k+1))
          ed(k) = ak(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k+1))*rinv
       enddo
!
! k=1, top-layer
!       
          u(i,1) = fd(1)
! from top to ground	  
       do k=2,km
          u(i,k) = ed(k)*u(i, k-1) + fd(k)
       enddo
    else                    !========== direction from surf to top
       ed(1) =  ck(1)/bk(1)
       fd(1)=   dk(1)/bk(1) 
       do k = 2, km, 1
          rinv = 1./(bk(k)-ak(k)*ed(k-1))
          ed(k) = ck(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k-1))*rinv
       enddo 
          u(i,km) = fd(km)     !   the bottom-layer      
       do k=km-1,1,-1
          u(i,k) = ed(k)*u(i,k+1) + fd(k)
       enddo                  
    endif    ! directional factorization: from [top=>bot] or [bot=>top]        
  enddo
!   
!         if(mvis_debug) call mpp_error(FATAL,"ERROR MDIF non positive value of plyr") 
   if((mvis_debug) .and. (is_master())) then
    write(77,*) ak
    write(77,*) bk  
    write(77,*) ck
    write(77,*) dk
    write(77,*) ed
    write(77,*) fd    
    write(77,*) U   
    write(77,*) dt     
    close(77)   
!    call mpp_error(FATAL,"STOP in get_molviscosity -DEBUG")                  
   endif  
 end subroutine  uvget_molviscosity
 
!=============================================================
! 
 subroutine  get_molviscosity(dt, vum, rhomol, dp2, pe, wgrav, is, ie, isd, ied, jsd, jed, km, u,strdis)
 	
  implicit none
  integer, intent(in)   ::  is, ie, isd,ied, jsd,jed, km
  real, intent(in)      ::  dt
  real, intent(inout)   ::  u(is :ie, km) 
  real, intent(in)      ::  vum(is:ie, km+1), dp2(is:ie, km), wgrav(is:ie, km+1), rhomol(is:ie, km+1)
  real, intent(in)      ::  pe(is:ie, km+1)
  character(len=*),  intent(in)      ::   strdis  
!
! compute coef for dU/dt=1/rho*d [(rho*K) dU/dz] /dz
!
!   rho*dz = -dp/g   dU/dt = g d[ K rho^2 *g (du/dp)]dp,  rho = Pe/(RT)

  integer, parameter :: dir_fac = 1
  
  integer  :: i, k
  real :: ak(km), bk(km),  ck(km), dk(km), ed(km), fd(km)
  real :: dta, rdt, dpc, dpf, gpc, gpf, rh2c, rh2a, tint, rinv  
  real :: uu(km)
  logical :: mvis_debug
  character(len=128) :: file_dis 
!  
!  
   mvis_debug = .false.
   dta = abs(dt)
   rdt = 1./dta
   
   if((mvis_debug) .and. (is_master())) then
!      print *, 'isd, ied, jsd, jed ', isd, ied, jsd, jed
!      print *, ' is, ie, km ',    is, ie, km 
     file_dis = trim(strdis//'_dedug.form') 
     open(unit=77, file=file_dis, status='unknown', form='formatted')
     write(77,*) is
     write(77,*) ie
     write(77,*) isd, ied, jsd, jed
     write(77,*) km
     write(77,*) vum
     write(77,*) rhomol
     write(77,*) dp2
     write(77,*) pe
     write(77,*) wgrav
     write(77,*) U    
   endif
  do i= is, ie
    ak(1) = 0.
    ck(km) =0.
    dpf = dp2(i,1)
    dpc = pe(i,1)
    gpf = .5*(wgrav(i,1)+wgrav(i,2))
    
    ck(1) = dta*vum(i,2)*gpf*wgrav(i,2)/dpf*rhomol(i,2)
    bk(1) = ck(1) + 1.
    
    dk(1) = u(i,1)
    do k =2, km-1
       ak(k) = ck(k-1)
       dk(k) = u(i,k)
       gpf = .5*(wgrav(i,k)+wgrav(i,k+1))*vum(i,k+1)      
       ck(k) = dta*gpf/dp2(i,k)*rhomol(i,k) 
       bk(k) = ak(k) + ck(k) + 1.
    enddo
    
       k= km
       ak(k) = ck(k-1)
       dk(k) = u(i,k)
       bk(k) = 1. +ak(k) 
    if (dir_fac == 1 ) then  
       ed(k) =  ak(k)/bk(k)
       fd(k)=   dk(k)/bk(k)
!
!   Y[k] = ed*Y[k-1] +fd[k]  ed = a(km)/b(km) fd = D(km)/b(km)
!       
       do k = km-1, 1, -1
          rinv = 1./(bk(k)-ck(k)*ed(k+1))
          ed(k) = ak(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k+1))*rinv
       enddo
!
! k=1, top-layer
!       
          u(i,1) = fd(1)
! from top to ground	  
       do k=2,km
          u(i,k) = ed(k)*u(i, k-1) + fd(k)
       enddo
    else                    !========== direction from surf to top
       ed(1) =  ck(1)/bk(1)
       fd(1)=   dk(1)/bk(1) 
       do k = 2, km, 1
          rinv = 1./(bk(k)-ak(k)*ed(k-1))
          ed(k) = ck(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k-1))*rinv
       enddo 
          u(i,km) = fd(km)     !   the bottom-layer      
       do k=km-1,1,-1
          u(i,k) = ed(k)*u(i,k+1) + fd(k)
       enddo                  
    endif    ! directional factorization: from [top=>bot] or [bot=>top]        
  enddo
!   
!         if(mvis_debug) call mpp_error(FATAL,"ERROR MDIF non positive value of plyr") 
   if((mvis_debug) .and. (is_master())) then
    write(77,*) ak
    write(77,*) bk  
    write(77,*) ck
    write(77,*) dk
    write(77,*) ed
    write(77,*) fd    
    write(77,*) U   
    write(77,*) dt     
    close(77)   
!    call mpp_error(FATAL,"STOP in get_molviscosity -DEBUG")                  
   endif  
 end subroutine  get_molviscosity
!=============================================================
!
! 
  subroutine get_molviscosity_ie1(dt, vum, rhomol, dp2, pe, wgrav, is, ie, isd, ied, jsd, jed, km, u, strdis) 
 implicit none
 integer, intent(in)   ::  is, ie, isd,ied, jsd,jed, km
 real, intent(in)      ::  dt
 real, intent(inout)   ::  u(is : ie+1, km) 
 real, intent(in)      ::  pe(is:ie+1, km+1), dp2(is:ie+1, km) 
 real, intent(in)      ::  vum(is:ie, km+1),  wgrav(is:ie, km+1), rhomol(is:ie, km+1)
 
 character(len=*),  intent(in)      ::   strdis  

  integer  :: i, k, ii
  real :: ak(km), bk(km),  ck(km), dk(km), ed(km), fd(km)
  real :: rdt, dpc, dpf, gpc, gpf, rh2c, rh2a, tint, rinv  
  real :: uu(km), dta
  logical :: mvis_debug
  integer, parameter :: dir_fac = 1
  character(len=128) :: file_dis 
!  
!  
   mvis_debug = .false.
   dta = abs(dt)
   rdt = 1./dta   
   if((mvis_debug) .and. (is_master())) then

     file_dis = trim(strdis//'_dedug.form') 
     open(unit=77, file=file_dis, status='unknown', form='formatted')
     write(77,*) is
     write(77,*) ie
     write(77,*) isd, ied, jsd, jed
     write(77,*) km
     write(77,*) vum
     write(77,*) rhomol
     write(77,*) dp2(is:ie, 1:km)  
     write(77,*) pe(is:ie, 1:km+1) 
     write(77,*) wgrav
     write(77,*) U(is:ie, 1:km)    
   endif
  do i= is, ie
    ak(1) = 0.
    ck(km) =0.
    dpf = dp2(i,1)
    dpc = pe(i,1)
    gpf = .5*(wgrav(i,1)+wgrav(i,2))
    
    ck(1) = dta*vum(i,2)*gpf*wgrav(i,2)/dpf*rhomol(i,2)
    bk(1) = ck(1) + 1.
    dk(1) = u(i,1)* dta
    
   
    do k =2, km-1
       ak(k) = ck(k-1)
       dk(k) = u(i,k)*dta
       gpf = .5*(wgrav(i,k)+wgrav(i,k+1))*vum(i,k+1)      
       ck(k) = dta*gpf/dp2(i,k)*rhomol(i,k) 
       bk(k) = ak(k) + ck(k) + 1.
    enddo
    
       k= km
       ak(k) = ck(k-1)
       dk(k) = u(i,k)*dta
       bk(k) = 1. +ak(k)
       
    if (dir_fac == 1 ) then  
       
       ed(k) =  ak(k)/bk(k)
       fd(k)=   dk(k)/bk(k)
!
!   Y[k] = ed*Y[k-1] +fd[k]  ed = a(km)/b(km) fd = D(km)/b(km)
!       
       do k = km-1, 1, -1
          rinv = 1./(bk(k)-ck(k)*ed(k+1))
          ed(k) = ak(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k+1))*rinv
       enddo    
          u(i,1) = fd(1)          !    k=1, top-layer
! from top to ground	  
       do k=2,km
          u(i,k) = ed(k)*u(i,k-1) + fd(k)
       enddo
!       
    else               !dirfac - 1
!    
       ed(1) =  ck(1)/bk(1)
       fd(1)=   dk(1)/bk(1) 
       do k = 2, km, 1
          rinv = 1./(bk(k)-ak(k)*ed(k-1))
          ed(k) = ck(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k-1))*rinv
       enddo 
          u(i,km) = fd(km)     !   the bottom-layer      
       do k=km-1,1,-1
          u(i,k) = ed(k)*u(i,k+1) + fd(k)
       enddo                  
    endif    ! directional factorization: from [top=>bot] or [bot=>top]  
  enddo
!   dt, vum, rhomol, dp2, pe, wgrav, is, ie, isd, ied, jsd, jed, km, u  
!         if(mvis_debug) call mpp_error(FATAL,"ERROR MDIF non positive value of plyr") 
   if((mvis_debug) .and. (is_master())) then
    write(77,*) ak
    write(77,*) bk  
    write(77,*) ck
    write(77,*) dk
    write(77,*) ed
    write(77,*) fd    
    write(77,*) U(is:ie,1:km)   
    write(77,*) dt     
    close(77)   
!    call mpp_error(FATAL,"STOP in get_molviscosity -DEBUG")                  
   endif  
!
! extra-point for V-winds
! 
   ii = ie  
   i = ie+1
   
    ak(1) = 0.
    ck(km) =0.
      
    dpf = dp2(i,1)
    dpc = pe(i,1)
    gpf = .5*(wgrav(ii,1)+wgrav(ii,2))
    
    ck(1) = dta*vum(ii,2)*gpf*wgrav(ii,2)/dpf*rhomol(ii,2)
    bk(1) = ck(1) + 1.
    dk(1) = u(i,1)* dta
    
    do k =2, km-1
       ak(k) = ck(k-1)
       dk(k) = u(i,k)*dta
       gpf = .5*(wgrav(ii,k)+wgrav(ii,k+1))*vum(ii,k+1)      
       ck(k) = dta*gpf/dp2(i,k)*rhomol(ii,k) 
       bk(k) = ak(k) + ck(k) + 1.
    enddo 
       
    k= km
       ak(k) = ck(k-1)
       dk(k) = u(i,k)* dta
       bk(k) = 1. +ak(k) 
  if (dir_fac == 1 ) then  
       ed(k) =  ak(k)/bk(k)
       fd(k)=   dk(k)/bk(k)         
    do k = km-1, 1, -1
          rinv = 1./(bk(k)-ck(k)*ed(k+1))
          ed(k) = ak(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k+1))*rinv
    enddo   
        u(i,1) = fd(1) 
    do k=2,km
          u(i,k) = ed(k)*u(i,k-1) + fd(k)
    enddo 
    
    else          ! directional factorization: from [top=>bot] or [bot=>top] 
       ed(1) =  ck(1)/bk(1)
       fd(1)=   dk(1)/bk(1) 
       do k = 2, km, 1
          rinv = 1./(bk(k)-ak(k)*ed(k-1))
          ed(k) = ck(k) *rinv 
	  fd(k) = (dk(k)+ck(k)*fd(k-1))*rinv
       enddo 
          u(i,km) = fd(km)     !   the bottom-layer      
       do k=km-1,1,-1
          u(i,k) = ed(k)*u(i,k+1) + fd(k)
       enddo        
    endif          
 end subroutine  get_molviscosity_ie1
! ========================================
!     
end module fv_mapz_mod
