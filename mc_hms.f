         subroutine mc_hms (p_spec, th_spec, dpp, x, y, z, dxdz, dydz,
     > x_fp, dx_fp, y_fp, dy_fp, tg_rad_len,m2,cos_ts,
     > sin_ts,spect, ms_flag, wcs_flag, decay_flag,
     > resmult, fry,frx,ok_hms)

C+______________________________________________________________________________
!
! Monte-Carlo of HMS spectrometer.
!  Note that we only pass on real*8 variables to the subroutine.
!  This will make life easier for portability.
!
! Author: David Potterveld, March-1993
!
! Modification History:
!
!  11-Aug-1993 	(D. Potterveld) Modified to use new transformation scheme in
!    which each transformation begins at the pivot.
!
!  19-AUG-1993  (D. Potterveld) Modified to use COSY INFINITY transformations.
!
!  15-SEP-1997  MODIFY stepping through spectrometer so that all drifts
!    use project.f (not transp.f), and so that project.f and
!    transp.f take care of decay.  Decay distances all assume
!    the pathlength for the CENTRAL RAY.
!
!  06-MAR-2008 (Anusha) Modified the Codes to use the Target Field Corrections.
!               Added extra subroutines track_from_tgt.f to project the particle
!               from the target to the field free region through the field (100 
!               cm from the target) and track_to_tgt.f to reconstruct the target 
!               quentities (together with simc_hms_recon.f)
C-______________________________________________________________________________

        implicit none

        include 'apertures.inc'
        include 'hms_mc.cmn'
        include 'track.inc'
        include 'constants.inc'

*G.&I.N. STUFF - for checking dipole exit apertures and vacuum pipe to HMS hut.
        include 'g_dump_all_events.inc'
 
        real*8 x_offset_pipes/0.d0/,y_offset_pipes/0.d0/
        real*8 x_offset_hut/3.50/,y_offset_hut/0.60/


! Spectrometer definitions

        integer*4 hms,sos
        parameter (hms = 1)
        parameter (sos = 2)

! Collimator (octagon) dimensions and offsets.

        real*8 h_entr,v_entr,h_exit,v_exit
        real*8 x_off,y_off,z_off,xtemp,ytemp

c        real*8 phi_e,phi_p
        real*8 bdl,xtgt,frx  ! frx beam right
        real*8 cos_ts,sin_ts
c        character*10 map 
        integer*4 spect

        ! Open values (no collimator)
c       parameter (h_entr = 20.0)
c       parameter (v_entr = 20.0)
c       parameter (h_exit = 20.0)
c       parameter (v_exit = 20.0)
 
        ! Old, 'large' or 'pion' collimator.
c        parameter (h_entr = 3.536)
c        parameter (v_entr = 9.003)
c        parameter (h_exit = 3.708)
c        parameter (v_exit = 9.444)


        ! New collimator for HMS-100 tune.
        parameter (h_entr = 4.560)
        parameter (v_entr = 11.646)
        parameter (h_exit = 4.759)
        parameter (v_exit = 12.114)        


c        parameter (x_off=0.057)        !+ve is slit DOWN      !!!! standard
        parameter (x_off=-0.043)        !+ve is slit DOWN      !!!! standard
        parameter (y_off=+0.030)        !+ve is slit LEFT (as seen from target)
        parameter (z_off=+40.17)      !HMS100 tune (dg 5/27/98)


! Math constants

        real*8 d_r,r_d,root
  
        parameter (d_r = pi/180.)
        parameter (r_d = 180./pi)
        parameter (root = 0.707106781)  !square root of 1/2


! The arguments

        real*8 x,y,z      !(cm)
        real*8 dpp        !delta p/p (%)
        real*8 dxdz,dydz  !X,Y slope in spectrometer
        real*8 x_fp,y_fp,dx_fp,dy_fp !Focal plane values to return
        real*8 p_spec,th_spec  !spectrometer setting
        real*8 tg_rad_len    !target length in r.l.
        real*8 fry           !vertical position@tgt (+y=up)
        logical ms_flag      !mult. scattering flag
        logical wcs_flag     !wire chamber smearing flag
        logical decay_flag   !check for particle decay
        logical ok_hms       !true if particle makes it

! Local declarations.

        integer*4  chan/1/,n_classes

cc       logical	first_time_hms

        real*8 dpp_recon,dth_recon,dph_recon !reconstructed quantities
        real*8 y_recon
        real*8 p,m2 !More kinematic variables.
        real*8 xt,yt !temporaries
        real*8 resmult !DC resolution factor

        logical dflag  !has particle decayed?
        logical ok

! Gaby's dipole shape stuff
        logical checkdip,checksieve
        logical check_dipole,check_sieve
        external check_dipole
        external check_sieve

        logical fieldon
        common /mag/fieldon

        save  !Remember it all!


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cc      first_time_hms = .true.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C================================ Executable Code =============================

! Initialize ok_hms to .false., reset decay flag

        ok_hms = .false.
        dflag = .false.   !particle has not decayed yet
        trials = trials + 1
c write(*,*)'trials',trials

! Save spectrometer coordinates.

        xs = x
        ys = y
        zs = z
        dxdzs = dxdz
        dydzs = dydz


c  write(6,*) x_offset_pipes,x_offset_hut


! particle momentum

        dpps = dpp
        p = p_spec*(1.+dpps/100.)

! Read in transport coefficients.

         if (first_time_hms) then   
           call transp_init(hms,n_classes)
          close (unit=chan)
          if (n_classes.ne.12) stop 'MC_HMS, wrong number of transport classes'
          first_time_hms = .false.
        endif

! Calculate multiple scattering in target
       
         if(ms_flag) call musc(m2,p,tg_rad_len,dydzs,dxdzs)

c if(ms_flag) call musc_ld(m2,p,6.0,12.0,0.41,dydzs,dxdzs)  !!!  check with C12

! Begin transporting particle.
! use target field corrections to the out going particle

**********************************************************************************
! code tracks particle to z=100 
! use the momentum in MeV

cc call track_from_tgt(xs,ys,zs,dxdzs,dydzs,p*1000.,m2,spect,ok_hms)
c  if (.not. ok_spec) then
c  xs = xs - zs * dxdzs
c  ys = ys - zs * dydzs
c  zs = 0.0d00
c
c  endif

c write(*,*) ' after ',xs,ys,zs,dxdzs,dydzs,p*1000.,m2,spect,ok_hms
c read(*,*) ok_hms
**********************************************************************************

! Do transformations, checking against apertures.      
! Check sieve slit holes if dosieve = true  !(1.68 m from the target)

          if(dosieve) then          
           xtemp = xs
           ytemp = ys
           call project(xtemp,ytemp,168.d0,decay_flag,dflag,m2,p) 
           checksieve = check_sieve(xtemp,ytemp)
           if(.not.checksieve) then
            where=25.
            x_stop=xtemp
            y_stop=ytemp
            goto 500
           endif
           endif

! Check front of fixed slit, at about 0.26 meter(actually 1.26 m from the target and 1.0 m consider as target field)

           call project(xs,ys,126.2d0+z_off,decay_flag,dflag,m2,p) !project and decay
           if (abs(ys-y_off).gt.h_entr) then
           slit_hor = slit_hor + 1
           where=1.
           x_stop=xs
           y_stop=ys
c write(*,*)'slit_hor',slit_hor,ys,y_off,ys-y_off,h_entr
           goto 500
c
c           if(where.EQ.1.)write(6,*) where

           endif
           if (abs(xs-x_off).gt.v_entr) then
           slit_vert = slit_vert + 1
           where=2.     
           x_stop=xs
           y_stop=ys
c write(*,*)'slit_vert',slit_vert
           goto 500
           endif
            if (abs(xs-x_off).gt.
     & (-(v_entr/h_entr*abs(ys-y_off))+3.d0*v_entr/2)) then
           slit_oct = slit_oct + 1
           where=3.
           x_stop=xs
           y_stop=ys
c write(*,*)'slit_oct',slit_oct
           goto 500

           endif
!Check back of fixed slit, at about 1.325 meter(1.26 m+0.06 m)
c  shut = shut + 1
c  write(*,*)'shut',shut

          call project(xs,ys,6.3d0,decay_flag,dflag,m2,p) !project and decay
          if (abs(ys-y_off).gt.h_exit) then
          slit_hor = slit_hor + 1
           where=4.
           x_stop=xs
           y_stop=ys
           goto 500
          endif
          if (abs(xs-x_off).gt.v_exit) then
          slit_vert = slit_vert + 1
           where=5.
           x_stop=xs
          y_stop=ys
          goto 500
          endif
          if (abs(xs-x_off).gt.
     &           (-(v_exit/h_exit*abs(ys-y_off))+3.d0*v_exit/2)) then
          slit_oct = slit_oct + 1
           where=6.
           x_stop=xs
           y_stop=ys
           goto 500
           endif

! Go to Q1 IN  mag bound.  Drift rather than using COSY matrices

          call project(xs,ys,(216.075d0-126.2d0-z_off-6.3d0),
     &               decay_flag,dflag,m2,p) !project and decay
          if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
           Q1_in = Q1_in + 1
           where=7.     
           x_stop=xs
           y_stop=ys
           goto 500
            endif

! Check aperture at 2/3 of Q1.

          call transp(hms,2,decay_flag,dflag,m2,p,126.0d0)
          if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
           Q1_mid = Q1_mid + 1
           where=8.
           x_stop=xs
           y_stop=ys
          goto 500
          endif

! Go to Q1 OUT mag boundary.

          call transp(hms,3,decay_flag,dflag,m2,p,63.0d0)
         if ((xs*xs + ys*ys).gt.r_Q1*r_Q1) then
           Q1_out = Q1_out + 1
           where=9.  
           x_stop=xs
           y_stop=ys
          goto 500
          endif

! Go to Q2 IN  mag bound.  Drift rather than using COSY matrices
!! call transp(hms,4,decay_flag,dflag,m2,p)

          call project(xs,ys,123.15d0,decay_flag,dflag,m2,p) !project and decay
          if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
           Q2_in = Q2_in + 1
           where=10.
           x_stop=xs
           y_stop=ys
           goto 500
           endif

! Check aperture at 2/3 of Q2.

          call transp(hms,5,decay_flag,dflag,m2,p,143.67d0)
          if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
           Q2_mid = Q2_mid + 1
           where=11.
           x_stop=xs
           y_stop=ys
          goto 500
          endif

! Go to Q2 OUT mag boundary.

          call transp(hms,6,decay_flag,dflag,m2,p,71.833d0)
          if ((xs*xs + ys*ys).gt.r_Q2*r_Q2) then
          Q2_out = Q2_out + 1
           where=12.
           x_stop=xs
           y_stop=ys
           goto 500
           endif

! Go to Q3 IN  mag bound.  Drift rather than using COSY matrices
!!  call transp(hms,7,decay_flag,dflag,m2,p)

          call project(xs,ys,94.225d0,decay_flag,dflag,m2,p) !project and decay
          if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
           Q3_in = Q3_in + 1
           where=13.
           x_stop=xs
           y_stop=ys
          goto 500
          endif

! Check aperture at 2/3 of Q3.

          call transp(hms,8,decay_flag,dflag,m2,p,145.7d0)
          if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
          Q3_mid = Q3_mid + 1
           where=14.
           x_stop=xs
           y_stop=ys
          goto 500
          endif

! Go to Q3 OUT mag boundary.

          call transp(hms,9,decay_flag,dflag,m2,p,72.9d0)
          if ((xs*xs + ys*ys).gt.r_Q3*r_Q3) then
           Q3_out = Q3_out + 1
           where=15.
           x_stop=xs
           y_stop=ys
           goto 500
          endif

! Go to D1 IN magnetic boundary, Find intersection with rotated aperture plane.
! Aperture has elliptical form.
c   call transp(hms,10,decay_flag,dflag,m2,p)
          call project(xs,ys,102.15d0,decay_flag,dflag,m2,p) !project and decay
          xt = xs      !  These were never filled before
          yt = ys      !  Why????????????
          call rotate_haxis(-6.0d0,xt,yt)
          checkdip = check_dipole(xt,yt)
         if (checkdip) then
          D1_in = D1_in + 1
           where=16.
           x_stop=xs    ! shouldn't these be xs and ys?
           y_stop=ys
            goto 500
         endif

! Go to D1 OUT magnetic boundary.
! Find intersection with rotated aperture plane.

            call transp(hms,11,decay_flag,dflag,m2,p,526.1d0)
          xt = xs
          yt = ys
c   call rotate_haxis(6.0d0,xt,yt)  ! why rotate back and not use xs,ys?
          checkdip = check_dipole(xt,yt)
          if (checkdip) then
           D1_out = D1_out + 1
           where=17.
           x_stop=xt
           y_stop=yt 
          goto 500
          endif

! Check a number of apertures in the vacuum pipes following the
! dipole.  First the odd piece interfacing with the dipole itself

          if ( (((xt-x_offset_pipes)**2+(yt-y_offset_pipes)**2)
     &         .gt.30.48**2).or. (abs((yt-y_offset_pipes))
     &         .gt.20.5232) ) then
           D1_out = D1_out + 1
           where=18.
           x_stop=xt
           y_stop=yt 
           goto 500
          endif

*
* now it gets tricky. first off, save xs and ys for safekeeping
*
*           xs_save=xs
*           ys_save=ys
*           xt_save=xt ! not needed anymore but just in case
*           yt_save=yt ! same as above

*
* put now xt in xs and yt in ys. The reason for this is that the 
* project routine we are about to call starts with xs/ys...
*


! Check the exit of the 26.65 inch pipe

          call project(xs,ys,64.77d0,decay_flag,dflag,m2,p) !project and decay
           if (((xs-x_offset_pipes-1.)**2+(ys-y_offset_pipes)**2).gt.1145.518)then
           D1_out = D1_out + 1
           where=19.
           x_stop=xt
          y_stop=yt
           goto 500
          endif

! check exit of long (117 inch) pipe (entrance is bigger than previous pipe)
! note: Rolf claims its 117.5 but the dravings say more like 116.x
! .. so i put 117 even.  Should be a 30.62 diameter pipe

!  Changed to 30.25 inch diameter - EC Nov,2001  !
!  Dropped the vacuum pipe by an extra 2 cm at the hut entrance  ! 
!  **  New Survey shows this pipe to be low by ~2.8 cm +/- .3    !
!  **  Nov. 2001                                                !

cc    call project(xs,ys,297.18d0,decay_flag,dflag,m2,p) !project and decay
cc!  Changed to 117.5 inch length and 30.25 diameter, EC June,2001  !
cc          if (((xs-x_offset_hut)**2
cc     &        +(ys-y_offset_hut)**2).gt.1475.90) then 
cc      D1_out = D1_out + 1
cc           where=20.
cc           x_stop=xt
cc           y_stop=yt
cc           goto 500
cc  endif

!!        Now project 60 cm further into the last pipe        !!
!!        This is optically the smallest section              !!  

          call project(xs,ys,298.58d0,decay_flag,dflag,m2,p) !project and decay
          if (((xs-x_offset_hut)**2
     &        +(ys-y_offset_hut)**2).gt.1475.90) then 
          D1_out = D1_out + 1
           where=20.
           x_stop=xt
           y_stop=yt
         goto 500
        endif

! lastly check the exit of the last piece of pipe. 45.5 inches long, 36.25 inch dia.
! Changed to 36.25 cm. diameter - EC June,2001 !

cc call project(xs,ys,+115.57d0,decay_flag,dflag,m2,p) !project and decay
cc          if (((xs-x_offset_hut)**2+
cc     &        (ys-y_offset_hut)**2).gt.2119.45)then
cc   D1_out = D1_out + 1
cc           where=21.
cc           x_stop=xt
cc           y_stop=yt
cc   goto 500
cc  endif

!!        Make up the 60cm from the last projection      !!

          call project(xs,ys,+114.17d0,decay_flag,dflag,m2,p) !project and decay
          if (((xs-x_offset_hut)**2+
     &        (ys-y_offset_hut)**2).gt.2119.45)then
           D1_out = D1_out + 1
           where=21.
           x_stop=xt
           y_stop=yt
          goto 500
         endif
*
* Now, if we passed all these tests restore the saved variables
* and go to the focal plane
*
*           xs=xs_save
*           ys=ys_save
*           xt=xt_save
*           yt=yt_save
*
*

! Note that we do NOT transport (project) to focal plane.  We will do this
! in mc_hms_hut.f so that it can take care of all of the decay, mult. scatt,
! and apertures.  Pass the current z position so that mc_hms_hut knows
! where to start.  Initial zposition for mc_hms_hut is -147.48 cm so that the
! sum of four drift lengths between pipe and focal plane is 625.0 cm
! (64.77+297.18+115.57+147.48=625)

! If we get this far, the particle is in the hut.

         shut = shut + 1
c	write(*,*)'shut',shut
! and track through the detector hut
c	  write(*,*) ' call mc_hms_hut'

        call mc_hms_hut(m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >  decay_flag,dflag,resmult,ok,-147.48d0)

c          write(*,*) ' return mc_hms_hut',ok
c          if(firsttime) write(6,*) ok,firsttime
  
           if (.not.ok) goto 500   


! replace xs,ys,... with 'tracked' quantities.
         xs=x_fp
          ys=y_fp
          dxdzs=dx_fp
          dydzs=dy_fp

! Reconstruct the target quantities. (Including the target field)
c       write(*,*) 'targ.mag.field mc_hms',fieldon
         call simc_hms_recon(dpp_recon,dth_recon,dph_recon,y_recon,
     >       -fry,xs,dxdzs,ys,dydzs)
c         write(*,*) ' call track to tgt'
         if ( fieldon) then
         p = p_spec*(1.+dpp_recon/100.)
c	     write(*,*) 'before go through the field',ok_hms	  
        call track_to_tgt(dpp_recon,y_recon,dph_recon,dth_recon,
     >       frx,fry,p*1000.,m2,cos_ts,sin_ts,spect,ok_hms,x_fp,
     >       dx_fp,y_fp,dy_fp,xtgt,bdl)
c	     write(*,*) 'after go through the field',ok_hms
         endif

! Fill output to return to main code
          dpp = dpp_recon
         dxdz = dph_recon
         dydz = dth_recon
         y = y_recon

         ok_hms = .true.
          where = 0.
          successes = successes + 1  
c       write(*,*) 'successes',successes

! We are done with this event, whether GOOD or BAD.

500      continue

        return
        end















