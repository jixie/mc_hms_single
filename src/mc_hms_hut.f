       subroutine mc_hms_hut (m2,p,x_fp,dx_fp,y_fp,dy_fp,ms_flag,wcs_flag,
     >           decay_flag,dflag,resmult,ok_hut,zinit)

C----------------------------------------------------------------------
C
C Monte-Carlo of HMS detector hut.
C Note that we only pass on real*8 variables to the subroutine.
C This will make life easier for portability.
C
C The particle is stepped through the detector (using project), and
C multiple scattering is applied for each detector or air gap.
C If particle decay in enabled, then project.f also checks for
C decay of particle.  The particle starts at z=zinit.  This
C needs to be before the first mult. scattering point (the exit window)
C or the decay distance is negative, and things are BAD.
C
C----------------------------------------------------------------------

        implicit none

        include 'hms_mc.cmn'
        include 'track.inc'

*
* G.&I.N. added
*
        include 'g_dump_all_events.inc'
***


C Math constants

        real*8 pi,d_r,r_d,root
        parameter (pi = 3.141592654)
        parameter (d_r = pi/180.0d0)
        parameter (r_d = 180.0d0/pi)
        parameter (root = 0.707106781) !square root of 1/2

        real*8 grnd,gauss1 !external functions

C----------------------------------------------------------------------
C HMS_MATERIALS
C CTP parameter file containing the materials of all the HMS detectors.
C For all materials AFTER the bend only have to do multiple scattering.
C     radlen = 1 radiation length (in cm)
C     thick  = thickness in cm
C In case a "+" is added to the comment, the thickness is a guess only.
C----------------------------------------------------------------------

C OLD:  spectrometer exit window, 15 mil Kevlar, 5 mil Mylar (use 20 mil,X0=53.3)
C NEW(2003-2005):                 20 mil Titanium, X0 = 3.56  -MEC

        real*8 hfoil_exit_radlen,hfoil_exit_thick
c       parameter (hfoil_exit_radlen = 53.3)
c        parameter (hfoil_exit_thick  = 0.020*2.54)
c        parameter (hfoil_exit_radlen = 3.56)
c 	parameter (hfoil_exit_thick  = 0.020*2.54)
        parameter (hfoil_exit_radlen = 8.9)
        parameter (hfoil_exit_thick  = 0.011*2.54)

C spectrometer air gaps
        real*8 hair_radlen
        parameter (hair_radlen = 30420.0d0)

C chamber entrance foil, 1 mil of Mylar
        real*8 hdc_entr_radlen,hdc_entr_thick
        parameter (hdc_entr_radlen = 28.7)
        parameter (hdc_entr_thick  = 0.001*2.54)

C chamber gas, 50/50 ethane/argon
        real*8 hdc_radlen,hdc_thick
        parameter (hdc_radlen = 16700.0)
        parameter (hdc_thick  = 1.8)

C effective wire plane, 25 micron W+Be/Cu gives <t>=pi*(.0025/2)**2
        real*8 hdc_wire_radlen,hdc_wire_thick
        parameter (hdc_wire_radlen = 0.35) !Assuming all Tungsten
        parameter (hdc_wire_thick  = 0.0000049)

C effective cathode plane, Be/Cu
        real*8 hdc_cath_radlen,hdc_cath_thick
        parameter (hdc_cath_radlen = 7.2) !'Ave' of Be/Cu
        parameter (hdc_cath_thick  = 0.000177)

C chamber exit foil, 1 mil of Mylar
        real*8 hdc_exit_radlen,hdc_exit_thick
        parameter (hdc_exit_radlen = 28.7)
        parameter (hdc_exit_thick  = 0.001*2.54)

C hodoscopes
        real*8 hscin_radlen
        parameter (hscin_radlen =  42.4)

C Cherenkov entrance foil, 40 mil of Al
        real*8 hcer_entr_radlen,hcer_entr_thick
        parameter (hcer_entr_radlen = 8.90)
        parameter (hcer_entr_thick  = 0.0394*2.54)

C Cerenkov gas.
        real*8 hcer_radlen
C 0.5 atm of CO2 for e/pi
C parameter (hcer_radlen = 36620.0d0)
C 0.5 atm of Freon for better e/pi
c  parameter (hcer_radlen = 9620.0d0)
C 2 atm of Freon for pi/p
C  parameter (hcer_radlen = 2405.0d0)
C 0.5 atm of C4F10
        parameter (hcer_radlen = 1105.0d0)    !!!  MEC - 4/28/05  !!!


C Cherenkov mirror, 75 mu plus 2 cm of Rotacell +
        real*8 hcer_mir_radlen,hcer_mir_thick
        parameter (hcer_mir_radlen = 400.0d0)
        parameter (hcer_mir_thick  = 2.0d0)

C Cherenkov exit foil
        real*8 hcer_exit_radlen,hcer_exit_thick
        parameter (hcer_exit_radlen = 8.90d0)
        parameter (hcer_exit_thick  = 0.0394*2.54)

C shower counter
        real*8 hcal_radlen
        parameter (hcal_radlen = 2.50)

C Wire chamber resolutions (sigma)
c	real*8 hdc_sigma(1:12)/	0.130,0.130,0.130,0.130,0.130,0.130,         !done
c     >				0.130,0.130,0.130,0.130,0.130,0.130/
        real*8 hdc_sigma(1:12)/ 0.090,0.090,0.090,0.090,0.090,0.090,
     > 0.090,0.090,0.090,0.090,0.090,0.090/
c	real*8 hdc_sigma(1:12)/	0.030,0.030,0.030,0.030,0.030,0.030,
c     >				0.030,0.030,0.030,0.030,0.030,0.030/



C Wire plane positions, construct hdc_zpos array using these parameters
        integer*4 hdc_nr_cham,hdc_nr_plan
        parameter (hdc_nr_cham = 2)                                          !done
        parameter (hdc_nr_plan = 6)

        real*8 hdc_1_zpos,hdc_1_left,hdc_1_right,hdc_1_top,hdc_1_bot
        real*8 hdc_1x_offset,hdc_1y_offset
        real*8 hdc_2_zpos,hdc_2_left,hdc_2_right,hdc_2_top,hdc_2_bot
        real*8 hdc_2x_offset,hdc_2y_offset
         real*8 hdc_del_plane

        parameter (hdc_1_zpos = -51.920)
        parameter (hdc_2_zpos =  29.291)
        parameter (hdc_del_plane = hdc_thick + hdc_wire_thick + hdc_cath_thick)  !above done
     
        parameter (hdc_1_left  =  26.0d0)       ! E99118
        parameter (hdc_1_right = -26.0d0)       ! E99118
c	 parameter (hdc_1y_offset = 1.443)     ! changed to match hdc.pos.nov95
        parameter (hdc_1y_offset = 1.2502)    ! E99118        
c	parameter (hdc_1y_offset = 0.00) ! test for ROSEN07
        parameter (hdc_1_top   = -56.5)
c        parameter (hdc_1_top   = -56.0d0)       ! E99118
        parameter (hdc_1_bot   =  56.5)
c        parameter (hdc_1_bot   =  56.0d0)       ! E99118
c	 parameter (hdc_1x_offset = 1.670)
        parameter (hdc_1x_offset = 1.6345)    ! E99118
c        parameter (hdc_1x_offset =0.00)    ! ROSEN07

        parameter (hdc_2_left  =  26.0d0)       ! E99118
        parameter (hdc_2_right = -26.0d0)       ! E99118 
c	 parameter (hdc_2y_offset = 2.753)     !changed as above,on 03/03,DD
        parameter (hdc_2y_offset = 2.651)     ! E99118
c        parameter (hdc_2y_offset = 0.00)     ! ROSEN07
        parameter (hdc_2_top   = -56.5)
c       parameter (hdc_2_top   = -56.0d0)       ! E99118
c	 parameter (hdc_2_bot   =  56.0d0)       ! E99118
        parameter (hdc_2_bot   =  56.5)
c	 parameter (hdc_2x_offset = 2.758) 
        parameter (hdc_2x_offset = 2.7825)    ! E99118
c        parameter (hdc_2x_offset = 0.00)    ! ROSEN07


C Scintillator positions and thicknesses                                         !done
        real*8 hscin_1x_zpos,hscin_1y_zpos
        real*8 hscin_1x_thick,hscin_1y_thick
        real*8 hscin_1x_left,hscin_1x_right,hscin_1x_offset
        real*8 hscin_1y_top,hscin_1y_bot,hscin_1y_offset
        real*8 hscin_2x_zpos,hscin_2y_zpos
        real*8 hscin_2x_thick,hscin_2y_thick
        real*8 hscin_2x_left,hscin_2x_right,hscin_2x_offset
        real*8 hscin_2y_top,hscin_2y_bot,hscin_2y_offset

        parameter (hscin_1x_zpos =  77.830)
        parameter (hscin_1y_zpos =  97.520)
        parameter (hscin_2x_zpos = 298.820)
        parameter (hscin_2y_zpos = 318.510)
        parameter (hscin_1x_thick = 1.067) !1 cm thick, 6.7% overlap
        parameter (hscin_1y_thick = 1.067)
        parameter (hscin_2x_thick = 1.067)
        parameter (hscin_2y_thick = 1.067)
        parameter (hscin_1x_left  =  37.75)
        parameter (hscin_1x_right = -37.75)
c parameter (hscin_1x_offset = -1.55)
        parameter (hscin_1x_offset = -1.30)    ! E99118
        parameter (hscin_1y_top   = -60.25)
        parameter (hscin_1y_bot   =  60.25)
c parameter (hscin_1y_offset = -0.37)
        parameter (hscin_1y_offset = -1.30)   ! E99118
        parameter (hscin_2x_left  =  37.75)
        parameter (hscin_2x_right = -37.75)
c parameter (hscin_2x_offset = -0.63)
        parameter (hscin_2x_offset = -0.60)   ! E99118
        parameter (hscin_2y_top   = -60.25)
        parameter (hscin_2y_bot   =  60.25)
c parameter (hscin_2y_offset = -1.46)
        parameter (hscin_2y_offset = -2.40)   ! E99118

C Cherenkov position
        real*8 hcer_zentrance,hcer_zmirror,hcer_zexit
        parameter (hcer_zentrance = 110.000)
        parameter (hcer_zmirror   = 230.000)
        parameter (hcer_zexit     = 265.000)

C Calorimeter position                                                                    !done
        real*8 hcal_1pr_zpos,hcal_2ta_zpos,hcal_3ta_zpos,hcal_4ta_zpos
        real*8 hcal_left,hcal_right,hcal_top,hcal_bottom
        parameter (hcal_1pr_zpos = 338.69)
        parameter (hcal_2ta_zpos = 349.69)
        parameter (hcal_3ta_zpos = 360.69)
        parameter (hcal_4ta_zpos = 371.69)
        parameter (hcal_left     =  35.00) !actual size of calorimeter
        parameter (hcal_right    = -35.00)
        parameter (hcal_top      = -70.40)     !was -69.66
        parameter (hcal_bottom   =  49.60)      !was 60.34
c parameter (hcal_top      = -69.66)
c parameter (hcal_bottom   =  60.34)

C The arguments

        real*8 p,m2 !momentum and mass of particle
        real*8 x_fp,y_fp,dx_fp,dy_fp !Focal plane values to return
        real*8 xcal,ycal !Position of track at calorimeter.
        real*8 zinit !Initial z-position (Not at F.P.)
        logical ms_flag !mult. scattering flag.
        logical wcs_flag !wire chamber smearing flag
        logical decay_flag !check for decay
        logical ok_hut !true if particle makes it

C Local declarations.

        integer*4 i,iplane,jchamber,npl_off
        integer*4 scintrig, scincount
        parameter (scintrig = 3) !set trigger to 3/4

        logical dflag !has particle decayed?

        real*8 resmult,tmpran,tmplim
        real*8 tmpran1,tmpran2 !temporary random numbers
        real*8 radw,drift

        real*8 nsig_max
        parameter(nsig_max=99.0d0) !max #/sigma for gaussian ran #s.

        real*8 dydzstemp,ztemp 

C These have to be real*4 for the CERNLIB lfit routine.
        real*4 badf !temporaries
        real*4 xfp4,yfp4,dxfp4,dyfp4 !real*4 versions of fp track.
        real*4 xdc(12),ydc(12),zdc(12) !positions at d.c. planes

C ================================ Executable Code =============================

C Initialize some variables
C These come from Doug's examination of events with >6 hits per chamber.
C In some fraction of events (~10% ?) there were >6 hits per chamber.  
C The position resolution for these event had approx. 3 times worse resolution.
C There was also some position (or delta) dependence.  The numbers below
C are based on examining one set of runs.  If you're going to use these
C to try and reproduce tails in the resolution, you should probably 
C see if the fraction with >6 hits, the resolution, and the delta depdendance
C are consistant with these numbers.  
C Note that if you play with the nominal DC resolutions, then you
C may want to reduce the resolution multiplier here.


! These are more conservative values (some mix of Kaon and NucPi values).
! Best to check what you see in your experiment.
        tmpran = grnd()
        tmplim = 0.10 !fraction of events with >6 hits.
        if (tmpran.lt.tmplim) then
          resmult = 1.25
        else
          resmult = 1.0
        endif

C Initialize ok_hut

        ok_hut = .false.

C Initialize the xdc and ydc arrays to zero

        do i=1,12
          xdc(i) = 0.0d0
          ydc(i) = 0.0d0
        enddo

C Initialize scincount to zero

        scincount = 0

C------------------------------------------------------------------------------C
C                           Top of loop through hut                            C
C------------------------------------------------------------------------------C

C Go to spectrometer exit foil, 25cm before DC1 (Drift forwards from zinit).
C As usual, neglect effect of nonzero dydzs and dxdzs on radw.


         drift = (hdc_1_zpos-25.0d0)-zinit !Need to adjust number if known

c        ztemp = zinit+drift 

c      write(6,*) hdc_1_zpos, zinit, drift   

         call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay
      
        ztemp = zinit+drift 

        radw = hfoil_exit_thick/hfoil_exit_radlen


        dydzstemp =  dydzs

        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
c      write(6,*) dydzs-dydzstemp

c      if(ms_flag) call musc_ld(m2,p,13.d0,27.d0,0.00314,dydzs,dxdzs)       
c       if(ms_flag) call musc_ld(m2,p,22.d0,47.9d0,0.115d0,dydzs,dxdzs)





C Go to first drift chamber set
C For simplicity, perform air mult. scattering (probably negligeable)
C before drifting instead of 1/2 way through.

         drift = (hdc_1_zpos - 0.5*hdc_nr_plan*hdc_del_plane) -
     > (hdc_1_zpos-25.0d0)
        radw = drift/hair_radlen
        if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
        call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay

        ztemp = ztemp + drift

        jchamber = 1
        radw = hdc_entr_thick/hdc_entr_radlen
        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
        npl_off = (jchamber-1)*hdc_nr_plan
        do iplane = 1,hdc_nr_plan
          radw = hdc_cath_thick/hdc_cath_radlen
          if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
          drift = 0.5*hdc_thick
         radw = drift/hdc_radlen
          if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
          drift = drift + hdc_cath_thick
          call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay

          ztemp = ztemp + drift

         radw = hdc_wire_thick/hdc_wire_radlen
          if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
          if(wcs_flag) then
           tmpran1 = gauss1(nsig_max) !gaussian, truncated at 99 sigma.
            tmpran2 = gauss1(nsig_max)
          else
           tmpran1 = 0.0d0
            tmpran2 = 0.0d0
          endif
          xdc(npl_off+iplane)=xs+hdc_sigma(npl_off+iplane)*tmpran1*resmult!+hdc_1x_offset !with DC offsets 
          ydc(npl_off+iplane)=ys+hdc_sigma(npl_off+iplane)*tmpran2*resmult!+hdc_1y_offset
          if (iplane.eq.2 .or. iplane.eq.5) then
           xdc(npl_off+iplane) = 0.0d0 !y plane, no x information   
          else
           ydc(npl_off+iplane) = 0.0d0 !x-like plane, no y info
          endif
          drift = 0.5*hdc_thick
          radw = drift/hdc_radlen
         if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
         drift = drift + hdc_wire_thick
          call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay

          ztemp = ztemp+drift

        enddo
        radw = hdc_exit_thick/hdc_exit_radlen
        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
        if ( xs.gt.(hdc_1_bot-hdc_1x_offset) .or.
     >       xs.lt.(hdc_1_top-hdc_1x_offset) .or.
     >       ys.gt.(hdc_1_left-hdc_1y_offset) .or.
     >       ys.lt.(hdc_1_right-hdc_1y_offset) ) then
         detec = detec + 1
          where=22.
          x_stop=xs
          y_stop=ys
         goto 500
         endif
        radw = hdc_cath_thick/hdc_cath_radlen
        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C at last cathode foil of first drift chamber set, drift to next

         drift = (hdc_2_zpos - 0.5*hdc_nr_plan*hdc_del_plane) -
     > (hdc_1_zpos + 0.5*hdc_nr_plan*hdc_del_plane)
        call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay

        ztemp = ztemp+drift

c       write(6,*) ztemp
 
        radw = drift/hair_radlen
        if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)

        jchamber = 2
        radw = hdc_entr_thick/hdc_entr_radlen
        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
        npl_off = (jchamber-1)*hdc_nr_plan
        do iplane = 1,hdc_nr_plan
          radw = hdc_cath_thick/hdc_cath_radlen
         if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
          drift = 0.5*hdc_thick
          radw = drift/hdc_radlen
          if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
          drift = drift + hdc_cath_thick
           call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay
          radw = hdc_wire_thick/hdc_wire_radlen
          if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
          if(wcs_flag) then
           tmpran1 = gauss1(nsig_max) !gaussian, truncated at 99 sigma.
           tmpran2 = gauss1(nsig_max)
          else
           tmpran1 = 0.0d0
           tmpran2 = 0.0d0
          endif
         xdc(npl_off+iplane)=xs+hdc_sigma(npl_off+iplane)*tmpran1*resmult!+hdc_2x_offset
         ydc(npl_off+iplane)=ys+hdc_sigma(npl_off+iplane)*tmpran2*resmult!+hdc_2y_offset
         if (iplane.eq.2 .or. iplane.eq.5) then
          xdc(npl_off+iplane) = 0.0d0   !y plane, no x information
         else
           ydc(npl_off+iplane) = 0.0d0   !x-like plane, no y info
         endif

          drift = 0.5*hdc_thick
          radw = drift/hdc_radlen
         if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
          drift = drift + hdc_wire_thick
         call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay

          ztemp = ztemp + drift

        enddo
        radw = hdc_exit_thick/hdc_exit_radlen
        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
        if ( xs.gt.(hdc_2_bot-hdc_2x_offset) .or.
     >       xs.lt.(hdc_2_top-hdc_2x_offset) .or.
     >       ys.gt.(hdc_2_left-hdc_2y_offset) .or. 
     >       ys.lt.(hdc_2_right-hdc_2y_offset) ) then
         detec = detec + 1
           where=23.
           x_stop=xs
           y_stop=ys
          goto 500
         endif
        radw = hdc_cath_thick/hdc_cath_radlen
        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C fit track to give new focal plane values, use LFIT from GENLIB
        do jchamber=1,hdc_nr_cham
          npl_off = (jchamber-1)*hdc_nr_plan
          do iplane=1,hdc_nr_plan
            if (jchamber.eq.1) zdc(npl_off+iplane) = hdc_1_zpos +
     > (iplane-0.5-0.5*hdc_nr_plan)*hdc_del_plane
            if (jchamber.eq.2) zdc(npl_off+iplane) = hdc_2_zpos +
     > (iplane-0.5-0.5*hdc_nr_plan)*hdc_del_plane
          enddo
         enddo




        call lfit(zdc,xdc,12,0,dxfp4,xfp4,badf)
        call lfit(zdc,ydc,12,0,dyfp4,yfp4,badf)

        x_fp = dble(xfp4)
        y_fp = dble(yfp4)
        dx_fp = dble(dxfp4)
        dy_fp = dble(dyfp4)

C at last cathode foil of second drift chamber set, drift to hodoscopes

        drift = hscin_1x_zpos - hdc_2_zpos - 0.5*hdc_nr_plan*hdc_del_plane
        radw = drift/hair_radlen
        if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
        call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay
        if (ys.lt.(hscin_1x_left+hscin_1y_offset) .and. 
     >      ys.gt.(hscin_1x_right+hscin_1y_offset) .and.
     >      xs.lt.(hscin_1y_bot+hscin_1x_offset) .and.
     >      xs.gt.(hscin_1y_top+hscin_1x_offset) ) scincount = 
     &             scincount + 1
        radw = hscin_1x_thick/hscin_radlen
        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

        drift = hscin_1y_zpos - hscin_1x_zpos
        radw = drift/hair_radlen
        if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
        call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay
        if (ys.lt.(hscin_1x_left+hscin_1y_offset) .and. 
     >      ys.gt.(hscin_1x_right+hscin_1y_offset) .and.
     >      xs.lt.(hscin_1y_bot+hscin_1x_offset) .and.
     >      xs.gt.(hscin_1y_top+hscin_1x_offset) ) scincount = 
     &             scincount + 1
        radw = hscin_1y_thick/hscin_radlen
        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C finished first hodoscope, drift to cherenkov

        drift = hcer_zentrance - hscin_1y_zpos
        radw = drift/hair_radlen
        if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
        call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay

        radw = hcer_entr_thick/hcer_entr_radlen
        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

        drift = hcer_zmirror - hcer_zentrance
        radw = drift/hcer_radlen
        if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
        call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay

        radw = hcer_mir_thick/hcer_mir_radlen
        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

        drift = hcer_zexit - hcer_zmirror
        radw = hcer_exit_thick/hcer_exit_radlen
        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)
        call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay

        radw = hcer_exit_thick/hcer_exit_radlen
        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C drift to second hodoscope

        drift = hscin_2x_zpos - hcer_zexit
        radw = drift/hair_radlen
        if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
        call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay
        if (ys.lt.(hscin_1x_left+hscin_1y_offset) .and. 
     >      ys.gt.(hscin_1x_right+hscin_1y_offset) .and.
     >      xs.lt.(hscin_1y_bot+hscin_1x_offset) .and.
     >      xs.gt.(hscin_1y_top+hscin_1x_offset) ) scincount = 
     &            scincount + 1
        radw = hscin_2x_thick/hscin_radlen
        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

        drift = hscin_2y_zpos - hscin_2x_zpos
        radw = drift/hair_radlen
        if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
        call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay
        if (ys.lt.(hscin_1x_left+hscin_1y_offset) .and. 
     >      ys.gt.(hscin_1x_right+hscin_1y_offset) .and.
     >      xs.lt.(hscin_1y_bot+hscin_1x_offset) .and.
     >      xs.gt.(hscin_1y_top+hscin_1x_offset) ) scincount = 
     &         scincount + 1
        radw = hscin_2y_thick/hscin_radlen
        if(ms_flag) call musc(m2,p,radw,dydzs,dxdzs)

C Test on scintillator trigger (DJG 8/18/98)

        if( scincount .lt. scintrig ) then
          detec = detec + 1
          where=24.
          x_stop=xs
          y_stop=ys
         goto 500
        endif

C Drift to calorimeter

        drift = hcal_4ta_zpos - hscin_2y_zpos
        radw = drift/hair_radlen
        if(ms_flag) call musc_ext(m2,p,radw,drift,dydzs,dxdzs,ys,xs)
        call project(xs,ys,drift,decay_flag,dflag,m2,p) !drift and decay

C Note that even with the standard PID trigger, the calorimeter is NOT
C required, since the trigger needs either the cerenkov OR the calorimeter.
C If you require the calorimeter in you analysis, then you need to make
C sure that the calorimeter is hit here, AND apply a seperate post-tracking
C fiducial cut (whatever is required by the engine or in your analysis).
C The fiducial cut comes later in the code.
C
C This cut simply requires that the event hit the calorimeter.  It does
C not insure that the event is within the good region of the calorimeter.
C The seperate fiducial cut is required to insure that the entire shower
C energy is contained in the calorimeter.  Or, if you like, you can require
C some distance between the track and the edge.  2-3 cm seems to be enough
C to get most or all of the energy in the calorimeter

!c           xif (ys.gt.(hcal_left-2.0d0) .or. ys.lt.(hcal_right+2.0d0) .or.
!c     >     xs.gt.(hcal_bottom-2.0d0) .or. xs.lt.(hcal_top+2.0d0)) then
!           if (ys.gt.hcal_left .or. ys.lt.hcal_right .or.
!     >      xs.gt.hcal_bottom .or. xs.lt.hcal_top) then
!            detec = detec + 1
!           goto 500
!           endif

C If you use a calorimeter cut in your analysis, the engine applied a
C a fiducial cut at the calorimeter.  This is based on the position of the
C TRACK at the calorimeter, not the real position of the event.  Go to the
C back of the calorimeter since engine uses a FID cut at the back.
C The standard fiducial cut is 5 cm from the edges of the block.

        xcal = x_fp + dx_fp * hcal_4ta_zpos
        ycal = y_fp + dy_fp * hcal_4ta_zpos
!        if (ycal.gt.(hcal_left-5.0d0) .or. ycal.lt.(hcal_right+5.0d0) .or.
!     >     xcal.gt.(hcal_bottom-5.0d0) .or. xcal.lt.(hcal_top+5.0d0)) then

        ok_hut = .true.
        where = 0.d0

C We are done with this event, whether GOOD or BAD.

500     continue

        return
        end
