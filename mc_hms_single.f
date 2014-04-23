        program mc_hms_single
C+______________________________________________________________________________
!
! Monte-Carlo of HMS spectrometer using uniform illumination.
!   This version uses TRANSPORT right-handed coordinate system.
!
! Author: David Potterveld, March-1993
!
! Modification History:
!
!  11-Aug-1993	(D. Potterveld) Modified to use new transformation scheme in
!which each transformation begins at the pivot.
!
!  19-AUG-1993  (D. Potterveld) Modified to use COSY INFINITY transformations.
!
!  06-MAR-2008  (Anusha) Modified the codes to use the target field corrections.
!               Used the subroutine trgInit.f to call the target field map and 
!               another couple of subroutines in the mc_hms.f.
C-______________________________________________________________________________

        implicit none

        include 'hms_mc.cmn'
        include 'hms_magnets.cmn'
        include 'track.inc'
        include 'constants.inc'
        include 'addmore.cmn'
*
* G.&I.N. added
*
        include 'g_dump_all_events.inc'
***

C Math constants
        real*8 root
        parameter (root = 0.707106781) !square root of 1/2
C HBOOK/NTUPLE common block and parameters.
        integer*4 pawc_size
        parameter (pawc_size = 4000000)
        common       /pawc/ hbdata(pawc_size)
        common       /quest/iquest
        integer        iquest(100)
         integer*4 hbdata
        character*8 hut_nt_names(19)/
     > 'XFOC', 'YFOC', 'DYDZ', 'DXDZ',
     > 'ZINI', 'DPPI', 'DTHI', 'DPHI',
     > 'ZREC', 'DPPR', 'DTHR', 'DPHR','FAIL_ID',
     > 'x_stop','y_stop','yini','yrec',
     > 'rast_x', 'rast_y'/
        real*4 hut(19)
C Local declarations.
        integer*4 i,j,
     > chan /1/,
     > n_trials,trial,
     > tmp_int
        integer*4       spect

        logical*4 iss
        integer*4       map
c	character*17 trg_field_map /'trg_field_map.dat'/
        real*8 theta_e,theta_p,phi_e,phi_p
C Histogram defining parameters.
        integer*4 ix,iy,iz,it           !Indices
        integer*4 h_nx,h_ny,h_nz,h_nt   !Number of dp/p,theta,phi,ytar bins
        real*8 h_xl,h_yl,h_zl,h_tl      !Hist. low edge in dp/p,theta,phi,y_tgt
        real*8 h_dx,h_dy,h_dz,h_dt      !Hist. bin width
        parameter (h_nx = 64)
        parameter (h_ny = 64)
        parameter (h_nz = 64)
        parameter (h_nt = 64)

        integer*4 cnts_ztgt(h_nt)
        integer*4 cnts_phi(h_nz)
        integer*4 cnts_th(h_ny)
        integer*4 cnts_p(h_nx)
        integer*4 cnts_pphi(h_nx,h_nz) !P-phi histogram
        integer*4 cnts_pth(h_nx,h_ny) !P-theta histogram
        integer*4 cnts_thphi(h_ny,h_nz) !theta-phi histogram
        integer*4 cnts_zth(h_nt,h_ny) !Z-theta histogram
        integer*2 cnts_pthph(h_nx,h_ny,h_nz) !P-theta-phi histogram

C Event limits, topdrawer limits, physics quantities
        real*8 gen_lim(6)     !M.C. phase space limits.
*
        real*8 gen_lim_up(3)
        real*8 gen_lim_down(3)
*
* G.&I.N. added up/down limits for dpp,dysz,dydz...
*
        real*8 dpp_hist,dth_hist,dph_hist,z_hist !Histogram limits for Topdrawer
        real*8 cut_dpp,cut_dth,cut_dph,cut_z !cuts on reconstructed quantities
        real*8 cos_ts,sin_ts      !cos and sin of spectrometer angle
        real*8 th_ev,cos_ev,sin_ev !cos and sin of event angle
        real*8 x,y,z,dpp,dxdz,dydz,t1,t2,t3,t4!temporaries
        real*8 xf,yf
        real*8 r2,th,r
        real*8 musc_targ_len !target length for multiple scattering
        real*8 tf_w                             !target window length for multiple scattering
        real*8 m2,p !particle mass squared.
        real*8 rad_len_cm !conversion r.l. to cm for target
        logical ok_hms !indicates whether event makes it in MC

C Initial and reconstructed track quantities.
        real*8 dpp_init,dth_init,dph_init,xtar_init
        real*8 ytar_init,ztar_init,xs_init,ys_init,zs_init
        real*8 dpp_recon,dth_recon,dph_recon,ztar_recon,ytar_recon
        real*8 x_fp,y_fp,dx_fp,dy_fp !at focal plane
        real*8 p_spec,th_spec          !spectrometer setting
        real*8 resmult

C CHANGE THESE FLAGS IN mc_hms.inp - NOT HERE!
        integer*4 p_flag !particle identification
        logical*4 ms_flag /.false./ !Flag for multiple scattering
        logical*4 wcs_flag /.false./ !Flag for wire-chamber smear.
        logical*4 fit_reverse /.false./ !Flag for fitting reverse coeff.
        logical*4 topd_flag /.false./
        logical*4 pth_flag /.false./
        logical*4 pphi_flag /.false./
        logical*4 thphi_flag /.false./
        logical*4 zth_flag /.false./
        logical*4 init_flag /.false./
        logical*4 hut_ntuple /.false./

        integer*4 ttl_indx
        real*8 dpp_var(2),dth_var(2),dph_var(2),ztg_var(2)
        real*8 t5,t6,rotang,xoff,yoff,zt,zoff
        real*8 hms_mispoint,t7,t8

        character*132 str_line
        character*9 TTL_name(0:1) /'RECON    ','INIT   '/

C Function definitions.
        
        real*8          tarlength
        integer*4 last_char
        logical*4 rd_real,rd_int
        logical         solid_tar
        logical         dummy_tar


        real*8 grnd,gauss1

        character*17 read_name
        character*70 inp_name
        character*70 ntp_name
        character*70 my_name
        character*10 inp_ext
        character*10 ntp_ext
        character*10 out_ext

        logical fieldon /.false./
        common /mag/fieldon
                                                                          
        save !Remember it all!

        inp_ext = '.inp'
        ntp_ext = '.rzdat'
        out_ext = '.dat'




C ================================ Executable Code =============================
       
        xoff = 0.0d0
        yoff = 0.0d0
        zoff = 0.0d0       

C Initialize

        trials = 0
        slit_hor = 0
        slit_vert = 0
        slit_oct = 0
        Q1_in = 0
        Q1_mid = 0
        Q1_out = 0
        Q2_in = 0
        Q2_mid = 0
        Q2_out = 0
        Q3_in = 0
        Q3_out = 0
        D1_in = 0
        D1_out = 0
        shut = 0
        incut = 0
        detec = 0
        successes= 0
        where = 0

        firsttime = .true.
c Adjust the beam_offset on the target X direction
C       beam_offset = 0.5
        beam_offset = 0.0
C Open setup file.


        write(6,*)'Give the name of the input file without extension (will be same for the output) '
        read(5,1968)read_name

 1968   format(a)

        inp_name = 'infiles/'//read_name(1:last_char(read_name))//inp_ext
        write(6,*)'input file = ',inp_name
        open (unit=chan,status='old',file=inp_name)
*       open (unit=10,status='old',name='c_kin.dat')
*	open (unit=chan,status='old',name='mc_hms.inp')

C Strip off header
! N_TRIALS:

C Strip off header
         read (chan,1001) str_line
        do while (str_line(1:1).eq.'!')
         read (chan,1001) str_line
        enddo
c
        write(*,*),str_line(1:last_char(str_line))
        iss = rd_int(str_line,n_trials)
        
        if (.not.iss) stop 'ERROR (ntrials) in setup!'

! Spectrometer momentum:
c	read (chan,*) p_spec
        read (chan,1001) str_line
        iss = rd_real(str_line,p_spec)
        if (.not.iss) stop 'ERROR (Spec momentum) in setup!'


! Spectrometer angle:
        read (chan,1001) str_line
        iss = rd_real(str_line,th_spec)
        if (.not.iss) stop 'ERROR (Spec theta) in setup!'
        th_spec = abs(th_spec) / degrad
        cos_ts = cos(th_spec)
        sin_ts = sin(th_spec)
 
! M.C. limits (half width's for dp,th,ph, full width's for x,y,z)
*
* G.&I.N.
*
         do i=1,3
          read (chan,1001) str_line
          write(6,*) str_line(1:last_char(str_line))
          iss = rd_real(str_line,gen_lim(i))
          if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
         gen_lim_down(i) = gen_lim(i)
          read (chan,1001) str_line
           write(6,*) str_line(1:last_char(str_line))
          iss = rd_real(str_line,gen_lim(i))
          if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
           gen_lim_up(i) = gen_lim(i)
         enddo
         do i=1,3
c	   write(6,*)'gen_lim_up/down = ',gen_lim_up(i),' '
c	1	,gen_lim_down(i)
          enddo
*
        do i = 4,6
         read (chan,1001) str_line
          write(6,*) str_line(1:last_char(str_line))
          iss = rd_real(str_line,gen_lim(i))
          if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
        enddo

! Cuts on reconstructed quantities
        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        if (.not.rd_real(str_line,cut_dpp)) stop 'ERROR (CUT_DPP) in setup!'

        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        if (.not.rd_real(str_line,cut_dth)) stop 'ERROR (CUT_DTH) in setup!'

        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        if (.not.rd_real(str_line,cut_dph)) stop 'ERROR (CUT_DPH) in setup!'

        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        if (.not.rd_real(str_line,cut_z)) stop 'ERROR (CUT_Z) in setup!'

! Limits defining topdrawer histograms.
        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        if (.not.rd_real(str_line,dpp_hist)) stop 'ERROR (DPP_HIST) in setup!'

        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        if (.not.rd_real(str_line,dth_hist)) stop 'ERROR (DTH_HIST) in setup!'

        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        if (.not.rd_real(str_line,dph_hist)) stop 'ERROR (DPH_HIST) in setup!'

        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        if (.not.rd_real(str_line,z_hist)) stop 'ERROR (Z_HIST) in setup!'

! Read in radiation length of target material in cm
        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        if (.not.rd_real(str_line,rad_len_cm)) stop 'ERROR (RAD_LEN_CM) in setup!'

! read in flag for particle type.

        read(chan,*) p_flag
        write(6,*) p_flag
cc	read (chan,1001) str_line
cc	write(6,*) str_line(1:last_char(str_line))
cc	if (.not.rd_int(str_line,p_flag)) stop 'ERROR: p_flag in setup file!'

! Read in flag for multiple scattering.

        read(chan,*) tmp_int
cc	read (chan,1001) str_line
cc	write(6,*) str_line(1:last_char(str_line))
cc if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: ms_flag in setup file!'
        if (tmp_int.eq.1) ms_flag = .true.
        write(6,*) ms_flag, "    =   Do multiple scattering:"

! Read in flag for wire chamber smearing.

        read(chan,*) tmp_int
cc	read (chan,1001) str_line
cc	write(6,*) str_line(1:last_char(str_line))
cc	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR: wcs_flag in setup file!'
        if (tmp_int.eq.1) wcs_flag = .true.
        write(6,*)  wcs_flag, "    =   Do wirechamber smearing?  "

! Read in flag FIT_REVERSE.
        
        read(chan,*) tmp_int
cc	read (chan,1001) str_line
cc	write(6,*) str_line(1:last_char(str_line))
cc	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR:fit_reverse in setup file!'
        if (tmp_int.eq.1) fit_reverse = .true.
        write(6,*)  fit_reverse,"    =   Fit reverse matrix elements?  "

! Read in flag for init/recon values.

        read(chan,*) tmp_int
cc	read (chan,1001) str_line
cc	write(6,*) str_line(1:last_char(str_line))
cc	if (.not.rd_int(str_line,tmp_int)) stop 'ERROR:init_flag in setup file!'
        if (tmp_int.eq.1) then
          init_flag = .true.
          ttl_indx = 1
        else
           ttl_indx = 0
        endif
        write(6,*) ttl_indx


! Read in flag for TOPDRAWER output.

        read(chan,*) tmp_int
        read (chan,1001) str_line
cc  write(6,*) str_line(1:last_char(str_line))
cc if (.not.rd_int(str_line,tmp_int)) stop 'ERROR:topd_flag in setup file!'
        if (tmp_int.eq.1) topd_flag = .true.
        write(6,*) topd_flag   

! Read in flag for HUT NTUPLE output.

        read(chan,*) tmp_int
cc read (chan,1001) str_line
cc write(6,*) str_line(1:last_char(str_line))
cc if (.not.rd_int(str_line,tmp_int)) stop 'ERROR:hut_ntuple in setup file!'
        if (tmp_int.eq.1) hut_ntuple = .true.
        write(6,*) hut_ntuple,"    =   output Hut ntuple"
* 
* N.G. added
*
! Read in flag for dumping ALL events into the HUT NTUPLE (note that
! the ...recon quantities will be ill defined but the FAIL_ID could be
! used to tell...


        tmp_int = 1
cc       read(chan,*) tmp_int
cc read (chan,1001) str_line
cc write(6,*) str_line(1:last_char(str_line))
cc if (.not.rd_int(str_line,tmp_int)) stop
cc 1    'ERROR:dump_all_in_ntuple in setup file!'
        if (tmp_int.eq.1) dump_all_in_ntuple = .true.
        write(6,*) dump_all_in_ntuple, "    =   Dumping all events?"

! Read in flag for doing sieve slit run.

        read(chan,*) tmp_int
cc read (chan,1001) str_line
cc write(6,*) str_line(1:last_char(str_line))
cc if (.not.rd_int(str_line,tmp_int)) stop 'ERROR:hut_ntuple in setup file!'
        if (tmp_int.eq.1) dosieve = .true.
        write(6,*) dosieve,"    =   Sieve Slit run? "
        

CCCCCCCC       Read in magnet corrections            CCCCCCCCC


        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        iss = rd_real(str_line,q1cor)
        if (.not.iss) stop 'ERROR (q1cor) in setup!'

        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        iss = rd_real(str_line,q2cor)
        if (.not.iss) stop 'ERROR (q2cor) in setup!'

        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        iss = rd_real(str_line,q3cor)
        if (.not.iss) stop 'ERROR (q3cor) in setup!'

        dcor = 0.

        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        iss = rd_real(str_line,dcor)
        if (.not.iss) stop 'ERROR (dcor) in setup!'

        if(dcor.LT.0..AND.p_spec.GT.3.5573) then      !  if input is < 0. then calculate dipole  !
         dcor = 100.*1.0755e-3*((p_spec-3.5573)**2)   !  phantom saturation effect               !      
         write(6,*)'Including phantom HMS saturation correction,
     &    dP(%)=',dcor
        else
         dcor = 0.d0
        endif

c        write(6,*) dcor, p_spec


CCCCCCCC       Beam position and target offsets         CCCCCCCCC


        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        iss = rd_real(str_line,xoff)
        if (.not.iss) stop 'ERROR (xoff) in setup!' 

        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        iss = rd_real(str_line,yoff)
        if (.not.iss) stop 'ERROR (yoff) in setup!'

        read (chan,1001) str_line
        write(6,*) str_line(1:last_char(str_line))
        iss = rd_real(str_line,zoff)
        if (.not.iss) stop 'ERROR (zoff) in setup!'

        read(chan,*) tmp_int
        if (tmp_int.eq.1) fieldon = .true.
        write(6,*) 'using the targ. mag. field',fieldon
       

**********
        theta_lab=th_spec*degrad
        sry_off=yoff
        srx_off=xoff !beam right

C Set particle masses.
        m2 = me2    !default to electron
       if(p_flag.eq.0) then
        m2 = me2
        spect=-1
        else if(p_flag.eq.1) then
          m2 = mp2
        spect=1
        else if(p_flag.eq.2) then
           m2 = md2
        else if(p_flag.eq.3) then
          m2 = mpi2
        else if(p_flag.eq.4) then
          m2 = mk2
        endif

C Set histogram limits:
        h_xl = -dpp_hist
        h_dx = 2*dpp_hist/h_nx

        h_yl = -dth_hist
        h_dy = 2*dth_hist/h_ny

        h_zl = -dph_hist
        h_dz = 2*dph_hist/h_nz

        h_tl = -z_hist
        h_dt = 2*z_hist/h_nt

C Initialize HBOOK/NTUPLE if used.
        if (hut_ntuple) then
         call hlimit(pawc_size)
         ntp_name = 'worksim/'//read_name(1:last_char(read_name))//ntp_ext
c	iquest(10)=65000
c        call hropen(30,'HUT',ntp_name,'NQ',4096,i)
          call hropen(30,'HUT',ntp_name,'N',1024,i)  !!  MEC  
        if (i.ne.0) then
           write(6,*),'HROPEN error: istat = ',i
           stop
         endif
         call hbookn(1,'HUT NTUPLE',19,'HUT',1000,hut_nt_names)
        endif
C Output results.
c write(*,*)'Give name of the output file'
c read(*,1968)my_name

         my_name = 'outfiles/'//read_name(1:last_char(read_name))//out_ext
C------------------------------------------------------------------------------C
C                           Top of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C

        if(dump_all_in_ntuple) then 
         write(6,*) "Dumping all events to ntuple"
        else
         write(6,*) "Only dumping good events to ntuple"
        endif

c        stime = 1.
c	stime = secnds(0.)


        first_time_hms = .true.
        firsttime = .true.

        if(spect.eq.-1) then ! electron beam left is negative
         theta_e = (th_spec*degrad)+80  ! perpendicular mag. field
         phi_e=180                      ! perpendicular mag. field
c    theta_e = 180 - (th_spec*degrad)! parallel mag. field
c    phi_e=0                         ! parallel mag. field
           theta_p = 0
         phi_p=0
        elseif(spect.eq.1) then ! proton
           theta_e = 0
          phi_e= 0
          theta_p = th_spec
          phi_p = 180
          endif
*        bdl=0

!! call the target magnetic field strength.******************************

        if ( fieldon) then
c write(*,*) 'theta_e,phi_e',theta_e,phi_e
        call trgInit(map,theta_e,phi_e,theta_p,phi_p) 
       endif
*************************************************************************
cc	p_spec_offset = -0.13 
c	p_spec = p_spec * ( 1. + p_spec_offset / 100. )

        do trial = 1,n_trials

c       if(first_time_hms) write(6,*) first_time_hms
c          if(firsttime) write(6,*) firsttime
        
* G.&I.N.
* reset the fail id
          where=0.d0
         x_stop=0.0
         y_stop=0.0
*
         if(mod(trial,10000).eq.0) write(*,*)'event #: ',trial


C Pick starting point within target. Z is picked uniformly, X and Y are
C chosen as truncated gaussians, using numbers picked above.
C Units are cm.


          x = gauss1(3.0d0)*0.003 ! intrinsic beam width
          y = gauss1(3.0d0)*0.003                        ! intrinsic beam height

c  x = gauss1(3.0d0)                        ! intrinsic beam width
c          y = gauss1(3.0d0)                        ! intrinsic beam height


CCCCC     Include  rastering - MEC     CCCCC

!!    Sinusoidal raster OFF    !!

          if(gen_lim(4).GE.0) then

           t5 = (2.*grnd() - 1.)*gen_lim(4)
            t6 = (2.*grnd() - 1.)*gen_lim(5)

             x = x + t5                  
             y = y + t6
            x = x + xoff !! X-beam (horizontal) +x is beam left  !! 
            y = y + yoff !! Y-beam (vertical)   +y is beam up     !!
 
          else

!!    Uniform sinusoidal raster ON   !!

            xf = (2.*grnd() - 1.)*0.05 
             yf = (2.*grnd() - 1.)*0.05
             t5 = (2.*grnd() - 1.)*gen_lim(4)
              t6 = (2.*grnd() - 1.)*gen_lim(5)
 
             r2 = grnd()*abs(gen_lim(4))
             th = 360.*grnd()

             r = sqrt(r2)

c     r = grnd()*gen_lim(4)
  
           t7=r*cos(th)
            t8=r*sin(th)
     
c      x = x + t5+t7                  
c      y = y + t6+t8 
c      x = x + xoff ! +beam_offset	!! X-beam (horizontal) +x is beam left  !! 
c      y = y + yoff	!! Y-beam (vertical)   +y is beam up     !!
 
c      x = x + t7                  
c      y = y + t8 
         x = t7                  
         y = t8 
         yrecs=t8 + yoff
          x = t7 + xoff+xf ! +beam_offset	!! X-beam (horizontal) +x is beam right  !! 
          y = t8 + yoff+yf                       !! +y is pointing up  !! 

          endif
c   open(2,file='raster.dat',status='unknown')
c   write (2,*) t7,xoff,x,t8,yoff,y



CCCCC        End rastering stuff        CCCCC

          solid_tar = .false.         !!  Initialize to false  !!
          dummy_tar = .false.
          
          if(gen_lim(6).LT.0.)then
           dummy_tar = .true.
           tarlength = abs(gen_lim(6))
          endif  
          if(.not.dummy_tar.AND.gen_lim(6).LE.3.0d0)  !!  If less than 3cm then assume solid target !!
     &             solid_tar = .true. 
CCCCC    Solid Target was rotated by 20.3 deg for E99-118   CCCCC


c          rotang = 20.3*pi/180.0d0

          rotang = 0.0d0  

          if(dosieve) then
           rotang = 0.0d0
          endif

          if(solid_tar) then
           zt = cos(rotang)*(grnd()-0.5)*gen_lim(6)           !!  along beam direction relative to target center !! 
           if(rotang.GT.0) zt = zt + (x*tan(rotang))          !!  target rotation                   !! 
           z = zt + zoff                                      !!  ztarget offset (+ is downstream)  !!
          elseif(dummy_tar) then                              !!  Currently assumes 4cm dummy       !! 
           z = (int(2.0d0*grnd())-0.5)*4.0d0
     &             +(grnd()-0.5)*tarlength/2.0d0
          else 
           zt = (1.0d0*grnd()-0.50)*gen_lim(6)          
           z = zt + zoff
          endif

c          write(6,*) (z/gen_lim(6))


***    The factor cos(rotang)  takes into                     ***
***    account the correct thickness due to rotation.  The    ***
***    z = z + (x*tan(rotang))  term takes into               ***
***    account that the central z position of the target      ***
***    walks with x beam.                                     ***


c       write(6,*) trial,x,y,z

C Pick scattering angles and DPP from independent, uniform distributions.
C dxdz and dydz in HMS TRANSPORT coordinates.


         dpp  = grnd()*(gen_lim_up(1)-gen_lim_down(1))
     &             + gen_lim_down(1)



          dydz = (grnd()*(gen_lim_up(2)-gen_lim_down(2))/1000.
     &           + gen_lim_down(2)/1000.)

          dxdz = grnd()*(gen_lim_up(3)-gen_lim_down(3))
     &          /1000.   + gen_lim_down(3)/1000.




C Transform from target to HMS (TRANSPORT) coordinates.
C Note that this assumes that HMS is on the right-hand side of the beam
C line (looking downstream).
c  The beam coordinate system is +y (vertical up) , +z ( downstream towards dump)
c    and +x should be beam left

          xs    = -y
         ys    = -(x * cos_ts) + z * sin_ts
           zs    = z * cos_ts - x * sin_ts

c          write(6,*) "1: ",ys
         ys = ys - hms_mispoint(th_spec*degrad)              !!  include hms mispointing  !!
c          write(6,*) "2: ",ys
         dpps  = dpp
         dxdzs = dxdz
          dydzs = dydz
c  write(*,*) 'theta corr',dxdzs,dydzs,dydzs - (0.025*dxdzs)
c  dydzs = dydzs + (0.025*dxdzs)

          xs_init = xs
          ys_init = ys
          zs_init = zs
c   write(*,*) ' xs = ',xs
c   write(*,*) ' ys = ',ys
c  write(*,*) ' zs = ',zs

C Save init values for later.
         xtar_init = x
         ytar_init = y
         ztar_init = z
         dpp_init = dpp
          dth_init = dydzs*1000. !mr
         dph_init = dxdzs*1000. !mr

C Calculate multiple scattering length of target
         cos_ev = (cos_ts+dydzs*sin_ts)/sqrt(1+dydzs**2+dxdzs**2)
         th_ev = acos(cos_ev)
          sin_ev = sin(th_ev)
  
C Case 1 : cryo target (2.65 inches diamater --> 3.37 cm radius)
C   Liquid + 5 mil Al (X0=8.89cm) beer can
          if (abs(gen_lim(6)).gt.3.) then
*   if (abs(gen_lim(6)/2. - z)/cos_ev .lt. 3.37/sin_ev) then !forward 
* N.G. added
             if (((gen_lim(6)/2. - z)/cos_ev .lt. 3.37/sin_ev).
     &           and.((gen_lim(6)/2. + z)/cos_ev .gt. 3.37/sin_ev)) then !forward or backward...
*	write(*,*)'forward/backward'
        	musc_targ_len = abs(gen_lim(6)/2. - z)/rad_len_cm/cos_ev

         	musc_targ_len = musc_targ_len + .005*2.54/8.89/cos_ev
           else !side
               musc_targ_len = 3.37/rad_len_cm/sin_ev
              musc_targ_len = musc_targ_len + .005*2.54/8.89/sin_ev
            endif


CCC  Next bit is to overide the target length seen by outgoing particles in the  CCC
CCC  tuna can target - MEC, 2005                                                 CCC


             call radlength_tuna(z,th_ev*degrad,musc_targ_len,tf_w)


             musc_targ_len = musc_targ_len/rad_len_cm + tf_w/8.89    !!! cryo + AL cellwall !!!


C Case 2 : solid target
          else
              musc_targ_len = abs(gen_lim(6)/2. - zt)/rad_len_cm/cos_ev
          endif


C Scattering before magnets:  Approximate all scattering as occuring AT TARGET.
C  16 mil Al scattering chamber window (X0=8.89cm)
C  15(?) cm air (X0=30420cm)  
C  4.5 mil Kevlar (X0=74.6cm)   !!!  2004  !!! 
C  2 x 2 mil Mylar  (X0=28.7cm) !!!  2004  !!!        Total of 0.60% rad. length.

          musc_targ_len = musc_targ_len + 0.016*2.54/8.89 +
     >          15./30420. + .0045*2.54/74.6 + 2*.002*2.54/28.7

C Begin transporting particle.  Note that the ".false." is the decay_flag
C (i.e. decay is always ignored), and "0.0" is fry (fast raster y position,
C which we should add at some point).

           ok_hms=.true.

*********** ADD SOME SUBROUTINES FOR THE TARGET FIELD CORRECTIONS TO THE FILE mc_hms.f ***************
c write(*,*) 'targ.mag.field',fieldon
          if ( fieldon) then
            p = p_spec*(1.+dpps/100.)
            call track_from_tgt(xs,ys,zs,dxdzs,dydzs,p*1000.,m2,spect,ok_hms)
c      write(*,*) 'after go through the field',ok_hms
          endif

          if (ok_hms) then
C drift back to zs = 0, the plane through the target center
c       
           xs = xs - zs * dxdzs
          ys = ys - zs * dydzs
             zs = 0.0d00
c     write(*,*) 'before go through the detec',ok_hms
          call mc_hms(p_spec,th_spec,dpps,xs,ys,zs,dxdzs,dydzs,
     >x_fp,dx_fp,y_fp,dy_fp,musc_targ_len,m2,cos_ts,
     >sin_ts,spect,ms_flag,wcs_flag,.false.,resmult,y,x,ok_hms)
c    write(*,*) 'after go through the detec',ok_hms
          endif
* G.&I.N. modified...
*

c          write(6,*) ' after call mc_hms'
          if (ok_hms.or.dump_all_in_ntuple) then !Success, increment arrays
             dpp_recon = dpps  
            dth_recon = dydzs*1000. !mr
           dph_recon = dxdzs*1000.

             ytar_recon = ys
            ztar_recon = (ytar_recon*cos_ts+xtar_init)
     >                /tan(th_spec-dth_recon)+ytar_recon*sin_ts

cc            write(6,*) dpp_recon,dth_recon,
cc     &                 dph_init,dph_recon


C Output NTUPLE entry.

           if (hut_ntuple) then
            hut(1) = x_fp
             hut(2) = y_fp
             hut(3) = dy_fp
               hut(4) = dx_fp
             hut(5) = ztar_init
             hut(6) = dpp_init
            hut(7) = dth_init
            hut(8) = dph_init
             hut(9) = ztar_recon
             hut(10)= dpp_recon
             hut(11)= dth_recon
             hut(12)= dph_recon
             hut(13)= where
              hut(14)= x_stop
              hut(15)= y_stop
              hut(16)= ys_init
              hut(17)= ytar_recon
              hut(18)= x
              hut(19)= y      
              call hfn(1,hut)
             endif

C Cut on reconstructed quantities.
           if ((abs(dpp_recon).gt.cut_dpp) .or.
     >(abs(dth_recon).gt.cut_dth) .or.
     >(abs(dph_recon).gt.cut_dph) .or.
     >(abs(ztar_recon).gt.cut_z)) then
               goto 500           !quit if failed
             endif

C Monte-Carlo trial is a 'GOOD' event, so histogram it.

            incut = incut + 1

C Get histogram indices. Use nint(x-.5) instead of int(x) to avoid
C doublecounts at x~=0, since int(-.999)=int(.999)=0.

            if (init_flag) then
             ix = nint( (dpp_init - h_xl)/h_dx - 0.5) + 1
               iy = nint( (dth_init - h_yl)/h_dy - 0.5) + 1
              iz = nint( (dph_init - h_zl)/h_dz - 0.5) + 1
             it = nint( (ztar_init - h_tl)/h_dt - 0.5) + 1
             else
              ix = nint( (dpp_recon - h_xl)/h_dx - 0.5) + 1
              iy = nint( (dth_recon - h_yl)/h_dy - 0.5) + 1
              iz = nint( (dph_recon - h_zl)/h_dz - 0.5) + 1
              it = nint( (ztar_recon - h_tl)/h_dt - 0.5) + 1
              endif

C Clip indices.  Note that this gives overflow bins.

            ix = min(max(ix,1),h_nx)
            iy = min(max(iy,1),h_ny)
            iz = min(max(iz,1),h_nz)
            it = min(max(it,1),h_nt)

C Increment histograms.
           cnts_p(ix) = cnts_p(ix) + 1
            cnts_th(iy) = cnts_th(iy) + 1
            cnts_phi(iz) = cnts_phi(iz) + 1
           cnts_ztgt(it) = cnts_ztgt(it) + 1
            cnts_zth(it,iy) = cnts_zth(it,iy) + 1
            cnts_pth(ix,iy) = cnts_pth(ix,iy) + 1
            cnts_pphi(ix,iz) = cnts_pphi(ix,iz) + 1
             cnts_thphi(iy,iz) = cnts_thphi(iy,iz) + 1
           cnts_pthph(ix,iy,iz) = cnts_pthph(ix,iy,iz) + 1

C Compute sums for calculating reconstruction variances.
           dpp_var(1) = dpp_var(1) + (dpp_recon - dpp_init)
            dth_var(1) = dth_var(1) + (dth_recon - dth_init)
           dph_var(1) = dph_var(1) + (dph_recon - dph_init) 
            ztg_var(1) = ztg_var(1) + (ztar_recon - ztar_init)

           dpp_var(2) = dpp_var(2) + (dpp_recon - dpp_init)**2
            dth_var(2) = dth_var(2) + (dth_recon - dth_init)**2
            dph_var(2) = dph_var(2) + (dph_recon - dph_init)**2
             ztg_var(2) = ztg_var(2) + (ztar_recon - ztar_init)**2
           endif !Incremented the arrays

C We are done with this event, whether GOOD or BAD.
C Loop for remainder of trials.

500        continue
c write(6,*) '1 event is over'
          enddo !End of M.C. loop

C------------------------------------------------------------------------------C
C                           End of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C


C        etime = 1

c	etime = secnds(stime)
c	type *,'Elapsed time = ',etime,' seconds'
c	type *,' '

C Close NTUPLE file.

        call hrout(1,i,' ')
        call hrend('HUT')
*
* G.&I.N. added
*
c open (unit=chan,status='unknown',name=my_name,
c     >	carriagecontrol='list')

        open (unit=chan,status='unknown',file=my_name)

* open (unit=chan,status='unknown',name='mc_hms.out',
*     >		carriagecontrol='list')

        write (chan,1002)
        write (chan,1003) p_spec,th_spec*degrad

!inp write (chan,1004) (gen_lim(i),i=1,6),
!inp     >	dpp_hist,dth_hist,dph_hist,h_xl,h_yl,h_zl,h_dx,h_dy,h_dz,
!inp     >  h_nx,h_ny,h_nz,h_entr,v_entr
        write (chan,1004) (gen_lim(i),i=1,6),
     >dpp_hist,dth_hist,dph_hist,h_xl,h_yl,h_zl,h_dx,h_dy,h_dz,
     > h_nx,h_ny,h_nz
       write (chan,1005) n_trials

C Indicate where particles are lost in spectrometer.

        write (chan,1015)
     >slit_hor,slit_vert,slit_oct,
     >Q1_in,Q1_mid,Q1_out,
     >Q2_in,Q2_mid,Q2_out,
     >Q3_in,Q3_out,
     >D1_in,D1_out

        write (chan,1006)
     >trials,shut,detec,successes,incut

         write(6,*)"check",trials,shut,detec,successes,incut

C Compute reconstruction resolutions.

        t1 = sqrt(max(0.,dpp_var(2)/incut - (dpp_var(1)/incut)**2))
        t2 = sqrt(max(0.,dth_var(2)/incut - (dth_var(1)/incut)**2))
        t3 = sqrt(max(0.,dph_var(2)/incut - (dph_var(1)/incut)**2))
        t4 = sqrt(max(0.,ztg_var(2)/incut - (ztg_var(1)/incut)**2))

        write (chan,1011) t1,t2,t3,t4

C Histogram dumps.

C$$$write (chan,1007) 'DPP-DTH'
C$$$write (chan,1008) ((cnts_pth(i,j),i=1,h_nx),j=1,h_ny)
C$$$write (chan,1007) 'DPP-PHI'
C$$$write (chan,1008) ((cnts_pphi(i,j),i=1,h_nx),j=1,h_nz)
C$$$write (chan,1007) 'DTH-PHI'
C$$$write (chan,1008) ((cnts_thphi(i,j),i=1,h_ny),j=1,h_nz)
!!!write (chan,1007) 'DPP-DTH-DPH'
!!!write (chan,1012) (((cnts_pthph(i,j,k),i=1,h_nx),j=1,h_ny),k=1,h_nz)
         close (unit=chan)

C Topdrawer plotfile output.

        if (topd_flag) then
c       open (unit=chan,status='unknown',name='mc_hms.top',
c     >		carriagecontrol='list')  
         open (unit=chan,status='unknown',file='mc_hms.top')
          write (chan,*)'set font duplex'

! Setup titles, etc. for DP/P histogram.

          write (chan,*) 'title top "DP/P histogram"'
          write (chan,*) 'title bottom "'//ttl_name(ttl_indx)

! Insert data points. Histogram, and begin new page.

          do i = 1,h_nx
          write (chan,*) h_xl+(i-0.5)*h_dx,cnts_p(i)
         enddo
c    write (chan,1001) 'HIST;NEW PAGE'

! Setup titles, etc. for DTH histogram.

          write (chan,*) 'title top "DTH  histogram"'
         write (chan,*) 'title bottom "'//ttl_name(ttl_indx)

! Insert data points. Histogram, and begin new page.

         do i = 1,h_ny
          write (chan,*) h_yl+(i-0.5)*h_dy,cnts_th(i)
          enddo
ccc  write (chan,1001) 'HIST;NEW PAGE'

! Setup titles, etc. for DPH histogram.

          write (chan,*) 'title top "DPHI  histogram"'
          write (chan,*) 'title bottom "'//ttl_name(ttl_indx)

! Insert data points. Histogram, and begin new page.

          do i = 1,h_nz
           write (chan,*) h_zl+(i-0.5)*h_dz,cnts_phi(i)
          enddo
c  write (chan,1001) 'HIST;NEW PAGE'

! Setup titles, etc. for ZTGT histogram.

          write (chan,*) 'title top "ZTGT  histogram"'
          write (chan,*) 'title bottom "'//ttl_name(ttl_indx)

! Insert data points. Histogram.

          do i = 1,h_nt
            write (chan,*) h_tl+(i-0.5)*h_dt,cnts_ztgt(i)
         enddo
           write (chan,1001) 'HIST;NEW PAGE'

! Setup titles, etc. for TH-PHI histogram.

          write (chan,*) 'title top "Theta-Phi  histogram"'
          write (chan,*) 'title bottom "'//ttl_name(ttl_indx)

! Insert data points. Histogram.

          write (chan,*) 'set storage x y z; set order x y z'
          do j = 1,h_ny
            do i = 1,h_nz
             if (cnts_thphi(i,j).gt.0)
     > write (chan,*) h_yl+(i-0.5)*h_dy,h_zl+(j-0.5)*h_dz,cnts_thphi(i,j)
          enddo
          enddo
ccc write (chan,1001) 'HIST'

! Topdrawer output finished.

         close (chan)
        endif

! Dump P-PHI histogram.

        if (pphi_flag) then
c  open (unit=chan,status='unknown',name='pphi.dat',carriagecontrol='list')
          open (unit=chan,status='unknown',file='pphi.dat')
        do i = 1,h_nx
          do j = 1,h_nz
             if (cnts_pphi(i,j).gt.0) write (chan,*)
     > h_xl+(i-0.5)*h_dx,h_zl+(j-0.5)*h_dz,cnts_pphi(i,j)
          enddo
         enddo
          close (chan)
       endif

! Dump Z-TH histogram.

        if (zth_flag) then
c open (unit=chan,status='unknown',name='zth.dat',carriagecontrol='list')
          open (unit=chan,status='unknown',file='zth.dat')
          do i = 1,h_nt
           do j = 1,h_ny
              if (cnts_zth(i,j).gt.0) write (chan,*)
     > h_tl+(i-0.5)*h_dt,h_yl+(j-0.5)*h_dy,cnts_zth(i,j)
            enddo
          enddo
          close (chan)
         endif

! Dump P-TH histogram.

        if (pth_flag) then
c	  open (unit=chan,status='unknown',name='pth.dat',carriagecontrol='list')
          open (unit=chan,status='unknown',file='pth.dat')
         do i = 1,h_nx
          do j = 1,h_ny
           if (cnts_pth(i,j).gt.0) write (chan,*)
     > h_xl+(i-0.5)*h_dx,h_yl+(j-0.5)*h_dy,cnts_pth(i,j)
          enddo
         enddo
          close (chan)
        endif

! Dump TH-PHI histogram.

        if (thphi_flag) then
c open (unit=chan,status='unknown',name='thphi.dat',carriagecontrol='list')
          open (unit=chan,status='unknown',file='thphi.dat')
          do i = 1,h_ny
           do j = 1,h_nz
            if (cnts_thphi(i,j).gt.0) write (chan,*)
     >  h_yl+(i-0.5)*h_dy,h_zl+(j-0.5)*h_dz,cnts_thphi(i,j)
            enddo
          enddo
          close (chan)
         endif

C ALL done!

c stop ' '
        stop

C =============================== Format Statements ============================

 1001   format(a)
 1002   format('!',/,'! Uniform illumination Monte-Carlo results')
 1003   format('!',/'! Spectrometer setting:',/,'!',/,
     >f8.3,' = P  spect (GeV)',/,
     >f8.3,' = TH spect (deg)')

 1004   format('!',/'! Monte-Carlo limits:',/,'!',/,
     >f8.3,' = GEN_LIM(1) - DP/P                    (half width, % )',/,
     >f8.3,' = GEN_LIM(2) - Theta                   (half width, mr)',/,
     >f8.3,' = GEN_LIM(3) - Phi                     (half width, mr)',/,
     >f8.3,' = GEN_LIM(4) - HORIZ (full width of 3 sigma cutoff, cm)',/,
     >f8.3,' = GEN_LIM(5) - VERT  (full width of 3 sigma cutoff, cm)',/,
     >f8.3,' = GEN_LIM(6) - Z                       (Full width, cm)',/,
     >'!',/,
     >'! Region Of Interest bounds (+/- these values):',/,'!',/,
     >f8.3,' = DPP_HIST (% )',/,
     >f8.3,' = DTH_HIST (mr)',/,
     >f8.3,' = DPH_HIST (mr)',/,
     >'!',/,
     >'! Histogram limits:',/,'!',/,
     >f8.3,' = Low edge  (DP, % )',/,
     >f8.3,' = Low edge  (DTH, mr)',/,
     >f8.3,' = Low edge  (DPH, mr)',/,
     >f8.3,' = Bin width (DP, % )',/,
     >f8.3,' = Bin width (DTH, mr)',/,
     >f8.3,' = Bin width (DPH, mr)',/,
     >i4,'    = Number of X bins',/,
     >i4,'    = Number of Y bins',/,
     >i4,'    = Number of Z bins',/,
     >'!',/
     >'! Fixed slit aperture:',/,'!')

!inp     >	,/,
!inp     >	g,' = Hor. 1/2 gap size (cm)',/,
!inp     >	g,' = Vert. 1/2 gap size (cm)')

 1005   format('!',/,'! Summary:',/,'!',/,i10,' Monte-Carlo trials:')

 1006   format(i10,' Initial Trials',/
     >i10,' Trials made it to the hut',/
     >i10,' Trials were cut in detectors',/
     >i10,' Trials made it thru the detectors and were reconstructed',/
     >i10,' Trials passed all cuts and were histogrammed.',/
     >)

c 1007 format('!',/,'! ',A,' histogram contents: ',/,'!')
c 1008 format(8i)
c 1009 format(1x,i4,g,i)
c 1010 format(a,i)
 1011 format(
     >'DPP resolution = ',f7.3,' %',/,
     >'DTH resolution = ',f7.3,' mr',/,
     >'DPH resolution = ',f7.3,' mr',/,
     >'YTG resolution = ',f7.3,' cm')

c 1012 format(1x,16i4)

 1015 format(/,
     >i8,' stopped in the FIXED SLIT HOR',/
     >i8,' stopped in the FIXED SLIT VERT',/
     >i8,' stopped in the FIXED SLIT OCTAGON',/
     >i8,' stopped in Q1 ENTRANCE',/
     >i8,' stopped in Q1 MIDPLANE',/
     >i8,' stopped in Q1 EXIT',/
     >i8,' stopped in Q2 ENTRANCE',/
     >i8,' stopped in Q2 MIDPLANE',/
     >i8,' stopped in Q2 EXIT',/
     >i8,' stopped in Q3 ENTRANCE',/
     >i8,' stopped in Q3 EXIT',/
     >i8,' stopped in D1 ENTRANCE',/
     >i8,' stopped in D1 EXIT',/
     >)
        end
