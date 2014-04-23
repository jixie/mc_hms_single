!      subroutine track_to_tgt(delta,y,dx,dy,frx,fry,mom,mass,ctheta,
!     >                        stheta,spect,ok,xfp,xpfp,yfp,ypfp)
!      subroutine track_to_tgt(delta,y,dx,dy,frx,fry,mom,mass,ctheta,
!     >                        stheta,spect,ok,xfp,xpfp,yfp,ypfp,xtgt)	! OR 4/04
      subroutine track_to_tgt(delta,y,dx,dy,frx,fry,mom,mass1,ctheta,
     >                   stheta,spect,ok,xfp,xpfp,yfp,ypfp,xtgt,bdl) ! OR 4/04

      implicit none

       real*8 delta,y,dx,dy! in: first guess reconstructed coords.
                                ! out: final reconstructed coords.
      real*8 frx                ! raster horizontal position (points right)
      real*8 fry                ! raster vertical position (points up)
      real*8 mom,one            ! momentum (MeV). (mom<0 for e-, mom>0 for p,d)
      real*8 mass1,mass             ! mass of particle (MeV)
      real*8 ctheta,stheta        ! cosine and sine of central spectrometer angle
      real*8 xfp,yfp,xpfp,ypfp
      integer spect
      logical ok

!      real*8  vT(6),vTx(6)
      real*8  vT(9),vTx(9) ! OR - 4/04
      real*8 xx,delx
      real*8 xxd
      real*8 save_delx,save_diff_delx
      integer*2 i,n
      real*8 vel,cc,eng,mom_0
      
      real*8 xtgt ! OR 4/04
      real*8 bdl ! OR - 4/04
      real*8 vtsave(6)
      integer ii
      logical forwd ! OR - 4/04
      common/fwd/forwd ! OR - 4/04
      integer flag_az ! OR - 7/04
      common /azimuth/flag_az ! OR 7-04	
      parameter (cc=29.9792458d00)
c        write(*,*) ' in track_to_tgt'
      forwd = .false. ! OR - 4/04 do the azimuthal angle correction for the backward
      flag_az = 1 ! OR - 7/04
      one=1.
! do reconstruction considering target field.  Taken from gen_track.f from 
! Markus Muehlbauer
!
        mass = sqrt(mass1)
       bdl = 0.0 ! OR - 4/04

c      write(*,*) 'to target',spect
c      call print_coord1('after first mc_hms_recon',y,delta,dx,dy,fry,0.)
c        write(*,*) 'call to track_tom_tgt'
c
c initialize xtgt = x_tar - OR 4/04
c
        xtgt = 0.d0
!
! copy vertical offset into another variable

      xx  = -fry

! use first call to mc_hms_recon as a first guess. Next trace back 100 cm 
! to enter field free region and calculate vector for field tracking program.
      if (spect.lt.0) mom = -mom
c      mom_0 = mom
      mom_0=  mom/(1.d0+delta/100.d0)
      vel = abs(mom)/sqrt(mom**2+mass**2)*cc
      eng = sign(one,mom)*sqrt(mom**2+mass**2)

c      mom_0 = mom/(1.d0+delta/100.d0)

c      write(*,*) 'vel,eng,mom_0 = ',vel,eng,mom_0,y,dx,dy
      vT(1) = -fry    + 100.d00*dx
      vT(2) = y + 100.d00*dy 
      vT(3) = 100.d00
      vT(6) = vel/SQRT(1+dy**2+dx**2)
      vT(4) = dx*vT(6)
      vT(5) = dy*vT(6) 

      do ii=1,6
         vtsave(ii)=vt(ii)
      enddo

c      call print_coord3('first track at z=100',vT)
      
! and track into the magnetic field to the beam plane (perp. to y)
c       write(*,*) ' Before 1st call trgtracktoplane ',vt
      ok = .true.
      CALL trgTrackToPlane (vT,eng,1.d00,0.d00,-ctheta,stheta,-frx,ok,
     &     spect)
c      CALL trgTrackToPlane (vT,eng,1.d00,0.d00,-ctheta,stheta,-100.d00,ok,
c     &     spect)

!write(*,*) ' Bdl - 1',vT(7),vT(8),vT(9),ok	!	OR - 4/04
      
c       write(*,*) ' 1st call trgtracktoplane ',-ctheta,stheta,frx,ok,vt
c      call print_coord3( 'first track on beam',vT)

      n  = 0
      delx = 1.
      save_diff_delx = 1.  
      save_delx=delx      
      DO WHILE ((delx .GT. .1) .AND. (n .LT. 10).and. (save_diff_delx .gt. 0) .AND. ok)
        delx = abs(-fry-vT(1))
c	 write(*,*) '----------------'
c	 write(*,*) 'do while: ',n,delx

        ! track to the z=0 plane to find a correction for the x-offset   

        vTx(1) = -fry 
        DO i=2,6
          vTx(i) = vT(i)
        ENDDO
    
c        call print_coord3( 'vT,  beam',vT)
c        call print_coord3( 'vTx, beam',vTx)
        CALL trgTrackToPlane (vT, eng,1.d00,0.d00,0.d00,1.d00,0.d00,
     &       ok,spect)

!write(*,*) ' Bdl - 2',vT(7),vT(8),vT(9),ok	!	OR - 4/04

        CALL trgTrackToPlane (vTx,eng,1.d00,0.d00,0.d00,1.d00,0.d00,
     &       ok,spect)

!write(*,*) ' Bdl - 3',vTx(7),vTx(8),vTx(9),ok	!	OR - 4/04

c        call print_coord3( 'vT,   z=0',vT)
c        call print_coord3( 'vTx,  z=0',vTx)

        xx = xx+min(1.,max(-1.,(vTx(1)-vT(1))))
        xxd = xx

        ! now find a better approximation 

c        call print_coord1('before mc_hms_recon',y,delta,dx,dy,-xxd,0.)

         call  simc_hms_recon (delta,dy,dx,y
     >                             ,xxd,xfp,xpfp,yfp,ypfp)
c           if (spect.lt.0) mom = -mom
c           mom_0 = mom
        mom = mom_0*(1.d0+delta/100.d0)
        vel = abs(mom)/sqrt(mom**2+mass**2)*cc
        eng = sign(one,mom)*sqrt(mom**2+mass**2)
        
c        call print_coord1('after mc_hms_recon',y,delta,dx,dy,-xxd,0.)

        ! drift to a field free region and calculate the velocities

        vT(1) = xx + 100.d00*dx
        vT(2) = y  + 100.d00*dy
        vT(3) = 100.d00
        vT(6) = vel/SQRT(1+dy**2+dx**2)
        vT(4) = dx*vT(6)
        vT(5) = dy*vT(6) 
     
        ! and track into the magnetic field to the beam plane (perp. to y)

c        call print_coord3( 'before last track',vT)
c
c 10/20/2003 change frx to -frx which gives better resolution.
c
c        write(*,*) n,' track to tgt ',vt
        CALL trgTrackToPlane (vT,eng,1.d00,0.d00,-ctheta,stheta,-frx,
     &       ok,spect)

        bdl = sqrt(vT(7)**2+vT(8)**2+vT(9)**2) !OR - 4/04
        delx = abs(-fry-vT(1))
        save_diff_delx = save_delx - delx
        if (save_diff_delx .lt. 0 .and.   n .ne. 0) then
           do ii=1,6
              vt(ii)=vtsave(ii)
           enddo
           ok = .false.
           if ( save_delx .le. 1.0) ok = .true.
        else
           n = n+1
           do ii=1,6
              vtsave(ii)=vt(ii)
           enddo
           save_delx = delx
        endif

!write(*,*) ' Bdl - 4',vT(7),vT(8),vT(9),ok,bdl !OR - 4/04

c        call print_coord3( 'after last track:',vT)
cc        n = n+1       
      ENDDO

      IF (n .ge. 10 ) ok = .FALSE.
      if (.not. ok) then
      endif      
cc      IF (delx .GT. .2) ok = .FALSE.
            
      ! calculate the result in HMS coordinates

      dy = vT(5)/vT(6)
      dx = vT(4)/vT(6)
      y  = vT(2)
      xtgt = vT(1) ! x_tar - OR 4/04   
      return
      end
