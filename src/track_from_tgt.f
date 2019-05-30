      subroutine track_from_tgt(x,y,z,dx,dy,mom,mass1,spect,ok)

C Given vertex coordinates, momentum and mass, tracks particle through a
C field to field-free region 100 cm from target.  It is assumed that the 
C code that calls this routine will reconstruct the track to z=0.
C
C x,y,z,dx,dy coordinates follow COSY a la Hall C:
C   z = into spectrometer
C   x = down (direction of increasing momentum)
C   y = z cross x.
C   dx = dx/dz 
C   dy = dy/dz
C
C coordinates for tracking are:
C   vT(1,2,3) are the position in X,Y,Z [cm]
C   vT(4,5,6) are the velocity in the X,Y,Z direction [cm/ns].


      implicit none


      real*8 x,y,z,dx,dy      ! in: initial coords.  out: coords of image track (in HMS cordinates)
      real*8 mom,one    ! momentum (MeV). (mom<0 for e-, mom>0 for p,d)
      real*8 mass1            ! mass squared
      real*8 mass             ! mass of particle (MeV)
      integer spect
      logical ok
 
      real*8 cc               ! speed of light in cm/ns
      parameter (cc = 29.9792458)
      real*8 vel              ! velocity of particle [cm/ns]
      real*8 eng              ! energy of particle
!      real*8 vT(6)
      real*8 vT(9)      ! OR - 4/04 (in HMS cordinates)
              
        logical forwd        ! OR - 4/04
        common/fwd/forwd        !        OR - 4/04
        integer flag_az        !        OR - 7/04
        common /azimuth/flag_az        !        OR 7-04
        
        logical b_corr
        real*8 field_off ! in cm
        real*8 x_init1,y_init1,z_init1,dx_init1,dy_init1
        real*8 x_init2,y_init2,z_init2,dx_init2,dy_init2
        real*8 x_final1,y_final1,z_final1,dx_final1,dy_final1
        real*8 x_final2,y_final2,z_final2,dx_final2,dy_final2
        real*8 projx,projy

        b_corr = .false.
        field_off = -2.0
c
c set forwd= false and flag_az=0
c  these were needed to be true and 1 
c   for RSS when target field pointing to 90deg.
c        forwd = .true.        
c        flag_az = 1  

        forwd = .false.        
        flag_az = 1 
        one=1.    
c

c        write(*,*) 'call to track_from_tgt'
c      write(*,*) 'from target',spect
c      write(*,*) mom,mass
c      call print_coord2('init track, beam: ',x,y,z,dx,dy)
      mass=sqrt(mass1)
      if (spect.lt.0) mom = -mom
      vel = abs(mom)/sqrt(mom**2+mass**2)*cc
      eng = sign(one,mom)*sqrt(mom**2+mass**2)
c
c
      vT(1) = x 
      vT(2) = y
      vT(3) = z 

      vT(6) = vel/sqrt(1+dx**2+dy**2)
      vT(4) = dx*vT(6)
      vT(5) = dy*vT(6)

! for debugging, run track first to z=0.
      ok = .true.
      call trgTrackToPlane(vT,eng,1.d00,0.d00,0.d00,1.d00,0.d00,ok,
     &     spect)
c      write(*,*) ' debug call to trgtracktoplane ' ,ok,spect
            
c      write(*,*) 'init track, z=0:',vt

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      x_init1 = vT(1) 
      y_init1 = vT(2)
      z_init1 = vT(3)
      dx_init1 = vT(4)/vT(6)
      dy_init1 = vT(5)/vt(6)

c      write(*,*) 'x,y,z,xp,yp at tgt1',x_init1,y_init1,z_init1,dx_init1,
c     &     dy_init1
c      write(*,*) 'vt 1',vT
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! track through magnetic field to z=100 cm plane
      ok = .true.
c      write(*,*) ' track_from_tgt = ',vt
c      write(*,*) 'x,y,z,xp,yp at trg',vT(1),vT(2),vT(3),
c     &     vT(4)/vT(6),vT(5)/vt(6)

      projx = vT(1)+100.0*vT(4)/vT(6)
      projy = vT(2)+100.0*vT(5)/vt(6)
c      write(*,*) 'proj. x/y at 100',projx,projy

      call trgTrackToPlane(vT,eng,1.d00,0.d00,0.d00,1.d00,-100.d00,
     &     ok,spect)
c      write(*,*) 'x,y,z,xp,yp at 100 cm',vT(1),vT(2),vT(3),
c     &     vT(4)/vT(6),vT(5)/vt(6)
c      write(*,*) ' after z=-100. call',ok,spect
c      write(*,*) 'init track, z=100:',vt
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccc add some b-field corrections ccccccccccccccccccccccccccccccccccc
c ys is for left in cm
      if(b_corr) then
         x_final1 = vT(1) 
         y_final1 = vT(2)
         z_final1 = vT(3)
         dx_final1 = vT(4)/vT(6)
         dy_final1 = vT(5)/vt(6)

c      write(*,*) 'x,y,z,xp,yp at 100cm for 1',x_final1,y_final1,
c     &        z_final1,dx_final1,dy_final1        

         vT(1) = x_init1+field_off
         vT(2) = y_init1
         vT(3) = z_init1
         vT(6) = vel/sqrt(1+dx_init1**2+dy_init1**2)
         vT(4) = dx_init1*vT(6)
         vT(5) = dy_init1*vT(6)      

         x_init2 = vT(1) 
         y_init2 = vT(2)
         z_init2 = vT(3)
         dx_init2 = vT(4)/vT(6)
         dy_init2 = vT(5)/vt(6) 
     
c      write(*,*) 'x,y,z,xp,yp at trg2',x_init2,y_init2,z_init2,
c     &        dx_init2,dy_init2
c      write(*,*) 'vt 2',vT

         ok = .true.
         call trgTrackToPlane(vT,eng,1.d00,0.d00,0.d00,1.d00,-100.d00,
     &        ok,spect)
         
         x_final2 = vT(1) 
         y_final2 = vT(2)
         z_final2 = vT(3)
         dx_final2 = vT(4)/vT(6)
         dy_final2 = vT(5)/vt(6)

c      write(*,*) 'x,y,z,xp,yp at 100cm for 2',x_final2,y_final2,
c     &        z_final2,dx_final2,dy_final2          

         x = x_final2 
         y = y_final2 - field_off
         z = z_init1

         dx = dx_final2
         dy = dy_final2   

c      write(*,*) 'pass',x,y,z,dx,dy
      else

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! translate back into SIMC variables

      x = vT(1)
      y = vT(2)
      z = vT(3)

      dx = vT(4)/vT(6)
      dy = vT(5)/vt(6)
      endif

      return
      end
     
