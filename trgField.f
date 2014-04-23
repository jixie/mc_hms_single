      SUBROUTINE trgField (x_,B_,spect)
      IMPLICIT NONE
      
      include 'addmore.cmn'
      
!      REAL*8 x_(3),B_(3)
      REAL*8 x_(6),B_(3)
* --  calculate actual field
*
*     Parameter:
*        x_   I : lab coordinates  (cm)
*        B_   O : B field in lab coordinates
*        spect I: id for spectrometer (-1=e-, +1=p)
*
*      Notes:
*      - 2-Dimensional Linear Interpolation:                               
*        Assumes uniform spacing of fieldmap in x,y        
*      - for performance reasons B_phi is always treated 0 
*
* GAW 99/12/06 modified to work with 2 spectrometers
       
      INTEGER    nz,nr 
      PARAMETER (nz = 51)
      PARAMETER (nr = 51)

      REAL*8    B_field_z(nz,nr),B_field_r(nz,nr),zz(nz),rr(nr)
      REAL*8    B_theta_e,B_stheta_e,B_ctheta_e,B_phi_e,B_sphi_e
      REAL*8    B_cphi_e 
      REAL*8    B_theta_p,B_stheta_p,B_ctheta_p,B_phi_p,B_sphi_p
      REAL*8    B_cphi_p 
       
      COMMON  /trgFieldStrength/ B_field_z,B_field_r,zz,rr
      COMMON  /trgFieldAngles_e/ B_theta_e,B_stheta_e,B_ctheta_e,
     >                           B_phi_e,  B_sphi_e,  B_cphi_e 
      COMMON  /trgFieldAngles_p/ B_theta_p,B_stheta_p,B_ctheta_p,
     >                           B_phi_p,  B_sphi_p,  B_cphi_p 

      REAL*8 B_tht, B_stht, B_ctht, B_ph, B_sph, B_cph

      INTEGER i,j,spect
      REAL*8    x(3),B(3),z,r,az,ar,a0,a1

        real*8 x_sur(2),y_sur(2),z_sur(2)        !        OR - 4/04
        real*8 azim,az0,az_corr,horiz       !        OR - 7/04
        real*8 B_scale,B_corr                        !        OR - 7/04        
        logical forwd        !        OR - 4/04 B offset forward only
        common/fwd/forwd        !        OR -4/04
        integer flag_az        !        OR - 7/04
        logical azcor
        common /azimuth/ flag_az        !        OR - 7/04       
        save B_scale
 
!        data flagaz /1/
c
c Survey system z=beam downstream, x=left of beam y= above beam
c        
!        data    x_sur /-.161d0,  .123d0/, 
!        1        y_sur / .147d0,  .108d0/, 
!        1        z_sur / .339d0, -.205d0/! Survey offsets cm, 1= perp
!!        1        z_sur / .197d0, -.205d0/! Survey offsets cm, 1= perp
        
        data    x_sur / .0d0,  .0d0/ 
        data    y_sur / .0d0,  .0d0/ 
        data    z_sur / .0d0,  .0d0/ ! Survey offsets cm, 1= perp

! set up angles for electron or proton spectrometer

!        if(flag_az.eq.1.and.forwd) then        !        OR - 7/04
c                 write(*,*) 'b_theta_e',b_theta_e,b_theta_e/180.*3.14159
c     &                ,b_theta_e-80
          horiz=atan(x_(5)/x_(6))
c          azim=atan(x_(4)/(x_(6)*cos(horiz+b_theta_e/180.*3.14159)))
          azim=atan(x_(4)/(x_(6)*cos(horiz+5.41/180.*3.14159)))


!          azim=x_(4)/x_(6)
!          flag_az = 0
!         end if
        
      if (spect.eq.-1) then
        B_tht    = B_theta_e
        B_stht   = B_stheta_e
        B_ctht   = B_ctheta_e
        B_ph     = B_phi_e
        B_sph    = B_sphi_e
        B_cph    = B_cphi_e
      else if (spect.eq.1) then
        B_tht    = B_theta_p
        B_stht   = B_stheta_p
        B_ctht   = B_ctheta_p
        B_ph     = B_phi_p
        B_sph    = B_sphi_p
        B_cph    = B_cphi_p
      endif

c          write(*,*) 'b_theta',B_tht,B_stht,B_ctht,B_ph,B_sph,B_cph
c
c Offset the center of the field: magnet center not on pivot
c

!        xx(1) = x_(1) - 0.147
!        xx(1) = x_(1) - 0.847
!        xx(2) = x_(2) - 0.112
!        xx(2) = x_(2) - 0.112
!        xx(3) = x_(3) + 0.228
!        xx(3) = x_(3) + 0.339

        
      ! rotate to coordinates with z' along field direction
c      write(*,*)'x(1) =', x(1),x_(1),x_(1)-0.43

c X(1) pointing down, X(2) pointing left and X(3) pointing beam down.
c      write(*,*)'theta_lab',theta_lab,y_off,x_off
c      write(*,*)'before x_(1),x_(2),x_(3)',x_(1),x_(2),x_(3)
c      x_(1)=x_(1)-sry_off
c      x_(2)=x_(2)-srx_off*sin(theta_lab)
c      x_(3)=x_(3)+srx_off*cos(theta_lab)
c      write(*,*)'after x_(1),x_(2),x_(3)',x_(1),x_(2),x_(3)      
! sry_off pointing up, srx_off pointing right
      x(1) =           x_(1)! + sry_off
      x(2) =  B_stht*x_(3) + B_ctht*x_(2) !- (beam_offset*B_ctheta_e)
      x(3) =  B_ctht*x_(3) - B_stht*x_(2) !+ srx_off !- (beam_offset*B_stheta_e) 

c      write(*,*)'x_(1),x_(1) + sry_off',x_(1),x_(1) + sry_off
c      write(*,*)'x(2)',x(2)-B_ctht*srx_off,x(2)
c      write(*,*)'x(3)',x(3)+B_stht*srx_off,x(3)
c      write(*,*)'x(1),x(2),x(3)',x(1),x(2),x(3)

c      x(1) =           x_(1)+0.24
c      x(2) =  B_stht*x_(3) + B_ctht*x_(2) + 0.24*B_ctheta_e
c      x(3) =  B_ctht*x_(3) - B_stht*x_(2) + 0.24*B_stheta_e 

 
!      x(1) =           xx(1)
!      x(2) =  B_stht*xx(3) + B_ctht*xx(2)
!      x(3) =  B_ctht*xx(3) - B_stht*xx(2)  

      ! compute zylinder coordinates

        if(forwd) then        ! to hut OR - 4/04 
            
         if(abs(B_tht).gt.70.d0) then        !        perp
!          xx(1) = x_(1) - y_sur(1)
!          xx(2) = x_(2) + x_sur(1)
!          xx(3) = x_(3) + z_sur(1)
          z = ABS  (x(3) + x_sur(1))
          r = SQRT ((x(1)+y_sur(1))**2 + (x(2)+z_sur(1))**2)
         else                        !        para
!          xx(1) = x_(1) - y_sur(2)
!          xx(2) = x_(2) + x_sur(2)
!          xx(3) = x_(3) + z_sur(2)
          z = ABS  (x(3)) + z_sur(2)
          r = SQRT ((x(1)+y_sur(2))**2 + (x(2)+x_sur(2))**2)
!          z = ABS  (x(3)) + x_sur(2)
!          r = SQRT ((x(1)+y_sur(2))**2 + (x(2)+Z_sur(2))**2)
         end if                      
        
        else        ! from hut
        
!         xx(1) = x_(1)
!         xx(2) = x_(2)
!         xx(3) = x_(3)
        z  = ABS  (x(3))
        r  = SQRT (x(1)**2 + x(2)**2)
                
        end if      
        
!        z  = ABS  (x(3))
!        r  = SQRT (x(1)**2 + x(2)**2)
        
        
      ! interpolate the field map 
      i = INT((z-zz(1))/(zz(2)-zz(1))) + 1                                          
   
      j = INT((r-rr(1))/(rr(2)-rr(1))) + 1                                          
   
      IF ((i+1 .GT. nz) .OR. (i .LT. 1) .OR. 
     >    (j+1 .GT. nr) .OR. (j .LT. 1)) THEN                
!        Missing initialization of field - inserted by OR 4/04
        B(1)=0.
        B(2)=0.
        B(3)=0.

        B_(1)=0.
        B_(2)=0.
        B_(3)=0.
      ELSE                                                                     
        ! calculate the Bz component 
        az = ((z-zz(i))/(zz(2)-zz(1))) 
        ar = ((r-rr(j))/(rr(2)-rr(1))) 
        a0=az*(B_field_z(i+1,j)  -B_field_z(i,j))  +B_field_z(i,j)                  
                        
        a1=az*(B_field_z(i+1,j+1)-B_field_z(i,j+1))+B_field_z(i,j+1)                
                          
        B(3) = (ar*(a1-a0)+a0)           
        IF (r .gt. 0.) THEN
          ! calculate the Bx,By components 
          a0=az*(B_field_r(i+1,j)  -B_field_r(i,j))  +B_field_r(i,j)                
                          
          a1=az*(B_field_r(i+1,j+1)-B_field_r(i,j+1))+B_field_r(i,j+1)              
                            
          B(2) = (ar*(a1-a0)+a0)/r
          IF (x(3) .LT. 0.) B(2)= -B(2)
          B(1) = B(2)*x(1)
          B(2) = B(2)*x(2)       
c
                   
c
c New azim corr
c
c        azcor = .true.
        azcor = .false.        

c        if(B_tht.lt.0) then        ! B_tht = -103.15 deg = perp
c        azcor = .true.        ! = -90 - 13.15 (theta_e)
c        else if(B_tht.gt.0) then ! B_tht = 166.85 deg = para
c        azcor = .false.            ! = 180 - 13.15
c        end if
        
!no fwd        if(azcor.and.forwd) then

        if(azcor) then
        az0 =0.06d0        !        OR - 7/04
        az_corr =-0.06d1
c        az0 =-0.75d0        !        OR - 7/04
c        az_corr =0.05d1
c        B_corr = (azim-az0)*az_corr
           B_corr = (azim-az0)*az_corr
           if(flag_az.eq.1.and.forwd) then
              B_scale = 5.003/(5.003+abs(B_corr))
!        B_scale = B(3)/(B(3)+B_corr)        ! B(3) = B_max = 5.003 the 1st. time
              flag_az = 0
           else
              B_scale = 1.d0
           end if
        
c     mkj temporary to not correct for forward azimuthal distortion
           if ( .not. forwd) then
              B_scale=1.0
              b_corr=0.0
           endif
c
c        if (forwd)   write(*,*) ' corr=',azim,az0,az_corr,b_corr,b_scale
c     >  ,b_theta_e,horiz
           B(3) = B(3)+B_corr
           B(3) = B(3)*B_scale
!        1        +(azim-az0)**2*az_corr2        ! quadratic
c           B(2) = B(2)+B_corr
c           B(2) = B(2)*B_scale
           B(1) = B(1)-B_corr   !B(1) is x component pointing down
           B(1) = B(1)*B_scale
        end if
!no fwd        end if
        
          ! transform B field to lab. system
          B_(1) =          B(1)
          B_(2) = (- B_stht)*B(3) + B_ctht*B(2)
          B_(3) =   B_ctht*B(3) + B_stht*B(2)  
        ELSE  
          B_(1) =   0.
          B_(2) = (- B_stht)*B(3)
          B_(3) =   B_ctht*B(3)
        ENDIF
      ENDIF           
c       write(*,*) r,z,x_,' b = ',b_
      RETURN
      END
