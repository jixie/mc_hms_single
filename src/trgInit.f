* load the field map and calculate the magnetic field strength  
* 
      SUBROUTINE trgInit(map,theta_e,phi_e,theta_p,phi_p)
      IMPLICIT NONE
c      CHARACTER*17 map
      INTEGER*4 map
c      CHARACTER map*(*)
      REAL*8      theta_e,phi_e,theta_p,phi_p
* --  read field map (for calculations in the LAB system)
*
*     Parameter:
*        map        I : filename of the fieldmap (=' ': uniform field test case)
*        theta_e,phi_e I : inplane(theta) & out of plane(phi) angle for e spect
*        theta_p,phi_p I : inplane(theta) & out of plane(phi) angle for p spect
*        
*        note: currently phi is always treated as 0
*
* GAW 99/12/06 modified to work with 2 spectrometers
c      INTEGER trial1
      INTEGER    nz,nr 
      PARAMETER (nz = 51)
      PARAMETER (nr = 51)

      REAL*8    B_field_z(nz,nr),B_field_r(nz,nr),zz(nz),rr(nr)
      REAL*8    B_theta_e,B_stheta_e,B_ctheta_e,B_phi_e,B_sphi_e
      REAL*8    B_cphi_e 
      REAL*8    B_theta_p,B_stheta_p,B_ctheta_p,B_phi_p,B_sphi_p
      Real*8    B_cphi_p 
       
      COMMON  /trgFieldStrength/ B_field_z,B_field_r,zz,rr
      COMMON  /trgFieldAngles_e/ B_theta_e,B_stheta_e,B_ctheta_e,
     >                           B_phi_e,  B_sphi_e,  B_cphi_e 
      COMMON  /trgFieldAngles_p/ B_theta_p,B_stheta_p,B_ctheta_p,
     >                           B_phi_p,  B_sphi_p,  B_cphi_p 
 
      REAL*8       pi180
      real*8 scale ! hard coded rescaling of B field for RSS
      parameter (scale=0.98104) ! B_RSS = 5.0033 T - OARA 4/12/04
c	parameter (scale=1.15) ! B_RSS = 5.0033 T - OARA 4/12/04
      PARAMETER (pi180 = 3.141592653d00/180.d00) 

      INTEGER ir,iz 
      REAL*8    xx
  
      B_theta_e  = theta_e
      B_stheta_e = SIN(theta_e*pi180)*cos(phi_e*pi180)  
      B_ctheta_e = COS(theta_e*pi180)

      B_theta_p  = theta_p
      B_stheta_p = SIN(theta_p*pi180)*cos(phi_p*pi180)  
      B_ctheta_p = COS(theta_p*pi180)

      ! Note: for performance reasons B_phi is always treated 0 in trgField
      B_phi_e    = phi_e
      B_sphi_e   = SIN(phi_e*pi180) 
      B_cphi_e   = COS(phi_e*pi180)

      B_phi_p    = phi_p
      B_sphi_p   = SIN(phi_p*pi180) 
      B_cphi_p   = COS(phi_p*pi180)
CGAW      write(*,*) 'trginit',theta_e,theta_p
! GAW 99/11/22: Add zero field option
   
c      map ='trg_field_map.dat'
      map = 1
   
c      trial1=0
c      IF (map.EQ.'0') THEN
      IF (map.eq.0) THEN
  
         DO ir=1,nr
            rr(ir) = 2.0*float(ir-1)
            zz(ir) = 2.0*float(ir-1)
            DO iz=1,nz
               B_field_r(iz,ir) = 0.
               B_field_z(iz,ir) = 0.
            ENDDO
         ENDDO
! GAW 99/11/22: End
    
c      ELSEIF (map.NE.' ') THEN !read in numerical field map
      ELSEIF (map.EQ.1) THEN !read in numerical field map
         OPEN (unit=1,file='trg_field_map.dat',status='old')
c         OPEN (unit=1,file=map,status='old')
         write(6,*) ' read the numerical field map trg_field_map.dat'
         DO ir=1,nr
            rr(ir) = 2.0*float(ir-1)
            zz(ir) = 2.0*float(ir-1)
            DO iz=1,nz
               READ (1,*)xx,xx,B_field_z(iz,ir),B_field_r(iz,ir),
     &              xx,xx,xx
c               write(6,*)'after reading the file'
c               write(6,*) iz,ir,B_field_z(iz,ir),B_field_r(iz,ir) 
c               trial1=trial1+1
c               write(6,*) trial1
   ! rescale field to desired value
               B_field_z(iz,ir) = B_field_z(iz,ir) * scale
               B_field_r(iz,ir) = B_field_r(iz,ir) * scale
c               write(6,*) B_field_z(iz,ir),B_field_r(iz,ir)
               
            ENDDO
         ENDDO
         CLOSE (unit=1)
      ELSE

! GAW 99/11/19: Must initialize rr and zz before going through loop since 
! do tests on them.

        DO ir=1,nr ! uniform 5T field over 26 cm in z
          rr(ir) = 2.0*float(ir-1) ! and 16 cm in r
          zz(ir) = 2.0*float(ir-1)
        ENDDO
        DO ir=1,nr ! uniform 5T field over 26 cm in z
          DO iz=1,nz
            B_field_r(iz,ir) = 0.
            IF (rr(ir) .LE. 16. .and. zz(iz) .LE. 26.) THEN
CGAW              B_field_z(iz,ir) = 0.0
              B_field_z(iz,ir) = 5.0
            ELSE
             B_field_z(iz,ir) = 0.0
           ENDIF
          ENDDO
        ENDDO
      ENDIF
 
      RETURN
      END
      
