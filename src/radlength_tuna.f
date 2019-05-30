      SUBROUTINE radlength_tuna(ztar,theta,tf_h2,dt)

CCCCC   This subroutine calculates the amount of target material in cm          CCCCC
CCCCC   that is transversed in the cryogen and cell wall of the Hall C          CCCCC
CCCCC   tuna can target for an electron detected at a scattering angle, theta.  CCCCC
CCCCC   Returns length of cryogen transversed (tf_h2) and amount of             CCCCC
CCCCC   cell wall transversed (dt) in cm.                                       CCCCC
CCCCC   This doesn't currently include the beam on target position.             CCCCC


      IMPLICIT none

      integer*8 i
      real*8 pi,mp,mp2,radcon,dens,diam,R,R1,thick,theta,trad
      real*8 ztar,zmax,zint,xint,zdiff,li,lf(2),tf_h2,dt
      real*8 a,b,c,ctan,ctan2

      pi = 3.141592654
      radcon = pi/180.
      mp = .9382727
      mp2 = mp*mp
      dens = 0.0723
      diam = 3.96     !!!  Tuna can inner diameter  !!!
      thick = 0.0125  !!!  Thickness of cell wall   !!!
      R1 = diam/2.
      
c      write(6,*) "theta = ",theta 
  
      zmax = 3.96/2 - 1.e-4
      if(ztar.GT.zmax) ztar = zmax  !!! don't go outside target wall  !!!

      if(theta.EQ.90) theta = 89.999
      trad = radcon*theta
      ctan = 1/tan(trad)
      ctan2 = ctan*ctan
      a = 1. + ctan2
      b = 2.*ztar*ctan
     
      do i=1,2      !!!  Loop over inner and outer radius of cell wall  !!!
       R = R1
       if(i.eq.2) R = R1 + thick
       c = ztar*ztar - R*R
       xint = ((-1.)*b + sqrt(b*b - 4.*a*c))/2./a !!! z intersection of line and circle     !!!   
       zint = sqrt(R*R - xint*xint)  !!! x intersection of line and circle     !!! 
       if(xint.LT.0) then   !!!  choose other solution  !!!
         xint = ((-1.)*b - sqrt(b*b - 4.*a*c))/2./a !!! z intersection of line and circle   !!! 
         zint = sqrt(R*R - zint*zint)  !!! x intersection of line and circle     !!! 
       endif
       zdiff = zint - ztar           !!! z distance transversed                !!!  
       lf(i) = sqrt(xint*xint + zdiff*zdiff)  !!!  total distance transversed  !!!
      enddo
      li = ztar + R1
      dt = lf(2) - lf(1)

      tf_h2 = lf(1)

    
c 920   format(5f13.4)

      return
      end






