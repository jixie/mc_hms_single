      SUBROUTINE radlength(ztar,theta,tf,tf_w)
CCCCC   This program calculates the average amount of target material    CCCCC
CCCCC   in radiation lengths that is transversed in the Hall C tuna can  CCCCC
CCCCC   target for an electron detected at a scattering angle, theta.    CCCCC

      IMPLICIT none

      integer*8 i,j
      real*8 pi,mp,mp2,radcon,dens,diam,R,R1,thick,theta,trad
      real*8 ztar,zint,xint,zdiff,xdiff,tan2,li,lf(2),dt
      real*8 ti,tf,a,b,c,sign,b_h,b_al,b_my,b_air   
      real*8 ti_w,ti_e,tf_w,tf_air,tf_spectwin,ti_tot 

      pi = 3.141592654
      radcon = pi/180.
      mp = .9382727
      mp2 = mp*mp
      dens = 0.0723
      diam = 3.96
      thick = 0.0125
      R1 = diam/2.
      sign = +1.
      
      if(theta.EQ.90) theta = 89.999
      if(theta.GT.90.0) sign = -1.
      trad = radcon*theta
      tan2 = tan(trad)*tan(trad)
      a = 1. + tan2
      b = -2.*ztar*tan2
     
      do i=1,2
       R = R1
       if(i.eq.2) R = R1 + thick
       c = ztar*ztar*tan2 - R*R
       zint = (-1.*b + sign*sqrt(b*b - 4.*a*c))/2./a            
       xint = sqrt(R*R - zint*zint) 
       zdiff = zint - ztar
       lf(i) = sqrt(xint*xint + zdiff*zdiff)
      enddo
      li = ztar + R1
      dt = lf(2) - lf(1)

      call bcalc(1.0,b_h)
      call bcalc(13.0,b_al)
      call bcalc(5.2,b_my)
      call bcalc(7.3,b_air)

      ti = li*0.0723/61.280               !!!  Hydrogen in       !!!
      tf = lf(1)*0.0723/61.280            !!!  Hydrogen out      !!!

      call bcalc(13,b_al) 
      ti_w = 0.0125*2.68/23.6311          !!!  Al cell wall in   !!!
      ti_e = 0.0406*2.68/23.6311          !!!  Entrance window   !!!

      tf_w = dt*2.68/23.6311              !!!  Al cell wall out  !!!
      tf_air = 30.0*0.00121/36.66
      tf_spectwin = 0.0127*1.39/39.95 +   !!!  Spect. window     !!!
     &              0.0432*0.75/55.20     !!!  Mylar + kevlar    !!!

      ti_tot = ti + ti_w + ti_e

c      write(6,*) ztar,R1,ztar + R1

c      write(6,*) xint,zint,zdiff,lf(1),dt
c       write(6,*) "ti = ",ti,"  tf = ",tf, " tf_aal = ",tf_w 
    
 920   format(5f13.4)

      return
      end






