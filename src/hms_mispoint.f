        real*8 function hms_mispoint(theta)
        implicit none

* Original version made on 04/28/05 by E. Christy
* to model the measured HMS mis-pointing.  Interpolate 
* using measured ztar mispointing of HMS

         real*8 theta,thetab(22),yb(22),thetalow,thetahi,dtheta
        integer i

c       data thetab/ 6., 10., 14., 16., 18., 20., 22., 25., 26., 30., 32., 35., 40.,  45.,  
c       >               50.,  55.,   65.,  70., 80., 90. /
c       data yb/   0.6, 0.55, 0.4, 0.2, 0.1, 0.3,0.21,0.65,0.62, 1.0,1.09, 1.3, 1.43, 1.6, 
c       >               1.55, 1.2,   1.2,  1.4, 2.0, 2.8 / 


        data thetab/ 6., 10., 14., 16., 18., 20., 22., 25., 26., 30., 32., 35., 40.,  45.,                    
     >               50.,  55.,58., 60.,  65.,  70., 80., 90. /
        data yb/   0.6, 0.55, 0.4, 0.2, 0.1, 0.3,0.21,0.65,0.62, 1.0,1.09, 1.3, 1.43, 1.6, 
     >               1.55, 1.2,0., 0., 1.2,  1.4, 2.0, 2.8 /  

        thetahi = 1000.
        thetalow = 0.
        i=0

c        write(6,*) "theta = ",theta

        dowhile(thetahi.GE.100)
          i = i+1
          if(thetab(i).GT.theta) then
            thetahi = thetab(i)
            thetalow = thetab(i-1)
          endif
        enddo
        dtheta = thetahi - thetalow
        hms_mispoint = (thetahi-theta)*yb(i-1) 
     &                     + (theta-thetalow)*yb(i)
        hms_mispoint = hms_mispoint/dtheta
        hms_mispoint = hms_mispoint/10.    !!  put in cm  !!    enddo 
     
c        write(6,*) theta,thetalow,thetahi,dtheta,hms_mispoint


        return
        end




