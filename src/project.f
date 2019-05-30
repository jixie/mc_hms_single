         subroutine project(x_new,y_new,z_drift,decay_flag,dflag,m2,ppi)

C+______________________________________________________________________________
!
! PROJECT - Calculate new track transverse coordinates after drifting in a field
!   free region for a distance z_drift from current track location.
!
! Added 9/15/98: Check for decay of particle. dflag is true of the particle
! has already decayed, so check for decay if dflag .eq. .false.  Take into
! account additional path length for rays not parallel to the z axis.
!
C-______________________________________________________________________________

        implicit none

        include 'track.inc'

        real*8 x_new,y_new,z_drift
        logical*4 decay_flag !check for decay
        logical*4 dflag !has particle decayed yet?

C Local declarations.
        real*8 ppi,m2
!dec    real*8 z_decay
!dec    real*8 rph,rth1,rth
!dec    real*8 beta,gamma,dlen
!dec    real*8 ef,pf,pxf,pyf,pzf,pxr,pyr,pzr
!dec    real*8 bx,by,bz,er,pr

!dec    real*8 grnd

C ============================= Executable Code ================================


C Check for decay of particle.

!dec    if (.not.decay_flag .or. dflag) then  !particle has already decayed, just drift

          x_new = x_new + dxdzs * z_drift
          y_new = y_new + dydzs * z_drift

!dec    else !check for decay
!dec
!dec     p_spec = ppi/(1.+dpps/100.)
!dec     beta = ppi/sqrt(ppi**2+m2)
!dec     gamma = 1./sqrt(1.-beta*beta)
!dec     dlen=ctau*beta*gamma   !1/e decay length (beta*c*tau*gamma)
!dec     if (z_drift.le.0) write(6,*) 'drift distance<0:  automatic decay! BAD!'
!dec
!dec     z_decay = -1.*dlen*log(1-grnd())
!dec
!decC Check if the decay is within the drift length (i.e. the drift lenght in
!decC z, times sqrt(1+dxdzs**2+dydzs**2) to correct for the true path of the ray.
!dec
!dec     if(z_decay .gt. z_drift*sqrt(1+dxdzs**2+dydzs**2)) then !no decay 
!dec
!dec      decdist(30)=decdist(30)+z_drift
!dec     x_new = x_new + dxdzs * z_drift
!dec     y_new = y_new + dydzs * z_drift
!dec
!dec   else		!DECAY.  Find out where and generate decay particle.
!dec
!dec     dflag = .true.
!dec     decdist(30)=decdist(30)+z_decay
!dec     x_new = x_new + dxdzs * z_decay
!dec     y_new = y_new + dydzs * z_decay
!dec
!decC Generate center/mass decay angles and momenta.
!dec     rph = grnd()*2.*pi
!dec     rth1 = grnd()*2.-1.
!dec     rth = acos(rth1)
!dec     er = 109.787
!dec     pr = 29.783
!dec     pxr = 29.783*sin(rth)*cos(rph)
!dec     pyr = 29.783*sin(rth)*sin(rph)
!dec     pzr = 29.783*cos(rth)
!dec     m2 = 105.67 **2	!need mass-squared for multiple scattering.
!dec     Mh2_final = m2	!for ntuple
!dec
!dec
!decC Boost to Lab frame, calculate new angles and momentum, finish drift
!dec
!dec     bx = beta * dxdzs / sqrt(1. + dxdzs**2 + dydzs**2)
!dec     by = beta * dydzs / sqrt(1. + dxdzs**2 + dydzs**2)
!dec     bz = beta *   1.  / sqrt(1. + dxdzs**2 + dydzs**2)
!dec     call loren(gamma,bx,by,bz,er,pxr,pyr,pzr,ef,pxf,pyf,pzf,pf)
!dec     dxdzs = pxf/pzf
!dec     dydzs = pyf/pzf
!dec     dpps = 100.*(pf/p_spec-1.)
!dec     ppi=pf
!dec     x_new = x_new + dxdzs * (z_drift-z_decay)
!dec     y_new = y_new + dydzs * (z_drift-z_decay)
!dec
!dec   endif	  !if decayed
!dec endif	!if need to check for decay

        return
        end
