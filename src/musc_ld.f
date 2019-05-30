          subroutine musc_ld(m2,p,Z,A,X,dth,dph)
C+_____________________________________________________________________
!
! MUSC_LD - Simulate multiple scattering of any particle.
!
! ASSUMPTIONS: DTH and DPH given in milli-radians.  The formula assumes a
!   gaussian distribution for the scattering angle. It is further assumed
!   that the angles DTH and DPH are located at an angle with respect 
!   to the beam direction that is large enough so that DTH-DPH space can
!   be approximated as a cartesian space.
!
! M.E. Christy - 04/29/05
!
! Formulas are from Lynch and Dahl, NIM B58 (1991) 6, which reproduce 
! Moliere distributions to ~2%.  
!
C-_____________________________________________________________________

        implicit none

        real*8 F,alpha
        parameter (F = 0.990)
        parameter (alpha = 1/137.03599)

        real*8 A,Z,X,Xc2,Xa2,v,omega,dth,dph
        real*8 beta, theta_sigma
        real*8 m2, p, ptemp,ddt

        real*8 nsig_max
        parameter(nsig_max=99.0d0)   !!!  max #/sigma for gaussian ran #s.  !!!

        real*8 gauss1

! Compute scattering angles, THETA_SCAT from a gaussian distribution,
! PHI_SCAT from uniform distribution.


        ptemp = p*1.0e3         ! change from GeV/c to MeV/c  !

        beta = ptemp/sqrt(m2+ptemp*ptemp)

        xc2 = 0.157*(Z*(Z+1.)/A)*X/(ptemp*beta)**2
        xa2 = 2.007e-5*(Z**(2./3.))*(1.+3.34*
     &            (Z*alpha/beta)**2)/ptemp/ptemp
        omega = xc2/xa2
       
c        write(6,*) x,z,a
c        write(6,*) m2,ptemp,p,beta,xc2,xa2,omega

        v = 0.5*omega/(1.-F) 

        theta_sigma = Xc2/(1.+F*F)*((1.+v)/v*log(1.+v)-1.)
        theta_sigma = sqrt(theta_sigma)


c       write(6,*) "LD: ",theta_sigma

! Compute new trajectory angles (units are rad)
       
c        write(6,*) "before: ",dth

        ddt = dth
        dth = dth + theta_sigma * gauss1(nsig_max)
        ddt = dth-ddt

c        write(6,*) "1: ",ddt

        ddt = dph
        dph = dph + theta_sigma * gauss1(nsig_max)
        ddt = dph-ddt

c        write(6,*) "2: ",ddt
c        write(6,*) theta_sigma*gauss1(nsig_max)

        return
        end

