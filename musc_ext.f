         subroutine musc_ext(m2,p,rad_len,x_len,dth,dph,y,x)
C+_____________________________________________________________________
!
! MUSC - Simulate multiple scattering of any particle.
!        Used for extended scatterers.
!
! According to Particle Data Booklet, July 1994
!
C-_____________________________________________________________________

        implicit none

        real*8 Es, epsilon
        parameter (Es = 13.6)  !!!  MeV  !!!
        parameter (epsilon = 0.088)

        real*8 rad_len, x_len, dth, dph, x, y
        real*8 beta, g1, g2, theta_sigma
        real*8 m2, p, ptemp

        real*8 nsig_max
        parameter(nsig_max=99.0d0)     !max #/sigma for gaussian ran #s.

        real*8 gauss1

        if (x_len.le.0.0d0 .or. rad_len.le.0.0d0) then
          write(6,*) 'x_len or rad_len < 0 in musc_ext.f'
c   write(6,*) 'This is bad.  Really bad.  Dont even ask how bad it is.'
        write(6,*) 'Just fix it now.'
          stop
        endif


        ptemp = p*1.0e3 ! change from GeV/c to MeV/c - MEC !

        beta = ptemp/sqrt(m2+ptemp*ptemp)

        theta_sigma = Es/ptemp/beta*sqrt(rad_len)*(1.+epsilon*log10
     &  (rad_len))

! Compute new trajectory angles and displacements (units are rad and cm)

        g1 = gauss1(nsig_max) ! gaussian, truncated at 10 sigma
        g2 = gauss1(nsig_max)
        dth = dth + theta_sigma*g1
        x   = x   + theta_sigma*x_len*g2/sqrt(12.0d0) + 
     &              theta_sigma*x_len*g1/2.0d0

        g1 = gauss1(nsig_max) ! gaussian, truncated at 10 sigma
        g2 = gauss1(nsig_max)
        dph = dph + theta_sigma*g1
        y   = y   + theta_sigma*x_len*g2/sqrt(12.0d0) + 
     &              theta_sigma*x_len*g1/2.0d0

        return
        end





