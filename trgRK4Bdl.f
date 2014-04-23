      SUBROUTINE trgRK4Bdl(u0,u1,h,spect)
      IMPLICIT NONE
      REAL*8     u0(9),u1(9),h
* --  Fourth-order Runge-Kutta from Numerical Recipes book
*     for tracking through the target field (incl. B/dl calculation)
*
*     Parameter:
*      u0  I  : input  coordinate vector
*      u1  O  : output coordinate vector
*                 u(1,2,3) : x, y, z
*                 u(4,5,6) : dx/dt, dy/dt, dz/dt 
*                 u(7,8,9) : integral Bxdx, Bydy, Bzdz   
*      h   I  : time step
*      spect I: -1 for e spectrometer, +1 for p spectrometer

      INTEGER i,spect
      REAL*8    ut(9),dudt(9),dut(9),dum(9),hh,h6

      hh=h*0.5
      h6=h/6.

      CALL trgDeriv(u0,dudt,spect)
      DO i=1,9
       ut(i) = u0(i) + hh*dudt(i)
      ENDDO

      CALL trgDeriv(ut,dut,spect)
      DO i=1,9
       ut(i) = u0(i) + hh*dut(i)
      ENDDO

      CALL trgDeriv(ut,dum,spect)
      DO i=1,9
      ut(i) = u0(i) +h*dum(i)
      dum(i)= dut(i)  +dum(i)
      ENDDO

      CALL trgDeriv(ut,dut,spect)
      DO i=1,9
        u1(i)=u0(i)+h6*(dudt(i)+dut(i)+2.*dum(i))
      ENDDO

      RETURN       
      END
