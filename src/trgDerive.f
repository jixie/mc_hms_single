* solve the differential equation of the particle  
*
      SUBROUTINE trgDeriv(u,dudt,spect)
      IMPLICIT NONE
      REAL*8 u(9),dudt(9)
* --  calculate the derivatives du(i)/dt for the runke kutta routine         
*
*     Parameter:
*       u     I : actual coordinate vector
*                   u(1,2,3)    I : x, y, z
*                   u(4,5,6)    I : dx/dt, dy/dt, dz/dt 
*                   u(7,8,9)    I : integral Bxdx, Bydy, Bzdz   
*       dudt  O : derivative du/dt
*                   dudt(1,2,3) : dx/dt, dy/dt, dz/dt 
*                   dudt(4,5,6) : d^2xdt^2, d^2ydt^2, d^2zdt^2
*                   dudt(7,8,9) : B x v
*       spect I : -1 for e spectrometer, +1 for p spectrometer

      REAL*8   factor
      COMMON /trgConversionFactor/factor
      INTEGER spect
      REAL*8   B(3)

      CALL trgField (u,B,spect)
c      write(*,*) 'B(1), B(2),B(3)',B(1),B(2),B(3)
      ! These are just the velocities
      dudt(1) = u(4)
      dudt(2) = u(5)
      dudt(3) = u(6)

      ! This is just (v_vec X B_vec)  
      dudt(7) = u(5)*B(3) - u(6)*B(2)
      dudt(8) = u(6)*B(1) - u(4)*B(3) 
      dudt(9) = u(4)*B(2) - u(5)*B(1)  

      ! This is just (v_vec X B_vec) * factor
      dudt(4) = dudt(7)*factor
      dudt(5) = dudt(8)*factor
      dudt(6) = dudt(9)*factor

      RETURN
      END
                                                       
