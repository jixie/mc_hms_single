      SUBROUTINE trgTrackToPlane (u,E,dl,a,b,c,d,ok,spect)
      IMPLICIT NONE
!      REAL*8    u(6),E,dl,a,b,c,d
      REAL*8    u(9),E,dl,a,b,c,d  !OR - 4/04
      INTEGER spect
      LOGICAL ok
* --  track a single particle with given start parameters
*     and find the intersection of the particle track with a given plane
*
*     Parameter:
*        u     IO : coordinate vector (initial/final)
*                     u0(1,2,3) : x, y, z [cm]
*                     u0(4,5,6) : dx/dt, dy/dt, dz/dt [cm/ns] 
*        E     I  : particle energy [MeV] * sign of particle charge
*                   (negative for electrons, positive for protons/deuterons)
*        dl    I  : step size [cm]
*        a..d  I  : parameter of the intersection plane 
*                   0 = a*x+b*y+c*z+d; 
*        ok    IO : status variable 
*                   - if false no action is taken 
*                   - set to false when no intersection point is found 
*                    
                 
      REAL*8   factor
      COMMON /trgConversionFactor/factor

!      REAL*8    ts,n,an,bn,cn,dn,maxdist,dist0,dist1,u0(6),u1(6)
      REAL*8    ts,n,an,bn,cn,dn,maxdist,dist0,dist1,u0(9),u1(9)  !  OR - 4/04
      REAL*8    one       
      INTEGER i,steps,max_steps
c
c      write(*,*) ' E = ',E,u
c For Bdl
        one=1.
        do i=7,9
        u(i)=0.0
        u0(i)=0.0
        u1(i)=0.0
        end do

      IF (.NOT. OK) RETURN   
c        write(*,*) ' u = ',u
c        write(*,*) ' track to plane a-d = ',a,b,c,d
        
      n  = 1.d00/SQRT (a*a+b*b+c*c)
c      n = -SIGN(one,d)*n
      an = a*n
      bn = b*n
      cn = c*n
      dn = d*n
      
    
      factor =  90.d00/E
      ts     = -dl/sqrt(u(4)**2+u(5)**2+u(6)**2)


      dist0   = u(1)*an + u(2)*bn + u(3)*cn + dn
      maxdist = max(ABS(dist0)*4.d00,1.0d00)
      
      ! check for the tracking direction 
!      CALL trgRK4(u,u1,ts,spect)
      CALL trgRK4Bdl(u,u1,ts,spect)
      dist1 = u1(1)*an + u1(2)*bn + u1(3)*cn + dn  
c      write(*,*) ' init =',dist0,dist1
      IF ((SIGN(one,dist0) .EQ. SIGN(one,dist1)) .AND.
     >    (ABS(dist0) .LT. ABS(dist1))) ts=-ts
c         write(*,*) ' ts = ',ts
         
      ! track through the intersection plane 
! GAW 99/11/22
! Previously, if dist1 and dist0 had different signs, it move the track one
! extra step so that interpolation in end is wrong.  The added if prevents that.
      steps = 0
      max_steps = int(max(dist0,10.*dl)/dl)*10
      IF (SIGN(one,dist0).eq.SIGN(one,dist1)) THEN
        dist1 = dist0   
        DO WHILE ((SIGN(one,dist0) .EQ. SIGN(one,dist1)) .AND. ok) 
!          CALL trgRK4(u1,u0,ts,spect)
        CALL trgRK4Bdl(u1,u0,ts,spect)
          dist0 = u0(1)*an + u0(2)*bn + u0(3)*cn + dn 
          IF (SIGN(one,dist0) .EQ. SIGN(one,dist1)) THEN
!            CALL trgRK4(u0,u1,ts,spect)
            CALL trgRK4Bdl(u0,u1,ts,spect)
            dist1 = u1(1)*an + u1(2)*bn + u1(3)*cn + dn  
          ENDIF
          ok = (ABS(dist1) .LT. maxdist).and.steps.lt.max_steps
c          write(*,*) dist0,dist1,ok
          steps = steps+1
        ENDDO        
      ELSE
        DO i=1,6
          u0(i) = u(i)
        ENDDO
      ENDIF
      
      IF (ok) THEN        
        ! calculate the intersection point
        DO i=1,6
          u(i) = u0(i) + (u1(i)-u0(i)) * dist0/(dist0-dist1)
        ENDDO

!Bdl

        do i=7,9
          u(i) = u0(i) + (u1(i)-u0(i)) * dist0/(dist0-dist1)
!u(i)=u0(i)
         end do

      ENDIF
           
      RETURN
      END
