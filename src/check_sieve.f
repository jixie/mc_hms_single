        logical function check_sieve(x,y)
        implicit none

* Original version made on 01/30/02 by M.E. Christy
* to more accurately model the size/shape of the HMS
* sieve slit (with much borrowed from check_dipole.f) ...

        include 'apertures.inc'
        real*8 x,y,xlocal,ylocal
        real*8 x_local,y_local,rad,rad1,rad2 
        integer i,j
        logical check0

        check_sieve=.false.
        check0 = .false.
        rad1 = 0.254/2.
        rad2 = 0.508/2.

* Let us observe first the obvious symmetry of the problem
* This helps reduce the checks to the first quadrant only...

        xlocal=abs(x)
        ylocal=abs(y)


* Now compare the current position and compare it with the different 
* apertures..
            
        if(.not.check_sieve) then
         do i=1,5
          x_local = xlocal - 2.540*(i-1)
          do j=1,4
           y_local = ylocal - 1.524*(j-1)
           rad = rad2 
           if(i.EQ.1.AND.j.EQ.1) rad = rad1
           if(i.EQ.2.AND.j.EQ.2.AND.x.GT.0.0.AND.y.GT.0.0) rad = 0.0
           if(i.EQ.3.AND.j.EQ.2.AND.x.LT.0.0.AND.y.LT.0.0) rad = 0.0
           
           check0 = ((x_local*x_local + y_local*y_local).LE.(rad*rad))
           if(check0) check_sieve=.true.
          enddo
         enddo
        endif

        return
        end




