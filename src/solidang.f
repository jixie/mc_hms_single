      program solidang

      integer*4 i,k
      real*4 domegav(20),accept(16,20),temp(16), de

      open(unit=10,file="taccept.dat.bak",status='old')

      de = 3.077*.01

      do i=1,20
       read(10,*) domegav(i)
      enddo

c      write(6,2000)domegav

      do i=1,20
       read(10,*) temp
       do j=1,16 
        accept(j,i) = temp(j)*domegav(i)
c        write(6,*) accept(j,i)
       enddo
      enddo

      do j=1,20
       do i=1,16
        temp(i) = accept(i,j)*de
c        write(6,*) accept(i,j)
       enddo
       write(6,2000) temp
      enddo

 2000 format(16f7.4)


      end
