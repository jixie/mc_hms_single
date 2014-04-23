          subroutine transp(spectr,class,decay_flag,dflag,m2,ppi,zd)
C+______________________________________________________________________________
!
! TRANSP - This subroutine transports a particle through the various
!  segments of a spectrometer. The passed variable CLASS determines
!   which transformation to use. Each transformation IS SEQUENTIAL,
!   and carries the particle from the last transformation to a particular
!   plane Z=const in the spectrometer.
!
!   NOTE: The coordinate system used here is the right handed "TRANSPORT"
!   coordinate system, in which +Z -> downstream, +X -> in median plane,
!   pointing in direction taken through bending magnet by high momentum rays,
!   +Y = transverse direction.
!
!   For an upward vertical bend spectrometer:
!
!   +X points towards the floor,
!   +Y points horizontally LEFT as one looks downstream.
!
! D. Potterveld - March 1993.
!
!   Modification History:
!
! August 1993.	(D. Potterveld) Modified to begin each transformation at
!  the pivot.
!
!  20-AUG-1993	(D. Potterveld) Modified to use COSY INFINITY forward maps.
!
!   03-MAR-1994	(D. Potterveld) Fixed bug to correctly compute cosy "A" and "B".
!
!  10-MAY-1995	(D. Potterveld) Switched to sequential transformations.
!
!  11-MAY-1995	(D. Potterveld) Added rejection of extremely small coeff.
!
!  30-OCT-1995	(D. Potterveld) Change to COSY-7 transport units.
!
!  16-SEP-1998  Check for decay of particle. dflag is true of the particle
!  has already decayed, so check for decay if dflag .eq. .false.
!  After pathlength is calculated, add extra check for decay
!  within the different path length.
!
!  13-NOV-2001  (E. Christy)  put in hack to allow field adjustment for magnets.
!               these are passed in hms_magnets.cmn
!
C-______________________________________________________________________________

        implicit none

        include 'track.inc'
        include 'hms_magnets.cmn'

C Arguments.

        integer*4 spectr  !HMS=1, SOS=2
        integer*4 class
        logical decay_flag !check for decay
        logical dflag !has particle decayed yet?
        real*8  m2,ppi,zd !decay variables.

C Parameters.

        integer*4 max_class
        parameter (max_class = 18) !maximum number of classes.
        real*8 coeff_min
        parameter(coeff_min = 1.0d-14) !coeff's smaller than this = 0

C Local declarations.

        integer*4 idummy

        real*8 delta_z !pathlength difference.
        integer*4 i,j,k,kk,chan,n_classes

        character*132 str_line
        character*132 file_name

C Cosy reconstruction matrix elements.

        integer*4 max_elements
        parameter (max_elements = 500)
        integer*4 nspectr
        parameter (nspectr=2)
        real*8  coeff(nspectr,5,max_elements,max_class)
        integer*2 expon(nspectr,5,max_elements,max_class)
        integer*4 n_terms(nspectr,max_class)
        real*8 sum(5),ray(5),term,temp



C Function definitions.

! integer*4 last_char
          logical*4 locforunt
!dec real*8 grnd

C No amnesia between calls!!!

          save

C ================================ Executable Code =============================

C A word from the sponsor:
C  a) Using sequential matrix elements, variables with the ending 's' should be
C     used, otherwise they should end on '_transp'
C  b) the sequential matrix elements have been made with COSY 7, i.e. the units
C      to be used are cm, mrad, not m, slopes

CCCCC    Eric's optics hack.  The correction back to the nominal     CCCCCC
CCCCC    are done at the end of this routine.                        CCCCCC


         if(class.EQ.2) dpps = dpps-q1cor   ! Q1 in   !
         if(class.EQ.5) dpps = dpps-q2cor    ! Q2 in   !
         if(class.EQ.8) dpps = dpps-q3cor    ! Q2 in   !
         if(class.EQ.11) dpps = dpps-dcor    ! D in    !

C Pack local copy of input coordinates.

        if(spectr.eq.1)then
         ray(1) = xs  !cm. ( "X" )
         ray(2) = dxdzs*1000. !mrad. ( "THETA" )
         ray(3) = ys !cm. ( "Y" )
         ray(4) = dydzs*1000. !mrad. ( "PHI" )
         ray(5) = dpps !Fractional "Delta P/P"
        else !SOS - uses sequential transformations.
          ray(1) = xs !cm.  ( "X" )
          ray(2) = dxdzs*1000. !mrad. ( "THETA" )
          ray(3) = ys !cm.  ( "Y" )
          ray(4) = dydzs*1000. !mrad. ( "PHI" )
          ray(5) = dpps !Fractional variation.
        endif 

C Reset COSY sums.

        do i = 1,5
          sum(i) = 0.
        enddo

C Compute COSY sums.

        k = class
        do i = 1,n_terms(spectr,k)
          term = 1.0
          do j = 1,5
            temp = 1.0
           if (expon(spectr,j,i,k).ne.0.) temp = ray(j)**expon
     &  (spectr,j,i,k) 
          term = term*temp
          enddo
          sum(1) = sum(1) + term*coeff(spectr,1,i,k) ! NEW "X"
          sum(2) = sum(2) + term*coeff(spectr,2,i,k) ! NEW "A"
          sum(3) = sum(3) + term*coeff(spectr,3,i,k) ! NEW "Y"
          sum(4) = sum(4) + term*coeff(spectr,4,i,k) ! NEW "B"
          sum(5) = sum(5) + term*coeff(spectr,5,i,k) ! NEW "dL"
        enddo
     
C Unpack output coordinates. Note that DPPS is unchanged by transformation.

        if(spectr.eq.1)then
          xs    = sum(1) !cm
          dxdzs = sum(2)/1000. !slope (mr)
          ys    = sum(3) !cm
          dydzs = sum(4)/1000. !slope (mr)
           delta_z = -sum(5) !deltaZ (cm)
        else
          xs    = sum(1) !cm
          dxdzs = sum(2)/1000. !slope (mr)
          ys    = sum(3) !cm
          dydzs = sum(4)/1000. !slope (mr)
          delta_z = -sum(5) !deltaZ (cm)
        endif


CCCCCC         Correct Magnets back to nominal            CCCCCC

         if(class.EQ.4) dpps = dpps+q1cor    ! Q1 out  !
         if(class.EQ.7) dpps = dpps+q2cor    ! Q2 out  !
         if(class.EQ.10) dpps = dpps+q3cor   ! Q2 out  !
         if(class.EQ.12) dpps = dpps+dcor    ! D out   !

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


         return

C###############################################################################

C Initialization entry points.

          entry transp_init(spectr,n_classes)

C Use default filename.
C    HMS: if spectr = 1
C    SOS: if spectr = 2

         if (spectr.eq.1) then
!!   write(6,*) 'TRANSP_INIT HMS!!!!!!'
          file_name='forward_cosy.dat'
         else
!!     write(6,*) 'TRANSP_INIT SOS!!!!!!'
          file_name='forward_cosy.dat'
         endif

C Open input file.

        if (.not.locforunt(chan))
     >  stop 'TRANSP_INIT: No I/O channels!'
        open(unit=chan,status='old',file=file_name)

C Strip away header.

        str_line = '!'
        do while (str_line(1:1).eq.'!')
        read (chan,1001) str_line
!!     if (str_line(1:1).eq.'!') write(6,*) str_line(1:min(79,last_char(str_line)))
        enddo

C Read in the transformation tables.

        n_classes = 0
        do i = 1,max_class
          n_terms(spectr,i) = 0
        enddo

        do while (.true.)
          kk = n_classes + 1

! If too many transformations, complain!

          if (kk.gt.max_class) stop 'TRANSP_INIT: too many 
     &  transformations!'

! Add data lines to table, looking for flag line.

          do while (str_line(1:4).ne.' ---')
           n_terms(spectr,kk) = n_terms(spectr,kk) + 1
            if (n_terms(spectr,kk).gt.max_elements)
     >      stop 'TRANSP_INIT: too many COSY terms!'

            read (str_line,1200)
     >(coeff(spectr,i,n_terms(spectr,kk),kk),i=1,5),
     >(expon(spectr,j,n_terms(spectr,kk),kk),j=1,4),idummy,
     > expon(spectr,5,n_terms(spectr,kk),kk)

! Ignore time-of-flight term.

           if (idummy.ne.0) then
           if (coeff(spectr,1,n_terms(spectr,kk),kk).ne.0.or.
     >coeff(spectr,2,n_terms(spectr,kk),kk).ne.0.or.
     >coeff(spectr,3,n_terms(spectr,kk),kk).ne.0.or.
     >coeff(spectr,4,n_terms(spectr,kk),kk).ne.0)
     >stop 'TRANSP_INIT: non-zero TOF terms!'
        n_terms(spectr,kk) = n_terms(spectr,kk) - 1
            endif

! Fetch next line from file.

            read (chan,1001) str_line
         enddo

! If flag line is seen, increment transformation counter.

          n_classes = kk
!!   write(6,*) 'TRANS, ORDER, TERMS =',kk,
!!     >       expon(spectr,5,n_terms(spectr,kk),kk),n_terms(spectr,kk)

! Read lines until a non-blank, non-comment non-terminal line is found.

150       read (chan,1001,end=200) str_line
         if (str_line(1:1).eq.'!'.or.str_line(1:4).eq.' ---'.or.
     >    str_line.eq.'    ') goto 150
       enddo

C Done with file.

200     close (unit=chan)

C Go home.

        return

C ============================== Format Statements =============================
1001    format(a)
1200    format(1x,5g14.7,1x,6i1)

         end









