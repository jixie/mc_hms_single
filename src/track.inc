C Common block with track variables.

        real*8 xs,ys,zs  !track x,y,z positions (cm)
        real*8 dxdzs,dydzs !track slop (dx/dz, dy/dz)
        real*8 dpps !delta dp/p (%).

        common  /track/ xs,ys,zs,dpps,dxdzs,dydzs

C xs,ys,zs,dxdzs,dydzs are positions (cm) and slopes (dx/dz and dy/dz) for
C the track.  Initially, you start with the pivot values.  As you go thru
C the magnets, these are the values output by transp.f.
