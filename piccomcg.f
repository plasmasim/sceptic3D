c
c SCEPTIC3D
c___________________________________________________________________________
c
c     This code is copyright (c)
c              Ian H Hutchinson    hutch@psfc.mit.edu.
c              Leonardo Patacchini patacchi@mit.edu
c
c     It may be used freely with the stipulation that any scientific or
c     scholarly publication concerning work that uses the code must give
c     an acknowledgement referring to the relevant papers
c
c     I.H. Hutchinson, Plasma Physics and Controlled Fusion, vol 44, p
c     1953 (2002), vol 45, p 1477 (2003).
c
c     L. Patacchini and I.H. Hutchinson, Plasma Physics and Controlled
c     Fusion, vol 49, p1193 (2007), vol 49, p 1719 (2007).
c
c     I.H. Hutchinson and L. Patacchini, Physics of Plasmas, vol 14,
c     p013505 (2007)
c
c     The code may not be redistributed except in its original package.
c
c     No warranty, explicit or implied, is given. If you choose to build
c     or run the code, you do so at your own risk.
c___________________________________________________________________________

c Here are stored the common declarations for the parallel bloc solver

      
c     Number of participating processors
      integer nprocs


c     Cartesian communicator
      integer icommcart

      parameter(ndims=3)
c Declared Dimensional structure of u must be (Li,Lj,Lk,..) passed using
c iLs, which embodies pointer steps. For 2d iLs=(1,Li,Li*Lj), 
c 3d (1,Li,Li*Lj,Li*Lj*Lk) etc 
      integer iLs(ndims+1)
c iuds used dimensions of u
      integer iuds(ndims)

c     dimension of each topology (nb of blocks)
      integer idims(ndims)

c     cartesian topology coords of this process (block)
      integer icoords(ndims)

c     structure of icoords (1,(idims(1)+1),(idims(1)+1)*(idims(2)+1),...)
      integer iLcoords(ndims+1)

c     origin of the blocks (in the total i,j frame), in a linear referencing
      parameter (norigmax=1000)
      integer iorig(norigmax)

c     dimension of my block (in the total i,j frame)
      integer myside(ndims)

c     myorig1 is just the origin of the 1st dimension
      integer myorig,myorig1,myorig2,myorig3

      integer ifull(ndims)

c     If we are on the outer or inner boundary, for the BC
      logical out,inn

c     cartid
      integer mycartid
      
c      common /sor_parallel/ nprocs,icommcart,iLs,iuds,idims,icoords
c     $     ,iorig,ifull,myside,mycartid,myorig,myorig1,myorig2
      common/cg_parallel/iLcoords,icommcart,mycartid,ifull,iorig,iuds
     $     ,idims,iLs,icoords,myside,myorig,myorig1,myorig2,myorig3,out
     $     ,inn
