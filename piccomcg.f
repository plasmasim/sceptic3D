
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
