c___________________________________________________________________________
c
c     SCEPTIC3D
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

c***********************************************************************
c General version allows choice of reinjection scheme.
c***********************************************************************
      subroutine reinject(i,dt,icolntype,bcr)

      integer bcr 

      if(bcr.eq.3) then
         call mcreinject(i)
      else
c      if(bcr.ne.0) then
         call maxreinject(i,dt,bcr)
c      elseif(icolntype.eq.1.or.icolntype.eq.5) then
c Injection from fv distribution at the boundary.
c         call fvreinject(i,dt,icolntype)
c      elseif(icolntype.eq.2.or.icolntype.eq.6)then
c Injection from a general gyrotropic distribution at infinity
c         call ogenreinject(i,dt)
c      else
c Injection from a shifted maxwellian at infinity
c         call oreinject(i,dt)
c      endif
      endif
      end
c***********************************************************************
      subroutine injinit(icolntype,bcr,colnwt)

      integer bcr

      if(bcr.eq.3) then
         call mcrinjinit(icolntype,colnwt)
      else
c      if(bcr.ne.0) then
c Injection from a maxwellian at boundary?
         call maxinjinit(bcr)
c      elseif(icolntype.eq.1.or.icolntype.eq.5) then
c Injection from fv distribution at the boundary.
c         call fvinjinit(icolntype)
c      elseif(icolntype.eq.2.or.icolntype.eq.6)then
c Injection from a general gyrotropic distribution at infinity
c         call ogeninjinit(icolntype)
c      else
c Injection from a shifted maxwellian at infinity
c         call oinjinit()
c      endif
      endif
      end
