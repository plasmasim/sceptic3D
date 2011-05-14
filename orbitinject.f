
c***********************************************************************
c General version allows choice of reinjection scheme.
c***********************************************************************
      subroutine reinject(i,dt,icolntype)

      integer bcr 

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
      end
c***********************************************************************
      subroutine injinit(icolntype,bcr)

      integer bcr

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
      end
