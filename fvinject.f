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
c Currently this is set up for a constant nu charge-exchange distrib.
c***********************************************************************
c Reinjection from a general ion distribution function.
      subroutine fvreinject(i,dt,icolntype)
      include 'piccom.f'
      include 'fvcom.f'
      include 'colncom.f'
      real qz(nzfvi:nzfva),vzp(nzfvi:nzfva)
      logical lpos, lpos1
      logical istrapped
      real eps
c diagnostics
      real vydist(nxfvi:nxfva),vzdist(nxfvi:nxfva)
      common /vinjdiag/vydist,vzdist

c A really small number
      eps=1.e-20
c___________________________________________________________________
c Pick a random th: pth
      idum=1
 1    y1=ran0(idum)*qthfv(nthfvsize)
      
      call f1invtfunc(qthfv,nthfvsize,y1,pth)
      
      ipth=pth
      fpth=pth-ipth
      costheta=(1.-fpth)*fvth(ipth)+ fpth*fvth(ipth+1)
c      write(*,'(a,2f10.5,i4,2f10.5)')'y,ymax,ipth,fpth,cos ',
c     $     y1,qthfv(nthfvsize),ipth,fpth,costheta

c___________________________________________________________________
c Pick a random vx: pxfv
 2    y2=ran0(idum)*((1.-fpth)*qxfv(nxfva,ipth)+fpth*qxfv(nxfva,ipth+1))
      call f2invtfunc(qxfv(nxfvi,ipth),qxfv(nxfvi,ipth+1)
     $     ,nxfva-nxfvi+1,y2,pxfv,(1.-fpth),fpth)
c      write(*,*) pxfv
      ixfv=pxfv
      fxfv=pxfv-ixfv
      ixfv=ixfv-1+nxfvi
c      write(*,*)(qxfv(j,ipth),j=nxfvi,nxfva)
c      write(*,*)(qxfv(j,ipth+1),j=nxfvi,nxfva)
      vx=(1.-fxfv)*vxfv(ixfv) + fxfv*vxfv(ixfv+1)
c      write(*,'(a,2f10.5,i4,2f10.5)')'y,ymax,ixfv,fxfv,vx  ',y2,
c     $     ((1.-fpth)*qxfv(nxfva,ipth)+fpth*qxfv(nxfva,ipth+1))
c     $     ,ixfv,fxfv,vx

c      write(*,*)(qfv(j,ixfv,ipth+1),j=nxfvi,nxfva)

c___________________________________________________________________
c Calculate the exact vztr: vz transition (where n.v=0),
c for this vx, costheta. If it is too extreme, try again.
      if(costheta.eq.0.)then
c         vztr=-1.e30
         goto 1
      else
         vztr=-vx*sqrt(1.-costheta**2)/costheta
      endif
      if(costheta.gt.0.)then
         if(vztr.gt.vzfv(nzfva))then
c            write(*,*)'Impossible vztr',vztr,costheta,vx,' trying again'
            goto 2
         elseif(vztr.lt.vzfv(nzfvi))then
            vztr=vzfv(nzfvi)+1.e-4
         endif
      else
         if(vztr.lt.vzfv(nzfvi))then
c            write(*,*)'Impossible vztr',vztr,costheta,vx,' trying again'
            goto 2
         elseif(vztr.gt.vzfv(nzfva))then
            vztr=vzfv(nzfva)-1.e-4
         endif         
      endif

c___________________________________________________________________
c Pick a random vz: pzfv
      s=(1.-fpth)*(1.-fxfv)
      t=(1.-fpth)*(fxfv)
      u=(fpth)*(1.-fxfv)
      v=(fpth)*(fxfv)
c     Total qz includes only positive or negative sections, depending on
c     whether th is positive or negative. The transition between positive
c     and negative occurs at a z-index ztrfv(vx). So if th is +ve, all
c     qfv values at ztr and below are zero, while if th is -ve all values
c     above ztr are zero (and the order is reversed).
      lpos=(fvth(ipth).ge.0.)
      lpos1=(fvth(ipth+1).ge.0.)
c zero crossing indices for the 4 interpolation points
      ztrs=ztrfv(ixfv,ipth)
      ztrt=ztrfv(ixfv+1,ipth)
      ztru=ztrfv(ixfv,ipth+1)
      ztrv=ztrfv(ixfv+1,ipth+1)
      xfv=float(ixfv)
c Zero the qz array
      do jz=nzfvi,nzfva
         qz(jz)=0.
      enddo
c Create an average qz that that interpolates between the 4 different
c cases of th and vx.
      do jz=nzfvi,nzfva
         zfv=jz
c Never allow non-zero cumulative distrib, for negative n.v.
c         if(zfv.gt.ztr)then
            if(lpos) then
               if(zfv.gt.ztrs) qz(jz)=qz(jz)+
     $              s*qfv(jz,ixfv,ipth)
               if(zfv.gt.ztrt) qz(jz)=qz(jz)+
     $              t*qfv(jz,ixfv+1,ipth)
            else
               if(zfv.le.ztrs) qz(nzfva-jz+nzfvi)=qz(nzfva-jz+nzfvi)+
     $              s*qfv(jz,ixfv,ipth)
               if(zfv.le.ztrt) qz(nzfva-jz+nzfvi)=qz(nzfva-jz+nzfvi)+
     $              t*qfv(jz,ixfv+1,ipth)
            endif
            if(lpos1) then
               if(zfv.gt.ztru) qz(jz)=qz(jz)+
     $              u*qfv(jz,ixfv,ipth+1)
               if(zfv.gt.ztrv) qz(jz)=qz(jz)+
     $              v*qfv(jz,ixfv+1,ipth+1)
            else
               if(zfv.le.ztru) qz(nzfva-jz+nzfvi)=qz(nzfva-jz+nzfvi)+
     $              u*qfv(jz,ixfv,ipth+1)
               if(zfv.le.ztrv) qz(nzfva-jz+nzfvi)=qz(nzfva-jz+nzfvi)+
     $              v*qfv(jz,ixfv+1,ipth+1)
            endif
c         endif
      enddo
      if(ldiaginj)then
         write(*,*)'stuv lpos lpos1',s,t,u,v,lpos,lpos1
         write(*,*)'ztrs-v=',ztrs,ztrt,ztru,ztrv
         write(*,*)'qfv for',ixfv,ipth
         write(*,'(8g10.3)') (qfv(j,ixfv,ipth)
     $        ,j=nzfvi,nzfva)
         write(*,*)'qfv for',ixfv+1,ipth
         write(*,'(8g10.3)') (qfv(j,ixfv+1,ipth),j=nzfvi,nzfva)
         write(*,*)'qfv for',ixfv,ipth+1
         write(*,'(8g10.3)') (qfv(j,ixfv,ipth+1),j=nzfvi,nzfva)
         write(*,*)'qfv for',ixfv+1,ipth+1
         write(*,'(8g10.3)') (qfv(j,ixfv+1,ipth+1),j=nzfvi,nzfva)
         write(*,*)'qz='
         write(*,'(8g10.3)') qz
c        write(*,'(8e10.3)') (qfv(j,ixfv,ipth),j=nzfvi,nzfva)
      endif
      y3=ran0(idum)*qz(nzfva)
      call f1invtfunc(qz,nzfva-nzfvi+1,y3,pzfv)

      izfv=pzfv
      fzfv=pzfv-izfv
c      write(*,*)'izfv,fzfv=',izfv,fzfv
c ----------------------------------------------
c Here we might to correct for non-linearity at the transition.
c As follows:
c              MAY BE WORKING.
c Solve to get exact crossing z-index ztrp, relative to initial index 1.
      call f1invtfunc(vzfv,nzfva-nzfvi+1,vztr,ztrp)
c If it does not solve which should never happen, then tell us 
      if(ztrp.eq.0.)then
         write(*,*)'Failure of f1invtfunc for vztr=',vztr
         stop
      endif
c If we have reversed the qfv, ie cos negative, then invert this ztrp.
      if(costheta.lt.0.)then
         ztrp=nzfva-nzfvi+2-ztrp
c         if(ztrp.eq.float(nzfva-nzfvi+1))ztrp=ztrp-1.e-4
      endif
c Solution is shifted only if it would give a negative projection.
c This is not exactly the same as linear interpolation and may not be best.
      if(pzfv.lt.ztrp)then 
c Shift the solution to correspond to ramp to first non-zero point,
c from ztrp, rather than from izfv as the zero.
         if(ldiaginj) write(*,*)'***** ztr correction on',pzfv,ztrp,vztr
         pzfv=fzfv*(int(ztrp)+1)+(1.-fzfv)*(ztrp)
         izfv=pzfv
         fzfv=pzfv-izfv
         if(ldiaginj) write(*,*)'changed to',pzfv,fzfv,vzrev
      endif
c End of non-linearity correction
c-----------------------------------------------
c Shift to nzfvi-based.
      izfv=izfv+nzfvi-1
c Translate back to positive and negative from all positive indices.
      if(costheta.lt.0.)then
         izfv=-izfv-1
         fzfv=1.-fzfv
c Plot backwards
         if(ldiaginj)then
            do j=nzfvi,nzfva
               vzp(j)=vzfv(nzfva+nzfvi-j)
            enddo
         endif
      else
         if(ldiaginj)then
            do j=nzfvi,nzfva
               vzp(j)=vzfv(j)
            enddo
         endif
      endif
      vz=(1.-fzfv)*vzfv(izfv) + fzfv*vzfv(izfv+1)
c Diagnostics assuming uniform vz arrays.
      izp=nint(nzfvi+(nzfva-nzfvi)*(vz-vzfv(nzfvi))
     $     /(vzfv(nzfva)-vzfv(nzfvi)))
c Obsolete     izp=nint(vz*nzfva/vzfv(nzfva))
      if(izp.lt.nzfvi) iyp=nzfvi
      if(izp.gt.nzfva) iyp=nzfva
      vzdist(izp)=vzdist(izp)+1

      sintheta=sqrt(1.-costheta**2)
c----------------------------------------------------------------------
c This should never be negative.
      vproj=costheta*vz+sintheta*vx
      if(vproj.lt.0.)then
         write(*,*) 'Negative projection'
c     $        , costheta,sintheta,vz,vx,
c     $      costheta*vz+sintheta*vx
         write(*,'(a,2f10.6,i4,2f10.6)')'y2,ymax,ixfv,fxfv,vx  ',y2,
     $     ((1.-fpth)*qxfv(nxfva,ipth)+fpth*qxfv(nxfva,ipth+1))
     $     ,ixfv,fxfv,vx
         write(*,'(a,3f10.6,i4,2f10.6)')'y3,ymax,izfv,fzfv,vz  '
     $        ,y3,qz(nzfva),y3/qz(nzfva),izfv,fzfv,vz
         write(*,*)'i,i+1,qz(i),qz(i+1)',izfv,izfv+1,
     $        qz(izfv),qz(izfv+1)
c Don't try to plot here. Just try again.
         goto 2
         ldiaginj=.true.
      endif
c End of solving for the reinjection th,vx,vz.
c___________________________________________________________________
      if(ldiaginj)then
c     Diagnostic plotting.
         write(*,'(a,4f10.4)')'Reinject costheta,vx,vz,v-projection='
     $        ,costheta,vx,vz,costheta*vz+sqrt(1.-costheta**2)*vx

         call autoplot(fvth(1),qthfv,nthfvsize)
         call polymark(fvth(1),qthfv,nthfvsize,1)
         call axlabels('cos(theta)','qthfv')
         call color(4)
         call vecw(-1.,y1,0)
         call vecw(1.,y1,1)
         call vecw(costheta,0.,0)
         call vecw(costheta,yn2yw(0.7),1)
         call color(15)
         call pltend()

         call autoplot(vxfv,qxfv(nxfvi,ipth),nxfva-nxfvi+1)
         call polymark(vxfv,qxfv(nxfvi,ipth),nxfva-nxfvi+1,1)
         call color(2)
         call polyline(vxfv,qxfv(nxfvi,ipth+1),nxfva-nxfvi+1)
         call polymark(vxfv,qxfv(nxfvi,ipth+1),nxfva-nxfvi+1,2)
         call axlabels('vx','qxfv')
         call color(4)
         call vecw(vxfv(nxfvi),y2,0)
         call vecw(vxfv(nxfva),y2,1)
         call vecw(vx,0.,0)
         call vecw(vx,yn2yw(0.7),1)
         call color(15)
         call pltend()

         call autoplot(vzfv,qfv(nzfvi,ixfv,ipth),nzfva-nzfvi+1)
         y4=0.3*yn2yw(0.7)
         write(*,*)'ztr interps',ztrs-nzfvi+1,ztrt-nzfvi+1,
     $        ztru-nzfvi+1,ztrv-nzfvi+1,ztrp,pzfv
         xs=finterp(vzfv,ztrs-nzfvi+1,nzfva-nzfvi+1)
         xt=finterp(vzfv,ztrt-nzfvi+1,nzfva-nzfvi+1)
         xu=finterp(vzfv,ztru-nzfvi+1,nzfva-nzfvi+1)
         xv=finterp(vzfv,ztrv-nzfvi+1,nzfva-nzfvi+1)
         x1=finterp(vzp,ztrp,nzfva-nzfvi+1)
         x2=finterp(vzp,pzfv,nzfva-nzfvi+1)
         call axlabels('vz','qfv, qz')
         call polymark(vzfv,qfv(nzfvi,ixfv,ipth),nzfva-nzfvi+1,1)
         call vecw(xs,0.,0)
         call vecw(xs,y4,1)
         call color(2)
         call polyline(vzfv,qfv(nzfvi,ixfv,ipth+1),nzfva-nzfvi+1)
         call polymark(vzfv,qfv(nzfvi,ixfv,ipth+1),nzfva-nzfvi+1,2)
         call vecw(xt,0.,0)
         call vecw(xt,y4,1)
         call color(3)
         call polyline(vzfv,qfv(nzfvi,ixfv+1,ipth),nzfva-nzfvi+1)
         call polymark(vzfv,qfv(nzfvi,ixfv+1,ipth),nzfva-nzfvi+1,3)
         call vecw(xv,0.,0)
         call vecw(xv,y4,1)
         call color(4)
         call polyline(vzfv,qfv(nzfvi,ixfv+1,ipth+1),nzfva-nzfvi+1)
         call polymark(vzfv,qfv(nzfvi,ixfv+1,ipth+1),nzfva-nzfvi+1,4)
         call vecw(xv,0.,0)
         call vecw(xv,y4,1)
         call color(5)
         call polyline(vzp,qz(nzfvi),nzfva-nzfvi+1)
         call polymark(vzp,qz(nzfvi),nzfva-nzfvi+1,5)
         call drcstr('qz')
         call vecw(vztr,0.,0)
         call vecw(vztr,y4,1)
         call color(8)
         call drcstr('vztr')
         write(*,*)'x1=',x1,' ztrp=',ztrp,' x2=',x2,' pzfv=',pzfv
         call vecw(x1,0.,0)
         call vecw(x1,y4,0)
         call color(6)
         call vecw(vzfv(nzfvi),y3,0)
         call vecw(vzfv(nzfva),y3,1)
         call vecw(vz,0.,0)
         call vecw(vz,yn2yw(0.7),1)
         call color(15)
         call pltend()
      endif
      if(vproj.lt.0.)then 
         ldiaginj=.false.
         goto 2
      endif
c_______________________________________________________________________
c Now pick an azimuthal angle of position. 
c The x-direction points at angle phiazim relative to cartesian 1-direc.
c The z-direction and 3-direction coincide.
c The x-z plane is the plane containing the surface normal.
      phiazim=pi*2.*ran0(idum)
      cosphi=cos(phiazim)
      sinphi=sin(phiazim)
c_______________________________________________________________________
c Pick the ignorable velocity in the y-direction. fqvxvz is normalized
c Pick a random vx: pxfv
      y=ran0(idum)
      call f2invtfunc(fqvxvz(nxfvi,izfv),fqvxvz(nxfvi,izfv+1)
     $     ,nxfva-nxfvi+1,y,pyfv,(1.-fzfv),fzfv)

c      write(*,*) pyfv,y,fzfv
      iyfv=pyfv
      fyfv=pyfv-iyfv
      iyfv=iyfv-1+nxfvi
      vy=(1.-fyfv)*vxfv(iyfv) + fyfv*vxfv(iyfv+1)
c Diagnostics
      iyp=nint(vy*nxfva/vxfv(nxfva))
      if(iyp.lt.nxfvi) iyp=nxfvi
      if(iyp.gt.nxfva) iyp=nxfva
      vydist(iyp)=vydist(iyp)+1
c___________________________________________________________________
c Install the reinjection velocity:
c Up till now, all velocities are in units of sqrt(2Ti/m_i).
c Convert to units of sqrt(T_e/m_i) using Ti value in units of Te.
      vscale=sqrt(2*Ti)
      xp(6,i)=vscale*vz
      xp(4,i)=vscale*(cosphi*vx-sinphi*vy)
      xp(5,i)=vscale*(sinphi*vx+cosphi*vy)
c The reinjection position. The surface normal for flux across
c surface element is the angle corresponding to costheta,
c which is therefore the inward normal.
      rs=-r(nr)*0.99999
      xp(3,i)=rs*costheta
      xp(2,i)=(rs*sintheta)*sinphi
      xp(1,i)=(rs*sintheta)*cosphi
c
c Obtain angle coordinate and map back to th for phihere.

      ct=-costheta

      call invtfunc(th(1),nth,ct,x)


      ic1=x
      ic2=ic1+1
      dc=x-ic1
c This expression should work for CIC And NGP.
      phihere=(phi(NRUSED,ic1)+phi(NRFULL,ic1))*0.5*(1.-dc)
     $        +(phi(NRUSED,ic2)+phi(NRFULL,ic2))*0.5*dc

c No longer do this for the new reinjection handling.      
c Increment the position by a random amount of the velocity.
c This is equivalent to the particle having started at an appropriately
c random position prior to reentering the domain.
c      xinc=ran0(idum)*dt
c      do j=1,3
c         xp(j,i)=xp(j,i)+xp(j+3,i)*xinc
c      enddo
c If the third bit (4) of icolntype is set this means we must add the 
c Eneutral acceleration for the effective prior step. 
c      if(mod(icolntype/4,2).eq.1)xp(6,i)=xp(6,i)+Eneutral*dt
c Do the outer flux accumulation.
      spotrein=spotrein+phihere
      nrein=nrein+1
c Reject particles that are already outside the mesh.
c With new reinjection, this should never happen, but did.
c Remove when satisfied:
c...........
      vp=xp(4,i)**2+xp(5,i)**2+xp(6,i)**2
      rp=xp(1,i)**2+xp(2,i)**2+xp(3,i)**2
      if(.not.rp.le.r(nr)*r(nr).or. rp.le.1.
     $     .or. .not.vp.lt.1.e8)then
         write(*,*)
         write(*,*)'Launch Error',nrein,sqrt(rp),xp(1,i),xp(2,i),xp(3,i)
         write(*,*)'velocity:',xp(4,i),xp(5,i),xp(6,i)
         write(*,*)'cosphi,sinphi,vx,vy,vz',cosphi,sinphi,vx,vy,vz
         write(*,*)'y,y3,iyfv,pyfv,fyfv,fzfv',y,y3,iyfv,pyfv,fyfv,fzfv
c        stop
c trying counting only once.
         nrein=nrein-1
c...........
         goto 1
      else
c Needs th to be initialized, I think.
         if(th(nth).eq.-1.)then
c The th array is initialized.
            if(istrapped(i))then
               ntrapre=ntrapre+1
c     v=sqrt(xp(4,i)**2+xp(5,i)**2+xp(6,i)**2)
c     write(*,*)'Trapped',vdx/rp,u,v,sqrt(u**2-2.*averein)
c     crt,czt,ceta,cosal
            endif
         endif
      endif

      end
c***********************************************************************
c***********************************************************************
      subroutine fvinjinit(icolntype)
c Initialize the needed data arrays.
      include 'piccom.f'
      include 'fvcom.f'
c Passing the drift velocity to fv.
      common /distfunc/ud
      parameter (vtrange=4.,udrange=8.)

c Velocity in this routine is normalized to a nominal ion thermal velocity
c which for a Maxwellian-related form is sqrt(2T_i/m).
      ud=vd/sqrt(2.*Ti)

c Decide what the positive and negative velocity ranges are and fit them
c to the mesh.
      vxfvi=-vtrange
      vxfva=vtrange
      vzfvi=min(-vtrange,-vtrange+udrange*ud)
      vzfva=max(vtrange,vtrange+udrange*ud)
c      write(*,*)'vzfvi,vzfva=',vzfvi,vzfva
      do i=0,nxfva
         vxfv(i)=vxfva*i/float(nxfva)
      enddo
      do i=-1,nxfvi,-1
         vxfv(i)=vxfvi*i/float(nxfvi)
      enddo
c Equal numbers of points for negative and positive. Inefficient.
c      do i=0,nzfva
c         vzfv(i)=vzfva*i/float(nzfva)
c      enddo
c      do i=-1,nzfvi,-1
c         vzfv(i)=vzfvi*i/float(nzfvi)
c      enddo
c Uniform vz-grid.
      do i=nzfvi,nzfva
         vzfv(i)=vzfvi+(vzfva-vzfvi)*(i-nzfvi)
     $        /(float(nzfva-nzfvi))
      enddo
c
      do i=1,nthfvsize
c theta array including poles
c Uniform in cos theta
         fvth(i)=1.-2.*(i-1)/(nthfvsize-1)
      enddo
c
      nthfv=(nthfvsize+1)/2
c Calculate qxfv:
      do ith=1,nthfvsize
c      do ith=1,nthfv
c For all theta positive?
         call calcqxfv(ith)
      enddo
c
c Integrate to give qthfv (th runs from 1 to -1 so we change sign).
      qthfv(1)=0.
      infty=nxfva
      do ith=2,nthfvsize
         qthfv(ith)=qthfv(ith-1) - (fvth(ith)-fvth(ith-1))*
     $        0.5*(qxfv(infty,ith)+qxfv(infty,ith-1))
      enddo
c The final result should be 1/sqrt(\pi)=.5641 for a Maxwellian. 
c This might be a good check.
      if(ud.eq.0. .and. myid.eq.0)
     $     write(*,*)'qthfv(nthfvsize)=',qthfv(nthfvsize),
     $     ' c.f. 0.5641 at ud=0'

c Calculate fqvxvz, the perpendicular cumulative distibution.
      call calcfqx()

      call srand(myid+1)
      end
c******************************************************************
c Calculate cumulative flux in direction given by 
c normalized direction vector (nx,nz): nz=th, nx=sqrt(1+n^2)
      subroutine calcqxfv(ith)
      include 'piccom.f'
      include 'fvcom.f'
      
c Cycle over vx calculating the rows.
      do ix=nxfvi,nxfva
         vx=vxfv(ix)
         call zintqfv(nzfvi,nzfva,vzfv(nzfvi),fvth(ith),vx,
     $        qfv(nzfvi,ix,ith),ztrfv(ix,ith))
c         ztrfv(ix,nthfvsize-ith+1)=ztrfv(ix,ith)
      enddo

      if(fvth(ith).ge.0.)then
c   Integrate to get qxfv only for positive angles.
c   Infinity values for p and m.
         infty =nzfva
         inftym=nzfvi
         qxfv(nxfvi,ith)=0.
         ithm=nthfvsize-ith+1
         qxfv(nxfvi,ithm)=0.
         do ix=nxfvi+1,nxfva
            qxfv(ix,ith)=qxfv(ix-1,ith)+ (vxfv(ix)-vxfv(ix-1))*
     $           0.5*(qfv(infty,ix,ith)+qfv(infty,ix-1,ith))
            ixm=(nxfva-ix+nxfvi)
c            if(ith.ne.ithm) Not sure about this?
            qxfv(ix,ithm)=qxfv(ix-1,ithm) -(vxfv(ixm)-vxfv(ixm+1))*
     $           0.5*(qfv(inftym,ixm,ith)+qfv(inftym,ixm+1,ith))
c     write(*,*)ix,ixm,qxfv(ix,ith),qxfv(ixm,ithm),
c     $        qfv(infty,ix,ith),qfv(inftym,ixm,ith)
         enddo
      endif
      end 

c******************************************************************
c Integrate along the vz direction, the n.v flux 
c of distribution fv(vx,vz) for specified vx.
      subroutine zintqfv(nzfvi,nzfva,vzfv1,th1,vx,qfv1,ztr1)
c Inputs
c The minimum and maximum vz indices
      integer nzfvi,nzfva
c The vz array
      real vzfv1(nzfvi:nzfva)
c The cosine(normal) and the vx velocity
      real th1,vx
c Output
c  the integral.
      real qfv1(nzfvi:nzfva)
c  the index of the zero-transition:
      real ztr1
c The trailing 1s remind us that this version is 1-d even if part of 
c a multi-d array is passed.

c A small number to prevent adopting exactly the end-points for ztr1
      eps=1.e-5
c  Direction normal
      zn=th1
      xn=sqrt(1-zn**2)
c For this integration we change the sign for negative cosines.
c Not necessary if doing abs(vdot).
      if(zn.lt.-1.e-8)then
         zn=-zn
         xn=-xn
      endif
      qfv10=0.
c Set initial values 
      qfv1(nzfvi)=0.
      ztr1=nzfvi+eps
      fvz=fv(vx,vzfv1(nzfvi))
      vdotn=xn*vx + zn*vzfv1(nzfvi)
      vdotn0=vdotn
      do iz=nzfvi+1,nzfva
         fvzm=fvz
         fvz=fv(vx,vzfv1(iz))
         if(fvz.lt.0.) then
            write(*,*)'zintqfv fv negative error',iz,fvz
         endif
         vdotnm=vdotn
         vdotn=xn*vx + zn*vzfv1(iz)
         if(vdotn*vdotnm.le.0 .and. vdotn.ne.0.) then
c Zero crossing index
            ztr1=((iz-1)*vdotn - iz*vdotnm)/(vdotn-vdotnm)
c Zero crossing velocity
            vz0=(vzfv1(iz-1)*vdotn - vzfv1(iz)*vdotnm)/(vdotn-vdotnm)
c Zero integral value
            qfv10=qfv1(iz-1) +
     $           (vz0-vzfv1(iz-1))*0.5*(fvzm*(vdotnm))
c            write(*,*)'Crossing',iz,vdotnm,vdotn,
c     $           fvz,fvzm,((iz-1)*vdotn-iz*vdotnm)/(vdotn-vdotnm)
c            write(*,*)'qm1,q,q10,diff=',qfv1(iz-1),qfv1(iz),qfv10,
c     $           qfv1(iz-1)-qfv10
            qfv1(iz)=qfv10 +
     $           (vzfv1(iz)-vz0)*0.5*(fvz*(vdotn))
         else
            qfv1(iz)=qfv1(iz-1) +
     $           (vzfv1(iz)-vzfv1(iz-1))*0.5*
     $           (fvz*(vdotn)+fvzm*(vdotnm))
         endif

      enddo
c Set qfv10 if it was not already set correctly.
      if(vdotn.le.0. .and. vdotn0.lt.0.) then
c         write(*,*) 'Adjusting right hand end',qfv1(nzfva),qfv10
         qfv10=qfv1(nzfva)
c We must avoid having exactly the maximum index (or minimum)
         ztr1=nzfva-eps
      endif
c Correct the zero level.
      do iz=nzfvi,nzfva
         qfv1(iz)=qfv1(iz)-qfv10
         if(qfv1(iz).lt.0.) then
            write(*,*)'zinitqfv q negative error at',iz,qfv1(iz)
         endif
      enddo

c      write(*,*)'th1=',th1,' vx=',vx,' ztr1=',ztr1
c 520  format('th1=',f7.4,' vx=',f7.3,' vz0=',f7.3,' q0=',e10.4,
c     $     ' q+,q-=',2e11.4)
c      write(*,520)th1,vx,vz0,qfv10,qfv1(nzfva),qfv1(nzfvi)
c      call autoplot(vzfv1,qfv1,nzfva-nzfvi+1)
c      call pltend()

      end
c***************************************************************
      subroutine calcfqx()
c integrate to get the cumulative vx distribution at fixed vz.
c This is the appropriate distribution to use for determining the
c ignorable coordinate only if the distribution fv(vx,vz) is independent
c of the third coordinate (y), i.e. separable in x,y. 
c But we have already implicitly assumed that the distribution is 
c isotropic in the x,y plane, and consequently a function of x^2+y^2.
c Therefore we require the distribution to be maxwellian exp(-vx^2-vy^2)
c Thus the apparent generality is not really justified.
c However, if distribution were not of this form, but still isotropic,
c then everwhere else, one could regard fv(vx,vz) as the distribution
c integrated over vy; while in the final ignorable coordinate choice
c we would have to specify both vx and vz, and the cumulative distrib
c we need would be three-dimensional.
      include 'fvcom.f'
      do j=nzfvi,nzfva
         fqvxvz(nxfvi,j)=0.
         fvkm=fv(vxfv(nxfvi),vzfv(j))
         do k=nxfvi+1,nxfva
            fvk=fv(vxfv(k),vzfv(j))
            fqvxvz(k,j)=fqvxvz(k-1,j)+
     $           (vxfv(k)-vxfv(k-1))*(fvk+fvkm)*0.5
            fvkm=fvk
         enddo
         if(fqvxvz(nxfva,j).ne.0.)then
            do k=nxfvi,nxfva
               fqvxvz(k,j)=fqvxvz(k,j)/fqvxvz(nxfva,j)
            enddo
         else
            write(*,*)'calcfqx: Warning, zero-integral inaccuracy!'
c hack a straight line
            do k=nxfvi,nxfva
               fqvxvz(k,j)=(k-nxfvi)/float(nxfva-nxfvi)
            enddo
         endif
      enddo
      end
c***************************************************************
c The distribution function, normalized so that integral over 
c normalized velocities is one. 2-d with y-coord ignorable.
      function fv(vx,vz)
c Velocities normalized to thermal sqrt(2T/m)
      real vx,vz
c Parameters of the distribution function
c At present this is a drift distribution corresponding to constant
c collision frequency charge-exchange, which has just one parameter,
c ud the normalized drift velocity.
c The perpendicular distribution is exp(-vx^2)/pi a Maxwellian.
      common /distfunc/ud
      fv=fvcx(vz,ud)*exp(-vx**2)/1.77245385
      end
c***************************************************************
c The distribution function, normalized so that integral over 
c normalized velocities is one. Gyrotropic
      function fvgyro(vx,vz)
c Velocities normalized to thermal sqrt(2T/m)
      real vx,vz
c Parameters of the distribution function
c At present this is a drift distribution corresponding to constant
c collision frequency charge-exchange, which has just one parameter,
c ud the normalized drift velocity.
c The perpendicular distribution is exp(-vx^2)/pi a Maxwellian.
      common /distfunc/ud
      fvgyro=fvcx(vz,ud)*exp(-vx**2)/3.1415926
      end
c****************************************************************
c FVCX function for 1-d drifting CX distribution.
      function fvcx(u,ud)
      real u,ud,v,vd,fvcx
c Return the normalized distribution function v_n f(v) for constant 
c cx collision frequency at a value of normalized velocity u=v/v_n,
c when the normalized drift velocity is ud= (a/\nu_c) /v_n,
c with v_n = sqrt(2T_n/m). a is acceleration, nu_c collision freq.
      if(ud.lt.0.) then
         v=-u
         vd=-ud
      else
         v=u
         vd=ud
      endif
      if(vd.eq.0.)then
         carg=20.
         earg=100
      else
         carg=0.5/vd-v
         earg=(0.5/vd)**2-v/vd
      endif
      if(carg.gt.10)then
c asymptotic form for large exp argument (small vd):
c  exp(-v^2)/[sqrt(\pi)(1-2 v_d v)]:
         fvcx=exp(-v**2)/1.77245385/(1.-2.*vd*v)
      elseif(carg.gt.-5.)then
         fvcx=exp(-v**2)*experfcc(carg)*0.5/vd
      else
c         fvcx=exp(earg)*erfcc(carg)*0.5/vd
         fvcx=exp(earg)/vd
      endif
c      write(*,*)'fvcx:vd,v,earg,fvcx',vd,v,earg,fvcx
      if(.not.fvcx.ge.0) then
         write(*,*)'fvcx error. u=',u,' ud=',ud,' f=',fvcx,carg
         fvcx=0.
      endif
      end
c****************************************************************
c (ERFCC is in randf.f) this is exp*erfc
      FUNCTION expERFCC(X)
      Z=ABS(X)      
      T=1./(1.+0.5*Z)
      expERFCC=T*EXP(-1.26551223+T*(1.00002368+T*(.37409196+
     *    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     *    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) expERFCC=2.*exp(z**2)-expERFCC
      END
c*******************************************************************

c********************************************************************
c Given a monotonic (increasing?) 
c function Q(x) on a 1-D grid x=1..nq, solve Q(x)=y for x.
c That is, invert Q to give x=Q^-1(y).
      subroutine f1invtfunc(Q,nq,y,x)
c Comment out this next declaration if you want Q to be a function
      real Q(nq)
      integer nq
      real y,x
c
      integer iqr,iql,iqx
      real Qx,Qr,Ql
      Ql=Q(1)
      Qr=Q(nq)
      iql=1
      iqr=nq
      if((y-Ql)*(y-Qr).gt.0.) then
c Value is outside the range.
         x=0
         return
      endif
 200  if(iqr-iql.eq.1)goto 210
      iqx=(iqr+iql)/2
      Qx=Q(iqx)
c      write(*,*)y,Ql,Qx,Qr,iql,iqr
c Formerly .lt. which is an error.
      if((Qx-y)*(Qr-y).le.0.) then
         Ql=Qx
         iql=iqx
      else
         Qr=Qx
         iqr=iqx
      endif
      goto 200
 210  continue
c Now iql and iqr, Ql and Qr bracket Q
      if(Qr-Ql.ne.0.)then
         x=(y-Ql)/(Qr-Ql)+iql
      else
         x=iql
         write(*,*)'****** Error!: finvtfunc coincident points'
      endif
      end
c**********************************************************************

c********************************************************************
c Given two monotonic functions (arrays) Q1(x) Q2(x)
c  on a 1-D grid x=1..nq, solve s*Q1(x)+t*Q2(x)=y for x.
c That is, invert Q=s Q1+ t Q2 to give x=Q^-1(y).
c Interpolating inside this routine reduces computational effort.
      subroutine f2invtfunc(Q1,Q2,nq,y,x,s,t)
      integer nq
c Comment out the next declaration if you want Q1, Q2 to be functions.
      real Q1(nq),Q2(nq)
      real y,x,s,t
c Statement function
      Q(j)=s*Q1(j)+t*Q2(j)
c Here on is just finvtfunc:
      integer iqr,iql,iqx
      real Qx,Qr,Ql
      Ql=Q(1)
      Qr=Q(nq)
      iql=1
      iqr=nq
      if((y-Ql)*(y-Qr).gt.0.) then
c Value is outside the range.
         x=0
         return
      endif
 200  if(iqr-iql.eq.1)goto 210
      iqx=(iqr+iql)/2
      Qx=Q(iqx)
c      write(*,*)y,Ql,Qx,Qr,iql,iqr
c Formerly .lt. which is an error.
      if((Qx-y)*(Qr-y).le.0.) then
         Ql=Qx
         iql=iqx
      else
         Qr=Qx
         iqr=iqx
      endif
      goto 200
 210  continue

c Now iql and iqr, Ql and Qr bracket Q
c Trap errors caused by flat sections.
      Qd=Qr-Ql
      if(Qd.eq.0.)then
         x=(iql+iqr)/2.
      else
         x=(y-Ql)/(Qr-Ql)+iql
      endif
      end
c**********************************************************************
c**********************************************************************
c Return the value of f(z) interpolated by index zi
      function finterp(f,zi,nf)
      integer nf
      real f(nf)
      real zi
      i=zi
      if(i.lt.1 .or. i.ge.nf)then 
         write(*,*)'***** Finterp error!',i,zi,nf
      endif
      fz=zi-i
      finterp=f(i)*(1.-fz)+f(i+1)*fz
      end
