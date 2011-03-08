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

      subroutine chargetomesh()
c Common data:
      include 'piccom.f'
      include 'errcom.f'

c      ninner=0

c     Zeroing the arrays take time if we zero everything all the
c     time. Just zero what needed. Cycle k,j,i rather than i,j,k for
c     cache issues
      
         
      do k=1,npsiused
         do j=1,nthused
            do i=1,nrused
               psum(i,j,k)=0.
            enddo
         enddo
      enddo


      if(diags) then
         do k=1,npsiused
            do j=1,nthused
               do i=1,nrused
                  vxsum(i,j,k)=0.
                  vysum(i,j,k)=0.
                  vzsum(i,j,k)=0.
               enddo
            enddo
         enddo
      endif
c Split the zeroying when diags is on, maybe faster with cache ??
      if(diags.or.samp) then
         do k=1,npsiused
            do j=1,nthused
               do i=1,nrused
                  vrsum(i,j,k)=0.
                  vtsum(i,j,k)=0.
                  vpsum(i,j,k)=0.
                  vr2sum(i,j,k)=0.
                  vt2sum(i,j,k)=0.
                  vp2sum(i,j,k)=0.
                  vrtsum(i,j,k)=0.
                  vrpsum(i,j,k)=0.
                  vtpsum(i,j,k)=0.
               enddo
            enddo
         enddo
      elseif(debyelen.eq.0) then
         do k=1,npsiused
            do j=1,nthused
               do i=1,2
                  vrsum(i,j,k)=0.
                  vr2sum(i,j,k)=0.
               enddo
            enddo
         enddo
      endif

c Perhaps this needs to be larger than npart for .not.lfixed.
c      write(*,*)'Starting chargetomesh',npart
      do i=1,iocprev
         if(ipf(i).gt.0)then
c         if(i.lt.10000)write(*,'(i6,$)')i
c Use fast ptomesh, half-quantities not needed.
            ih=0
            hf=99.
           
            call ptomesh(i,irl,rf,ithl,thf,ipl,pf,st,ct,sp,cp,rp
     $           ,zetap,ih,hf)
c           If lgotooutput tripped, skip to end
            if (lgotooutput) then
               write (*,*) 'lgotooutput tripped in chargetomesh.'
               goto 404
            endif

            if(rf.lt.0..or.rf.gt.1.)then
               rp=sqrt(xp(1,i)**2+xp(2,i)**2+xp(3,i)**2)
               write(*,*)'Outside mesh, rf error in chargetomesh',
     $              rf,irl,i,rp
            else
               call chargeassign(i,irl,rf,ithl,thf,
     $              ipl,pf,st,ct,sp,cp,rp)

            endif
         endif
      enddo

c      do iw=1,nrused
c         write(*,*) iw
c         write(*,*) ((vrsum(iw,jw,kw),jw=1,nthused),kw=1 ,npsiused)
c      enddo
 404  continue
      end
c***********************************************************************
c Accumulate particle charge into rho mesh and other diagnostics.
      subroutine chargeassign(i,irl,rf,ithl,thf,ipl,pf,st,ct,sp,cp,rp)
c      implicit none
      integer i
c Common data:
      include 'piccom.f'
      include 'errcom.f'

c Cyclic ipl in the poloidal direction
      integer iplp
      if(ipl.eq.npsi) then
         iplp=1
      else 
         iplp=ipl+1
      endif

c Assign as if cube for now. Volume weighting might be better.
c Charge summation.

      
      psum(irl,ithl,ipl)=psum(irl,ithl,ipl) + (1.-rf)*(1.-thf)*(1-pf)
      psum(irl+1,ithl,ipl)=psum(irl+1,ithl,ipl) + rf*(1.-thf)*(1-pf)
      psum(irl,ithl+1,ipl)=psum(irl,ithl+1,ipl) + (1.-rf)*thf*(1-pf)
      psum(irl+1,ithl+1,ipl)=psum(irl+1,ithl+1,ipl) + rf*thf*(1-pf)

      psum(irl,ithl,iplp)=psum(irl,ithl,iplp) + (1.-rf)*(1.-thf)*pf
      psum(irl+1,ithl,iplp)=psum(irl+1,ithl,iplp) + rf*(1.-thf)*pf
      psum(irl,ithl+1,iplp)=psum(irl,ithl+1,iplp) + (1.-rf)*thf*pf
      psum(irl+1,ithl+1,iplp)=psum(irl+1,ithl+1,iplp) + rf*thf*pf

c     calculate the sums only if needed, i.e. at the probe edge if Lde=0
c     or if we plot the diagrams (i.e. not -g)

      if(diags.or.samp.or.(debyelen.eq.0.and.irl.le.2)) then

      vxy=xp(4,i)*cp + xp(5,i)*sp
      vr=vxy*st + xp(6,i)*ct

      vrsum(irl,ithl,ipl)=vrsum(irl,ithl,ipl) + 
     $     (1.-rf)*(1.-thf)*(1-pf)*vr
      vrsum(irl+1,ithl,ipl)=vrsum(irl+1,ithl,ipl) + 
     $     rf*(1.-thf)*(1-pf)*vr
      vrsum(irl,ithl+1,ipl)=vrsum(irl,ithl+1,ipl) + 
     $     (1.-rf)*thf*(1-pf)*vr
      vrsum(irl+1,ithl+1,ipl)=vrsum(irl+1,ithl+1,ipl) + 
     $     rf*thf*(1-pf)*vr
      vrsum(irl,ithl,iplp)=vrsum(irl,ithl,iplp) + 
     $     (1.-rf)*(1.-thf)*pf*vr
      vrsum(irl+1,ithl,iplp)=vrsum(irl+1,ithl,iplp) + 
     $     rf*(1.-thf)*pf*vr
      vrsum(irl,ithl+1,iplp)=vrsum(irl,ithl+1,iplp) + 
     $     (1.-rf)*thf*pf*vr
      vrsum(irl+1,ithl+1,iplp)=vrsum(irl+1,ithl+1,iplp) + 
     $     rf*thf*pf*vr

      vr2=vr*vr

      vr2sum(irl,ithl,ipl)=vr2sum(irl,ithl,ipl) + 
     $     (1.-rf)*(1.-thf)*(1-pf)*vr2
      vr2sum(irl+1,ithl,ipl)=vr2sum(irl+1,ithl,ipl) + 
     $     rf*(1.-thf)*(1-pf)*vr2
      vr2sum(irl,ithl+1,ipl)=vr2sum(irl,ithl+1,ipl) + 
     $     (1.-rf)*thf*(1-pf)*vr2
      vr2sum(irl+1,ithl+1,ipl)=vr2sum(irl+1,ithl+1,ipl) + 
     $     rf*thf*(1-pf)*vr2
      vr2sum(irl,ithl,iplp)=vr2sum(irl,ithl,iplp) + 
     $     (1.-rf)*(1.-thf)*pf*vr2
      vr2sum(irl+1,ithl,iplp)=vr2sum(irl+1,ithl,iplp) + 
     $     rf*(1.-thf)*pf*vr2
      vr2sum(irl,ithl+1,iplp)=vr2sum(irl,ithl+1,iplp) + 
     $     (1.-rf)*thf*pf*vr2
      vr2sum(irl+1,ithl+1,iplp)=vr2sum(irl+1,ithl+1,iplp) + 
     $     rf*thf*pf*vr2

      endif

      if(diags.or.samp) then

      vt= vxy*ct - xp(6,i)*st

      vtsum(irl,ithl,ipl)=vtsum(irl,ithl,ipl) + 
     $     (1.-rf)*(1.-thf)*(1-pf)*vt
      vtsum(irl+1,ithl,ipl)=vtsum(irl+1,ithl,ipl) + 
     $     rf*(1.-thf)*(1-pf)*vt
      vtsum(irl,ithl+1,ipl)=vtsum(irl,ithl+1,ipl) + 
     $     (1.-rf)*thf*(1-pf)*vt
      vtsum(irl+1,ithl+1,ipl)=vtsum(irl+1,ithl+1,ipl) + 
     $     rf*thf*(1-pf)*vt
      vtsum(irl,ithl,iplp)=vtsum(irl,ithl,iplp) + 
     $     (1.-rf)*(1.-thf)*pf*vt
      vtsum(irl+1,ithl,iplp)=vtsum(irl+1,ithl,iplp) + 
     $     rf*(1.-thf)*pf*vt
      vtsum(irl,ithl+1,iplp)=vtsum(irl,ithl+1,iplp) + 
     $     (1.-rf)*thf*pf*vt
      vtsum(irl+1,ithl+1,iplp)=vtsum(irl+1,ithl+1,iplp) + 
     $     rf*thf*pf*vt


      vp=-xp(4,i)*sp  + xp(5,i)*cp

      vpsum(irl,ithl,ipl)=vpsum(irl,ithl,ipl) + 
     $     (1.-rf)*(1.-thf)*(1-pf)*vp
      vpsum(irl+1,ithl,ipl)=vpsum(irl+1,ithl,ipl) + 
     $     rf*(1.-thf)*(1-pf)*vp
      vpsum(irl,ithl+1,ipl)=vpsum(irl,ithl+1,ipl) + 
     $     (1.-rf)*thf*(1-pf)*vp
      vpsum(irl+1,ithl+1,ipl)=vpsum(irl+1,ithl+1,ipl) + 
     $     rf*thf*(1-pf)*vp
      vpsum(irl,ithl,iplp)=vpsum(irl,ithl,iplp) + 
     $     (1.-rf)*(1.-thf)*pf*vp
      vpsum(irl+1,ithl,iplp)=vpsum(irl+1,ithl,iplp) + 
     $     rf*(1.-thf)*pf*vp
      vpsum(irl,ithl+1,iplp)=vpsum(irl,ithl+1,iplp) + 
     $     (1.-rf)*thf*pf*vp
      vpsum(irl+1,ithl+1,iplp)=vpsum(irl+1,ithl+1,iplp) + 
     $     rf*thf*pf*vp

      vt2=vt*vt

      vt2sum(irl,ithl,ipl)=vt2sum(irl,ithl,ipl) + 
     $     (1.-rf)*(1.-thf)*(1-pf)*vt2
      vt2sum(irl+1,ithl,ipl)=vt2sum(irl+1,ithl,ipl) + 
     $     rf*(1.-thf)*(1-pf)*vt2
      vt2sum(irl,ithl+1,ipl)=vt2sum(irl,ithl+1,ipl) + 
     $     (1.-rf)*thf*(1-pf)*vt2
      vt2sum(irl+1,ithl+1,ipl)=vt2sum(irl+1,ithl+1,ipl) + 
     $     rf*thf*(1-pf)*vt2
      vt2sum(irl,ithl,iplp)=vt2sum(irl,ithl,iplp) + 
     $     (1.-rf)*(1.-thf)*pf*vt2
      vt2sum(irl+1,ithl,iplp)=vt2sum(irl+1,ithl,iplp) + 
     $     rf*(1.-thf)*pf*vt2
      vt2sum(irl,ithl+1,iplp)=vt2sum(irl,ithl+1,iplp) + 
     $     (1.-rf)*thf*pf*vt2
      vt2sum(irl+1,ithl+1,iplp)=vt2sum(irl+1,ithl+1,iplp) + 
     $     rf*thf*pf*vt2

      vp2=vp*vp

      vp2sum(irl,ithl,ipl)=vp2sum(irl,ithl,ipl) + 
     $     (1.-rf)*(1.-thf)*(1-pf)*vp2
      vp2sum(irl+1,ithl,ipl)=vp2sum(irl+1,ithl,ipl) + 
     $     rf*(1.-thf)*(1-pf)*vp2
      vp2sum(irl,ithl+1,ipl)=vp2sum(irl,ithl+1,ipl) + 
     $     (1.-rf)*thf*(1-pf)*vp2
      vp2sum(irl+1,ithl+1,ipl)=vp2sum(irl+1,ithl+1,ipl) + 
     $     rf*thf*(1-pf)*vp2
      vp2sum(irl,ithl,iplp)=vp2sum(irl,ithl,iplp) + 
     $     (1.-rf)*(1.-thf)*pf*vp2
      vp2sum(irl+1,ithl,iplp)=vp2sum(irl+1,ithl,iplp) + 
     $     rf*(1.-thf)*pf*vp2
      vp2sum(irl,ithl+1,iplp)=vp2sum(irl,ithl+1,iplp) + 
     $     (1.-rf)*thf*pf*vp2
      vp2sum(irl+1,ithl+1,iplp)=vp2sum(irl+1,ithl+1,iplp) + 
     $     rf*thf*pf*vp2

      vrt=vr*vt

      vrtsum(irl,ithl,ipl)=vrtsum(irl,ithl,ipl) + 
     $     (1.-rf)*(1.-thf)*(1-pf)*vrt
      vrtsum(irl+1,ithl,ipl)=vrtsum(irl+1,ithl,ipl) + 
     $     rf*(1.-thf)*(1-pf)*vrt
      vrtsum(irl,ithl+1,ipl)=vrtsum(irl,ithl+1,ipl) + 
     $     (1.-rf)*thf*(1-pf)*vrt
      vrtsum(irl+1,ithl+1,ipl)=vrtsum(irl+1,ithl+1,ipl) + 
     $     rf*thf*(1-pf)*vrt
      vrtsum(irl,ithl,iplp)=vrtsum(irl,ithl,iplp) + 
     $     (1.-rf)*(1.-thf)*pf*vrt
      vrtsum(irl+1,ithl,iplp)=vrtsum(irl+1,ithl,iplp) + 
     $     rf*(1.-thf)*pf*vrt
      vrtsum(irl,ithl+1,iplp)=vrtsum(irl,ithl+1,iplp) + 
     $     (1.-rf)*thf*pf*vrt
      vrtsum(irl+1,ithl+1,iplp)=vrtsum(irl+1,ithl+1,iplp) + 
     $     rf*thf*pf*vrt

      vrp=vr*vp

      vrpsum(irl,ithl,ipl)=vrpsum(irl,ithl,ipl) + 
     $     (1.-rf)*(1.-thf)*(1-pf)*vrp
      vrpsum(irl+1,ithl,ipl)=vrpsum(irl+1,ithl,ipl) + 
     $     rf*(1.-thf)*(1-pf)*vrp
      vrpsum(irl,ithl+1,ipl)=vrpsum(irl,ithl+1,ipl) + 
     $     (1.-rf)*thf*(1-pf)*vrp
      vrpsum(irl+1,ithl+1,ipl)=vrpsum(irl+1,ithl+1,ipl) + 
     $     rf*thf*(1-pf)*vrp
      vrpsum(irl,ithl,iplp)=vrpsum(irl,ithl,iplp) + 
     $     (1.-rf)*(1.-thf)*pf*vrp
      vrpsum(irl+1,ithl,iplp)=vrpsum(irl+1,ithl,iplp) + 
     $     rf*(1.-thf)*pf*vrp
      vrpsum(irl,ithl+1,iplp)=vrpsum(irl,ithl+1,iplp) + 
     $     (1.-rf)*thf*pf*vrp
      vrpsum(irl+1,ithl+1,iplp)=vrpsum(irl+1,ithl+1,iplp) + 
     $     rf*thf*pf*vrp

      vtp=vt*vp

      vtpsum(irl,ithl,ipl)=vtpsum(irl,ithl,ipl) + 
     $     (1.-rf)*(1.-thf)*(1-pf)*vtp
      vtpsum(irl+1,ithl,ipl)=vtpsum(irl+1,ithl,ipl) + 
     $     rf*(1.-thf)*(1-pf)*vtp
      vtpsum(irl,ithl+1,ipl)=vtpsum(irl,ithl+1,ipl) + 
     $     (1.-rf)*thf*(1-pf)*vtp
      vtpsum(irl+1,ithl+1,ipl)=vtpsum(irl+1,ithl+1,ipl) + 
     $     rf*thf*(1-pf)*vtp
      vtpsum(irl,ithl,iplp)=vtpsum(irl,ithl,iplp) + 
     $     (1.-rf)*(1.-thf)*pf*vtp
      vtpsum(irl+1,ithl,iplp)=vtpsum(irl+1,ithl,iplp) + 
     $     rf*(1.-thf)*pf*vtp
      vtpsum(irl,ithl+1,iplp)=vtpsum(irl,ithl+1,iplp) + 
     $     (1.-rf)*thf*pf*vtp
      vtpsum(irl+1,ithl+1,iplp)=vtpsum(irl+1,ithl+1,iplp) + 
     $     rf*thf*pf*vtp

      endif

      if(ldistf) then

c     Distribution function diagnostics

c     Find velocity bins
      call invtarray(vxbins,2,nvx+1,xp(4,i),ivx)
      ivx = ivx-1
      call invtarray(vybins,2,nvy+1,xp(5,i),ivy)
      ivy = ivy-1
      call invtarray(vzbins,2,nvz+1,xp(6,i),ivz)
      ivz = ivz-1

c     Find spatial bins
      ix = min(nx,max(1,1+(xp(1,i)/rcc(nrused)+1)/2*nx))
      iy = min(nx,max(1,1+(xp(2,i)/rcc(nrused)+1)/2*ny))
      iz = min(nx,max(1,1+(xp(3,i)/rcc(nrused)+1)/2*nz))

c     Increment counter for corresponding bin
      vxdistf(ivx,ix,iy,iz) = vxdistf(ivx,ix,iy,iz) + 1;
      vydistf(ivy,ix,iy,iz) = vydistf(ivy,ix,iy,iz) + 1;
      vzdistf(ivz,ix,iy,iz) = vzdistf(ivz,ix,iy,iz) + 1;


      endif

      if(diags)then

      vx=xp(4,i)
      vxsum(irl,ithl,ipl)=vxsum(irl,ithl,ipl) + 
     $     (1.-rf)*(1.-thf)*(1-pf)*vx
      vxsum(irl+1,ithl,ipl)=vxsum(irl+1,ithl,ipl) + 
     $     rf*(1.-thf)*(1-pf)*vx
      vxsum(irl,ithl+1,ipl)=vxsum(irl,ithl+1,ipl) + 
     $     (1.-rf)*thf*(1-pf)*vx
      vxsum(irl+1,ithl+1,ipl)=vxsum(irl+1,ithl+1,ipl) + 
     $     rf*thf*(1-pf)*vx
      vxsum(irl,ithl,iplp)=vxsum(irl,ithl,iplp) + 
     $     (1.-rf)*(1.-thf)*pf*vx
      vxsum(irl+1,ithl,iplp)=vxsum(irl+1,ithl,iplp) + 
     $     rf*(1.-thf)*pf*vx
      vxsum(irl,ithl+1,iplp)=vxsum(irl,ithl+1,iplp) + 
     $     (1.-rf)*thf*pf*vx
      vxsum(irl+1,ithl+1,iplp)=vxsum(irl+1,ithl+1,iplp) + 
     $     rf*thf*pf*vx

      vy=xp(5,i)
      vysum(irl,ithl,ipl)=vysum(irl,ithl,ipl) + 
     $     (1.-rf)*(1.-thf)*(1-pf)*vy
      vysum(irl+1,ithl,ipl)=vysum(irl+1,ithl,ipl) + 
     $     rf*(1.-thf)*(1-pf)*vy
      vysum(irl,ithl+1,ipl)=vysum(irl,ithl+1,ipl) + 
     $     (1.-rf)*thf*(1-pf)*vy
      vysum(irl+1,ithl+1,ipl)=vysum(irl+1,ithl+1,ipl) + 
     $     rf*thf*(1-pf)*vy
      vysum(irl,ithl,iplp)=vysum(irl,ithl,iplp) + 
     $     (1.-rf)*(1.-thf)*pf*vy
      vysum(irl+1,ithl,iplp)=vysum(irl+1,ithl,iplp) + 
     $     rf*(1.-thf)*pf*vy
      vysum(irl,ithl+1,iplp)=vysum(irl,ithl+1,iplp) + 
     $     (1.-rf)*thf*pf*vy
      vysum(irl+1,ithl+1,iplp)=vysum(irl+1,ithl+1,iplp) + 
     $     rf*thf*pf*vy

      vz=xp(6,i)
      vzsum(irl,ithl,ipl)=vzsum(irl,ithl,ipl) + 
     $     (1.-rf)*(1.-thf)*(1-pf)*vz
      vzsum(irl+1,ithl,ipl)=vzsum(irl+1,ithl,ipl) + 
     $     rf*(1.-thf)*(1-pf)*vz
      vzsum(irl,ithl+1,ipl)=vzsum(irl,ithl+1,ipl) + 
     $     (1.-rf)*thf*(1-pf)*vz
      vzsum(irl+1,ithl+1,ipl)=vzsum(irl+1,ithl+1,ipl) + 
     $     rf*thf*(1-pf)*vz
      vzsum(irl,ithl,iplp)=vzsum(irl,ithl,iplp) + 
     $     (1.-rf)*(1.-thf)*pf*vz
      vzsum(irl+1,ithl,iplp)=vzsum(irl+1,ithl,iplp) + 
     $     rf*(1.-thf)*pf*vz
      vzsum(irl,ithl+1,iplp)=vzsum(irl,ithl+1,iplp) + 
     $     (1.-rf)*thf*pf*vz
      vzsum(irl+1,ithl+1,iplp)=vzsum(irl+1,ithl+1,iplp) + 
     $     rf*thf*pf*vz

      endif
      
      end
c***********************************************************************
c Calculate potential phi from rho.
      subroutine fcalc(dt)
c Common data:
      include 'piccom.f'
      include 'errcom.f'
      real phi0mphi1(0:nthsize,0:npsisize),delphi0(0:nthsize,0:npsisize)
      parameter (nthp1=nthsize+1)
      parameter (npsip1=npsisize+1)
      parameter (ntot=nthp1*npsip1)
c      real phi1ave
      real bcifac,bcpfac,bci,bcp,bvf
      real relax
      real cs(nthsize,npsisize),csd(nthsize,npsisize)
      real vs,vsd(nthsize,npsisize)
      real psum_r(2,nthsize,npsisize)
      real vrsum_r(2,nthsize,npsisize)
      real vr2sum_r(2,nthsize,npsisize)
      real ncs
      logical first
      integer kk1,kk2
      data relax/1./
c Reduced with respect to 2D
      data bcifac/.05/bcpfac/.025/
      data bvf/1.2071/
      data first/.true./
      data phi0mphi1/ntot*0./
      data delphi0/ntot*0./
      data ncs/50./
      real dr
      save
      cerr=0.

      dr=rcc(2)-rcc(1)

c Calculate the potential on the grid according to phi=log(rho)
      do k=1,npsiused
         do j=1,nthused
            relax=1.
            do i=3,nrused
               if(rho(i,j,k).le.0.)then
                  write(*,*)'rho=0',i,j,k
c                  stop
c                 Trigger go to output
                  lgotooutput = .true.
               endif
c Relaxed Boltzmann scheme.
               delta=phi(i,j,k)-log(rho(i,j,k))
               if(abs(delta).le.cerr)cerr=abs(delta)
               phi(i,j,k)=phi(i,j,k)-relax*delta
            enddo

            do i=1,2
               if(rho(i,j,k).le.0.)then
                  write(*,*)'rho=0',i,j,k
c                  stop
c                 Trigger go to output
                  lgotooutput = .true.
               endif
               relax=1.               
c Relaxed Boltzmann scheme.
               delta=phi(i,j,k)-log(rho(i,j,k))
               if(abs(delta).le.cerr)cerr=abs(delta)
               phi(i,j,k)=phi(i,j,k)-relax*delta
               
c want at least 40 particles for average, and relax<0.25
               relax=(40/psum(i,j,k)-1)*dr/sqrt(1+Ti)/dt
               relax=1/max(4.,relax)

               delta=psum_r(i,j,k)-psum(i,j,k)
               psum_r(i,j,k)=psum_r(i,j,k)-relax*delta
               delta=vrsum_r(i,j,k)-vrsum(i,j,k)
               vrsum_r(i,j,k)=vrsum_r(i,j,k)-relax*delta
               delta=vr2sum_r(i,j,k)-vr2sum(i,j,k)
               vr2sum_r(i,j,k)=vr2sum_r(i,j,k)-relax*delta
            enddo

         enddo
      enddo

c Probe boundary condition.
      if(first)then
         do j=1,nthused
            do k=1,npsiused
               cs(j,k)=-sqrt(1.+Ti)
            enddo
         enddo
         first=.false.
      endif


c If we impose the bohm condition
      if(bohm) then


      do j=1,nthused
         do k=1,npsiused
c     calculate the exact density at the sheath edge at this precise
c     time-step. Indeed, the sheath density calculated in rhocalc is at
c     the previous step.
            vs=vrincellave(j,k)/(1e-6+fincellave(j,k))
            fluxofangle=fincellave(j,k)*(nthused-1)*npsiused/(4*pi
     $           *dt*r(1)**2)
            if(fluxofangle.ne.0) then
               s0=-fluxofangle/vs
               rho(0,j,k)=s0/rhoinf
               if(j.eq.1.or.j.eq.nthused) rho(0,j,k)=2*rho(0,j,k)
               rho(1,j,k)=rho(0,j,k)
               s0=s0/((nthused-1)*npsiused)
               p0=((vr2incellave(j,k)*fincellave(j,k)-vrincellave(j,k)
     $              **2)/(fincellave(j,k)**2))*s0
            else
               s0=0
               p0=0
            endif
            

c            s1=psum_r(1,j,k)
c            p1=vr2sum_r(1,j,k)*s1-vrsum_r(1,j,k)**2
            s2=psum_r(2,j,k)
            p2=vr2sum_r(2,j,k)*s2-vrsum_r(2,j,k)**2

c            if(p1.ne.0)then
c               if(s1.gt.0)then
c                  p1=(p1/s1)*volinv(1)
c                  v1=vrsum(1,j,k)/(psum(1,j,k)+1e-6)
c                  s1=s1*volinv(1)
c               else
c                  write(*,*)'s1=0'
c                  stop
c               endif
c            endif
            if(p2.ne.0)then
               if(s2.gt.0)then
                  p2=(p2/s2)*volinv(2)
                  v2=vrsum(2,j,k)/(psum(2,j,k)+1e-6)
                  s2=s2*volinv(2)
               else
                  write(*,*)'s2=0'
c                  stop
c                  Trigger go to output
                   lgotooutput = .true.
               endif
            endif

c Old way
c            vs=(1.+bvf)*v1-bvf*v2
c            csd(j,k)=(p2 - p1)/(s2-s1)+1.

c New way
            vs=vrincellave(j,k)/(1e-6+fincellave(j,k))
            csd(j,k)=(p2 - p0)/(s2-s0)+1.

c     There is a strong sensitivity issue when s1~s2. So use the fact
c     that cs must be larger than csiso, and smaller than csadia (gamma
c     is between 1 (isotherm) and 3 (1D adiabatic))
            csiso=p0/(s0+1e-6) + 1.
            csadia=3.*p0/(s0+1e-6) + 1.     
            if(csd(j,k).lt.csiso)then
               csd(j,k)=csiso
            elseif(csd(j,k).gt.csadia)then
               csd(j,k)=csadia
            elseif(abs(s2-s0).le.1e-6)then
               csd(j,k)=0.5*(csio+csadia)
            endif
               
            csd(j,k)=-sqrt(csd(j,k))

c     Clip the excursion symmetrically. (Be careful that if csd ever
c     goes to zero, it will stay to zero forever)
            if(csd(j,k).lt.2.*cs(j,k))  csd(j,k)=2.*cs(j,k)

c Average the sound-speed value over ncs steps.
            cs(j,k)=((ncs-1.)*cs(j,k)+csd(j,k))/ncs
            if(cs(j,k).gt.0)then
               write(*,*)'cs positive',j,cs(j,k),csd(j,k)
     $              ,p0,p2,psum(1,j,k),psum(2,j,k),ncs
c               stop
c           Trigger go to output
            lgotooutput = .true.
            endif

c I am not sure I understand the cs**2 (except to put the sign back to positive
            bci=-bcifac*cs(j,k)**2*dt/dr
            bcp=-bcpfac*dr/(dt*cs(j,k))
            delphinew=(vs - cs(j,k))*bci

            phidrop=phi0mphi1(j,k)+delphinew +(delphinew -delphi0(j,k))
     $           *bcp

            phi(0,j,k)=phi(1,j,k)+phidrop

c     The following line says that the artificial potential of the first
c     cell can not be higher than what is calculated using
c     log(rho). Hence sometimes (almost always in fact), the
c     potential plot at rcc=1 is capped at log(rho(rcc=1)).

c     Physically, this is because if the velocity is already higher than
c     the sound speed, Bohm condition is already satisfied
            if(phidrop.gt.0)then
               phi(0,j,k)=phi(0,j,k)-phidrop
            endif

            
c     the sheath entrance potential must be bounded negatively,
c     otherwise in the trailing with E*B field we can get accel
c     overflow. Physically, there is no reason for the sheath potential
c     to be much more negative than ~-(Te+Ek_inf+|vd*Bz|*rp)
            if(phi(0,j,k).le.-3*(1+0.5*vd**2+abs(vd*Bz)))then
               phi(0,j,k)=-3*(1+0.5*vd**2+abs(vd*Bz))
            endif

            delphi0(j,k)=delphinew
            phi0mphi1(j,k)=phi(0,j,k)-phi(1,j,k)

c Adjusting the potential of the first cell.
            phi(1,j,k)=phi(0,j,k)

         enddo
      enddo

      endif

c     We must average the potential of each psi-cell at theta=0 or
c     theta=pi. The boundary conditions assume the cell center is on
c     axis, while what is indeed calculated is the potential at the
c     center of the first/last theta-cells.

      do i=2,nrused
         psiave1=0.
         psiave2=0.
         do k=1,npsiused
            psiave1=psiave1+phi(i,1,k)
            psiave2=psiave2+phi(i,nthused,k)
         enddo
         psiave1=psiave1/npsiused
         psiave2=psiave2/npsiused
         do k=1,npsiused
            phi(i,1,k)=psiave1
            phi(i,nthused,k)=psiave2
         enddo
      enddo


c We must set the potential of the shadow theta-cells to the
c physical value of the potential at psi+pi. The zero-derivative
c condition at theta=0,pi is not valid anymore in the 3D case
      
      do k=1,npsiused
c kk1=kk2 if npsi is even
         kk1=mod(k+3*npsi/2-1,npsi)+1
         kk2=mod(k+(3*npsi+1)/2-1,npsi)+1
         do i=1,nrused
            phi(i,0,k)=0.5*(phi(i,2,kk1)+phi(i,2,kk2))
            phi(i,nthused+1,k)=
     $           0.5*(phi(i,nthused-1,kk1)+phi(i,nthused-1,kk2))
         enddo
      enddo
c Set the theta shadow cells to their proper value to ensure periodicity
      do i=1,nrused
         do j=0,nthused+1
            phi(i,j,npsiused+1)=phi(i,j,1)
            phi(i,j,0)=phi(i,j,npsiused)
         enddo
      enddo



c      write(*,*)'At end cs=',(cs(kk),kk=1,nth)
c      write(*,'(10f7.4)')((phi(i,ih),i=0,9),ih=0,9)
      cerr=1.

      end
c***********************************************************************
c cic version.
      subroutine getaccel(i,accel,il,rf,ith,tf,ipl,pf,
     $     st,ct,sp,cp,rp,zetap,ih,hf)
c Evaluate the cartesian acceleration into accel. Using half-mesh
c parameters.
c accel is minus the gradient of phi for the ith particle.
c Be careful with variables in this routine.

      implicit none
      integer i
      real accel(3)
c Common data:
      include 'piccom.f'
      include 'errcom.f'
c Radial, Theta and Phi accelerations
      real ar,at,ap
      real ct,st,cp,sp,rp
      real zetap,hf
      integer ih
      real dp,dth,dpinv,dthinv
      integer ii,jj
      integer ipl,ith
      integer il,ir,ithp1,ithp2,ithm1,ilm1
      integer iplp1,iplp2,iplm1
      real rlm1,rf,tf,pf,rr,rl,dpsi,tflin

      real philm1tt,philm1pt,philm1tp,philm1pp
      real philm1t,philm1p

      real phihp1tt,phihp1pt,phihp1mt,phihp12t
      real phihp1tp,phihp1pp,phihp1mp,phihp12p
      real phihp1tm,phihp1pm
      real phihp1t2,phihp1p2

      real phihp1t,phihp1p,phihp1m,phihp12
      real phihp1tX,phihp1pX,phihp1mX,phihp12X

      real phih1tt,phih1tp,phih1tm, phih1t2
      real phih1pt,phih1pp,phih1pm,phih1p2
      real phih1mt,phih1mp
      real phih12t,phih12p

      real phih1t,phih1p,phih1m,phih12
      real phih1tX,phih1pX,phih1mX,phih12X

      real sinhere,sinplus

      data dp/0./


c Don't use parameter statements. Just set values first time
      if(dp.eq.0.)then
         dp=2.*pi/np
         dth=pi/(nth-1)
         dpinv=np/2./pi
         dthinv=(nth-1)/pi
      endif

      rl=r(ih)
      ir=ih+1
      ilm1=ih-1
c      if(hf.gt.1. .or. hf.lt.0.)write(*,*)'hf=',hf

c     Theta indexes
      ithp1=ith+1
      ithp2=ith+2
      ithm1=ith-1

c     Psi indexes
      iplp1=mod(ipl,npsiused)+1
      iplp2=mod(ipl+1,npsiused)+1
      iplm1=mod(ipl+npsiused-2,npsiused)+1


c     phihXXXX: phi at i=ih
c     phihpXXXX : phi at i=ih+1
c     philmXXXX: phi at i=ih-1

c     The first 1 reminds us that we have a linear interpolation in the angles
c     Then the letter (t,p,m,2) is for theta
c     Then the letter (t,p,m,2) is for psi


c     Potential al i=ih
      phih1tt=phi(ih,ith,ipl)
      phih1tp=phi(ih,ith,iplp1)
      phih1tm=phi(ih,ith,iplm1)
      phih1t2=phi(ih,ith,iplp2)
      phih1pt=phi(ih,ithp1,ipl)
      phih1pp=phi(ih,ithp1,iplp1)
      phih1pm=phi(ih,ithp1,iplm1)
      phih1p2=phi(ih,ithp1,iplp2)
      phih1mt=phi(ih,ithm1,ipl)
      phih1mp=phi(ih,ithm1,iplp1)
      phih12t=phi(ih,ithp2,ipl)
      phih12p=phi(ih,ithp2,iplp1)
      
c     Weighted in the psi direction
      phih1t=phih1tp*pf+phih1tt*(1-pf)
      phih1p=phih1pp*pf+phih1pt*(1-pf)
      phih1m=phih1mp*pf+phih1mt*(1-pf)
      phih12=phih12p*pf+phih12t*(1-pf)

c     Theta weighting tf is the fractional theta position interpolated
c     linearly in cos(theta). But to calculate ap correctly (in
c     particular close to the axis), we need tflin, fractional theta
c     position interpolated linearly in theta.
      tflin=(acos(ct)-thang(ith))/(thang(ith+1)-thang(ith))

      phih1tX=phih1pt*tflin +phih1tt*(1-tflin)
      phih1pX=phih1pp*tflin +phih1tp*(1-tflin)
      phih1mX=phih1pm*tflin +phih1tm*(1-tflin)
      phih12X=phih1p2*tflin +phih1t2*(1-tflin)


c Potentials at i=ih=il+1 for neighboring theta and psi
      if(ih.eq.nr)then
c Deal with r-boundary conditions.
c Constant slope at r-boundary. (Zero second derivative).
         phihp1tt=2*phih1tt-phi(ih-1,ith,ipl)
         phihp1pt=2*phih1pt-phi(ih-1,ithp1,ipl)
         phihp1mt=2*phih1mt-phi(ih-1,ithm1,ipl)
         phihp12t=2*phih12t-phi(ih-1,ithp2,ipl)

         phihp1tp=2*phih1tp-phi(ih-1,ith,iplp1)
         phihp1pp=2*phih1pp-phi(ih-1,ithp1,iplp1)
         phihp1mp=2*phih1mp-phi(ih-1,ithm1,iplp1)
         phihp12p=2*phih12p-phi(ih-1,ithp2,iplp1)

         phihp1tm=2*phih1tm-phi(ih-1,ith,iplm1)
         phihp1pm=2*phih1pm-phi(ih-1,ithp1,iplm1)

         phihp1t2=2*phih1t2-phi(ih-1,ith,iplp2)
         phihp1p2=2*phih1p2-phi(ih-1,ithp1,iplp2)

         rr=2*rl-r(ih-1)

      else
         phihp1tt=phi(ir,ith,ipl)
         phihp1tp=phi(ir,ith,iplp1)
         phihp1tm=phi(ir,ith,iplm1)
         phihp1t2=phi(ir,ith,iplp2)
         phihp1pt=phi(ir,ithp1,ipl)
         phihp1pp=phi(ir,ithp1,iplp1)
         phihp1pm=phi(ir,ithp1,iplm1)
         phihp1p2=phi(ir,ithp1,iplp2)
         phihp1mt=phi(ir,ithm1,ipl)
         phihp1mp=phi(ir,ithm1,iplp1)
         phihp12t=phi(ir,ithp2,ipl)
         phihp12p=phi(ir,ithp2,iplp1)

         rr=r(ir)
      endif

c     For the radial and theta acceleration, use values of the potential
c     already weighted in the psi direction
      phihp1t=phihp1tp*pf+phihp1tt*(1-pf)
      phihp1p=phihp1pp*pf+phihp1pt*(1-pf)
      phihp1m=phihp1mp*pf+phihp1mt*(1-pf)
      phihp12=phihp12p*pf+phihp12t*(1-pf)

c     Theta weighting
      phihp1tX=phihp1pt*tflin +phihp1tt*(1-tflin)
      phihp1pX=phihp1pp*tflin +phihp1tp*(1-tflin)
      phihp1mX=phihp1pm*tflin +phihp1tm*(1-tflin)
      phihp12X=phihp1p2*tflin +phihp1t2*(1-tflin)



c------------------------------
c Deal with radial acceleration
c -----------------------------
      
c Here we control whether we use the zeta or r interpolation.
      if(debyelen.lt.1.e-2)then
c Linear approx to sqrt form at boundary.
         if(ih.eq.1)then

c This is the way it was done in SCEPTIC2D. I don't like the arbirary bdyfc
c            philm1tt=phi(il,ith,ipl)-bdyfc*sqrt(2.*(rr-rl))*0.25
c            philm1pt=phi(il,ithp1,ipl)-bdyfc*sqrt(2.*(rr-rl))*0.25
c            philm1tp=phi(il,ith,iplp1)-bdyfc*sqrt(2.*(rr-rl))*0.25
c            philm1pp=phi(il,ithp1,iplp1)-bdyfc*sqrt(2.*(rr-rl))*0.25

c     Instead, mirror phi(2,j,k). We can not use phi(0,j,k) because in
c     the quasineutral regime, it represents a dummy variable used in
c     fcalc.
            philm1tt=2*phi(il,ith,ipl)-phi(ir,ith,ipl)
            philm1pt=2*phi(il,ithp1,ipl)-phi(ir,ithp1,ipl)
            philm1tp=2*phi(il,ith,iplp1)-phi(ir,ith,iplp1)
            philm1pp=2*phi(il,ithp1,iplp1)-phi(ir,ithp1,iplp1)

            rlm1=2.*rl - rr
c     Constant slope
         else
            philm1tt=phi(ilm1,ith,ipl)
            philm1pt=phi(ilm1,ithp1,ipl)
            philm1tp=phi(ilm1,ith,iplp1)
            philm1pp=phi(ilm1,ithp1,iplp1)
            rlm1=r(ilm1)
         endif
c     Psi weighting
         philm1t=philm1tp*pf+philm1tt*(1-pf)
         philm1p=philm1pp*pf+philm1pt*(1-pf)

c Uniform interpolation of d\phi/d\zeta times 1/\zeta
         if(zetap.le.1.e-2)zetap=1.e-2
         ar=(  ( (phihp1t-phih1t)/(zeta(ir)-zeta(ih))*hf +
     $        (phih1t-philm1t)/(zeta(ih)-zeta(ilm1))*(1.-hf)
     $           )*(1.-tflin)
     $        +( (phihp1p-phih1p)/(zeta(ir)-zeta(ih))*hf +
     $        (phih1p-philm1p)/(zeta(ih)-zeta(ilm1))*(1.-hf)
     $         )*tflin )/zetap
      else
c Interpolate in r.
         philm1tt=phi(ilm1,ith,ipl)
         philm1pt=phi(ilm1,ithp1,ipl)
         philm1tp=phi(ilm1,ith,iplp1)
         philm1pp=phi(ilm1,ithp1,iplp1)
         rlm1=r(ilm1)
c     Psi weighting
         philm1t=philm1tp*pf+philm1tt*(1-pf)
         philm1p=philm1pp*pf+philm1pt*(1-pf)

         ar=(  ( (phihp1t-phih1t)/(rr-rl)*hf +
     $        (phih1t-philm1t)/(rl-rlm1)*(1.-hf)
     $           )*(1.-tflin)
     $        +( (phihp1p-phih1p)/(rr-rl)*hf +
     $        (phih1p-philm1p)/(rl-rlm1)*(1.-hf)
     $         )*tflin )

      endif

c The minus sign is because E=-grad Phi         
      ar=-ar


c------------------------------
c Deal with theta acceleration
c -----------------------------

c     In sceptic2D, we always divide by rl. It think it is wrong. One
c     should divide by rl or rr depending on the radius at which the
c     potential derivative is calculated

      if(tflin.le.0.5)then
         at= ( (phih1p-phih1t)*(tflin)*2.
     $          /(rl*(thang(ithp1)-thang(ith)))
     $        +(phih1p-phih1m)*(0.5-tflin)*2.
     $          /(rl*(thang(ithp1)-thang(ithm1))) ) * (1.-rf)
     $        + ( (phihp1p-phihp1t)*(tflin)*2.
     $          /(rr*(thang(ithp1)-thang(ith)))
     $        +(phihp1p-phihp1m)*(0.5-tflin)*2.
     $          /(rr*(thang(ithp1)-thang(ithm1))) ) * rf
      else
         at= ( (phih12-phih1t)*(tflin-0.5)*2.
     $          /(rl*(thang(ithp2)-thang(ith)))
     $        +(phih1p-phih1t)*(1.-tflin)*2.
     $          /(rl*(thang(ithp1)-thang(ith))) ) * (1.-rf)
     $        + ( (phihp12-phihp1t)*(tflin-0.5)*2.
     $          /(rr*(thang(ithp2)-thang(ith)))
     $        +(phihp1p-phihp1t)*(1.-tflin)*2.
     $          /(rr*(thang(ithp1)-thang(ith))) ) * rf
      endif

      
      at=-at
      if(lat0)at=0.

c------------------------------
c Deal with psi acceleration
c -----------------------------


c Uniform Psi spacing
      dpsi=pcc(2)-pcc(1)     

      if(pf.le.0.5)then
         ap= ( (phih1pX-phih1tX)*(pf)*2./rl +0.5*(phih1pX-phih1mX)*(0.5
     $        -pf)*2./rl ) * (1.-rf) + ( (phihp1pX-phihp1tX)*(pf)*2./rr
     $        +0.5*(phihp1pX-phihp1mX)*(0.5-pf)*2./rr ) * rf
      else
         ap= ( 0.5*(phih12X-phih1tX)*(pf-0.5)*2. /rl +(phih1pX-phih1tX)
     $        *(1.-pf)*2./rl ) * (1.-rf) + ( 0.5*(phihp12X-phihp1tX)*(pf
     $        -0.5)*2./rr +(phihp1pX-phihp1tX)*(1.-pf)*2./rr ) * rf
      endif



c     divide by dpsi*st, because dl=r*sin(theta)*dpsi, and r has already
c     been taken into account in the previous ap calculation. The minus
c     sign is because E=-grad(phi)
      ap=-ap/(st*dpsi+1e-7)
      if(lap0)ap=0.

 501  format(a,6f10.4)


c 3D acceleration
      accel(3)=ar*ct - at*st
      accel(2)=(ar*st+ at*ct)*sp+ap*cp
      accel(1)=(ar*st+ at*ct)*cp-ap*sp
c Trap errors.

      if(.not.abs(accel(1)).lt.1.e5)then
         write(*,*) 'i: ',i,' x: ',xp(1,i),' y: ',xp(2,i),' z: ',xp(3,i)
         write(*,*)'Accel Excessive: ar,at,ap,ct,sp,cp,'
         write(*,*) ar,at,ap,ct,sp,cp
         write(*,*) 'phi at ith and ithp1:'
         write(*,*) philm1t,phih1t,phihp1t
         write(*,*) philm1p,phih1p,phihp1p
         write(*,*) 'zetap=',zetap,'  bdyfc=',bdyfc
         write(*,*) 'zeta '
         write(*,*) zeta(ilm1),zeta(ih),zeta(ir)
         write(*,*) 'ih,ith,ipl,iplp1,hf,rf,tflin,pf'
         write(*,*) ih,ith,ipl,iplp1,hf,rf,tflin,pf
         write(*,'(10f8.4)')((phi(ii,jj,ipl),ii=1,10),jj=1,10)
c         stop
c        Trigger go to output
         lgotooutput = .true.
      endif



      end
c***********************************************************************

c**********************************************************************
      subroutine esforce(ir,qp,fz,epz,fbz,fx,epx,fbx,fy,epy,fby)
      include 'piccom.f'
      include 'errcom.f'
c 3D version of esforce
c Return the charge qp, esforce fz,x,y , and electron pressure force epz,x,y.


c Coefficients for fz      
      real ercoefZ(nthsize),etcoefZ(nthsize),epcoefZ(nthsize)
     $     ,ertcoefZ(nthsize)
c Coefficients for fx,y
      real ercoefX(nthsize),etcoefX(nthsize),epcoefX(nthsize)
     $     ,ertcoefX(nthsize) ,erpcoefX(nthsize)
c Coefficients for the charge qp
      real qpcoef(nthsize)

c Psi angular spacing
      real dpsi
c sin theta between four cell centers, interpolated in theta and not cos(theta)
      real shalf

      real sp,cp,spp,cpp,phihere,sd,sb
      real vx,vy,vz,partsum,frac

      logical lnotinit
      data lnotinit/.true./
      
      save

      dpsi=pcc(2)-pcc(1)

      if(lnotinit)then
c Initialize coefficient arrays
         do j=1,nthused-1
        
            qpcoef(j)=-(th(j+1)-th(j))

            ercoefZ(j)=-0.25*(th(j+1)**2-th(j)**2)
            etcoefZ(j)=0.5*(th(j+1)**2/2.-th(j+1)**4/4. -th(j)**2/2.
     $           +th(j)**4/4.)/(th(j+1)-th(j))**2
            epcoefZ(j)=-ercoefZ(j)/dpsi**2
            ertcoefZ(j)=-(1.-(th(j+1)*(th(j+1)+th(j))+th(j)**2)/3.)

            ercoefX(j)=0.25*((thang(j+1)-th(j+1)*sqrt(1-th(j+1)**2)
     $           )-(thang(j)-th(j)*sqrt(1-th(j)**2)))       
            etcoefX(j)=-0.125*((3.*thang(j+1)/2-th(j+1)*sqrt(1-th(j+1
     $           )**2)*(3./2+1.-th(j+1)**2))- (3.*thang(j)/2-th(j)
     $           *sqrt(1-th(j)**2)*(3./2+1.-th(j)**2)))/(th(j+1)-th(j))
     $           **2
            epcoefX(j)=-ercoefX(j)/dpsi**2

c Old ertcoefX. Only valid to first order in ntheta            
            ertcoefX(j)=-1/(3.*(th(j+1)-th(j))) *((1-th(j+1)**2)**(3./2)
     $           -(1-th(j)**2)**(3./2))
c     Alternative version of ertcoefX (correct to second order)
c            shalf=sin(0.5*(thang(j+1)+thang(j)))
c            ertcoefX(j)=-1/(2.*(th(j+1)-th(j)))*shalf*((1-th(j+1)**2)
c     $           -(1-th(j)**2))

            erpcoefX(j)=-1./dpsi*(th(j+1)-th(j))


         enddo

c Set the coeffs at nthused to zero
         j=nthused

         qpcoef(j)=0.
         ercoefZ(j)=0.
         etcoefZ(j)=0.
         epcoefZ(j)=0.
         ertcoefZ(j)=0.
         ercoefX(j)=0.
         etcoefX(j)=0.
         epcoefX(j)=0.
         ertcoefX(j)=0.
         erpcoefX(j)=0.
         
         lnotinit=.false.
      endif

c Calculate charge and electrostatic force on probe surface
c We do this at a specified radius node.
      k=ir
      if(k.gt.nrused)then
         write(*,*)'esforce radial node number too large. Reset.'
         k=nrused
      endif
      delr=(r(2)-r(1))

      
      if(k.eq.1)then
         rkp=0.5*(rcc(k)+rcc(k+1))
         rkp2=0.5*(rcc(k+1)+rcc(k+2))
      elseif(k.eq.nrused)then
         rkm=0.5*(rcc(k)+rcc(k-1))
         rkm2=0.5*(rcc(k-1)+rcc(k-2))
      else
         rkp=0.5*(rcc(k)+rcc(k+1))
         rkm=0.5*(rcc(k)+rcc(k-1))
      endif

      qp=0.
      epz=0.
      epx=0.
      epy=0.
      fz=0.
      fx=0
      fy=0.

      do l=1,npsiused
         sp=sin(pcc(l))
         cp=cos(pcc(l))
         spp=sin(pcc(l+1))
         cpp=cos(pcc(l+1))
         if(l.eq.npsiused) then
            spp=sin(pcc(1))
            cpp=cos(pcc(1))
         endif

         do j=1,nthused-1
c     The potential at psi=npsiused (and psi=0) already accounts for the
c     periodicity
            phihere1= 0.5*(phi(k,j,l)+phi(k,j,l+1))
            phihere2= 0.5*(phi(k,j+1,l)+phi(k,j+1,l+1))
            phihere=0.5*(phihere1+phihere2)
                
            if(k.eq.1)then
               phiherep1=0.5*(phi(k+1,j,l)+phi(k+1,j,l+1))
               phiherep2=0.5*(phi(k+1,j+1,l)+phi(k+1,j+1,l+1))
               phiherepp1=0.5*(phi(k+2,j,l)+phi(k+2,j,l+1))
               phiherepp2=0.5*(phi(k+2,j+1,l)+phi(k+2,j+1,l+1))

               er1=-(rkp*rkp*(phiherep1-phihere1)*(1+.5)- .5
     $              *rkp2*rkp2*(phiherepp1-phiherep1))/delr
               er2=-(rkp*rkp*(phiherep2-phihere2)*(1+.5)- .5
     $              *rkp2*rkp2*(phiherepp2-phiherep2))/delr
               
c Old version (assumes er=a+b*cos(theta))
c                er=0.5*(er1+er2)
c New version (assumes er=a+b*theta)
               er=(er1*thang(j+1)-er2*thang(j))/(thang(j+1)-thang(j))
     $              +(er1-er2)/((thang(j+1)-thang(j))*(th(j+1)-th(j)))
     $              *(sqrt(1-th(j+1)**2)-sqrt(1-th(j)**2)-th(j+1)
     $              *thang(j+1)+th(j)*thang(j))
               
            elseif(k.eq.nrused)then
               phiherem1=0.5*(phi(k-1,j,l)+phi(k-1,j,l+1))
               phiherem2=0.5*(phi(k-1,j+1,l)+phi(k-1,j+1,l+1))
               phiheremm1=0.5*(phi(k-2,j,l)+phi(k-2,j,l+1))
               phiheremm2=0.5*(phi(k-2,j+1,l)+phi(k-2,j+1,l+1))

               er1=-(rkm*rkm*(phihere1-phiherem1)*(1+.5)-
     $              .5*rkm2*rkm2*(phiherem1-phiheremm1))/delr
               er2=-(rkm*rkm*(phihere2-phiherem2)*(1+.5)-
     $              .5*rkm2*rkm2*(phiherem2-phiheremm2))/delr

c Old version (assumes er=a+b*cos(theta))
c                er=0.5*(er1+er2)
c New version (assumes er=a+b*theta)
               er=(er1*thang(j+1)-er2*thang(j))/(thang(j+1)-thang(j))
     $              +(er1-er2)/((thang(j+1)-thang(j))*(th(j+1)-th(j)))
     $              *(sqrt(1-th(j+1)**2)-sqrt(1-th(j)**2)-th(j+1)
     $              *thang(j+1)+th(j)*thang(j))

            else
               phiherem= 0.25*(phi(k-1,j,l)+phi(k-1,j+1,l)
     $              +phi(k-1,j,l+1) +phi(k-1,j+1,l+1)) 
               phiherep= 0.25*(phi(k+1,j,l)+phi(k+1,j+1,l)
     $              +phi(k+1,j,l+1) +phi(k+1,j+1,l+1)) 
               er=-0.5*(rkp*rkp*(phiherep-phihere)+
     $              rkm*rkm*(phihere-phiherem))/delr
            endif
            
            epz=epz+2*ercoefZ(j)*exp(phihere)*dpsi
            epx=epx+2*ercoefX(j)*exp(phihere)*(spp-sp)
            epy=epy+2*ercoefX(j)*exp(phihere)*(cp-cpp)

            qp=qp+qpcoef(j)*er*dpsi

c     E=-gradPhi
            eth=-0.5*(phi(k,j+1,l)+phi(k,j+1,l+1)-phi(k,j,l)-phi(k,j,l+1
     $           ))

         epsi1=-(phi(k,j,l+1)-phi(k,j,l))/(sqrt(1-tcc(j)**2)+1e-7)
         epsi2=-(phi(k,j+1,l+1)-phi(k,j+1,l))/(sqrt(1-tcc(j+1)**2)+1e-7)
c     Assume (assumes epsi=a+b*theta) epsi is actually dphi/dpsi/sin(theta)
c            epsi=(epsi1*thang(j+1)-epsi2*thang(j))/(thang(j+1)-thang(j))
c     $           +(epsi1-epsi2)/((thang(j+1)-thang(j))*(th(j+1)-th(j)))
c     $           *(sqrt(1-th(j+1)**2)-sqrt(1-th(j)**2)-th(j+1) *thang(j
c     $           +1)+th(j)*thang(j))

c     Assumes epsi=a+b*cos(theta)). Seems to work better that way
            epsi=0.5*(epsi1+epsi2)


            fz=fz+ (ercoefZ(j)*er**2/rcc(k)**2
     $           +etcoefZ(j)*eth*eth
     $           +epcoefZ(j)*epsi*epsi
     $           +ertcoefZ(j)*er*eth/rcc(k) )*dpsi
            
            fx=fx+(ercoefX(j)*er**2/rcc(k)**2
     $           +etcoefX(j)*eth*eth
     $           +epcoefX(j)*epsi*epsi
     $           +ertcoefX(j)*er*eth/rcc(k) )*(spp-sp)
     $           +erpcoefX(j)*er*epsi/rcc(k)*(cpp-cp)
            
            fy=fy+(ercoefX(j)*er**2/rcc(k)**2
     $           +etcoefX(j)*eth*eth
     $           +epcoefX(j)*epsi*epsi
     $           +ertcoefX(j)*er*eth/rcc(k) )*(cp-cpp)
     $           +erpcoefX(j)*er*epsi/rcc(k)*(spp-sp)

c            if(ir.eq.nrused)then
c               write(*,*) j,ercoefX(j),er
c            endif

         enddo
      enddo


      qp=qp
      epz=-epz*r(k)**2
      epx=-epx*r(k)**2
      epy=-epy*r(k)**2



c     Calculate the Lorentz force on the ions inside the control surface
c     in a frame where the plasma at infinity is at rest     
      if(ir.ne.1) then

c Old way
c            vx=0.
c            vy=0.
c            vz=0.
c            partsum=0.
c            sd=sqrt(1-cd**2)
c            sb=sqrt(1-cb**2)
            
c            do k=1,npsiused
c               do j=1,nthused
c                  do i=1,ir-1
c     Old way, not accurate enough
c     vx=vx+(vrsum(i,j,k)*sqrt(1-tcc(j)**2)+vtsum(i,j,k)
c     $                 *tcc(j))*cos(pcc(k))-vpsum(i,j,k)*sin(pcc(k))
c     vy=vy+(vrsum(i,j,k)*sqrt(1-tcc(j)**2)+vtsum(i,j,k)
c     $                 *tcc(j))*sin(pcc(k))+vpsum(i,j,k)*cos(pcc(k))
c     vz=vz+vrsum(i,j,k)*tcc(j)
c     $                 -vtsum(i,j,k)*sqrt(1-tcc(j)**2)             
c                     vx=vx+vxsum(i,j,k)
c                     vy=vy+vysum(i,j,k)
c                     vz=vz+vzsum(i,j,k)
c                     partsum=partsum+psum(i,j,k)
c                  enddo
c               enddo
c            enddo
c            i=ir
c            if(i.ne.nrused) then
c     We only sum the inner half of the last radial cell.  frac\sim 0.5
c     is the volume ratio of the first half of a cell (radial direction)
c     over the full cell (first order in delr/rcc(i))
c               frac=(rcc(i)-0.5*delr)/(2*rcc(i))
c            else
c     if ir.eq.nrused, the cell is already half the size
c               frac=1.
c            endif
c            do k=1,npsiused
c               do j=1,nthused
c                  vx=vx+frac*vxsum(i,j,k)
c                  vy=vy+frac*vysum(i,j,k)
c                  vz=vz+frac*vzsum(i,j,k)
c                  partsum=partsum+frac*psum(i,j,k)
c               enddo
c            enddo
            
c         else
         vx=curr(1)
         vy=curr(2)
         vz=curr(3)
         partsum=curr(4)
c         endif

         sd=sqrt(1-cd**2)
         sb=sqrt(1-cb**2)
         

         vx=vx
         vy=(vy-partsum*vd*sd)
         vz=(vz-partsum*vd*cd)

         
         fbx=Bz*(vy*cb-vz*sb)
         fby=-Bz*cb*vx
         fbz=Bz*sb*vx

      else
         fbx=0.
         fby=0.
         fbz=0.
      endif

      end
