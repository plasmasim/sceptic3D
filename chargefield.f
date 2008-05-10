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


      if(diags.or.savelor) then
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
      if(diags) then
         do k=1,npsiused
            do j=1,nthused
               do i=1,nrused
                  vrsum(i,j,k)=0.
                  vtsum(i,j,k)=0.
                  vpsum(i,j,k)=0.
                  v2sum(i,j,k)=0.
                  vr2sum(i,j,k)=0.
                  vtp2sum(i,j,k)=0.
               enddo
            enddo
         enddo
      elseif(debyelen.le.0.01) then
         do j=1,nthused
            do k=1,npsiused
               do i=1,2
                  vrsum(i,j,k)=0.
                  vtsum(i,j,k)=0.
                  vpsum(i,j,k)=0.
                  v2sum(i,j,k)=0.
                  vr2sum(i,j,k)=0.
                  vtp2sum(i,j,k)=0.
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
      end
c***********************************************************************
c Accumulate particle charge into rho mesh and other diagnostics.
      subroutine chargeassign(i,irl,rf,ithl,thf,ipl,pf,st,ct,sp,cp,rp)
c      implicit none
      integer i
c Common data:
      include 'piccom.f'

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

      if(diags.or.(debyelen.le.0.01.and.irl.le.2)) then

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

      v2=(xp(4,i)*xp(4,i) +xp(5,i)*xp(5,i) +xp(6,i)*xp(6,i))
      v2sum(irl,ithl,ipl)=v2sum(irl,ithl,ipl) + 
     $     (1.-rf)*(1.-thf)*(1-pf)*v2
      v2sum(irl+1,ithl,ipl)=v2sum(irl+1,ithl,ipl) + 
     $     rf*(1.-thf)*(1-pf)*v2
      v2sum(irl,ithl+1,ipl)=v2sum(irl,ithl+1,ipl) + 
     $     (1.-rf)*thf*(1-pf)*v2
      v2sum(irl+1,ithl+1,ipl)=v2sum(irl+1,ithl+1,ipl) + 
     $     rf*thf*(1-pf)*v2
      v2sum(irl,ithl,iplp)=v2sum(irl,ithl,iplp) + 
     $     (1.-rf)*(1.-thf)*pf*v2
      v2sum(irl+1,ithl,iplp)=v2sum(irl+1,ithl,iplp) + 
     $     rf*(1.-thf)*pf*v2
      v2sum(irl,ithl+1,iplp)=v2sum(irl,ithl+1,iplp) + 
     $     (1.-rf)*thf*pf*v2
      v2sum(irl+1,ithl+1,iplp)=v2sum(irl+1,ithl+1,iplp) + 
     $     rf*thf*pf*v2

      vtp2=v2-vr2
      vtp2sum(irl,ithl,ipl)=vtp2sum(irl,ithl,ipl) + 
     $     (1.-rf)*(1.-thf)*(1-pf)*vtp2
      vtp2sum(irl+1,ithl,ipl)=vtp2sum(irl+1,ithl,ipl) + 
     $     rf*(1.-thf)*(1-pf)*vtp2
      vtp2sum(irl,ithl+1,ipl)=vtp2sum(irl,ithl+1,ipl) + 
     $     (1.-rf)*thf*(1-pf)*vtp2
      vtp2sum(irl+1,ithl+1,ipl)=vtp2sum(irl+1,ithl+1,ipl) + 
     $     rf*thf*(1-pf)*vtp2
      vtp2sum(irl,ithl,iplp)=vtp2sum(irl,ithl,iplp) + 
     $     (1.-rf)*(1.-thf)*pf*vtp2
      vtp2sum(irl+1,ithl,iplp)=vtp2sum(irl+1,ithl,iplp) + 
     $     rf*(1.-thf)*pf*vtp2
      vtp2sum(irl,ithl+1,iplp)=vtp2sum(irl,ithl+1,iplp) + 
     $     (1.-rf)*thf*pf*vtp2
      vtp2sum(irl+1,ithl+1,iplp)=vtp2sum(irl+1,ithl+1,iplp) + 
     $     rf*thf*pf*vtp2
      endif

      if(diags.or.savelor) then

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
      real phi0mphi1(0:nthsize,0:npsisize),delphi0(0:nthsize,0:npsisize)
      parameter (nthp1=nthsize+1)
      parameter (npsip1=npsisize+1)
      parameter (ntot=nthp1*npsip1)
c      real phi1ave
      real bcifac,bcpfac,bci,bcp,bvf
      real relax
      real cs(nthsize,npsisize),csd(nthsize,npsisize)
      real ncs
      logical first
      integer kk1,kk2
c     stencil is the excursion around a given cell for the sound speed
c     averaging
      integer stencil
      parameter (stencil=1)
      data relax/1./
      data bcifac/.2/bcpfac/.1/
      data bvf/1.2071/
      data first/.true./
      data phi0mphi1/ntot*0./
      data delphi0/ntot*0./
      data ncs/50./
      save
      cerr=0.

      do k=1,npsiused
         do j=1,nthused
            do i=1,nrused
               if(rho(i,j,k).le.0.)then
                  write(*,*)'rho=0',i,j,k
                  stop
               endif
c Simplistic Boltzmann scheme. May need relaxation.
               delta=phi(i,j,k)-log(rho(i,j,k))
               if(abs(delta).le.cerr)cerr=abs(delta)
               phi(i,j,k)=phi(i,j,k)-relax*delta
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
c      write(*,*)'p1,p2,v1,v2,csd,cs,vs,phi0,phi1,delphinew'
c      write(*,*)'cs=',(cs(kk),kk=1,nth)
      do j=1,nthused
         do k=1,npsiused
            p1=0.
            v1=0.
            s1=0.
            p2=0.
            v2=0.
            s2=0.
            do l=-stencil,stencil
               ja=j+l
               if(ja.eq.0.or.ja.eq.nthused+1) ja=j
               do m=-stencil,stencil
                  ka=k+m
                  if(ka.eq.0) then
                     ka=npsiused
                  elseif(ka.eq.npsiused+1) then
                     ka=1
                  endif
                  p1=p1+vr2sum(1,ja,ka)*psum(1,ja,ka)-vrsum(1,ja,ka)**2
                  v1=v1+vrsum(1,ja,ka)
                  s1=s1+psum(1,ja,ka)
                  p2=p2+vr2sum(2,ja,ka)*psum(2,ja,ka)-vrsum(2,ja,ka)**2
                  v2=v2+vrsum(2,ja,ka)
                  s2=s2+psum(2,ja,ka)
               enddo
            enddo

            if(p1.ne.0)then
               if(s2.gt.0)then
                  p1=p1/s1
                  v1=v1/s1
               else
                  write(*,*)'s1=0'
                  stop
               endif
            endif
            if(p2.ne.0)then
               if(s2.gt.0)then
                  p2=p2/s2
                  v2=v2/s2
               else
                  write(*,*)'s2=0'
                  stop
               endif
            endif
            vs=(1.+bvf)*v1-bvf*v2
c Fix the psum difference so it can't be zero.
c Original derivative form
            csd(j,k)=(p2 - p1)/(s2-s1+.5)+1.

c Gamma Ti form with gamma=3.
c         csd(j)=3.*p1/(psum(1,j)+.5) + 1.
            if(csd(j,k).lt.0.)csd(j,k)=0.
            csd(j,k)=-sqrt(csd(j,k))
c Clip the excursion symmetrically.
            if(csd(j,k).lt.2.*cs(j,k))  csd(j,k)=2.*cs(j,k)
c Average the sound-speed value over ncs steps.
            cs(j,k)=((ncs-1.)*cs(j,k)+csd(j,k))/ncs
            if(cs(j,k).gt.0)then
               write(*,*)'cs positive',j,cs(j,k),csd(j,k)
     $              ,p1,p2,psum(1,j,k),psum(2,j,k),ncs
               stop
            endif
            bci=-bcifac*cs(j,k)**2*dt/(r(2)-r(1))
            bcp=-bcpfac*(r(2)-r(1))/(dt*cs(j,k))
            delphinew=(vs - cs(j,k))*bci
            phi(0,j,k)=phi(1,j,k)+phi0mphi1(j,k)+delphinew+
     $           (delphinew-delphi0(j,k))*bcp
            delphi0(j,k)=delphinew
            if(phi(0,j,k).gt.phi(1,j,k))phi(0,j,k)=phi(1,j,k)
            phi0mphi1(j,k)=phi(0,j,k)-phi(1,j,k)
c         write(*,'(10f8.3)')p1,p2,v1,v2,csd(j),cs(j),vs,
c     $        phi(0,j),phi(1,j),delphinew
c Adjusting the potential of the first cell.

            phi(1,j,k)=phi(0,j,k)
            if(.not.abs(phi(1,j,k)).lt.1.e20)then
               write(*,*)'phi1 overflow', phi(1,j,k),bcp,cs(j,k),ncs
     $              ,delphinew ,p1,p2,psum(1,j,k),psum(2,j,k),csd(j,k)
               stop
            endif
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
c Radial, Theta and Phi accelerations
      real ar,at,ap
      real ct,st,cp,sp,rp
      real zetap,hf
      integer ih
      real dp,dth,dpinv,dthinv
      integer ipl,ith
      integer il,ir,ithp1,ithp2,ithm1,ilm1
      integer iplp1,iplp2,iplm1
      real rlm1,rf,tf,pf,rr,rl,dpsi

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

c     Theta weighting
      phih1tX=phih1pt*tf +phih1tt*(1-tf)
      phih1pX=phih1pp*tf +phih1tp*(1-tf)
      phih1mX=phih1pm*tf +phih1tm*(1-tf)
      phih12X=phih1p2*tf +phih1t2*(1-tf)


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
      phihp1tX=phihp1pt*tf +phihp1tt*(1-tf)
      phihp1pX=phihp1pp*tf +phihp1tp*(1-tf)
      phihp1mX=phihp1pm*tf +phihp1tm*(1-tf)
      phihp12X=phihp1p2*tf +phihp1t2*(1-tf)




c------------------------------
c Deal with radial acceleration
c -----------------------------
      
c Here we control whether we use the zeta or r interpolation.
      if(debyelen.lt.1.e-2)then
c Linear approx to sqrt form at boundary.
         if(ih.eq.1)then
c            rr=r(ir)
            philm1tt=phi(il,ith,ipl)-bdyfc*sqrt(2.*(rr-rl))*0.25
            philm1pt=phi(il,ithp1,ipl)-bdyfc*sqrt(2.*(rr-rl))*0.25
            philm1tp=phi(il,ith,iplp1)-bdyfc*sqrt(2.*(rr-rl))*0.25
            philm1pp=phi(il,ithp1,iplp1)-bdyfc*sqrt(2.*(rr-rl))*0.25
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
     $           )*(1.-tf)
     $        +( (phihp1p-phih1p)/(zeta(ir)-zeta(ih))*hf +
     $        (phih1p-philm1p)/(zeta(ih)-zeta(ilm1))*(1.-hf)
     $         )*tf )/zetap
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
     $           )*(1.-tf)
     $        +( (phihp1p-phih1p)/(rr-rl)*hf +
     $        (phih1p-philm1p)/(rl-rlm1)*(1.-hf)
     $         )*tf )

      endif

c The minus sign is because E=-grad Phi         
      ar=-ar


c------------------------------
c Deal with theta acceleration
c -----------------------------

c     In sceptic2D, we always divide by rl. It think it is wrong. One
c     should divide by rl or rr depending on the radius at which the
c     potential derivative is calculated

      if(tf.le.0.5)then
         at= ( (phih1p-phih1t)*(tf)*2.
     $          /(rl*(thang(ithp1)-thang(ith)))
     $        +(phih1p-phih1m)*(0.5-tf)*2.
     $          /(rl*(thang(ithp1)-thang(ithm1))) ) * (1.-rf)
     $        + ( (phihp1p-phihp1t)*(tf)*2.
     $          /(rr*(thang(ithp1)-thang(ith)))
     $        +(phihp1p-phihp1m)*(0.5-tf)*2.
     $          /(rr*(thang(ithp1)-thang(ithm1))) ) * rf
      else
         at= ( (phih12-phih1t)*(tf-0.5)*2.
     $          /(rl*(thang(ithp2)-thang(ith)))
     $        +(phih1p-phih1t)*(1.-tf)*2.
     $          /(rl*(thang(ithp1)-thang(ith))) ) * (1.-rf)
     $        + ( (phihp12-phihp1t)*(tf-0.5)*2.
     $          /(rr*(thang(ithp2)-thang(ith)))
     $        +(phihp1p-phihp1t)*(1.-tf)*2.
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
      ap=-ap/(dpsi*st+1e-6)
      if(lap0)ap=0.

 501  format(a,6f10.4)


c 3D acceleration
      accel(3)=ar*ct - at*st
      accel(2)=(ar*st+ at*ct)*sp+ap*cp
      accel(1)=(ar*st+ at*ct)*cp-ap*sp
c Trap errors.


      if(.not.accel(1).lt.1.e5)then
         write(*,*) 'i: ',i,' x: ',xp(1,i),' y: ',xp(2,i),' z: ',xp(3,i)
         write(*,*)'Accel Excessive: ar,at,ap,ct,sp,cp,'
         write(*,*) ar,at,ap,ct,sp,cp
         write(*,*) 'phi at ith and ithp1:'
         write(*,*) philm1t,phi(ih,ith,ipl),phihp1t
         write(*,*) philm1p,phi(ih,ithp1,ipl),phihp1p
         write(*,*) 'zetap=',zetap,'  bdyfc=',bdyfc
         write(*,*) 'zeta '
         write(*,*) zeta(ilm1),zeta(ih),zeta(ir)
         write(*,'(10f7.4)')((phi(i,ih,ipl),i=1,10),ih=1,10)
         stop
      endif
      end
c***********************************************************************

c**********************************************************************
      subroutine esforce(ir,qp,fz,epz,fbz,fx,epx,fbx,fy,epy,fby)
      include 'piccom.f'
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
c sin theta at cell boundary
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

            shalf=sin(acos(0.5*(th(j+1)+th(j))))
c     The 3D force calculation uses a cariant of ercoef with respect to
c     the 2D routine.
            ercoefZ(j)=-0.25*(th(j+1)**2-th(j)**2)
            
            etcoefZ(j)=0.5*(th(j+1)**2/2.-th(j+1)**4/4. -th(j)**2/2.
     $           +th(j)**4/4.)/(th(j+1)-th(j))**2
            
            epcoefZ(j)=-ercoefZ(j)/shalf**2/dpsi**2
            
            ertcoefZ(j)=-(1.-(th(j+1)*(th(j+1)+th(j))+th(j)**2)/3.)
            
            ercoefX(j)=0.25*((thang(j+1)-th(j+1)*sqrt(1-th(j+1)**2)
     $           )-(thang(j)-th(j)*sqrt(1-th(j)**2)))
            
            etcoefX(j)=-0.125*((3.*thang(j+1)/2-th(j+1)*sqrt(1-th(j+1
     $           )**2)*(3./2+1.-th(j+1)**2))- (3.*thang(j)/2-th(j)
     $           *sqrt(1-th(j)**2)*(3./2+1.-th(j)**2)))/(th(j+1)-th(j))
     $           **2
            
            epcoefX(j)=-ercoefX(j)/shalf**2/dpsi**2
            
            ertcoefX(j)=-1/(3.*(th(j+1)-th(j))) *((1-th(j+1)**2)**(3./2)
     $           -(1-th(j)**2)**(3./2))
            
            erpcoefX(j)=-1./dpsi*(th(j+1)-th(j))/shalf


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
         
c         write(*,*)'j, qpcoef,  ercoefZ,  etcoefZ,  epcoefZ, ertcoefZ'
c         write(*,'(i3,5f10.5)')(j,qpcoef(j),ercoefZ(j),etcoefZ(j)
c     $        ,epcoefZ(j), ertcoefZ(j),j=1,nthused)

c         write(*,*) ""

c         write(*,*)'j, ercoefX,  etcoefX,  epcoefX,  ertcoefX, erpcoefX'
c         write(*,'(i3,5f10.5)')(j,ercoefX(j),etcoefX(j)
c     $        ,epcoefX(j), ertcoefX(j),erpcoefX(j),j=1,nthused)

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
            phihere= 0.25*(phi(k,j,l)+phi(k,j+1,l)+phi(k,j,l+1)+phi(k,j
     $           +1,l+1))
         
            if(k.eq.1)then
               phiherep= 0.25*(phi(k+1,j,l)+phi(k+1,j+1,l)
     $              +phi(k+1,j,l+1) +phi(k+1,j+1,l+1)) 
               phiherepp= 0.25*(phi(k+2,j,l)+phi(k+2,j+1,l)
     $              +phi(k+2,j,l+1) +phi(k+2,j+1,l+1)) 
               er=-(rkp*rkp*(phiherep-phihere)*(1+.5)- .5
     $              *rkp2*rkp2*(phiherepp-phiherep))/delr
            elseif(k.eq.nrused)then
               phiherem= 0.25*(phi(k-1,j,l)+phi(k-1,j+1,l)
     $              +phi(k-1,j,l+1) +phi(k-1,j+1,l+1)) 
               phiheremm= 0.25*(phi(k-2,j,l)+phi(k-2,j+1,l)
     $              +phi(k-2,j,l+1) +phi(k-2,j+1,l+1)) 
               er=-(rkm*rkm*(phihere-phiherem)*(1+.5)-
     $              .5*rkm2*rkm2*(phiherem-phiheremm))/delr
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
            epsi=-0.5*(phi(k,j,l+1)+phi(k,j+1,l+1)-phi(k,j,l)-phi(k,j+1
     $           ,l))

            
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
            
            
c     write(*,*)j,erp,eth,ercoef(j+1)*erp*erp,etcoef(j)*eth*eth,
c     $           - ertcoef(j)*0.5*(erp+er)*eth
c     write(*,*)fz
         enddo
      enddo
      qp=qp
      epz=-epz*r(k)**2
      epx=-epx*r(k)**2
      epy=-epy*r(k)**2

c     Calculate the Lorentz force on the ions inside the control surface
c     in a frame where the plasma at infinity is at rest

      
      if(ir.ge.2.and.savelor) then
         vx=0.
         vy=0.
         vz=0.
         partsum=0.
         sd=sqrt(1-cd**2)
         sb=sqrt(1-cb**2)

         do k=1,npsiused
            do j=1,nthused
               do i=1,ir-1
c Old way, not accurate enough
c                  vx=vx+(vrsum(i,j,k)*sqrt(1-tcc(j)**2)+vtsum(i,j,k)
c     $                 *tcc(j))*cos(pcc(k))-vpsum(i,j,k)*sin(pcc(k))
c                  vy=vy+(vrsum(i,j,k)*sqrt(1-tcc(j)**2)+vtsum(i,j,k)
c     $                 *tcc(j))*sin(pcc(k))+vpsum(i,j,k)*cos(pcc(k))
c                  vz=vz+vrsum(i,j,k)*tcc(j)
c     $                 -vtsum(i,j,k)*sqrt(1-tcc(j)**2)
                  
                  vx=vx+vxsum(i,j,k)
                  vy=vy+vysum(i,j,k)
                  vz=vz+vzsum(i,j,k)
                  partsum=partsum+psum(i,j,k)
               enddo
            enddo
         enddo
         i=ir
         if(i.ne.nrused) then
c     We only sum the inner half of the last radial cell.  frac\sim 0.5
c     is the volume ratio of the first half of a cell (radial direction)
c     over the full cell (first order in delr/rcc(i))
            frac=(rcc(i)-0.5*delr)/(2*rcc(i))
         else
c if ir.eq.nrused, the cell is already half the size
            frac=1.
         endif
         do k=1,npsiused
            do j=1,nthused
               vx=vx+frac*vxsum(i,j,k)
               vy=vy+frac*vysum(i,j,k)
               vz=vz+frac*vzsum(i,j,k)
               partsum=partsum+frac*psum(i,j,k)
            enddo
         enddo
         
               
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
