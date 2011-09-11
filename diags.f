
      subroutine chargediag(dt,istep,icolntype)
      real dt
      integer istep,icolntype
c Common data:
      include 'piccom.f'
      include 'errcom.f'
      include 'fvcom.f'
      real rhoplot(nrsize),rho1theta(nthsize),rhomidtheta(nthsize)
      real rhomidave(nthsize),rho1ave(nthsize)
      real phiave(nrsize)
      real phitemp
c      real riave
      save
c Calculate rhoplot,diagphi,diagrho,rho1theta,rhomidtheta
      do i=1,nr
 510     format(10f8.1)
         rhoplot(i)=0.
         phiave(i)=0.
         nrp=0
      enddo


c Assume the unperturbed region is upstream outside the magnetic shadow      
      do k=1,NPSIUSED
         do j=1,NTHUSED

c     If (very, i.e. in the first 0.25% of the domain) upstream (the
c     "*vd" at the end is just to capture vd's sign)
            if((cd*tcc(j)+sqrt(1-cd**2)*sqrt(1-tcc(j)**2)*sin(pcc(k)))
     $           *(vd+1e-7)/(abs(vd)+1e-7).le.-0.5) then

c      If outside the magnetic shadow (the term (0.5)**2 means that we
c      consider the magnetic shadow to be rcc(nrused)/2)
               if(Bz.eq.0.or.((sqrt(1-cB**2) *tcc(j)-cB *sqrt(1-tcc(j)
     $              **2)*sin(pcc(k)))**2+(-sqrt(1-tcc(j)**2) *cos(pcc(k)
     $              )*sqrt(1 -cB**2))**2+(sqrt(1-tcc(j)**2) *cos(pcc(k))
     $              *cB)**2.ge.(0.5)**2)) then
                  nrp=nrp+1
                  do i=1,nr
c     rhoplot is unnormalized. All others are normalized.
                     rhoplot(i)=rhoplot(i)+rho(i,j,k)
                     phiave(i)=phiave(i)+phi(i,j,k)
                  enddo
               endif
            endif
         enddo
      enddo
      do i=1,nr
         rhoplot(i)=rhoplot(i)/nrp
         phiave(i)=phiave(i)/nrp
         diagphi(i)=(diagphi(i)*(nstepsave-1)+
     $        phiave(i))/nstepsave
         diagrho(i)=(diagrho(i)*(nstepsave-1) + (rhoplot(i)))/nstepsave
      enddo

c Calculate diagchi (outer potential as a function of nth normalized
c to the ion thermal velocity)
c Necessary for the fortran diagnostics
      phiout=0.
      do j=1,NTHUSED
         phitemp=0.
         do k=1,NPSIUSED
            phitemp=phitemp+phi(NRUSED,j,k)
         enddo
         diagchi(j)=(diagchi(j)*(nstepsave-1)+phitemp/Ti/NPSIUSED)
     $        /nstepsave
         phiout=phiout+phitemp
      enddo
      phiout=phiout/NTHUSED

c  Calcualte the averaged data on collected particles
      do k=1,npsiused
         do j=1,nthused
            delta=fincellave(j,k)-nincellstep(j,k,istep-1)
            fincellave(j,k)=fincellave(j,k)-delta/nstepsave
            delta=vrincellave(j,k)-vrincellstep(j,k,istep-1)
            vrincellave(j,k)=vrincellave(j,k)-delta/nstepsave
            delta=vr2incellave(j,k)-vr2incellstep(j,k,istep-1)
            vr2incellave(j,k)=vr2incellave(j,k)-delta/nstepsave
         enddo
      enddo


c Calculate rhoinf if it is not being done in main.
c      call rhoinfcalc(dt,icolntype,colnwt)

      end

c********************************************************************
      real function smaxflux(uc,chi)
c     Return the total flux to a unit radius sphere from a unit density
c     maxwellian distribution shifted by velocity
      real uc
c     normalized to sqrt(2T/m), in a spherically symmetric potential
c     having a value on the sphere normalized to Ti of minus
      real chi

      real eps,pi
      data eps/1.e-3/pi/3.1415927/

      erf=1.-erfcc(uc)
      sqpi=sqrt(pi)
      if(abs(uc).lt.eps) then
         erfbyu=(2./sqpi)*(1.-uc**2 /3.)
      else
         erfbyu=erf/uc
      endif

      smaxflux=pi*sqrt(2.)*(uc*erf +(0.5+chi)*erfbyu + exp(-uc**2)/sqpi)
      end
c*******************************************************************
      real function smaxflux2(chi)
c     Returns the total flux to a unit radius sphere from a unit density
c     Maxwellian accelerated adiabatiquely in the z direction.

      real flux,chi
      real pi
      data pi/3.1415927/
      sqpi=sqrt(pi)

      flux=1-sqpi/2*sqrt(chi)*exp(chi)*erfcc(sqrt(chi))
c      write(*,*) "flux",flux
c    Seems not to converge if the potential is not prescribed in advance


      smaxflux2=sqpi*2*flux
      end

c***********************************************************************
c Contouring of the charge density, rho, on distorted mesh., averaged over psi
      subroutine rhodisplay(ir,rhomin,rhomax)
      integer ir
c Common data:
      include 'piccom.f'
      include 'errcom.f'
c      save
      character cworka(nrsize*nthsize+13)
      integer ncont
      parameter (ncont=9)
      real zclv(ncont)
      real xrho(nrsize+1,0:nthsize),zrho(nrsize+1,0:nthsize)
      real rhoave(0:nrsize,0:nthsize)
      logical lfirst
      data lfirst/.true./
      save xrho,zrho,lfirst,first
      
      rmi=rhomin
      rma=rhomax
      do i=1,NRUSED
         do j=1,NTHUSED
            zrho(i,j)=rcc(i)*tcc(j)
            xrho(i,j)=rcc(i)*sqrt(1.-tcc(j)**2)
            rhoave(i,j)=0.
            do k=1,NPSIUSED
               rhoave(i,j)=rhoave(i,j)+rho(i,j,k)
            enddo
c     rhoave is the density averaged in the psi direction
            rhoave(i,j)=rhoave(i,j)/NPSIUSED
         enddo
         if(.not.LCIC)then
            zrho(i,0)=rcc(i)
            xrho(i,0)=0.
            zrho(i,nth)=-rcc(i)
            xrho(i,nth)=0.
         endif
      enddo

      if(rma.eq.0)then
         call minmax(rhoave(1,1),NRUSED,rmi,rma)
         call minmax(rhoave(1,NTHUSED),NRUSED,rn2,rx2)
         rmi=min(rmi,rn2)
         rma=max(rma,rx2)
         if(rma.lt.1.2)rma=1.
         if(rma.gt.5.)then
            rma=5.
            rmi=0.
         endif
         if(Ezext.ne.0. .and. rma.gt.2.)then
            rma=2.
            rmi=0.
         endif
         if(rmi.lt.0.7)rmi=0.
      endif
      call fitrange(rmi,rma,ncont-1,ipow,fac10,delta,first,xlast)
c      write(*,*)'rhomin/max=',rmi,rma,first,xlast,delta
      ncl=0
      do j=1,ncont
         zclv(j)=(j-1)*delta + first
         if(zclv(j).gt. xlast*1.0001) goto 440
         ncl=ncl+1
      enddo
 440  continue
      icl=-ncl
      call pltinaspect(-r(nr),r(nr),0.,r(nr))
      if(lfirst)then
         call accisgradinit(-50000,10000,-10000,130000,65000,65000)
c         call accisgradinit(-25000,00000,25000,130000,65000,130000)
         lfirst=.false.
      endif
c     There are signs that the following is overflowing bounds somehow.
      ntype=2+16+32
c      call contourl(rho(1,1),cworka,nr,nr,nth,zclv,icl,xrho,zrho,ntype)
      if(LCIC)then
         call contourl(rhoave(1,1),cworka,nrsize+1,NRUSED,nth,
     $        zclv,icl,zrho(1,1),xrho(1,1),ntype)
      else
         call contourl(rhoave(1,0),cworka,nrsize+1,NRUSED,nth+1,
     $        zclv,icl,zrho,xrho,ntype)
      endif
c      call color(15)
      call gradlegend(zclv(1),zclv(ncl),
     $     .1,1.25,.9,1.25,0.02,.true.)
      call axis()
      call axlabels('z     density color, velocity arrows','r sin!Aq!@')
c      call legendline(.1,-.3,258,'density color-contours')
      call color(12)
      if(ir.le.0.or.ir.ge.100) ir=10
      j0=1
      if(LCIC)j0=2
      do j=j0,NTHUSED-1,NTHUSED/5+1
         do i=1,NRUSED,max(nr/ir,1)
            vri=0.
            vti=0.
            do k=1,NPSIUSED
               vri=vri+vrsum(i,j,k)/(psum(i,j,k)+1.e-5)
               vti=vti+vtsum(i,j,k)/(psum(i,j,k)+1.e-5)
            enddo
            vri=vri/NPSIUSED
            vti=vti/NPSIUSED
            size=0.02*sqrt(vri**2+vti**2)
            angle=atan2(vti,vri)+acos(tcc(j))
            call charsize(size,0.5*size)
            call charangl(180.*angle/3.141593)
            call jdrwstr(wx2nx(zrho(i,j)),wy2ny(xrho(i,j)),
     $           char(ichar('_')+128)//char(0),0.)
         enddo
      enddo
      call charangl(0.)
      call charsize(0.,0.)
      call color(15)
c      call pltend()
      end

c***********************************************************************
c Contouring of the charge density at theta=pi/2 (r,psi)
      subroutine rhodisplayPsi(ir,rhomin,rhomax)
      integer ir
c Common data:
      include 'piccom.f'
      include 'errcom.f'
c      save
      character cworka(nrsize*nthsize+13)
      integer ncont
      parameter (ncont=9)
      real zclv(ncont)
      real xrho(nrsize+1,0:npsisize),yrho(nrsize+1,0:npsisize)
      integer nthhalf
      real rhoave(0:nrsize,0:npsisize)
      logical lfirst
      data lfirst/.true./
      save xrho,yrho,lfirst,first
      
      nthhalf=nint(0.5*nthused)
      rmi=rhomin
      rma=rhomax
      do i=1,NRUSED
         do j=1,NPSIUSED
            xrho(i,j)=rcc(i)*cos(pcc(j))
            yrho(i,j)=rcc(i)*sin(pcc(j))
            rhoave(i,j)=rho(i,nthhalf,j)
         enddo
         xrho(i,NPSIFULL)=xrho(i,1)
         yrho(i,NPSIFULL)=yrho(i,1)
         rhoave(i,NPSIFULL)=rhoave(i,1)
      enddo

      if(rma.eq.0)then
         call minmax(rhoave(1,1),NRUSED,rmi,rma)
         call minmax(rhoave(1,nint(0.5*NPSIUSED)),NRUSED,rn2,rx2)
         rmi=min(rmi,rn2)
         rma=max(rma,rx2)
         if(rma.lt.1.2)rma=1.
         if(rma.gt.5.)then
            rma=5.
            rmi=0.
         endif
         if(Ezext.ne.0. .and. rma.gt.2.)then
            rma=2.
            rmi=0.
         endif
         if(rmi.lt.0.7)rmi=0.
      endif
      call fitrange(rmi,rma,ncont-1,ipow,fac10,delta,first,xlast)
c      write(*,*)'rhomin/max=',rmi,rma,first,xlast,delta
      ncl=0
      do j=1,ncont
         zclv(j)=(j-1)*delta + first
         if(zclv(j).gt. xlast*1.0001) goto 440
         ncl=ncl+1
      enddo
 440  continue
      icl=-ncl
      call pltinaspect(-r(nr),r(nr),-r(nr),r(nr))
      if(lfirst)then
         call accisgradinit(-50000,10000,-10000,130000,65000,65000)
c         call accisgradinit(-25000,00000,25000,130000,65000,130000)
         lfirst=.false.
      endif
c     There are signs that the following is overflowing bounds somehow.
      ntype=2+16+32
c      call contourl(rho(1,1),cworka,nr,nr,nth,zclv,icl,xrho,yrho,ntype)
      if(LCIC)then
         call contourl(rhoave(1,1),cworka,nrsize+1,NRUSED,npsifull,
     $        zclv,icl,xrho(1,1),yrho(1,1),ntype)
      else
         call contourl(rhoave(1,0),cworka,nrsize+1,NRUSED,npsifull+1,
     $        zclv,icl,xrho,yrho,ntype)
      endif
c      call color(15)
      call gradlegend(zclv(1),zclv(ncl),
     $     1.15,.1,1.15,.9,0.02,.false.)
      call axis()
      call axlabels('x density color, velocity arrows','y')
c      call legendline(.1,-.3,258,'density color-contours')
      call color(12)
      if(ir.le.0.or.ir.ge.100) ir=10
      j0=1
      do j=j0,NPSIUSED-1,NPSIUSED/10+1
         do i=1,NRUSED,max(nr/ir,1)
            vri=vrsum(i,nthhalf,j)/(psum(i,nthhalf,j)+1.e-5)
            vpi=vpsum(i,nthhalf,j)/(psum(i,nthhalf,j)+1.e-5)
            size=0.02*sqrt(vri**2+vpi**2)
            angle=atan2(vpi,vri)+pcc(j)
            call charsize(size,0.5*size)
            call charangl(180.*angle/3.141593)
            call jdrwstr(wx2nx(xrho(i,j)),wy2ny(yrho(i,j)),
     $           char(ichar('_')+128)//char(0),0.)
         enddo
      enddo
      call charangl(0.)
      call charsize(0.,0.)
      call color(15)
c      call pltend()
      end

c*******************************************************************

c Overplot orbits on existing plot.
      subroutine plotorbits
      include 'piccom.f'
      include 'errcom.f'
      call winset(.true.)
      do k=1,norbits
         call color(7)
         call polyline(zorbit(1,k),rorbit(1,k),iorbitlen(k))
         call color(15)
         call charsize(0.01,0.01)
         call accircle(wx2nx(zorbit(iorbitlen(k),k)),
     $        wy2ny(rorbit(iorbitlen(k),k)))
         call charsize(0.0,0.0)
      enddo
      call winset(.false.)
      end

c*********************************************************************
      subroutine slices(jstepth,rhomin,rhomax,time)
c Slices, averaged over psi

      include 'piccom.f'
      include 'errcom.f'
      character charin*100,ch2*100
      real rccleft(nrsize),phip1(nrsize)
      real rhoave(0:nrsize,0:nthsize),potave(0:nrsize,0:nthsize)
      logical sfirst
      data sfirst/.true./
      save rccleft,sfirst

      if(sfirst)then
         do j=1,NRUSED
            rccleft(j)=-rcc(j)
         enddo
         sfirst=.false.
      endif
      
      do i=1,NRUSED
         do j=1,NTHUSED
            rhoave(i,j)=0.
            potave(i,j)=0.
            do k=1,NPSIUSED
               rhoave(i,j)=rhoave(i,j)+rho(i,j,k)
               potave(i,j)=potave(i,j)+phi(i,j,k)
            enddo
            rhoave(i,j)=rhoave(i,j)/NPSIUSED
            potave(i,j)=potave(i,j)/NPSIUSED
c            if(rhoave(i,j).le.0) write(*,*)i,j,rhoave(i,j)
         enddo
      enddo

      jstepth=NTHUSED/2-1
      if(rhomax.eq.0.)then

c     I had to change the second argument to NRSIZE+1 instead of NRUSED
c     in the 2D version of SCEPTIC. I have no clue why
         call minmax2(rhoave(1,1),NRSIZE+1,NRUSED,NTHUSED,rhomin,rhopm)
         rhopm=max(1.2,.1*nint(10.*rhopm))
      else
         rhopm=rhomax
      endif
      if(rhopm.gt.4.)then
         phi0rho=rhopm*.9
      else
         phi0rho=1.2
      endif
      if(vprobe.ne.0.) then
         phiscale=phi0rho*abs(1./vprobe)
      else
         phiscale=phi0rho*0.25
      endif
      call pltinit(-rcc(NRUSED),rcc(NRUSED),
     $     -phi0rho/phiscale,(rhopm-phi0rho)/phiscale)
      call axptset(1.,1.)
      call jdrwstr(wx2nx(rcc(NRUSED)),wy2ny(-0.05/phiscale),
     $     '!Af!@',-1.5)
      call ticrev()
      call yaxis(0.,0.)
      call ticrev()
      call axptset(0.,0.)
      call scalewn(-rcc(NRUSED),rcc(NRUSED),0.,rhopm,
     $     .false.,.false.)
      call axis()
      call boxtitle('Density and Potential')
      call axlabels('Radial position','n/n!A!d;!d!@')
c     call polyline(rcc,rhoright,NRUSED)
c     call polyline(rccleft,rholeft,NRUSED)
  

      do j=1,NTHUSED/2,jstepth
         call color(mod(j,15)+1)
         call winset(.true.)
 
         call polyline(rccleft,rhoave(1,NTHUSED-j+1),NRUSED)
         call polyline(rcc(1),rhoave(1,j),NRUSED)
 
         call dashset(2)
         do k=1,NRUSED
            phip1(k)=phiscale*potave(k,NTHUSED-j+1)+phi0rho
         enddo
         call polyline(rccleft,phip1,NRUSED)
         do k=1,NRUSED
            phip1(k)=phiscale*potave(k,j)+phi0rho
         enddo
         call polyline(rcc(1),phip1,NRUSED)
         call dashset(0)
         call winset(.false.)
         write(charin,'(f4.0)')thang(j)*180./3.1415927
         call legendline(0.7,0.08*(j-1)/jstepth+.1,
     $        0,charin(1:4)//char(0))
      enddo
      call color(15)
      call vecw(-rcc(NRUSED),1.,0)
      call vecw(rcc(NRUSED),1.,1)
      call vecw(-rcc(NRUSED),phi0rho,0)
      call vecw(rcc(NRUSED),phi0rho,1)
c      call fwrite(vd,iwdth,2,charin)
c      call fwrite(time,iw2,2,ch2)
c      call jdrwstr(0.03,0.37,'v!dd!d='//charin(1:iwdth)//' t='
c     $     //ch2(1:iw2),1.)

      end

c*********************************************************************
      subroutine angularPsi(jstepth,rhomin,rhomax,time)
c     Angular distribution of density at the probe edge as a function of
c     theta for various Psi

      include 'piccom.f'
      include 'errcom.f'
      character charin*100,ch2*100
      real rccleft(nrsize),phip1(nrsize)
      real rhoave(0:nthsize,0:npsisize)

      
      do i=1,NTHUSED
         do j=1,NPSIFULL
            rhoave(i,j)=rho(1,i,j)
c            if(rhoave(i,j).le.0) write(*,*)i,j,rhoave(i,j)
         enddo
      enddo

      jsteppsi=NPSIUSED/3-1
      if(rhomax.eq.0.)then

c     I had to change the second argument to NRSIZE+1 instead of NRUSED
c     in the 2D version of SCEPTIC. I have no clue why
         call minmax2(rhoave(1,1),NTHSIZE+1,NTHUSED,NPSIFULL,rhomin
     $        ,rhopm)
         rhopm=max(1.2,.1*nint(10.*rhopm))
      else
         rhopm=rhomax
      endif
      
      call pltinit(-1.,1., rhomin,rhopm)
c      call axptset(1.,1.)

      call axptset(0.,0.)
      call scalewn(-1.,1.,0.,rhopm, .false.,.false.)

      call axis()
c      call boxtitle('Density')
      call axlabels('cos!Aq!@','n/n!A!d;!d!@')

 
      do j=1,NPSIUSED,jsteppsi

         call color(mod(j,15)+1)
         call winset(.true.)
 
         call polyline(tcc(1),rhoave(1,j),NTHUSED)
 
         call winset(.false.)
         write(charin,'(f4.0)')pcc(j)*180./3.1415927
         call legendline(0.7,0.08*(j-1)/jsteppsi+.1,
     $        0,charin(1:4)//char(0))
      enddo
      call color(15)
      call vecw(-1.,1.,0)
      call vecw(1.,1.,1)
      call vecw(-1.,phi0rho,0)
      call vecw(1.,phi0rho,1)

      end
