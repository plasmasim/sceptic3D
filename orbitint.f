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

      program orbitint
      include 'piccom.f'
      real cflux(nthsize)
      real vdist(nvel),vcd(nvel)
      real tdist(nthsize),ftdist(nthsize)
c crdist angular distribution of launches at infinity
c cidist theoretical ditto
c tdist angular distribution at Boundary.
c cflux coulflux angular distrib at boundary if inverse-square.
c vdist velocity distribution of launches at infinity.
c vdanal analytic theory of ditto.
c fventer velocity distribution entering at Boundary.
      real crdist(nthsize),cidist(nthsize)
      common/rtest/crdist,cidist,tdist,vdist
      real vdanal(nvel),vdnumer(nvel),cndist(nQth)
      real Iu2
      real venter(nvel),fventer(nvel)
      logical istrapped
      character*100 string,charin

      call pfset(2)
c Full size arrays by default. Can be changed later by switches.
      nr=nrsize
      nth=nthsize
      if(LCIC)nth=nthsize-1
      if(LCIC)then
         if(nth.gt.nthsize-1)then
            write(*,*)'Too many theta points:',nth,'  Set to',nthsize-1
            nth=nthsize-1
         endif
         NTHUSED=nth
         NTHFULL=nth+1
      else
         if(nth.gt.nthsize)then
            write(*,*)'Too many theta points:',nth,'  Set to',nthsize
            nth=nthsize
         endif
         NTHUSED=nth-1
         NTHFULL=nth
      endif
c Using more than 200000 seems to make the code hang. Don't know why.
      ninjects=100000
      Ti=1.
      vd=.1
c Don't make much smaller. Messes up normalization:
      dt=.0001
      rmax=5.
      r(1)=0.
      vprobe=-5.
      nrein=0
      ninner=0
      cftot=0
      debyelen=1.
      adeficit=0.
      icolntype=0
      bcr=0
      Pc(nQth,nvel)=0.

c Initialize the mesh and poisson coefficients
      if(LCIC)then
         call meshinitcic(rmax)
c         call poisinitcic()
      else
         call meshinitngp(rmax)
c         call poisinitngp()
      endif
c Signal that we have not read a file.
      diagrho(1)=0.
c Deal with arguments
      if(iargc().eq.0) goto 3
      do 1 i=1,iargc()
         call getarg(i,string)
         if(string(1:1) .eq. ' ') then
            goto 3
         endif
         if(string(1:2) .eq. '-n') read(string(3:),*)ninjects
         if(string(1:2) .eq. '-a') read(string(3:),*)adeficit
         if(string(1:2) .eq. '-t') read(string(3:),*)Ti
c         if(string(1:2) .eq. '-x') read(string(3:),*)rmax
         if(string(1:2) .eq. '-v') read(string(3:),*)vd
         if(string(1:2) .eq. '-p') read(string(3:),*)vprobe
         if(string(1:2) .eq. '-l') read(string(3:),*)debyelen
         if(string(1:3) .eq. '-kt') read(string(4:),*)icolntype
         if(string(1:1) .ne. '-' ) then
c Try to read file:
            write(*,*)'Trying to read file',string
            open(10,file=string,status='old',err=101)
c     Line for nothing.
            read(10,*)charin
            read(10,'(2f8.5,f8.4,i6,f12.4,f12.6,f7.4,2f14.5)',err=201)
     $           dtf,vd,Ti,istepsf,rhoinff,phiinff,
     $           favef,debyelen,vprobe
 201        continue
c            write(*,*)'dt,vd,Ti,isteps,rhoinf,phiinf,fave,debyelen,Vp'
c            write(*,*)dt,vd,Ti,isteps,rhoinf,phiinf,fave,debyelen,vprobe
            read(10,*,err=101)nrTi
c Overload ninner to transmit the number.
            ninner=nrTi
            if(ninner.gt.nr)then
               write(*,*)'r-array to long for storage'
               call exit()
            endif
c     write(*,*)'nrTi=',nrTi
            do ir=1,nrTi
               read(10,*,err=101)rcc(ir),diagphi(ir)
               write(*,*)rcc(ir),diagphi(ir)
            enddo

c Continue reading until we have the angular distribution that sceptic 
c calculated.
            read(10,*)charin
c     write(*,*)charin(1:78)
            read(10,*,err=101)nsteps
c     write(*,*)nsteps
            if(nsteps.gt.nstepmax) then
               write(*,*)'Number of steps',nsteps,
     $              ' exceeds allocation',nstepmax
               call exit
            endif
            read(10,*)(fluxprobe(j),j=1,nsteps)
c     Read theta cells
            read(10,*)charin
c     write(*,*)charin(1:78)
            read(10,*)nthhere,nsteps
c     write(*,*)nthhere,nsteps
            do ir=1,nsteps
               read(10,*)(ninthstep(j,ir),j=1,nthhere)
            enddo
            
            read(10,*)charin
c     write(*,*)charin(1:78)
            read(10,*)nastep
c     write(*,*)'nastep',nastep
            nread=nthhere
            read(10,*)(ninth(j),j=1,nread)
            write(*,*)'nastep',nastep,' ninth:'
            write(*,*)( ninth(j),j=1,nthhere)
            close(10)

            diagrho(1)=999.
 101        write(*,*)'End of file reading'
c            stop
         endif

 1    continue
 3    continue
c r(nr) is needed by reinject.
      r(nr)=rmax
      write(*,*)'Initializing injection, icolntype=',icolntype
c Testing new injinit
      call injinit(icolntype,bcr)


      call autoplot(vcom,pu1,nvel)
      call polyline(vcom,pu2,nvel)
      call boxtitle('Velocity distribution functions pu1, pu2')
      call pltend()
c Testing.
c      stop
c
      call finit()
      if(diagrho(1).eq.999.)then
c         initialize to set averein
         call reinject(1,dt,icolntype,bcr)
      else
         averein=vprobe/rmax
      endif
      venr=sqrt(2.*Ti)*Vcom(nvel)-2.*averein
      if(icolntype.eq.2.)then
         smaxcomp=sqrt(2.)*2.*3.141593**2*(pu1(1)-pu2(1)*averein/Ti)
         write(*,*)'Ti,averein,vd=',Ti,averein,vd
         write(*,*)'The following numerical vs analytic values ',
     $        'should agree for zero drift'
         write(*,*)'smaxcomp=',smaxcomp
         write(*,*)'smaxflux=',smaxflux(vd/sqrt(2.*Ti),(-averein/Ti))
      endif

      do i=1,nth
         tdist(i)=0.
         crdist(i)=0.
         if(i.le.nth-1)then
c            write(*,*)acos(th(i)),(-averein/Ti),vd/sqrt(Ti)
            cflux(i)=coulflux(acos(th(i)),(-averein/Ti),vd/sqrt(Ti))
            cftot=cftot+cflux(i)
         endif
      enddo
      do i=1,nth
         cflux(i)=cflux(i)/cftot
      enddo
      do i=1,nvel
c         vdist(i)=0.
         fventer(i)=0.
         venter(i)=i*venr/nvel
      enddo
      ntrapped=0
      write(*,*)'Doing ',ninjects,' injections. Type=',icolntype
      do k=1,ninjects
         kp=1+mod(k-1,npartmax)
         call reinject(kp,dt,icolntype,bcr)
c         write(*,501)(xp(j,kp),j=1,6)
c 501  format('Pos=',3f10.4,' Vel=',3f10.4)
         ct=xp(3,kp)/sqrt(xp(1,kp)**2+xp(2,kp)**2+xp(3,kp)**2)
         ith=(1.-ct)*0.5*(nth-1) + 1
         tdist(ith)=tdist(ith)+1
         v=sqrt(xp(4,kp)**2+xp(5,kp)**2+xp(6,kp)**2)
         ivdiag=1+max(0,nint(nvel*(v/venr)))
         if(ivdiag.gt.nvel) ivdiag=nvel
         fventer(ivdiag)=fventer(ivdiag)+1.
         if(istrapped(kp)) ntrapped=ntrapped+1
      enddo
c      write(*,*)'venter,   fventer'
c      write(*,'(f10.3,f10.1)')(venter(j),fventer(j),j=1,nvel)
      ts=nthused/2.
      do i=1,NTHUSED
         cidist(i)=ts*Iu2(vd/sqrt(2.*Ti),(-th(i)),(-averein/Ti))/cftot
         tdist(i)=tdist(i)/ninjects
         crdist(i)=ts*crdist(i)/nrein
      enddo
      if(LCIC)then
         crdist(NTHUSED)=2.*crdist(NTHUSED)
         crdist(1)=2.*crdist(1)
      endif
      write(*,*)'averein=',averein,'  vd=',vd,'  Ti=',Ti,
     $     '  adeficit=',adeficit
      write(*,*)'nrein=',nrein,'  ninject=',ninjects,
     $     '  vdpositive=',ninner,
     $     '  ntrapped=',ntrapped
c      write(*,'(3f10.4)')(th(j),tdist(j),cflux(j),j=1,nth-1)
c      write(*,*)cidist
c      write(*,'(2f10.4)')(cidist(j),crdist(j),j=1,nth-1)

      if(LCIC)then
c adjust the tcc mesh to account for proper averaging.
         tcc(1)=tcc(1)+0.25*(tcc(2)-tcc(1))
         tcc(NTHUSED)=tcc(NTHUSED)+0.25*(tcc(NTHUSED-1)-tcc(NTHUSED))
      endif
      call automark(tcc(1),crdist,NTHUSED,2)
      if(Pc(nQth,nvel).ne.0. .and. icolntype.eq.2)then
c We are using the numerical integrations. Get the angular distrib
c from the place they are stored, in Gcom (1 and 2)
         do j=1,nQth
            cndist(j)=(Gcom(1,j)-Gcom(2,j)*(averein/Ti))/
     $           (Gcom(3,1)-Gcom(3,2)*(averein/Ti))
c Adjust plot position Not any more:
c            if(j.lt.nQth) Qcom(j)=(Qcom(j)+Qcom(j+1))/2.
         enddo
c          write(*,*)'Qcom(1:nQth)=',Qcom(1),Qcom(nQth)
         call polyline(Qcom,cndist,nQth)
c         call polymark(Qcom,cndist,nQth,3)
      else
         write(*,*)'Comparing with cidist'
         call polyline(th(1),cidist,NTHUSED)
      endif
      call boxtitle('Angular distribution at infinity.')
      call axlabels('cos!Aq!@','')
      call pltend()

      Uc=abs(vd)/sqrt(2.*Ti)
      do i=1,nvel
         vcd(i)=Vcom(i)+0.5*(Vcom(2)-Vcom(1))
c Old approach to normalization.
c         vdist(i)=vdist(i)*((nvel-1)/(5.+abs(uc)))*
c     $        (pu1(1)-(averein/Ti)*pu2(1))
c     $        /float(nrein)
         vdanal(i)=(exp(-(Uc-Vcom(i))**2)-exp(-(Uc+Vcom(i))**2))*
     $        (Vcom(i)**2 - averein/Ti)
     $        /(pu1(1)-(averein/Ti)*pu2(1))
c Normalization so that \int vdist dv = 1 :
         vdist(i)=vdist(i)/Vcom(2)/float(nrein)
c This is not the analytical approx but the numerical:
         vdnumer(i)=Pc(nQth,i)*(Vcom(i)**2-averein/Ti)*Vcom(i)
     $        /(pu1(1)-(averein/Ti)*pu2(1))

      enddo
      call automark(Vcd,vdist,nvel,1)
      if(Pc(nQth,nvel).ne.0.)then
         call polyline(Vcom,vdnumer,nvel)
         call legendline(.6,.2,0,' numerical')
      else
         call polyline(Vcom,vdanal,nvel)
         call legendline(.6,.2,0,' analytical')
      endif
      call legendline(.6,.15,1,' reinjects')
      call boxtitle('Velocity distribution at infinity.')
      call pltend()

c      write(*,*)'tdist',tdist
      call automark(th(1),tdist,nth-1,1)
      call boxtitle('Angular distribution on BC radius & coulflux')
      call axlabels('cos!Aq!@','')
      call winset(.true.)
      call polyline(th(1),cflux,nth-1)

      if(nthhere.eq.nthused)then
c We can do a comparison with file.
         nisum=0
         do j=1,nthhere
            nisum=nisum+ninth(j)
         enddo
         if(LCIC)nisum=nisum+ninth(1)+ninth(nthhere)
         do j=1,nthhere
            ftdist(j)=cftot*float(ninth(j))/nisum
c            write(*,*)j,ftdist(j)
         enddo
         if(LCIC)then
            ftdist(1)=2.*ftdist(1)
            ftdist(nthhere)=2.*ftdist(nthhere)
         endif
      call polymark(tcc(1),ftdist,nthhere,3)
      endif
      call legendline(.4,.05,3,' SCEPTIC')
      call legendline(.4,.1,1,' orbit integral')
      call winset(.false.)
      call jdrwstr(.05,.03,'(Agreement only when lambda is large)',1.)
      call pltend()

      call automark(venter,fventer,nvel,1)
c      call polyline(venter,fventer,nvel)
c Don't know why these avereins don't get scaled by /Ti.
c I guess the venr is in standard, not ion, units.
      call vecw(sqrt(-averein*2.),0,0)
      call vecn(wx2nx(sqrt(-averein*2.)),.7,1)
      call boxtitle('Flux distribution at boundary.')
      call axlabels('v','')
      call pltend()

      end
c**********************************************************************
c Return angular distribution of flux collected from a shifted maxwellian
c distribution at infinity, as a function of angle (at infinity)
c cos(theta)=c  of the velocity to the drift velocity, U. 
c Normalized to sqrt(2T/m). Collecting potential is -chi.
c When U=0 its value is 1+chi.
      real function Iu2(U,c,chi)
      real U,c
      if(c.lt.-1. .or. c.gt.+1.) then
         Iu2=0.
         write(*,*)'Error Iu2. cosine out of range'
      endif
      uc=U*c
      uc2=uc**2
      temp=sqrt(3.141593)*0.5*erfcc((-uc))*(uc*(3.+2.*uc2)+ 2.*uc*chi)
      Iu2=exp(-U**2)*(uc2+1.+chi) + exp(-U**2+uc2)*temp
c      Iu2=exp(-U**2)*(uc2+1.+chi) + exp(-U**2*(1.-c**2))*temp
      end
