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


c Main program for pic code
      program sceptic3D

c     For now, put maxsteps in piccom
      integer ipstep
      logical success,finaldiags,ieradset
      real rmax
      character*100 string
      character*10 cfinal
      logical readpart,writepart
      logical lcolcont,lpstore
      logical lsmoothT,lsmoothP
      integer m2,rshield

      logical ltmp

c Communicator and id for the conjugate gradient communicator
c In the case of cgparallel=false, myid2=0 only
c does the potential calculation
      integer cg_comm,myid2

c Common storage
      include 'piccom.f'
      include 'errcom.f'
      include 'colncom.f'
      include 'fvcom.f'
      data rmax/5./
      data dtf/0.025/bdt/1./
      data success/.false./
      data readpart/.false./
      data writepart/.false./

c     Variables used for looping over CG for nonlin Poisson solver
      real dconverge,dconvergenl,maxdphi
      integer niter,totiter
      logical nonlin
      data nonlin/.false./

c     Variables used for calculating matrix A for debugging
      real inputvect(nrsize-1,0:nthsize,0:npsisize),
     $  outputvect(nrsize-1,0:nthsize ,0:npsisize)
      integer n1,n2,n3,m,n,o


c Parallel processing MPI options.
#ifdef MPI
      include 'mpif.h'
      double precision cgtime
      integer rc
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
#else
      myid=0
      numprocs=1
      cg_comm=0
#endif

c      write(*,*)'Starting',myid

      
c Variables for debugging
      lsavemat = .false.
      lsavephi = .false.
      lfirsttime = .true.
      stepcount = 0
      saveatstep = 1
c     Use biconjugate gradient method as default solver
      lbcg=.true.

c Full size arrays by default. Can be changed later by switches.
      nr=nrsize
      nth=nthsize
      if(LCIC) nth=nthsize-1
      npsi=npsisize-1
c Note: npsi should not be set lower than 4 since this causes indexing
c   problems in the parallel solver, and nth cannot be lower than 3
c   (determined during debugging).

c Default MC injection values
      mcrntheta = 5
      mcrnpsi = 5
      mcrnpart = 100000
c Set cout of injected particles to zero
      mcrninjd = 0
      

c Common defaults. Avoid block data.
      vd=.000
      cd=1.0
      cnd=1.0
      cB=1.
      psind=pi/2.
      debyelen=.1
      vprobe=-4.
      Ezext=0.
      bdyfc=4.
      Ti=1.
      Tneutral=1.
      diags=.true.
      samp=.false.
      diagsamp=4
      finaldiags=.true.
      ldist=.false.
      avelim=0.6
      lsmoothT=.false.
      lsmoothP=.false.
      lcolcont=.true.
      ircell=1
      itcell=1
      norbits=0
      npart=npartmax
      ninjcomp0=npartmax
      ninjcomp=ninjcomp0
      linsulate=.false.
      lfloat=.false.
      rmtoz=1.
      rhomin=0.
      rhomax=0.
c     lat0 and lap0 mean that we set the polar or poloidal acceleration
c     artificially to zero
      lat0=.false.  
      lap0=.false.

c     Don't go to output unless triggered
      lgotooutput = .false.

      ieradset=.false.

      cgparallel=.false.
      collcic=.true.
      maxsteps=500

      colnwt=0.
      icolntype=0
      localinj=.false.
      ipfsw=3
      call pfset(0)
      lfixedn=.true.
      ipstep=1
      vneutral=0.
c     bcphi : BC on potential
      bcphi=0
c     bcr : BC on reinjection 1 by default. In fact sofar I have only
c     set up bcr1 to work on SCEPTIC3D. This corresponds to a standard
c     drifting Maxwellian reinjection.
      bcr=1
      Bz=0.
      infdbl=.false.
      orbinit=.false.
      lsubcycle=.false.
      lfulloutput=.false.
      verlet=.false.
      bohm=.false.
c Signal that fvcom is not initialized. After initialization it is .ne.0
      qthfv(nthfvsize)=0.

      

c Deal with arguments
      if(iargc().eq.0) goto 51
      do 1 i=1,iargc()
         call getarg(i,string)
         if(string(1:1) .eq. ' ') then
            goto 3
         endif
         if(string(1:8) .eq. '--rhoinf')then
            read(string(9:),*)rhoinf
            lfixedn=.false.
         endif
         if(string(1:2) .eq. '-r') readpart=.true.
         if(string(1:2) .eq. '-w') writepart=.true.
         if(string(1:3) .eq. '-gc') then
            lcolcont=.false.
         elseif(string(1:3) .eq. '-ge') then
         elseif(string(1:2) .eq. '-g')then
            if(string(3:3).ne.' ')then
               read(string(3:),*,err=259,end=259)ipstep
            else
               diags=.false.
            endif
            goto 260
 259        ipstep=1
 260        continue
         endif
         if(string(1:8) .eq. '--subcyc')then
            lsubcycle=.true.
         endif
         if(string(1:12) .eq. '--fulloutput')then
            lfulloutput=.true.
         endif

         if(string(1:2) .eq. '-f') finaldiags=.false.
         if(string(1:3) .eq. '-er') then
            ieradset=.true.
            read(string(4:),*)ierad
         elseif(string(1:2) .eq. '-e') then
            ldist=.true.
c try to read the cell from the rest of the string.
            read(string(3:),*,err=261,end=261)ircell,itcell
            goto 262
c restore defaults on read error.
 261        ircell=1
            itcell=1
 262        continue
         endif
         if(string(1:4) .eq. '--bz') then 
            read(string(5:),*)Bz
         elseif(string(1:5) .eq. '-bohm') then 
            bohm=.true.
         endif
#ifdef MPI
         if(string(1:4) .eq. '--sp') then
            cgparallel=.true.
         endif
#endif
c        For debugging, allow use of minimum residual method
         if(string(1:8) .eq. '--minres') then
            lbcg=.false.
         endif
c        For debugging, give option to save A and/or phi
         if(string(1:9) .eq. '--savemat') then
            lsavemat=.true.
         elseif(string(1:9) .eq. '--savephi') then
            lsavephi=.true.
         endif
         
         if(string(1:2) .eq. '-s') then
            read(string(3:),*)maxsteps
         endif
         if(string(1:2) .eq. '-b') read(string(3:),*)bdyfc
         if(string(1:2) .eq. '-t') read(string(3:),*)Ti
         if(string(1:3) .eq. '-tn') then
            read(string(4:),*)Tneutral
         else
            Tneutral = Ti
         endif
         if(string(1:2) .eq. '-d') read(string(3:),*)dtf
         if(string(1:5) .eq. '-cmin')then
            read(string(6:),*)rhomin
         elseif(string(1:5) .eq. '-cmax')then
            read(string(6:),*)rhomax
         elseif(string(1:3) .eq. '-cd')then
            read(string(4:),*)cd
         elseif(string(1:3) .eq. '-cB')then
            read(string(4:),*)cB
         elseif(string(1:2) .eq. '-c') then
            read(string(3:),*)bdt
         endif
         if(string(1:2) .eq. '-x') read(string(3:),*)rmax
         if(string(1:3) .eq. '-vn') then
            read(string(4:),*)vneutral
         elseif(string(1:4) .eq. '-cnd')then
            read(string(5:),*)cnd
         elseif(string(1:6) .eq. '-psind')then
            read(string(7:),*)psind
         elseif(string(1:4) .eq. '-ver') then
            verlet=.true.            
         elseif(string(1:2) .eq. '-v')then
            read(string(3:),*)vd
         endif
         if(string(1:3) .eq. '-pf') then
            lfloat=.true.
         elseif(string(1:3) .eq. '-pi') then
            linsulate=.true.
         elseif(string(1:2) .eq. '-p') then
            read(string(3:),*)vprobe
         endif
         if(string(1:2) .eq. '-z') read(string(3:),*)Ezext
         if(string(1:2) .eq. '-l') then
            read(string(3:),*)debyelen
         endif
         if(string(1:3) .eq. '-nr') read(string(4:),*)nr
         if(string(1:3) .eq. '-nt') read(string(4:),*)nth
         if(string(1:3) .eq. '-np') read(string(4:),*)npsi
         if(string(1:3) .eq. '-ni') read(string(4:),*)npart
         if(string(1:2) .eq. '-m') read(string(3:),*)rmtoz
         if(string(1:2) .eq. '-a') then
            read(string(3:),*,err=251,end=251)ipfsw
            call pfset(ipfsw)
            goto 252
 251        ipfsw=0
 252        continue
         endif
         if(string(1:5) .eq. '--bcr') read(string(6:),*)bcr
         if(string(1:7) .eq. '--bcphi') read(string(8:),*)bcphi
         if(string(1:7) .eq. '-infdbl') infdbl=.true.
         
         if(string(1:3) .eq. '-oi') then
            read(string(4:),*)norbits
            orbinit=.true.
         elseif(string(1:2) .eq. '-o') then
            read(string(3:),*)norbits
         endif
         if(string(1:3) .eq. '-kt')then
            read(string(4:),*)icolntype
         elseif(string(1:2) .eq. '-k') then
            read(string(3:),*)colnwt
         endif
         if(string(1:2) .eq. '-?') goto 51
         if(string(1:9) .eq. '--smoothP') then 
            lsmoothP=.true.
         elseif(string(1:9) .eq. '--smoothT') then 
            lsmoothT=.true.
         endif
         if(string(1:5) .eq. '--at0') lat0=.true.
         if(string(1:5) .eq. '--ap0') lap0=.true.
         if(string(1:10) .eq. '--localinj') localinj=.true.
         if(string(1:8) .eq. '--nonlin') nonlin=.true.
        
c     MC Reinject variables
         if(string(1:7) .eq. '--mcrnt') read(string(8:),*)mcrntheta
         if(string(1:7) .eq. '--mcrnp') read(string(8:),*)mcrnpsi
         if(string(1:7) .eq. '--mcrni') read(string(8:),*)mcrnpart


 1    continue
 3    continue

c Check ion velocity angle and set velocity vector
      if(cd.lt.-1) then
         write(*,*)'Illegal ion velocity theta, set to cd=-1'
         cd=-1
      elseif(cd.gt.1) then
         write(*,*)'Illegal ion velocity theta, set to cd=1'
         cd=1.
      endif

c Drift vector taken to be in y-z plane
      drvect(3)=vd*cd
      drvect(2)=(vd*sqrt(1-cd**2))*1.
      drvect(1)=(vd*sqrt(1-cd**2))*0.

c Check neutral velocity angles and set velocity vector
      if(cnd.lt.-1) then
         write(*,*)'Illegal neutral velocity theta, set to cnd=-1'
         cnd=-1
      elseif(cnd.gt.1) then
         write(*,*)'Illegal neutral velocity theta, set to cnd=1'
         cnd=1.
      endif

      if((psind.lt.0).or.(psind.ge.(2*pi))) then
         psind = acos(cos(psind))
      endif

      vneut(3)=vneutral*cnd
      vneut(2)=(vneutral*sqrt(1-cnd**2))*sin(psind)
      vneut(1)=(vneutral*sqrt(1-cnd**2))*cos(psind)

c Set relative drift
      do i=1,3
c        Relative drift taken to be what is specified with vd
         reldrift(i) = drvect(i)
cc         reldrift(i) = drvect(i) - vneut(i)
      enddo

c Check magnetic field angle and set magnetic field vector
      if(cB.lt.-1) then
         write(*,*)'Illegal magnetic field theta, set to cB=-1'
         cB=-1
      elseif(cB.gt.1) then
         write(*,*)'Illegal magnetic field theta, set to cB=1'
         cB=1.
      endif

c Magnetic field taken to be in y-z plane
      Bvect(3)=Bz*cB
      Bvect(2)=(Bz*sqrt(1-cB**2))*1.
      Bvect(1)=(Bz*sqrt(1-cB**2))*0.
      if (Bz.ne.0.) then
c        Set direction of magnetic field
         do i=1,3
            magdir(i) = Bvect(i)/Bz
         enddo
c        Find ExB drift
         do i=1,3
            ecbdrift(i) = reldrift(i) -
     $        dot(magdir(1),reldrift(1),3)*magdir(i)
         enddo
      else
         do i=1,3
            magdir(i) = 0.
            ecbdrift(i) = 0.
         enddo
      endif


c Set Array sizes, allowed variable.
      if(nr.gt.nrsize-1)then
         write(*,*)'Too many radial points:',nr,'  Set to',nrsize-1
         nr=nrsize-1
      endif
      if(npsi.gt.npsisize-1)then
         write(*,*)'Too many poloidal points:',npsi,' Set to',npsisize-1
         npsi=npsisize-1
      endif
      if(npart.gt.npartmax)then
         write(*,*)'Too many ions:',npart,'  Set to',npartmax
         npart=npartmax
      endif
      if(nth.gt.nthsize-1)then
         write(*,*)'Too many theta points:',nth,'  Set to',nthsize-1
         nth=nthsize-1
      endif
      NRUSED=nr
      NTHUSED=nth
      NRFULL=nr+1
      NTHFULL=nth+1
      NPSIUSED=npsi
      NPSIFULL=npsi+1

c Check that MC Reinjection Array sizes are within limits
      if(mcrntheta.gt.mcrthetasize)then
         write(*,*)'Too many theta points for MC reinjection:',
     $     mcrntheta,'  Set to',mcrthetasize
         mcrntheta=mcrthetasize
      endif
      if(mcrnpsi.gt.mcrpsisize)then
         write(*,*)'Too many poloidal points in MC reinjection:',
     $     mcrnpsi,' Set to',mcrpsisize
         mcrnpsi=mcrpsisize
      endif
      if(mcrnpart.gt.mcrpartsize)then
         write(*,*)'Too many ions in MC reinjection:',
     $     mcrnpart,'  Set to',mcrpartsize
         mcrnpart=mcrpartsize
      endif
      

      if(.not.ieradset) then
         ierad=nrused
      endif
c Set when we start the averagings for forces, fluxes ... (m2)
c Averages over the last 1/4 of the steps
      m2=3*maxsteps/4

c Orbit tracking setting. This is when we want to follow orbits of particles
c having a thermal speed directed in the z direction, but different impact
c parameters
      trackinit=maxsteps/2
      
      if(myid.eq.0)then
c         write(*,*)'Command line and defaults:'
         write(*,505)nr,nth,npsi,npart,Bz
         write(*,502)debyelen,Ti,maxsteps,dtf,rmax,vd,vprobe
 505     format(' Mesh: ',i3,' x',i3,' x',i3,'  Particles:',i7,' Bz:',f6
     $        .2)
 502     format('l_d=',f10.3,' Ti=',f7.3,' Steps=',i4,
     $        ' dt=',f6.4,' rmax=',f5.1,' vd=',f6.3,' vp=',f8.4)
         if(lsubcycle)write(*,*)'Subcycling on!'
         if(nonlin)write(*,*)'Nonlinear Poisson solver used!'
      endif
      if(maxsteps.gt.nstepmax-1) Stop 'Too many steps requested.'
      lplot=diags.and.(myid.eq.0)
      k=0

c     Parallel solver stuff

c Don't need the parallel solver for 0 Debye length
      if ((debyelen.eq.0).or.infdbl) cgparallel=.false.

#ifdef MPI
      if (cgparallel) then
         call cgparinit(myid2,cg_comm)
      elseif (myid.eq.0) then
         myid2=0
      else
         myid2=-1
      endif
#else
      myid2=0
#endif
c Initialize the mesh and poisson coefficients


      call meshinitcic(rmax)
      call poisinitcic()


c     For debugging, don't do this since clutters output
c      if(nr.eq.10)then
c         write(*,*)'apc'
c         write(*,521)apc
c         write(*,*)'bpc'
c         write(*,521)bpc
c         write(*,*)'cpc'
c         write(*,521)cpc
c         write(*,*)'dpc'
c         write(*,521)epc
c         write(*,*)'epc'
c         write(*,521)dpc
c         write(*,*)'fpc'
c         write(*,521)((fpc(i,j),i=1,10),j=0,NTHFULL)
c 521     format(10g8.2)
c      endif

c Set seed for rand(), which is called through the reshuffler ran0()
      call srand(myid+1)
c Initialize collisions and related fields.
      call colninit(colnwt,icolntype)
c Initialize the fields. (Changed the order of finit and injinit)
      call finit()
c Initialize the random functions for reinjection.
      call injinit(icolntype,bcr,colnwt)
c Initialize velocity diagnostics
      vrange=8.*sqrt(Ti)+1.4*abs(vd)+1.4*sqrt(abs(vprobe))
      do kk=1,nvmax
         nvdiag(kk)=0
         vrdiagin(kk)=0
         vtdiagin(kk)=0
         vdiag(kk)=vrange*(float(kk-1)/nvmax-0.499)
      enddo

      if(.not.lfixedn)then
         ninjcomp0=dtf*rhoinf*sqrt(Ti)*
     $        smaxflux(vd/sqrt(2.*Ti),0.)*r(NRUSED)**2
         ninjcomp=ninjcomp0
         if(myid.eq.0) write(*,*)'Initial:  rhoinf=',rhoinf
     $        ,'  ninjcomp=',ninjcomp
      endif
      iocprev=npart

c Read in the previous particle distribution and averages.
      if(readpart) call partrd(success)
c     For debugging, notify if read particles
      if(success) write(*,*)'Read particles.'
c Or Initialize (load) particles
      

      if(.not.success) call pinit(colnwt)
c Need to reset i to avoid problems at RhoDiag ... initialization
      i=0
      samp=.true.
c Assign charge to mesh
      call chargetomesh()      
c Collect the partial sums of moments of distribution.
      call sumreduce()

      call rhocalc(lsmoothT,lsmoothP,i,dt)

c Don't forget to put back samp to false.
      samp=.false.

c For debugging, don't do this since clutters output
cc Print out diagnostics for tiny meshes.
c      if(nr.eq.10) then
c         write(*,*)'rho,phi,psum'
c         write(*,504)(((rho(iw,jw,kw),iw=1,nr),jw=1,NTHUSED),
c     $        kw=1,NPSIUSED)
c         write(*,504)(((phi(iw,jw,kw),iw=1,nr),jw=0,NTHFULL),
c     $        kw=1,NPSIFULL)
c         write(*,503)(((psum(iw,jw,kw),iw=1,nr),jw=0,NTHFULL),
c     $        kw=1,NPSIFULL)
c      endif
      if(lplot)then
         if(ldist) then
            call multiframe(2,3,3)
         else
            call multiframe(2,2,3)
         endif
      endif

      time=0.
c Save the permanent plot switch.
      lpstore=lplot
      maccel=maxsteps/3
      
c Must start nstepsave at 1, not 0
      nstepsave=1
      nsamax=min(199,maxsteps/20)+1
      if(success)nstepsave=nsamax

c     Timing
#ifdef MPI
      cgtime=MPI_WTIME()
#endif
      
      if (myid.eq.0) then
         write(*,*) "Maxsteps : ",maxsteps
      endif

c Main Stepping loop.
      do i=1,maxsteps
         

c Plot at some subset of steps.
c         write(*,*)mod(i-1,ipstep),i,ipstep
         lplot=lpstore.and.(mod(i-1,ipstep).eq.0)
c Acceleration for first 1/3rd of steps.
c The following two expressions are sensibly similar
c        bdtnow=max(1.,(bdt-1.)*(maccel-i+2)/(maccel+1.)+1.)
         if (maccel .ne. 1.) then
            bdtnow=max(1.,bdt*(maccel-i)/(maccel-1.))
         else
            bdtnow = 1.
         endif
         dt=bdtnow*dtf
c         write(*,*) colnwt,Eneutral
c         write(*,*)'bdtnow,btd,maccel,dt,dtf=',bdtnow,bdt,maccel,dt,dtf
         ninjcomp=bdtnow*ninjcomp0
c     Here there is a problem if we have ipstep!=1 because it is
c     essential that right after changing pfset one should call
c     pltinit. Otherwise initialization may be incorrect. This problem
c     is avoided by ensuring that the lplot is not changed back to true,
c     on leaving the loop, until we are about to reinitialize plotting.
c     As a result the last plot is only saved if it is on the maxsteps step.

c     If the flag orbinit is on, then we force the tracked particles to
c     have the following properties.
         if (myid.eq.0.and.i.eq.trackinit.and.orbinit) then 
            h=1./2*rmax
            do k=1,norbits
               iorbitlen(k)=0
               xp(3,k)=-4*rmax/5;
               xp(1,k)=0
               xp(2,k)=h*k/norbits*sin(2*pi*k/norbits)
               xp(1,k)=h*k/norbits
               xp(2,k)=0
               xp(6,k)=vd+0*sqrt(2*Ti)*1.
               xp(4,k)=0.
               xp(5,k)=0.
            enddo
         endif

         if(i.eq.maxsteps)call pfset(ipfsw)

c Assign charge to mesh
         call chargetomesh()

c Collect the partial sums of moments of distribution.
#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
      
c Check if any node tripped lgotooutput, and if so, leave do loop
#ifdef MPI
         ltmp = lgotooutput
         call MPI_ALLREDUCE(lgotooutput,ltmp,1,MPI_INTEGER,MPI_LOR
     $        ,MPI_COMM_WORLD,ierr)
         lgotooutput = ltmp
#endif
         if (lgotooutput) goto 50

         call sumreduce()
         call aveupstep()

        
c     One step prior to when the momenta are needed, set samp=true so
c     that chargediag calculates them.
         if(i.ge.(maxsteps-diagsamp*nstepsave)-diagsamp) then
            if(mod(i-(maxsteps-diagsamp*nstepsave)+diagsamp,diagsamp).eq
     $           .(diagsamp-1))then
               samp=.true.
            else
               samp=.false.
            endif
         endif

c myid.eq.0 for the master node.
         if(myid.eq.0)then
            
c     Plot density contours.
            if(lplot)then
               if(lcolcont)then
                  call rhodisplay(10,rhomin,rhomax)
                  if(norbits.gt.0) call plotorbits()
                  call rhodisplayPsi(10,rhomin,rhomax)
               else
                  call pltinit(0.,1.,0.,1.)
               endif
            endif            
            call rhocalc(lsmoothT,lsmoothP,i,dt)
         endif
         
#ifdef MPI
c     Broadcast back to the slaves.
         if (cgparallel) then
            call MPI_BCAST(rho,(nrsize+1)*(nthsize+1)*(npsiused+2),
     $           MPI_REAL,0,MPI_COMM_WORLD,ierr)
         endif
#endif

         if(myid.eq.0) then          
c Plot the slices through the probe.
            if(lplot) then 
               call slices(jstepth,rhomin,rhomax,time)
               call angularPsi(jsteppsi,rhomin,rhomax,time)
            endif
               
c Document charge accumulated; calculate diagrho, finthave.
            call chargediag(dt,i,icolntype)
c Calculate rhoinf
            call rhoinfcalc(dt,icolntype,colnwt)
         endif

 
c If we use floating potential, the slaves may need rhoinf and fluxes
#ifdef MPI
c     Needs to bcast rhoinf in any case to get the last line of the
c     particles diagnostic file right, ie if -w
         if (cgparallel.and.(linsulate.or.lfloat)) then
            call MPI_BCAST(rhoinf,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(fincellave,nthsize*npsiused,MPI_REAL,0,
     $           MPI_COMM_WORLD,ierr)
         endif
#endif
         if(myid.eq.0) then
c  Write step information.
            cfinal=' ReinjFrac'
            if(icolntype.gt.0)cfinal=' ColnFreq'
            if(linsulate.or.lfloat)cfinal=' Vprobe'
            if(mod(i,5).eq.1) then
               if(i.gt.1)then
                  mtrap=mtrapped()
                  if(linsulate.or.lfloat)then
                     yvpre=vprobe
                  else
                     if(icolntype.gt.0 .and. debyelen.ne.0.) then
                        yvpre=NCneutral/(dt*npart)
                     elseif(lfixedn)then
c                        yvpre=fluxrein/nrein
                        yvpre=float(nrein)/npart
                     else
                        yvpre=float(nrein)/ninjcomp
                     endif
                  endif
                  
c     Potential diagnostics averaged over psi at the outer boundary. We
c     used diagchi (averaged over the last nstepsave) here, while
c     SCEPTIC2D uses the potential at the current time step.
                  write(*,'(2f7.3,i7,f7.3)') fluxprobe(i-1)/(4.*pi*r(1)
     $                 **2)/rhoinf/dtf, Ti*diagchi(NTHUSED/2), mtrap
     $                 ,yvpre
               endif
               if(.not.lfixedn) write(*,601)ninjcomp,nrein,rhoinf,
     $              bdtnow,iocprev
 601           format('ninjcomp=',i6,'  nrein=',i6,'  rhoinf=',f8.2,
     $              '  bdtnow',f7.4
     $              ,'  iocprev=',i6)
               if(debyelen.ne.0. .and. mod(i,100).eq.1)then
                  write(*,*)
     $              ' Step:CG  Step:CG  Step:CG',
     $              '  Step:CG  Step:CG',
     $              '  Flux, Phiedge, Trap,', cfinal
               endif
            endif
            write(*,'(i5,$)')i
         endif

c        For debugging: Save Poisson matrices if flag set and small grid
         stepcount = stepcount+1
         if (lsavemat .and.
     $     nr .le. nrsizesave .and.
     $     nth .le. nthsizesave .and.
     $     npsi .le. npsisizesave .and.
     $     i .le. nstepssave) then
            lsavemat = .true.
         elseif (i .le. nstepssave) then
            lsavemat = .false.
         endif


c     Calculate the potential field, unless it is fixed.
         if (myid2.eq.0.and.infdbl) then
            call fcalc_infdbl(dt)
         elseif (myid2.eq.0.and.debyelen.eq.0) then
            call fcalc(dt)
         else
            totiter=0
            maxdphi=1.
c            dconverge=1.e-5
c           For debugging, change convergence crit.
            dconverge=1.e-5
            dconvergenl = dconverge
            call shielding_bc(dt,rshield,icolntype,colnwt)
            if(.not.cgparallel.and.myid2.eq.0) then
               if(nonlin) then
 620              if((maxdphi.gt.dconverge.or.dconvergenl.gt.dconverge)
     $              .and. totiter.lt.999) then
c                    For debugging, limit total number of iterations (above)
                     dconvergenl =
     $                 max((maxdphi/nrused/debyelen)**2,dconverge)
                     call shielding3D(dt,rshield,dconvergenl,niter,
     $                    maxdphi)
                     totiter=totiter+niter
cc                    For debugging, indicate each loop
c                     write(*,'(a,$)') "."
                     goto 620
                  endif
               else
                  call shielding3D(dt,rshield,dconverge,niter,maxdphi)
                  totiter=niter
               endif
               write(*,'('':'',i3,$)')totiter
#ifdef MPI
            elseif(cgparallel.and.myid2.ge.0) then
cc           For debugging, broadcast potential and density to slaves
c               call MPI_BCAST(gpc,(nthsize+1)*(npsisize+1)*5,
c     $           MPI_REAL,0,cg_comm,ierr)
c               call MPI_BCAST(rho,(nrsize+1)*(nthsize+1)*(npsiused+2),
c     $           MPI_REAL,0,cg_comm,ierr)
c               call MPI_BCAST(phi,(nrsize+1)*(nthsize+1)*(npsiused+2),
c     $           MPI_REAL,0,cg_comm,ierr)
c               call MPI_BCAST(phiaxis,(nrsize+1)*2*(npsiused+2),
c     $           MPI_REAL,0,cg_comm,ierr)
               if(nonlin) then
 629              if((maxdphi.gt.dconverge.or.dconvergenl.gt.dconverge)
     $              .and. totiter.lt.999) then
c                    For debugging, limit total number of iterations (above)
                     dconvergenl =
     $                 max((maxdphi/nrused/debyelen)**2,dconverge)
                     call shielding3D_par(dt,rshield,cg_comm,myid2,
     $                    dconvergenl,niter,maxdphi)
                     totiter=totiter+niter
c     Broadcast phi and phiaxis for next CG loop
                     call MPI_BCAST(phi,
     $                    (nrsize+1)*(nthsize+1)*(npsiused+2),
     $                    MPI_REAL,0,cg_comm,ierr)
                     call MPI_BCAST(phiaxis,(nrsize+1)*2*(npsiused+2),
     $                    MPI_REAL,0,cg_comm,ierr)
                     call MPI_BCAST(maxdphi,1,
     $                    MPI_REAL,0,cg_comm,ierr)
                     call MPI_BCAST(totiter,1,
     $                    MPI_INTEGER,0,cg_comm,ierr)
                     call MPI_BARRIER(cg_comm,ierr)
cc                    For debugging, indicate each loop
c                     if(myid2.eq.0) write(*,'(a,$)') "."
                     goto 629
                  endif
               else
                  call shielding3D_par(dt,rshield,cg_comm,myid2,
     $                 dconverge,niter,maxdphi)
                  totiter=niter
               endif
c     Output the total number of CG steps
               if(myid2.eq.0) write(*,'('':'',i3,$)')totiter
#endif
            endif
         endif  
         
  

c     Broadcast back to the slaves.
#ifdef MPI
         call MPI_BCAST(phi,(nrsize+1)*(nthsize+1)*(npsiused+2),
     $        MPI_REAL,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(phiaxis,(nrsize+1)*2*(npsiused+2),
     $        MPI_REAL,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(adeficit,1, MPI_REAL,0,MPI_COMM_WORLD,ierr)
c     Averein is always broadcasted from 0 since it is calculated
c     in diags, called only by myid=0, not myid2=0
         call MPI_BCAST(averein,1, MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif

c        For debugging: Save phi for first nstepssave if small grid
         if (lsavephi .and.
     $     nr .le. nrsizesave .and.
     $     nth .le. nthsizesave .and.
     $     npsi .le. npsisizesave .and.
     $     i .le. nstepssave) then
            do j=0,npsi
               do k=0,nth
                  do l=0,nr
                     phisave(l,k,j,i) = phi(l,k,j)
                     phiaxissave(l,1,j,i) = phiaxis(l,1,j)
                     phiaxissave(l,2,j,i) = phiaxis(l,2,j)
                  enddo
               enddo
            enddo
         elseif (i .le. nstepssave) then
            lsavephi = .false.
         endif


         if (myid.eq.0) then

c             For debugging, save matrix A and its transpose,
c               the vector b, and the result x
                  if (lsavemat .and. stepcount .eq. saveatstep) then
                     n1 = rshieldingsave
                     n2 = nthused
                     n3 = npsiused
                     do j=1,n3
                        do k=1,n2
                           do l=1,n1
                              jkl = l + (k-1)*n1 + (j-1)*n2*n1
                              bsavevect(jkl) = bsave(l,k,j)
                              xsavevect(jkl) = xsave(l,k,j)
                              do m=1,n3
                                 do n=1,n2
                                    do o=1,n1
                                       mno = o + (n-1)*n1 + (m-1)*n2*n1;
                                       Amat(mno,jkl)
     $                                   = Asave(o,n,m,l,k,j)
                                       Atmat(mno,jkl)
     $                                   = Atsave(o,n,m,l,k,j)
                                    enddo
                                 enddo
                              enddo
                           enddo
                        enddo
                     enddo
                  endif

            
c     Evaluate total charge and z-forces. Not yet implemented for NGP.
c     In sceptic3D, ierad is the outer radius of evaluation, and by
c     default we also evaluate at the probe surface

            call esforce(1,qp,fz,epz,fbz,fx,epx,fbx,fy,epy,fby)
            zmom(i,enccharge,1)=qp
            zmom(i,fieldz,1)=fz
            zmom(i,epressz,1)=epz
            zmom(i,lorentz,1)=fbz
            xmom(i,fieldz,1)=fx
            xmom(i,epressz,1)=epx
            xmom(i,lorentz,1)=fbx
            ymom(i,fieldz,1)=fy
            ymom(i,epressz,1)=epy
            ymom(i,lorentz,1)=fby


            call esforce(ierad,qp,fz,epz,fbz,fx,epx,fbx,fy,epy,fby)
            zmom(i,enccharge,2)=qp
            zmom(i,fieldz,2)=fz
            zmom(i,epressz,2)=epz
            zmom(i,lorentz,2)=fbz
            xmom(i,fieldz,2)=fx
            xmom(i,epressz,2)=epx
            xmom(i,lorentz,2)=fbx
            ymom(i,fieldz,2)=fy
            ymom(i,epressz,2)=epy
            ymom(i,lorentz,2)=fby
          
         endif

c Check if any node tripped lgotooutput, and if so, leave do loop
#ifdef MPI
         ltmp = lgotooutput
         call MPI_ALLREDUCE(lgotooutput,ltmp,1,MPI_INTEGER,MPI_LOR
     $        ,MPI_COMM_WORLD,ierr)
         lgotooutput = ltmp
#endif
         if (lgotooutput) goto 50


c     Main particle advance, including collisions.
         
         call padvnc(dt,icolntype,colnwt,i,maccel,ierad)


c     Reduce back the flux and distribution data from the particle advance.
         call partreduce(i)


c     Adjust to the flux and momenta that would have occurred for
c     standard step size.
         if(myid.eq.0) then
            fluxprobe(i)=fluxprobe(i)/bdtnow
            zmom(i,partz,1)=zmom(i,partz,1)/dt
            zmom(i,partz,2)=zmom(i,partz,2)/dt
            xmom(i,partz,1)=xmom(i,partz,1)/dt
            xmom(i,partz,2)=xmom(i,partz,2)/dt
            ymom(i,partz,1)=ymom(i,partz,1)/dt
            ymom(i,partz,2)=ymom(i,partz,2)/dt
            enertot(i)=enertot(i)/dt

c For debugging, don't do this since clutters output
cc Small mesh written diagnostics.
c            if(nr.eq.10) then
c               write(*,*)'rho,phi,psum'
c               write(*,504)(((rho(iw,jw,kw),iw=1,nr),jw=1,NTHUSED),
c     $              kw=1,NPSIUSED)
c               write(*,504)(((phi(iw,jw,kw),iw=1,nr),jw=0,NTHFULL),
c     $              kw=1,NPSIFULL)
c            endif
         endif

         time=time+dt
 503     format(10f8.1)
 504     format(10f8.3) 

c Check if any node tripped lgotooutput, and if so, leave do loop
#ifdef MPI
         ltmp = lgotooutput
         call MPI_ALLREDUCE(lgotooutput,ltmp,1,MPI_INTEGER,MPI_LOR
     $        ,MPI_COMM_WORLD,ierr)
         lgotooutput = ltmp
#endif
         if (lgotooutput) goto 50

      enddo

 50   continue
c End of particle stepping section.

      itotsteps=maxsteps
      if(myid.eq.0)then
c         write(*,'(/,a)')'Exhausted maxsteps'
         write(*,*)
         write(*,505)nr,nth,npsi,npart,Bz
         write(*,502)debyelen,Ti,maxsteps,dtf,rmax,vd,vprobe
         if(iocprev.ne.npart)write(*,*)'iocprev=',iocprev
         if(rmtoz.ne.1.) write(*,'(''rmtoz='',f10.4)')rmtoz
         if(colnwt.ne.0) write(*,
     $        '(''Collisions: type='',i4,'' density='',f10.5)')
     $        icolntype,colnwt

c Average various fluxes and fields:

         call avefluxes(itotsteps,m2,dt,fave, zmomave,fezave
     $        ,zmoutave,xmomave,fexave,xmoutave,ymomave, feyave,ymoutave
     $        ,qprobeave,epzave,epxave,epyave)

         

         if(lplot.and.finaldiags) then
            write(string,'(''Probe Flux Density='',f8.4)')fave
            call jdrwstr(0.05,0.7,string,1.)
         endif
      endif 

      diags=.true.

      call chargetomesh()
  
      call sumreduce()

      if(myid.eq.0)then
c Write the output files.
         call output(dt,i-1,fave,icolntype,colnwt)
#ifdef HDF
         call outputhdf(dt,i-1,fave,icolntype,colnwt)
#endif
         if (norbits.ge.1) call orbitoutput()
      endif
      if(writepart) call partwrt()

      if(lplot) call pltend()
c Restore the permanent plotting switch.
      lplot=lpstore

c Draw some final diagnostic figures, if flag set.    
      if(finaldiags.and.myid.eq.0)then
         call multiframe(0,0,0)
         call yautoplot(fluxprobe,itotsteps)
         call axlabels('step','Particles to probe in step')
         call pltend()
         call multiframe(1,1,0)
         call rhodisplay(15,rhomin,rhomax)
         if(norbits.gt.0) call plotorbits()
         call pltend()
      endif

#ifdef MPI
      if (myid2.eq.0) then
         write(*,*) "Time : ",MPI_WTIME()-cgtime
      endif
      call MPI_FINALIZE(rc)
#endif

c Use call exit to avoid silly fortran message.
      call exit(0)

 51   continue
c Help section
      write(*,*)
     $ 'Usage: sceptic3D [-t.. -s.. -d.. -c..',
     $     ' -x.. -v.. -p.. -l.. -e[..] -r -w -g -?]'
      write(*,*)' e.g. ./sceptic -v0.2 -nr50 -nt30 -d.05'
      write(*,*)' Switch arguments (defaults):[.ff=float, nnn=integer]'
      write(*,*)' -t.ff,-tn.ff ion and neutr. temp (Ti=1, Tn=Ti);',
     $     ' -snnn max time steps (500);',
     $     ' --bz.ff Larmor frequency'
      write(*,*)' -d.ff time step dt (.025);',
     $     ' -c.ff initial accel factor (1)'
      write(*,*)' -cmax.ff -cmin.ff limits of rho plotting (auto)'
      write(*,*)' -v.ff ion drift vel. rel. to neut. due to E (0.0);',
     $     ' -x.ff max r of domain (5)'
      write(*,*)' -p.ff probe potential (-4.) -pf float -pi insulate;',
     $     ' -b.ff Boundary factor (4)'
      write(*,*) ' -l.ff Debye length (.1) ; -z.ff External Ez (0)'
      write(*,*) ' -onnn Rand Orbits to plot (0) ; -oinnn Aligned',
     $     'orbits (0)'
      write(*,*)' -e[ir,it] distribution function diagnostics',
     $     ' [in cell ir,it](1,1);'
      write(*,*)' -r read particle data (no), -w write ptcl. data (no)'
      write(*,*)' -g<nnn> diag plots only on nnn th step;',
     $     ' -a<n> save plots to disk (pfset n)'
      write(*,*)' -g no diags, -f no final diags -? Print this help.'
      write(*,*)' -nrnnn, -ntnnn, -ninnn, set radial, angle mesh-size',
     $     ', particle number.'
      write(*,*)' -ernnn radius at which to calculate q,E-force.',
     $     ' -mfff ratio of mass to Z (1.)'
      write(*,*)' -ktnnn collision type (0: none, 1 direct, 2 remote).'
     $     ,' -k.ff collision freq.'
      write(*,*)' -vn.ff neutral drift velocity for collisions (0).'
      write(*,*)' --at0 set theta acceleration artificially to zero.'
      write(*,*)' --ap0 set poloidal acc. artificially to zero.'
      write(*,*)' --localinj use local potential for injection.',
     $     ' --smooth symmetrize in angle.'
      write(*,*)' -cd.ff,-cnd.ff,-cB.ff cos(th) ion & neutr dr., and B'
      write(*,*)' -psind.ff psi of neutral drift (pi/2)'
      write(*,*)' --rhoinf.ff specify constant injection rate, rather',
     $     ' than particle number.'
      write(*,*)' --sp use parallel bloc solver, '
      write(*,*)' --nonlin use non-linear Poisson solver, '
      write(*,*)' --bcphi(0) BC potl (0:spherical sym, ',
     $     '1 : Quasi neutrality on outer 15%,',
     $     '    2: Phiout=0, 3: dPhiout/dz=0,',
     $     '    4: dPhiout/dr=-Phiout/r. use with vd=0 and k>0)'
      write(*,*)' --bcr(0) BC reinject (0:spherical sym, 1: Simple',
     $     'Maxwellian, 2: Adiabatic, 3: Monte Carlo)'
      write(*,*)' --mcrnt<nnn>(5), --mcrnp<nnn>(5), MC inj. mesh size.'
      write(*,*)' --mcrni<nnn>(100000) MC inj. number of particles.'
      write(*,*)' --subcyc use step subcycling near probe.'
      write(*,*)' -ver Old verlet integrator.',  
     $     '-bohm Impose Bohm condition when LDe=0.'
      write(*,*)' --fulloutput do not suppress output of large arrays.'

      end

c***********************************************************************
c End of Main Program.
c***********************************************************************

c Calculate rho from the psum.
      subroutine rhocalc(lsmoothT,lsmoothP,i,dt)

c If lsmooth, then symmetrize the density, by averaging.
      logical lsmoothT,lsmoothP
c Current time-step
      integer i
c Common data:
      include 'piccom.f'
      include 'errcom.f'
 501  format(10f8.1)

c Now we have added up fractional charges assigned to each mesh point. 
c We need to divide by the volume corresponding to each.
      if(.not.rhoinf.gt.1.e-4) then
c Initialize rhoinf approximately:
         write(*,*)'Rhoinf in rhocalc too small',rhoinf
         rhoinf=numprocs*npart/(4.*pi*r(NRUSED)**3/3.)
         write(*,*)'Rhoinf in rhocalc approximated as',rhoinf
      endif
      do ipsi=1,NPSIUSED
         do ith=1,NTHUSED
            do ir=1,NRUSED
c     Volumes are now calculated in meshinit.
               if(psum(ir,ith,ipsi).eq.0. .and. debyelen.eq.0.)then
c     write(*,*)'Rhocalc fixup: zero particles in cell',ir,ith,ipsi
                  rho(ir,ith,ipsi)=
     $                 0.01*volinv(ir)*(nthused-1.)*npsiused*np/rhoinf
               else
                  rho(ir,ith,ipsi)=psum(ir,ith,ipsi)
     $                 *volinv(ir)*(nthused -1.)*npsiused*np/rhoinf
               endif
            enddo
         enddo
      enddo
      do ipsi=1,NPSIUSED
         do ir=1,NRUSED
c Fix the theta boundaries
            rho(ir,1,ipsi)=2.*rho(ir,1,ipsi)
            rho(ir,nthused,ipsi)=2.*rho(ir,nthused,ipsi)

c     Smoothing in theta angle
            if(lsmoothT) then
               rhotot=0.
               do ith=1,NTHUSED
                  rhotot=rhotot+rho(ir,ith,ipsi)
               enddo
               rhotot=rhotot/nthused
               do ith=1,NTHUSED
                  rho(ir,ith,ipsi)=rhotot
               enddo
            endif
         enddo
      enddo
c     Smoothing in psi angle
      if(lsmoothP) then
         do ith=1,NTHUSED
            do ir=1,NRUSED
               rhotot=0.
               do ipsi=1,NPSIUSED
                  rhotot=rhotot+rho(ir,ith,ipsi)
               enddo
               rhotot=rhotot/npsiused
               do ipsi=1,NPSIUSED
                  rho(ir,ith,ipsi)=rhotot
               enddo
            enddo
         enddo
      endif
         
c     If vs<0 (which is always the case unless there is exceptional
c     noise), get the density at r=1 by dividing the ion flux by the
c     average velocity
         do ipsi=1,npsi
            do ith=1,nth
               vs=vrincellave(ith,ipsi)/(1e-6+fincellave(ith,ipsi))
               fluxofangle=fincellave(ith,ipsi)*(nth-1)*npsi/(4*pi
     $              *dt*r(1)**2)
               if(vs.lt.0) then
                  rho(1,ith,ipsi)=-fluxofangle/(vs*rhoinf)+1e-6
                  if(ith.eq.1.or.ith.eq.nth) then
                     rho(1,ith,ipsi)=2*rho(1,ith,ipsi)
                  endif
               endif
            enddo
         enddo
         
c     Storing the averaged diagnostic density.  Because densities at
c     following time-steps are correlated, just save every 4
c     time-steps. This in effect multiplies by 4 the length of
c     averaging.

      if(i.ge.maxsteps-diagsamp*nstepsave) then
         if(mod(i-(maxsteps-diagsamp*nstepsave),diagsamp).eq.0)then 
c     past=0 for its first calculation. This complicated averaging
c     simply ensures that at the end, every rho at each time step has
c     equal weight
            past=1.*(i-maxsteps)/diagsamp+nstepsave
            do i3=1,npsiused
               do i2=1,nthused
                  do i1=1,nrused

                     rhoDiag(i1,i2,i3)=(rhoDiag(i1,i2,i3)*past +rho(i1
     $                    ,i2,i3))/(past+1)
                     pDiag(i1,i2,i3)=(pDiag(i1,i2,i3)*past +psum(i1,i2
     $                    ,i3))/(past+1)
                     vrDiag(i1,i2,i3)=(vrDiag(i1,i2,i3)*past +vrsum(i1
     $                    ,i2,i3))/(past+1)
                     vtDiag(i1,i2,i3)=(vtDiag(i1,i2,i3)*past +vtsum(i1
     $                    ,i2,i3))/(past+1)
                     vpDiag(i1,i2,i3)=(vpDiag(i1,i2,i3)*past +vpsum(i1
     $                    ,i2,i3))/(past+1)
                     vr2Diag(i1,i2,i3)=(vr2Diag(i1,i2,i3)*past
     $                    +vr2sum(i1,i2,i3))/(past+1)
                     vt2Diag(i1,i2,i3)=(vt2Diag(i1,i2,i3)*past
     $                    +vt2sum(i1,i2,i3))/(past+1)
                     vp2Diag(i1,i2,i3)=(vp2Diag(i1,i2,i3)*past
     $                    +vp2sum(i1,i2,i3))/(past+1)
                     vrtDiag(i1,i2,i3)=(vrtDiag(i1,i2,i3)*past
     $                    +vrtsum(i1,i2,i3))/(past+1)
                     vrpDiag(i1,i2,i3)=(vrpDiag(i1,i2,i3)*past
     $                    +vrpsum(i1,i2,i3))/(past+1)
                     vtpDiag(i1,i2,i3)=(vtpDiag(i1,i2,i3)*past
     $                    +vtpsum(i1,i2,i3))/(past+1)
                     
                  enddo
               enddo
            enddo
         endif
      endif

      end
 
c*********************************************************************

      subroutine sumreduce()
#ifdef MPI
c The particle number total
      include 'piccom.f'
      include 'errcom.f'
      include 'mpif.h'
      logical first

      character*11 fil

      real ptot(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vrtot(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vttot(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vptot(1:nrsize-1,1:nthsize-1,1:npsisize-1)

      real vr2tot(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vt2tot(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vp2tot(1:nrsize-1,1:nthsize-1,1:npsisize-1)

      real vrttot(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vrptot(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vtptot(1:nrsize-1,1:nthsize-1,1:npsisize-1)

      real vxtot(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vytot(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vztot(1:nrsize-1,1:nthsize-1,1:npsisize-1)

      real currtot(4)


c     psum is always required for the density calculation. Reduce the
c     full array except in the third dimension (psi)

 
      call MPI_REDUCE(psum,ptot,(nrsize-1)*(nthsize-1)*npsiused,
     $     MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

c     curr is the ion current in the domain, necessary to calculate the
c     Lorentz force. Altough it is only needed for i>m2 (last part of
c     the run), it is only a 1*4 array so we reduce it at each timestep
      call MPI_REDUCE(curr,currtot,4,
     $     MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

c     If diags is off, we sumreduce the information needed to calculate
c     the sound speed at the probe edge, needed only in the zero debye
c     length case. vpsum and vtsum are needed for the Lorentz force
c     calculation


      if(diags) then
c         write(*,*) vxsum,vysum,vzsum
         call MPI_REDUCE(vxsum,vxtot,(nrsize-1)*(nthsize-1)*npsiused
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(vysum,vytot,(nrsize-1)*(nthsize-1)*npsiused
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)         
         call MPI_REDUCE(vzsum,vztot,(nrsize-1)*(nthsize-1)*npsiused
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      endif


      if(diags.or.samp.or.debyelen.eq.0) then
         call MPI_REDUCE(vrsum,vrtot,(nrsize-1)*(nthsize-1)*npsiused
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(vr2sum,vr2tot,(nrsize-1)*(nthsize-1)*npsiused
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      endif


      if(diags.or.samp)then
         call MPI_REDUCE(vpsum,vptot,(nrsize-1)*(nthsize-1)*npsiused
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(vtsum,vttot,(nrsize-1)*(nthsize-1)*npsiused
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(vt2sum,vt2tot,(nrsize-1)*(nthsize-1)*npsiused
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(vp2sum,vp2tot,(nrsize-1)*(nthsize-1)*npsiused
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(vrtsum,vrttot,(nrsize-1)*(nthsize-1)*npsiused
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(vrpsum,vrptot,(nrsize-1)*(nthsize-1)*npsiused
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(vtpsum,vtptot,(nrsize-1)*(nthsize-1)*npsiused
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)        
      endif
      if(myid.eq.0)then

         do i1=1,4
            curr(i1)=currtot(i1)
         enddo

         do i3=1,npsiused
            do i2=1,nthused
               do i1=1,nrused
                  psum(i1,i2,i3)=ptot(i1,i2,i3)
               enddo
            enddo
         enddo

         if(diags.or.samp.or.debyelen.eq.0) then
            do i3=1,npsiused
               do i2=1,nthused
                  do i1=1,nrused
                     vrsum(i1,i2,i3)=vrtot(i1,i2,i3)
                     vr2sum(i1,i2,i3)=vr2tot(i1,i2,i3)
                  enddo
               enddo
            enddo
         endif
         if(diags) then
            do i3=1,npsiused
               do i2=1,nthused
                  do i1=1,nrused
                     vxsum(i1,i2,i3)=vxtot(i1,i2,i3)
                     vysum(i1,i2,i3)=vytot(i1,i2,i3)
                     vzsum(i1,i2,i3)=vztot(i1,i2,i3)
                  enddo
               enddo
            enddo
         endif
         if(diags.or.samp)then
            do i3=1,npsiused
               do i2=1,nthused
                  do i1=1,nrused
                     vt2sum(i1,i2,i3)=vt2tot(i1,i2,i3)
                     vp2sum(i1,i2,i3)=vp2tot(i1,i2,i3)
                     vrtsum(i1,i2,i3)=vrttot(i1,i2,i3)
                     vrpsum(i1,i2,i3)=vrptot(i1,i2,i3)
                     vtpsum(i1,i2,i3)=vtptot(i1,i2,i3)
                     vpsum(i1,i2,i3)=vptot(i1,i2,i3)
                     vtsum(i1,i2,i3)=vttot(i1,i2,i3)
                  enddo
               enddo
            enddo
         endif

      endif

#endif
      end
c***********************************************************************
      subroutine partreduce(i)
      integer i
      include 'piccom.f'
      include 'errcom.f'
#ifdef MPI
      include 'mpif.h'
      real nvdiagtot(nvmax)
         call MPI_REDUCE(nrein,nreintot,1,MPI_INTEGER,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(nreintry,nreintrytot,1,MPI_INTEGER,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(spotrein,spotreintot,1,MPI_REAL,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(fluxrein,fluxreintot,1,MPI_REAL,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(ninner,nintot,1,MPI_INTEGER,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(nincell,nincellstep(1,1,i),nthsize*npsisize
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(vrincell,vrincellstep(1,1,i),nthsize*npsisize
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(vr2incell,vr2incellstep(1,1,i),nthsize*npsisize
     $        ,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(zmomprobe,zmom(i,partz,1),1,MPI_REAL,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(xmomprobe,xmom(i,partz,1),1,MPI_REAL,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(ymomprobe,ymom(i,partz,1),1,MPI_REAL,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(zmout,zmom(i,partz,2),1,MPI_REAL,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(xmout,xmom(i,partz,2),1,MPI_REAL,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(ymout,ymom(i,partz,2),1,MPI_REAL,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(enerprobe,enertot(i),1,MPI_REAL,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         if(diags)then
            call MPI_REDUCE(nvdiag,nvdiagtot,nvmax,MPI_REAL,
     $           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if(myid.eq.0) then
               do kk=1,nvmax
                  nvdiag(kk)=nvdiagtot(kk)
               enddo
            endif
         endif
         if(ldist)then
            call MPI_REDUCE(vrdiagin,nvdiagtot,nvmax,MPI_REAL,
     $           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if(myid.eq.0) then
               do kk=1,nvmax
                  vrdiagin(kk)=nvdiagtot(kk)
               enddo
            endif
            call MPI_REDUCE(vtdiagin,nvdiagtot,nvmax,MPI_REAL,
     $           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if(myid.eq.0) then
               do kk=1,nvmax
                  vtdiagin(kk)=nvdiagtot(kk)
               enddo
            endif
         endif
#else
c This shuffle is necessary to accommodate the reduce.
         nreintot=nrein
         nreintrytot=nreintry
         spotreintot=spotrein
         fluxreintot=fluxrein
         nintot=ninner
         do j=1,nth
            do k=1,npsi
               nincellstep(j,k,i)=nincell(j,k)
               vrincellstep(j,k,i)=vrincell(j,k)
               vr2incellstep(j,k,i)=vr2incell(j,k)
            enddo
         enddo
         zmom(i,partz,1)=zmomprobe
         zmom(i,partz,2)=zmout
         xmom(i,partz,1)=xmomprobe
         xmom(i,partz,2)=xmout
         ymom(i,partz,1)=ymomprobe
         ymom(i,partz,2)=ymout
         enertot(i)=enerprobe
#endif
c         write(*,*)'zmomstep=',zmomstep(i),'  zoutstep=',zmoutstep(i)
         nrein=nreintot
         nreintry=nreintrytot
         spotrein=spotreintot
         fluxrein=fluxreintot
         fluxprobe(i)=nintot
         ninner=nintot
         end
c***************************************************************
      subroutine aveupstep()
      include 'piccom.f'
      include 'errcom.f'
c Increment the averaging processes only here.
      if(nstepsave.lt.nsamax)nstepsave=nstepsave+1
      do j=1,nvmax
         nvdiag(j)=nvdiag(j)*(nstepsave-1.)/nstepsave
      enddo
      do j=1,nvmax
         vtdiagin(j)=vtdiagin(j)*(nstepsave-1.)/nstepsave
         vrdiagin(j)=vrdiagin(j)*(nstepsave-1.)/nstepsave
      enddo
      end
c*****************************************************************
      subroutine avefluxes(itotsteps,m2,dt,fave, zmomave,fezave
     $     ,zmoutave,xmomave,fexave,xmoutave,ymomave, feyave,ymoutave
     $     ,qprobeave,epzave,epxave,epyave)
      include 'piccom.f'
      include 'errcom.f'
      
      integer m2
c     Average the flux to the probe from the few steps (given by m2).
c     Also average the momenta components and place them in
c     zmom(nstepmax,*,*). maxsteps is limited to be nstepmax-1.  The
c     average is what is output.
      fave=0.
      zmomave=0.
      zmoutave=0.
      zmoutvar=0.
      xmomave=0.
      xmoutave=0.
      xmoutvar=0.
      ymomave=0.
      ymoutave=0.
      ymoutvar=0.
      fezave=0.
      fexave=0.
      feyave=0.
      qprobeave=0.
      epzave=0.
      epxave=0.
      epyave=0.

      do j=2,5
         zmom(nstepmax,j,1)=0.
         zmom(nstepmax,j,2)=0.
         xmom(nstepmax,j,1)=0.
         xmom(nstepmax,j,2)=0.
         ymom(nstepmax,j,1)=0.
         ymom(nstepmax,j,2)=0.
      enddo
      zmom(nstepmax,1,1)=0.
      zmom(nstepmax,1,2)=0.


      enertot(nstepmax)=0.
      if(itotsteps.ge.2)then
         do k=m2,itotsteps
            fave=fave+fluxprobe(k)
            do j=2,5
               zmom(nstepmax,j,1)=zmom(nstepmax,j,1)+zmom(k,j,1)
               zmom(nstepmax,j,2)=zmom(nstepmax,j,2)+zmom(k,j,2)
               xmom(nstepmax,j,1)=xmom(nstepmax,j,1)+xmom(k,j,1)
               xmom(nstepmax,j,2)=xmom(nstepmax,j,2)+xmom(k,j,2)
               ymom(nstepmax,j,1)=ymom(nstepmax,j,1)+ymom(k,j,1)
               ymom(nstepmax,j,2)=ymom(nstepmax,j,2)+ymom(k,j,2)
            enddo
            zmom(nstepmax,1,1)=zmom(nstepmax,1,1)+zmom(k,1,1)
            zmom(nstepmax,1,2)=zmom(nstepmax,1,2)+zmom(k,1,2)

            enertot(nstepmax)=enertot(nstepmax)+enertot(k)

            zmomave=zmomave+zmom(k,partz,1)
            fezave=fezave+zmom(k,fieldz,1)
            zmoutave=zmoutave+zmom(k,partz,2)
            zmoutvar=zmoutvar+zmom(k,partz,2)**2
            qprobeave=qprobeave+zmom(k,enccharge,1)
            epzave=epzave+zmom(k,epressz,2)

            xmomave=xmomave+xmom(k,partz,1)
            fexave=fexave+xmom(k,fieldz,1)
            xmoutave=xmoutave+xmom(k,partz,2)
            xmoutvar=xmoutvar+xmom(k,partz,2)**2
            epxave=epxave+xmom(k,epressz,2)

            ymomave=ymomave+ymom(k,partz,1)
            feyave=feyave+ymom(k,fieldz,1)
            ymoutave=ymoutave+ymom(k,partz,2)
            ymoutvar=ymoutvar+ymom(k,partz,2)**2
            epyave=epyave+ymom(k,epressz,2)
         enddo
         fave=fave/(itotsteps-m2+1)
c     normalize to rhoinf and the probe area to give normalized flux density
         fave=fave/(4.*pi*r(1)**2)/rhoinf/dt

         do j=2,5
            zmom(nstepmax,j,1)=zmom(nstepmax,j,1)/(itotsteps-m2+1)
            zmom(nstepmax,j,2)=zmom(nstepmax,j,2)/(itotsteps-m2+1)
            xmom(nstepmax,j,1)=xmom(nstepmax,j,1)/(itotsteps-m2+1)
            xmom(nstepmax,j,2)=xmom(nstepmax,j,2)/(itotsteps-m2+1)
            ymom(nstepmax,j,1)=ymom(nstepmax,j,1)/(itotsteps-m2+1)
            ymom(nstepmax,j,2)=ymom(nstepmax,j,2)/(itotsteps-m2+1)
         enddo
         zmom(nstepmax,1,1)=zmom(nstepmax,1,1)/(itotsteps-m2+1)
         zmom(nstepmax,1,2)=zmom(nstepmax,1,2)/(itotsteps-m2+1)

         j=lorentz

         zmom(nstepmax,j,1)=zmom(nstepmax,j,1)/rhoinf
         zmom(nstepmax,j,2)=zmom(nstepmax,j,2)/rhoinf
         xmom(nstepmax,j,1)=xmom(nstepmax,j,1)/rhoinf
         xmom(nstepmax,j,2)=xmom(nstepmax,j,2)/rhoinf
         ymom(nstepmax,j,1)=ymom(nstepmax,j,1)/rhoinf
         ymom(nstepmax,j,2)=ymom(nstepmax,j,2)/rhoinf

         j=partz
         
         zmom(nstepmax,j,1)=zmom(nstepmax,j,1)/rhoinf
         zmom(nstepmax,j,2)=zmom(nstepmax,j,2)/rhoinf
         xmom(nstepmax,j,1)=xmom(nstepmax,j,1)/rhoinf
         xmom(nstepmax,j,2)=xmom(nstepmax,j,2)/rhoinf
         ymom(nstepmax,j,1)=ymom(nstepmax,j,1)/rhoinf
         ymom(nstepmax,j,2)=ymom(nstepmax,j,2)/rhoinf

         zmomave=zmomave/(itotsteps-m2+1)/(rhoinf)
         fezave=fezave/(itotsteps-m2+1)
         xmomave=xmomave/(itotsteps-m2+1)/(rhoinf)
         fexave=fexave/(itotsteps-m2+1)
         ymomave=ymomave/(itotsteps-m2+1)/(rhoinf)
         feyave=feyave/(itotsteps-m2+1)

         qprobeave=qprobeave/(itotsteps-m2+1)
         zmoutave=zmoutave/(itotsteps-m2+1)/(rhoinf)
         zmoutvar=zmoutvar/(itotsteps-m2+1)/(rhoinf**2)-zmoutave**2
         epzave=epzave/(itotsteps-m2+1)
         xmoutave=xmoutave/(itotsteps-m2+1)/(rhoinf)
         xmoutvar=xmoutvar/(itotsteps-m2+1)/(rhoinf**2)-xmoutave**2
         epxave=epxave/(itotsteps-m2+1)
         ymoutave=ymoutave/(itotsteps-m2+1)/(rhoinf)
         ymoutvar=ymoutvar/(itotsteps-m2+1)/(rhoinf**2)-ymoutave**2
         epyave=epyave/(itotsteps-m2+1)

         enertot(nstepmax)=
     $        enertot(nstepmax)/(itotsteps-m2+1)/rhoinf

         write(*,501)fave
 501     format('Probe flux density=',f8.4)
         write(*,503)enertot(nstepmax)
 503     format('Probe energy flux=',f8.4)
         ztotal1=zmom(nstepmax,fieldz,1)*debyelen**2
     $        +zmom(nstepmax,epressz,1)+zmom(nstepmax,partz,1)
     $        +zmom(nstepmax,lorentz,1)
         ztotal2=zmom(nstepmax,fieldz,2)*debyelen**2
     $        +zmom(nstepmax,epressz,2)+zmom(nstepmax,partz,2)
     $        +zmom(nstepmax,lorentz,2)
         xtotal1=xmom(nstepmax,fieldz,1)*debyelen**2
     $        +xmom(nstepmax,epressz,1)+xmom(nstepmax,partz,1)
     $        +xmom(nstepmax,lorentz,1)
         xtotal2=xmom(nstepmax,fieldz,2)*debyelen**2
     $        +xmom(nstepmax,epressz,2)+xmom(nstepmax,partz,2)
     $        +xmom(nstepmax,lorentz,2)
         ytotal1=ymom(nstepmax,fieldz,1)*debyelen**2
     $        +ymom(nstepmax,epressz,1)+ymom(nstepmax,partz,1)
     $        +ymom(nstepmax,lorentz,1)
         ytotal2=ymom(nstepmax,fieldz,2)*debyelen**2
     $        +ymom(nstepmax,epressz,2)+ymom(nstepmax,partz,2)
     $        +ymom(nstepmax,lorentz,2)
         write(*,*)'  z E-field   z Electrons',
     $        '   z Ions     z Lorentz    z Total    Charge'
         write(*,506)(zmom(nstepmax,j,1),j=2,5)
     $        ,ztotal1,zmom(nstepmax,1,1)
         write(*,506)(zmom(nstepmax,j,2),j=2,5)
     $        ,ztotal2,zmom(nstepmax,1,2)
         write(*,*) ""
         write(*,*)'  x E-field   x Electrons',
     $        '   x Ions     x Lorentz    x Total'
         write(*,507)(xmom(nstepmax,j,1),j=2,5),xtotal1
         write(*,507)(xmom(nstepmax,j,2),j=2,5),xtotal2

         write(*,*) ""
         write(*,*)'  y E-field   y Electrons',
     $        '   y Ions     y Lorentz    y Total'
         write(*,507)(ymom(nstepmax,j,1),j=2,5),ytotal1
         write(*,507)(ymom(nstepmax,j,2),j=2,5),ytotal2

 506     format(6f12.6)
 507     format(5f12.6)

         phiout=0.
         do j=1,NTHUSED
            do k=1,NPSIUSED
               phiout=phi(NRUSED,j,k)+phiout
            enddo
         enddo
         phiout=phiout/NTHUSED/NPSIUSED
 
         write(*,*) 'rhoinf', rhoinf
         write(*,*) 'Outer potential',phiout


c tests
c         write(*,*)' Steps zmominner,outer'
c         do k=1,itotsteps
c            write(*,*)k,zmom(k,partz,1)/rhoinf,zmom(k,partz,2)/rhoinf
c         enddo
         write(*,'(''Ion Forces: Inner, outer, uncertainty,'',2f10.4,
     $        ''+-'',f7.4)')
     $        zmomave,zmoutave,sqrt(zmoutvar/(itotsteps-m2+1))
      endif
      end



