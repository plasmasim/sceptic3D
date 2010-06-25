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


c Advance the particles
      subroutine padvnc(dtin,icolntype,colnwt,step,maccel,ierad)

      integer step,maccel
      real dt,dtin
c Common data:
      include 'piccom.f'
      include 'errcom.f'
      include 'colncom.f'

      real accel(3)
      real rn
      real cosomdt,sinomdt
c temp data:
      real temp
      integer idum
      integer isubcycle
      logical lcollide,lcstep
      real ctc,spsi,cpsi,rad,sd,sB


c Choose the collision cycle here and set tau appropriately:
c icycle of 1 costs about 20% extra cf false. (Mostly alog, I'd guess).
      if(colnwt.gt.0.)then
c         icycle=1./(20.*colnwt*dtin)
c         icycle=1./(50.*colnwt*dtin)
         icycle=1
         if(.not.icycle.ge.1) icycle=1
         ichoose=0
c         ichoose=ran0(idum)*icycle
         tau=1./(colnwt*icycle)
         lcollide=.true.
      else
         lcollide=.false.
         icycle=1
         ichoose=1
         tau=1.e20
      endif

      idum=1

c Xp is the three x-coordinates followed by the 3 v coordinates.
c Use a leapfrog scheme, so interpret the v-coords as half a step
c behind the x-coords. 
      tisq=sqrt(Ti)
c If lsubcycle, use multiple fractional steps near inner boundary.
      dt=dtin
      rp2=r(1)**2
c Zero the sums.
      ncollide=0
      nrein=0
      nreintry=0
      spotrein=0.
      ninner=0
      fluxrein=0.
      ntrapre=0
      zmomprobe=0.
      xmomprobe=0.
      ymomprobe=0.
      enerprobe=0.
      zmout=0.
      xmout=0.
      ymout=0.
      iocthis=0.

c Velocity and magnetic axis angle sinus
      sd=sqrt(1-cd**2)
      sB=sqrt(1-cB**2)

      do j=1,nth
         do k=1,npsi
            nincell(j,k)=0
            vrincell(j,k)=0
            vr2incell(j,k)=0
         enddo
      enddo
 
c Set the sum of particles velocities to zero
      do k=1,4
         curr(k)=0;
      enddo

      ido=npart

c      write(*,*)'colnwt,tau,Eneutral,icycle',colnwt,tau,Eneutral,icycle
c End of setup
c------------------ Iterate over particles --------------------------
c No-subcycle default. Never gets changed w/o subcycling.
      dts=dtin
      isubcycle=1
      do i=1,ido

         if(ipf(i).gt.0) then
c ````````````````````````````````````````` Treatment of active slot.
c     Find the mesh position and the trigonometry.
c     Here we do need half quantities.
            ih=1
            hf=88.
            call ptomesh(i,il,rf,ith,tf,ipl,pf,st,ct,sp,cp,rp
     $           ,zetap,ih,hf)
c           If lgotooutput tripped, skip to end
            if (lgotooutput) then
               write (*,*) 'lgotooutput tripped in padvnc.'
               goto 401
            endif

            
c .................... Subcycle Loop .................
            remdt=dtin
            ic=0
            lcstep=.false.
c            do 81 ic=1,isubcycle    Obsolete.
c Here is the start of the modified loop, now explicit.
c We iterate till we have used up the whole time step dtin (remdt=0).
c Steps may be shortened by subcycling and collisions.
 80         ic=ic+1

c  Now we know where we are in radius rp. 
c  We decide the level of subcycling.
            if(lsubcycle) then

c This seems not to be useful
c               if(step.gt.maccel)then
c Require the time-step to be smaller than 10% of the Larmor period
c                  isubcycle=max(r(nrfull)/rp,dtin/(0.628/(Bz+1e-6)))
c               else
                  isubcycle=r(nrfull)/rp
c               endif

c          if(mod(i,1000).eq.0) write(*,'(i1,$)')isubcycle

               dts=dtin/isubcycle*1.00001
               
            endif

c If prior step was ended by a collision, restart the particle velocity.

            if(lcstep)then
               call postcollide(i,tisq)
               lcstep=.false.
c     Because postcollide selects the velocity and the position at the
c     same time, we need to set dtprec to zero, in order to offset v and
c     x by half a time step properly.
               dtprec(i)=0
            endif
            dt=min(dts,remdt)

            if(lcollide .and. mod(i,icycle).eq.ichoose)then
c Here we calculate the time to next collision: cdt
c Based on random number draw and poisson distribution.
               cdt= -alog(ran0(idum))*tau
c Using this approximation instead saves negligible time:
c                  cdt= ran0(idum)*tau
c So I conclude that the only loss of time is initialization.
               if(cdt.lt.dt)then
c Collision at the end of cdt step.
                  dt=cdt
                  lcstep=.true.
                  ncollide=ncollide+1
               endif
            endif
            if(.not.dt.lt.1000.)then
c     Error trap 
               write(*,*)'dt error: dtin, dt, cdt, dts, remdt:',
     $              dtin, dt, cdt, dts, remdt
               write(*,*)'isubcycle, r(nrfull), rp:',
     $              isubcycle, r(nrfull), rp
c              Trigger go to output
               lgotooutput = .true.
               goto 401
            endif
            remdt=remdt-dt
            
c Except for the first time, find new position.
            if(ic.ne.1)then 
               ih=1
               hf=77.

               
               
               call ptomesh(i,il,rf,ith,tf,ipl,pf,st,ct,
     $              sp,cp,rp,zetap,ih,hf)
c              If lgotooutput tripped, skip to end
               if (lgotooutput) then
                  write (*,*) 'lgotooutput tripped in padvnc.'
                  goto 401
               endif
            endif            
            call getaccel(i,accel,il,rf,ith,tf,ipl,pf,st,ct,
     $           sp,cp,rp,zetap,ih,hf)
c           If lgotooutput tripped, skip to end
            if (lgotooutput) then
               write (*,*) 'lgotooutput tripped in padvnc.'
               goto 401
            endif


c For acceleration, when dt is changing, use the average of prior and
c present values: dtnow.
c               if(dtprec(i).eq.0.)dtprec(i)=dt

            dtnow=0.5*(dt+dtprec(i))

c     The first way to deal with E^B drift is with a convective
c     field. Use an other approach
c     accel(1)=accel(1)-vd*sd*Bz

            if(.not.verlet) then 

c     Getaccel returns the accel based on the charge-field calculation.
               accel(3)=accel(3)+Eneut(3)
               accel(2)=accel(2)+Eneut(2)
               accel(1)=accel(1)+Eneut(1)

c     Kick
               do j=4,6
                  xp(j,i)=xp(j,i)+accel(j-3)*dtnow
               enddo
c     Drift with Bz.ne.0 (Cyclotronic integrator) Perpendicular direction
               if(Bz.ne.0) then
c ------

c Account for the E*B drift by working in a frame where Econvective=0
                  xp(5,i)=xp(5,i)-vd*sd
                  xp(6,i)=xp(6,i)-vd*cd
c B is not aligned with the z-axis
                  if(cB.lt.0.999) then     
                     temp=xp(2,i)
                     xp(2,i)=temp*cB-xp(3,i)*sB
                     xp(3,i)=xp(3,i)*cB+temp*sB 
                     temp=xp(5,i)
                     xp(5,i)=temp*cB-xp(6,i)*sB
                     xp(6,i)=xp(6,i)*cB+temp*sB 
                  endif
c ------
                  cosomdt=cos(Bz*dt)
                  sinomdt=sin(Bz*dt)
                  xp(1,i)=xp(1,i)+
     $                 (xp(5,i)*(1-cosomdt)+xp(4,i)*sinomdt)/Bz
                  xp(2,i)=xp(2,i)+
     $                 (xp(4,i)*(cosomdt-1)+xp(5,i)*sinomdt)/Bz
                  
                  temp=xp(4,i)
                  xp(4,i)=temp*cosomdt+xp(5,i)*sinomdt
                  xp(5,i)=xp(5,i)*cosomdt-temp*sinomdt

                  xp(3,i)=xp(3,i)+xp(6,i)*dt

c ------
c B is not aligned with the z-axis (Rotate back)
                  if(cB.lt.0.999) then
                     temp=xp(2,i)
                     xp(2,i)=temp*cB+xp(3,i)*sB
                     xp(3,i)=xp(3,i)*cB-temp*sB 
                     temp=xp(5,i)
                     xp(5,i)=temp*cB+xp(6,i)*sB
                     xp(6,i)=xp(6,i)*cB-temp*sB 
                  endif
c Account for the E*B drift (Transform back)
                  xp(2,i)=xp(2,i)+vd*sd*dt
                  xp(3,i)=xp(3,i)+vd*cd*dt
                  xp(5,i)=xp(5,i)+vd*sd
                  xp(6,i)=xp(6,i)+vd*cd
c ------

               else
                  do j=1,3
                     xp(j,i)=xp(j,i)+xp(j+3,i)*dt
                  enddo
               endif
               
            else
c     Old Verlet integrator

c     Getaccel returns the accel based on the charge-field calculation.
               accel(3)=accel(3)+Eneut(3)
               accel(2)=accel(2)+Eneut(2)
               accel(1)=accel(1)+Eneut(1)
               
               if(Bz.eq.0.)then
c     Don't use split steps if Bz=0, for speed gain of 9%.
                  do j=4,6
                     xp(j,i)=xp(j,i)+accel(j-3)*dtnow
                  enddo
                  do j=1,3
                     xp(j,i)=xp(j,i)+xp(j+3,i)*dt
                  enddo
                  
               else
c Old Boris integrator

c First half of velocity advance:    AccelPhi/2+AccelBz+AccelPhi/2
                  do j=4,6
                     xp(j,i)=xp(j,i)+accel(j-3)*dtnow/2
                  enddo
c B-field rotation
c ------

c Account for the E*B drift
                  xp(5,i)=xp(5,i)-vd*sd
                  xp(6,i)=xp(6,i)-vd*cd
c B is not aligned with the z-accis
                  if(cB.lt.0.999) then   
                     temp=xp(5,i)
                     xp(5,i)=temp*cB-xp(6,i)*sB
                     xp(6,i)=xp(6,i)*cB+temp*sB 
                  endif
c ------
                  cosomdt=cos(Bz*dtnow)
                  sinomdt=sin(Bz*dtnow)         
                  temp=xp(4,i)
                  xp(4,i)=temp*cosomdt+xp(5,i)*sinomdt
                  xp(5,i)=xp(5,i)*cosomdt-temp*sinomdt
c ------
c B is not aligned with the z-accis (Rotate back)
                  if(cB.lt.0.999) then
                     temp=xp(5,i)
                     xp(5,i)=temp*cB+xp(6,i)*sB
                     xp(6,i)=xp(6,i)*cB-temp*sB 
                  endif
c Account for the E*B drift (Transform back)
                  xp(5,i)=xp(5,i)+vd*sd
                  xp(6,i)=xp(6,i)+vd*cd
c ------

c Second half of velocity advance
                  do j=4,6
                     xp(j,i)=xp(j,i)+accel(j-3)*dtnow/2
                  enddo
                     
                  do j=1,3
                     xp(j,i)=xp(j,i)+xp(j+3,i)*dt
                  enddo
               endif
               
            endif

            dtprec(i)=dt
            rn2=0.
            xdv=0.
            v2=0.

c Position advance
            do j=1,3
               rn2=rn2+xp(j,i)**2
               xdv=xdv+xp(j,i)*xp(j+3,i)
               v2=v2+xp(j+3,i)**2
            enddo
           
c The time prior to step end of closest approach
            tm=xdv/v2
            rn=sqrt(rn2)

c The following does not make a difference almost)

c  Test if we went through the probe and came back out.
            if((0..lt.tm .and. tm.lt.dt .and.
     $           (rn2 - tm**2*v2).lt.rp2))then
c For a long time this had an error: used  tm**2/v2 erroneously. 
c Corrected 9 Apr 07.
               if(rn.gt.r(1))then
c     write(*,*)'Through probe',tm,(rn2 - tm**2*v2)
                  rn=0.
               endif
            endif

c-----------------------------------------------------------------               
c Handling boundaries :

            if(rn.le.r(1).or.rn.ge.r(nr))then
c     We left the computational domain
               
c     Rewind time to get the right velocities at the domain exit. Assume
c     that the particle left the domain at a random time during the last
c     time-step.
               dtl=-ran0(idum)*dt



c     Don't need the electrostatic part of the correction, since it
c     averages to zero (the first dt/2 is to finish the step since
c     leapfrog is offset half a timestep)

c     do j=4,6
c         xp(j,i)=xp(j,i)+accel(j-3)*(dt/2+dtl)
c     enddo

               if(Bz.ne.0) then

c     For strong magnetic fields, the convective Efield may be strong.
c     In order to get accurate force calculations, we must account for
c     the fact that the particle may have left the domain between 0 and
c     dt before now. 
               
              
c     Old verlet integrator here because we only really care about the
c     magnetic field effect

c Account for the E*B drift
                  xp(5,i)=xp(5,i)-vd*sd
                  xp(6,i)=xp(6,i)-vd*cd
c B is not aligned with the z-accis
                  if(cB.lt.0.999) then   
                     temp=xp(5,i)
                     xp(5,i)=temp*cB-xp(6,i)*sB
                     xp(6,i)=xp(6,i)*cB+temp*sB 
                  endif
c     ------
                  cosomdt=cos(Bz*dtl)
                  sinomdt=sin(Bz*dtl)         
                  temp=xp(4,i)
                  xp(4,i)=temp*cosomdt+xp(5,i)*sinomdt
                  xp(5,i)=xp(5,i)*cosomdt-temp*sinomdt
c     ------
c B is not aligned with the z-accis (Rotate back)
                  if(cB.lt.0.999) then
                     temp=xp(5,i)
                     xp(5,i)=temp*cB+xp(6,i)*sB
                     xp(6,i)=xp(6,i)*cB-temp*sB 
                  endif
c Account for the E*B drift (Transform back)
                  xp(5,i)=xp(5,i)+vd*sd
                  xp(6,i)=xp(6,i)+vd*cd

               endif

            else
c     Do nothing
            endif

            if(rn.le.r(1)) then

               ninner=ninner+1

c Collision point
               xc=xp(1,i)
               yc=xp(2,i)
               zc=xp(3,i)

               rad=xc**2+yc**2
c cos(theta) of collision
               ctc=zc/sqrt(rad+zc**2)
c sin(psi) of collision
               spsi=yc/sqrt(rad)
c cos(psi) of collision
               cpsi=xc/sqrt(rad)

c     Compute the ion collection the old way, i.e. full count on the
c     collection cell

               vxy=xp(4,i)*cpsi + xp(5,i)*spsi
               vr=vxy*sqrt(1-ctc**2) + xp(6,i)*ctc

c     When the time-step is too large, it can happen that a particle
c     crosses half the sphere and vr becomes positive. In that rare
c     case, put vr=abs(vr), otherwise problem in how to calculate the
c     density on the first radial cells.
               if(vr.gt.0) vr=-vr

               if(.not.collcic) then
c     Interpolate onto the mesh as in ptomesh

                  ithc=interpth(ctc,tfc)
                  if(LCIC)then
                     icell=nint(ithc+tfc)
                  else
                     icell=ithc
                  endif
                  jpsic=interppsi(spsi,cpsi,pfc)
                  jcell=nint(jpsic+pfc)
                  if(jcell.eq.npsiused+1)jcell=1
                  nincell(icell,jcell)=nincell(icell,jcell)+1
                  vrincell(icell,jcell)=vrincell(icell,jcell)+vr
                  vr2incell(icell,jcell)=vr2incell(icell,jcell)+vr**2
c     Compute the ion collection by linear extrapolation
               else
                  ithc=interpth(ctc,tfc)
                  jpsic=interppsi(spsi,cpsi,pfc)
                  nincell(ithc,jpsic)=nincell(ithc,jpsic)
     $                 +(1-tfc)*(1-pfc)
                  nincell(ithc+1,jpsic)=nincell(ithc+1,jpsic)
     $                 +tfc*(1-pfc)
                  vrincell(ithc,jpsic)=vrincell(ithc,jpsic)
     $                 +vr*(1-tfc)*(1-pfc)
                  vrincell(ithc+1,jpsic)=vrincell(ithc+1,jpsic)
     $                 +vr*tfc*(1-pfc)
                  vr2incell(ithc,jpsic)=vr2incell(ithc,jpsic)
     $                 +vr**2*(1-tfc)*(1-pfc)
                  vr2incell(ithc+1,jpsic)=vr2incell(ithc+1,jpsic)
     $                 +vr**2*tfc*(1-pfc)
                  if(jpsic.eq.npsiused)jpsic=0
                  nincell(ithc,jpsic+1)=nincell(ithc,jpsic+1)
     $                 +(1-tfc)*pfc
                  nincell(ithc+1,jpsic+1)=nincell(ithc+1,jpsic+1)
     $                 +tfc*pfc
                  vrincell(ithc,jpsic+1)=vrincell(ithc,jpsic+1)
     $                 +vr*(1-tfc)*pfc
                  vrincell(ithc+1,jpsic+1)=vrincell(ithc+1,jpsic+1)
     $                 +vr*tfc*pfc
                  vr2incell(ithc,jpsic+1)=vr2incell(ithc,jpsic+1)
     $                 +vr**2*(1-tfc)*pfc
                  vr2incell(ithc+1,jpsic+1)=vr2incell(ithc+1,jpsic+1)
     $                 +vr**2*tfc*pfc
               endif

c     Collected momentum and energy

               zmomprobe=zmomprobe+xp(6,i)
               xmomprobe=xmomprobe+xp(4,i)
               ymomprobe=ymomprobe+xp(5,i)
               enerprobe=enerprobe+0.5*v2


            elseif(rn.ge.r(nr))then
c     Left the grid outer boundary.
               zmout=zmout-xp(6,i)
               xmout=xmout-xp(4,i)
               ymout=ymout-xp(5,i)

c     Did not leave the grid. Jump to subcycle end.
            else
               goto 81
            endif


c We left. If we haven't exhausted complement, restart particle i.
            if(nrein.lt.ninjcomp) then


               call reinject(i,dtin,icolntype,bcr)
c              Keep record of injected particles for testing
               if (mcrninjd.lt.npartmax) then
                  mcrninjd = mcrninjd + 1
                  mcrxpinjd(1,mcrninjd) = xp(1,i)
                  mcrxpinjd(2,mcrninjd) = xp(2,i)
                  mcrxpinjd(3,mcrninjd) = xp(3,i)
                  mcrxpinjd(4,mcrninjd) = xp(4,i)
                  mcrxpinjd(5,mcrninjd) = xp(5,i)
                  mcrxpinjd(6,mcrninjd) = xp(6,i)
               endif

               ipf(i)=1
               zmout=zmout+xp(6,i)
               xmout=xmout+xp(4,i)
               ymout=ymout+xp(5,i)
               
               vzinit(i)=xp(6,i)
               if(i.le.norbits) then
                  if (.not.(orbinit))
     $                 iorbitlen(i)=0
               endif
c If an ion is reinjected but was to collide outside, it should not collide
c after reinjection !
               lcstep=.false.


c     New reinjection handling. Simply use the rest of the time step with
c     the new particle starting just at the edge. Get new position:
               ih=1
               hf=77.
               call ptomesh(i,il,rf,ith,tf,ipl,pf,st,ct,
     $              sp,cp,rp,zetap,ih,hf)
c              If lgotooutput tripped, skip to end
               if (lgotooutput) then
                  write (*,*) 'lgotooutput tripped in padvnc.'
                  goto 401
               endif
c     Set the external step length, (isubcycle=1).
               dts=dtin
c     The ion is reinjected with v and x synchronized. We set dtprec to
c     zero to offset v and x by half a timestep
               dtprec(i)=0
c     Call the timestep fraction-remaining random.
               remdt=dtin*ran0(idum)

c     Try to use the real remaining time for remdt (remember dtl is negative)
c     I think it's wrong
c               remdt=remdt+dtl


c     Jump to subcycle end.
               goto 81
            else
               ipf(i)=0
            endif
c Break from subcycles after dealing with a particle that left.
            goto 82  
           
 81         continue
c     Explicit cycle controlled by remaining time in step: (Sometimes
c     problems with roundings when adding multiple subcycles to remdt,
c     so put a cat at 10^-8)
            if(remdt.gt.1e-8) goto 80
c     .................... End of Subcycle Loop .................
c     Break jump point:
 82         continue
c -----------------------------------------------------------
               
            rn=sqrt(xp(1,i)**2+xp(2,i)**2+xp(3,i)**2)

            if(ldist) then
c Start of Various distribution diagnostics.
c     Diagnostics of f_r(rmax):
               if(rn.gt.r(nr-1))then
                  v=(xp(4,i)*xp(1,i)+xp(5,i)*xp(2,i)+xp(6,i)*xp(3,i))/rn
                  ivdiag=1+max(0,nint(nvmax*(v/vrange + .499)))
                  if(ivdiag.gt.nvmax) ivdiag=nvmax
                  nvdiag(ivdiag)=nvdiag(ivdiag)+1
               elseif(rn.gt.r(ircell).and.rn.le.r(ircell+1))then
c     Inner distribution Diagnostics: Assumes reinject never gets here.
                  ctc=xp(3,i)/rn
                  ithc=interpth(ctc,thc)
                  if(ithc.eq.itcell)then
                     vz=xp(6,i)
                     vxy=(xp(4,i)*xp(1,i)+xp(5,i)*xp(2,i))/
     $                    sqrt(xp(1,i)**2+ xp(2,i)**2)
                     vr=vz*ct+vxy*st
                     vt=-vz*st+vxy*ct
c     Radial
                     ivdiag=1+max(0,nint(nvmax*(vr/vrange + .499)))
                     if(ivdiag.gt.nvmax) ivdiag=nvmax
                     vrdiagin(ivdiag)=vrdiagin(ivdiag)+1
c     Angular
                     ivdiag=1+max(0,nint(nvmax*(vt/vrange + .499)))
                     if(ivdiag.gt.nvmax) ivdiag=nvmax
                     vtdiagin(ivdiag)=vtdiagin(ivdiag)+1
c     write(*,502)rn,ithc,vr
                  endif
               endif
            endif
c     Orbit diagnostics
            if(i.le.norbits) then
               iorbitlen(i)=iorbitlen(i)+1
               xorbit(iorbitlen(i),i)=xp(1,i)
               yorbit(iorbitlen(i),i)=xp(2,i)
               rorbit(iorbitlen(i),i)=sqrt(xp(1,i)**2+xp(2,i)**2)
               zorbit(iorbitlen(i),i)=xp(3,i)
               vxorbit(iorbitlen(i),i)=xp(4,i)
               vyorbit(iorbitlen(i),i)=xp(5,i)
               vzorbit(iorbitlen(i),i)=xp(6,i)
c     write(*,503)i,iorbitlen(i),xorbit(iorbitlen(i),i)
c     $           ,yorbit(iorbitlen(i),i),zorbit(iorbitlen(i),i)
c     $           ,rorbit(iorbitlen(i),i)
            endif
c------------------------End distribution diagnostics ---------------
            if(ipf(i).gt.0)iocthis=i
            
         elseif(nrein.lt.ninjcomp)then

c ```````````````````````````````````````` Treatment of INactive slot.
c Case for ipf(i) le 0 (empty slot) but still wanting to inject. 
c We should not come here unless .not.lfixedn.
c            write(*,*)'Reinjecting empty slot',i


c Still need to work on here.

            call reinject(i,dtin,icolntype,bcr)
            dtprec(i)=dtin
            ipf(i)=1
            iocthis=i
         elseif(i.ge.iocprev)then
c     Break if inactive slot and we have exhausted the complement of
c     injections.  And we have reached the maximum occupied slot of
c     previous run.
            goto 401
         endif
c---------------- End of padvnc particle iteration ------------------


c     add the current particle velocity synchronized with its current
c     position to the total current. Only first order in accel. In
c     theory we should need to recalculate accel, but we don't do it
c     until the next timestep

         if(rn.le.rcc(ierad)) then
            curr(1)=curr(1)+xp(4,i)+0.5*dt*accel(1)
            curr(2)=curr(2)+xp(5,i)+0.5*dt*accel(2)
            curr(3)=curr(3)+xp(6,i)+0.5*dt*accel(3)
            curr(4)=curr(4)+1
         endif

      enddo


 401  continue

      NCneutral=ncollide
c      write(*,*)'ncollide=',ncollide,' icycle=',icycle
c We just want the diagnostics with the true particles for now
c      iocthis=min(iocthis,npartmax)

      iocprev=iocthis
c     if(.not.lfixedn)write(*,504)ninjcomp,nrein,i,iocprev
 504  format('  ninjcomp=',i6,'  nrein=',i6,'  i=',i6,
     $     '  iocprev=',i6)
 503  format('Orbit',i3,' length=',i5,' position=',4f7.3)
 502  format('Distrib. rn=',f6.3,' ithc=',i4,' vr=',f6.3)
 501  format('accel=',3f11.4,' xp=',3f11.4)


      end
c***********************************************************************
c***********************************************************************
c Version using precalculated functions. About 30% faster.
      subroutine ptomesh(i,irl,rf,ithl,thf,ipl,pf,st,ct,sp,cp,rp
     $     ,zetap,ih,hf)


c Return the left hand mesh point and the fractional mesh distance of the
c position of particle i, in irl,rf,itl,tf,ipl,pf 
c Return the sines and cosines of theta and psi in st,ct,sp,cp
c If ih.ne.0 on entry, calculate the half-mesh postion, zetap,ih,hf.
      implicit none
      integer i
      integer irl,ithl,ipl
      real rf,thf,pf
      real ct,st,cp,sp,rp
      real zetap,hf
      integer ih
      real hp
c      common /angles/ct,st,cp,sp,rp
c Common data:
      include 'piccom.f'
      include 'errcom.f'
      real rsp
      real x,y,z
      external interpth,interppsi
      integer interpth,interppsi

c Testing
      if(.not. xp(1,i).le.400.)then
         write(*,*)'Ptomesh particle overflow on entry'
         write(*,*)i,(xp(ipl,i),ipl=1,6)
         write(*,*)i,irl,rf,ithl,thf,ipl,pf,st,ct,sp,cp,rp,zetap,ih,hf
c         stop
c        Trigger go to output
         lgotooutput = .true.
         goto 404
      endif

C Find the cell and cell fraction we are at.
      x=xp(1,i)
      y=xp(2,i)
      z=xp(3,i)

      rsp=x**2+y**2
c The square roots here cost perhaps 1/3 of this routine. 
      rp=sqrt(rsp+z**2)
c 

      if(.not. rp.le.r(nr))then
         write(*,*)'Ptomesh particle outside on entry'
         write(*,*)'xp:',(xp(ipl,i),ipl=1,6)
         write(*,*)'i,r(nr),rp,zetap,ih,hf'
         write(*,*)i,r(nr),rp,zetap,ih,hf
         write(*,*) 'x:',x,' y:',y,' z:',z
c         stop
c        Trigger go to output
         lgotooutput = .true.
         goto 404
      endif

c psi sin/cos
      rsp=sqrt(rsp)
      if(rsp .gt. 1.e-9) then
         cp=x/rsp
         sp=y/rsp
      else
         cp=1.
         sp=0.
      endif

      ipl=interppsi(sp,cp,pf)
  
      if(pf.lt.0.) then
         write(*,*)"Negative pf from ippre. i,ipl,pcc(ipl),sp,cp,pf="
     $        ,i,ipl,pcc(ipl),sp,cp,pf
         write(*,*) 'r: ',xp(1,i)**2+xp(2,i)**2+xp(3,i)**2
      elseif(pf.gt.1.) then
         write(*,*)"pf>1 from ippre. i,ipl,pcc(ipl),sp,cp,pf="
     $        ,i,ipl,pcc(ipl),sp,cp,pf
         write(*,*) 'r: ',xp(1,i)**2+xp(2,i)**2+xp(3,i)**2
         write(*,*)
      endif
      
c theta sin/cos
      st=rsp/rp
      ct=z/rp

      if(abs(1+int((ct-th(1))*tfac)).gt.ntpre)then
         write(*,*)'ptomesh overflow. Probably particle NAN'
         write(*,*)'i,irl,rf,ithl,thf',i,irl,rf,ithl,thf
         write(*,*)'ct,th(1),tfac,z,rp',ct,th(1),tfac,z,rp
         write(*,*)'xp',xp(1,i),xp(2,i),xp(3,i),xp(4,i),xp(5,i),xp(6,i)
         write(*,*)'x,y,z',x,y,z
c         stop
c        Trigger go to output
         lgotooutput = .true.
         goto 404
      endif
      ithl=interpth(ct,thf)

      irl=irpre(1+int((rp-r(1))*rfac))
      rf=(rp-r(irl))/(r(irl+1)-r(irl))
      if(rf.lt.0.) then
         write(*,*)"Negative rf from irpre. i,ih,irl,rf,rp="
     $        ,i,ih,irl,rf,rp
         write(*,*) 'r: ',xp(1,i)**2+xp(2,i)**2+xp(3,i)**2
      endif
c "While not"      
 402  if(rf.le.1.)goto 401
      if(irl.eq.nr)then
         write(*,*)'ptomesh rf gt 1 error:',rf,irl
c         stop
c        Trigger go to output
         lgotooutput = .true.
         goto 404
      else
         irl=irl+1
         rf=(rp-r(irl))/(r(irl+1)-r(irl))
      endif
      goto 402
 401  continue
c      return
c New section for halfmesh quantities. Adds about 10% to time.
c Now we have identified the whole mesh position. The half mesh is very
c near it, either irl or irl+1.
      if(ih.ne.0)then
         ih=irl+1
         hp=rp-r(1)
         zetap=sqrt(2.*hp)
         hf=zetap-zetahalf(ih)
         if(hf.lt.0.)ih=ih-1
c     This is the halfmesh fraction 'x'
         hf=(zetap-zetahalf(ih))/(zetahalf(ih+1)-zetahalf(ih))
         
         if(hf.gt.1.or.hf.lt.0.or.zetap.lt.0..or.ih.le.0
     $        )then
c     $        .or. ih.eq.NRFULL)then
            write(*,*)'hf error, ih,irl,rf,zetahalf',ih,irl,rf,
     $           zetahalf(ih),zetahalf(ih+1)
            write(*,*)'zetap,zeta(ih),zeta(ih+1),hf',
     $           zetap,zeta(ih),zeta(ih+1),hf
         endif
      endif
 404  continue
      end

c****************************************************************** 
c     Set the finite volumes coefficients for the outer boundary, as
c     well as the probe potential.
      subroutine shielding_bc(dt,n1,icolntype,colnwt)

      include 'piccom.f'
      include 'errcom.f'
c      include 'colncom.f'
      real dt
      integer n1
      real phislopeconst,phislopefac

cc     Position vector and associated quantities
c      real rpos(3), rs, ct, st, cp, sp

c Set appropriate probe potential.
      call innerbc(1,dt)


c Potential calculation in the shielding region

c  bcp=1 -> Quasineutrality on the 15% outer crone
      if(bcphi.eq.1) then
         n1=nint(NRUSED*.85)-1

         do k=1,npsiused
            do j=1,nthused
               do i=n1+1,nrused
                  phi(i,j,k)=log(rho(i,j,k))
               enddo
               gpc(j,k,4)=phi(n1+1,j,k)
            enddo
         enddo


c bcp=2 -> Phiout=0
         
      elseif(bcphi.eq.2) then
         n1=nrused-1

         do k=0,npsiused+1
            do j=0,nthused+1
               phi(n1+1,j,k)=0.
cc              When allowing different parallel neutral and ion drifts,
cc                there needs to be a corresponding electric field at the boundary
c               rs = rcc(n1+1)
c               ct = tcc(j)
c               st = sqrt(1.-ct**2)
c               cp = cos(pcc(k))
c               sp = sin(pcc(k))
c               rpos(3)=rs*ct
c               rpos(2)=(rs*st)*sp
c               rpos(1)=(rs*st)*cp
c               phi(n1+1,j,k) = phi(n1+1,j,k) + dot(rpos, Eneut, 3)

               gpc(j,k,4)=phi(n1+1,j,k)
            enddo
         enddo

c bcp=4 -> dPhi/drout=-Phiout/r
      elseif(bcphi.eq.4) then
         n1=nrused

c bcp=3 -> dPhi/dzout=0
      elseif(bcphi.eq.3) then
         n1=nrused

c bcp=0 -> Use the spherical symmetry approximation (Hutch paper2)
      elseif (bcphi.eq.0) then
         n1=nrused
         redge=rcc(n1)
         delredge=rcc(n1)-rcc(n1-1)
c     Screening k-number combines electrons and ions.
         if(debyelen.gt. 1.e-10) then
            el2=(1.+1./Ti)/debyelen**2
         else
            el2=2.e20
         endif
         el=sqrt(el2)
         afactor=0.02
         alpha=1./(1.+(afactor*redge/debyelen)**2)
         rxl=el*redge
         expE1=(alog(1.+1./rxl) - 0.56/(1.+4.1*rxl+0.9*rxl**2))
         rindex=alpha*(redge*el+1.)+ (1.-alpha)*2.
c At high collisionality reduce the debye gradient term
         if(icolntype.eq.1 .or. icolntype.eq.2)then
            rindex=(rindex-1.)/(1.+(colnwt*redge)**2/Ti)+1.
         endif
         adeficit=0
c     Boundary slope factor calculations:
         do k=1,npsiused
            do j=1,nthused
c     Current fractional ion deficit due to collection.
c     Coefficient of 1/r^2 in modified shielding equation is
c     a = deficitj * r_edge^2 / \lambda_De^2
               deficitj=1-phi(n1,j,k)/Ti -rho(n1,j,k)
c Reduce the deficit term when collisionality is significant.
c Because it no longer applies. (Perhaps ought to account for vd).
               deficitj=deficitj/(1.+(colnwt*redge)**2/Ti)
               blfac1=(deficitj/debyelen**2) * redge
               adeficit=adeficit+blfac1
c     BC modification is (a/r_edge)[exp(EL*r) E_1(El*r)] given by approx.
               blfac=blfac1*expE1
               blfac=alpha*blfac
               phislopeconst=blfac*redge*delredge/
     $              (redge+delredge*rindex*0.5)
               phislopefac=(redge-delredge*rindex*0.5)/
     $              (redge+delredge*rindex*0.5)
c     Set gpc array
               gpc(j,k,1)=1+2*(phislopefac-1)
               gpc(j,k,2)=0.
               gpc(j,k,3)=0.
               gpc(j,k,4)=-2*phislopeconst
               gpc(j,k,5)=0.
            enddo
         enddo
c     Actual a factor averaged over angles:
         adeficit=adeficit*redge/NTHUSED
         if(adeficit.lt.0.)then
c     write(*,*)'Negative adeficit',adeficit,' set to zero'
            adeficit=0.
         endif

      endif

      end

c***********************************************************************
c Initialization for Collisions
      subroutine colninit(colnwt,icolntype)
      real colnwt
      integer icolntype
      include 'piccom.f'
      include 'errcom.f'
      include 'colncom.f'
      integer i
c      write(*,*)'Initialized collisions',colnwt,icolntype
      if(icolntype.eq.1 .or. icolntype.eq.2
     $     .or. icolntype.eq.5 .or. icolntype.eq.6)then
         if (Bz.ne.0.) then
c        Constant nu collisions. Eneutral must be consistent with reldrift:
c          Taking Eneutral and Eneut to just give parallel relative drift
         Eneutral = colnwt*dot(magdir(1),reldrift(1),3)
c           Start with component of Eneut needed to give right ExB dr.
cc            call cross(magdir(1),ecbdrift(1),Eneut(1))
cc            do i=1,3
cc               Eneut(i) = Eneut(i)/Bz
cc            enddo
c           Add component parallel to B to give right average parallel drift
            do i=1,3
c               Eneut(i) = Eneut(i) +
               Eneut(i) =
     $           dot(magdir(1),reldrift(1),3)*colnwt*magdir(i)
            enddo
         else
c           No magnetic field, so set Eneut to give right average drift
            Eneutral = colnwt*sqrt(dot(reldrift(1),reldrift(1),3))
            do i=1,3
               Eneut(i) = colnwt*reldrift(i)
            enddo
         endif
c Testing
c         Eneutral=0.
         if(myid .eq.0) write(*,*)'colnwt,vd,vneutral=',colnwt,vd
     $        ,vneutral,' Eneutral=',Eneutral
     $        
      elseif(icolntype.eq.0)then
c Need more code here for other types. Not yet implemented.
c Must set Eneutral to zero by default.

         Eneutral=0.
      else
         write(*,*)'Incorrect icolntype',icolntype
         stop
      endif
      end
c*******************************************************************
      subroutine fcalc_infdbl(dt)
      include 'piccom.f'
      include 'errcom.f'
      real dt,rmax
      integer imin,kk1,kk2
      imin=1.
      rmax=rcc(nrused)
c      call innerbc(imin,dt)
      decay=debyelen
c     max(debyelen,0.01)/sqrt(1+1/(Ti+vd**2))

      sB=sqrt(1-cB**2)
      sd=sqrt(1-cd**2)

      Exext=-vd*Bz*(cB*sd-sB*cd)

      do k=1,npsiused
         do j=1,nthused
            do i=1,nrused
c     Debye Huckel initialization.
c               phi(i,j,k)=vprobe*r(1)/r(i)*exp(-(r(i)-r(1))/decay)
c     $              +Exext*cos(pcc(k))*sqrt(1-tcc(j)**2)* (r(1)/r(i))**2
c     $              *((r(i)+decay)/(1+decay))*exp(-(r(i)-r(1))/decay)

c  Enclosed Coulomb initialization
c               phi(i,j,k)=vprobe*r(1)/r(i)*(rmax-r(i))/(rmax-r(1))+Exext
c     $              *cos(pcc(k))*sqrt(1-tcc(j)**2)*(r(1)/r(i))**2
c     $              *(rmax**3-r(i)**3)/(rmax**3-r(1)**3)

c Pseudo vanishing DH potential

               phi(i,j,k)=vprobe/r(i)*(exp(-(r(i)-r(1))/decay)
     $              /(1-exp(-2*(rmax-1)/decay))+exp((r(i)-r(1))
     $              /decay)/(1-exp(2*(rmax-1)/decay)))
     $              +Exext*cos(pcc(k))*sqrt(1-tcc(j)**2)* (r(1)/r(i))**2
     $              *((r(i)+decay)/(1+decay))*exp(-(r(i)-r(1))/decay)


            enddo
            phi(0,j,k)=2.5*phi(1,j,k)-2*phi(2,j,k)+0.5*phi(3,j,k)
         enddo
      enddo

      do k=1,npsiused
         kk1=mod(k+3*npsi/2-1,npsi)+1
         kk2=mod(k+(3*npsi+1)/2-1,npsi)+1
         do i=1,nrused
            phi(i,0,k)=0.5*(phi(i,1+imin,kk1)+phi(i,1+imin,kk2))
            phi(i,nthused+imin,k)= 0.5*(phi(i,nthused-imin,kk1)
     $           +phi(i ,nthused-imin,kk2))
         enddo
      enddo
c Set the theta shadow cells to their proper value to ensure periodicity
      do i=0,nrused
         do j=0,nthused+1
            phi(i,j,npsiused+1)=phi(i,j,1)
            phi(i,j,0)=phi(i,j,npsiused)
         enddo
      enddo
c There seems to be a problem on Loki with the following expression
c      write(*,'($)')" "

      end

c*******************************************************************
      subroutine innerbc(imin,dt)
      include 'piccom.f'
      include 'errcom.f'
      real flogfac
      real fluxofangle(nthsize,npsisize)

c parameters to find the floating potentail in presence of Bz
      integer nPhi
      data nPhi/300/
      real Irep(1:300)
      real ncs
      real LS(nthsize,npsisize)
      real z,iota,dpdr,Tau,eta,beta_e,beta_i
      data ncs/25./
c phispan is for floating potential if Bz.ne.0
      phispan=5
      beta_e=0
      iota=0


      if(linsulate.or.lfloat) then

         if(Bz.ne.0) then
            beta_i=Bz*sqrt(2/(Ti*pi))
            beta_e=beta_i*sqrt(Ti)*sqrt(rmtoz*1837.)
            z=beta_e/(1+beta_e)
            iota=1-0.0946*z-0.305*z**2+0.95*z**3-2.2*z**4+1.15*z**5
            
            do j=1,nthused
               do k=1,npsiused
               
                  dpdr=1/phi(1,j,k)*(-phi(3,j,k)+4*phi(2,j,k)
     $                 -3*phi(1,j,k)) /((rcc(3)-rcc(1)))
                  LS(j,k)=-1/(min(dpdr,-1.01)+1)

c               LS(j,k)=debyelen/sqrt(1+1/(Ti+rmtoz*vd**2))
c               LS(j,k)=debyelen/sqrt(1+1/Ti)
c               LS(j,k)=LS(j,k)+debyelen*log(1+1/debyelen)
               enddo
            enddo
         endif



         flogfac=0.5*alog(2.*pi/(rmtoz*1837.))
         totflux=0.
         do j=1,nthused
            do k=1,npsiused
c     Calculate the flux to each angular cell
               if(lcic)then
                  fluxofangle(j,k)=fincellave(j,k)*(nthused-1.)/
     $              (4.*pi*rhoinf*dt*r(1)**2)
                  if(j.eq.1 .or. j.eq.nthused)
     $                 fluxofangle(j,k)=fluxofangle(j,k)*2.
               else
                  fluxofangle(j,k)=fincellave(j,k)*(nthused)/
     $              (4.*pi*rhoinf*dt*r(1)**2)
               endif
               totflux=totflux+fluxofangle(j,k)


               if(linsulate)then

                  if(Bz.ne.0) then
                     do l=1,nPhi
                        Qth=1.-2*(j-1.)/(nthused-1.)
                        eta=phispan*k/nPhi/beta_e*
     $                    (1+beta_e/4*(1-exp(-4/(LS(j,k)*beta_e))))
                        Tau=eta/(1+eta)
                        A=0.678*Tau+1.543*Tau**2-1.212*Tau**3
                        Irep(l)=exp(-phispan*l/nPhi) *((2*(A+(1-A)*iota)
     $                       -1)+2*(1-(A+(1-A)*iota))*abs(Qth))
                     enddo
                  endif
                  if(fluxofangle(j,k).gt.0.)then
                     if(Bz.ne.0) then
                        call invtfunc(Irep,nPhi,exp(flogfac)
     $                       *fluxofangle(j,k),x)
                        phi(imin,j,k)= (-x*phispan/nPhi+(ncs-1.)
     $                       *phi(imin,j,k))/ncs
                     else
                        phi(imin,j,k)=(alog(fluxofangle(j,k))+flogfac
     $                       +(ncs-1.)*phi(imin,j,k))/ncs
                     endif
                  else
                     phi(imin,j,k)=phi(imin,j,k)
                  endif
               endif
            enddo
         enddo

         if(totflux.gt.0.)then
            if(Bz.eq.0) then
               vprobe=(alog(totflux/nthused)+
     $              flogfac+(ncs-1.)*vprobe)/ncs
c     comparison with old version
c               vprobe=alog(totflux/nthused)+flogfac
            else
c     calculate the total e- current to the sphere -> Irep
               do l=1,nPhi
                  Irep(l)=0
                  do j=1,nthused
                     do k=1,npsiused
                        axis=1
                        if (j.eq.1.or.j.eq.nthused) axis=0.5
                        Qth=1.-2*(j-1.)/(nthused-1.)
                        eta=phispan*l/nPhi/beta_e*
     $                       (1+beta_e/4*(1-exp(-4/(LS(j,k)*beta_e))))
                        Tau=eta/(1+eta)
                        A=0.678*Tau+1.543*Tau**2-1.212*Tau**3
                        Irep(l)=Irep(l)+axis*exp(-phispan*l/nPhi) *((2
     $                       *(A+(1-A)*iota)-1)+2*(1-(A+(1-A)*iota))
     $                       *abs(Qth))
                     enddo
                  enddo
                  Irep(l)=Irep(l)/(nthused-1)/npsiused
               enddo
               call invtfunc(Irep,nPhi,exp(flogfac)
     $              *totflux/nthused,x)
               vprobe=(-x*phispan/nPhi+(ncs-1.)*vprobe)/ncs
            endif
         endif
         
c         write(*,*) vprobe,alog(totflux/nthused)+flogfac
         if(lfloat)then
            do j=1,nthused
               do k=1,npsiused
                  phi(imin,j,k)=vprobe
               enddo
            enddo
         endif
c     write(*,*)
c     write(*,*)'fluxofangle=',(fluxofangle(j),j=1,NTHUSED)
c     write(*,*)'phi=',(phi(imin,j),j=1,NTHUSED)



c If prespecified probe potential
c Account for the convective Electric field due to Bz
      else
         sB=sqrt(1-cB**2)
         sd=sqrt(1-cd**2)
         Exext=-vd*Bz*(cB*sd-sB*cd)
         do j=1,nthused
            do k=1,npsiused
c               phi(imin,j,k)=vprobe
               phi(imin,j,k)=vprobe+Exext*cos(pcc(k))*sqrt(1-tcc(j)**2)
            enddo
         enddo     
      endif

      end



c********************************************************************
      subroutine postcollide(i,tisq)
      include 'piccom.f'
      include 'errcom.f'
      include 'colncom.f'
c Get new velocity; reflects neutral maxwellian shifted by vneut and v_ExB.
      xp(4,i)=tisq*gasdev(idum) + vneut(1) + ecbdrift(1)
      xp(5,i)=tisq*gasdev(idum) + vneut(2) + ecbdrift(2)
      xp(6,i)=tisq*gasdev(idum) + vneut(3) + ecbdrift(3)
      end
