

c***********************************************************************
c Orbit injection for general gyrotropic distribution function.
c***********************************************************************
c Other versions are in other source files.
      subroutine ogenreinject(i,dt)
      integer i
      real dt
c Common data:
      include 'piccom.f'
      parameter (eup=1.e-10)
      external pugen
      logical istrapped
c Testing
      real vdist(nvel)
      real tdist(nthsize)
      real crdist(nthsize),cidist(nthsize)
      common/rtest/crdist,cidist,tdist,vdist

c In this routine we work in velocity units relative to ion thermal till end.
      vscale=sqrt(2.*Ti)
      idum=1
      ilaunch=0
      ilocalfail=0
 1    continue
      ilaunch=ilaunch+1
c Section for dealing with multiple launch failures
      if(mod(ilaunch,1000).eq.0)then
         write(*,*)'ilaunch excessive. averein=',averein,' brcsq=',
     $        brcsq,' ierr=',ierr,' ilocalfail=',ilocalfail
         write(*,*)'Setting averein to zero artificially.'
         averein=0.
         if(ilaunch .gt. 10000)then
c In extremis just dump the particle in the region.
            write(*,*)'Desperation drop of particle in region'
            xp(1,i)=r(nr/2)
            xp(2,i)=0.
            xp(3,i)=0.
            xp(4,i)=0.
            xp(5,i)=0.
            xp(6,i)=0.
            return
         endif
      endif
c Pick normal velocity from cumulative Pu
      y=ran0(idum)
      call finvtfunc(pugen,nvel,y,u)
      iv=u
      dv=u-iv
      u=dv*Vcom(iv+1)+(1.-dv)*Vcom(iv)
      if(dv.gt.1)write(*,*)'Error in u calculation',iv,dv
      vdist(iv)=vdist(iv)+1.
c Pick angle from cumulative Pc based on numerical distributions.
 2    y=ran0(idum)
      y=y*( (1.-dv)*Pc(nQth,iv) + dv*Pc(nQth,iv+1))
      call f2invtfunc(Pc(1,iv),Pc(1,iv+1),nQth,y,xc,(1.-dv),dv)
      ixc=xc
      if(xc.lt.1.) then
         write(*,*) 'Theta choice error xc,y,iv,dv=',
     $        xc,y,iv,dv
         write(*,*)'Pc(1,iv),Pc(1,iv+1),Pc(nQth,iv),Pc(nQth,iv+1)',
     $        Pc(1,iv),Pc(1,iv+1),Pc(nQth,iv),Pc(nQth,iv+1)
c         call autoplot(Qcom,Pc(1,iv),nQth)
c         call polyline(Qcom,Pc(1,iv+1),nQth)
c         call pltend()
         goto 2
      endif
      fxc=xc-ixc
      crt=(1.-fxc)*Qcom(ixc)+ fxc*Qcom(ixc+1)
c Testing angular distribution.
      if(LCIC)then
         icr=(1.+crt)*0.5*(NTHUSED-1) + 1.5
      else
         icr=(1.+crt)*0.5*(nth-1) + 1
      endif
      crdist(icr)=crdist(icr)+1.
c End of testing distribution monitor.
      srt=sqrt(1.- crt**2)
c Pick angle zt of poloidal impact and angle eta of impact parameter
      zt=ran0(idum)*2.*pi
      czt=cos(zt)
      szt=sin(zt)
      eta=ran0(idum)*2.*pi
      ceta=cos(eta)
      seta=sin(eta)
c Choose impact parameter, preventing overflow.
      chium2=-averein/Ti/(u+eup)**2
      if(chium2.le.-1.) then
         write(*,*)'Impossible chium2=',chium2,' averein=', averein,
     $        ' u=',u,' iv=',iv
c         stop
      endif
c      if(.not.lfixedn)chium2=0.
      brcsq=ran0(idum)*(1.+chium2)
c Reject a particle that will not reach boundary.
      if(brcsq.lt.0.) then
         goto 1
      endif
      brc=sqrt(brcsq)
c Get cosine and sine of impact angle relative to distant position.
c Based on integration.
      p2=brcsq*2.*Ti*u**2
      ierr=0
      if(debyelen.gt..001)then
c Orbit integration angle calculation.
c There is an overflow with this at zero debyelen. Ought to be properly fixed.
         call alphaint(p2,brcsq,cosal,sinal,ierr)
         if(ierr.ne.0)goto 1
c      write(*,'(4f9.4)')cosal-alcos(brc,chium2),sinal-alsin(brc,chium2)
c     Now ilaunch is the number of launches at infinity it took to get
c     one that reached the boundary.
      else
c Alternative based on analytic orbit calculation.
c Used for low debyelen, but really assumes negligible boundary potential.
         call alcossin(brc,chium2,cosal,sinal)
         cosal=alcos(brc,chium2)
         sinal=alsin(brc,chium2)
      endif
c Install reinjection position
      a1=crt*ceta*sinal+srt*cosal
      rs=r(nr)*0.99999
      xp(1,i)=rs*(czt*a1+szt*seta*sinal)
      xp(2,i)=rs*(-szt*a1+czt*seta*sinal)
      xp(3,i)=rs*(-srt*ceta*sinal + crt*cosal)

c Obtain angle coordinate and map back to th for phihere.
      ct=xp(3,i)/rs
      call invtfunc(th(1),nth,ct,x)
      ic1=x
      ic2=ic1+1
      dc=x-ic1

      phihere=phi(NRUSED,ic1)*(1.-dc)+phi(NRUSED,ic2)*dc
c Section to correct the injection velocity and direction (but not the
c position) to account for local potential. 26 July 2004.
      if(localinj)then
         brcsq=(brcsq*(1.-phihere/Ti/(u+eup)**2)/(1.+chium2))
         if(brcsq.lt. 0.) then
c     This launch cannot penetrate at this angle. But it would have done
c     if the potential were equal here to averein. Thus it probably
c     should not be counted as part of the launch effort. So
            ilaunch=ilaunch-1
            ilocalfail=ilocalfail+1
            goto 1
         endif
         chium2=-phihere/Ti/(u+eup)**2
         brc=sqrt(brcsq)
      endif
c Injection velocity components normalized in the rotated frame:
      ua1=-brc*cosal -sqrt(1.+chium2-brcsq)*sinal
      ua3=brc*sinal - sqrt(1.+chium2-brcsq)*cosal
      ua=crt*ceta*ua1+srt*ua3
c Install reinjection velocity in Te-scaled units
      u=u*vscale
      xp(4,i)=u*(czt*ua+szt*seta*ua1)
      xp(5,i)=u*(-szt*ua+czt*seta*ua1)
      xp(6,i)=u*(-srt*ceta*ua1 + crt*ua3)
c Remove with new advancing code:
c Increment the position by a random amount of the velocity.
c This is equivalent to the particle having started at an appropriately
c random position prior to reentering the domain.
c      xinc=ran0(idum)*dt
cc      xinc=0.
c      vdx=0.
c      do j=1,3
c         vdx=vdx+xp(j,i)*xp(j+3,i)
c         xp(j,i)=xp(j,i)+xp(j+3,i)*xinc
c      enddo
c      if(vdx.gt.0.)then
c         write(*,*)'Positive projection. u,phi=',u,phihere
c 601     format(a,5G10.5)
c      endif
c      rcyl=xp(1,i)**2+xp(2,i)**2
c      rp=rcyl+xp(3,i)**2
      rp=rs
c Reject particles that are already outside the mesh.
      if(.not.rp.lt.r(nr)*r(nr))then
c      if(.not.rp.le.r(nr)*r(nr))then
c         write(*,*)'Relaunch',rp,xp(1,i),xp(2,i),xp(3,i)
         goto 1
      else

c Do the outer flux accumulation.
c In order to accumulate the number of launches at infinity, rather than
c just the number of reinjections, we weight this by ilaunch
         spotrein=spotrein+phihere*ilaunch
         nrein=nrein+ilaunch
         fluxrein=fluxrein+1.
         if(istrapped(i))then
            ntrapre=ntrapre+ilaunch
c            v=sqrt(xp(4,i)**2+xp(5,i)**2+xp(6,i)**2)
c            write(*,*)'Trapped',vdx/rp,u,v,sqrt(u**2-2.*averein)
c crt,czt,ceta,cosal
         endif
      endif
      end
c********************************************************************
c***********************************************************************
c Calculate the cumulative probability for velocity index iu such that
c         u= vspread*(iu-1.)/(nvel-1.)   as per injinit
      real function pugen(iu)
      integer iu
c     averein is the average potential of reinjected particles, which is
c     used as an estimate of the potential at the reinjection boundary.
c     It is expressed in units of Te so needs to be scaled to Ti.
      include 'piccom.f'
      pudenom=pu1(1)-pu2(1)*averein/Ti
      pugen=1.- (pu1(iu)-pu2(iu)*averein/Ti)/pudenom
      end
c********************************************************************
c Initialize the distributions describing reinjected particles
      subroutine ogeninjinit(icolntype)
      integer icolntype
c Common data:
      include 'piccom.f'
c Passing the drift velocity to fv.
      common /distfunc/ud
c Velocity in this routine is normalized to a nominal ion thermal velocity
c which for a Maxwellian-related form is sqrt(2T_i/m).
      ud=vd/sqrt(2.*Ti)
c Range of velocities permitted for injection.
      vspread=3.+5.*abs(ud)

      do i=1,nQth
c Qcom is here used as 
c the cosine angle of the ith angle interpolation position. 
c         Qcom(i)=1.-2.*(i-1.)/(nQth-1.)
c Perhaps we want this going from -1 to +1 not +1 to -1.
c         Qcom(i)=-1.+2.*(i-1.)/(nQth-1.)
c But it is better for this to be uniform in theta not cos(theta)
         Qcom(i)=cos(3.141593*(1-(i-1.)/(nQth-1.)))
      enddo
c As a function of radial velocity index j
      do j=1,nvel
c on the mesh Vcom
         Vcom(j)=vspread*(j-1.)/(nvel-1.)
         Pc(1,j)=0.
c Integrate fv with respect to costheta.
         do i=2,nQth
            ci=0.5*(Qcom(i)+Qcom(i-1))
            si=sqrt(1-ci**2)
            vx=Vcom(j)*si
            vz=Vcom(j)*ci
            Pc(i,j)=Pc(i-1,j)+(Qcom(i)-Qcom(i-1))*
     $           fvgyro(vx,vz)
c Pc is the integral in cos(angle) Qcom(i) at velocity Vcom(j) of fv
         enddo
      enddo
c Initialize integrations
      pu1(1)=0.
      pu2(1)=0.
      do j=2,nvel
c Integrate along the velocity to get the interpolation functions.
         du=(Vcom(j)-Vcom(j-1))
         p1=(Pc(nQth,j)*Vcom(j)+Pc(nQth,j-1)*Vcom(j-1))*0.5
         p3=(Pc(nQth,j)*Vcom(j)**3+Pc(nQth,j-1)*Vcom(j-1)**3)*0.5
         pu1(j)=pu1(j-1)+du*p3
         pu2(j)=pu2(j-1)+du*p1
      enddo
      do j=1,nvel
c Make pu1,2 monotonically decreasing to zero for orbitinject:
         pu1(j)=pu1(nvel) -pu1(j)
         pu2(j)=pu2(nvel) -pu2(j)
      enddo
c Now pu1(1) = \int u   f(u,c) dc du and
c     pu2(1) = \int u^3 f(u,c) dc du.
c So the flux at infinity is 2\pi*(pu1(1)*u^2-pu2(1)*\chi_b).
c
c For angular comparisons only, not used for actual reinjection here,
c we want the integral over velocity of the flux at angle Qcom.
c We store this in Gcom's first two rows. 
c (The second is to be chi_b weighted).
      do i=1,nQth
         Gcom(1,i)=0.
         Gcom(2,i)=0.
         ci=-Qcom(i)
         si=sqrt(1-ci**2)
         fvp=fvgyro(0.,0.)
         do j=2,nvel
            du=(Vcom(j)-Vcom(j-1))
            fvn=fvgyro(Vcom(j)*si,Vcom(j)*ci)
            Gcom(1,i)=Gcom(1,i)+du*(Vcom(j)**3*fvn+Vcom(j-1)**3*fvp)/2.
            Gcom(2,i)=Gcom(2,i)+du*(Vcom(j)*fvn+Vcom(j-1)*fvp)/2.
            fvp=fvn
         enddo
      enddo
c We store the integral over theta of Gcoms in G(3,1 and 2)
      Gcom(3,1)=0.
      Gcom(3,2)=0.
      do i=2,nQth
         dth=Qcom(i)-Qcom(i-1)
         Gcom(3,1)=Gcom(3,1)+dth*(Gcom(1,i)+Gcom(1,i-1))/2.
         Gcom(3,2)=Gcom(3,2)+dth*(Gcom(2,i)+Gcom(2,i-1))/2.
      enddo
      if(myid.eq.0)then 
         write(*,*)'pu1(1),pu2(1)      ',pu1(1),pu2(1)
         write(*,*)'Gcom(3,1),Gcom(3,2)',Gcom(3,1),Gcom(3,2)
      endif
      call srand(myid+1)
      end
