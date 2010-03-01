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
      subroutine reinject(i,dt)
      integer i
      real dt
c Common data:
      include 'piccom.f'
      parameter (eup=1.e-7)
      external pu
      logical istrapped
c Testing
      real vdist(nvel)
      real tdist(nthsize)
      real crdist(nthsize),cidist(nthsize)
      common/rtest/crdist,cidist,tdist,vdist

c In this routine we work in velocity units relative to ion thermal till end.
      vscale=sqrt(2.*Ti)
      idum=1
 1    continue
c Pick normal velocity from cumulative Pu
      y=ran0(idum)
      call finvtfunc(pu,nvel,y,u)
      iv=u
      dv=u-iv
      u=dv*Vcom(iv+1)+(1.-dv)*Vcom(iv)
      if(dv.gt.1)write(*,*)'Error in u calculation',iv,dv
      vdist(iv)=vdist(iv)+1.
c Pick angle from cumulative Pc.
      y=ran0(idum)
c Here the drift velocity is scaled to the ion temperature.
      Uc=vd/vscale
      uu2=2.*Uc*u
      if(uu2.gt.50.) then
         crt=1.+alog(y)/uu2
      elseif(uu2.lt.-50.) then
         crt=-1.+alog(1-y)/uu2
      elseif(abs(uu2).lt.1.e-5)then
         crt=2.*y -1.
      else
         expuu2=exp(uu2)
c This expression is evaluated very inaccurately if expuu2 is nearly 1.
c That is why such a large limit on abs(uu2) is adopted.
         crt=alog(y*expuu2 + (1-y)/expuu2)/uu2
c The following do not do any better at solving this problem.
c         crt=alog( (y*expuu2 + (1-y)/expuu2)**(1./uu2))
c         crt=-1.+alog(1.+(expuu2**2-1.)*y)/uu2
      endif
      if(.not. abs(crt).le.1)then
c         write(*,*)'Strange crt:',crt,y,expuu2,uu2
c It seems impossible to avoid occasional strange results when uu2 is small.
         crt=2.*y-1.
      endif
c Testing
      icr=(1.+crt)*0.5*(nth-1) + 1
      crdist(icr)=crdist(icr)+1.
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
      brcsq=ran0(idum)*(1.+chium2)
c Reject a particle that will not reach boundary.
      if(brcsq.lt.0.)goto 1
      brc=sqrt(brcsq)
c Get cosine and sine of impact angle relative to distant position.
c Based on orbit calculation.
      call alcossin(brc,chium2,cosal,sinal)
C      cosal=alcos(brc,chium2)
c      sinal=alsin(brc,chium2)
c Install reinjection position
      a1=crt*ceta*sinal+srt*cosal
      rs=r(nr)
      xp(1,i)=rs*(czt*a1+szt*seta*sinal)
      xp(2,i)=rs*(-szt*a1+czt*seta*sinal)
      xp(3,i)=rs*(-srt*ceta*sinal + crt*cosal)

c Obtain angle coordinate and map back to th for phihere.
      ct=xp(3,i)/rs
      call invtfunc(th(1),nth,ct,x)
      ic1=x
      ic2=ic1+1
      dc=x-ic1
c This expression should work for CIC And NGP.
      phihere=phi(NRUSED,ic1)(1.-dc)+phi(NRUSED,ic2)*dc
c Injection velocity components normalized in the rotated frame:
      ua1=-brc*cosal -sqrt(1.+chium2-brcsq)*sinal
      ua3=brc*sinal - sqrt(1.+chium2-brcsq)*cosal
      ua=crt*ceta*ua1+srt*ua3
c Install reinjection velocity in Te-scaled units
      u=u*vscale
      xp(4,i)=u*(czt*ua+szt*seta*ua1)
      xp(5,i)=u*(-szt*ua+czt*seta*ua1)
      xp(6,i)=u*(-srt*ceta*ua1 + crt*ua3)
c Increment the position by a random amount of the velocity.
c This is equivalent to the particle having started at an appropriately
c random position prior to reentering the domain.
      xinc=ran0(idum)*dt
c      xinc=0.
      vdx=0.
      do j=1,3
         vdx=vdx+xp(j,i)*xp(j+3,i)
         xp(j,i)=xp(j,i)+xp(j+3,i)*xinc
      enddo
      if(vdx.gt.0.)then
         write(*,*)'Positive projection. u,phi=',u,phihere
 601     format(a,5G10.5)
      endif
      rp=xp(1,i)**2+xp(2,i)**2+xp(3,i)**2
c Reject particles that are already outside the mesh.
      if(.not.rp.le.r(nr)*r(nr))then
c         write(*,*)'Relaunch',rp,xp(1,i),xp(2,i),xp(3,i)
         goto 1
      else
c Do the outer flux accumulation.
         spotrein=spotrein+phihere
         nrein=nrein+1
         if(istrapped(I))then
            ntrapre=ntrapre+1
            v=sqrt(xp(4,i)**2+xp(5,i)**2+xp(6,i)**2)
c            write(*,*)'Trapped',vdx/rp,u,v,sqrt(u**2-2.*averein)
c crt,czt,ceta,cosal
         endif
      endif
      end
c********************************************************************
c Initialize the distributions describing reinjected particles
      subroutine injinit()
c Common data:
      include 'piccom.f'

c Here the drift velocity is scaled to the ion temperature.
c And U's are in units of sqrt(2T/m), unlike vd.
      Uc=abs(vd)/sqrt(2.*Ti)
c Can't use these formulas for Uc exactly equal to zero.
      if(abs(Uc).lt.1.e-4)then
         Uc=1.e-4
      endif

c Range of velocities (times (Ti/m_i)^(1/2)) permitted for injection.
      vspread=5.+abs(Uc)

      do i=1,nvel
         u0= vspread*(i-1.)/(nvel-1.)
         Vcom(i)=u0
         uplus=u0+Uc
         uminus=u0-Uc
         pu2(i)=0.5*sqrt(pi)*(erfcc(uminus)-erfcc(uplus))
         pu1(i)=0.5*(-uminus*exp(-uplus**2)+uplus*exp(-uminus**2))
     $        +(Uc**2 +0.5)*pu2(i)
      enddo

      call srand(myid+1)
      end
c***********************************************************************
c Calculate the cumulative probability for velocity index iu such that
c         u= vspread*(iu-1.)/(nvel-1.)   as per injinit
      real function pu(iu)
      integer iu
c     averein is the average potential of reinjected particles, which is
c     used as an estimate of the potential at the reinjection boundary.
c     It is expressed in units of Te so needs to be scaled to Ti.
      include 'piccom.f'
      pudenom=pu1(1)-pu2(1)*averein/Ti
      pu=1.- (pu1(iu)-pu2(iu)*averein/Ti)/pudenom
      end
c********************************************************************
c Given a monotonic (increasing?) 
c function Q(x) on a 1-D grid x=1..nq, solve Q(x)=y for x.
c That is, invert Q to give x=Q^-1(y).
      subroutine finvtfunc(Q,nq,y,x)
c Somehow this breaks the passing of a function reference.
c      implicit none
c      real external Q
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
      x=(y-Ql)/(Qr-Ql)+iql
      end
c**********************************************************************
cInverse square law (phi\propto 1/r) injection functions:
c**********************************************************************
      subroutine alcossin(s,c,cosal,sinal)
      real s,c,cosal,sinal
      cosal=alcos(s,c)
      sinal=alsin(s,c)
      end
c**********************************************************************
      real function alcos(s,c)
      real s,c
      if(s.le.1.e-12*c)then
         alcos=-1.
         return
      else
         r=c/(2.*s)
         alcos=-(sqrt(1.+c-s**2)-(s-r)*r)/(1+r**2)
      endif
      end
c**********************************************************************
      real function alsin(s,c)
      real s,c
      if(s.le.1.e-12*c)then
         alsin=0.
         return
      else
         r=c/(2.*s)
         alsin=(sqrt(1.+c-s**2)*r+(s-r))/(1+r**2)
      endif
      end


