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

c To use for simple Maxwellian reinjection
c The adiabatic reinjection is considerent as a restriction of the 3D one,
c So that we only have 1 maxinjinit for the 3D case.


      subroutine maxreinject(i,dt)

      integer i
      real dt,sd,temp
c Common data:
      include 'piccom.f'

      vscale=sqrt(Ti)
      vdi=vd/vscale
      cerr=0.
      idum=1
 1    continue
      y=ran0(idum)
c Quick fixing to prevent some errors
      if(y.gt.1-1e-6) then
c         write(*,*) y,nrein
         y=1-1e-6
      endif

c Pick angle from cumulative Q.
      call invtfunc(Qcom,nQth,y,x)
      ic1=x
      if(x.lt.1. .or. x.ge.float(nQth))then
         write(*,*)  'REINJECT Q-Error'
         write(*,*)'y,x,nQth=',y,x,nQth
         write(*,*)'Qcom=',Qcom
         goto 1
      endif
c      if(ic1.ge.nQth)ic1=ic1-1
      ic2=ic1+1
      dc=x-ic1
 2    continue
      x=dc+ic1
      yy=ran0(idum)
c Quick fixing to prevent some errors
      if(yy.gt.1-1e-6) then
c         write(*,*) yy,nrein
         yy=1-1e-6
      endif
c Pick normal velocity from cumulative G.
      call invtfunc(Gcom(1,ic1),nvel,yy,v1)
      call invtfunc(Gcom(1,ic2),nvel,yy,v2)
      vr=dc*v2+(1.-dc)*v1

      if(vr.lt.1. .or. vr.ge.nvel) then
         write(*,*) 'REINJECT V-Error'
         write(*,*) yy,v1,v2,ic1,ic2,nvel,vr,nQth
         goto 2
      endif
      iv=vr
      dv=vr-iv
      vr=dv*Vcom(iv+1)+(1.-dv)*Vcom(iv)
c New angle interpolation.
      ct=1.-2.*(x-1.)/(nQth-1)
c Map back to th for phihere.
      call invtfunc(th(1),nth,ct,x)
      ic1h=x
      ic2h=ic1h+1
      dch=x-ic1h
c Old version used th() directly.
c      ct=th(ic1)*(1.-dc)+th(ic2)*dc
c      write(*,*)'ic1,ic2,dc,ct',ic1,ic2,dc,ct
c ct is cosine of the angle of the velocity -- opposite to the radius.      
      st=sqrt(1.- ct**2)
c Now we have cosine theta=c and normal velocity normalized to v_ti.
c Theta and phi velocities are (shifted) Maxwellians but we are working
c in units of vti.
      vt=gasdev(idum)- st*vdi
      vp=gasdev(idum)
c All velocities now.
      p=2.*pi*ran0(idum)
      cp=cos(p)
      sp=sin(p)
c      write(*,*)ct,st,cp,sp


c use averein to reject particles with too low an energy
c      phihere=averein
c To avoid issues at low Ti, cap phihere
c      phihere=max(phihere,-3*Ti)
      vv2=(vt**2 + vr**2 + vp**2)*vscale**2

c     It is better not to use any acceleration. averein is just used in
c     rhoinfcalc.
      phihere=0

c     Reject particles that have too low an energy. forget about the
c     adiabatic reinjection for now
      if(.not.vv2.gt.-2.*phihere) goto 2


c If velocity is normalized to sqrt(Te/mi), and Ti is Ti/Te really,
c then a distribution with standard deviation sqrt(Ti/mi) is obtained
c from the unit variance random distribution multiplied by sqrt(Ti)=vscale.

c Here we consider the drift to be in the z-direction
      xp(6,i)=(vr*ct - vt*st)*vscale
      xp(5,i)=((vr*st+ vt*ct)*sp + vp*cp)*vscale
      xp(4,i)=((vr*st+ vt*ct)*cp - vp*sp)*vscale

      rs=-r(nr)*0.99999
      xp(3,i)=rs*ct
      xp(2,i)=(rs*st)*sp
      xp(1,i)=(rs*st)*cp

c Rotate the drift/position to the right direction
c      write(*,*) (xp(iw,i),iw=1,6)
      sd=sqrt(1-cd**2)
      temp=xp(2,i)
      xp(2,i)=temp*cd+xp(3,i)*sd
      xp(3,i)=-temp*sd+xp(3,i)*cd
      temp=xp(5,i)
      xp(5,i)=temp*cd+xp(6,i)*sd
      xp(6,i)=-temp*sd+xp(6,i)*cd
      
      rcyl=xp(1,i)**2+xp(2,i)**2
      rp=rcyl+xp(3,i)**2

c Do the outer flux accumulation.
      spotrein=spotrein+phihere
      nrein=nrein+1

c Direct ic1 usage
      end

c********************************************************************
c Initialize the distributions describing reinjected particles
      subroutine maxinjinit()
c Common data:
      include 'piccom.f'
      real chi
      integer*2 idum
      real gam(nQth)
c      character*1 work(nvel,nth)

c Range of velocities (times (Ti/m_i)^(1/2)) permitted for injection.
      vspread=5.+abs(vd)/sqrt(Ti)
c Random interpolates
      sq2pi=1./sqrt(2.*pi)
      sq2=1./sqrt(2.)
      Qcom(1)=0.
      dqp=0.

      do i=1,nQth
         t=NTHUSED*i/nQth
c Depending on the reinjection, we have a flux depending on chi or not
         chi=diagchi(int(t)+1)*(t-int(t))+diagchi(int(t))*(int(t)+1-t)

c Qth is the cosine angle of the ith angle interpolation position. 
c Qsin is the corresponding sinus
c     We used to used th(i). This is the equivalent definition.
           Qth=1.-2.*(i-1.)/(nQth-1.)
           Qsin=sqrt(1-Qth**2)
c Here the drift velocity is scaled to the ion temperature.
           vdr=vd*Qth/sqrt(Ti)
           dqn=sq2pi*exp(-0.5*vdr**2)+.5*vdr*erfcc(-sq2*vdr)
           if(i.gt.1) then
              Qcom(i)=Qcom(i-1) + dqp +dqn
           endif
   
c     Gamma is the total flux over all velocities at this angle. 
c But it is never used 
         gam(i)=dqn
         dqp=dqn
c At this angle,
         do j=1,nvel
c Contruct the Cumulative distribution of radial velocity
c On the mesh Vcom
            Vcom(j)=vspread*(j-1.)/(nvel-1.)
c The cumulative distribution is Gcom, 
            Gcom(j,i)=(dqn - sq2pi*exp(-0.5*(Vcom(j)-vdr)**2)
     $           - .5*vdr*erfcc(sq2*(Vcom(j)-vdr)) )/dqn
         enddo
         Gcom(1,i)=0.
         Gcom(nvel,i)=1.
      enddo
      do i=2,nQth
         Qcom(i)=Qcom(i)/Qcom(nQth)
      enddo
c Now Gcom(*,i) is the cumulative distribution of radial velocity at cos(Qth)
c normalized to the ion thermal velocity, not sqrt(T_e/m_i).
c And Qcom() is the cumulative distribution in cosine angles Qth
c      work(1,1)=' '
     
 501  format(a,11f8.4)
      idum=4
      call srand(myid+1)
      end

