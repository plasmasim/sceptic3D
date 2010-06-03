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
c Interpolate onto the theta mesh. Return nearest index, fraction in thf.
      integer function interpth(ct,thf)
      include 'piccom.f'
      include 'errcom.f'
      ithl=itpre(1+int((ct-th(1))*tfac))
      thf=(ct-th(ithl))/(th(ithl+1)-th(ithl))
      if(thf.gt.1.)then
         if(ithl+2.le.NTHFULL)then
            ithl=ithl+1
            thf=(ct-th(ithl))/(th(ithl+1)-th(ithl))
         else
            write(*,*)'INTERPTH error. ithl, thf incorrect'
            write(*,*)ithl,thf,ct
         endif
      endif
      interpth=ithl
      end
c***********************************************************************
c Interpolate onto the psi mesh. Return nearest index, fraction in pf.
      integer function interppsi(sp,cp,pf)
      include 'piccom.f'
      include 'errcom.f'
      real psi

      psi=atan2(sp,cp)
      if(psi.lt.0) psi=psi+2*pi

      ipl=ippre(1+int((psi-pcc(1))*pfac))
      pf=(psi-pcc(ipl))/(pcc(ipl+1)-pcc(ipl))

      if(pf.gt.1.)then
         if(ipl+2.le.NPSIFULL)then
            ipl=ipl+1
            pf=(psi-pcc(ipl))/(pcc(ipl+1)-pcc(ipl))
         else
c            write(*,*)'INTERPPSI error. ipl, psi incorrect'
c            write(*,*)ipl,psi

c     There must be a very slight problem with the number of decimals
c     used for pi. Not a problem
            pf=1.
         endif
      endif
      interppsi=ipl
      end


c***********************************************************************
      subroutine meshinitcic(rmax)
      real rmax
c Common data:
      include 'piccom.f'
      include 'errcom.f'

      r(0)=1.-(rmax-1.)/(NRUSED-1)
      rcc(0)=r(0)
      r(NRFULL)=rmax+(rmax-1.)/(NRUSED-1)
      rcc(NRFULL)=r(NRFULL)

      do i=1,NRUSED
c linear r mesh   
         r(i)=1.+(i-1)*(rmax-1.)/(NRUSED-1)
         rcc(i)=r(i)
c distance from the probe surface, called \rho in notes.
         hr(i)=r(i)-r(1)
         zeta(i)=sqrt(2.*hr(i))
      enddo
c Uniform r-mesh extrapolation.
      zeta(0)=-zeta(2)
      zetahalf(0)=-0.5*(zeta(2)+zeta(3))
      zeta(nr+1)=sqrt(2.*(2.*r(nr)-r(nr-1)-r(1)))
      do i=1,NRUSED
c Half-mesh quantities
         if(i.eq.1)then
            rat=1.
         elseif(i.eq.NRUSED)then
            rat=(sqrt(sqrt(2.*(2.*r(NRUSED)-r(NRUSED-1))))
     $           -sqrt(zeta(i)))/
     $        (sqrt(zeta(i))-sqrt(zeta(i-1)))
         else
            rat=(sqrt(zeta(i+1))-sqrt(zeta(i)))/
     $        (sqrt(zeta(i))-sqrt(zeta(i-1)))
         endif
         zetahalf(i)=0.5*(zeta(i)+zeta(i-1))
         cminus(i)=(rat-2./rat +1.)/6.
         cmid(i)=(rat+1./rat -2.)/6.
         cplus(i)=(1./rat -2.*rat +1)/6.
      enddo
      zetahalf(nr+1)=0.5*(zeta(nr)+zeta(nr+1))
      do i=1,nth
c theta array including poles
c Uniform in cos theta
         th(i)=1.-2.*(i-1)/(nth-1)
      enddo
c Additional angle positions are given past the ends for the purposes of the
c boundary conditions. They are a distance beyond the ends equal to the
c last step.
      th(0)=2.*th(1)-th(2)
      th(NTHUSED+1)=2.*th(NTHUSED)-th(NTHUSED-1)
      do i=1,nth
c     Cic version
         tcc(i)=th(i)
         thang(i)=acos(th(i))
      enddo
      tcc(NTHUSED+1)=2.*tcc(NTHUSED)-tcc(NTHUSED-1)
      thang(NTHUSED+1)=2.*thang(NTHUSED)-thang(NTHUSED-1)
      thang(0)=2.*thang(1)-thang(2)

c Poloidal angles, uniform psi-spacing
      do i=1,npsi
         pcc(i)=0+(i-1)*2*pi/npsi
      enddo
      pcc(0)=2*pcc(1)-pcc(2)
      pcc(NPSIUSED+1)=2*pcc(NPSIUSED)-pcc(NPSIUSED-1)


c For debugging, don't do this since crowding output
c      if(NRUSED.le.10 .and. nth.le.10) then
c         write(*,*)'r,rcc,th,tcc,thang'
c         write(*,*)(r(j),j=0,nrfull)
c         write(*,*)(rcc(j),j=0,nrfull)
c         write(*,*)(th(j),j=0,nthfull)
c         write(*,*)(tcc(j),j=0,nthfull)
c         write(*,*)(thang(j),j=0,nthfull)
c         write(*,*)(pcc(j),j=0,npsifull)
c      endif
cc      write(*,*)'th=',th



c Calculate the mesh volumes
      rim=0.
      rm2=0.
      rm3=0.
      voltot=0.
c Silence warnings. Not otherwise necessary.
      ri=r(1)
      ri1=r(2)

c
      do i=1,NRUSED
         if(i.lt.NRUSED)then
            ri=r(i)
            ri1=r(i+1)
            rs3=ri1**3+ri1**2*ri+ri1*ri**2+ri**3
            rs2=ri1**2+ri1*ri+ri**2
         else
            rs3=0.
            rs2=0.
         endif
         vol= ri1*rs2-0.75*rs3 + 0.75*rm3-rim*rm2
         volinv(i)=3./(4.*pi*vol)
         rim=ri
         rm2=rs2
         rm3=rs3
         voltot=voltot+vol
      enddo
c      if(myid.eq.0) write(*,*)'Voltot',voltot,'   Ratio to actual',
c     $     voltot/(r(NRFULL)**3-r(1)**3)
c      write(*,*)'Volinv',volinv
c Zero the ninth storage.
c in the 3D version, remove the ninthstep storage
c      do k=1,nstepmax
c         do j=1,nth
c            ninthstep(j,k)=0
c         enddo
c      enddo
      call precalc()
      end
c***********************************************************************
c Initializing particles.
      subroutine pinit(colnwt)
c Input variables
c     Collision frequency
      real colnwt
c Common data:
      include 'piccom.f'
      include 'errcom.f'
      real sd

c drift velocity angle
      sd=sqrt(1-cd**2)

c For now use the whole array.
      ntries=0
      ntrapped=0
      rmax=r(NRUSED)
      rmax2=rmax*rmax
      idum=1
      if(rmax2.le.1.) stop 'Error: rmax is less than 1.'

c     We initialize the 'true' particles'
c     Set velocities if using bcr=3
      if (bcr.eq.3) then
         call mcgenpart(xp,6,npart,4,colnwt)
      endif
      do i=1,npart
c        Set positions
         ipf(i)=1
 1       continue
         ntries=ntries+1
         xp(1,i)=rmax*(2.*ran0(idum)-1.)
         xp(2,i)=rmax*(2.*ran0(idum)-1.)
         xp(3,i)=rmax*(2.*ran0(idum)-1.)
         rc=0.
         do j=1,3
            rc=rc+xp(j,i)**2
         enddo
c        If we are not in the plasma region, try again.
         if(rc.ge.rmax2 .or. rc.le.1.) goto 1

c        Set velocities if using bcr!=3
         if (bcr.ne.3) then
            Ti0=Ti
            tisq=sqrt(Ti0)
            xp(4,i)=tisq*gasdev(idum)
            xp(5,i)=tisq*gasdev(idum) + vd*sd
            xp(6,i)=tisq*gasdev(idum) + vd*cd
         endif
   

c         if(istrapped(i).eq.1)then
c            ntrapped=ntrapped+1
c If this goto is included then trapped particles are rejected.
c But that tends to deplete the region close to the probe.
c            goto 1
c         endif
c         if(bcr.eq.2) then
c     Remove particles with too low vz
c            if (istrapped2(i)) then
c               ntrapped=ntrapped+1
c               goto 1
c            endif
c         endif

c     vzinit is the z momentum a particle had when reinjected. allows to
c     get the usual Fc (collection force)
         vzinit(i)=xp(6,i)

      enddo


c Set flag of unused slots to 0
c      do i=npart+1,npartmax
c         ipf(i)=0
c      enddo
      
c      write(*,*)'Initialized ','id=',myid,
c     $     '  n=',npart,'  ntries=',ntries,'  ntrapped=',ntrapped
c Initialize rhoinf:
      rhoinf=numprocs*npart/(4.*pi*r(NRUSED)**3/3.)
c Initialize orbit tracking
      do ko=1,norbits
         iorbitlen(ko)=0
      enddo


      end
c***********************************************************************
c     Initializing the fields Remove the section to read an external
c     potential, since reading the external particle position is enough
c     (potential straightforwardly obtained by poisson's equation
      subroutine finit()
      include 'piccom.f'
      include 'errcom.f'
      real decay
      
c     When the simulation starts, the ion density is uniform, so start
c     with a Debye Huckel form only accounting for the electron response
      decay=debyelen
c max(debyelen,0.01)/sqrt(1+1/(Ti+vd**2))

      sB=sqrt(1-cB**2)
      sd=sqrt(1-cd**2)
      Exext=-vd*Bz*(cB*sd-sB*cd)
        
      do k=1,npsiused
         do j=1,nthused
            do i=1,nrused
c     Debye Huckel initialization.
               phi(i,j,k)=vprobe*r(1)/r(i)*exp(-(r(i)-r(1))/decay)
     $              +Exext*cos(pcc(k))*sqrt(1-tcc(j)**2)* (r(1)/r(i))**2
     $              *((r(i)+decay)/(1+decay))*exp(-(r(i)-r(1))/decay)
            enddo
            phi(0,j,k)=2.5*phi(1,j,k)-2*phi(2,j,k)+0.5*phi(3,j,k)
         enddo
      enddo

c     We must average the potential of each psi-cell at theta=0 or
c     theta=pi. The boundary conditions assume the cell center is on
c     axis, while what is indeed calculated is the potential at the
c     center of the first/last theta-cells.

      do i=0,nrused
         psiave1=0.
         psiave2=0.
         do k=1,npsiused
            psiave1=psiave1+phi(i,1,k)
            psiave2=psiave2+phi(i,nthused,k)
         enddo
         psiave1=psiave1/npsiused
         psiave2=psiave2/npsiused
         do k=1,npsiused
            phiaxis(i,1,k)=phi(i,1,k)
            phiaxis(i,2,k)=phi(i,nthused,k)
            phi(i,1,k)=psiave1
            phi(i,nthused,k)=psiave2
         enddo
      enddo

c We must set the potential of the shadow theta-cells to the
c physical value of the potential at psi+pi. The zero-derivative
c condition is not valid anymore in the 3D case
      do k=1,npsiused
         kk1=mod(k+3*npsi/2-1,npsi)+1
         kk2=mod(k+(3*npsi+1)/2-1,npsi)+1
         do i=0,nrused
            phi(i,0,k)=0.5*(phi(i,2,kk1)+phi(i,2,kk2))
            phi(i,nthused+1,k)=
     $           0.5*(phi(i,nthused-1,kk1)+phi(i,nthused-1,kk2))
         enddo
      enddo
         
c Set the psi shadow cells to their proper value to ensure periodicity
      do j=0,nthused+1
         do i=0,nrused
            phi(i,j,npsiused+1)=phi(i,j,1)
            phi(i,j,0)=phi(i,j,npsiused)
         enddo
      enddo

c      call fcalc_infdbl(0.1)

      end
c************************************************************************
      subroutine precalc()
c     Precalculation functions
      include 'piccom.f'
      include 'errcom.f'

      
      rfac=(nrpre-1.)/(r(NRUSED)-r(1))
      tfac=(ntpre-1.)/(th(nth)-th(1))
c     Since the psimesh is periodic, unlike the theta mesh, we need to
c     go until npsifull=npsiused+1. psi(npsifull)=psi(1)+2*pi
      pfac=(nppre-1.)/(pcc(npsifull)-pcc(1))


      do j=1,ntpre
c     finding the theta precalculated mesh.
         thp=(j-1.)/tfac+th(1)
         thl=th(1)
         thr=th(nth)
         itl=1
         itr=nth
 200     if(itr-itl.le.1)goto 210
         itx=(itr+itl)/2
         thx=th(itx)
         if(thx.ge.thp) then
            thl=thx
            itl=itx
         else
            thr=thx
            itr=itx
         endif
         goto 200
 210     continue
         itpre(j)=itl
      enddo

c     finding the psi precalculated mesh. 

      do j=1,nppre
         psip=(j-1.)/pfac+pcc(1)
         psil=pcc(1)
         psir=pcc(npsifull)
         ipl=1
         ipr=npsifull
c     Opposite direction than for theta
 202     if(ipr-ipl.le.1) goto 211
         ipx=(ipr+ipl)/2
         psix=pcc(ipx)
         if(psix.le.psip) then
            psil=psix
            ipl=ipx
         else
            psir=psix
            ipr=ipx
         endif
         goto 202
 211     continue
         ippre(j)=ipl
      enddo

c     r grid may be nonlinear find grid position by bisection.
      do j=1,nrpre
         rp=(j-1.)/rfac+r(1)
         rl=r(1)
         rr=r(NRUSED)
         il=1
         ir=NRUSED
 201     ix=(ir+il)/2
         rx=r(ix)
         if(rx.le.rp) then
            rl=rx
            il=ix
         else
            rr=rx
            ir=ix
         endif
         if(ir-il.gt.1)goto 201
c     Now il and ir, rl and rr bracket the radius.
         irpre(j)=il
      enddo
c      write(*,*)'Precalculated the r and theta mesh lookups.'
c      write(*,*)itpre,tfac

c Now irpre(1+int(rp-r(1))*rfac)) is the irl except for rounding etc.
c The irpre spacing must be small enough that the maximum increment of
c irpre from j to j+1 is 1. Then it is possible that the downward rounding
c causes irpre to be at most 1 too small.
c The same applies to itpre

      end


c***********************************************************************
c Initialize the poison iteration coefficients. Must be done after
c mesh initiation.
      subroutine poisinitcic()

c Common data:
      include 'piccom.f'
      include 'errcom.f'
      real rave,dpsi
c Initiate rave to avoid warnings at compilation
      rave=1.

c the psi cells are evenly spaced
      dpsi=pcc(2)-pcc(1)

c Start by calculating the stiffness matrix coefficients for the inner nodes
      do i=2,NRUSED
         ri=  (r(i)+r(i-1))/2.
         rip1=(r(i+1)+r(i))/2.
c
         rave=(rip1**2+2.*rip1*ri+ri**2)/4.
         apc(i)=debyelen**2 *rip1**2/rave/(rcc(i+1)-rcc(i))
     $        /(rip1-ri)
         bpc(i)=debyelen**2 *ri**2/rave/(rcc(i)-rcc(i-1))
     $        /(rip1-ri)

c        For debugging, check for negative elements
         if (apc(i).lt.0 .or. bpc(i).lt.0) then
            write (*,*)'apc or bpc negative: ',apc(i),bpc(i),i
         endif


         do j=1,NTHUSED
            cpc(i,j)=debyelen**2/rave
     $           *(0.5*(nth-1)*sin(acos((th(j+1)+th(j))/2.)))**2
            dpc(i,j)=debyelen**2/rave
     $           *(0.5*(nth-1)*sin(acos((th(j)+th(j-1))/2.)))**2

            if(j.eq.1)then
               cpc(i,j)=2.*cpc(i,j)
               dpc(i,j)=0.
            elseif(j.eq.NTHUSED)then
               dpc(i,j)=2.*dpc(i,j)
               cpc(i,j)=0.
            endif
            epc(i,j)=debyelen**2/rave/sin(acos(th(j)))**2/dpsi**2

            if(j.eq.1)then
               epc(i,j)=debyelen**2/rave
     $              /sin(acos(0.25*(3*th(j)+th(j+1))))**2/dpsi**2
               epc(i,j)=epc(i,j)
            elseif(j.eq.NTHUSED)then
               epc(i,j)=debyelen**2/rave
     $              /sin(acos(0.25*(3*th(j)+th(j-1))))**2/dpsi**2
               epc(i,j)=epc(i,j)
            endif

            fpc(i,j)=apc(i)+bpc(i)+cpc(i,j)+dpc(i,j)+2*epc(i,j)

c           Set factor to multiply each equation by to ensure symmetry of A
            if (lmultpc) then
               multpc(i,j)=rave
               if(j.eq.1)then
                  multpc(i,j)=rave/2.
               elseif(j.eq.NTHUSED)then
                  multpc(i,j)=rave/2.
               endif
            else
               multpc(i,j)=1.
c              Note that i=1 remains zero, so if that element is ever
c                used this needs to be modified
c              For debuggin, also set i=0 elements
               multpc(1,j)=1.
            endif

         enddo
      enddo


c Now, calculate the coefficients for the outer boundary, depending
c on the chosen boundary conditions      

c bcp=1 -> Quasineutrality on the 15% outer crone
      if(bcphi.eq.1) then
         do j=1,nthused
            do k=1,npsiused
c gpc(j,k,4) needs to be recomputed at each iteration
               gpc(j,k,1)=0.
               gpc(j,k,2)=0.
               gpc(j,k,3)=0.
               gpc(j,k,4)=0.
               gpc(j,k,5)=0.
            enddo
         enddo

c bcp=2 -> Phiout=0
      elseif(bcphi.eq.2) then
         do j=1,nthused
            do k=1,npsiused
c Zero for Dirichlet condition
               gpc(j,k,1)=0.
               gpc(j,k,2)=0.
               gpc(j,k,3)=0.
               gpc(j,k,4)=0.
               gpc(j,k,5)=0.
            enddo
         enddo

c bcp=4 -> dPhi/drout=-Phiout/r
      elseif(bcphi.eq.4) then
         delredge=rcc(nrused+1)-rcc(nrused-1)
         do j=1,nthused
            do k=1,npsiused
               gpc(j,k,1)=1
               gpc(j,k,2)=0
               gpc(j,k,3)=0
               gpc(j,k,4)=0
               gpc(j,k,5)=-delredge/rcc(nrused)
            enddo
         enddo

c bcp=3 -> dPhi/dzout=0
      elseif(bcphi.eq.3) then
         delredge=rcc(nrused+1)-rcc(nrused-1)
         delcosth=2./(nthused-1.)
         do j=1,nthused
            do k=1,npsiused
               gpc(j,k,1)=1
               gpc(j,k,2)=-delredge/(2*delcosth*rcc(nrused))*tcc(j)
     $              *(1-tcc(j)**2)
               gpc(j,k,3)=-gpc(j,k,2)
               gpc(j,k,4)=0
               gpc(j,k,5)=-delredge*((1-tcc(j)**2)/(2.*rcc(nrused))
     $              +(1-tcc(j)**2)**(1.5)/debyelen)
            enddo
         enddo

c bcp=0 -> Use the spherical symmetry approximation (Hutch paper2)
      elseif (bcphi.eq.0) then
c In this case the gpc coefficients can not be precalculated. We
c only initialize them to zero.
         do j=1,nthused
            do k=1,npsiused
               gpc(j,k,1)=0.
               gpc(j,k,2)=0.
               gpc(j,k,3)=0.
               gpc(j,k,4)=0.
               gpc(j,k,5)=0.
            enddo
         enddo

      endif

      end
c***********************************************************************
      integer function istrapped(i)

c     Return as logical whether the particle i is trapped or not.  It is
c     considered trapped if the energy available for radial velocity at
c     the outer boundary and at the probe, conserving angular momentum
c     and energy is negative. The assumption of angular momentum
c     conservation is false for non-symmetric situations.

      include 'piccom.f'
      include 'errcom.f'

      ih=0
      hf=66.
      call ptomesh(i,il,rf,ith,tf,ipl,pf,st,ct,sp,cp,rp
     $     ,zetap,ih,hf)

      rn=sqrt(xp(1,i)**2+xp(2,i)**2+xp(3,i)**2)

c The interpolation here might not be correct for both schemes.
      phin=(phi(il,ith,ipl)*(1.-tf)+phi(il,ith+1,ipl)*tf)*(1.-rf) +
     $     (phi(il+1,ith,ipl)*(1.-tf)+phi(il+1,ith+1,ipl)*tf)*rf

c Definition as being that the particle does not leave the domain
c      phie= phi(NRUSED,ith)
c Definition that the particle does not reach infinity.
      phie=0.
      phip= vprobe

      vr2=(xp(4,i)*xp(1,i)+xp(5,i)*xp(2,i)+xp(6,i)*xp(3,i))**2
     $     /rn
      v2=xp(4,i)**2+xp(5,i)**2+xp(6,i)**2
      vt2=v2-vr2

c Domain definition.
c      vte2=(vt2*(rn/rcc(NRUSED))**2)
      vte2=0.
c conservation of angular momentum imposes the particule to have the
c following vphi as minimum on the probe
      vtp2=vt2*rn**2

      if( (vte2 .gt. v2 + 2.*(phin-phie) ) .and.
     $    (vtp2 .gt. v2 + 2.*(phin-phip) ) ) then
         istrapped=1
      else
         istrapped=0
      endif
c      write(*,*)'phin=',phin,'  phie=',phie,'  phip=',phip
c      write(*,*)'vte=',vte,' vtp=',vtp

      end

c***********************************************************************
      logical function istrapped2(i)

c     Return as logical whether the particle i is trapped or not.  It is
c     considered trapped if the energy available for radial velocity at
c     the outer boundary and at the probe, conserving angular momentum
c     and energy is negative. The assumption of angular momentum
c     conservation is false for non-symmetric situations.

      include 'piccom.f'
      include 'errcom.f'

      ih=0
      hf=66.
      call ptomesh(i,il,rf,ith,tf,ipl,pf,st,ct,sp,cp,rp
     $     ,zetap,ih,hf)

      phihere=(phi(il,ith,ipl)*(1.-tf)+phi(il,ith+1,ipl)*tf)*(1.-rf) +
     $     (phi(il+1,ith,ipl)*(1.-tf)+phi(il+1,ith+1,ipl)*tf)*rf

      vv=xp(4,i)**2+xp(5,i)**2+xp(6,i)**2
      vz2=xp(6,i)**2
      if (vz2.le.-2*phihere) then
         istrapped2=.true.
      else
         istrapped2=.false.
      endif
c      write(*,*)'phin=',phin,'  phie=',phie,'  phip=',phip
c      write(*,*)'vte=',vte,' vtp=',vtp

      end
c********************************************************************




      function mtrapped()
      include 'piccom.f'
      include 'errcom.f'
      integer istrapped

      mtrapped=0
      do i=1,iocprev
         if(ipf(i).gt.0) then
            if(istrapped(i).eq.1) mtrapped=mtrapped+1
         endif
      enddo

      end
