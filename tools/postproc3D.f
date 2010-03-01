C Version 3.0 SCEPTIC3D Jan 2008

c*******************************************************************


      character*100 string,filename
      include 'piccompost.f'
c      include 'cic/piccompost.f'

      real rholocal(0:NRFULL,0:NTHFULL,0:NPSIFULL)
c      real thanglocal(0:NTHFULL)
      common /forces/ charge1,ffield1,felec1,fion1,ftot1,
     $     charge2,ffield2,felec2,fion2,ftot2

      real phipic(1000),rhopic(1000)
      real rpic(1000),rpicleft(1000),phicos(1000)
      integer nti0
      parameter (nti0=100)
      integer jstepth
      logical lpcic,lphip,lreaddiag,lgraph,larrows,lconline
      logical ledge,ldens
      real canghere,af

      data lpcic/.false./
      data lconline/.false./
      data lphip/.false./
      data lreaddiag/.false./
      data lgraph/.true./
      data larrows/.false./
      data ledge/.false./
      data ldens/.false./
      data jstepth/1/


c Deal with arguments
      do 1 i=1,iargc()
         call getarg(i,string)
C         write(*,*)'Argument:',string(1:40)
         if(string(1:1) .eq. ' ') then
            goto 3
         endif
         if(string(1:1) .eq. '-') then
         if(string(1:2) .eq. '-p') lgraph=.false. 
         if(string(1:2) .eq. '-f') lphip=.true.
         if(string(1:2) .eq. '-n') ldens=.true.
         if(string(1:2) .eq. '-?') goto 51
         else
            filename=string
         endif
 1    continue
 3    continue
      if(i.eq.1)goto 51

      do i=1,len(filename)-2
         if(filename(i:i+1).eq.'Sp')then
            filename(i:i+1)='Ti'
            write(*,*) 'Using Ti file ',filename
            goto 110
         endif
      enddo
 110  continue


 

c Read the outputfile
      call readoutput(lreaddiag,lpcic,ledge,
     $     filename,rholocal,nrhere,nthhere,npsihere,
     $     phipic,rhopic,rpic,rpicleft,phicos,
     $     rhomax,rhomin,
     $     nrti,phiinf,nastep,nsteps,
     $     dt,rmax,fave,debyelen,vprobe,
     $     icolntype,colnwt,Eneutral,vneutral,Tneutral,
     $     ierr)
      if(ierr.eq.101) goto 101

      goto 102
 101  write(*,*) 'Does not seem to be a Ti... file here'
 102  continue


      do k=1,npsihere
         pcc(k)=2*3.141592*(k-1)/npsihere
      enddo

      

c Set arrow scale
      v1=max(1.,vd)
c Start of Plotting:

      call pfset(3)

      if(lgraph)then
c Time trace of the particle flux


         call yautoplot(fluxprobe,nsteps)
c         call yautoplot(pcc,npsihere)
         call axlabels('step','Particles to probe')

c         pause

         call pltend()

c         pause

         call multiframe(2,1,3)
         call autoplot(rpic,rhopic,nrhere)
         call axlabels('r','angle averaged density')
         call winset(.true.)
         call vecw(rpic(1),1.,0)
         call vecw(rpic(nrhere),1.,1)
         call winset(.false.)
         call autoplot(rpic,phipic,nrhere)
         call axlabels('r','angle averaged potential')
         call winset(.true.)
         call vecw(rpic(1),1.,0)
         call vecw(rpic(nrhere),1.,1)
         call winset(.false.)
         call pltend()
         call multiframe(0,0,0)

      endif


      if(lphip) then
         call condisphi2D(ir,jstepth,nrhere,nthhere,npsihere,v1 ,larrows
     $        ,lconline,lpcic,ledge)

         call pltend()

      endif


      if(ldens) then
         call conrho2D(ir,jstepth,nrhere,nthhere,npsihere,v1 ,larrows
     $        ,lconline,lpcic,ledge)

         call pltend()

         call conrhoPsi(ir,jstepth,nrhere,nthhere,npsihere,v1 ,larrows
     $        ,lconline,lpcic,ledge)

         call pltend()

      endif
      
c     calculate the flux as a function of the azimuthal angle along the
c     drift velocity direction
      do j=1,angsize
         cang(j)=1-2.*(j-1)/(angsize-1)
         fluxang(j)=0.
         aweight(j)=0.
      enddo

      call precalc()
      do j=1,nthhere
         do k=1,npsihere
            canghere=sqrt(1-tcc(j)**2)*sin(pcc(k))*sqrt(1-cd**2)+tcc(j)
     $           *cd
            ial=interpang(canghere,af)
            fluxang(ial)=fluxang(ial)+nincell(j,k)*(1-af)
            fluxang(ial+1)=fluxang(ial+1)+nincell(j,k)*af
            aweight(ial)=aweight(ial)+(1-af)
            aweight(ial+1)=aweight(ial+1)+af
c     Double the flux in the half cells
            if(j.eq.1.or.j.eq.nthhere)then
               aweight(ial)=aweight(ial)-(1-af)/2.
               aweight(ial+1)=aweight(ial+1)-af/2.
            endif
         enddo
      enddo

      do j=1,angsize
         fluxang(j)=fluxang(j)*npsihere*(nthhere-1)/(aweight(j)+1e-5)/(4
     $        .*pi*rhoinf*dt*nastep)
      enddo

      call multiframe(0,0,0)
      call autoplot(cang(1),fluxang,angsize)
      call pltend()
      

      call exit
 51   write(*,*)"Usage: postproc3D [-f -x ... ] filename"
      write(*,*)" -p No profiles",
     $     " -f Potential color plot psi-average",
     $     " -n Density color plot psi-average and x-y (z=0)"
      call exit

      end


c***************************************************************************
c Contouring of the potential averaged over Psi, on distorted mesh.
      subroutine condisphi2D(ir,it,nrhere,nthhere,npsihere,v1,larrows
     $     ,lconline,lpcic,ledge)
      integer ir,it
      real rhomax,rhomin,v1
      logical larrows,lconline,lpcic,ledge
c Common data:
      include 'piccompost.f'
c      include 'cic/piccompost.f'
c      save
c      character*20 cstring
      character cworka(nr*(NTHFULL+1+1))
      integer ncont
      parameter (ncont=12)
      real zclv(2*ncont)
      real xrho(NRFULL+1,0:NTHFULL+1),zrho(NRFULL+1,0:NTHFULL+1)
      real x2rho(NRFULL+1,0:NTHFULL+1),z2rho(NRFULL+1,0:NTHFULL+1)
      real potave(NRFULL+1,0:NTHFULL+1)
c      real Tr(NRFULL+1,0:NTHFULL+1),Ttp(NRFULL+1,0:NTHFULL+1)
      save xrho,zrho
      real basesize
      parameter (basesize=.02)

      if(nthhere.gt.NTHFULL)then
         write(*,*)' Condisplay error. Mesh required:',nrhere,nthhere,
     $        ' Exceeds allocated:',NRFULL,NTHFULL
         stop
      endif

      do i=1,nrhere
         do j=0,nthhere+1
            potave(i,j)=0.
            do k=1,npsihere
               potave(i,j)=potave(i,j)+phi(i,j,k)
            enddo
            potave(i,j)=potave(i,j)/npsihere
         enddo
      enddo      

      call minmax2(potave(1,1),nr+1,nrhere,nthhere,rhomax,rhomin)

      do i=1,nrhere
          do j=1,nthhere
            zrho(i,j)=rcc(i)*tcc(j)
            xrho(i,j)=rcc(i)*sqrt(1.-tcc(j)**2)
            z2rho(i,j)=zrho(i,j)
            x2rho(i,j)=xrho(i,j)            
         enddo
         zrho(i,0)=rcc(i)
         xrho(i,0)=0.
         zrho(i,nthhere+1)=-rcc(i)
         xrho(i,nthhere+1)=0.
         if(lpcic)then
            z2rho(i,1)=zrho(i,0)
            x2rho(i,1)=xrho(i,0)
            z2rho(i,nthhere)=zrho(i,nthhere+1)
            x2rho(i,nthhere)=xrho(i,nthhere+1)
         endif
      enddo
c Adjust the outside angle centers if necessary.
c Should never be necessary.
c      if(tcc(1).eq.1)then
c         tcc(1)=0.25*(3.+tcc(2))
c         tcc(0)=1.
c         tcc(nthhere)=0.25*(-3.+tcc(nthhere-1))
c         tcc(nthhere+1)=-1.
c      endif

c ledge fixing
      if (ledge) then
         rhomax=rhomax/2.
      endif

      zscale=log(1000.)/(ncont-1.)
      do j=1,ncont
c Logarithmic
         zclv(j)=rhomax + (rhomin-rhomax)*
     $        exp(zscale*float(ncont-j))/exp(zscale*float(ncont-1))
c Linear
         zclv(j)=(rhomax-rhomin)*(1.*(j-1)/float(ncont-1))+ rhomin
c         write(*,*) "zcle : ",zclv(j)
      enddo
c      write(*,*) "min,max",rhomin,rhomax

c      write(*,*)'Contours=',zclv
      icl=-ncont

c      rpmax=2.*rcc(nrhere)-rcc(nrhere-1)
      rpmax=rcc(nrhere)
c      write(*,*)rcc(nrhere),nrhere,nthhere,rpmax
c      call multiframe(2,2,3)

      call ticnumset(8)
      call pltinaspect(-rpmax,rpmax,0.,rpmax)
      call accisgradinit(-25000,00000,35000,130000,65000,130000)
      ntype=2+16+32
      call contourl(potave(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv,icl,zrho,xrho,ntype)
      call gradlegend(zclv(1),zclv(abs(icl)),
     $     .1,1.25,.9,1.25,0.02,.true.)

c Call a second time for contours, without the highest.
         if(lconline)then
c If this is not very bipolar
            if(abs(rhomax)-abs(rhomin).gt.0.2*abs(rhomax-rhomin))then
c Rational logarithmic contours
               do indx=1,ncont-2,3
                  base=10.**(-2+(indx-1)/3)
                  zclv(ncont+indx)=-base
                  zclv(ncont+indx+1)=-base*2.
                  zclv(ncont+indx+2)=-base*5.
c     New positive contours.
                  zclv(ncont-indx)=base
                  zclv(ncont-(indx+1))=base*2.
                  zclv(ncont-(indx+2))=base*5.
               enddo
               zclv(ncont)=0.
            endif
            write(*,'(a,30f6.2)')'Contours=',zclv
            ntype=2
            icl=(ncont-1)
c            write(*,'(2f8.4)')
c     $           (zrho(nrhere,k),z2rho(nrhere,k),k=0,nthhere+1)
c            call contourl(phi(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
c     $           zclv(ncont),icl,zrho,xrho,ntype)
c Fixed contour plot avoiding the hacked coordinates:
            call color(igreen())
            call contourl(potave(1,1),cworka,NRFULL+1,nrhere,nthhere,
     $           zclv(ncont),icl,z2rho(1,1),x2rho(1,1),ntype)
            call color(iskyblue())
c            call contourl(phi(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
c     $           zclv(1),icl,zrho,xrho,ntype)
            call contourl(potave(1,1),cworka,NRFULL+1,nrhere,nthhere,
     $           zclv(1),icl,z2rho(1,1),x2rho(1,1),ntype)
            call color(15)
         endif
c      endif
      call legendline(0.47,1.07,258,'!Af!@'//char(0))
      call axis()
      call axlabels('z','r sin!Aq!@')
      
c     call pltend()
      end


c***************************************************************************
c Contouring of the density averaged over Psi, on distorted mesh.
      subroutine conrho2D(ir,it,nrhere,nthhere,npsihere,v1,larrows
     $     ,lconline,lpcic,ledge)
      integer ir,it
      real rhomax,rhomin,v1
      logical larrows,lconline,lpcic,ledge
c Common data:
      include 'piccompost.f'
c      include 'cic/piccompost.f'
c      save
c      character*20 cstring
      character*30 tstring
      character cworka(nr*(NTHFULL+1+1))
      integer ncont
      parameter (ncont=12)
      real zclv(2*ncont)
      real xrho(NRFULL+1,0:NTHFULL+1),zrho(NRFULL+1,0:NTHFULL+1)
      real x2rho(NRFULL+1,0:NTHFULL+1),z2rho(NRFULL+1,0:NTHFULL+1)
      real rhoave(NRFULL+1,0:NTHFULL+1)
c      real Tr(NRFULL+1,0:NTHFULL+1),Ttp(NRFULL+1,0:NTHFULL+1)
      save xrho,zrho
      real basesize
      parameter (basesize=.02)

      if(nthhere.gt.NTHFULL)then
         write(*,*)' Condisplay error. Mesh required:',nrhere,nthhere,
     $        ' Exceeds allocated:',NRFULL,NTHFULL
         stop
      endif

      do i=1,nrhere
         do j=0,nthhere+1
            rhoave(i,j)=0.
            do k=1,npsihere
               rhoave(i,j)=rhoave(i,j)+rho(i,j,k)
            enddo
            rhoave(i,j)=rhoave(i,j)/npsihere
         enddo
      enddo      

      call minmax2(rhoave(1,1),nr+1,nrhere,nthhere,rhomin,rhomax)

      do i=1,nrhere
          do j=1,nthhere
            zrho(i,j)=rcc(i)*tcc(j)
            xrho(i,j)=rcc(i)*sqrt(1.-tcc(j)**2)
            z2rho(i,j)=zrho(i,j)
            x2rho(i,j)=xrho(i,j)            
         enddo
         zrho(i,0)=rcc(i)
         xrho(i,0)=0.
         zrho(i,nthhere+1)=-rcc(i)
         xrho(i,nthhere+1)=0.
         if(lpcic)then
            z2rho(i,1)=zrho(i,0)
            x2rho(i,1)=xrho(i,0)
            z2rho(i,nthhere)=zrho(i,nthhere+1)
            x2rho(i,nthhere)=xrho(i,nthhere+1)
         endif
      enddo
c Adjust the outside angle centers if necessary.
c Should never be necessary.
c      if(tcc(1).eq.1)then
c         tcc(1)=0.25*(3.+tcc(2))
c         tcc(0)=1.
c         tcc(nthhere)=0.25*(-3.+tcc(nthhere-1))
c         tcc(nthhere+1)=-1.
c      endif

c ledge fixing
      if (ledge) then
         rhomax=rhomax/2.
      endif

      zscale=log(1000.)/(ncont-1.)
      do j=1,ncont
c Logarithmic
         zclv(j)=rhomax + (rhomin-rhomax)*
     $        exp(zscale*float(ncont-j))/exp(zscale*float(ncont-1))
c Linear
         zclv(j)=(rhomax-rhomin)*(1.*(j-1)/float(ncont-1))+ rhomin
c         write(*,*) "zcle : ",zclv(j)
      enddo
c      write(*,*) "min,max",rhomin,rhomax

c      write(*,*)'Contours=',zclv
      icl=-ncont

c      rpmax=2.*rcc(nrhere)-rcc(nrhere-1)
      rpmax=rcc(nrhere)
c      write(*,*)rcc(nrhere),nrhere,nthhere,rpmax
c      call multiframe(2,2,3)

      call ticnumset(8)
      call pltinaspect(-rpmax,rpmax,0.,rpmax)
      call accisgradinit(-25000,00000,35000,130000,65000,130000)
      ntype=2+16+32
      call contourl(rhoave(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv,icl,zrho,xrho,ntype)
      call gradlegend(zclv(1),zclv(abs(icl)),
     $     .1,1.25,.9,1.25,0.02,.true.)
c Call a second time for contours, without the highest.
      tstring(1:1)=char(0)
         if(lconline)then
c If this is not very bipolar
            if(abs(rhomax)-abs(rhomin).gt.0.2*abs(rhomax-rhomin))then
c Rational logarithmic contours
               do indx=1,ncont-2,3
                  base=10.**(-2+(indx-1)/3)
                  zclv(ncont+indx)=-base
                  zclv(ncont+indx+1)=-base*2.
                  zclv(ncont+indx+2)=-base*5.
c     New positive contours.
                  zclv(ncont-indx)=base
                  zclv(ncont-(indx+1))=base*2.
                  zclv(ncont-(indx+2))=base*5.
               enddo
               zclv(ncont)=0.
            endif
            write(*,'(a,30f6.2)')'Contours=',zclv
            ntype=2
            icl=(ncont-1)
c            write(*,'(2f8.4)')
c     $           (zrho(nrhere,k),z2rho(nrhere,k),k=0,nthhere+1)
c            call contourl(phi(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
c     $           zclv(ncont),icl,zrho,xrho,ntype)
c Fixed contour plot avoiding the hacked coordinates:
            call color(igreen())
            call contourl(rhoave(1,1),cworka,NRFULL+1,nrhere,nthhere,
     $           zclv(ncont),icl,z2rho(1,1),x2rho(1,1),ntype)
            call color(iskyblue())
c            call contourl(phi(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
c     $           zclv(1),icl,zrho,xrho,ntype)
            call contourl(rhoave(1,1),cworka,NRFULL+1,nrhere,nthhere,
     $           zclv(1),icl,z2rho(1,1),x2rho(1,1),ntype)
            call color(15)
         endif
c      endif
         call legendline(0.47,1.07,258,'n/n!A!d;!d!@'//tstring)
         call axis()
         call axlabels('z','r sin!Aq!@')
      
c     call pltend()
      end
c***********************************************************************
c Contouring of the charge density at theta=pi/2 (r,psi)
      subroutine conrhoPsi(ir,jstepth,nrhere,nthhere,npsihere,v1
     $     ,larrows,lconline,lpcic,ledge)

      integer ir
      real rhomax,rhomin,v1
      logical larrows,lconline,lpcic,ledge
c Common data:
      include 'piccompost.f'
c      include 'cic/piccompost.f'
c      save
c      character*20 cstring
      character*30 tstring
      character cworka(NRFULL*NPSIFULL+13)
      integer ncont
      parameter (ncont=12)
      real zclv(2*ncont)
      real xrho(NRFULL+1,0:NPSIFULL),yrho(NRFULL+1,0:NPSIFULL)
      real x2rho(NRFULL+1,0:NPSIFULL),y2rho(NRFULL+1,0:NPSIFULL)
      integer nthhalf
      real rhoave(0:NRFULL,0:NPSIFULL)
      save xrho,yrho
      real basesize
      parameter (basesize=.02)

      
      nthhalf=nint(0.5*nthhere)

      do i=1,nrhere
         do j=1,npsihere
            xrho(i,j)=rcc(i)*cos(pcc(j))
            yrho(i,j)=rcc(i)*sin(pcc(j))
            rhoave(i,j)=rho(i,nthhalf,j)
         enddo
         xrho(i,npsihere+1)=xrho(i,1)
         yrho(i,npsihere+1)=yrho(i,1)
         rhoave(i,npsihere+1)=rhoave(i,1)
      enddo

      call minmax(rhoave(1,1),nrhere,rmi,rma)
      call minmax(rhoave(1,nint(0.5*npsihere)),nrhere,rn2,rx2)
      rmi=min(rmi,rn2)
      rma=max(rma,rx2)

      if(rma.gt.5.)then
         rma=5.
         rmi=0.
      endif
      if(Ezext.ne.0. .and. rma.gt.2.)then
         rma=2.
         rmi=0.
      endif
      if(rmi.lt.0.5)rmi=0.5
      rhomin=rmi
      rhomax=rma

c ledge fixing
      if(ledge) then
         rhomax=rhomax/2.
      endif


       zscale=log(1000.)/(ncont-1.)
       do j=1,ncont
c     Logarithmic
          zclv(j)=rhomax + (rhomin-rhomax)*
     $         exp(zscale*float(ncont-j))/exp(zscale*float(ncont-1))
c     Linear
          zclv(j)=(rhomax-rhomin)*(1.*(j-1)/float(ncont-1))+ rhomin
c     write(*,*) "zcle : ",zclv(j)
       enddo
      icl=-ncont
      rpmax=rcc(nrhere)

      call ticnumset(8)
      call pltinaspect(-rpmax,rpmax,-rpmax,rpmax)
c      call pltinit(-rpmax,rpmax,-rpmax,rpmax)
      call accisgradinit(-25000,00000,-35000,130000,65000,130000)
      ntype=2+16+32
      call contourl(rhoave(1,1),cworka,NRFULL+1,nrhere,npsihere+1,
     $        zclv,icl,xrho(1,1),yrho(1,1),ntype)

      call gradlegend(zclv(1),zclv(abs(icl)),
     $     .1,1.25,.9,1.25,0.02,.true.)



c Call a second time for contours, without the highest.
      tstring(1:1)=char(0)
         if(lconline)then
c     If this is not very bipolar
            if(abs(rhomax)-abs(rhomin).gt.0.2*abs(rhomax-rhomin))then
c     Rational logarithmic contours
               do indx=1,ncont-2,3
                  base=10.**(-2+(indx-1)/3)
                  zclv(ncont+indx)=-base
                  zclv(ncont+indx+1)=-base*2.
                  zclv(ncont+indx+2)=-base*5.
c     New positive contours.
                  zclv(ncont-indx)=base
                  zclv(ncont-(indx+1))=base*2.
                  zclv(ncont-(indx+2))=base*5.
               enddo
               zclv(ncont)=0.
            endif
            write(*,'(a,30f6.2)')'Contours=',zclv
            ntype=2
            icl=(ncont-1)
c            write(*,'(2f8.4)')
c     $           (zrho(nrhere,k),z2rho(nrhere,k),k=0,nthhere+1)
c            call contourl(phi(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
c     $           zclv(ncont),icl,zrho,xrho,ntype)
c Fixed contour plot avoiding the hacked coordinates:
            call color(igreen())
            call contourl(rhoave(1,1),cworka,NRFULL+1,nrhere,npsihere,
     $           zclv(ncont),icl,x2rho(1,1),y2rho(1,1),ntype)
            call color(iskyblue())
c            call contourl(phi(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
c     $           zclv(1),icl,zrho,xrho,ntype)
            call contourl(rhoave(1,1),cworka,NRFULL+1,nrhere,npsihere,
     $           zclv(1),icl,x2rho(1,1),y2rho(1,1),ntype)
            call color(15)
         endif
c      endif
         call legendline(0.47,1.07,258,'n/n!A!d;!d!@'//tstring)
         call axis()
         call axlabels('x','y')

      call color(15)
c      call pltend()
      end

c ******************************************************************
c subroutine to compute the precalculated table for fluxofangle

      subroutine precalc()
c     Precalculation functions
      include 'piccompost.f'

      afac=(napre-1.)/(cang(angsize)-cang(1))

      do j=1,napre
c     finding the theta precalculated mesh.
         ap=(j-1.)/afac+cang(1)
         angl=cang(1)
         angr=cang(angsize)
         ial=1
         iar=angsize
 200     if(iar-ial.le.1)goto 210
         iax=(iar+ial)/2
         angx=cang(iax)
         if(angx.ge.ap) then
            angl=angx
            ial=iax
         else
            angr=angx
            iar=iax
         endif
         goto 200
 210     continue
         iapre(j)=ial
      enddo

      end

c***********************************************************************
c Interpolate onto the angle mesh. Return nearest index, fraction in angf.
      integer function interpang(ang,angf)
      include 'piccompost.f'
      real ang,angf
      ial=iapre(1+int((ang-cang(1))*afac))
      angf=(ang-cang(ial))/(cang(ial+1)-cang(ial))
            
      if(angf.gt.1.)then
         if(ial+1.le.napre)then
            ial=ial+1
            angf=(ang-cang(ial))/(cang(ial+1)-cang(ial))
         else
            write(*,*)'INTERPANG error. ial, cang incorrect'
            write(*,*)ial,angf,ang
         endif
      endif
      interpang=ial
      end

c***************************************************************************
c Data reading subroutine
c Version 3.0, Reading outputs from SCEPTIC3D

      subroutine readoutput(lreaddiag,lpcic,ledge,
     $     filename,rholocal,nrhere,nthhere,npsihere,
     $     phipic,rhopic,rpic,rpicleft,phicos,
     $     rhomax,rhomin,
     $     nrti,phiinf,nastep,nsteps,
     $     dt,rmax,fave,debyelen,vprobe,
     $     icolntype,colnwt,Eneutral,vneutral,Tneutral,
     $     ierr)
      logical lreaddiag,lpcic,ledge
      character*100 string,filename
      real phipic(1000),rhopic(1000)
      real rpic(1000),rpicleft(1000),phicos(1000)
      include 'piccompost.f'
      real rholocal(0:NRFULL,0:NTHFULL,0:NPSIFULL)
      character*256 charin
      common /forces/ charge1,ffield1,felec1,fion1,ftot1,
     $     charge2,ffield2,felec2,fion2,ftot2

      ierr=0
c Read the data file:
c__________________________________________________________________

      open(10,file=filename,status='old',err=101)
c Line for nothing.
      read(10,*)charin
      read(10,'(a)')charin
c      write(*,*)charin
      read(charin,*,err=201,end=201) dt,vd,cd,Ti,isteps,rhoinf,phiinf
     $     ,fave,debyelen,vprobe,Bz ,icolntype,colnwt
 201  continue
      write(*,'(a,a)')'  dt    vd     cd    Ti     steps  rhoinf ' ,
     $       'phiinf  fave  debyelen Vp   Bz...'
      write(*
     $     ,'(2f7.4,f7.3,f7.3,i5,f8.1,f7.3,f8.4,f8.3,f8.3, f7.3,$)'
     $     ) dt,vd,cd,Ti,isteps,rhoinf,phiinf,fave,debyelen,vprobe
     $     ,Bz
      if(icolntype.gt.0)then
         write(*,'(i2,f7.3)',err=212)icolntype,colnwt
      else
         write(*,*)
      endif

  
 212  read(10,*,err=202)nrTi
      nrhere=nrTi
c      write(*,*)'nrTi=',nrTi
      do i=1,nrTi
         read(10,*,err=203)rpic(i),phipic(i),rhopic(i)
      enddo
      read(10,*)charin
c      write(*,*)charin(1:78)
      read(10,*,err=204)nsteps
c      write(*,*)nsteps
      if(nsteps.gt.nstepmax) then
         write(*,*)'Number of steps',nsteps,
     $        ' exceeds allocation',nstepmax
         call exit
      endif

      read(10,*)(fluxprobe(j),j=1,nsteps)



c      write(*,*)(fluxprobe(j),j=1,nsteps)
c Read theta cells
      read(10,*)charin

c      write(*,*)charin(1:78)
      read(10,*)nthhere,npsihere,nsteps
c      write(*,*) nthhere,npsihere,nsteps
      
c      write(*,*)nthhere,npsihere,nsteps
c      do i=1,nsteps
c         read(10,*)(ninthstep(j,i),j=1,nthhere)
c      enddo

      read(10,*)charin
c      write(*,*)charin(1:78)
      read(10,*)nastep
c      write(*,*)'nastep',nastep
      nread=nthhere
c This is not necessary here and indeed breaks the combined read because
c lpcic has not yet been set.
c      if(.not.lpcic)nread=nread+1

      read(10,*)((nincell(j,l),j=1,nread),l=1,npsihere)
      if(lreaddiag)then
         write(*,*)'nastep',nastep,' nincell:'
         do l=1,npsihere
            write(*,*)(nincell(j,l),j=1,nthhere)
         enddo
      endif
      read(10,*)charin
      do j=1,nrhere
         read(10,'(10f8.3)')((phi(j,k,l),k=1,nthhere),l=1,npsihere)
      enddo
      read(10,*)charin
      do j=1,nrhere
         read(10,'(10f8.3)')((rholocal(j,k,l),k=1,nthhere),l=1,npsihere)
      enddo
      if(lreaddiag)then
         write(*,*)' rhoprobe:'
         do l=1,npsihere
            write(*,*)(rholocal(1,j,l),j=1,nthhere)
         enddo
      endif

      read(10,*)charin
      read(10,*)(volinv(k),k=1,nrhere)
c We don't use open and close for combined files.
      if(filename(1:2).eq.'Ti')then
c But this must be using the split version.
         close(10)
c Read in  summed results.
         filename(1:2)='Sp'
         open(10,file=filename,err=210,status='old')
      endif
      read(10,'(a)')string
      read(10,*,err=200)
     $     dt,vd,Ti,i,rmax,rhoinf,debyelen,vprobe
 200  continue
      if(rhoinf.lt.1.)rhoinf=exp(phiinf)
      read(10,*)nrhere,nthhere,npsihere
      if(nrhere.gt.NRUSED .or. nthhere.gt.NTHUSED)then
         write(*,*)'Required dimensions: nr',nrhere,' nth',nthhere
         write(*,*)'are too large for the allocated values:'
     $        ,NRUSED,NTHUSED
         stop
      endif

      

      read(10,*)string
      read(10,*)(((psum(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      read(10,*)string
      read(10,*)(((vrsum(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      read(10,*)string
      read(10,*)(((vtsum(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      read(10,*)string
      read(10,*)(((vpsum(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      read(10,*)string
      read(10,*)(((v2sum(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      read(10,*)string
      read(10,*)(((vr2sum(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      read(10,*)string
      read(10,*)(((vtp2sum(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      read(10,*)string
      read(10,*)(((vxsum(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      read(10,*)string
      read(10,*)(((vysum(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      read(10,*)string
      read(10,*)(((vzsum(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      read(10,*)string

      read(10,*)(rcc(k1),k1=1,nrhere)
      read(10,*)string
      read(10,*)(volinv(k1),k1=1,nrhere)
      read(10,*)string
      read(10,*)(tcc(k2),k2=1,nthhere)      
      read(10,*,err=402,end=402)string
      read(10,*)charge1,ffield1,felec1,fion1,ftot1
      read(10,*)charge2,ffield2,felec2,fion2,ftot2
      read(10,*,err=402,end=402)string
      read(10,*,err=402,end=402)
     $     icolntype,colwt,Eneutral,vneutral,Tneutral
 402  close(10)
      write(*,*)'nrhere,nthhere,npsihere,icolntype,colwt'
      write(*,*)nrhere,nthhere,npsihere,icolntype,colwt
      if(lreaddiag)then
         write(*,*)'Finished reading'
         write(*,*)'vrsum(1)'
         write(*,501)((vrsum(1,k2,k3), k2=1,nthhere),k3=1,npsihere)
         write(*,*)'psum(1)'
         write(*,501)((psum(1,k2,k3), k2=1,nthhere),k3=1,npsihere)
         write(*,*)'vr(1)'
         write(*,501)((vrsum(1,k2,k3)/psum(1,k2,k3), k2=1,nthhere),k3=1
     $        ,npsihere)
         write(*,*)'rcc'
         write(*,*)(rcc(k1),k1=1,nrhere)
         write(*,*)'volinv'
         write(*,*)(volinv(k1),k1=1,nrhere)
         write(*,*)'tcc'
         write(*,*)(tcc(k2),k2=1,nthhere)
      endif
 501  format(10f8.3)
c__________________________________________________________________
c End of reading the data section

c__________________________________________________________________
c Fix up data 
      if(tcc(1).eq.1)lpcic=.true.
c Correct the outside angle centers if necessary.
      if(lpcic)then
         tcc(1)=0.25*(3.+tcc(2))
         tcc(0)=1.
         tcc(nthhere)=0.25*(-3.+tcc(nthhere-1))
         tcc(nthhere+1)=-1.

         thcells=nthhere-1.
      else
         thcells=nthhere
      endif


c In postproc, th is used as the weighting of the cells, essentially 
c the delta of cos theta that corresponds to each cell.
      do i=1,nthhere
         th(i)=2./(thcells)
      enddo
c The ends are half-weighted for cic.
      if(lpcic)then
         th(1)=0.5*th(1)
         th(nthhere)=0.5*th(nthhere)
      endif


c Calculate rho(i,j) from the psum, volinv. Normalized by rhoinf.
      rhomax=0.
      do k3=1,npsihere
         do k1=1,nrhere
            do k2=1,nthhere
               rho(k1,k2,k3)=psum(k1,k2,k3)*volinv(k1)*thcells*npsihere
     $              /rhoinf
c     For now, use rholocal to fix --ds problem, maybe wrong ?
               if(rholocal(k1,k2,k3).gt.rhomax)rhomax=rholocal(k1,k2,k3)
            enddo
            if(lpcic)then
c fix up rho as double on boundary.
               rho(k1,1,k3)=2.*rho(k1,1,k3)
               rho(k1,nthhere,k3)=2.*rho(k1,nthhere,k3)
            endif
c fix angle ends of rho and psi
            rho(k1,0,k3)=rho(k1,1,k3)
            rho(k1,nthhere+1,k3)=rho(k1,nthhere,k3)
            phi(k1,0,k3)=phi(k1,1,k3)
            phi(k1,nthhere+1,k3)=phi(k1,nthhere,k3)
            rholocal(k1,0,k3)=rholocal(k1,1,k3)
            rholocal(k1,nthhere+1,k3)=rholocal(k1,nthhere,k3)
         enddo
      enddo
      ir=10.
      if(ledge)then
         rhomax=min(rhomax,1.5)
         rhomin=.5
      else
         rhomin=0.
      endif

      if(lreaddiag)then
         write(*,*)'rho   ','rholocal',
     $     ' ;  rho is from psum, rholocal from Ti file'
         do i=1,nrhere
            do k=1,npsihere
               write(*,*)rho(i,1,k),rholocal(i,1,k),rho(i,1,k)
     $              /rholocal(i,1,k)
            enddo
         enddo
         write(*,*)'End of rho comparison'
      endif

      phiinf=0.

      jmin=1
      jmax=nthhere
      if(lreaddiag)write(*,*)'jmin,jmax',jmin,jmax
      do i=1,nrhere
         rhopic(i)=0.
         phicos(i)=0.
c         write(*,*)'th   tcc   phi'
         do j=jmin,jmax
            do k=1,npsihere
               rhopic(i)=rhopic(i)+rholocal(i,j,k)
c     rhopic(i)=rhopic(i)+rho(i,j)
c \int cos(\theta) \phi(\theta) d\cos(\theta)
               phicos(i)=phicos(i)+th(j)*tcc(j)*phi(i,j,k)
c            write(*,'(4f10.4)')th(j),tcc(j),phi(i,j),phicos(i)
            enddo
         enddo
         rhopic(i)=rhopic(i)/float(jmax-jmin+1)/npsihere
      enddo

c     rescale rho; but usually this is the identity transformation.
      do i=1,nrhere
         rpicleft(i)=-rpic(i)
      enddo

      return
c End of data fix-up section


c__________________________________________________________________
 202  write(*,*)"nr error"
      call exit
 203  write(*,*)"rpicphipic error"
      call exit
 204  write(*,*)"nsteps error"
      call exit
 210  write(*,*)'Error opening file: ',filename(1:50)
      call exit
 101  ierr=101
      end

