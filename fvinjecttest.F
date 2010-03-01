      program fvinjecttest
c To test fvinject
c Grabbed initially from sceptic.F
      logical success,finaldiags
      real rmax
      character*100 string
      logical readpart,writepart
      logical lsmooth,lcolcont
c Common storage:
      include 'piccom.f'
      include 'fvcom.f'
      real ftest(nxfvi:nxfva,nzfvi:nzfva)

c diagnostics
      real vydist(nxfvi:nxfva),fmaxwell(nxfvi:nxfva),vzdist(nxfvi:nxfva)
      common /vinjdiag/vydist,vzdist
c Local diagnostics
      integer ndist
      parameter (ndist=40)
      real vphidist(-ndist:ndist),vphi(-ndist:ndist)
     $     ,fmdist(-ndist:ndist)
      real phiangle(npartmax), pindex(npartmax),tanvel(npartmax)
     $     ,cosarr(npartmax)
      parameter (ncdist=20)
      real cosdist(ncdist),cdangle(ncdist)
      real cosangle(ncdist),acdist(ncdist),ratc(ncdist)
      real acint(ncdist)
      real qthtest(nthfvsize),fvthtest(nthfvsize)
c      real worka(0:nxfva-nxfvi+2,0:nzfva-nzfvi+2)
c Data
      data rmax/5./
      data readpart/.false./writepart/.false./

c Diagnostic array
      vmax=4.
      do j=-ndist,ndist
         vphi(j)=j*vmax/ndist
         vphidist(j)=0.
      enddo

c Full size arrays by default. Can be changed later by switches.
      nr=nrsize
      nth=nthfvsize
c Common defaults. Avoid block data.
      success=.false.
      vd=.000
      debyelen=.1
      vprobe=-4.
      Ezext=0.
      bdyfc=4.
      Ti=1.
      diags=.false.
      finaldiags=.false.
      ldist=.false.
      avelim=0.6
      nstepsave=0
      nsamax=40
      lsmooth=.false.
      lcolcont=.true.
      ircell=1
      itcell=1
      norbits=0
      npart=1
      ninjcomp0=npartmax
      ninjcomp=ninjcomp0
      linsulate=.false.
      lfloat=.false.
      rmtoz=1.
      ierad=1
      rhomin=0.
      rhomax=0.
      lat0=.false.
      lfext=.false.
      colnwt=0.
      icolntype=0
      vneutral=0.
      localinj=.false.
      damplen=0.
      ipfsw=3
      call pfset(0)
      lfixedn=.true.
      ipstep=1

c Deal with arguments. Lots of these are spurious carry-overs.
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
         if(string(1:2) .eq. '-f') finaldiags=.true.
         if(string(1:2) .eq. '-s') diags=.true.
         if(string(1:3) .eq. '-er') then
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
         if(string(1:2) .eq. '-b') read(string(3:),*)bdyfc
         if(string(1:2) .eq. '-t') read(string(3:),*)Ti
         if(string(1:2) .eq. '-x') read(string(3:),*)rmax
         if(string(1:2) .eq. '-v') read(string(3:),*)vd
         if(string(1:3) .eq. '-pf') then
            lfloat=.true.
         elseif(string(1:3) .eq. '-pi') then
            linsulate=.true.
         elseif(string(1:2) .eq. '-p') then
            read(string(3:),*)vprobe
         endif
         if(string(1:2) .eq. '-z') read(string(3:),*)Ezext
         if(string(1:2) .eq. '-l') read(string(3:),*)debyelen
         if(string(1:3) .eq. '-nr') read(string(4:),*)nr
c         if(string(1:3) .eq. '-nt') read(string(4:),*)nth
c Can't change here.
         if(string(1:3) .eq. '-ni') read(string(4:),*)npart
         if(string(1:3) .eq. '-kt') then
            read(string(4:),*)icolntype
         elseif(string(1:2) .eq. '-k') then
            read(string(3:),*)colnwt
         endif
         if(string(1:2) .eq. '-m') read(string(3:),*)myid
         if(string(1:2) .eq. '-a') then
            read(string(3:),*,err=251,end=251)ipfsw
            call pfset(ipfsw)
            goto 252
 251        ipfsw=0
 252        continue
         endif
         if(string(1:2) .eq. '-?') goto 51
 1    continue
 3    continue

c      write(*,*)linsulate,diags

c Set Array sizes, allowed variable.
      if(nr.gt.nrsize)then
         write(*,*)'Too many radial points:',nr,'  Set to',nrsize
         nr=nrsize
      endif
      if(npart.gt.npartmax)then
         write(*,*)'Too many ions:',npart,'  Set to',npartmax
         npart=npartmax
      endif
      write(*,505)nr,nthfvsize,npart
 505  format(' Mesh: ',i3,' x',i3,'  Particles:',i7)

      r(nr)=rmax
      Eneutral=colnwt*(vd-vneutral)
      write(*,*)'vd,icolntype,colnwt,Eneutral',
     $     vd,icolntype,colnwt,Eneutral
c********************
c Real init testing section.
c
      ud=vd/sqrt(2.*Ti)
      actot=0.
      acdist(1)=0.
      acint(1)=0.
      cosangle(1)=-1.
      do j=2,ncdist
         cosdist(j)=0.
         cdangle(j) =-1.+2.*(j-1.5)/(ncdist-1)
         acdist(j)=canglecx(-cdangle(j),ud)
         actot=actot+acdist(j)
         cosangle(j)=-1.+2.*(j-1)/(ncdist-1)
c         acint(j)=acint(j-1)+0.5*(acdist(j)+acdist(j-1))
         acint(j)=acint(j-1)+acdist(j)
     $        *(cosangle(j)-cosangle(j-1))
      enddo
c Renormalize to sum equals 1.
      do j=2,ncdist
         acdist(j)=0.5*(ncdist-1)*acdist(j)/actot
      enddo
c Now actot is \int canglecx dcostheta.
      ldiaginj=diags
      write(*,*) 'Testing fvinject. Drift velocity=',vd
      call fvinjinit()
      write(*,*) 'Finished initialization on grid',
     $     nthfvsize,nxfva,nzfva
      write(*,*)'acint,qthfv(nthfvsize)',acint(ncdist),qthfv(nthfvsize)


 510  format(10f8.3)

      if(nthfvsize.le.10)then
         write(string,'(''('',i3,''f8.3)'')')nthfvsize
         write(*,*)'qxfv(ix,ith)='
         write(*,string)((qxfv(i,j),j=1,nthfvsize),i=nxfvi,nxfva)
         write(*,*)'nthfvsize=',nthfvsize,' vx=',nxfvi,nxfva
         write(*,*)'     th     qxfv(nxfva,ith)   qthfv    qfv(1,1,ith)'
         write(*,'(4f12.4)')(fvth(i),qxfv(nxfva,i),qthfv(i),
     $        qfv(1,1,i),
     $        i=1,nthfvsize)
      endif
c         write(*,*)'ztrfv(ix,ith)='
c         write(*,string)((qxfv(i,j),j=1,nthfvsize),i=nxfvi,nxfva)

c Test the function evaluation
      do i=nxfvi,nxfva
         do j=nzfvi,nzfva
            ftest(i,j)=fv(vxfv(i),vzfv(j))
c            write(*,'(a,3f10.5)')'f-value at vx,vy',
c     $           ftest(i,j),vxfv(i),vzfv(j)
         enddo
      enddo
c Plot the function.
      if(lfloat) then
 100     call pltinit(0.,1.,0.,1.)
         j=1 + 256*10 + 256*256*7
         call hidweb(vxfv(nxfvi),vzfv(nzfvi),ftest,
     $        nxfva-nxfvi+1,nxfva-nxfvi+1,nzfva-nzfvi+1,j)
         call boxtitle('f(vx,vz)')
         call ax3labels('vx','vz','f(vx,vz)')
         call axident3()
         if(ieye3d().ne.0) goto 100
c plot ztrfv
         call pltinit(vxfv(nxfvi),vxfv(nxfva),float(nzfvi),float(nzfva))
         call axis()
         call axlabels('vxfv','ztrfv')
         do i=1,nthfvsize
            call polyline(vxfv,ztrfv(nxfvi,i),nxfva-nxfvi+1)
            call polymark(vxfv,ztrfv(nxfvi,i),nxfva-nxfvi+1,mod(i,16))
         enddo
         call pltend()
      endif

c Plot the integrals.
      if(linsulate) then
 600     call pltinit(0.,1.,0.,1.)
         jsw=1 + 256*10 + 256*256*7
         call hidweb(vxfv,fvth(1),qxfv,
     $        nxfva-nxfvi+1,nxfva-nxfvi+1,nthfvsize,jsw)
         call boxtitle('qxfv(vx,th)')
         call axident3()
         call ax3labels('vx','th','qxfv(vx,th)')
         if(ieye3d().ne.0) goto 600
         do j=1,nthfvsize
c qfv for all the thetas  
            call minmax2(qfv(nzfvi,nxfvi,j),nzfva-nzfvi+1,nzfva-nzfvi+1,
     $           nzfva-nzfvi+1,qmin,qmax)
            if(qmin.ne.qmax)then
               write(*,*)'j,th',j,fvth(j)
 601           call pltinit(0.,1.,0.,1.)
               jsw=1 + 256*10 + 256*256*7
               call hidweb(vzfv,vxfv,qfv(nzfvi,nxfvi,j),
     $              nzfva-nzfvi+1,nzfva-nzfvi+1,nxfva-nxfvi+1,jsw)
               call boxtitle('qfv(vz,vx)')
               call axident3()
               call ax3labels('vz','vx','qfv(vz,vx)')
               if(ieye3d().ne.0) goto 601
            else
               write(*,*)'Flat qfv',j,qmin,qmax
            endif
         enddo

c         endif
 620     call pltinit(0.,1.,0.,1.)
         jsw=1 + 256*10 + 256*256*7
         call hidweb(vxfv,vzfv,fqvxvz,
     $        nxfva-nxfvi+1,nxfva-nxfvi+1,nzfva-nzfvi+1,jsw)
         call boxtitle('fqvxvz(vx,vz)'//
     $        ' tangential cumulative distribution')
         call axident3()
         call ax3labels('vx','vz','fqvxvz(vx,vz)')
         if(ieye3d().ne.0) goto 620

      endif

      if(npart.le.100)
     $     write(*,*)'Reinjection particle No, position, velocity'
      dt=.0001
      do i=1,npart
         call fvreinject(i,dt,icolntype)
         if(npart.le.100) write(*,701)i,(xp(k,i),k=1,6)
 701     format(i6,6f10.4)
         phiangle(i)=atan2(xp(5,i),xp(4,i))
         pindex(i)=i
         vp=(xp(1,i)*xp(5,i)-xp(2,i)*xp(4,i))/
     $        sqrt(xp(1,i)**2+xp(2,i)**2)
         tanvel(i)=vp
         ip=nint(ndist*vp/vmax)
         if(ip.lt.-ndist)ip=-ndist
         if(ip.gt.ndist)ip=ndist
         vphidist(ip)=vphidist(ip)+1.
         ca=xp(3,i)/sqrt(xp(1,i)**2+xp(2,i)**2+xp(3,i)**2)
         cosarr(i)=ca
         icdist=nint(1.5+(ncdist-1)*(ca+1)/2.000001)
         cosdist(icdist)=cosdist(icdist)+1
      enddo
      itotal=0
      aver=0
      am1=0
      do j=2,ncdist
         itotal=itotal+cosdist(j)
         cosdist(j)=0.5*(ncdist-1)*cosdist(j)/npart
         ratc(j)=cosdist(j)/acdist(j)
         aver=aver+ratc(j)
         am1=am1+ratc(j)*cosangle(j)
      enddo
      aver=aver/ncdist
      am1=am1/ncdist
      write(*,*)'npart,itotal,average,1st moment',npart,itotal,aver,am1

c Plot some diagnostics      
c The fvinject theta distribution function. Includes the poles.
      do j=1,nthfvsize
         qthtest(j)=qxfv(nxfva,j)
         fvthtest(j)=-fvth(j)
      enddo
c Integrate:
      total=0.
      do j=1,nthfvsize-1
         total=total+
     $        (fvthtest(j+1)-fvthtest(j))*(qthtest(j+1)+qthtest(j))*0.5
      enddo
      do j=1,nthfvsize
         qthtest(j)=qthtest(j)/abs(total)
      enddo
c This then seems to give excellent agreement with analytic form.
c And provided dt is small enough, so does the reinjection.

      call autoplot(cosangle,acint,ncdist)
      call axlabels('cos(!Aq!@)','Cumulative angle distribution')
      call color(4)
      call dashset(1)
      call polyline(fvthtest,qthfv,nthfvsize)
      call pltend()

      write(*,'(2f10.5)')(cdangle(j),acdist(j),j=1,ncdist)

      call automark(cdangle(2),cosdist(2),ncdist-1,1)
      call axlabels('cos(!Aq!@)','Injection angle distribution')
      call polyline(cdangle(2),acdist(2),ncdist-1)
      call color(4)
      call dashset(1)
      call polyline(fvthtest,qthtest,nthfvsize)
      call pltend()

      call automark(cosangle,ratc,ncdist,3)
      call axlabels('cos(!Aq!@)','Distribution ratio')
      call pltend()

c      call automark(pindex,phiangle,npart,10)
c      call axlabels('index','phiangle')
c      call pltend()
c      call automark(pindex,tanvel,npart,10)
c      call axlabels('index','azimuthal velocity')
c      call pltend()
      call autoplot(vphi,vphidist,ndist+ndist+1)
      call axlabels('v','Azimuthal velocity distribution')
      call minmax(vphidist,ndist+ndist+1,vphimin,vphimax)
      do j=-ndist,ndist
         v=vphi(j)
c velocity is in Te units; convert to those.
         fmdist(j)=exp(-v**2/(2.*Ti))*vphimax
      enddo
      call dashset(1)
      call polyline(vphi,fmdist,ndist+ndist+1)
      call dashset(0)
      call pltend()
c
      call autoplot(vxfv,vydist,nxfva-nxfvi+1)
      call axlabels('v','y-velocity distribution')
      call minmax(vydist,nxfva-nxfvi+1,vymin,vymax)
      do j=nxfvi,nxfva
         v=vxfv(j)
         fmaxwell(j)=exp(-v**2)*vymax
      enddo
      call dashset(1)
      call polyline(vxfv,fmaxwell,nxfva-nxfvi+1)
      call dashset(0)
      call pltend()
c
      call autoplot(vzfv,vzdist,nxfva-nxfvi+1)
      call axlabels('v','z-velocity distribution')
      call minmax(vzdist,nxfva-nxfvi+1,vzmin,vzmax)
      do j=nxfvi,nxfva
         v=vzfv(j)
         fmaxwell(j)=fvcx(v,vd/sqrt(2.*Ti))*npart*(vzfv(2)-vzfv(1))
      enddo
      call dashset(1)
      call polyline(vzfv,fmaxwell,nxfva-nxfvi+1)
      call dashset(0)
      call jdrwstr(.02,.02,
     $'Comparing injection flux to f(vz) distrib. Should not agree',
     $     1.)
      call pltend()

      call exit(0)
 51   write(*,*) 'Usage: fvinjecttest -vfff.ff [other options]'
      write(*,*) ' Test the fvinjection routine specifying drift vd'
      write(*,*) ' -pf plot f(vx,vz). -pi plot integrated f (qxfv).'
      write(*,*) ' -nrii -ntii set r- and theta- array size ii.'
      write(*,*) ' -aii set plotfile ii. -mii set myid (seed)'
      write(*,*) ' -s give solving diagnostics. -ni injections.'
      end
c******************************************************************

c Notes
c 8 Oct 05. Got it apparently working. Doing time tests.
c With nrxnt=90x90, and vx x vz =50 x 50
c  Initialization approximately 0.6s
c  1000000 reinjections approximately 4s. Actually may be just 200000.
c  Probably small enough to be negligible.

c***************************************************************
      function canglecx(ct,ud)
c Return the angle distribution function for a cx constant collision
c frequency distribution when the normalized velocity is ud, at angle
c cos(theta)=ct
      real ct,ud
      
      x=abs(2.*ct*ud)
      if(x.lt.1.e-10)x=1.e-10
      canglecx=0.5*(ct*ud*(1.+ sign(experfcc(1/x),ct*ud))
     $     + 1./sqrt(3.1415926))
      end
