
c*********************************************************************
c     Writes the main output file
      subroutine output(dt,i,fave,icolntype,colnwt)
c Common data:
      include 'piccom.f'
      include 'colncom.f'
      character*55 filename
c      integer iti,it2
c Construct a filename that contains many parameters
c Using the routines in strings_names.f
      filename=' '
      call nameappendexp(filename,'T',Ti,1)
      call nameappendint(filename,'v',nint(100*vd),3)
      if(cd.ne.1)then
         call nameappendint(filename,'a',nint(10*cd),2)
      endif
      if(cB.ne.1)then
         call nameappendint(filename,'b',nint(10*cB),2)
      endif
      if(r(nr).ge.100)then
         call nameappendint(filename,'R',ifix(r(nr)/10.),2)
      else
         call nameappendint(filename,'r',ifix(r(nr)),2)
      endif
      call nameappendint(filename,'P',ifix(abs(Vprobe)),3)
      if (infdbl) then
         call nameappendexp(filename,'LI',debyelen,1)
      else
         call nameappendexp(filename,'L',debyelen,1)
      endif
      if(Bz.ne.0.) call nameappendexp(filename,'B',Bz,2)
      if(icolntype.eq.1.or.icolntype.eq.5)
     $     call nameappendexp(filename,'c',colnwt,1)
      if(icolntype.eq.2.or.icolntype.eq.6)
     $     call nameappendexp(filename,'C',colnwt,1)
     
      if(vneutral.ne.0)
     $     call nameappendint(filename,'N',nint(100*vneutral),3)

c Turn those on when we test convergence with changing grid size
c      call nameappendint(filename,'Nr',nrused,3)
c      call nameappendint(filename,'Nt',nthused,3)
c      call nameappendint(filename,'Np',npsiused,3)

      idf=nbcat(filename,'.dat')
c Write out averaged results.
      open(10,file=filename)
      write(10,'(a,a)')'  dt    vd    cd      Ti     steps  rhoinf ' ,
     $       'phiinf  fave  debyelen Vp   Bz  cB ...'
      write(10,'(2f7.4,f7.3,f7.3,i5,f8.1,f7.3,f8.4,f12.4,f9.3,
     $     f7.3,f7.3,$)') dt,vd,cd,Ti,i,rhoinf,log(rhoinf),fave,debyelen
     $     ,vprobe,Bz,cB

      if(icolntype.gt.0)then
         write(10,'(i2,f8.4)')icolntype,colnwt
      else
         write(10,*)
      endif
      write(10,*)NRUSED
      do j=1,NRUSED
         write(10,*)rcc(j),diagphi(j),diagrho(j)
c -log(rhoinf)
      enddo
      write(10,'(a)')'Number of steps, Particles to probe each step'
      write(10,*)i
      write(10,*)(fluxprobe(j),j=1,i)
      write(10,'(a,a)')'Number of theta cells, Number of psi cells,'
     $     ,' Number of steps'
      write(10,*)NTHUSED,NPSIUSED,i
      do j=1,NTHUSED
         do k=1,NPSIUSED
            nincell(j,k)=0
            vrincell(j,k)=0
            vr2incell(j,k)=0
         enddo
      enddo
      nastep=0
      do k=1,i
c     For the 3D version, do not write the history of particle
c     collection, because it would take too much space
c         write(10,*)(nincellstep(j,k),j=1,NTHUSED)
c     Just save the last quarter for the average
         if(k.gt.3*i/4)then
            nastep=nastep+1
            do j=1,NTHUSED
               do l=1,NPSIUSED
                  nincell(j,l)=nincell(j,l)+nincellstep(j,l,k)
                  vrincell(j,l)=vrincell(j,l)+vrincellstep(j,l,k)
                  vr2incell(j,l)=vr2incell(j,l)+vr2incellstep(j,l,k)
               enddo
            enddo
         endif
      enddo
      write(10,'(a,a)')'Particle angular distrib summed over last'
     $     ,' quarter of steps, numbering:'
      write(10,*)nastep
      write(10,*)((nincell(j,l),j=1,NTHUSED),l=1,NPSIUSED)
      write(10,'(a,a)')'Particle angular distrib*vr summed over last'
     $     ,' quarter of steps, numbering:'
      write(10,*)((vrincell(j,l),j=1,NTHUSED),l=1,NPSIUSED)
      write(10,'(a,a)')'Particle angular distrib*vr^2 summed over last'
     $     ,' quarter of steps, numbering:'
      write(10,*)((vr2incell(j,l),j=1,NTHUSED),l=1,NPSIUSED)

      write(10,'(a,i4,i4,i4)')'Mesh potential. Grid',NRUSED,NTHUSED
     $     ,NPSIUSED
c      call minmax2(phi(1,1),nrsize+1,NRUSED,NTHUSED,phimin,phimax)
c      if(max(abs(phimin),abs(phimax)).lt.1.)then
c         do j=1,NRUSED
c            write(10,'(10f8.4)')(phi(j,k),k=1,NTHUSED)
c         enddo
c      else

         do j=1,NRUSED
            write(10,'(10f8.3)')((phi(j,k,l),k=1,NTHUSED),l=1,NPSIUSED)
         enddo

c      endif
         write(10,'(a,i4,i4,i4)')'Mesh density/infinity. Grid',NRUSED
     $        ,NTHUSED,NPSIUSED
      do j=1,NRUSED
         write(10,'(10f8.3)')((rhoDiag(j,k,l),k=1,NTHUSED),l=1,NPSIUSED)
      enddo
      write(10,'(a,i4,i4)')'Volinv. Grid',NRUSED
      write(10,'(10f8.3)')(volinv(k),k=1,NRUSED)
      write(10,*)'rhoDiag'
c      do j=1,NRUSED
c         write(10,'(10f8.3)')((rhoDiag(j,k,l),k=1,NTHUSED),l=1,NPSIUSED)
c      enddo

      call outsums(dt,i+1)

c Output time-averages of z-force components stored in zmom(nstepmax,*,*).
c Particle units nTr^2, Electric nT lambda_D^2.
      ztotal1=zmom(nstepmax,fieldz,1)*debyelen**2
     $     +zmom(nstepmax,epressz,1)+zmom(nstepmax,partz,1)
     $     +zmom(nstepmax,lorentz,1)
      ztotal2=zmom(nstepmax,fieldz,2)*debyelen**2
     $     +zmom(nstepmax,epressz,2)+zmom(nstepmax,partz,2)
     $     +zmom(nstepmax,lorentz,2)
      xtotal1=xmom(nstepmax,fieldz,1)*debyelen**2
     $     +xmom(nstepmax,epressz,1)+xmom(nstepmax,partz,1)
     $     +xmom(nstepmax,lorentz,1)
      xtotal2=xmom(nstepmax,fieldz,2)*debyelen**2
     $     +xmom(nstepmax,epressz,2)+xmom(nstepmax,partz,2)
     $     +xmom(nstepmax,lorentz,2)
      ytotal1=ymom(nstepmax,fieldz,1)*debyelen**2
     $     +ymom(nstepmax,epressz,1)+ymom(nstepmax,partz,1)
     $     +ymom(nstepmax,lorentz,1)
      ytotal2=ymom(nstepmax,fieldz,2)*debyelen**2
     $     +ymom(nstepmax,epressz,2)+ymom(nstepmax,partz,2)
     $     +ymom(nstepmax,lorentz,2)
      
      write(10,*)'Charge     z E-field     z  Electrons',
     $     '     z Ions    z Lorentz  z Total'
      write(10,*)(zmom(nstepmax,j,1),j=1,5),ztotal1
      write(10,*)(zmom(nstepmax,j,2),j=1,5),ztotal2

      write(10,*)'x E-field     x  Electrons',
     $     '     x Ions    x Lorentz  x Total'
      write(10,*)(xmom(nstepmax,j,1),j=2,5),xtotal1
      write(10,*)(xmom(nstepmax,j,2),j=2,5),xtotal2
      
      write(10,*)'y E-field     y  Electrons',
     $     '     y Ions    y Lorentz  y Total'
      write(10,*)(ymom(nstepmax,j,1),j=2,5),ytotal1
      write(10,*)(ymom(nstepmax,j,2),j=2,5),ytotal2

      if(rmtoz.ne.1.) write(10,'(''rmtoz='',f10.4)')rmtoz
      write(10,*)'Collisions: Type,Weight,Eneutral,vneutral,Tneutral'
      if(icolntype.ne.0) write(10,701) icolntype,colnwt,Eneutral
     $     ,vneutral,Tneutral
      write(10,*) 'Energy flux to the probe'
      write(10,*) enertot(nstepmax)
      
 701  format(10x,i3,4f10.5)
c     701  format('Collisions: type=',i4,' weight=',f8.4,' Eneutral=',
c     $     f10.5,' vneutral=',f8.4,' Tneutral=',f8.4)

c End of output file.
      close(10)

c     For the 3D version, remove the outforce file. This has never been
c     used anyways

c      innm=lentrim(filename)
c      filename(innm-3:innm)='.frc'
c      call outforce(filename,i)
      end
c************************************************************************
c     Writes a txt file with the orbits of the traced particles
      subroutine orbitoutput()
c     Common data
      include 'piccom.f'

      character*30 filename
      integer iti,it2

c Construct a filename that contains many parameters

      write(filename,'(a)')'T'
      iti=nint(alog10(Ti)-0.49)
      it2=nint(Ti/10.**iti)
      write(filename(2:2),'(i1.1)')it2
      if(iti.lt.0) then
         filename(3:3)='m'
         iti=-iti
      else
         filename(3:3)='e'
      endif
      write(filename(4:4),'(i1.1)')iti

      filename(5:5)='v'
      write(filename(6:8),'(i3.3)')nint(100*vd)
      filename(9:9)='r'
      write(filename(10:11),'(i2.2)')ifix(r(nr))
      filename(12:12)='P'
      write(filename(13:14),'(i2.2)')ifix(abs(Vprobe))

      filename(15:15)='L'
      if(debyelen.gt.1.e-10)then
         iti=nint(alog10(debyelen)-0.49)
         it2=nint(debyelen/10.**iti)
      else
         it2=0
         iti=0
      endif
      write(filename(16:16),'(i1.1)')it2
      if(iti.lt.0) then
         filename(17:17)='m'
         iti=-iti
      else
         filename(17:17)='e'
      endif
      write(filename(18:18),'(i1.1)')iti

      filename(19:19)='B'
      if(Bz.gt.1.e-10)then
         iti=nint(alog10(Bz)-0.49)-1
         it2=nint(Bz/10.**iti)
      else
         it2=0
         iti=0
      endif
      write(filename(20:21),'(i2.2)')it2
      if(iti.lt.0) then
         filename(22:22)='m'
         iti=-iti
      else
         filename(22:22)='e'
      endif
      write(filename(23:23),'(i1.1)')iti

      filename(24:27)='.orb'

      open(15,file=filename)
      write(15,*) 'Number of orbits'
      write(15,*) norbits
      do k=1,norbits
         write(15,*) k,'th orbit'
         write(15,*) iorbitlen(k)
         do i=1,iorbitlen(k)
            write(15,590) xorbit(i,k),yorbit(i,k),zorbit(i,k),
     $      vxorbit(i,k),vyorbit(i,k),vzorbit(i,k)
         enddo
      enddo
 590  format(6f9.4)
      close(15)
      end
c*********************************************************************
c Write a second file with the force data as a function of step
      subroutine outforce(filename,istepmax)
      character*(*) filename
      include 'piccom.f'
      real zmn(nstepmax,4,2)
c Apply normalization factors but don't change the zmom.
c Perhaps this extra storage is unnecessary.
      do i=1,istepmax
         do k=1,2
            zmn(i,enccharge,k)=zmom(i,enccharge,k)
            zmn(i,partz,k)=zmom(i,partz,k)/rhoinf
            zmn(i,epressz,k)=zmom(i,epressz,k)
            zmn(i,fieldz,k)=zmom(i,fieldz,k)*debyelen**2
         enddo
      enddo
      open(9,file=filename)
      write(9,*)istepmax
      write(9,*)'Step    Charge     E-field      Electrons',
     $        '      Ions   Total Force'
      write(9,'(i5,5f12.5)')((i,(zmn(i,j,k),j=1,4)
     $     ,zmn(i,partz,k)+zmn(i,fieldz,k)+zmn(i,epressz,k)
     $     ,k=1,2),i=1,istepmax)
      close(9)
      end
c**********************************************************************
c Write out the particle data.
      subroutine partwrt()
c Common data:
      include 'piccom.f'
      character*11 filename

      write(filename,'(''part'',i3.3,''.dat'')')myid
c Delete the file first to help with nfs problems.
      open(11,file=filename,status='unknown')
      close(11,status='delete')
c
      open(11,file=filename,status='unknown')
      write(11,*)npartmax,npart,nr,nth,ndim,np
      write(11,*)xp
      write(11,*)rhoinf,spotrein,averein
      write(*,*)'rhoinf,spotrein,averein',rhoinf,spotrein,averein
      close(11)
      end

c**********************************************************************
c Read in the particle data.
      subroutine partrd(success)
      logical success
c Common data:
      include 'piccom.f'
      character*11 filename

      write(filename,'(''part'',i3.3,''.dat'')')myid
      success=.false.
      open(11,file=filename,status='old',err=101)
      read(11,*,err=100,end=100)ipartmax,ipart,ir,ith,idim,ip
      if(ipartmax.eq.npartmax .and. ipart.eq.npart
     $     )then
c     $ .and. ir.eq.nr .and. ith.eq.nth )then
         write(*,*)'Using saved particle data.'
         read(11,*,err=100,end=100)xp
         read(11,*,err=100,end=100)rhoinf,spotrein,averein
      write(*,*)'rhoinf,spotrein,averein',rhoinf,spotrein,averein
         success=.true.
      else
         write(*,*)'Particle data mismatch',ipartmax,npartmax,ipart,
     $        npart,ir,nr,ith,nth
      endif
      close(11)
      return
 100  close(11)
      write(*,*) 'Error reading pardata.dat'
      return
 101  write(*,*) 'No particle file to read.'
      end
c**********************************************************************
c Get the average and slope over the rmesh range i1,i2.
      subroutine slopegen(phi,r,nr,i1,i2,slope,average)
      integer nr
      real phi(nr),r(nr)

c Assume r-mesh is linear
      rmom0=0.
      rmom1=0.
      rmid=(r(i2)+r(i1))/2.
      do i=i1,i2
         rmom0=rmom0+phi(i)
         rmom1=rmom1+(r(i)-rmid)*phi(i)
      enddo
      average=rmom0/(i2-i1+1)
c      rave=rmom1/(i2-i1+1)
      rlen=r(i2)-r(i1)
      slope=12.*(rmom1)/(rlen*rlen)/(i2-i1+1)
c      write(*,*)rmom0,rmom1,r(i1),r(i2),rmid,rlen
      end

c**********************************************************************
      subroutine outsums(dt,i)
c Common data:
      include 'piccom.f'
c Write out summed results.
      nrhere=NRUSED
      nthhere=NTHUSED
      npsihere=NPSIUSED
      nphere=np
c Combined files. Don't have to open.
c      open(10,file=filename)
      write(10,'(a,a)')
     $     '    dt       vd       Ti  steps  rmax',
     $     ' rhoinf debyelen Vp /nr,nth,np; sums'
      write(10,'(2f8.5,f8.4,i6,f8.3,f12.3,2f14.5)')
     $     dt,vd,Ti,i,r(nr),rhoinf,debyelen,vprobe
      write(10,*)nrhere,nthhere,npsihere
      write(10,*)'pDiag'
      write(10,*)(((pDiag(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      write(10,*)'vrDiag'
      write(10,*)(((vrDiag(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      write(10,*)'vtDiag'
      write(10,*)(((vtDiag(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      write(10,*)'vpDiag'
      write(10,*)(((vpDiag(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      write(10,*)'vr2Diag'
      write(10,*)(((vr2Diag(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      write(10,*)'vt2Diag'
      write(10,*)(((vt2Diag(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      write(10,*)'vp2Diag'
      write(10,*)(((vp2Diag(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      write(10,*)'vrtDiag'
      write(10,*)(((vrtDiag(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
      write(10,*)'vrpDiag'
      write(10,*)(((vrpDiag(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)
       write(10,*)'vtpDiag'
      write(10,*)(((vtpDiag(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
     $     ,npsihere)

c Don't save cvx, vy,vzsum
c      write(10,*)'vxsum'
c      write(10,*)(((vxsum(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
c     $     ,npsihere)
c      write(10,*)'vysum'
c      write(10,*)(((vysum(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
c     $     ,npsihere)
c      write(10,*)'vzsum'
c      write(10,*)(((vzsum(k1,k2,k3),k1=1,nrhere),k2=1,nthhere),k3=1
c     $     ,npsihere)

      write(10,*)'r[cc]'
      write(10,*)(rcc(k1),k1=1,nrhere)
      write(10,*)'volinv'
      write(10,*)(volinv(k1),k1=1,nrhere)
      write(10,*)'t[cc]'
      write(10,*)(tcc(k2),k2=1,nthhere)
      write(10,*) 'p[cc]'
      write(10,*)(pcc(k3),k3=1,npsihere)
      end




