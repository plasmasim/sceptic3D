c Tektronix 4014 driver.
C********************************************************************
      blockdata scrndat
      include 'plotcom.h'
c      data scrxpix,scrypix,ncolor,vmode/1024,779,15,4010/
      data scrxpix,scrypix,ncolor,vmode/4096,3120,15,4014/
      end
C********************************************************************
c Switch to graphics mode Tek 4010/14.
      subroutine svga(scrxpix,scrypix,vmode,ncolor)
      integer scrxpix,scrypix,vmode,ncolor
c Enter Tek mode. Modified for Xterm.
      write(*,'(1x,a)')char(27)//'[?38h'
      write(*,*)'                                 '
c      write(*,'(1x,a)')char(27)//char(12)
c      write(*,*)'                                 '
c Clear screen twice seems to do the trick for Kermit.
      write(*,'(1x,a)')char(27)//char(12)
      write(*,*)'                                 '
c Extra clear screen for non Kermit.
      write(*,'(1x,a)')char(12)
      return
      end
C********************************************************************
      subroutine txtmode
      character resp
c Write padding, so we don't lose characters on switch back.
      read(*,'(a1)')resp
      write(*,*)char(24),'                                         '
      write(*,*)char(24),'                                         '
c Fix for xterm. Did not work.
c      write(*,'(1x,a)')char(27)//'[?38l'
      end
c*********************************************************************
c 4014 vector driver.
      subroutine vec(x,y,iud)
      integer x,y,iud
      include 'plotcom.h'
      integer oylst,oyhi,oxlow,oxhi
      integer ylow,yhi,xlow,xhi,i,xlst,ylst
      integer hlab,ylab,xlab,istart
      parameter (hlab=32,ylab=96,xlab=64)
      character*6 start
c      parameter (start=char(29)//char(hlab)//char(ylab)//char(ylab)
c     $  //char(hlab)//char(xlab))
      parameter (start=' '''' @')
      character*80 outchr
      save
      data i,istart/6,1/
      data outchr(1:6)/start/
c Crunch on to the screen.
      if(x.lt.0)x=0
      if(x.gt.scrxpix)x=scrxpix
      if(y.lt.0)y=0
      if(y.gt.scrypix)y=scrypix
c Separate the coordinate 'nibbles'.
      xhi=x/128
      xlst=x-128*xhi
      xlow=xlst/4
      xlst=xlst-(xlow)*4
      xhi=xhi+hlab
      xlow=xlow+xlab
      ylst=3119-y
      yhi=ylst/128
      ylst=ylst-128*yhi
      ylow=ylst/4
      ylst=ylst-ylow*4
      yhi=yhi+hlab
      ylow=ylow+ylab
      ylst=ylab+4*ylst+xlst
      if(iud.gt.0)then
c Continuing vector. Send only the necessary parts.
	 if(yhi.ne.oyhi)then
	    i=i+1
	    outchr(i:i)=char(yhi)
	 endif
	 if(ylst.ne.oylst)then
	    i=i+1
	    outchr(i:i)=char(ylst)
	 endif
	 i=i+1
	 outchr(i:i)=char(ylow)
	 if(xhi.ne.oxhi)then
	    i=i+1
	    outchr(i:i)=char(xhi)
	 endif
	 i=i+1
	 outchr(i:i)=char(xlow)
         istart=0
      else
c Start vector.
         istart=1
      endif
      if(i.ge.74.or.istart.eq.1)then
c Finish draw and start again, if we are longer than a line or starting.
         write(*,999)outchr(1:i)
  999	 format(1x,a)
	 outchr(1:6)=char(29)//char(yhi)//char(ylst)
     $	      //char(ylow)//char(xhi)//char(xlow)
	 i=6
      endif
c Update the vector.
      oxhi=xhi
      oxlow=xlow
      oyhi=yhi
      oylst=ylst
      return
      end
c*****************************************************************
      subroutine scolor(li)
      include 'plotcom.h'
      integer li,mask
c Scale the color so that (0..15) -> (0...ncpre-1) if ncpre >16.
c But a nonzero color always makes some mark because mask>0. check this!
      if(li.eq.0)then
         mask=0
      elseif(li.eq.8)then
	 mask=7
      else
         mask=mod(li,8)
      endif
c Terminate the vector so as to flush buffer
      call vecn(crsrx,crsry,0)
c Set the color.
c ANSI sequences for kermit, assumed if we are a 4014. Flush buffer first.
c      write(*,*)
c      write(*,'(1x,a1,''['',i2,''m'')')char(27),30+mask
      return
      end
c*********************************************************************



