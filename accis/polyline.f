C********************************************************************
      subroutine polyline(x,y,npts)
c Dashed line version.
      real x(npts),y(npts)
      integer npts
      integer i
      include 'plotcom.h'
c Dashed line code
      real nx,ny
      real vlen,dx,dy,cx,cy,plen,flen,dlen
      real wx2nx, wy2ny
      integer cud

c dashlen is the arc length in normalized units of the the ith line
c segment. dashdist is the starting fractional part of the ith arc.
c Segments alternate pen down, pen up.
      logical ldash
      real dashlen,dashdist
      integer MASKNO,dashmask,jmask
      parameter (MASKNO=4)
      dimension dashmask(MASKNO),dashlen(MASKNO)
      common/dashline/ldash,dashlen,dashdist,dashmask,jmask

      if(npts.le.0) return
      call vecw(x(1),y(1),0)
      do 3 i=2,npts
	 if(.not.ldash) then
	    call vecw(x(i),y(i),1)
	 else
c We shall bypass vecw and go straight to normal.
	    nx=wx2nx(x(i))
	    ny=wy2ny(y(i))
c Lengths of total vector:
	    cx=nx
	    cy=ny
	    dx=nx- crsrx
	    dy=ny- crsry
	    vlen=sqrt(DX*DX+DY*DY)
c Partial length remaining:
	    plen=vlen
	    if(vlen.eq.0)return
c Distance to end of segment
    1       dlen=(dashlen(jmask)-dashdist)
	    if(plen.gt.dlen)then
c Vector longer than this segment. Draw segment and iterate.
	       dashdist=0
	       plen=plen-dlen
	       flen=dlen/vlen
	       nx= crsrx+dx*flen
	       ny= crsry+dy*flen
	       cud=dashmask(jmask)
c	       call optvecn(nx,ny,cud)
	       call vecn(nx,ny,cud)
	       jmask=mod(jmask,MASKNO)+1
	       goto 1
	    else
c Vector ends before segment. Draw to end of vector and quit.
	       dashdist=plen+dashdist
	       nx=cx
	       ny=cy
	       cud=dashmask(jmask)
c	       call optvecn(nx,ny,cud)
	       call vecn(nx,ny,cud)
	    endif
	 endif
    3 continue
      end
C********************************************************************
      block data dashdata

      logical ldash
      real dashlen,dashdist
      integer MASKNO,dashmask,jmask
      parameter (MASKNO=4)
      dimension dashmask(MASKNO),dashlen(MASKNO)
      common/dashline/ldash,dashlen,dashdist,dashmask,jmask
      data dashmask/1,0,1,0/
      data ldash/.false./
      data dashlen/0.03,0.01,0.01,0.01/
      data dashdist/1.e-5/
      data jmask/1/
      end
C********************************************************************
      subroutine dashset(i)
c Set some line styles. 0 solid (dashing off).
      integer i

      logical ldash
      real dashlen,dashdist
      integer MASKNO,dashmask,jmask
      parameter (MASKNO=4)
      dimension dashmask(MASKNO),dashlen(MASKNO)
      common/dashline/ldash,dashlen,dashdist,dashmask,jmask
      save
      if(i.eq.0)then
	 ldash=.false.
	 return
      elseif(i.eq.1)then
c Long dashes.
	 dashlen(1)=.03
	 dashlen(2)=.03
	 dashlen(3)=.03
	 dashlen(4)=.03
      elseif(i.eq.2)then
c Short dashes.
	 dashlen(1)=.01
	 dashlen(2)=.01
	 dashlen(3)=.01
	 dashlen(4)=.01
      elseif(i.eq.3)then
c Long/Short.
	 dashlen(1)=.03
	 dashlen(2)=.01
	 dashlen(3)=.01
	 dashlen(4)=.01
      elseif(i.eq.4)then
c 'Dots'.
	 dashlen(1)=.001
	 dashlen(2)=.01
	 dashlen(3)=.001
	 dashlen(4)=.01
      elseif(i.eq.5)then
c 'Medium/dot'.
	 dashlen(1)=.02
	 dashlen(2)=.01
	 dashlen(3)=.001
	 dashlen(4)=.01
      elseif(i.eq.6)then
c 'Long Dashes short breaks'.
	 dashlen(1)=.03
	 dashlen(2)=.01
	 dashlen(3)=.03
	 dashlen(4)=.01
      endif
      dashdist=1.e-6
      jmask=1
      ldash=.true.
      return
      end
C********************************************************************
      subroutine polymark(x,y,nx,nmark)
      integer nx,nmark,i
      character*4 mark
      real x(*),y(*),xp,yp,wx2nx,wy2ny
      include 'plotcom.h'
      ipf=pfPS
c      if(nmark.eq.3)then
c         mark='!A3'//char(0)
c  That does not quite align well. Better not to mix thinkgs up.
      if(nmark.lt.10)then
         pfPS=0
	 mark=char(nmark+176)//char(0)
      elseif(nmark.eq.10)then
	 mark='+'//char(0)
      elseif(nmark.eq.11)then
	 mark='!AX'//char(0)
      elseif(nmark.eq.12)then
	 mark='!A*'//char(0)
      elseif(nmark.eq.13)then
	 mark='!A-'//char(0)
      elseif(nmark.eq.15)then
	 mark='!A'//char(48)//char(0)
      else
	 mark=char(nmark)//char(0)
      endif
      do 1 i=1,nx
	 xp=wx2nx(x(i))
	 yp=wy2ny(y(i))
      if(nmark.eq.14)then
         call actrid(xp,yp)
      else
	 call jdrwstr(xp,yp,mark,0.)
      endif
    1 continue
      pfPS=ipf
      end
C********************************************************************
c Plot error bars from y to y+err.
      subroutine polyerr(x,y,err,nx)
      integer nx,i
      real x(*),y(*),err(1)
      do 1 i=1,nx
	 call vecw(x(i),y(i),0)
	 call vecw(x(i),y(i)+err(i),1)
    1 continue
      end
c******************************************************************
      subroutine legendline(xg,yg,nsym,string)
      real xg,yg
      integer nsym
      character*(*) string
c Draw a legendline
c At fractional position in plot box: xg, yg 
c Use symbol nsym if positive and <=256. 
c If nsym<0 use both line and symbol. If nsym=0 put only line.
c If nsym>256 use nothing but put the string at the usual place.
      include 'plotcom.h'

      real xd(2),yd(2)
      xleglen=.07
      xn=wx2nx(wxmin)*(1.-xg)+wx2nx(wxmax)*xg
      yn=wy2ny(wymin)*(1.-yg)+wy2ny(wymax)*yg
      xd(1)=xn2xw(xn)
      xd(2)=xn2xw(xn+xleglen)
      yd(1)=yn2yw(yn)
      yd(2)=yn2yw(yn)
      if(nsym.eq.0)then
         call polyline(xd,yd,2)
      elseif(nsym.eq.257)then
         call vecn(xn+xleglen,yn,0)
      elseif(nsym.eq.258)then
         call vecn(xn,yn,0)
      elseif(nsym.gt.0)then
         call polymark(xd(2),yd(2),1,nsym)
      elseif(nsym.lt.0)then
         call polymark(xd,yd,2,abs(nsym))
         call polyline(xd,yd,2)
      endif
      call drcstr(string)
c      write(*,*)xg,yg,xn,yn,xd(1),yd(1),nsym,string
      end


C********************************************************************
      subroutine polydraw(x,y,nx,drawfn)
      integer nx,i
      real x(*),y(*),xp,yp,wx2nx,wy2ny
      external drawfn
      do 1 i=1,nx
	 xp=wx2nx(x(i))
	 yp=wy2ny(y(i))
         call drawfn(xp,yp)
    1 continue
      end
c********************************************************************
      subroutine accircle(x,y)
      include 'plotcom.h'
      real x,y
      integer nangle
      parameter (nangle=40)
      iud=0
      do i=1,nangle
         angle=2.*3.14159*(i-1)/(nangle-1)
         xp=x+0.5*chrswdth*cos(angle)
         yp=y+0.5*chrshght*sin(angle)
         call vecn(xp,yp,iud)
         iud=1
      enddo
      end
c********************************************************************
      subroutine acx(x,y)
      include 'plotcom.h'
      real x,y
         xp=x+0.5*chrswdth
         yp=y+0.5*chrshght
         call vecn(xp,yp,0)
         xp=x-0.5*chrswdth
         yp=y-0.5*chrshght
         call vecn(xp,yp,1)
         xp=x+0.5*chrswdth
         yp=y-0.5*chrshght
         call vecn(xp,yp,0)
         xp=x-0.5*chrswdth
         yp=y+0.5*chrshght
         call vecn(xp,yp,1)
      end
c********************************************************************
      subroutine acplus(x,y)
      include 'plotcom.h'
      real x,y
         xp=x+0.5*chrswdth
         yp=y
         call vecn(xp,yp,0)
         xp=x-0.5*chrswdth
         yp=y
         call vecn(xp,yp,1)
         xp=x
         yp=y-0.5*chrshght
         call vecn(xp,yp,0)
         xp=x
         yp=y+0.5*chrshght
         call vecn(xp,yp,1)
      end
c********************************************************************
      subroutine acast(x,y)
      include 'plotcom.h'
      real x,y
      wtemp=chrswdth
      htemp=chrshght
      chrswdth=wtemp*0.707
      chrshght=htemp*0.707
      call acx(x,y)
      chrswdth=wtemp
      chrshght=htemp
      call acplus(x,y)
      end

c********************************************************************
      subroutine actrid(x,y)
      include 'plotcom.h'
      real x,y
c Inverted triangle
      call vecn(x,y-.33*sqrt(3.)*chrshght,0)
      call vecn(x+0.5*chrswdth,y+.17*sqrt(3.)*chrshght,1)
      call vecn(x-0.5*chrswdth,y+.17*sqrt(3.)*chrshght,1)
      call vecn(x,y-.33*sqrt(3.)*chrshght,1)
c      call pathfill()
      end
c********************************************************************
c General mechanism: Call acgset first then pass acgen.
      subroutine acgen(x,y)
      include 'plotcom.h'
      real x,y
c acgen data:
      integer nacg,nacgmax
      parameter (nacgmax=100)
      real xacg(nacgmax),yacg(nacgmax)
      common /acgcom/nacg,xacg,yacg,iacfill

      iud=0
      do i=1,nacg
         call vecn(x+xacg(i)*chrswdth,y+yacg(i)*chrshght,iud)
         iud=1
      enddo
      if(iacfill.ne.0)call pathfill()
      end
c********************************************************************
c Set the general accis symbol to the data in x,y, length n, scale chrs.
c if ifill ne 0 then the symbol is filled.
      subroutine acgset(x,y,n,ifill)
      real x(*),y(*)
      integer n
c acgen data:
      integer nacg,nacgmax
      parameter (nacgmax=100)
      real xacg(nacgmax),yacg(nacgmax)
      common /acgcom/nacg,xacg,yacg,iacfill
      nacg=min(n,nacgmax)
      iacfill=ifill
      do i=1,nacg
         xacg(i)=x(i)
         yacg(i)=y(i)
      enddo
      end



