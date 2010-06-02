      program datwrt
c Write a block data program for filling a big character array.
      character*1 ch
c Characters. Block from plotih.
c Link with getfont to be able to get the font data.
      include 'plotcom.h'
      integer ip,i,j,ii,ip5
      character*72 fline
      write(*,*)'      block data fonts'
      write(*,*)'      character*19000 ca'
      write(*,*)'      integer*2 a(384)'
      write(*,*)'      include ''plotcom.h'''
      write(*,*)'      character*2 crlf'
      write(*,*)'      parameter (crlf=char(13)//char(10))'
      write(*,*)'      equivalence (ca,chrsfont)'
      write(*,*)'      equivalence (a,chrsaddr)'

      call getfont('cgi.fnt     ')
c write the index data
      do 21 i=1,64
	 ii=6*(i-1)+1
	 ip5=ii+5
	 write(*,'(1x,''      data (a(j),j='',i3,'',''i3'')/'',
     $	      i5,'','',i5,'','',i5,'','',i5,'','',i5,'','',i5,''/'')')
     $	      ii,ip5, (chrsaddr(j),j=ii,ip5)
   21 continue

c Write the font data.
      fline(1:14)='      data ca('
      fline(20:20)=':'
      fline(26:28)=')/'//char(39)
      j=28
      ip=1
      do 1 i=1,BUFFER
	 ch=chrsfont(i)
	 if(ch.eq.char(10))then
c	    if(j.gt.28) call lnwrt(fline,j,ip,i)
c	    write(*,'(1x,
c     $        ''      data ca('',i5.5,'':'',i5.5,'')/10/'')')i,i
c	    ip=i+1
	 elseif(ch.eq.char(13))then
	    if(j.gt.28) call lnwrt(fline,j,ip,i)
c	    write(*,'(1x,
c     $        ''      data ca('',i5.5,'':'',i5.5,'')/13/'')')i,i
	    write(*,'(1x,
     $        ''      data ca('',i5.5,'':'',i5.5,'')/crlf/'')')i,i+1
	    ip=i+2
	 elseif(ch.eq.char(26))then
	    if(j.gt.28) call lnwrt(fline,j,ip,i)
	    goto 2
	 else
c Changed to 60 to try to fix f2c. Didn't work.
	    if(j.ge.70) call lnwrt(fline,j,ip,i)
	    call lnadd(fline,j,ch)
	 endif
    1 continue
      if(j.gt.28)call lnwrt(fline,j,ip,i)
    2 write(*,*)'      end'
      end

c*******************************************************************
      subroutine lnadd(fline,j,ch)
      character ch
      character*72 fline
      integer j
      j=j+1
      fline(j:j)=ch
      end
c*******************************************************************
      subroutine lnwrt(fline,j,ip,i)
      character*72 fline
      integer ip,i,j
      write(fline(15:19),'(i5.5)')ip
      write(fline(21:25),'(i5.5)')i-1
      fline(j+1:j+2)=char(39)//'/'
      write(*,*)fline(1:j+2)
      j=28
      ip=i
      end
