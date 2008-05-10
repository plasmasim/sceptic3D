c New version
      program fontshow
c Plot aligned sets of the three installed fonts.
      character*256 entered
      character*129 str1
      real js,px,py
      integer i,j
c
c      write(*,*)' Plot to file? (0:no,1:hp,2:ps,3:eps)'
c      read(*,*)i
c      call pfset(i)
      call pltinit(0.,1.,0.,1.)
      call charsize(.1,.1)
      do 8 k=1,1
	 if(k.eq.3)then
	    noff=1
	 else
	    noff=0
	 endif
	 do 9 j=1,3
	    js=0.
	    do 10 i=1,6
	       py=.7-(4*(k-1)+j)*0.043*4.
	       px=0.02+.03*i*4.
               if(char(31+32*(j-1)+i).eq.'\'
     $              .or.char(31+32*(j-1)+i).eq.'!')then 
                  str1="!"//char(63+k)//'!'//
     $                 char(31+32*(j-1)+i)//"!@"//char(0)
               else
                  str1="!"//char(63+k)//char(31+32*(j-1)+i)
     $                 //"!@"//char(0)
               endif
	       call drwstr(px,py,str1)
   10	    continue
    9	 continue
    8 continue
      call pltend()
      stop
      end

