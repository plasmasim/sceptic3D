c Worksurface common:
      integer nrows,ncolumns,nframe,multype
      real naxmin,naxmax,naymin,naymax,naxpt,naypt
      common/wksrfc/naxmin,naxmax,naymin,naymax,naxpt,naypt,
     $	   nrows,ncolumns,nframe,multype

c Worlds and labels common:
      real xticlen,yticlen,xticoff,yticoff
      integer nxlabw,nxlabp,nylabw,nylabp,ticnum
      logical lxlog,lylog
      integer nxpow,nypow
      real xpow,ypow,wxmin,wxmax,wymin,wymax,w2nx,w2ny,n2sy
      logical ltlog
      common/worlds/
     $	xticlen,yticlen,xticoff,yticoff,
     $	nxlabw,nxlabp,nylabw,nylabp,ticnum,
     $	lxlog,lylog,nxpow,nypow,xpow,ypow,
     $  wxmin,wxmax,wymin,wymax,w2nx,w2ny,n2sy,ltlog
c
c  Cursor position
      real crsrx,crsry
      integer*2 updown
      common/crsrxy/crsrx,crsry,updown
c
c  Truncation monitoring
      real trcxmi,trcxma,trcymi,trcyma
      common/trncat/ trcxmi,trcxma,trcymi,trcyma
c
c  Characters
      integer BUFFER,NOCHARS
      parameter (BUFFER=19000,NOCHARS=384)
      real chrscos,chrssin,chrswdth,chrshght,chrsslnt
      integer*2 chrsaddr
      character*1 chrsfont
      common/chrcter/chrscos,chrssin,chrswdth,chrshght,chrsslnt
     &  ,chrsaddr(NOCHARS)
      common/chrfnt/ chrsfont(BUFFER)
c
c  Plot-to-file control
      integer pfsw
      integer pfilno,pfnextsw,pfPS,psini
      common/pltfil/pfsw,pfilno,pfnextsw,pfPS,psini
c
c  Screen Parameters
      integer*2 scrxpix,scrypix,ncolor,vmode
      real yoverx
      common/screenp/scrxpix,scrypix,ncolor,vmode,yoverx

