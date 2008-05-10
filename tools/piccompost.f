c Version 3.0 SCEPTIC3D Jan 2008

      integer npartmax,npart,nr,nth,ndim,np
      parameter (npartmax=200000,nr=121,nth=81,npsi=80,np=1,ndim=6)
c NGP definitions
      logical LCIC
      integer NRUSED,NTHUSED,NPSIUSED,NRFULL,NTHFULL,NPSIFULL
      parameter (LCIC=.false.)
      parameter (NRUSED=nr-1,NTHUSED=nth-1,NPSIUSED=npsi,NRFULL=nr
     $     ,NTHFULL=nth,NPSIFULL=npsi)
c Positions and velocities of particles (6-d phase-space).
      real xp(ndim,npartmax)
c The particle number
      real psum(0:NRFULL,0:NTHFULL,0:NPSIFULL)
c The sum of particle radial velocities
      real vrsum(0:NRFULL,0:NTHFULL,0:NPSIFULL)
c Sum of theta velocities
      real vtsum(0:NRFULL,0:NTHFULL,0:NPSIFULL)
c Sum of psi velocities
      real vpsum(0:NRFULL,0:NTHFULL,0:NPSIFULL)
c The sum of particle velocities squared
      real v2sum(0:NRFULL,0:NTHFULL,0:NPSIFULL)
c The sum of radial particle velocities squared
      real vr2sum(0:NRFULL,0:NTHFULL,0:NPSIFULL)
c The sum of non-radial particle velocities squared
      real vtp2sum(0:NRFULL,0:NTHFULL,0:NPSIFULL)
c The sum of particle z-velocities
      real vxsum(0:NRFULL,0:NTHFULL,0:NPSIFULL)
c The sum of particle z-velocities
      real vysum(0:NRFULL,0:NTHFULL,0:NPSIFULL)
c The sum of particle z-velocities
      real vzsum(0:NRFULL,0:NTHFULL,0:NPSIFULL)
c The potential normalized to Te/e
      real phi(0:NRFULL,0:NTHFULL,0:NPSIFULL)
c Charge density
      real rho(0:NRFULL,0:NTHFULL,0:NPSIFULL)
c Mag field
      real Bz

      real pi
      parameter (pi=3.1415927)
      real cerr,bdyfc,Ti,vd,cd,damplen
      logical diags,lplot
      integer myid,numprocs
      common /piccom/xp,npart, psum,vrsum,vtsum,vpsum,v2sum,vr2sum
     $     ,vtp2sum,vxsum,vysum,vzsum, phi,rho,cerr,bdyfc,Ti,vd,cd,diags
     $     ,lplot ,myid,numprocs,Bz,damplen

c Radius mesh
      real r(0:NRFULL),rcc(0:NRFULL)
c Theta angle mesh
      real th(0:NTHFULL),tcc(0:NTHFULL)
c Psi angle mesh (In radians)
      real pcc(0:NPSIFULL)
c Theta mesh radians
      real thang(0:NTHFULL)
c Inverse of volume of radial shells.
      real volinv(0:NRFULL)
c Precalculation functions
      integer nrpre,ntpre,nppre
      parameter (nrpre=1000,ntpre=100,nppre=100)
      integer irpre(nrpre),itpre(ntpre),ippre(nppre)
      real rfac,tfac
c Non-uniform handling quantities.
      real hr(0:nr+1),zeta(0:nr+1),zetahalf(0:nr+1),
     $     cminus(nr),cmid(nr),cplus(nr)
c Lower limit of averaging range. 0.6 by default
      real avelim

      common /meshcom/r,rcc,th,tcc,pcc,thang,volinv,irpre,itpre,ippre
     $     ,rfac,tfac,hr,zeta,zetahalf,cminus,cmid,cplus,avelim

      integer nvel
      parameter (nvel=50)
      real Qcom(nth)
      real Gcom(nvel,nth)
      real Vcom(nvel)
      common /rancom/Gcom,Vcom,Qcom

      integer nvmax,nrein,ninner,nstepmax
      parameter (nvmax=60,nstepmax=10000)
      real nvdiag(nvmax),nvdiagave(nvmax),vdiag(nvmax)
      real vrange
      real diagrho(nr),diagphi(nr)
c Total particle flux to probe
      real fluxprobe(nstepmax)
c Number of particles striking probe in theta cell
      real nincell(nth,npsi)

c Number of particles reinjected per theta cell.
      integer noutrein(nth),ivoutrein(nth)
      real rhoinf

c Angular distribution with respect to the drift velocity
      integer angsize
c It is better that angsize be odd for the interpolation
      parameter (angsize=21)
c      real canghere
      real cang(angsize)
      real fluxang(angsize)
      real aweight(angsize)
      integer napre
      parameter(napre=4*angsize)
      integer iapre(napre)
      real afac

      common /diagcom/nvdiag,nvdiagave,vdiag,vrange,diagrho,diagphi,
     $     nrein,ninner,fluxprobe,nincell,rhoinf,noutrein ,ivoutrein
     $     ,cang,fluxang,aweight,iapre,afac



