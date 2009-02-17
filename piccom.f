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


   
      integer npartmax,npart,nr,nth,npsi,ndim,np
c Number of particles: npartmax, radial and theta mesh size: nr, nth.
c Don't change anything else.
      parameter (npartmax=200000,np=1,ndim=6)
c Use of particle advance subcycling in inner regions for accuracy.
      logical lsubcycle
c Integrator type. True=old, False=new symplectic schemes
      logical verlet
c CIC definitions
      logical LCIC,collcic
      integer NRUSED,NTHUSED,NPSIUSED,NRFULL,NTHFULL,NPSIFULL
      parameter (LCIC=.true.)
      integer nrsize,nthsize,npsisize
c These correspond to nrfull, nthfull.and npsifull
      parameter (nrsize=126,nthsize=36,npsisize=26)
c Positions and velocities of particles (6-d phase-space).
      real xp(ndim,npartmax)
      real vzinit(npartmax)

      real dtprec(npartmax)
c Flag of particle slot status (e.g. in use or not)
      integer ipf(npartmax)
c The potential normalized to Te/e
      real phi(0:nrsize,0:nthsize,0:npsisize)
c The potential on axis (cos(theta)=+-1) before averaging
      real phiaxis(0:nrsize,2,0:npsisize)
c Charge density
      real rho(0:nrsize,0:nthsize,0:npsisize)
      real rhoDiag(0:nrsize,0:nthsize,0:npsisize)
c Injection complement. How many particles to reinject each step
      integer ninjcomp
c Highest occupied particle slot.
      integer iocprev

      real pi
      parameter (pi=3.1415927)
      real cerr,bdyfc,Ti,vd,cd,cB,Bz
      logical diags,lplot,ldist,linsulate,lfloat,lat0,lap0,localinj
      logical lfixedn
      integer myid,numprocs
      real rmtoz
      common /piccom/xp,npart,vzinit,dtprec,phi,rho,rhoDiag,cerr,bdyfc
     $     ,Ti,vd,cd,cB,diags,ninjcomp,lplot,ldist,linsulate,lfloat
     $     ,lat0,lap0 ,localinj,lfixedn,myid,numprocs,rmtoz,ipf,iocprev
     $     ,Bz,lsubcycle,verlet,collcic,phiaxis


c *******************************************************************
c Momenta of the distribution function

c The particle number
      real psum(1:nrsize-1,1:nthsize-1,1:npsisize-1)
c The sum of particle radial velocities
      real vrsum(1:nrsize-1,1:nthsize-1,1:npsisize-1)
c Sum of theta velocities
      real vtsum(1:nrsize-1,1:nthsize-1,1:npsisize-1)
c Sum of phi velocities
      real vpsum(1:nrsize-1,1:nthsize-1,1:npsisize-1)
c The sum of particle velocities squared
      real vr2sum(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vt2sum(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vp2sum(1:nrsize-1,1:nthsize-1,1:npsisize-1)
c Off diagonal terms
      real vrtsum(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vrpsum(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vtpsum(1:nrsize-1,1:nthsize-1,1:npsisize-1)

c The sum of particle xyz-velocities
      real vxsum(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vysum(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vzsum(1:nrsize-1,1:nthsize-1,1:npsisize-1)

c Diagnostic sums
      real pDiag(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vrDiag(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vtDiag(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vpDiag(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vr2Diag(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vt2Diag(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vp2Diag(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vrtDiag(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vrpDiag(1:nrsize-1,1:nthsize-1,1:npsisize-1)
      real vtpDiag(1:nrsize-1,1:nthsize-1,1:npsisize-1)

c     Total sum of particle xyz-velocities, i.e. total current curr(4)
c     is the particle sum
      real curr(4)

      common /momcom/psum,vrsum,vtsum,vpsum,vr2sum,vt2sum,vp2sum ,vrtsum
     $     ,vrpsum,vtpsum,vzsum,vxsum,vysum,curr,pDiag,vrDiag,vtDiag
     $     ,vpDiag,vr2Diag,vt2Diag,vp2Diag ,vrtDiag,vrpDiag,vtpDiag
c*********************************************************************
c Radius mesh
      real r(0:nrsize),rcc(0:nrsize)
c Theta angle mesh
      real th(0:nthsize),tcc(0:nthsize)
c Theta mesh radians
      real thang(0:nthsize)
c Poloidal (Psi) Mesh. Only cell center are useful
      real pcc(0:npsisize)
c Inverse of volume of radial shells.
      real volinv(0:nrsize)
c Precalculation functions
      integer nrpre,ntpre,nppre
      parameter (nrpre=4*nrsize,ntpre=4*nthsize,nppre=4*npsisize)
      integer irpre(nrpre),itpre(ntpre),ippre(nppre)
      real rfac,tfac,pfac
c Non-uniform handling quantities.
      real hr(0:nrsize+1),zeta(0:nrsize+1),zetahalf(0:nrsize+1),
     $     cminus(nrsize),cmid(nrsize),cplus(nrsize)
c Lower limit of averaging range. 0.6 by default
      real avelim
c Parallel or serial solving
      logical cgparallel
c Parallel bloc solver arguments
      integer idim1,idim2,idim3

      common /meshcom/r,rcc,th,tcc,thang,volinv,irpre,itpre,rfac,tfac,
     $     pcc,ippre,pfac, hr,zeta,zetahalf,cminus,cmid,cplus ,avelim
     $     ,nr,NRFULL,NRUSED,NPSIFULL,NPSIUSED,nth,npsi,NTHFULL,NTHUSED
     $     ,cgparallel,idim1,idim2,idim3
c********************************************************************
c Random interpolate data.
      integer nvel,nQth
      parameter (nvel=50)
      parameter (nQth=200)
      real Qcom(nQth) 
      real Gcom(nvel,nQth)
      real Vcom(nvel)
      real pu1(nvel),pu2(nvel)
      real Pc(nQth,nvel)
c New BC
      integer bcphi,bcr
      logical infdbl
c Reinjection flux as a function of cos(theta) (line) and chi (column,
c from 0 to 9)
      common /rancom/Gcom,Vcom,Qcom,pu1,pu2,Pc,infdbl,bcphi,bcr
c********************************************************************
c diagnostic data
      integer nvmax,nrein,nreintry,ninner,nstepmax
      parameter (nvmax=60,nstepmax=10001)
      real nvdiag(nvmax),nvdiagave(nvmax),vdiag(nvmax)
      real vrdiagin(nvmax),vtdiagin(nvmax)
      real vrange
      real diagrho(nrsize),diagphi(nrsize)
      real diagchi(0:nthsize)
      real phiout
      integer partz,fieldz,epressz,enccharge,lorentz
      parameter(enccharge=1,fieldz=2,epressz=3,partz=4,lorentz=5)
c Total particle flux to probe
      real fluxprobe(nstepmax)
c Total momentum flux to probe
      real zmomprobe,xmomprobe,ymomprobe
c Total energy collected
      real enerprobe
c Z-momentum flux across outer boundary.
      real zmout,xmout,ymout
c Combined zmom data: fields, electron pressure, ion momentum.
c For inner 1, outer 2. zmom also carries the probe charge
      real zmom(nstepmax,5,2),xmom(nstepmax,2:5,2),ymom(nstepmax,2:5,2)
c enertot is the reduced enerprobe for each time-step
      real enertot(nstepmax)
c Number of particles striking probe in theta/psi cell
      real nincellstep(nthsize,npsisize,0:nstepmax)
      real vrincellstep(nthsize,npsisize,0:nstepmax)
      real vr2incellstep(nthsize,npsisize,0:nstepmax)
      real nincell(nthsize,npsisize)
      real vrincell(nthsize,npsisize)
      real vr2incell(nthsize,npsisize)
c Ave flux
      real fincellave(nthsize,npsisize)
c Ave radial mom flux
      real vrincellave(nthsize,npsisize)
c Ave radial vr2 flux
      real vr2incellave(nthsize,npsisize)
c Number of particles reinjected per theta cell.
c      integer noutrein(nth),ivoutrein(nth)
c Sum and average of potentials at which particles were reinjected.
      real spotrein,averein
c Flux of reinjection.
      real fluxrein
c Density at infinity.
      real rhoinf
c Number of trapped reinjections
      integer ntrapre
c Coefficient of density deficit, for external solution
      real adeficit
c Cell in which to accumulate distribution functions
      integer ircell,itcell
      common /diagcom/nvdiag,nvdiagave,vdiag,vrange,diagrho,diagphi,
     $     diagchi,phiout,nrein,nreintry,ninner,fluxprobe,nincellstep
     $     ,vrincellstep,vr2incellstep,nincell,vrincell,vr2incell,rhoinf
     $     ,vrdiagin,vtdiagin, spotrein,averein ,fluxrein ,ntrapre
     $     ,adeficit, ircell,itcell ,zmout,xmout,ymout ,zmomprobe
     $     ,ymomprobe,xmomprobe,fincellave ,vrincellave,vr2incellave
     $     ,zmom,xmom,ymom ,enerprobe ,enertot
c*********************************************************************
c Poisson coefficients for iterative solution, etc.

      real debyelen,vprobe,Ezext
      real apc(0:nrsize),bpc(0:nrsize)
      real cpc(0:nrsize,0:nthsize),dpc(0:nrsize,0:nthsize)
      real epc(0:nrsize,0:nthsize)
      real fpc(0:nrsize,0:nthsize)
      real gpc(0:nthsize,0:npsisize,1:5)
      common /poisson/debyelen,vprobe,Ezext,apc,bpc,cpc,dpc,fpc,epc,gpc
c*********************************************************************
c Smoothing steps
      integer nstepsave,nsamax,diagsamp
      logical samp
      common /stepave/nstepsave,nsamax,diagsamp,samp
c*********************************************************************
c Orbit plotting storage for tracking the first norbits orbits.
      integer nobsmax,norbits
      parameter (nobsmax=100)
      real xorbit(nstepmax,nobsmax),yorbit(nstepmax,nobsmax),
     $     zorbit(nstepmax,nobsmax)
      real vxorbit(nstepmax,nobsmax),vyorbit(nstepmax,nobsmax),
     $     vzorbit(nstepmax,nobsmax),rorbit(nstepmax,nobsmax)
      integer iorbitlen(nobsmax)
      common /orbits/norbits,iorbitlen,xorbit,yorbit,zorbit,rorbit,
     $     vxorbit,vyorbit,vzorbit

c*********************************************************************
c Data necessary for the orbit tracking
      logical orbinit
      integer maxsteps,trackinit
      common /orbtrack/orbinit,maxsteps,trackinit
