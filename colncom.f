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

c Collision data common.
c Drift velocity of neutral distribution
      real vneutral
      real vneut(3)
c Cos theta and psi of neutral drift velocity
      real cnd, psind
c Temperature of neutral distribution
      real Tneutral
c Required acceleration to give rise to specified ion average drift
c   (assuming no magnetic field)
      real Eneutral
      real Eneut(3)
c Number of collisions this time step
      integer NCneutral
      common /colcom/vneutral,Tneutral,Eneutral,NCneutral
     $  ,vneut,cnd,psind,Eneut
