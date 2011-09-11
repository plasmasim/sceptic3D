c Collision data common.
c Drift velocity of neutral distribution
      real vneutral
      real vneut(3)
c Cos theta and psi of neutral drift velocity
      real cnd, psind
c Temperature of neutral distribution
      real Tneutral
c Required acceleration to give rise to specified ion average drift
      real Eneutral
      real Eneut(3)
c Relative drift of ions compared to neutrals
      real reldrift(3)
c Number of collisions this time step
      integer NCneutral
      common /colcom/vneutral,Tneutral,Eneutral,NCneutral
     $  ,vneut,cnd,psind,Eneut,reldrift
