
c Collision data common.
c Drift velocity of neutral distribution
      real vneutral
c Temperature of neutral distribution
      real Tneutral
c Required acceleration to give rise to specified ion average drift.
      real Eneutral
c Number of collisions this time step
      integer NCneutral
      common /colcom/vneutral,Tneutral,Eneutral,NCneutral
