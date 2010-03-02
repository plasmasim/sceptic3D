c*********************************************************************
c Calculate the dot product of vectors a and b with dimension n
      real function dot(a, b, n)
      implicit none
      real a(*), b(*)
      integer i,n
      dot = 0.
      do i=1,n
         dot = dot + a(i)*b(i)
      enddo
      return
      end

c*********************************************************************
c Calculate the cross product c of 3D vectors a and b
      subroutine cross(a, b, c)
      implicit none
      real a(3), b(3), c(3), tmp(3)
      integer i
      tmp(1) = a(2)*b(3) - a(3)*b(2)
      tmp(2) = -( a(1)*b(3) - a(3)*b(1) )
      tmp(3) = a(1)*b(2) - a(2)*b(1)
c     Using tmp in case a or b is same array as c
      c(1) = tmp(1)
      c(2) = tmp(2)
      c(3) = tmp(3)
      end

