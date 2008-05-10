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

c******************************************************************
      subroutine shielding3D_par(dt,n1,cg_comm,myid2)

      include 'piccom.f'
c cg_comm is the subset of MPI_COMM_WORLD communicator used for the
c bloc conjugate gradient
      integer cg_comm,myid2
      real dt,dconverge
      integer maxits
      integer n1
      integer kk1,kk2


      maxits=2*(nrused*nthused*npsiused)**0.333
      dconverge=1.e-5


c     In the serial solver, one can treat the potential with size
c     (nrsize+1)*(ntsize+1)*(npsize+1) and the density as
c     (nrsize-1)*(ntsize-1)*(npsize-1). With the paralle version,
c     because of the indexes, it is much simpler to take all the arrays
c     with size (nrsize+1)*(ntsize+1)*(npsize+1)
      call fcalc3Dpar(nrsize+1,nthsize+1,npsisize+1,n1+1,nthused+2
     $     ,npsiused+2,maxits,dconverge,iter,cg_comm,myid2)
      

c Output the number of iterations
      if(myid2.eq.0)  then
         write(*,'('':'',i3,$)')iter

c We must set the potential of the shadow theta-cells to the
c physical value of the potential at psi+pi. The zero-derivative
c condition is not valid anymore in the 3D case

         do k=1,npsiused
            kk1=mod(k+3*npsi/2-1,npsi)+1
            kk2=mod(k+(3*npsi+1)/2-1,npsi)+1
            do i=1,nrused
               phi(i,0,k)=0.5*(phi(i,2,kk1)+phi(i,2,kk2))
               phi(i,nthused+1,k)=
     $              0.5*(phi(i,nthused-1,kk1)+phi(i,nthused-1,kk2))
            enddo
         enddo
         
c Set the theta shadow cells to their proper value to ensure periodicity
         do j=0,nthused+1
            do i=1,nrused
               phi(i,j,npsiused+1)=phi(i,j,1)
               phi(i,j,0)=phi(i,j,npsiused)
            enddo
         enddo
         
      endif

      end
