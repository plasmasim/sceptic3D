
c******************************************************************
      subroutine shielding3D_par(dt,n1,cg_comm,myid2)

      include 'piccom.f'
      include 'errcom.f'
c cg_comm is the subset of MPI_COMM_WORLD communicator used for the
c bloc conjugate gradient
      integer cg_comm,myid2
      real dt,dconverge
      integer maxits
      integer n1
      integer kk1,kk2


      maxits=2*(nrused*nthused*npsiused)**0.333
      dconverge=1.e-5

c Set the potential on axis to what the solver found at the previous timestep
      do i=2,n1
         do k=1,npsiused
            phi(i,1,k)=phiaxis(i,1,k)
            phi(i,nthused,k)=phiaxis(i,2,k)
         enddo
      enddo

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

c     We set the potential on the inner shadow cell by second order
c     extrapolation from the potential at i=1,2,3
      do j=1,nthused
         do k=1,npsiused
c            phi(0,j,k)=2*phi(1,j,k)-phi(2,j,k) Only 1st order
            phi(0,j,k)=2.5*phi(1,j,k)-2*phi(2,j,k)+0.5*phi(3,j,k)
         enddo
      enddo


c     We must average the potential of each psi-cell at theta=0 or
c     theta=pi. The boundary conditions assume the cell center is on
c     axis, while what is indeed calculated is the potential at the
c     center of the first/last theta-cells.

      do i=0,nrused
         psiave1=0.
         psiave2=0.
         do k=1,npsiused
            psiave1=psiave1+phi(i,1,k)
            psiave2=psiave2+phi(i,nthused,k)
         enddo
         psiave1=psiave1/npsiused
         psiave2=psiave2/npsiused
         do k=1,npsiused
            phiaxis(i,1,k)=phi(i,1,k)
            phiaxis(i,2,k)=phi(i,nthused,k)
            phi(i,1,k)=psiave1
            phi(i,nthused,k)=psiave2
         enddo
      enddo


c We must set the potential of the shadow theta-cells to the
c physical value of the potential at psi+pi. The zero-derivative
c condition is not valid anymore in the 3D case

         do k=1,npsiused
            kk1=mod(k+3*npsi/2-1,npsi)+1
            kk2=mod(k+(3*npsi+1)/2-1,npsi)+1
            do i=0,nrused
               phi(i,0,k)=0.5*(phi(i,2,kk1)+phi(i,2,kk2))
               phi(i,nthused+1,k)=
     $              0.5*(phi(i,nthused-1,kk1)+phi(i,nthused-1,kk2))
            enddo
         enddo
         
c Set the psi shadow cells to their proper value to ensure periodicity
         do j=0,nthused+1
            do i=0,nrused
               phi(i,j,npsiused+1)=phi(i,j,1)
               phi(i,j,0)=phi(i,j,npsiused)
            enddo
         enddo
         
      endif

      end
