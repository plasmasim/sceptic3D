
      subroutine cgparinit(myid2,cg_comm)

c     Creates a new communicator, cg_comm, who contains a subset of
c     MPI_COMM_WORLD, because we don't need all the processes for the
c     bloc CG (Conjugate gradient)

c     Not to be confused with mpicommcart, which is a reorganisation of
c     cg_comm

c     myid2 is the process id in the new communicator. If myid2 is negative,
c     the process does not belong to the communicator

      integer myid2,nproccg

c     Make sure the size is large enough when using many
c     processors. With Loki, can't be larger than 512.
      integer members(0:511)

      integer group_world,cg_group,cg_comm


      include 'piccom.f'
      include 'errcom.f'
      include 'mpif.h'

      idim1=nint((numprocs*NRUSED**2/ (NTHUSED*NPSIUSED*1.))**(1./3))
      idim2=nint((numprocs*NTHUSED**2/ (NRUSED*NPSIUSED*1.))**(1./3))
      idim3=(numprocs/(idim1*idim2)*1.)

      
      if(idim1.le.1) then
         idim1=2
         if(idim2.ge.4) then
            idim2=idim2/2
         else
            idim3=idim3/2
         endif
      endif
      if(idim2.le.1) then
         idim2=2
         if(idim3.ge.4) then
            idim3=idim3/2
         else
            idim1=idim1/2
         endif
      endif
      if(idim3.le.1) then
         idim3=2
         if(idim2.ge.4) then
            idim2=idim2/2
         else
            idim1=idim1/2
         endif
      endif
      
      nproccg=idim1*idim2*idim3
      if(nproccg.gt.numprocs) idim1=idim1/2
      nproccg=idim1*idim2*idim3

c     Quick fix for only 8 processors 
      if(numprocs.eq.8) then
         nproccg=8
         idim1=2
         idim2=2
         idim3=2
      endif

      do i=0,nproccg-1
         members(i)=i
      enddo

      call MPI_COMM_GROUP(MPI_COMM_WORLD,group_world,ierr)
      call MPI_GROUP_INCL(group_world,nproccg,members,cg_group,ierr)
      call MPI_GROUP_RANK(cg_group,myid2,ierr)
      call MPI_COMM_CREATE(MPI_COMM_WORLD,cg_group,cg_comm,ierr)
      if (myid2.ge.0) then
         call MPI_COMM_RANK(cg_comm,myid2,ierr)
      endif


      if (myid2.eq.0) then
c         write(*,*) "cgp : ",idim1,idim2,nproccg,numprocs
         write(*,510) idim1,idim2,idim3,nproccg
 510     format('BlocX=',i2,'  BlocY=',i2,'  BlocZ=',i2,'  CgProc=',i2)
      endif

      end

c ------------------------------------------------------------------------

c     Solve an elliptical problem in 3d by conjugate gradient
c     A(u)=q(x,y,z) where A is a second order symmetric linear operator
c     represented by a difference stencil of specified coefficients, and
c     q is the "charge density".
c
      subroutine cg3dmpi(cg_comm,Li,Lj,Lk,ni,nj,nk,bcphi,u,q
     $     ,ictl,ierr,mpiid,idim1,idim2,idim3,apc,bpc,cpc,dpc,epc,fpc
     $     ,gpc,b,x,p,res,z,pp,resr,zz,lbcg)

      integer cg_comm,mpiid
c     The number of dimensions, 3 here.
      integer nd
      parameter (nd=3,nd2=nd*2)
c     Li leading dimension of i 
c     ni active dimension of i
      integer Li,Lj,Lk,ni,nj,nk
c     In this mpi routine we must use linear addressing for cij,u,q
c     so that the pointers can be used for each block.
      real apc(*),bpc(*),cpc(*),dpc(*),epc(*),fpc(*)
      real gpc(0:Lj-1,0:Lk-1,1:5)
c     u  potential to be solved for (initialized on entry).
c        its boundaries are at 1,ni; 1,nj; 
      real u(*)
c     q  "charge density" input
      real q(*)
c     mpiid Returns my MPI process number for this mpi version.
      integer ictl,ierr,kc
c     idim1,2,3   The number of blocks in dimensions 1, 2 and 3
      integer idim1,idim2,idim3
     

      common /cg3dctl/icg_mi,cg_eps,cg_del,icg_k

c      Other things that we might want control over include the maximum
c     number of iterations and the convergence size.
      real delta,deltamax
      logical lconverged


c     icg_prec is the number of iterations at the previous
c     timestep. Necessary for an estimation of when to start checking
c     for convergence.
      real icg_prec

c     BC
      integer bcphi
c     lflag, decide if we call bbdy for the first time or not
      logical lflag

      
      include 'piccomcg.f'
      include 'errcom.f'

c The origin of blocks structure may be considered
c      integer iorig(idim1+1,idim2+1)
c we declare it as a 1-d array. This is the first time here:
c It must be of dimension greater than the number of processes (blocks)
c      parameter (norigmax=1000)
c      integer iorig(norigmax)
c      parameter (ndims=3)
c bbdydecl declares most things for bbdy, using parameter ndims.
c      include 'bbdydecl.f'
c MPI parameters

      include 'mpif.h'
      logical lperiod(ndims)
      data lperiod(1)/.false./lperiod(2)/.false./lperiod(3)/.true./

c     Index for array storage
      integer index
      real bknum,bknumR,bkden,akden,akdenR
c     Temporary arrays for the cg solver
      real b(*),x(*),p(*),res(*),z(*),pp(*),resr(*),zz(*)
      
c     Flag to use biconjugate gradient method (not minimum residual)
      logical lbcg

 
c-------------------------------------------------------------------

      save lflag

c Required iterations at the previous step
      icg_prec=icg_k

c Initialize denominators to avoid compilation warnings
      akden=0.
      akdenR=0.
      bknum=0.
      bkden=1.
      bknumR=0.

      idims(1)=idim1
      idims(2)=idim2
      idims(3)=idim3
      iuds(1)=ni
      iuds(2)=nj
      iuds(3)=nk
      ifull(1)=Li
      ifull(2)=Lj
      ifull(3)=nk
      if((idims(1)+1)*(idims(2)+1)*idims(3).gt.norigmax)then
         write(*,*)'Too many processes',idims(1),' x ',idims(2),' x '
     $        ,idims(3),' for norigmax=',norigmax
         stop
      endif


c Define mpi block structure.
      call bbdydefine(ndims,idims,ifull,iuds,iorig,iLs)
c End of control functions.
c-------------------------------------------------------------------
      ierr=0
      lconverged=.false.

      

c Do block boundary communications, returns block info icoords...myid.
      call bbdy(cg_comm,iLs,iuds,u,icg_k,iorig,ndims,idims,lperiod,
     $     icoords,iLcoords,myside,myorig,myorig1,myorig2,myorig3,
     $     icommcart,mycartid,mpiid,lflag,out,inn)



c For debugging, if lAdebug flag set, only multiply by A; do not solve
      if (lAdebug) then
c        Communicate x on boundaries (needed for psi periodicity)
         call bbdy(cg_comm,iLs,iuds,x,icg_k,iorig,ndims,idims,lperiod,
     $     icoords,iLcoords,myside,myorig,myorig1,myorig2,myorig3,
     $     icommcart,mycartid,mpiid,lflag,out,inn)
c        Do matrix multiplication
         call atimesmpi(myside(1),myside(2),myside(3),Li,Lj,Lk,x(myorig)
     $     ,res(myorig),u(myorig),apc(myorig1),bpc(myorig1) ,cpc(myorig1
     $     +Li*myorig2) ,dpc(myorig1+Li*myorig2),epc(myorig1+Li*myorig2)
     $     ,fpc(myorig1+Li*myorig2)
     $     ,gpc(myorig2,myorig3,1), out, lAtranspose)
c        Do the final mpi_gather
         kc=-1
         call bbdy(cg_comm,iLs,iuds,res,kc,iorig,ndims,idims,lperiod,
     $     icoords,iLcoords,myside,myorig,myorig1,myorig2,myorig3,
     $     icommcart,mycartid,mpiid,lflag,out,inn)
c        Don't do anything else
         return
      endif



c Initialize right hand side and x, the updated potential

         

      do k=1,myside(3)
         do j=1,myside(2)
            do i=1,myside(1)
               index=myorig+(i-1)*iLs(1)+(j-1)*iLs(2)+(k-1)*iLs(3)
               x(index)=u(index)
               b(index)=exp(u(index))*(1-u(index))-q(index)
               if(out) then
                  if(i.eq.myside(1)-1) then
                     b(index)=b(index) -apc(myorig1+i-1)*gpc(myorig2+j-1
     $                    ,myorig3+k-1,4)
                  elseif(i.eq.myside(1))then
                     x(index)=0.
                  endif
               endif
               if(inn) then
                  if(i.eq.2) then
                     b(index)=b(index)-bpc(myorig1+i-1) *u(index-iLs(1))
c     Formally set the potential at the probe edge to zero, since the
c     inner boundary condition lies in the right hand side of the equation (b)
                  elseif(i.eq.1)then
                     x(index)=0.
c                    For debugging, set b to zero as well
                     b(index)=0.
                  endif
               endif
            enddo
         enddo
      enddo


c Outputs Ax, where A is the cg matrix
      call atimesmpi(myside(1),myside(2),myside(3),Li,Lj,Lk,x(myorig)
     $     ,res(myorig),u(myorig),apc(myorig1),bpc(myorig1) ,cpc(myorig1
     $     +Li*myorig2) ,dpc(myorig1+Li*myorig2),epc(myorig1+Li*myorig2)
     $     ,fpc(myorig1+Li *myorig2) ,gpc(myorig2,myorig3,1),out,
     $     .false.)

      
c     Calculate the initial residual r=b-Ax, where x is the potential at
c     the previous time-step. Whith the conjugate gradient method, the
c     first search direction is the first residual
      

      
      do k=2,myside(3)-1
         do j=2,myside(2)-1
            do i=1,myside(1)-1
               index=myorig+(i-1)*iLs(1)+(j-1)*iLs(2)+(k-1)*iLs(3)
               res(index)=b(index)-res(index)
               if(inn.and.i.eq.1) res(index)=0.
c              The following line is required for the bcg method
               resr(index)=res(index)
            enddo
         enddo
      enddo


c For the next atimesmpi, need res also on the shadow cells
      call bbdy(cg_comm,iLs,iuds,res,icg_k,iorig,ndims,idims,lperiod,
     $     icoords,iLcoords,myside,myorig,myorig1,myorig2,myorig3,
     $     icommcart,mycartid,mpiid,lflag,out,inn)


      

c     Following call used for minimum residual method
      if (.not. lbcg) then
         call atimesmpi(myside(1),myside(2),myside(3),Li,Lj,Lk
     $     ,res(myorig)
     $     ,resr(myorig),u(myorig),apc(myorig1) ,bpc(myorig1)
     $     ,cpc(myorig1+Li*myorig2) ,dpc(myorig1+Li *myorig2)
     $     ,epc(myorig1+Li*myorig2) ,fpc(myorig1+Li*myorig2)
     $     ,gpc(myorig2,myorig3 ,1),out,
     $     .false.)
      endif


      call asolvempi(myside(1),myside(2),myside(3),Li,Lj,Lk
     $     ,res(myorig),z(myorig),u(myorig),apc(myorig1),fpc(myorig1
     $     +Li*myorig2),gpc(myorig2,myorig3 ,1),out)


c     Start Main iteration  

      
      


      do icg_k=1,icg_mi


c Do block boundary communications

         call asolvempi(myside(1),myside(2),myside(3),Li,Lj,Lk
     $        ,resr(myorig),zz(myorig),u(myorig),apc(myorig1)
     $        ,fpc(myorig1 +Li*myorig2),gpc(myorig2,myorig3 ,1),out)

         
         
         bknum=0.
         
         do k=2,myside(3)-1
            do j=2,myside(2)-1
               do i=2,myside(1)-1
                  index=myorig+(i-1)*iLs(1)+(j-1)*iLs(2)+(k-1)*iLs(3)
                  bknum=bknum+z(index)*resr(index)
               enddo
            enddo
         enddo
         
c     Sum bknum over all the participating nodes
         call MPI_ALLREDUCE(bknum,bknumR,1,MPI_REAL, MPI_SUM
     $        ,icommcart,ierr)
         bknum=bknumR
         
         if(icg_k.eq.1) then
            do k=2,myside(3)-1
               do j=2,myside(2)-1
                  do i=2,myside(1)-1
                     index=myorig+(i-1)*iLs(1)+(j-1)*iLs(2)+(k-1)*iLs(3)
                     p(index)=z(index)
                     pp(index)=zz(index)
                  enddo
               enddo
            enddo
            i=1
            if(inn) then
               do k=2,myside(3)-1
                  do j=2,myside(2)-1
                     index=myorig+(i-1)*iLs(1)+(j-1)*iLs(2)+(k-1)*iLs(3)
                     p(index)=0.
                     pp(index)=0.
                  enddo
               enddo
            endif
                  
         else
            bk=bknum/bkden
            do k=2,myside(3)-1
               do j=2,myside(2)-1
                  do i=2,myside(1)-1
                     index=myorig+(i-1)*iLs(1)+(j-1)*iLs(2)+(k-1)*iLs(3)
                     p(index)=bk*p(index)+z(index)
                     pp(index)=bk*pp(index)+zz(index)
                  enddo
               enddo
            enddo
            i=1
            if(inn) then
               do k=2,myside(3)-1
                  do j=2,myside(2)-1
                     index=myorig+(i-1)*iLs(1)+(j-1)*iLs(2)+(k-1)*iLs(3)
                     p(index)=0.
                     pp(index)=0.
                  enddo
               enddo
            endif
                     
         endif

     

         bkden=bknum

         call bbdy(cg_comm,iLs,iuds,p,icg_k,iorig,ndims,idims,lperiod,
     $        icoords,iLcoords,myside,myorig,myorig1,myorig2,myorig3,
     $        icommcart,mycartid,mpiid,lflag,out,inn)


         call atimesmpi(myside(1),myside(2),myside(3),Li,Lj,Lk,p(myorig)
     $        ,z(myorig),u(myorig),apc(myorig1),bpc(myorig1),cpc(myorig1
     $        +Li*myorig2) ,dpc(myorig1+Li*myorig2),epc(myorig1+Li
     $        *myorig2) ,fpc(myorig1+Li *myorig2) ,gpc(myorig2,myorig3,1
     $        ) ,out,
     $     .false.)

         
         akden=0.
         
         do k=2,myside(3)-1
            do j=2,myside(2)-1
               do i=2,myside(1)-1
                  index=myorig+(i-1)*iLs(1)+(j-1)*iLs(2)+(k-1)*iLs(3)
                  akden=akden+z(index)*pp(index)
               enddo
            enddo
         enddo
c Sum akden over all the participating nodes
         call MPI_ALLREDUCE(akden,akdenR,1,MPI_REAL, MPI_SUM
     $        ,icommcart,ierr)
         akden=akdenR

         ak=bknum/akden
         
         call bbdy(cg_comm,iLs,iuds,pp,icg_k,iorig,ndims,idims,lperiod,
     $        icoords,iLcoords,myside,myorig,myorig1,myorig2,myorig3,
     $        icommcart,mycartid,mpiid,lflag,out,inn)

c        Implement bcg option by using lbcg as transpose flag
         call atimesmpi(myside(1),myside(2),myside(3),Li,Lj,Lk,
     $     pp(myorig),zz(myorig),u(myorig),apc(myorig1),bpc(myorig1)
     $     ,cpc(myorig1+Li*myorig2) ,dpc(myorig1+Li*myorig2)
     $     ,epc(myorig1+Li*myorig2) ,fpc(myorig1+Li *myorig2)
     $     ,gpc(myorig2,myorig3,1), out,
     $     lbcg)

         deltamax=0.
         
         do k=2,myside(3)-1
            do j=2,myside(2)-1
               do i=2,myside(1)-1
                  index=myorig+(i-1)*iLs(1)+(j-1)*iLs(2)+(k-1)*iLs(3)
                  delta=ak*p(index)
                  x(index)=x(index)+delta
                  if(abs(delta).gt.deltamax)deltamax=abs(delta)
                  res(index)=res(index)-ak*z(index)
                  resr(index)=resr(index)-ak*zz(index)
               enddo
            enddo
         enddo



c     Test convergence (Only every 5 steps and after a reasonable amount
c     of iterations)
         if (mod(icg_k+2,5).eq.0.and.icg_k.gt.0.8*icg_prec-5) then
            call testifconverged(cg_eps,deltamax,
     $           lconverged,cg_comm)
c            if(mpiid.eq.1)   write(*,*) icg_k,deltamax,lconverged
         endif

         if(lconverged.and.icg_k.ge.2)goto 11


         call asolvempi(myside(1),myside(2),myside(3),Li,Lj,Lk
     $        ,res(myorig),z(myorig),u(myorig),apc(myorig1) ,fpc(myorig1
     $        +Li*myorig2),gpc(myorig2,myorig3 ,1),out)

         
      enddo
c     After a loop is finished, fortran seems to increase the counter one more
      icg_k=icg_k-1

c-------------------------------------------------------------------
 11   continue

c Do the final mpi_gather [or allgather if all processes need
c the result].
      kc=-1

      call bbdy(cg_comm,iLs,iuds,x,kc,iorig,ndims,idims,lperiod,
     $        icoords,iLcoords,myside,myorig,myorig1,myorig2,myorig3,
     $        icommcart,mycartid,mpiid,lflag,out,inn)

c For debugging, also gather b for saving
      call bbdy(cg_comm,iLs,iuds,b,kc,iorig,ndims,idims,lperiod,
     $        icoords,iLcoords,myside,myorig,myorig1,myorig2,myorig3,
     $        icommcart,mycartid,mpiid,lflag,out,inn)


      cg_del=deltamax
      ierr=icg_k


c Indirection is needed here because otherwise the finalize call
c seems to cause the return to fail. Probably unnecessary.
c      mpiid=myid

      end


c***********************************************************************
c Outputs res=Ax, where A is the finite volumes stiffness matrix
      subroutine atimesmpi(ni,nj,nk,Li,Lj,Lk,x,res,u,a,b,c,d,e,f,g,out,
     $  ltrnsp)

      
c      integer ni,nj,nk
c      integer Li,Lj,Lk

      logical out
      real x(Li,Lj,nk),res(Li,Lj,nk),u(Li,Lj,nk)
      real a(ni),b(ni),c(Li,nj),d(Li,nj),e(Li,nj),f(Li,nj),g(Lj,Lk,5)
      logical ltrnsp


c The numbering starts at 1 and ends at n, but 1 and n are the ghost cells

      if (ltrnsp) then


c Matrix multiplication
      do k=2,nk-1
         do j=2,nj-1
            do i=2,ni-1
               res(i,j,k) = b(i+1)*x(i+1,j,k)
     $           + a(i-1)*x(i-1,j,k)
     $           + d(i,j+1)*x(i,j+1,k)
     $           + c(i,j-1)*x(i,j-1,k)
     $           + e(i,j)*(x(i,j,k+1)+x(i,j,k-1))
     $           - (f(i,j)+exp(u(i,j,k)))*x(i,j,k)
            enddo
         enddo
      enddo

c Take care of outer boundary condition
      if(out) then
         do k=2,nk-1
            do j=2,nj-1
               i=ni-2
               res(i,j,k) = res(i,j,k)
     $           + g(j,k,1)*a(i+1)*x(i+1,j,k)
               i=ni-1
               res(i,j,k) = res(i,j,k)
     $           + g(j+1,k,2)*a(i)*x(i,j+1,k)
     $           + g(j-1,k,3)*a(i)*x(i,j-1,k)
     $           + g(j,k,5)*a(i)*x(i,j,k)
            enddo
         enddo
      endif



      else

c Matrix multiplication
      do k=2,nk-1
         do j=2,nj-1
            do i=2,ni-2
               res(i,j,k) = a(i)*x(i+1,j,k)
     $           + b(i)*x(i-1,j,k)
     $           + c(i,j)*x(i,j+1,k)
     $           + d(i,j)*x(i,j-1,k)
     $           + e(i,j)*(x(i,j,k+1)+x(i,j,k-1))
     $           - (f(i,j)+exp(u(i,j,k)))*x(i ,j,k)
            enddo
         enddo
      enddo


c Take care of outer boundary condition
      i=ni-1
      if(out) then
         do k=2,nk-1
            do j=2,nj-1
               x(i+1,j,k) = g(j,k,1)*x(i-1,j,k)
     $           + g(j,k,2)*x(i,j-1,k)
     $           + g(j,k,3)*x(i,j+1,k)
     $           + 0*g(j,k,4)
     $           + g(j,k,5)*x(i,j,k)

               res(i,j,k) = a(i)*x(i+1,j,k)
     $           + b(i)*x(i-1,j,k)
     $           + c(i,j)*x(i,j+1,k)
     $           + d(i,j)*x(i,j-1,k)
     $           + e(i,j)*(x(i,j,k+1)+x(i,j,k-1))
     $           - (f(i,j)+exp(u(i,j,k)))*x(i ,j,k)
            enddo
         enddo
      else
         do k=2,nk-1
            do j=2,nj-1
               res(i,j,k) = a(i)*x(i+1,j,k)
     $           + b(i)*x(i-1,j,k)
     $           + c(i,j)*x(i,j+1,k)
     $           + d(i,j)*x(i,j-1,k)
     $           + e(i,j)*(x(i,j,k+1)+x(i,j,k-1))
     $           - (f(i,j)+exp(u(i,j,k)))*x(i ,j,k)
            enddo
         enddo
      endif
      
      endif

      
      end

c***********************************************************************
c     Preconditioning subroutine. If Atilde is the preconditioning
c     matrix, returns z=Atilde^-1*b
      subroutine asolvempi(ni,nj,nk,Li,Lj,Lk,b,z,u,a,f,g,out)

      
c      integer ni,nj,nk
c      integer Li,Lj,Lk
      logical out
      real b(Li,Lj,nk),z(Li,Lj,nk),u(Li,Lj,nk)
      real a(ni),f(Li,nj),g(Lj,Lk,5)

c Matrix multiplication
      
      
      do k=2,nk-1
         do j=2,nj-1
            do i=2,ni-1
               z(i,j,k)=b(i,j,k)/(-f(i,j)-exp(u(i,j,k)))
            enddo
         enddo
      enddo

c Outer boundary
      i=ni-1
      if(out) then
         do k=2,nk-1               
            do j=2,nj-1
               z(i,j,k)=b(i,j,k)/(-f(i,j)-exp(u(i,j,k))+a(i)*g(j,k,5)) 
            enddo
         enddo
      else
         do k=2,nk-1
            do j=2,nj-1
               z(i,j,k)=b(i,j,k)/(-f(i,j)-exp(u(i,j,k)))     
            enddo
         enddo
      endif
         
      end



c***********************************************************************
c The challenge here is to ensure that all processes decide to end
c at the same time. If not then a process will hang waiting for message.
c So we have to universalize the convergence test. All block must be
c converged, but the total spread depends on multiple blocks.
      subroutine testifconverged(cg_eps,delta,lconverged,cg_comm)
      include 'mpif.h'
      logical lconverged
      real delta,deltaR,cg_eps
      integer cg_comm
      
c Here we need to allreduce the data, selecting the maximum values,

c  !!!!! not doing it in place
c Loki's version of mpich seems not to like MPI_IN_PLACE

      call MPI_ALLREDUCE(delta,deltaR,1,MPI_REAL, MPI_MAX,cg_comm,ierr)

      delta=deltaR

      if(delta.lt.cg_eps) then
         lconverged=.true.
      else
         lconverged=.false.
      endif

      end

c***********************************************************************
c Cut here and throw the rest away for the basic routines
c***********************************************************************
      subroutine fcalc3Dpar(Li,Lj,Lk,ni,nj,nk,mi,eps,k,cg_comm,myid2)


      include 'piccom.f'
      include 'errcom.f'
      integer Li,Lj,Lk,ni,nj,nk
      integer myid2
      integer cg_comm
      integer nd
      parameter (nd=3,nd2=nd*2)

      real cg_eps,cg_del,eps
      integer icg_k,icg_mi,mi
      

      common /cg3dctl/icg_mi,cg_eps,cg_del,icg_k

      real b(0:nrsize,0:nthsize,0:npsisize),x(0:nrsize,0:nthsize
     $     ,0:npsisize),p(0:nrsize,0:nthsize,0:npsisize),res(0:nrsize
     $     ,0:nthsize ,0:npsisize),z(0:nrsize,0:nthsize,0:npsisize)
     $     ,pp(0:nrsize,0:nthsize ,0:npsisize),resr(0:nrsize,0:nthsize
     $     ,0:npsisize),zz(0:nrsize ,0:nthsize,0:npsisize)

c     Variables used for calculating matrix A for debugging
      real inputvect(0:nrsize,0:nthsize,0:npsisize),
     $  outputvect(0:nrsize,0:nthsize,0:npsisize)
      integer n1,n2,n3,j,k,l,m,n,o,jkl,mno


c testing arrays
      integer iuds(nd),ifull(nd)

c     Initialize variables for debugging
      lAdebug = .false.


      icg_mi=mi
      cg_jac=jac
      cg_eps=eps


      ifull(1)=Li
      ifull(2)=Lj
      ifull(3)=Lk
      iuds(1)=ni
      iuds(2)=nj
      iuds(3)=nk
      
      
      ictl=1
  
      call cg3dmpi(cg_comm,Li,Lj,Lk,ni,nj,nk,bcphi, phi(1,0,0)
     $     ,rho(1,0,0),ictl,ierr,myid,idim1,idim2,idim3,apc(1),bpc(1)
     $     ,cpc(1,0),dpc(1,0),epc(1,0),fpc(1,0),gpc(0,0,1)
     $     ,b(1,0,0),x(1,0,0)
     $     ,p(1,0,0) ,res(1,0,0),z(1,0,0) ,pp(1,0,0),resr(1,0,0) ,zz(1,0
     $     ,0),lbcg)

      
      
c For debugging, save matrix A and its transpose, as well as b and x
      if (lsavemat .and. stepcount.eq.saveatstep) then
c        Set flag for cg3dmpi to only multiply by A
         lAdebug = .true.
         n1 = ni-1
         rshieldingsave = n1
         n2 = nthused
         n3 = npsiused
c        First, initialize bsave and xsave just in case
         do k=0,npsisizesave
            do j=0,nthsizesave
               do i=1,nrsizesave-1
                  bsave(i,j,k) = 0.
                  xsave(i,j,k) = 0.
               enddo
            enddo
         enddo
         do k=0,n3
            do j=0,n2
               do i=1,n1
                  bsave(i,j,k) = b(i,j,k)
                  xsave(i,j,k) = x(i,j,k)
               enddo
            enddo
         enddo
         do j=1,n3
            do k=1,n2
               do l=1,n1
c                 Pass unit vectors to atimes to build A
                  inputvect(l,k,j) = 1.
                  lAtranspose = .false.
c                 Note that all the arrays passed in the x(*) form
c                   are offset by one element so that when they are
c                   passed from within that routine as x(0) it is the
c                   actual beginning of the array that is passed.
                  call cg3dmpi(cg_comm,Li,Lj,Lk,ni,nj,nk,bcphi
     $              ,phi(1,0,0),rho(1,0,0),ictl,ierr,myid,idim1,idim2
     $              ,idim3,apc(1),bpc(1),cpc(1,0),dpc(1,0),epc(1,0)
     $              ,fpc(1,0),gpc(0,0,1),b(1,0,0)
     $              ,inputvect(1,0,0),p(1,0,0),outputvect(1,0,0)
     $              ,z(1,0,0),pp(1,0,0),resr(1,0,0),zz(1,0,0),lbcg)
                  do m=1,n3
                     do n=1,n2
                        do o=1,n1
                           Asave(o,n,m,l,k,j) = outputvect(o,n,m)
c                          atimes may change input vector, so reset
                           inputvect(o,n,m) = 0.
c                          reset output to to be safe
                           outputvect(o,n,m) = 0.
                        enddo
                     enddo
                  enddo
c                 Pass unit vectors to atimes to build A'
                  inputvect(l,k,j) = 1.
                  lAtranspose = .true.
                  call cg3dmpi(cg_comm,Li,Lj,Lk,ni,nj,nk,bcphi
     $              ,phi(1,0,0),rho(1,0,0),ictl,ierr,myid,idim1,idim2
     $              ,idim3,apc(1),bpc(1),cpc(1,0),dpc(1,0),epc(1,0)
     $              ,fpc(1,0),gpc(0,0,1),b(1,0,0)
     $              ,inputvect(1,0,0),p(1,0,0),outputvect(1,0,0)
     $              ,z(1,0,0),pp(1,0,0),resr(1,0,0),zz(1,0,0),lbcg)
                  do m=1,n3
                     do n=1,n2
                        do o=1,n1
                           Atsave(o,n,m,l,k,j) = outputvect(o,n,m)
c                          atimes may change input vector, so reset
                           inputvect(o,n,m) = 0.
c                          reset output to to be safe
                           outputvect(o,n,m) = 0.
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

c Write x, the temporary potential file, to phi, and find maxchange
      if(myid2.eq.0) then
         do k=1,npsiused
            do j=1,nthused
               do i=2,ni-1
                  phi(i,j,k)=x(i,j,k)
               enddo
            enddo
         enddo
      endif


      k=icg_k
      
      end
