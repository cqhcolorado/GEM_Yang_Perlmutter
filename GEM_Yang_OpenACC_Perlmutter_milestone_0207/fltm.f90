      subroutine fltm(u)   
      use gem_com
      use equil
      use fft_wrapper
      implicit none
      INTEGER :: i,j,k,k1,l,m,n,nky,ifirst,numpol,deltam
      real(8) :: u(0:imx,0:jmx,0:1)
      complex(8) :: cdum,cdum1
      complex(8) :: umn(0:imx,0:50),v(0:imx,0:50)
      real :: dum,r,s,zeta,ydum,udum,wy0,wy1,dzeta
      integer,dimension(:,:),allocatable :: mpolpf
      real,dimension(:),allocatable :: dthf
      complex(8),dimension(:,:,:),allocatable :: efac1 
      complex(8),dimension(:,:,:),allocatable :: efac2         

      save ifirst,mpolpf,dthf,efac1,efac2
      nky = 1
      n = nky*lyfra             !toroidal mode number
      deltam = 16
      numpol = 2*deltam+1
      dzeta = pi2/lyfra/jmx

      if(ifirst .ne. -99)then
         allocate(dthf(1:imx-1),mpolpf(1:imx-1,0:numpol-1),efac1(1:imx-1,0:jmx-1,0:1),efac2(1:imx-1,0:1,0:numpol-1))
         do i = 0,imx
            dthf(i) = (gthf(i,1)-gthf(i,0))*isgnq
         end do
         do i = 0,imx
            do m = 0, numpol-1
               mpolpf(i,m) = floor(-n*gsf(i)-deltam+m) !m+nq in [-deltam,deltam]
            end do
         end do

         do i = 1,imx-1
            do j = 0,jmx-1
               do k = 0,1
                  zeta = modulo(gsf(i)*gthf(i,k)-yg(j)*q0/r0,pi2/lyfra) !zeta at (x,y,z)
                  dum = zeta-gsf(i)*gthf(i,k)
                  efac1(i,j,k) = exp(IU*n*dum)
               end do
            
               do k = 0,1
                  do m = 0,numpol-1
                     efac2(i,k,m) = exp(IU*(n*gsf(i)+mpolpf(i,m))*gthf(i,k))
                  end do
               end do
            end do
         end do
         ifirst = -99
      end if

      umn = 0.
      do i = 1,imx-1
         do m = 0,numpol-1
            cdum = 0.
            cdum1 = 0.
            do j = 0,jmx-1
               udum = u(i,j,0) 
               cdum = cdum+udum/efac1(i,j,0)*dzeta 

               udum = u(i,j,1) 
               cdum1 = cdum1+udum/efac1(i,j,1)*dzeta                
            end do
            cdum = cdum*lyfra !integration over 2*pi
            cdum1 = cdum1*lyfra
            umn(i,m) = (cdum/efac2(i,0,m)+cdum1/efac2(i,1,m))*dthf(i)/2
         end do
      end do


      cnt = (imx-1)*numpol
      call mpi_allreduce(umn(1:imx-1,0:numpol-1),v(1:imx-1,0:numpol-1),cnt,MPI_DOUBLE_COMPLEX,mpi_sum, &
                              tube_comm,ierr)
      v = v/pi2**2


!inverse transform
      u = 0.
      do i = 1,imx-1
         do j = 0,jmx
            do k = 0,1
               zeta = modulo(gsf(i)*gthf(i,k)-yg(j)*q0/r0,pi2/lyfra) !zeta at (x,y,z)
               dum = zeta-gsf(i)*gthf(i,k)
               do m = 0,numpol-1
                  cdum = v(i,m)*efac2(i,k,m)*efac1(i,j,k)
                  u(i,j,k) = u(i,j,k)+cdum+dconjg(cdum)
               end do
            end do
         end do
      end do

      call enfxy(u(:,:,:))
      call enfz(u)

      return
      end
