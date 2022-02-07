      subroutine fltmp(u)   
      use gem_com
      use equil
      use fft_wrapper
      implicit none
      INTEGER :: i,j,k,k1,l,m,n,nky,deltam
      real(8) :: u(0:imx,0:jmx,0:1)
      complex(8) :: cdum,cdum1
      complex(8) :: umn(0:imx,0:50),v(0:imx,0:50)
      real(8) :: dum,dum1,r,s,zeta,ydum,ydum1,udum,udum1,wy0,wy1,dzeta
      real(8) :: dthf(0:imx)
      integer :: mpolpf(0:imx,0:50)

      nky = 1
      n = nky*lyfra !toroidal mode number
      deltam = 8
      dzeta = pi2/lyfra/jmx
      do i = 0,imx
         dthf(i) = gthf(i,1)-gthf(i,0)
      end do
      do i = 0,imx
         do m = 0, 2*deltam-1
            mpolpf(i,m) = floor(-n*gsf(i)-deltam+m)  !m+nq in [-deltam,deltam]
         end do
         mpolpf(i,2*deltam) = mpolpf(i,2*deltam-1)+1
      end do

      umn = 0.
      do i = 1,imx-1
         do m = 0,2*deltam
            cdum = 0.
            cdum1 = 0.
            do j = 0,jmx-1
               zeta = j*dzeta

               ydum = modulo(r0/q0*(gsf(i)*gthf(i,0)-zeta),ly)
               l = int(ydum/dy)
               wy0=float(l+1)-ydum/dy
               wy1=1.-wy0
               udum = wy0*u(i,l,0)+wy1*u(i,l+1,0)
               dum = -(gsf(i)*gthf(i,0)-zeta)  !zeta-gsf(i)*gthf(i,0)
               cdum = cdum+udum*exp(-IU*n*dum)*dzeta  !integration over a wedge

               ydum1 = modulo(r0/q0*(gsf(i)*gthf(i,1)-zeta),ly)
               l = int(ydum1/dy)
               wy0=float(l+1)-ydum1/dy
               wy1=1.-wy0
               udum1 = wy0*u(i,l,1)+wy1*u(i,l+1,1)
               dum1 = -(gsf(i)*gthf(i,1)-zeta)  !zeta-gsf(i)*gthf(i,0)
               cdum1 = cdum1+udum1*exp(-IU*n*dum1)*dzeta  !integration over a wedge

               
            end do
            cdum = cdum*lyfra !integration over 2*pi
            cdum1 = cdum1*lyfra !integration over 2*pi
            umn(i,m) = (cdum*exp(-IU*(n*gsf(i)+mpolpf(i,m))*gthf(i,0))*dthf(i)+cdum1*exp(-IU*(n*gsf(i)+mpolpf(i,m))*gthf(i,1))*dthf(i))/2
         end do
      end do


      cnt = (imx-1)*(2*deltam+1)
      call mpi_allreduce(umn(1:imx-1,0:2*deltam),v(1:imx-1,0:2*deltam),cnt,MPI_DOUBLE_COMPLEX,mpi_sum, &
                              tube_comm,ierr)
      v = v/pi2**2


!inverse transform
      u = 0.
      do i = 1,imx-1
         do j = 0,jmx
            do k = 0,1
               zeta = modulo(gsf(i)*gthf(i,k)-yg(j)*q0/r0,pi2/lyfra) !zeta at (x,y,z)
               dum = zeta-gsf(i)*gthf(i,k)
               do m = 0,2*deltam
                  cdum = v(i,m)*exp(IU*(n*gsf(i)+mpolpf(i,m))*gthf(i,k)+IU*n*dum)
                  u(i,j,k) = u(i,j,k)+cdum+dconjg(cdum)
               end do
            end do
         end do
      end do

      call enfxy(u(:,:,:))
      call enfz(u)

      return
      end
