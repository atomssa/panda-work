      program pbarp_pippim_cos
      implicit real*8(a-h,o-z)
      
      open(unit=6,file='pbarp_pippim_cos.dat',status='unknown',
     &     form='formatted')
      
      xmin=-0.94d0
      xmax= 0.94d0
      dx=0.04d0
      xmax=xmax+dx/10.d0
      
      do x=xmin,xmax,dx
        write(6,9100) x
      end do
      
 9100 format(1x,f6.2)
      end
