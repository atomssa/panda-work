      program pbarp_pi0pi0_be
      implicit real*8(a-h,o-z)      
c
c
c     section  efficace differentielle pbar p ---> pi0 pi0  basse energie
c     Ref:R.S; Dulude et al.  Phys. Lett. 79B (1978) p: 329-334 et 335-339
c
c            
      character*4 chisqs
      character*130 chfich 
       
      common/allmasse2c/rmassea,rmasseb,rmasse1,rmasse2 
      common/cqp_labab12/qpa_lab(0:3),qpb_lab(0:3),qp1_lab(0:3),
     &                   qp2_lab(0:3) 
     
              
      dimension sqstab(10),a0tab(10),a2tab(10),a4tab(10),
     &    a6tab(10),a8tab(10),sigtot(10) 
         
      dimension tetatab(1000)

      dimension bboost_labvcm(3),bboost_cmvlab(3) 
                       
      dimension qpaps_up(0:3),qpbps_up(0:3),qp1ps_up(0:3),
     &  qp2ps_up(0:3)      
      
c      dimension qpaps_do(0:3),qpbps_do(0:3),qp1ps_do(0:3),
c     &  qp2ps_do(0:3)

c      dimension  qpacm_up(0:3),qpbcm_up(0:3)
      dimension  qp1cm_up(0:3),qp2cm_up(0:3)

                  
      
      open(unit=6,file='pbarp_pi0pi0_be.ou',status='unknown',
     &    form='formatted')
     
    
      pi=acos(-1.d0)
      fd=180.d0/pi

      rmassea=0.93827231d0
      rmasseb=0.93827231d0
      rmasse1=0.13498d0
      rmasse2=0.13498d0


      
      nsqs=8
      call atab_fill(nsqs,sqstab,a0tab,a2tab,a4tab,a6tab,a8tab,sigtot)
      
      do i=1,nsqs
        write(6,9100) sqstab(i),a0tab(i),a2tab(i),a4tab(i),a6tab(i),
     &      a8tab(i),sigtot(i)
      end do
c
c     section efficace en nanobarn/sr
c
      fact=1000.d0
      do isqs=1,nsqs
        a0tab(isqs)=a0tab(isqs)*1000.d0
        a2tab(isqs)=a2tab(isqs)*1000.d0	
        a4tab(isqs)=a4tab(isqs)*1000.d0	
        a6tab(isqs)=a6tab(isqs)*1000.d0	
        a8tab(isqs)=a8tab(isqs)*1000.d0
	sigtot(isqs)=sigtot(isqs)*1000.d0		
      end do

c
c
c      
      call prep_ang_lab(nteta,tetatab)
      
      write(6,*) '   '
      
      write(6,*) ' nteta=',nteta
      do iteta=1,nteta
        tetad=tetatab(iteta)*fd
	write(6,*) ' iteta tetatab',iteta,tetad
      end do            

c
c     calcul dans le centre de masse
c
      
      do 2000 isqs=1,nsqs     
        sqs=sqstab(isqs)
        q2=sqs**2
	s_mandel=q2
	write(6,*) '  sqs (GeV)=',sqs
	
	iisqs=nint(sqs*1000.d0)
	write(chisqs,'(i4)') iisqs		
	iunit=8
	chfich='xsect_pbarp_pi0pi0_be_cm_'//chisqs//'.ou'
c	write(6,9900) chfich
	open(unit=iunit,file=chfich,status='unknown',
     &       form='formatted')
	do 1000 iteta=1,nteta
	  teta1=tetatab(iteta)
	  teta1cm=teta1
	  a0=a0tab(isqs)
	  a2=a2tab(isqs)	  
	  a4=a4tab(isqs)	  
	  a6=a6tab(isqs)	  
	  a8=a8tab(isqs)	  
	  x=cos(teta1cm)
	  call xsect_pi0pi0_leg(x,a0,a2,a4,a6,a8,xsect)
          write(6,*) ' teta1cm x  xsect',teta1cm*fd,x,xsect
	  write(iunit,9200) teta1cm*fd,x,xsect
 1000	continue
      close  (iunit)
 2000 continue

      
c
c     calcul dans le laboratoire
c
      
      do 4000 isqs=1,nsqs     
        sqs=sqstab(isqs)
        q2=sqs**2
	s_mandel=q2
	write(6,*) '  sqs (GeV)=',sqs
	
        s_mandel=q2
      
        t_cin=(q2-4.d0*rmassea**2)/(2.d0*rmassea)     

        write(6,*) ' q**2=',q2
        write(6,*) ' energie cinetique de p bar (GeV)=',t_cin
        write(6,*) '   '	
	
	
        do mu=0,3
          qpaps_up(mu)=0.d0
	  qpbps_up(mu)=0.d0
          qp1ps_up(mu)=0.d0
          qp2ps_up(mu)=0.d0
        end do
        qpaps_up(0)=rmassea+t_cin
        ea=qpaps_up(0)
        pa2=ea**2 - rmassea**2
        pa=sqrt(pa2)	
        qpaps_up(3)=sqrt(pa2)
        qpbps_up(0)=rmasseb
	
	
		
c
c         transformation de lorentz pure
c

        do i=1,3
          bboost_labvcm(i)=qpaps_up(i)/(qpaps_up(0)+rmasseb)
          bboost_cmvlab(i)=-bboost_labvcm(i)
        end do	
	
        e12cm=s_mandel/4.d0
        e1cm=sqrt(e12cm)
        p12cm=e12cm-rmasse1**2
        p1cm=sqrt(p12cm)
        e2cm=e1cm
        p2cm=p1cm	

        write(6,*)  ' p1cm',p1cm
	
	
	iisqs=nint(sqs*1000.d0)
	write(chisqs,'(i4)') iisqs		
	iunit=8
	chfich='xsect_pbarp_pi0pi0_be_lab_'//chisqs//'.ou'
c	write(6,9900) chfich
	open(unit=iunit,file=chfich,status='unknown',
     &       form='formatted')
        phi1=0.d0
	do 3000 iteta=1,nteta
	  teta1=tetatab(iteta)
	  write(6,*) ' teta1 ',teta1*fd
          call cin2cppbarpi0pi0(t_cin,teta1,phi1)	  
	  
          do i=0,3
            qqdif=qpa_lab(i) + qpb_lab(i) - qp1_lab(i) - qp2_lab(i)
	    if(qqdif.gt.1.d-8) then
              write(6,*) ' i qqdif ',i,qqdif
	    end if
          end do
	  	  
          do i=0,3
            qpaps_up(i)=qpa_lab(i)
            qpbps_up(i)=qpb_lab(i)
            qp1ps_up(i)=qp1_lab(i)
            qp2ps_up(i)=qp2_lab(i)
          end do	  
	  
	  	  	  
          call blorentz(qp1ps_up,bboost_labvcm,qp1cm_up)
          call blorentz(qp2ps_up,bboost_labvcm,qp2cm_up)	  	  
      
          p1xcm=qp1cm_up(1)
          p1ycm=qp1cm_up(2)	  
          p1zcm=qp1cm_up(3)	  

          p12cmp=p1xcm**2 + p1ycm**2 +p1zcm**2
	  p1cmp=sqrt(p12cmp)
      
          call SPHER(p1Xcm,p1Ycm,p1Zcm,pp1cmp,TETA1cm,PHI1cm)      
      	  write(6,*) ' teta1cm  phi1cm',teta1cm*fd,phi1cm*fd
	  write(6,*) ' p1cmp  pp1cmp',p1cmp,pp1cmp
	  
	  
	  a0=a0tab(isqs)
	  a2=a2tab(isqs)	  
	  a4=a4tab(isqs)	  
	  a6=a6tab(isqs)	  
	  a8=a8tab(isqs)	  
	  x=cos(teta1cm)
	  call xsect_pi0pi0_leg(x,a0,a2,a4,a6,a8,xsect)

	  xsect_cm=xsect
          write(6,*) '   x  xsect_cm',x,xsect_cm
          pa2=qpaps_up(1)**2 + qpaps_up(2)**2 + qpaps_up(3)**2 
          p12=qp1ps_up(1)**2 + qp1ps_up(2)**2 + qp1ps_up(3)**2
          pa=sqrt(pa2)
          p1=sqrt(p12)
          ea=qpaps_up(0)
          e1=qp1ps_up(0)
          e2=qp2ps_up(0)
          p2=sqrt(e2**2 - rmasse2**2)

          fac_xsect= 389.39d0 * 1000.d0
          fac_xsect=fac_xsect/2.d0


          xlab=abs( e2 + e1 * (1.d0 - pa/p1 * cos(teta1)) )
          den_xsect=16.d0 * (2.d0*pi)**2 * rmasseb * pa * xlab
          rnum_xsect= fac_xsect * p1        
	  	  	  
	           
          xcm=sqrt(s_mandel)
          den_xsect_cm=16.d0 * (2.d0*pi)**2 * rmasseb * pa * xcm	  
          rnum_xsect_cm= fac_xsect * p1cm	  	  	  

	  fjac_labscm= rnum_xsect/den_xsect	  
          fjac_labscm=fjac_labscm /(rnum_xsect_cm/den_xsect_cm) 	  
	  
	  xsect=xsect_cm * fjac_labscm	  
	  write(6,*)  ' fjac_labscm',fjac_labscm
          write(6,*) ' teta1  xsect',teta1*fd,xsect
	  if(xsect.lt.0.d0) stop
	  	  
	  write(iunit,9200) teta1*fd,xsect
 3000	continue
      close  (iunit)
 4000 continue
      
 9100 format(1x,10(1x,f5.2))
 9200 format(1x,f10.4,10(1x,e13.6))      
 9900 format(a80)     
      end
************************************************************************* 
************************************************************************       
      subroutine xsect_pi0pi0_leg(x,a0,a2,a4,a6,a8,xsect)      
      implicit real*8(a-h,o-z)      
      
      dimension plt(0:8)
      
      do L=0,8,2
        LL=L
        call PL0(LL,X,PL)
	plt(LL)=pl
      end do
      
      xsect=a0*plt(0)+a2*plt(2)+a4*plt(4)+a6*plt(6)+a8*plt(8)

      kpass=0
      if(xsect.lt.0.d0) kpass=1
      if(kpass.ne.0) then
        write(6,*) ' a0,a2,a4,a6,a8',a0,a2,a4,a6,a8
        write(6,*) ' x plt(0,2,4,6,8)',x,plt(0),plt(2),plt(4),
     &   plt(6),plt(8)
      end if	            
      return
      end
************************************************************************* 
************************************************************************      
      SUBROUTINE PL0(L,X,PL)
      IMPLICIT REAL*8(A-H,O-Z)
      PL=1.D0
      IF(L.EQ.0) GO TO 50
      PL=X
      IF(L.EQ.1) GO TO 50
      P0=1.D0
      P1=X
      LDEB=2
      DO 20 N=LDEB,L
      RN=DBLE(N)
      P2=2.D0*X*P1-P0-(X*P1-P0)/RN
      P0=P1
      P1=P2
 20   CONTINUE
      PL=P2
 50   CONTINUE
      RETURN
      END      
************************************************************************* 
************************************************************************      
      subroutine prep_ang_lab(nteta,tetatab)      
      implicit real*8(a-h,o-z)      
      dimension tetatab(*)       
      pi=acos(-1.d0) 
      fd=180.d0/pi
      
      dteta=0.5d0/fd
      
      ii=0
      tetamin=0.d0/fd     
      do iteta=1,11
        ii=ii+1
        tetatab(ii)=tetamin + dble(iteta-1)*dteta
      end do

      
     


c      write(6,*) ii      
      tetamin=tetatab(ii)            
      dteta=1./fd            
      do iteta=1,5
        ii=ii+1
        tetatab(ii)=tetamin + dble(iteta)*dteta
      end do      


c      write(6,*) ii       
      tetamin=tetatab(ii)            
      dteta=2./fd            
      do iteta=1,80
        ii=ii+1
        tetatab(ii)=tetamin + dble(iteta)*dteta
      end do


c      write(6,*) ii      
      tetamin=tetatab(ii)            
      dteta=1./fd            
      do iteta=1,9
        ii=ii+1
        tetatab(ii)=tetamin + dble(iteta)*dteta
      end do


 

c      ii=ii+1
      tetatab(ii)=180.d0/fd

      nteta=ii
      
c      write(6,*) ' nteta=',nteta
c      do iteta=1,nteta
c        tetad=tetatab(iteta)*fd
c	write(6,*) ' iteta tetatab',iteta,tetad
c      end do

      
      
      return
      end
************************************************************************* 
************************************************************************
      subroutine atab_fill(n,sqstab,a0tab,a2tab,a4tab,a6tab,a8tab,
     &                       sigtot)
      implicit real*8(a-h,o-z)
c
c     coefficients pour obtenir la section efficace differentielle
c     en microbarn/sr
c      
      dimension sqstab(n),a0tab(n),a2tab(n),a4tab(n),
     &    a6tab(n),a8tab(n),sigtot(n)


c
c     sqrt(s) en GeV
c      
      pi=acos(-1.d0)

      sqstab(1)=2.10d0
      sqstab(2)=2.14d0
      sqstab(3)=2.20d0
      sqstab(4)=2.24d0
      sqstab(5)=2.30d0
      sqstab(6)=2.34d0
      sqstab(7)=2.40d0
      sqstab(8)=2.45d0
      
      a0tab(1)=9.60d0
      a0tab(2)=7.00d0      
      a0tab(3)=4.20d0      
      a0tab(4)=2.80d0      
      a0tab(5)=2.00d0      
      a0tab(6)=1.60d0      
      a0tab(7)=1.40d0
      a0tab(8)=1.10d0



c      a2tab(1)=10.d0
c      a2tab(2)= 5.d0
c      a2tab(3)= 0.5d0     
c      a2tab(4)=-1.0d0      
c      a2tab(5)=-1.0d0      
c      a2tab(6)=-0.4d0      
c      a2tab(7)= 0.5d0
c      a2tab(8)= 1.0d0


      a2tab(1)=10.d0
      a2tab(2)= 5.d0
      a2tab(3)= 0.5d0     
      a2tab(4)=-0.5d0      
      a2tab(5)=-1.0d0      
      a2tab(6)=-0.4d0      
      a2tab(7)= 0.5d0
      a2tab(8)= 1.0d0

      
      a4tab(1)=13.0d0      
      a4tab(2)= 7.7d0
      a4tab(3)= 3.5d0
      a4tab(4)= 1.5d0
      a4tab(5)= 1.0d0
      a4tab(6)= 1.2d0
      a4tab(7)= 1.7d0
      a4tab(8)= 1.7d0
      
c      a6tab(1)= -6.0d0
c      a6tab(2)= -8.5d0     
c      a6tab(3)= -9.0d0
c      a6tab(4)= -8.0d0     
c      a6tab(5)= -3.0d0
c      a6tab(6)= -1.0d0
c      a6tab(7)=  0.0d0
c      a6tab(8)=  0.5d0


      a6tab(1)= -6.0d0
      a6tab(2)= -8.0d0     
      a6tab(3)= -8.0d0
      a6tab(4)= -6.0d0     
      a6tab(5)= -3.0d0
      a6tab(6)= -1.0d0
      a6tab(7)=  0.0d0
      a6tab(8)=  0.5d0

      
      a8tab(1)= 0.5d0
      a8tab(2)= 1.0d0
      a8tab(3)= 1.8d0     
      a8tab(4)= 2.5d0     
      a8tab(5)= 3.2d0 
      a8tab(6)= 3.7d0     
      a8tab(7)= 2.5d0     
      a8tab(8)= 1.5d0     


      do i=1,n
        sigtot(i)=2.d0*pi*a0tab(i)
      end do
      
      
          
      return
      end
************************************************************************* 
************************************************************************
      subroutine cin2cppbarpi0pi0(t_cin,teta1,phi1)
      implicit real*8(a-h,o-z)

      common/allmasse2c/rmassea,rmasseb,rmasse1,rmasse2
      common/cqp_labab12/qpa_lab(0:3),qpb_lab(0:3),qp1_lab(0:3),
     &                   qp2_lab(0:3)

      common/cipass/ipass

      
      dimension qpa(0:3),qpb(0:3),qp1(0:3),qq(0:3)
      dimension qp2(0:3)

      dimension itab(2),dift(2)
      dimension t1t(2)


      pi=acos(-1.d0)
      epsil_masse=1.d-4
      fd=180.d0/pi
      
      rma=rmassea
      rmb=rmasseb
      rm1=rmasse1
      rm2=rmasse2


      q=rma+rmb-(rm1+rm2)
      ta=t_cin
      
c      write(6,*) ' T_cin=',ta
      if(ipass.eq.1) then
        if(rm1.lt.epsil_masse) then
          write(6,*) ' calcul effectue avec rm1=rm2=0.'
        else
          write(6,*) ' calcul effectue avec rm1 et rm2 ne 0'
        end if
      end if
      

c      write(6,*) ' rma rmb ',rma,rmb
c      write(6,*) ' rm1 rm2 ',rm1,rm2
      
      pa2= ta* ( ta + 2.d0 * rma )
      pa=sqrt(pa2)
      ei=rma+ta+rmb
      ea=rma+ta
      pax=0.d0
      pay=0.d0
      paz=pa


      a=q*q + 2.d0 * q * ( rm2 + ta ) + 2.d0 * ta * (rm2 - rma )
      b= - 2.d0 * ( q + ta + rm1 + rm2 )

      gc=a*a


      qpa(0)=rma+ta
      qpa(1)=0.d0
      qpa(2)=0.d0
      qpa(3)=pa

      qpb(0)=rmb
      qpb(1)=0.d0
      qpb(2)=0.d0
      qpb(3)=0.d0

      teta=teta1
      phi=phi1
      
c      write(6,*) ' teta=',teta*fd
      
      cs=cos(teta)
      cs2=cs*cs
      gb=2.d0 * ( a*b -4.d0*rm1*pa2*cs2)
      ga=b*b-4.d0*pa2*cs2
      delta=gb*gb-4.d0*ga*gc
      
c      write(6,*)  ' delta cin=',delta
      
      
      if(abs(cs).lt.1.d-10) delta=1.d-12


      if(rm1.lt.epsil_masse) then
c        write(6,*) ' a  b',a,b
c        write(6,*) ' pa  cs ',pa,cs
        t1= -a/(2.d0*pa*cs + b)
c       write(6,*) ' - 3 t1 ',t1
      end if

      if(rm1.ge.epsil_masse) then
        if(abs(cs).lt.1.d-10) go to 10
        if(abs(delta).lt.1.d-10) then
	  delta=1.d-10
	  go to 10
	end if
	if(delta.lt.0.d0) then
	  write(6,*) ' delta negatif '
	  stop
	end if  
 10     continue      
      
      
      
        t1t(1)=(-gb+sqrt(delta))/(2.d0*ga)
        t1t(2)=(-gb-sqrt(delta))/(2.d0*ga)

        nsol=0

        do i=1,2
          itab(i)=0
        end do

        e1=t1t(1)+rm1
        if(e1.gt.rm1) then
          nsol=nsol+1
          itab(1)=1
        end if

        e1=t1t(2)+rm1
        if(e1.gt.rm1) then
          nsol=nsol+1
          itab(2)=1
        end if
        if(nsol.eq.0) then
          write(6,*) ' pas de solution pour t1 '
          stop
        end if

        if(nsol.eq.1) then
          do j=1,2
            if(itab(j).ne.0) j0=j
          end do
          t1=t1t(j0)
c          write(6,*) ' j0 t1',j0,t1
        end if

        if(nsol.eq.2) then
          do isol=1,2
            e1=t1t(isol)+rm1
            p12=e1*e1-rm1**2
            p1=sqrt(p12)
            aa1=a+b*t1t(isol)
            bb1=-2.d0*pa*p1*cs
            dift(isol)=abs(aa1-bb1)
          end do
          aa1=dift(1)
          aa2=dift(2)
c          write(6,*) ' dift(1)=',dift(1)
c          write(6,*) ' dift(2)=',dift(2)
          if(aa1.le.aa2) then
            t1=t1t(1)
	  else
            t1=t1t(2)
          end if
c          write(6,*) ' - 2 sol:t1=',t1
        end if
      end if
        
      e1=rm1+t1
      qp1(0)=e1
      p1=sqrt(t1*(t1+2.d0*rm1))
      p1x=p1*sin(teta)*cos(phi)
      p1y=p1*sin(teta)*sin(phi)
      p1z=p1*cs
      qp1(1)=p1x
      qp1(2)=p1y
      qp1(3)=p1z
      
      p2x=-p1x
      p2y=-p1y
      p2z=paz-p1z
      qp2(1)=p2x
      qp2(2)=p2y
      qp2(3)=p2z
      
      p22=p2x**2 + p2y**2 + p2z**2
      p2=sqrt(p22)
      
      e22=rm2**2 + p22
      e2=sqrt(e22)
      qp2(0)=e2
      t2=e2-rm2
      ef=e1+e2
      edif=ef-ei
      call SPHER(p2X,p2Y,p2Z,pp2,TETA2,PHI2)

      do i=0,3
        qdif=qpa(i)+qpb(i)-qp1(i)-qp2(i)
        aqdif=abs(qdif)
        if(aqdif.gt.1.d-12) then
          write(6,*) ' i qdif ',i,qdif
        end if
      end do
      

      
      do i=0,3
        qq(i)=qpa(i)+qpb(i)
      end do
      q2=qq(0)**2 - ( qq(1)**2 + qq(2)**2 + qq(3)**2 )


      do i=0,3
        qpa_lab(i)=qpa(i)
        qpb_lab(i)=qpb(i)
        qp1_lab(i)=qp1(i)
        qp2_lab(i)=qp2(i)
      end do
      
      return
      end
***************************************************************************
***************************************************************************
      SUBROUTINE SPHER(X,Y,Z,R,TETA,PHI)
      IMPLICIT REAL*8(A-H,O-Z)
      PI=acos(-1.d0)
      R2=X*X+Y*Y+Z*Z
      R=DSQRT(R2)
      if(R.lt.1.d-37) then
        teta=0.d0
        phi=0.d0
        return
      end if

C
      TETA=DACOS(Z/R)
C
      ATETA=ABS(TETA)
      PITETA=ABS(PI-TETA)
      IF(ATETA.LT.1.D-36.OR.PITETA.LT.1.D-36) THEN
      PHI=0.D0
      GO TO 500
      ENDIF
      AX=ABS(X)
      IF(AX.LT.1.D-36) THEN
             IF(Y.LT.0.D0) THEN
             PHI=1.5D0*PI
             ELSE
             PHI=PI/2.D0
             ENDIF
      GO TO 500
      ENDIF
C
      AY=ABS(Y)
      IF(AY.LT.1.D-36) Y=1.D-36
      IF(X.GT.0.D0.AND.Y.GT.0.D0)THEN
           PHI=DATAN(Y/X)
           ELSE IF(X.GT.0.D0.AND.Y.LT.0.D0) THEN
           PHI=DATAN(Y/X)+2.D0*PI
           ELSE IF(X.LT.0.D0.AND.Y.GT.0.D0) THEN
           PHI=DATAN(Y/X)+PI
           ELSE
           PHI=DATAN(Y/X)+PI
      ENDIF
 500  CONTINUE
      RETURN
      END
***************************************************************************
***************************************************************************
      subroutine blorentz(qv,bboost,qvp)
      implicit real*8(a-h,o-z)
      dimension qv(0:3),qvp(0:3),bboost(3)

      beta2=bboost(1)**2 + bboost(2)**2 + bboost(3)**2
      gamma=sqrt(1.d0/(1.d0-beta2))
      bscal=bboost(1)*qv(1) + bboost(2)*qv(2) + bboost(3)*qv(3)
      t=bscal*gamma/(gamma+1.d0) - qv(0)
      qvp(0)= gamma * ( qv(0) - bscal )
      do i=1,3
        qvp(i)=qv(i) + bboost(i) * gamma * t
      end do
      return
      end
***************************************************************************
***************************************************************************



















      
