!    According to Wikipedia formalism.
!  This code solves Rayleight-Plesset eq. using 4th order Runge-Kutta
!  method.
!
!   Define dR/dt=Z, then we solve 2nd order ODE by:
!   dR/dt=Z
!   dZ/dt=f . With f=-[1.5*Z^2+4v*Z/R+2*gamm/rho/R+DeltaP(t)/rho]/R .
! where DeltaP(t) is given in functional form, all other
! constants:v,gamm,rho are known (value will be read in).
      implicit none
      integer i,j,k,n
      double precision gamm,v,rho,w  ! constants to be read in
      double precision DeltaP,f    ! functions   
      double precision h,Ra,Rb,Za,Zb,t,Z0,R0,tt,rr  ! our variables    
      double precision j1,j2,j3,j4,k1,k2,k3,k4 ! RK intermediate quantities
      integer, parameter :: funt=23,funt0=24,funt1=25,funt2=26,funt3=27
      character (len=110) :: fname, fname5, fname1,fname2,fname3,fname4&
     &,fname6,fname7,fname8,fname9,fname10,fname11,fname12,fname13,&
     &fname14,fname15,fname16,fname17,fname18,fname19,fname20,fname21

!! First, read in gamm,v,rho, h from fort.10 , h is the time step size, tt is
!! the total time we want to evaluate.      
      read(10,*) gamm,v,rho,h,tt
!! also read in initial values, i.e., R(t=0), and dR/dt(t=0).      
      read(10,*) R0,Z0
!! read in fort.99 the oscillator frequency w, which is used in function
!! DeltaP at the end of the code.     
      read(99,*) w
      print*,'gamm,v,rho,time-step',gamm,v,rho,h
      n=int(tt/h)
      print*,'number of time steps',n
      print*,'initial R, dR/dt=',R0,Z0

      Za=Z0
      Ra=R0
      t=0.d0
      write(100,*) t*1000,Ra/R0,Ra,Za

      write(fname,'(i4)') int(w)
       fname = adjustl(fname)
      write(fname1,'(i3)') int(gamm)
       fname1 = adjustl(fname1)
      write(fname2,'(i3)') int(v)
       fname2 = adjustl(fname2)
      write(fname3,'(i3)') int(R0)
       fname3 = adjustl(fname3)
       rr=(R0-int(R0))*100.
      write(fname6,'(i3)') int(rr)
       fname6 = adjustl(fname6)
 
      fname4='bubble_Ro'//trim(fname3)//'.'//trim(fname6)//'_w_&
     &'//trim(fname)//'_gamm_'//trim(fname1)//'_v_'//trim(fname2)
      open (funt, file=fname4, status='replace')
      fname5='wave_w_'//trim(fname)
      open (funt1, file=fname5, status='replace')


      do i=1,n+1
              
      j1=h*Za
      j2=h*(Za+j1*0.5d0)
      j3=h*(za+j2*0.5d0)
      j4=h*(za+j3)
      Rb=Ra+j1/6.d0+j2/3.d0+j3/3.d0+j4/6.d0

      k1=h*f(t,Za,Ra,gamm,v,rho,R0,w)
      k2=h*f(t+0.5d0*h,Za+k1*0.5d0,Ra+j1*0.5d0,gamm,v,rho,R0,w)
      k3=h*f(t+0.5d0*h,Za+k2*0.5d0,Ra+j2*0.5d0,gamm,v,rho,R0,w)
      k4=h*f(t+h,Za+k3,Ra+j3,gamm,v,rho,R0,w)
      Zb=Za+k1/6.d0+k2/3.d0+k3/3.d0+k4/6.d0

      t=t+h
!!! For R<0, the bubble will blow up      
      if(Rb.lt.0.d0) then
      print*,'blow up at t=',t
      goto 100
      endif

      write(funt,*) t*1000,Rb/R0,Rb,Zb
      write(funt1,*) t*1000,DeltaP(t,w)
      Ra=Rb
      Za=Zb
      enddo

100      end

!!!!!!!
      double precision function f(ti,Z,R,gamm,v,rho,R0,w)     
      implicit none
      double precision ti,Z,R,gamm,v,rho,DeltaP,pgo,ka,R0,w
      ka=1.d0  ! ka is 1 for isothermal case
      pgo=101.!101. ! pgo is the inner pressure inside bubble at t=0 
!      f=-(1.5d0*Z*Z+4.d0*v*Z/R+2.d0*gamm/rho/R+DeltaP(ti)/rho)/R !for Wiki
      f=-(1.5d0*Z*Z+4.d0*v*Z/R+2.d0*gamm/rho/R+(DeltaP(ti,w)-&
     &pgo*(R0/R)**(3*ka))/rho)/R !for hal-00265882.
      return
      end 

!!!!!!!!
      double precision function DeltaP(ta,w)
      implicit none
      double precision ta,pi,pinf,w
!      pinf=101.3d0 ! pinf is the external pressure infinitely far from the bubble

      pi=3.1415926535897932384626d0
!      DeltaP=(101.-50.5*sin(2.*pi*55.1138*ta))-102.36
      DeltaP=(101.-200.*sin(2.*pi*w*ta))
      return
      end      
