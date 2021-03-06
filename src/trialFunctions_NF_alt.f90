program SincInterpolation !( nx,nz,rho,lam,mu,dx,dz,dt )
  implicit none
  integer i
  integer ix,iz
  integer jx,jz
  integer imx,imz,inx,inz ! 'small' coordinates inside M and N 
  integer mx,mz,nx,nz ! for trial functions
  integer mxmax,mxmin,mzmax,mzmin ! the max and min values of mx,mz,nx,nz 
  double precision xm,zm
  integer ndis,ngrid,npTF
  logical, parameter :: sincfunction = .false.
  double precision x,xx,z,zz
  !double precision rho(nx+1,nz+1),lam(nx+1,nz+1),mu(nx+1,nz+1)
  double precision, allocatable :: lam(:,:),mu(:,:),rho(:,:)
  double precision dx,dz,dt,smalldx,smalldz
  double precision dx2,dz2,dxdz,dt2
  double precision, allocatable :: phix(:),phiz(:) ! non-zero values for phix, phiz
  double precision, allocatable :: phixderiv(:),phizderiv(:) ! and theirs derivatives
  double precision, dimension (:,:,:,:), allocatable :: T0,H11,H13,H31,H33
  double precision, parameter :: pi = 3.141592653589793238462643383
  double precision, parameter :: eps = 1.d-20 ! if the denominator is smaller than eps, we do not perform any division
 

 ! dt = 1.d0
  dx = 1.d0
  dz = 1.d0

  
  !npTF defines points in scheme (3, 5, 7)
  npTF = 3
  ngrid = (npTF-1)/2
  ndis = 100


  smalldx = dx/dble(ndis)
  smalldz = dz/dble(ndis)


  mxmin = 1
  mxmax = npTF

  mzmin = 1
  mzmax = npTF



  allocate(phix(-ngrid*ndis:ngrid*ndis))
  allocate(phiz(-ngrid*ndis:ngrid*ndis))
  allocate(phixderiv(-ngrid*ndis:ngrid*ndis))
  allocate(phizderiv(-ngrid*ndis:ngrid*ndis))
  allocate(lam(-ngrid*ndis:ngrid*ndis,-ngrid*ndis:ngrid*ndis))
  


  allocate(T0(2:2,2:2,mxmin:mxmax,mzmin:mzmax))
  allocate(H11(2:2,2:2,mxmin:mxmax,mzmin:mzmax))
  allocate(H13(2:2,2:2,mxmin:mxmax,mzmin:mzmax))
  allocate(H31(2:2,2:2,mxmin:mxmax,mzmin:mzmax))
  allocate(H33(2:2,2:2,mxmin:mxmax,mzmin:mzmax))

 ! lam =1.d0

 ! dt2 = dt * dt
 ! dx2 = dx * dx
 ! dz2 = dz * dz
 ! dxdz = dx * dz
  
  xm = 0.d0
  zm = 0.d0


  ! I put mx, mz to be the centre (2,2)
  mx = 2
  mz = 2
  

  ! then mx, mz have mxmin:mxmax and mzmin:mzmax


  ! Initialising vectors

  phix = 0.d0
  phiz = 0.d0
  
  phixderiv = 0.d0
  phizderiv = 0.d0



  !trialfunction decides on sinc(true) or linear(false) interpolation
  !sincfunction = .true.



  ! Trial function calculations for sinc/spline
  
  if (sincfunction) then
     
     do i=-ngrid*ndis, ngrid*ndis

        xx =dble(i/ndis)*dx
        if(abs(xx)<eps) then
           phix(i) = 1.d0
           phixderiv(i) = 0.d0
        elseif(abs(xx)>eps) then
           phix(i) = sin(pi*xx)/(pi*xx)
           phixderiv(i) = pi*cos(pi*xx)/(pi*xx) - sin(pi*xx)/(pi*xx*xx)
        endif
        phiz(i) = phix(i)
        phizderiv(i) =phixderiv(i)
        

        
          

       ! open(unit=8,file="phix.dat",form="formatted"&
       !      ,status="replace",action="write")
        
       ! write(8,*)'phix', phix(ix,iz)f
        
       ! close(8)
        
        
     enddo
     
     
     
  else
     

     ! B-splines (for the moment only for 3 points so it's not correct for 5-, 7- point schemes)

     do i=-ngrid*ndis,0
        
        x = xm + dble(i/ndis)*dx
        phix(i+1+ndis) = (x+dx)/dx
        phixderiv(i+1+ndis) = 1.d0/dx

        z = zm + dble(i/ndis)*dz
        phiz(i+1+ndis) = (z+dz)/dz
        phizderiv(i+1+ndis) = 1.d0/dz
        
        
     
     enddo
     
     do i=0,ngrid*ndis
      
        x = xm + dble(i/ndis)*dx
        phix(i+1+ndis) = (-x+dx)/dx
        phixderiv(i+1+ndis) = -1.d0/dx

        z = z+ dble(i/ndis)*dz
        phiz(i+1+ndis) = (-z+dz)/dz
        phizderiv(i+1+ndis) = -1.d0/dz
        
        
     enddo


   
  endif

  


  
  ! writing phix and phixderiv in 'fort.12'

  do i = -ngrid*ndis,ngrid*ndis
     write(12,*) dble(i/ndis)*dx,phix(i),phixderiv(i)
  enddo
  
  ! end

  mx = 2
  mz = 2

  
  T0(mx,mz,nx,mz) = 0.d0
  H11(mx,mz,nx,mz) = 0.d0
  H13(mx,mz,nx,mz) = 0.d0
  H31(mx,mz,nx,mz) = 0.d0
  H33(mx,mz,nx,mz) = 0.d0



  do nx = mxmin,mxmax
     do nz = mzmin,mzmax
        

        
        do inx = -ngrid*ndis,ngrid*ndis
           

           imx = inx-(mx-nx)*ndis

           
           if((imx-(-ngrid*ndis))*(imx-(ngrid*ndis-1))<=0) then
              ! the integrand of two phix functions has non-zero value if-and-only-if the 'small' coordinates for M and N (imx, inx) 
              ! have values inside -ngrid*ndis:ngrid:ndis (otherwise phix is not defined)



              do inz = -ngrid*ndis, ngrid*ndis
                 imz = inz-(mz-nz)*ndis
                 
                 if((imz-(-ngrid*ndis))*(imz-(ngrid*ndis-1))<=0) then
                     ! same story for phiz

!                    T0(mx,mz,nx,nz) = T0(mx,mz,nx,nz)  &
!                         + 5.d-1*(phix(inx)+phix(inx+1))*5.d-1*(phix(imx)+phix(imx+1))   &
!                         *5.d-1*(phiz(inz)+phiz(inz+1))*5.d-1*(phiz(imz)+phiz(imz+1)) &
!                         *smalldx*smalldz

                    ! NF : T0 should be doubled if the same phix and phiz are used for x and z component trial functions

!                    H11(mx,mz,nx,nz) = H11(mx,mz,nx,nz) + 5.d-1*(phixderiv(inx)+phixderiv(inx+1)) * &
!                         5.d-1*(phixderiv(imx)+phixderiv(imx+1))*5.d-1*(phiz(inz)+phiz(inz+1))* &
!                         5.d-1*(phiz(imz)+phiz(imz+1))*smalldx*smalldz

!                    H13(mx,mz,nx,nz) = H13(mx,mz,nx,nz) + 5.d-1*(phizderiv(inz)+phizderiv(inz+1))* &
!                         5.d-1*(phixderiv(imx)+phixderiv(imx+1))*5.d-1*(phix(inx)+phix(inx+1))* &
!                         5.d-1*(phiz(imz)+phiz(imz+1))*smalldx*smalldz

!                    H31(mx,mz,nx,nz) = H31(mx,mz,nx,nz) + 5.d-1*(phixderiv(inx)+phixderiv(inx+1))* &
!                         5.d-1*(phizderiv(imz)+phizderiv(imz+1))*5.d-1*(phiz(inz)+phiz(inz+1))* &
!                         5.d-1*(phix(imx)+phix(imx+1))*smalldx*smalldz

!                    H33(mx,mz,nx,nz) = H33(mx,mz,nx,nz) + 5.d-1*(phizderiv(inz)+phizderiv(inx+1))* &
!                         5.d-1*(phizderiv(imz)+phizderiv(imz+1))*5.d-1*(phix(inx)+phix(inx+1)) &
!                         *5.d-1*(phix(imx)+phix(imx+1))*smalldx*smalldz



                    T0(mx,mz,nx,nz) = T0(mx,mz,nx,nz)+phix(inx+1+ngrid*ndis)*phix(imx+1+ngrid*ndis)&
                         *phiz(inz+1+ngrid*ndis)*phiz(imz+1+ngrid*ndis)*smalldx*smalldz


                    H11(mx,mz,nx,nz) = H11(mx,mz,nx,nz)+phixderiv(inx+1+ngrid*ndis)*phixderiv(imx+1+ngrid*ndis)&
                         *phiz(inz+1+ngrid*ndis)*phiz(imz+1+ngrid*ndis)*smalldx*smalldz 

                    H13(mx,mz,nx,nz) = H13(mx,mz,nx,nz)+phizderiv(inz+1+ngrid*ndis)*phixderiv(imx+1+ngrid*ndis)&
                         *phix(inx+1+ngrid*ndis)*phiz(imz+1+ngrid*ndis)*smalldx*smalldz

                    H31(mx,mz,nx,nz) = H31(mx,mz,nx,nz)+phixderiv(inx+1+ngrid*ndis)*phizderiv(imz+1+ngrid*ndis)&
                         *phiz(inz+1+ngrid*ndis)*phix(imx+1+ngrid*ndis)*smalldx*smalldz

                    H33(mx,mz,nx,nz) = H33(mx,mz,nx,nz)+phizderiv(inz+1+ngrid*ndis)*phizderiv(imz+1+ngrid*ndis)&
                         *phix(inx+1+ngrid*ndis)*phix(imx+1+ngrid*ndis)*smalldx*smalldz

!                    print *,'mx,mz,nx,nz', mx,mz,nx,nz,'T0', T0(mx,mz,nx,nz)
!                    print *,'mx,mz,nx,nz', mx,mz,nx,nz,'H11', H11(mx,mz,nx,nz), 'H13', H13(mx,mz,nx,nz), 'H31', H31(mx,mz,nx,nz)&
!                           , 'H33', H33(mx,mz,nx,nz)
                    
                 endif
              enddo
           endif

           
        enddo



        
   
      

     
     enddo


  enddo
  



  open(unit=8,file="coeffs.dat",form="formatted"&
       ,status="replace",action="write")
  
  
  mx = 2
  mz = 2

  do nx = mxmin,mxmax
     do nz = mzmin,mzmax
        write(8,*) 'mx,mz,nx,nz',mx,mz,nx,nz 
        write(8,*)'T0', T0(mx,mz,nx,nz)
        write(8,*)'H11', H11(mx,mz,nx,nz)
        write(8,*)'H13', H13(mx,mz,nx,nz)
        write(8,*)'H31', H31(mx,mz,nx,nz)
        write(8,*)'H33', H33(mx,mz,nx,nz)
        
     enddo
  enddo

     

  close(8)



end program SincInterpolation
