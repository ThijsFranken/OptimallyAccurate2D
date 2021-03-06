


subroutine cales( nx,nz,rho,lam,mu,dt,dx,dz,e1, e2, e3, e4, e5, e6, e7, e8,&
     e13,e14,e15,e16,e17,e18,e19,e20, &
     f1, f2, f3, f4, f5, f6, f7, f8, &
     f13,f14,f15,f16,f17,f18,f19,f20 )
  implicit none
  integer nx,nz
  double precision rho(nx+1,nz+1),lam(nx+1,nz+1),mu(nx+1,nz+1)
  double precision dt,dx,dz
  double precision  e1(nx+1,nz+1), e2(nx+1,nz+1), e3(nx+1,nz+1)
  double precision  e4(nx+1,nz+1), e5(nx+1,nx+1), e6(nx+1,nz+1)
  double precision  e7(nx+1,nz+1), e8(nx+1,nz+1)
  double precision e13(nx+1,nz+1),e14(nx+1,nz+1),e15(nx+1,nz+1)
  double precision e16(nx+1,nz+1),e17(nx+1,nz+1),e18(nx+1,nz+1)
  double precision e19(nx+1,nz+1),e20(nx+1,nz+1)
  double precision  f1(nx+1,nz+1), f2(nx+1,nz+1), f3(nx+1,nz+1)
  double precision  f4(nx+1,nz+1), f5(nx+1,nz+1), f6(nx+1,nz+1)
  double precision  f7(nx+1,nz+1), f8(nx+1,nz+1)
  double precision f13(nx+1,nz+1),f14(nx+1,nz+1),f15(nx+1,nz+1)
  double precision f16(nx+1,nz+1),f17(nx+1,nz+1),f18(nx+1,nz+1)
  double precision f19(nx+1,nz+1),f20(nx+1,nz+1)
  integer ix,iz
  double precision dt2,dx2,dz2,dxdz

  ! for smoothed part of the model :
  !  we use Zahradnik operators and optimally accurate operators

  
  dt2 = dt * dt
  dx2 = dx * dx
  dz2 = dz * dz
  dxdz = dx * dz
  
  do iz=2,nz
     do ix=2,nx
        e1(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix-1,iz) + lam(ix,iz) ) &
             + 2.d0 * ( mu(ix-1,iz) + mu(ix,iz) ) ) &
             / ( 2.d0 * dx2 )
        e2(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix,iz) + lam(ix+1,iz) ) &
             + 2.d0 * ( mu(ix,iz) + mu(ix+1,iz) ) ) &
             / ( 2.d0 * dx2 )
        e3(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix,iz-1) + mu(ix,iz) ) &
             / ( 2.d0 * dz2 )
        e4(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix,iz) + mu(ix,iz+1) ) &
             / ( 2.d0 * dz2 )
        e5(ix,iz) = dt2 / rho(ix,iz) * lam(ix-1,iz) &
             / ( 4.d0 * dxdz )
        e6(ix,iz) = dt2 / rho(ix,iz) * lam(ix+1,iz) &
             / ( 4.d0 * dxdz )
        e7(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz-1) &
             / ( 4.d0 * dxdz )
        e8(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+1) &
             / ( 4.d0 * dxdz )
        e13(ix,iz) = dt2 / rho(ix,iz) * lam(ix-1,iz) &
             * ( -5.d0 ) / ( 1728.d0 * dxdz )
        e14(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz) &
             * ( -3.d0 ) / ( 1728.d0 * dxdz )
        e15(ix,iz) = dt2 / rho(ix,iz) * lam(ix+1,iz) &
           * ( +9.d0 ) / ( 1728.d0 * dxdz )
        if ( ix+2.le.nx+1 ) then
           e16(ix,iz) = dt2 / rho(ix,iz) * lam(ix+2,iz) &
                * ( -1.d0 ) / ( 1728.d0 * dxdz )
        else
           e16(ix,iz) = 0.d0
        endif
        e17(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz-1) &
             * ( -5.d0 ) / ( 1728.d0 * dxdz )
        e18(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz) &
             * ( -3.d0 ) / ( 1728.d0 * dxdz )
        e19(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+1) &
             * ( +9.d0 ) / ( 1728.d0 * dxdz )
        if ( iz+2.le.nz+1) then
           e20(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz+2) &
                * ( -1.d0 ) / ( 1728.d0 * dxdz )
        else
           e20(ix,iz) = 0.d0
        endif
        f1(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix-1,iz) + mu(ix,iz) ) &
             / ( 2.d0 * dx2 )
        f2(ix,iz) = dt2 / rho(ix,iz) &
             * ( mu(ix,iz) + mu(ix+1,iz) ) &
             / ( 2.d0 * dx2 )
        f3(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix,iz-1) + lam(ix,iz) ) &
             + 2.d0 * ( mu(ix,iz-1) + mu(ix,iz) ) ) &
             / ( 2.d0 * dz2 )
        f4(ix,iz) = dt2 / rho(ix,iz) &
             * ( ( lam(ix,iz) + lam(ix,iz+1) ) &
           + 2.d0 * ( mu(ix,iz) + mu(ix,iz+1) ) ) &
           / ( 2.d0 * dz2 )
        f5(ix,iz) = dt2 / rho(ix,iz) * mu(ix-1,iz) &
             / ( 4.d0 * dxdz )
        f6(ix,iz) = dt2 / rho(ix,iz) * mu(ix+1,iz) &
             / ( 4.d0 * dxdz )
        f7(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-1) &
             / ( 4.d0 * dxdz )
        f8(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz+1) &
             / ( 4.d0 * dxdz )
        if ( ix-2.ge.1 ) then
           f13(ix,iz) = dt2 / rho(ix,iz) * mu(ix-2,iz) &
                * (  1.d0 ) / ( 1728.d0 * dxdz )
        else
           f13(ix,iz) = 0.d0
        endif
        f14(ix,iz) = dt2 / rho(ix,iz) * mu(ix-1,iz) &
             * ( -9.d0 ) / ( 1728.d0 * dxdz )
        f15(ix,iz) = dt2 / rho(ix,iz) * mu(ix,iz) &
             * (  3.d0 ) / ( 1728.d0 * dxdz )
        f16(ix,iz) = dt2 / rho(ix,iz) * mu(ix+1,iz) &
             * (  5.d0 ) / ( 1728.d0 * dxdz )
        if ( iz-2.ge.1 ) then
           f17(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-2) &
                * (  1.d0 ) / ( 1728.d0 * dxdz )
        else
           f17(ix,iz) = 0.d0
        endif
        f18(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz-1) &
             * ( -9.d0 ) / ( 1728.d0 * dxdz )
        f19(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz) &
             * (  3.d0 ) / ( 1728.d0 * dxdz )
        f20(ix,iz) = dt2 / rho(ix,iz) * lam(ix,iz+1) &
             * (  5.d0 ) / ( 1728.d0 * dxdz )

        
     enddo
  enddo
  
  open(8,file='MultiOperators.out',status='new')

  ix=nx
  iz=nz
  
  write (8,*)'e1', e1(ix,iz)
  write (8,*)'e2', e2(ix,iz)
  write (8,*)'e3', e3(ix,iz)
  write (8,*)'e4', e4(ix,iz)
  write (8,*)'e5', e5(ix,iz)
  write (8,*)'e6', e6(ix,iz)
  write (8,*)'e7', e7(ix,iz)
  write (8,*)'e8', e8(ix,iz)
  write (8,*)'e13', e13(ix,iz)
  write (8,*)'e14', e14(ix,iz)
  write (8,*)'e15', e15(ix,iz)
  write (8,*)'e16', e16(ix,iz)
  write (8,*)'e17', e17(ix,iz)
  write (8,*)'e18', e18(ix,iz)
  write (8,*)'e19', e19(ix,iz)
  write (8,*)'e20', e20(ix,iz)

  write (8,*)'f1', f1(ix,iz)
  write (8,*)'f2', f2(ix,iz)
  write (8,*)'f3', f3(ix,iz)
  write (8,*)'f4', f4(ix,iz)
  write (8,*)'f5', f5(ix,iz)
  write (8,*)'f6', f6(ix,iz)
  write (8,*)'f7', f7(ix,iz)
  write (8,*)'f8', f8(ix,iz)
  write (8,*)'f13', f13(ix,iz)
  write (8,*)'f14', f14(ix,iz)
  write (8,*)'f15', f15(ix,iz)
  write (8,*)'f16', f16(ix,iz)
  write (8,*)'f17', f17(ix,iz)
  write (8,*)'f18', f18(ix,iz)
  write (8,*)'f19', f19(ix,iz)
  write (8,*)'f20', f20(ix,iz)
  
  close(8,status='keep')

  stop
  return
end subroutine cales
