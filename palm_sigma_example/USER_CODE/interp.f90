!> @file interp.f90
!--------------------------------------------------------------------------------------------------!
! This file is part of PALM-Sigma.
!
! PALM-Sigma consists of a series of user-defined modules and modified PALM modules, developed for
! wave-phase-resolved LES. PALM-Sigma should be used as USER_CODE based on the framework of PALM.
!
! PALM-Sigma is free software: you can redistribute it and/or modify it under the terms of the GNU
! General Public License as published by the Free Software Foundation, either version 3 of the
! License, or (at your option) any later version.
!
! PALM-Sigma is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
! even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
! General Public License for more details.
!
! Copyright 2025 Xu Ning
!--------------------------------------------------------------------------------------------------!
!
!
! Description:
! ------------
!> A user-defined module for linear interpolation.
!--------------------------------------------------------------------------------------------------!

MODULE interp

IMPLICIT NONE

 CONTAINS

 INTEGER FUNCTION locate(x,x0)

 ! Locate a value in low-to-high sorted array in a cyclic way (if x0 is out of the range of x, it will be add or substract the
 ! integer times of (x(n)-x(1)) until it is within the range).
 ! x should be a low-to-high sorted array
 ! x0 can be any real number

 REAL(8), INTENT(IN) :: x(:), x0
 INTEGER :: n,jl,jm,jr

 n = SIZE(x)

 jl = 1
 jr = n+1

 DO
    IF (jr - jl <= 1) exit
    jm = (jr + jl)/2
    IF (x0 > x(jm)) THEN
       jl = jm
    ELSE
       jr = jm
    END IF
 END DO

 locate = jl

 IF ( x0 < x(1) .OR. x0 > x(n) ) THEN
    WRITE(*,*) 'Err: out of range'
    locate = -1
 ENDIF

 !IF (x0 == x(1)) THEN
 !   locate = 1
 !ELSE IF (x0 == x(n)) THEN
 !   locate = n-1
 !ELSE IF (x0 > x(n) .OR. x0 < x(1)) THEN
 !   locate = -1
 !ELSE
 !   locate = jl
 !END IF

 END FUNCTION locate

 REAL FUNCTION interp3d(x,y,z,v,x0,y0,z0)

     REAL(8), INTENT(IN) :: x(:), y(:), z(:), v(:,:,:), x0, y0, z0

     REAL(8) :: x0_, y0_, z0_

     INTEGER :: i1, i2, j1, j2, k1, k2, nx, ny, nz

     REAL :: xd, yd, zd, v000, v100, v001, v101, v010, v110, v011, v111, v00, v01, v10, v11, v0, v1

     IF(size(x).NE.size(v,1)) WRITE(*,*) "size of x does not match array of v"
     IF(size(y).NE.size(v,2)) WRITE(*,*) "size of y does not match array of v"
     IF(size(z).NE.size(v,3)) WRITE(*,*) "size of z does not match array of v"

     x0_ = x0 ! in order not to change x0
     y0_ = y0 ! in order not to change y0
     z0_ = z0 ! in order not to change z0
     nx = SIZE(x)
     ny = SIZE(y)
     nz = SIZE(z)

     IF ( x0_ < x(1) ) THEN
        DO
           x0_ = x0_ + ( x(nx) - x(1) )
           IF ( x0_ >= x(1) ) EXIT
        ENDDO
     ELSEIF ( x0_ > x(nx) ) THEN
        DO
           x0_ = x0_ - ( x(nx) - x(1) )
           IF ( x0_ <= x(nx) ) EXIT
        ENDDO
     ENDIF

     IF ( y0_ < y(1) ) THEN
        DO
           y0_ = y0_ + ( y(ny) - y(1) )
           IF ( y0_ >= y(1) ) EXIT
        ENDDO
     ELSEIF ( y0_ > y(ny) ) THEN
        DO
           y0_ = y0_ - ( y(ny) - y(1) )
           IF ( y0_ <= y(ny) ) EXIT
        ENDDO
     ENDIF

     IF ( z0_ < z(1) ) THEN
        DO
           z0_ = z0_ + ( z(nz) - z(1) )
           IF ( z0_ >= z(1) ) EXIT
        ENDDO
     ELSEIF ( z0_ > z(nz) ) THEN
        DO
           z0_ = z0_ - ( z(nz) - z(1) )
           IF ( z0_ <= z(nz) ) EXIT
        ENDDO
     ENDIF

     i1 = locate(x,x0_)
     IF ( i1 == SIZE(x) ) THEN
        i2 = 1
     ELSE
        i2 = i1 + 1
     ENDIF

     j1 = locate(y,y0_)
     IF ( j1 == SIZE(y) ) THEN
        j2 = 1
     ELSE
        j2 = j1 +1
     ENDIF

     k1 = locate(z,z0_)
     IF ( k1 == SIZE(z) ) THEN
        k2 = 1
     ELSE
        k2 = k1 + 1
     ENDIF

     !IF(i1==-1) THEN
     !    WRITE(*,*) 'ERROR: x0 out of bound.'
     !END IF

     !IF(j1==-1) THEN
     !    WRITE(*,*) 'ERROR: y0 out of bound.'
     !END IF

     !IF(k1==-1) THEN
     !    WRITE(*,*) 'ERROR: z0 out of bound.'
     !END IF

     xd = (x0_ - x(i1)) / (x(i2) - x(i1))
     yd = (y0_ - y(j1)) / (y(j2) - y(j1))
     zd = (z0_ - z(k1)) / (z(k2) - z(k1))

     v000 = v(i1,j1,k1)
     v100 = v(i2,j1,k1)
     v001 = v(i1,j1,k2)
     v101 = v(i2,j1,k2)
     v010 = v(i1,j2,k1)
     v110 = v(i2,j2,k1)
     v011 = v(i1,j2,k2)
     v111 = v(i2,j2,k2)

     v00 = v000*(1 - xd) + v100*xd
     v01 = v001*(1 - xd) + v101*xd
     v10 = v010*(1 - xd) + v110*xd
     v11 = v011*(1 - xd) + v111*xd

     v0 = v00*(1 - yd) + v10*yd
     v1 = v01*(1 - yd) + v11*yd

     interp3d = v0*(1 - zd) + v1*zd

 END FUNCTION interp3d


END MODULE interp
