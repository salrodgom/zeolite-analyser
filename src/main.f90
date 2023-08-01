! ZEOLYSIS code:
! --------------
! By Salvador R.G. Balestra, 2020,1) A. R. Ruiz-Salvador, 1), and Rocio Semino,2)
!    1) Departamento de Sistemas Físicos, Químicos y Naturales, Universidad Pablo de Olavide, Seville, Spain
!    2) Institute Charles Gerhardt Montpellier, Montpellier, France
!
! This code calculates angle distributions in zeolitic structures
! The code analyses CIF files. The CIF File is provided by input a list file
! with the path to the listed CIFFiles

! Some details:
! 1. The CIF file must have this format for coordinates:
!
!loop_
!_atom_site_label
!_atom_site_fract_x
!_atom_site_fract_y
!_atom_site_fract_z
!_atom_site_occupancy
!  Si          0.86398      0.51151      0.63128  1.0000
!  O           0.30362      0.86548      0.69030  1.0000
!
module functions
! used functions to manage strings in fortran
 implicit none
 private
 public clen_trim, clen, print_help
 contains
 pure integer function clen(s)       ! returns same result as len unless:
  character(*),intent(in) :: s       ! last non-blank char is null
  integer :: i
  clen = len(s)
  i = len_trim(s)
  if (s(i:i) == char(0)) clen = i-1  ! len of c string
  return
 end function clen
!
 pure integer function clen_trim(s)  ! returns same result as len_trim unless:
  character(*),intent(in) :: s       ! last char non-blank is null, if true:
  integer :: i                       ! then len of c string is returned, note:
                                     ! ctrim is only user of this function
  i = len_trim(s) ; clen_trim = i
  if (s(i:i) == char(0)) clen_trim = clen(s)   ! len of c string
 end function clen_trim
 !
 subroutine print_help()
  write(6,'(a)')'Usage: ./zeolite_analyser.exe --all'
 end subroutine
end module functions
!
module vector_module
 implicit none
 public
 type  :: vector
  sequence
  real :: x
  real :: y
  real :: z
 end type vector
 contains
!
 pure type(vector) function Array2Vector(a)
  implicit none
  real,intent(in) :: a(1:3)
  array2vector%x = a(1)
  array2vector%y = a(2)
  array2vector%z = a(3)
 end function Array2Vector
!
 pure function Vector2Array(v) result(a)
  implicit none
  type(vector),intent(in) :: v
  real :: a(1:3)
  a(1) = v%x
  a(2) = v%y
  a(3) = v%z
 end function Vector2Array
!
 pure type(vector) function PointCoordinatesAtDistanceFromOriginAndDirection(o,u,d)
  implicit none
  type(vector),intent(in) :: o,u
  real,intent(in)         :: d
  PointCoordinatesAtDistanceFromOriginAndDirection%x = o%x + u%x*d
  PointCoordinatesAtDistanceFromOriginAndDirection%y = o%y + u%y*d
  PointCoordinatesAtDistanceFromOriginAndDirection%z = o%z + u%z*d
 end function PointCoordinatesAtDistanceFromOriginAndDirection
!
 pure type (vector) function vector_unitary(v1)
  implicit none
  type(vector), intent(in) :: v1
  real                     :: norma
  norma=absvec(v1)
  vector_unitary%x = v1%x/norma
  vector_unitary%y = v1%y/norma
  vector_unitary%z = v1%z/norma
 end function vector_unitary
!
 pure type (vector) function vector_scale(v1, r )
  implicit none
  type(vector), intent(in) :: v1
  real,intent(in)          :: r
  vector_scale%x = v1%x*r
  vector_scale%y = v1%y*r
  vector_scale%z = v1%z*r
 end function vector_scale
!
 pure type (vector) function vector_add(v1, v2)
  implicit none
  type (vector), intent(in) :: v1, v2
  vector_add%x = v1%x + v2%x
  vector_add%y = v1%y + v2%y
  vector_add%z = v1%z + v2%z
 end function vector_add
!
 pure type (vector) function vector_sub(v1, v2)
  implicit none
  type (vector), intent(in) :: v1, v2
  vector_sub%x = v1%x - v2%x
  vector_sub%y = v1%y - v2%y
  vector_sub%z = v1%z - v2%z
 end function vector_sub
!
 pure type(vector) function cross(a, b)
  implicit none
  type(vector),intent(in) :: a, b
  cross%x = a%y * b%z - a%z * b%y
  cross%y = a%z * b%x - a%x * b%z
  cross%z = a%x * b%y - a%y * b%x
 end function cross
!
 pure real function dot(a, b)
  implicit none
  type(vector), intent(in) :: a, b
  dot = a%x * b%x + a%y*b%y + a%z*b%z
 end function dot
!
 pure real function absvec(a)
  type (vector), intent (in) :: a
  absvec = sqrt(a%x*a%x + a%y*a%y + a%z*a%z)
 end function absvec
!
 pure real function Dist2Vectors(a, b)
  type(vector), intent(in) :: a,b
  dist2vectors = absvec( vector_sub(b,a) )
 end function  Dist2Vectors
!
 pure real function angle2vectors(a, b)
  ! vectors are the  r_ij vector, (r_ij, r_jk)
  type(vector), intent(in) :: a, b
  angle2vectors = acos( dot(a,b)/(absvec(a)*absvec(b)) )
 end function  angle2vectors
!
 pure real function angle3vectors_jik(r_i, r_j, r_k)
  type(vector), intent(in) :: r_i, r_j, r_k
  type(vector)             :: r_ij, r_ik !, r_jk
  real                     :: x
  r_ij = vector_sub(r_i,r_j)
  !r_jk = vector_sub(r_j,r_k)
  r_ik = vector_sub(r_i,r_k)
  x = dot(r_ij,r_ik)/(absvec(r_ij)*absvec(r_ik))
  if (x > 1.0) then
   angle3vectors_jik = acos(1.0)
  elseif ( x < -1.0 ) then
   angle3vectors_jik = acos(-1.0)
  else
   angle3vectors_jik = acos(x)
  end if
  return
 end function  angle3vectors_jik
!
 pure real function dihedral_angle4vectors_ijkl(r_i, r_j, r_k, r_l)
  type(vector), intent(in) :: r_i, r_j, r_k, r_l
  type(vector)             :: r_ij, r_jk, r_kl
  real                     :: x
  r_ij = vector_sub(r_i,r_j)
  r_jk = vector_sub(r_j,r_k)
  r_kl = vector_sub(r_k,r_l)
  x = dot(cross(r_ij,r_jk),cross(r_jk,r_kl))/(absvec(cross(r_ij,r_jk))*absvec(cross(r_jk,r_kl)))
  if (x > 1.0) then
   dihedral_angle4vectors_ijkl = acos(1.0)
  elseif ( x < -1.0 ) then
   dihedral_angle4vectors_ijkl = acos(-1.0)
  else
   dihedral_angle4vectors_ijkl = acos(x)
  end if
  return
 end function  dihedral_angle4vectors_ijkl
!
 pure real function dihedral_angle4vectors_ijkm(r_i, r_j, r_k, r_m)
  type(vector), intent(in) :: r_i, r_j, r_k, r_m
  type(vector)             :: b_0,b_1,b_2
  real                     :: x
  b_0 = vector_sub(r_j,r_i)
  b_1 = vector_sub(r_k,r_i)
  b_2 = vector_sub(r_m,r_i)
  x = dot(cross(b_0,b_1),cross(b_0,b_2))/(absvec(cross(b_0,b_1))*absvec(cross(b_0,b_2)))
  if (x > 1.0) then
   dihedral_angle4vectors_ijkm = acos(1.0)
  elseif ( x < -1.0 ) then
   dihedral_angle4vectors_ijkm = acos(-1.0)
  else
   dihedral_angle4vectors_ijkm = acos(x)
  end if
  return
 end function  dihedral_angle4vectors_ijkm
end module vector_module
!
module types
 use vector_module
 implicit none
 integer, parameter :: max_atom_number = 1000
 integer, parameter :: max_n_componets = 2
 integer, parameter :: max_n_componets_in_Geometric_Thing = 100000
 ! Structure for atoms, molecules, ligands, or nodes.
 type :: particle
  integer           :: element
  character(len=2)  :: label_element
  character(len=4)  :: label
  character(len=6)  :: label_from_CIFFile="Xxxxxx"
  integer           :: degree
  integer           :: N_SiSi
  real              :: Sum_SiSi_d
  real              :: Sum_SiOSi_angle
  real              :: Sum_SiOSi_dot_SiO_product
  real              :: Sum_rho
  real              :: charge
  real              :: Q
  real              :: radius
  real              :: mass
  type(vector)      :: uCoordinate ! x,y,z
  type(vector)      :: xCoordinate ! x,y,z
 end type
 type :: CollectiveVariable
  real               :: k1,k2,k3,k4
  real               :: max_,min_
 end type
 ! Structure for Structure Files
 type :: CIFFile
  character(len=100) :: filename
  integer            :: n_atoms
  real               :: rv(1:3,1:3)
  real               :: vr(1:3,1:3)
  real               :: cell_0(1:6)
  real               :: energy
  type(particle), allocatable, dimension(:) :: atom
 end type
end module types
!
module GeneralVariables
 use types
 type(CIFfile),allocatable                    :: CIFFiles(:)
 integer                                      :: n_files = 0
 logical                                      :: generate_list_file = .false.
 logical                                      :: preoptimization_flag = .false.
 logical                                      :: opti_flag = .true.
 logical                                      :: shellonly_flag = .false.
 logical                                      :: restart_flag = .false.
 real                                         :: infinite = 3.4028e38
 real, parameter                              :: k_B = 8.617332478e-5
 real, parameter                              :: r_min_criteria_connectivity = 0.15
end module GeneralVariables
!
module histograms
 implicit none
 contains
 subroutine histogram(data,n_datos,n_boxs,ave,adev,sdev,var,skew,curt,suma, max_,min_)
   implicit none
   INTEGER,intent(in) :: n_datos,n_boxs
   REAL,   intent(in) :: data(1:n_datos)
   REAL               :: bound(0:n_boxs), histo(0:n_boxs)
   REAL,   intent(out):: ave,adev,sdev,var,skew,curt,suma, max_,min_
   INTEGER            :: i
   suma = 0.0
   max_ = maxval(data)
   min_ = minval(data)
   bound(0) = min_
   do i = 1,n_boxs
    bound(i) = bound(i-1) - (max_ - min_)/real(n_boxs)
   end do
   do i = 1,n_boxs
    histo(i) = count( data <= bound(i) .and. data >= bound(i-1))
   enddo
   suma = SUM(histo)
   histo(0:n_boxs) = histo(0:n_boxs) / suma
   do i=1,n_boxs
    histo(i)=histo(i)/ave
   end do
   call moment(data,n_datos,ave,adev,sdev,var,skew,curt)
   return
 END SUBROUTINE histogram
! 
 SUBROUTINE moment(data,n,ave,adev,sdev,var,skew,curt)
  IMPLICIT NONE
! Numerical Recipes (Fortran 90), pp 607-608.
! Given an array of data(1:n), its returns its mean ave, average deviation adev,
! standar deviation sdev, variance var, skewness skew, and kurtosis curt.
  INTEGER :: n,j
  REAL    :: adev,ave,curt,sdev,skew,var,data(n)
  REAL    :: p,ep
  REAL    :: s = 0.0
  IF (n<=1) PRINT*,'n must be at least 2 in moment'
  !DO j=1,n
  ! s=s+data(j)
  !END DO
  ave  = sum(data(1:n))/real(n)
  adev = sum(abs(data(1:n)-ave))/real(n)
  var  = 0.0
  skew = 0.0
  curt = 0.0
  ep   = 0.0
  storage: DO j=1,n
   s    = data(j) - ave
   ep   = ep + s
   p    = s*s
   var  = var + p
   p    = p*s
   skew = skew + p
   p    = p*s
   curt = curt + p
  END DO storage
  var  = (var-ep*ep/real(n))/(n-1)
  sdev = sqrt(var)
  IF(var/=0.0)THEN
   skew = skew/(n*adev**3)
   curt = curt/(n*var*var)-3.0
  ELSE
   skew=0.0
   curt=0.0
  end if
  RETURN
 END SUBROUTINE
end module histograms 
!
module GeometricProperties
 use GeneralVariables
 implicit none
 public
 contains
 !
 pure function Crystal2BoxCoordinates(rv, r_c) result (r_b)
  implicit none
  real,intent(in)  :: r_c(1:3), rv(1:3,1:3)
  real             :: r_b(1:3)
  integer          :: j
  forall ( j=1:3 )
   r_b(j) = rv(j,1)*r_c(1)  + rv(j,2)*r_c(2)  + rv(j,3)*r_c(3)
  end forall
 end function Crystal2BoxCoordinates
!
 pure function Box2CrystalCoordinates(vr,r_b) result (r_c)
  implicit none
  real,intent(in)  :: r_b(1:3), vr(1:3,1:3)
  real             :: r_c(1:3)
  integer          :: j
  forall ( j=1:3 )
   r_c(j) = mod(vr(j,1)*r_b(1)  + vr(j,2)*r_b(2)  + vr(j,3)*r_b(3) + 100.0,1.0)
  end forall
 end function Box2CrystalCoordinates
!
 subroutine ThreeCoordinatesInSameSpace(rv,r1,r2,r3,v1,v2,v3)
  implicit none
  real,intent(in)  :: r1(1:3),r2(1:3),r3(1:3),rv(1:3,1:3)!
  real,intent(out) :: v1(1:3),v2(1:3),v3(1:3)
  real             :: s
  v2 = r2
  call make_distances(.true.,r1,v2,rv,v1,s)
  call make_distances(.true.,r3,v2,rv,v3,s)
  return
 end subroutine ThreeCoordinatesInSameSpace
!
 subroutine FourCoordinatesInSameSpace(rv,r1,r2,r3,r4,v1,v2,v3,v4)
  implicit none
  real,intent(in)  :: r1(1:3),r2(1:3),r3(1:3),r4(1:3),rv(1:3,1:3)!
  real,intent(out) :: v1(1:3),v2(1:3),v3(1:3),v4(1:3)
  real             :: s
  v1 = r1 ! First "pivote" V1 -> V2 -> V3 -> V4
  call make_distances(.true.,r2,v1,rv,v2,s)
  call make_distances(.true.,r3,v2,rv,v3,s)
  call make_distances(.true.,r4,v3,rv,v4,s)
  return
 end subroutine FourCoordinatesInSameSpace
!
 subroutine Cell(rv,vr,cell_0)
 implicit none
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real, parameter   :: pi = acos(-1.0)
 real :: alp, bet
 real :: cosa, cosb, cosg
 real :: gam, sing
 real :: degtorad
 degtorad=pi/180.0
 if(cell_0(4) == 90.0) then
  cosa = 0.0
 else
  alp=cell_0(4)*degtorad
  cosa=cos(alp)
 endif
 if(cell_0(5) == 90.0) then
  cosb = 0.0
 else
  bet = cell_0(5)*degtorad
  cosb = cos(bet)
 endif
 if(cell_0(6) == 90.0) then
  sing = 1.0
  cosg = 0.0
 else
  gam = cell_0(6)*degtorad
  sing = sin(gam)
  cosg = cos(gam)
 endif
 rv(1,1) = cell_0(1)
 rv(1,2) = cell_0(2)*cosg
 rv(1,3) = cell_0(3)*cosb
 rv(2,1) = 0.0
 rv(2,2) = cell_0(2)*sing
 rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
 rv(3,1) = 0.0
 rv(3,2) = 0.0
 rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3))
 call inverse(rv,vr,3)
 return
 end subroutine cell
!
 subroutine uncell(rv,cell_0)
  implicit none
  real,intent(out)   :: cell_0(6)
  real,intent(in)    :: rv(3,3)
  integer            :: i,j
  real               :: temp(6)
  real, parameter    :: pi = acos(-1.0)
  real, parameter    :: radtodeg = 180.0/pi
  do i = 1,3
   temp(i) = 0.0
   do j = 1,3
    temp(i) = temp(i) + rv(j,i)*rv(j,i)
   end do
   temp(i) = sqrt(temp(i))
  enddo
  cell_0(1) = abs(temp(1))
  cell_0(2) = abs(temp(2))
  cell_0(3) = abs(temp(3))
  do i = 1,3
   temp(3+i) = 0.0
  enddo
  do j = 1,3
   temp(4) = temp(4) + rv(j,2)*rv(j,3)
   temp(5) = temp(5) + rv(j,1)*rv(j,3)
   temp(6) = temp(6) + rv(j,1)*rv(j,2)
  enddo
  temp(4) = temp(4)/(temp(2)*temp(3))
  temp(5) = temp(5)/(temp(1)*temp(3))
  temp(6) = temp(6)/(temp(1)*temp(2))
  cell_0(4) = radtodeg*acos(temp(4))
  cell_0(5) = radtodeg*acos(temp(5))
  cell_0(6) = radtodeg*acos(temp(6))
  do i=4,6
   if (abs(cell_0(i) - 90.0 ) < 0.00001) cell_0(i) = 90.0
   if (abs(cell_0(i) - 120.0) < 0.00001) cell_0(i) = 120.0
  end do
  return
 end subroutine uncell
!
 subroutine Inverse(a,c,n)
  implicit none
  integer    :: n
  real       :: a(n,n), c(n,n)
  real       :: l(n,n), u(n,n), b(n), d(n), x(n)
  real       :: coeff
  integer    :: i,j,k
  l=0.0 ; u=0.0 ; b=0.0
  do k=1,n-1
   do i=k+1,n
    coeff=a(i,k)/a(k,k)
    l(i,k) = coeff
    do j=k+1,n
     a(i,j) = a(i,j) - coeff*a(k,j)
    end do
   end do
  end do
  do i=1,n
   l(i,i) = 1.0
  end do
  do j=1,n
   do i=1,j
     u(i,j) = a(i,j)
   end do
  end do
  do k=1,n
   b(k)=1.0
   d(1) = b(1)
   do i=2,n
     d(i)=b(i)
     do j=1,i-1
       d(i) = d(i) - l(i,j)*d(j)
     end do
   end do
   x(n)=d(n)/u(n,n)
   do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
     x(i)=x(i)-u(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
   end do
   do i=1,n
    c(i,k) = x(i)
   end do
   b(k)=0.0
  end do
  return
 end subroutine Inverse
!
subroutine make_distances(flag,r2,r1,rv,r3,dist)
 use vector_module, only: dist2vectors, array2vector
 implicit none
 real,    intent(in)  :: r1(3),r2(3),rv(3,3)    ! coordenadas y matriz de cambio
 real,    intent(out) :: dist,r3(1:3)                     ! matriz de distancias n x n
 real                 :: d_image(1:27),image(3,27)        ! array de distancias
 real                 :: distance,rcm(3),phi
 integer              :: k,l,m,n,o,i,j                    ! variables mudas
 real                 :: atom(3),ouratom(3)               ! coordenadas preparadas
 logical              :: flag                             ! out the coordinate of the atom
! {{ calculamos la matriz de distancias }}
  k=0
  do l=-1,1
   do m=-1,1
      do n=-1,1
         k = k + 1
         ouratom = r1
         atom(1) = r2(1) + real(l)
         atom(2) = r2(2) + real(m)
         atom(3) = r2(3) + real(n)
         d_image(k) = dist2vectors( array2vector(Crystal2BoxCoordinates(rv,ouratom)),&
                                    array2vector(Crystal2BoxCoordinates(rv,atom   )) )
         !d_image(k) = distance(atom,ouratom,rv)
         forall ( i=1:3)       ! r3 sera la imagen de atom a la menor distancia
          image(i,k) = atom(i) ! entre todas las imagenes.
         end forall
     enddo
   enddo
  enddo
  dist=minval(d_image)
  if(flag)then
   phi=9999999.999999 ! initial infinite distance
   k=1                ! selected image
   do l=1,27
    if(d_image(l)<=phi)then
      phi=d_image(l) ! seleccionamos la imagen con la menor distancia
      k=l            !
    endif
   enddo
   forall ( l=1:3)
     r3(l)=image(l,k)
   end forall
  else
   r3(1:3)=0.0
  end if
 return
end subroutine make_distances
!
pure real function volume(rv)
  implicit none
  real, intent(in)  :: rv(3,3)
  real       :: r1x
  real       :: r1y
  real       :: r1z
  real       :: r2x
  real       :: r2y
  real       :: r2z
  real       :: r3x
  real       :: r3y
  real       :: r3z
  real       :: vol
!
  r1x = rv(1,1)
  r1y = rv(2,1)
  r1z = rv(3,1)
  r2x = rv(1,2)
  r2y = rv(2,2)
  r2z = rv(3,2)
  r3x = rv(1,3)
  r3y = rv(2,3)
  r3z = rv(3,3)
  vol = r1x*(r2y*r3z - r2z*r3y) + r1y*(r3x*r2z - r3z*r2x) + r1z*(r2x*r3y - r2y*r3x)
  volume = abs(vol)
  return
end function
!
 pure real function distance(atom,ouratom,rv)
  implicit none
  integer :: j
  real,intent(in) :: atom(3), ouratom(3), rv(3,3)
  real            :: dist(3),o_atom(3),o_ouratom(3)
  forall ( j=1:3 )
   o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
   o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
   dist(j) = o_ouratom(j) - o_atom(j)
  end forall
  distance = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
 end function
!
 subroutine output_gulp(ciffiles)
  use types
  use functions
  implicit none
  type(ciffile),intent(inout)     :: ciffiles
  character(len=100)              :: gulpfilename, gulpfilenameout
  integer                         :: u
  integer                         :: i,k,ierr
  real                            :: mmm,rrr
  integer                         :: zzz
  character(len=2)                :: zlz
  character(len=4)                :: extension=".gin"
  character(len=500)              :: command
!
  gulpfilename=ciffiles%filename(1:clen_trim(ciffiles%filename)-4)//extension
  gulpfilenameout=adjustl(ciffiles%filename(1:clen_trim(ciffiles%filename)-4)//".gout")
  gulpfilename=adjustl(gulpfilename)
  open(newunit=u,file=gulpfilename,iostat=ierr)
  if (ierr/=0) stop 'Error'
  if (shellonly_flag) then
   write(6,'(a)') 'Single-point calculation (optimising shells):' 
   write(u,'(a)') 'opti conv nosymm shellonly'
  else
   if(opti_flag) then
    write(6,'(a)') 'Optimising structure:'
    write(u,'(a)') 'opti conp nosymm'
   end if
  end if
  write(u,'(a,5x,a)')'name',ciffiles%filename(1:clen_trim(ciffiles%filename)-4)
  write(u,'(a)') 'switch_min rfo gnorm 0.08'//new_line('a')//&
   'cuts 1.1'//new_line('a')//&
   'stepmx opt 0.1'//new_line('a')//&
   'switch_stepmx 0.3 gnorm 0.005'//new_line('a')//&
   !'ftol 1e-005'//new_line('a')//&
   !'gtol 0.0001'//new_line('a')//&
   !'xtol 1e-005'//new_line('a')//&
   'maxcyc 500'//new_line('a')//&
   'cell'
  write(u,'(6(f9.5,1x))') (ciffiles%cell_0(i) , i=1,6)
  write(u,'(a)')'fractional'
  do i=1,ciffiles%n_atoms
   write(u,'(a4,1x,3(f14.7,1x),1x,f14.7)')ciffiles%atom(i)%label,&
    ciffiles%atom(i)%uCoordinate%x,&
    ciffiles%atom(i)%uCoordinate%y,&
    ciffiles%atom(i)%uCoordinate%z,&
    ciffiles%atom(i)%charge
  end do
  write(u,'(a)')'space'
  write(u,'(a)')' 1'
  write(u,'(a)')'supercell 1 1 1'
  write(u,'(a)')'species'
  write(u,'(a)')'O  core  O  core'
  write(u,'(a)')'O  shel  O  shel'
  write(u,'(a)')'Si core  Si core'
  write(u,'(a)')'library catlow nodump'
  if(opti_flag)then
   write(u,'(a)')'output cif '//ciffiles%filename(1:clen_trim(ciffiles%filename)-4)//'_opt'
  end if
  write(u,'(a)')'dump every 1 '//ciffiles%filename(1:clen_trim(ciffiles%filename)-4)//'.grs'
  close(u)
  if (opti_flag) then
   if(restart_flag)then
    command = "export OMP_NUM_THREADS=1; gulp < "//&
              ciffiles%filename(1:clen_trim(ciffiles%filename)-4)//'.grs > '//gulpfilenameout(1:clen(gulpfilenameout))
    call system(command)
   else
    command = "export OMP_NUM_THREADS=1; gulp < "//&
              gulpfilename(1:clen(gulpfilename))//" > "//gulpfilenameout(1:clen(gulpfilenameout))
    call system(command)
   end if
   ! **** Optimisation achieved ****
   if(optimisation_achieved(gulpfilenameout(1:clen(gulpfilenameout))))then
    write(6,'(a)') 'Optimisation achieved'
   else
    if(shellonly_flag.eqv..false.)then
     write(6,'(a)') 'Try again with optimisation'
     ! Remove shell particles
     command = "sed -i '/O     shel/d' "//ciffiles%filename(1:clen_trim(ciffiles%filename)-4)//".grs"
     call system(command)
     command = "export OMP_NUM_THREADS=1; gulp < "//&
              ciffiles%filename(1:clen_trim(ciffiles%filename)-4)//'.grs > '//gulpfilenameout(1:clen(gulpfilenameout))
     call system(command)
     if(optimisation_achieved(gulpfilenameout(1:clen(gulpfilenameout)))) then
      write(6,'(a)')'Optimisation achieved'
     else
      write(6,'(a)') 'Error in optimisation: '//ciffiles%filename(1:clen_trim(ciffiles%filename)-4)
     end if
    end if
   end if
   command = "grep 'l e' "//gulpfilenameout(1:clen(gulpfilenameout))//"  | awk '{print $4}' > tmp"
   call system(command)
   open(555, file="tmp")
   read(555, *) CIFFiles%energy
   close(555)
   call system("rm -rf tmp")
  else 
   CIFFiles%energy = 0.0
  end if
  deallocate(CIFFiles%atom)
  return
 end subroutine output_gulp
!
 logical function optimisation_achieved(abc)
  use functions
  implicit none
  character(len=100),intent(in)   :: abc
  character(len=100)              :: command
  integer                         :: i
  optimisation_achieved = .false.
  call system ('x=$(grep "**** Optimisation achieved ****" '//abc(1:clen(abc))//&
          ' | wc -l); if [ "$x" -eq "0" ]; then echo ".false" > tmp ; else echo ".true." > tmp ; fi')
  open(555, file="tmp"); read(555, *) optimisation_achieved; close(555) ; call system("rm -rf tmp")
  return
 end function optimisation_achieved
!
 subroutine output_extended_xyz(CIF)
  use types
  use functions
  implicit NONE
  type(CIFFile),intent(in) :: CIF
  character(len=100)       :: extended_xyz_filename
  integer                  :: u, i
  extended_xyz_filename=CIF%filename(1:clen_trim(CIF%filename)-4)//".exyz"
  extended_xyz_filename=adjustl(extended_xyz_filename)
  open(newunit=u, file=extended_xyz_filename)
  write(u,'(i0)') CIF%n_atoms
  write(u,'(a4)') "%PBC"
  do i = 1, CIF%n_atoms
   write(u,'(a4,1x,4(f14.7,1x))')CIF%atom(i)%label,&
    CIF%atom(i)%xCoordinate%x,CIF%atom(i)%xCoordinate%y,CIF%atom(i)%xCoordinate%z,CIF%atom(i)%Q
  end do
  write(u,'(a)')" "
  write(u,'(a7,1x,3(f14.7))')'Vector1',( CIF%rv(1,i), i=1,3)
  write(u,'(a7,1x,3(f14.7))')'Vector2',( CIF%rv(2,i), i=1,3)
  write(u,'(a7,1x,3(f14.7))')'Vector3',( CIF%rv(3,i), i=1,3)
  write(u,'(a7,1x,3(f14.7))')'Offset ',( 0.0, i=1,3)
  close(u)
 end subroutine output_extended_xyz
!
end module geometricproperties
!
module GetStructures
 use types
 use functions
 use GeometricProperties
 use vector_module
 implicit none
 private
 public GenerateCIFFileList, ReadListOfCIFFiles, AnalysisCIFFileList
 contains
!
 subroutine GenerateCIFFileList()
   implicit none
   character(len=1000)  :: string = " "
   write(string,'(a,1x,a)')"if [ -f list ] ; then rm list ; touch list ; fi",&
    "; for file in ReferenceDataBase/*.cif ; do echo $file ; done > list"
   call system(string)
 end subroutine GenerateCIFFileList
!
 subroutine ReadListOfCIFFiles(n_files)
   implicit none
   character(len=100)  :: line = " "
   integer             :: ierr = 0
   integer,intent(out) :: n_files
   n_files = 0
   open(111,file="list",iostat=ierr)
   read_list: do
    read(111,'(a)',iostat=ierr) line
    if(ierr/=0) exit read_list
    n_files=n_files+1
   end do read_list
   close(111)
   return
 end subroutine
!
 subroutine timestamp ( )
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!    31 May 2001   9:45:54.872 AM
  implicit none
  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  call date_and_time ( values = values )
  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)
  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if
  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
  return
 end subroutine
! ..............................................................................
 subroutine AnalysisCIFFileList(n_files,CIFFiles)
   implicit none
   real                        :: energy = 0.0
   real                        :: lene = 0.0
   real                        :: FD = 0.0
   real,parameter              :: Quartz = -128.703350846667
   real,parameter              :: eV2kJmol = 96.48530749925793
   integer                     :: i,k,u
   integer                     :: ierr = 0
   integer,intent(in)          :: n_files
   type(CIFfile),intent(inout) :: CIFFiles(n_files)
   type(CollectiveVariable)    :: Descriptor(1:8)
   character(len=100)          :: line = " "
   open(newunit=u,file="list",iostat=ierr)
   if (ierr/=0) stop "[ERROR] Problem with list file"
   read_CIFFiles: do i=1,n_files
    read(u,'(a)')line
    read(line(1:49),'(a)') CIFFiles(i)%filename
    write(6,'(a,1x,a)')"File name:",trim(CIFFiles(i)%filename)
    call ReadCIFFile(CIFFiles(i),"rasp",.false.)
    call output_gulp(CIFFiles(i))
    if(opti_flag)then
     call system("cp "//trim(CIFFiles(i)%filename)//" "//trim(CIFFiles(i)%filename)//".back")
     call system("mv "//adjustl(CIFFiles(i)%filename(1:clen_trim(CIFFiles(i)%filename)-4))//"_opt.cif "//&
          trim(CIFFiles(i)%filename))
    end if
    call ReadCIFFile(CIFFiles(i),"rasp",.true.)
    !call output_extended_xyz(CIFFiles(i))
    call TopologicalAnalysis(CIFFiles(i), Descriptor)
    energy = CIFFiles(i)%energy/(CIFFiles(i)%n_atoms/3.0)
    lene = LEnergy((1000.0*CIFFiles(i)%n_atoms/3.0)/volume(CIFFiles(i)%rv),Descriptor)
    FD = (1000.0*CIFFiles(i)%n_atoms/3.0)/volume(CIFFiles(i)%rv)
    write(6,'(a,3(f20.10,1x,a))')'Energy: ', energy-Quartz,'[eV/T], Estimated:', lene-Quartz, &
      '[eV/T], Absolute Error:', abs(energy-lene),'[eV/T]'
    write(6,'(a,3(f20.10,1x,a))')'Energy: ',eV2kJmol*(energy-Quartz),'[kJ/mol], Estimated:',&
      eV2kJmol*(lene-Quartz), '[kJ/mol], Absolute Error:', eV2kJmol*abs(energy-lene),'[kJ/mol]'
    write(6,'(a,f20.10,1x,a)')'Density:', FD,'[A^3/1000T]'
    write(6,'(a,f20.10,1x,a)')'Energy density: ',1000.0*(energy-Quartz)/FD, '[eV/A^3]'
    write(6,'(30a,5a,30a)')("=",k=1,30),' End ',("=",k=1,30)
   end do read_CIFFiles
   close(u)
   return
 end subroutine AnalysisCIFFileList
! 
 subroutine ReadCIFFile(CIF,mode,visual)
   implicit none
   character(len=4),intent(in) :: mode
   logical,intent(in)          :: visual
   integer                     :: i,j,k,u,uu
   integer                     :: ierr = 0
   type(CIFfile),intent(inout) :: CIF
   character(len=100)          :: line = " "
   character(len=100)          :: filename = " "
   character(len=20)           :: spam
   character(len=80)           :: string_stop_head
   select case(mode)
    case("gulp")
     string_stop_head = "_atom_site_occupancy"
    case default
     string_stop_head = "_atom_site_charge"
   end select
   open(newunit=uu,file=trim(CIF%filename),status='old',iostat=ierr)
   if(ierr/=0) stop 'CIFFile does not found'   
   read_cif: do
    read(uu,'(a)',iostat=ierr) line
    if (ierr/=0) exit read_cif
    if (line(1:14)=="_cell_length_a") then
     read(line,*)spam,CIF%cell_0(1)
     cycle read_cif
    end if
    if (line(1:14)=="_cell_length_b") then
     read(line,*)spam,CIF%cell_0(2)
     cycle read_cif
    end if
    if (line(1:14)=="_cell_length_c") then
     read(line,*)spam,CIF%cell_0(3)
     cycle read_cif
    end if
    if (line(1:17)=="_cell_angle_alpha") then
     read(line,*)spam,CIF%cell_0(4)
     cycle read_cif
    end if
    if (line(1:16)=="_cell_angle_beta") then
     read(line,*)spam,CIF%cell_0(5)
     cycle read_cif
    end if
    if (line(1:17)=="_cell_angle_gamma") then
     read(line,*)spam,CIF%cell_0(6)
     cycle read_cif
    end if
    if (line(1:)==string_stop_head) exit read_cif
   end do read_cif
   call cell(CIF%rv,CIF%vr,CIF%cell_0)
   CIF%n_atoms=0
   read_natoms: do
     read(uu,'(a)',iostat=ierr) line
     if ( ierr /= 0 ) exit read_natoms
     if ( trim(line) .eq. '' ) cycle read_natoms ! avoid blank lines
     CIF%n_atoms=CIF%n_atoms+1
   end do read_natoms
   rewind(uu)
   allocate( CIF%atom(1:CIF%n_atoms) )
   if(visual)write(6,'(a,1x,i4,1x,a)')'Atoms:', CIF%n_atoms, CIF%filename
   if(visual)write(6,'(a)')'transformation matrix:'
   if(visual)write(6,'(3(f14.7,1x),a,3(f14.7,1x))')(CIF%rv(1,j),j=1,3),' : ',(CIF%vr(1,j),j=1,3)
   if(visual)write(6,'(3(f14.7,1x),a,3(f14.7,1x))')(CIF%rv(2,j),j=1,3),' : ',(CIF%vr(2,j),j=1,3)
   if(visual)write(6,'(3(f14.7,1x),a,3(f14.7,1x))')(CIF%rv(3,j),j=1,3),' : ',(CIF%vr(3,j),j=1,3)
   do
    read(uu,'(a)') line
    if (line(1:)==string_stop_head) exit
   end do
   if(visual)write(6,'(a)')'Atom coordinates:'
   read_xyz: do j=1,CIF%n_atoms
    read(uu,'(a)', iostat=ierr) line
    if ( trim(line) .eq. '' ) cycle read_xyz
    if ( ierr /= 0 ) exit read_xyz
    select case(mode)
      case("gulp")
      ! mode: "gulp"
      read(line,*) CIF%atom(j)%label,&
       CIF%atom(j)%uCoordinate%x,&
       CIF%atom(j)%uCoordinate%y,&
       CIF%atom(j)%uCoordinate%z,&
       CIF%atom(j)%charge
      CIF%atom(j)%charge = 0.0
      CIF%atom(j)%label_from_CIFFile = CIF%atom(j)%label
      case default
      ! mode: default
      read(line,*) CIF%atom(j)%label_from_CIFFile,CIF%atom(j)%label, &
       CIF%atom(j)%uCoordinate%x,&
       CIF%atom(j)%uCoordinate%y,&
       CIF%atom(j)%uCoordinate%z,&
       CIF%atom(j)%charge
    end select
!    
    CIF%atom(j)%xCoordinate = Array2Vector( Crystal2BoxCoordinates(CIF%rv,&
                                            vector2array(CIF%atom(j)%uCoordinate)))
    call checkatom(CIF%atom(j)%label,  &
                   CIF%atom(j)%mass,   &
                   CIF%atom(j)%radius, &
                   CIF%atom(j)%element,&
                   CIF%atom(j)%label_element)
!
    if(j < 10 .or. j > CIF%n_atoms-10 )then
     if(visual)write(6,'(a4,1x,a2,1x,i2,1x,7(f14.10,1x))') CIF%atom(j)%label,                   &
         CIF%atom(j)%label_element,CIF%atom(j)%element, CIF%atom(j)%mass,             &
         CIF%atom(j)%uCoordinate%x,CIF%atom(j)%uCoordinate%y,CIF%atom(j)%uCoordinate%z,&
         CIF%atom(j)%xCoordinate%x,CIF%atom(j)%xCoordinate%y,CIF%atom(j)%xCoordinate%z
    elseif ( j==10 ) then
     if(visual)write(6,'(a)')'[...]'
    end if
   end do read_xyz
   close(uu)
   !deallocate(CIF%atom)
   return
 end subroutine ReadCIFFile
!
 logical function visited_angle(check_visited,i,j,k,n_bends,bend)
  implicit NONE
  logical, intent(in)             :: check_visited
  integer,intent(in)              :: i,j,k
  integer,intent(inout)           :: n_bends
  character(len=21),intent(inout) :: bend(0:n_bends-1)
  integer                         :: l,ii,jj,kk
  visited_angle = .false.
  if (.not.check_visited) return
  if (n_bends /= 0) then
   check_bend: do l = 1, n_bends-1
    read(bend(l)(1:21),'(3i7)') ii,jj,kk
    if((ii==i.and.jj==j.and.kk==k).or.&
       (ii==k.and.jj==j.and.kk==i)) then
     visited_angle = .true.
     exit check_bend
    end if
   end do check_bend
  end if
  return
 end function visited_angle
!
 character(len=40) function add_angle(i,j,k,angle)
  implicit none
  integer,intent(in) :: i,j,k
  real,intent(in)    :: angle
  write(add_angle(1:40),'(3i7,1x,a,1x,f14.7)') i,j,k,'#',angle
  return
 end function add_angle
!
 subroutine TopologicalAnalysis(CIF,Descriptor)
  implicit none
  type(CIFfile),intent(inout) :: CIF
  type(vector)             :: v1, v2, v3, v4
  type(CollectiveVariable),intent(out) :: Descriptor(1:8)
  integer                  :: i, j, k, l, m, n, o, N_Q
  integer                  :: sio_unit, osio_unit, siosi_unit, ooo_unit, sisisi_unit, sisi_unit, q_unit,NMR_SiSi_unit
  integer                  :: n_bends, ios, bends_unit, staggering_unit
  real                     :: DistanceMatrix(CIF%n_atoms,CIF%n_atoms)
  logical                  :: ConnectedAtoms(CIF%n_atoms,CIF%n_atoms)
  logical                  :: flag =.false.
  real                     :: r, r1(1:3), r2(1:3), r3(1:3), r4(1:3), angle, angle_tmp, theta, rho, dd, r_rho
  character(len=20)        :: abc
  character(len=100)       :: sio_file_name,  sisisi_file_name,NMR_SiSi_filename
  character(len=100)       :: osio_file_name, siosi_file_name, ooo_file_name
  real, parameter          :: r_min_criteria_connectivity = 0.15, degtorad = acos(-1.0)/180.0
  integer, parameter       :: n_bends_max = 10000000
  logical                  :: check_visited = .false.
  character(len=100)       :: bends_filename, staggering_filename, sisi_filename, q_filename
  character(len=50)        :: bends(0:n_bends_max)
  integer                  :: n_staggering_angles, n_staggering_angles_2
  real                     :: staggering_angle(1:9), staggering_angle_2(1:9), Q, q_l,q_m, cos3theta
  logical                  :: visited_staggering(CIF%n_atoms,CIF%n_atoms)
  character(len=40)        :: label
! ...
  n_bends = 0 ; ios = 0
  bends(0:n_bends_max)(1:40) = " "
  staggering_angle(1:9) = 999.99
  staggering_angle_2(1:9) = 999.99
! File names
! Descriptor: density, TO, OTO, TOT, TTT, TT, Q, Energy
! -----------------------------------------------------------
  sio_file_name=CIF%filename(1:clen_trim(CIF%filename)-4)//"_sio.dat"
  osio_file_name=CIF%filename(1:clen_trim(CIF%filename)-4)//"_osio.dat"
  siosi_file_name=CIF%filename(1:clen_trim(CIF%filename)-4)//"_siosi.dat"
  ooo_file_name=CIF%filename(1:clen_trim(CIF%filename)-4)//"_ooo.dat"
  sisisi_file_name=CIF%filename(1:clen_trim(CIF%filename)-4)//"_sisisi.dat"
  bends_filename=CIF%filename(1:clen_trim(CIF%filename)-4)//"_Bonds.dat"
  sisi_filename=CIF%filename(1:clen_trim(CIF%filename)-4)//"_sisi.dat"
  q_filename=CIF%filename(1:clen_trim(CIF%filename)-4)//"_Q.dat"
  staggering_filename=CIF%filename(1:clen_trim(CIF%filename)-4)//"_Staggering.dat"
  NMR_SiSi_filename=CIF%filename(1:clen_trim(CIF%filename)-4)//"_29SiNMR.dat"
! File units:
  osio_unit = 900 ; siosi_unit = 901 ; ooo_unit = 902        ; sisisi_unit = 903 ; NMR_SiSi_unit = 908
  q_unit    = 904 ; bends_unit = 905 ; staggering_unit = 906 ; sisi_unit   = 907 ; sio_unit = 777
!
  open(unit=sio_unit, file=sio_file_name, iostat=ios)
  if ( ios /= 0 ) stop "Error opening sio file "
  open(unit=osio_unit, file=osio_file_name, iostat=ios)
  if ( ios /= 0 ) stop "Error opening osio file "
  open(unit=siosi_unit, file=siosi_file_name, iostat=ios)
  if ( ios /= 0 ) stop "Error opening siosi file "
  open(unit=ooo_unit, file=ooo_file_name, iostat=ios)
  if ( ios /= 0 ) stop "Error opening ooo file "
  open(unit=sisisi_unit, file=sisisi_file_name, iostat=ios)
  if ( ios /= 0 ) stop "Error opening sisisi file "
  open(unit=staggering_unit, file=staggering_filename, iostat=ios)
  if ( ios /= 0 ) stop "Error opening staggering file "
  open(unit=sisi_unit, file=sisi_filename, iostat=ios)
  if ( ios /= 0 ) stop "Error Si-Si file "
  open(unit=q_unit, file=q_filename, iostat=ios)
  if( ios /= 0 ) stop "Error opening file q_filename"
  open(unit=NMR_SiSi_unit, file=NMR_SiSi_filename, iostat=ios)
  if( ios /= 0 ) stop "Error opening file q_filename"
!
  write(6,'(a)')'Detecting natural bonds [...]'
  !
  DistanceMatrix = 0.0
  ConnectedAtoms = .false.
  CIF%atom(:)%degree = 0.0
  do_i_bond: do i = 1, CIF%n_atoms
   do_j_bond: do j = i + 1, CIF%n_atoms
    call make_distances(.false.,vector2array(CIF%atom(j)%uCoordinate), &
                        vector2array(CIF%atom(i)%uCoordinate),CIF%rv,r3,r)
    if(r > 0.1 .and. r <= CIF%atom(i)%radius + CIF%atom(j)%radius + r_min_criteria_connectivity)then
     CIF%atom(i)%degree = CIF%atom(i)%degree + 1
     CIF%atom(j)%degree = CIF%atom(j)%degree + 1
     DistanceMatrix(i,j) = r
     DistanceMatrix(j,i) = DistanceMatrix(i,j)
     ConnectedAtoms(i,j) = .true.
     ConnectedAtoms(j,i) = .true.
    else
     if(r<5.and.CIF%atom(i)%element==8.and.CIF%atom(j)%element==16)then
      write(6,*) r, CIF%atom(i)%radius, CIF%atom(j)%radius, CIF%atom(i)%radius+CIF%atom(j)%radius + r_min_criteria_connectivity
     end if
     DistanceMatrix(i,j) = r
     DistanceMatrix(j,i) = DistanceMatrix(i,j)
    end if
   end do do_j_bond
  end do do_i_bond
!
  write(6,'(a)') 'Connectivity array for each atom:'
  write(abc,'("(",i0,"(i2,1x))")') CIF%n_atoms
  write(6,trim(abc)) (CIF%atom(i)%degree, i=1,CIF%n_atoms)
!
  write(6,'(a)')'Searching for O-(Si)-O-(Si)-O angles:'
  do_i_search_OOO: do i = 1, CIF%n_atoms                                    ! O,  i
   if(CIF%atom(i)%element == 8) then
    do_j_search_OOO: do j=1,CIF%n_atoms                                     ! Si, j
     if(CIF%atom(j)%element == 14.and.ConnectedAtoms(i,j)) then
      !write(6,*) '{DEV}', DistanceMatrix(i,j), i, j
      write(sio_unit,*) DistanceMatrix(i,j)
      do_k_search_OOO: do k=1,CIF%n_atoms                                   ! O,  k
       if(CIF%atom(k)%element == 8.and.ConnectedAtoms(j,k).and.k/=i) then
!       O - Si - O angle: i,j,k
        if (.not.visited_angle(check_visited,i,j,k,n_bends,bends(0:n_bends-1)(1:21))) then
         call ThreeCoordinatesInSameSpace( CIF%rv,&
          Vector2Array(CIF%atom(i)%uCoordinate),&
          Vector2Array(CIF%atom(j)%uCoordinate),&
          Vector2Array(CIF%atom(k)%uCoordinate),&
          r1,r2,r3 )
!
         v1 = Array2Vector(Crystal2BoxCoordinates(CIF%rv,r1))
         v2 = Array2Vector(Crystal2BoxCoordinates(CIF%rv,r2))
         v3 = Array2Vector(Crystal2BoxCoordinates(CIF%rv,r3))
!
         angle = angle3vectors_jik(v2,v1,v3)/degtorad
         n_bends = n_bends + 1
         write(bends(n_bends),'(3i7,1x,f14.7,1x,i4,1x,3(a2,1x))') i,j,k,angle,n_bends,&
          CIF%atom(i)%label_element,CIF%atom(j)%label_element,CIF%atom(k)%label_element
         write(osio_unit,*) angle
        end if
!
        do_l_search_OOO: do l=1,CIF%n_atoms                                     ! Si, l
         if(CIF%atom(l)%element == 14.and.ConnectedAtoms(k,l).and.l/=j) then
!         Si - O - Si angle: j,k,l
          if (.not.visited_angle(check_visited,j,k,l,n_bends,bends(0:n_bends-1)(1:21))) then
           call ThreeCoordinatesInSameSpace( CIF%rv,&
            Vector2Array(CIF%atom(j)%uCoordinate),&
            Vector2Array(CIF%atom(k)%uCoordinate),&
            Vector2Array(CIF%atom(l)%uCoordinate),&
            r1,r2,r3 )
           v1=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r1))
           v2=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r2))
           v3=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r3))
!
           angle = angle3vectors_jik(v2,v1,v3)/degtorad
           n_bends = n_bends + 1
           write(bends(n_bends),'(3i7,1x,f14.7,1x,i4,1x,3(a2,1x))') j,k,l,angle,n_bends,&
            CIF%atom(j)%label_element,CIF%atom(k)%label_element,CIF%atom(l)%label_element
           write(siosi_unit,*) angle
          end if
!
          do_m_search_OOO: do m=1,CIF%n_atoms                                   ! O,  m
           if(CIF%atom(m)%element == 8.and.ConnectedAtoms(m,l).and.m/=k.and.m/=i) then
!           O - O - O angle : i,k,m
            if (.not.visited_angle(check_visited,i,k,m,n_bends,bends(0:n_bends-1)(1:21))) then
             call ThreeCoordinatesInSameSpace( CIF%rv,&
              Vector2Array(CIF%atom(i)%uCoordinate),&
              Vector2Array(CIF%atom(k)%uCoordinate),&
              Vector2Array(CIF%atom(m)%uCoordinate),&
              r1,r2,r3 )
             v1=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r1))
             v2=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r2))
             v3=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r3))
             !
             angle = angle3vectors_jik(v2,v1,v3)/degtorad
             n_bends = n_bends + 1
             write(bends(n_bends),'(3i7,1x,f14.7,1x,i4,1x,3(a2,1x))') i,k,m,angle,n_bends,&
              CIF%atom(i)%label_element,CIF%atom(k)%label_element,CIF%atom(m)%label_element
             write(ooo_unit,*) angle
            end if
!
            do_n_search_OOO: do n=1,CIF%n_atoms                                 ! Si, n
             if(CIF%atom(n)%element == 14.and.ConnectedAtoms(m,n).and.n/=l.and.n/=j) then
!            Si - Si - Si angle: j,l,n
              if (.not.visited_angle(check_visited,j,l,n,n_bends,bends(0:n_bends-1)(1:21))) then
               call ThreeCoordinatesInSameSpace( CIF%rv,&
                Vector2Array(CIF%atom(j)%uCoordinate),&
                Vector2Array(CIF%atom(l)%uCoordinate),&
                Vector2Array(CIF%atom(n)%uCoordinate),&
                r1,r2,r3 )
               v1=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r1))
               v2=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r2))
               v3=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r3))
               !
               angle = angle3vectors_jik(v2,v1,v3)/degtorad
               n_bends = n_bends + 1
               write(bends(n_bends),'(3i7,1x,f14.7,1x,i4,1x,3(a2,1x))') j,l,n,angle,n_bends,&
                CIF%atom(j)%label_element,CIF%atom(l)%label_element,CIF%atom(n)%label_element
               write(sisisi_unit,*) angle
              end if
             end if
            end do do_n_search_OOO
           end if
          end do do_m_search_OOO
         end if
        end do do_l_search_OOO
       end if
      end do do_k_search_OOO
     end if
    end do do_j_search_OOO
   end if
  end do do_i_search_OOO
  open(unit=bends_unit, file=bends_filename, iostat=ios)
  if ( ios /= 0 ) stop "Error opening file bends_filename"
  do i=1,n_bends
   write(bends_unit,'(a)') bends(i)
  end do
!
  Q=0.0
  N_Q=0
  CIF%atom(:)%Q = 0.0
! Tetrahedral order parameter
! and 29Si NMR chemical Shift
! https://pubs.acs.org/doi/suppl/10.1021/jacs.5b08098/suppl_file/ja5b08098_si_001.pdf
  do_i_search_q: do i=1,CIF%n_atoms
   if(CIF%atom(i)%element == 14) then
    N_Q=N_Q+1
    q_l = 0.0
    l = 0
    do_j_search_q: do j=1,CIF%n_atoms
     if( CIF%atom(j)%element == 8.and.ConnectedAtoms(j,i)) then
     do_k_search_q: do k=1,CIF%n_atoms
      if (CIF%atom(k)%element == 8.and.ConnectedAtoms(i,k).and.j/=k) then
       q_m = 0.0 
       do_m_search_q: do m=1,CIF%n_atoms
        if( CIF%atom(m)%element == 8.and.ConnectedAtoms(i,m).and.m/=j.and.m/=k) then
         l = l + 1
         call FourCoordinatesInSameSpace(CIF%rv,&
          Vector2Array(CIF%atom(i)%uCoordinate),&
          Vector2Array(CIF%atom(j)%uCoordinate),&
          Vector2Array(CIF%atom(k)%uCoordinate),&
          Vector2Array(CIF%atom(m)%uCoordinate),&
          r1,r2,r3,r4 )
!
          v1=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r1))   ! i 
          v2=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r2))   ! j
          v3=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r3))   ! k
          v4=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r4))   ! m   
!
          cos3theta = abs(cos(3.0*dihedral_angle4vectors_ijkm(v1,v2,v3,v4)))
          angle = angle3vectors_jik(v1,v2,v4)               ! angle i-j-m
          q_m = q_m + cos3theta*exp(-0.5*(angle-acos(-1.0/3.0))**2/(12.0*degtorad)**2)
        end if
       end do do_m_search_q
!      i,j,k
       angle = angle3vectors_jik(v1,v2,v3)
       q_l = q_l + exp(-0.5*(angle-acos(-1.0/3.0))**2/(12.0*degtorad)**2)*q_m/24.0
      end if
     end do do_k_search_q
     end if
    end do do_j_search_q
    write(q_unit,*) q_l
    CIF%atom(i)%Q = q_l
    !Q= Q + q_l
   end if
  end do do_i_search_q
! 
  visited_staggering(1:CIF%n_atoms,1:CIF%n_atoms) = .false.
!
  do_j_search_staggering: do j = 1, CIF%n_atoms                            ! Si, j
   if(CIF%atom(j)%element == 14) then
    do_k_search_staggering: do k=1,CIF%n_atoms                             ! O,  k
    if(CIF%atom(k)%element == 8.and.ConnectedAtoms(j,k)) then
     do_l_search_staggering: do l=1,CIF%n_atoms                            ! Si, l
      if(CIF%atom(l)%element == 14.and.ConnectedAtoms(k,l).and.j/=l) then
      if(visited_staggering(j,l).eqv..false.)then
! ---- Si - Si distance:
       write(sisi_unit,*)DistanceMatrix(j,l)
       CIF%atom(j)%N_SiSi=CIF%atom(j)%N_SiSi+1
       CIF%atom(l)%N_SiSi=CIF%atom(l)%N_SiSi+1
       ! NMR
       !----
       CIF%atom(j)%Sum_SiSi_d=CIF%atom(j)%Sum_SiSi_d+DistanceMatrix(j,l)
       CIF%atom(l)%Sum_SiSi_d=CIF%atom(l)%Sum_SiSi_d+DistanceMatrix(l,j)
       !----
       call ThreeCoordinatesInSameSpace( CIF%rv,&
        Vector2Array(CIF%atom(j)%uCoordinate),&
        Vector2Array(CIF%atom(k)%uCoordinate),&
        Vector2Array(CIF%atom(l)%uCoordinate),&
        r1,r2,r3 )
       v1=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r1))
       v2=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r2))
       v3=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r3))
!
       angle = angle3vectors_jik(v2,v1,v3)/degtorad
       CIF%atom(j)%Sum_SiOSi_angle=CIF%atom(j)%Sum_SiOSi_angle+angle
       CIF%atom(l)%Sum_SiOSi_angle=CIF%atom(l)%Sum_SiOSi_angle+angle
! angle SiOSi:
       CIF%atom(j)%Sum_rho = CIF%atom(j)%Sum_rho + cos(angle/degtorad)/(cos(angle/degtorad)+1)
       CIF%atom(l)%Sum_rho = CIF%atom(l)%Sum_rho + cos(angle/degtorad)/(cos(angle/degtorad)+1)
! Sum_SiOSi_dot_SiO_product
       rho = cos(degtorad*angle)/(cos(degtorad*angle)+1)
       CIF%atom(j)%Sum_SiOSi_dot_SiO_product=CIF%atom(j)%Sum_SiOSi_dot_SiO_product+&
              rho * DistanceMatrix(j,k)
       CIF%atom(l)%Sum_SiOSi_dot_SiO_product=CIF%atom(l)%Sum_SiOSi_dot_SiO_product+&
              rho * DistanceMatrix(l,k) 
! ----  Staggering Angle:
! {{ Dihedral angle O-(si)-(Si)-O: i, j, l, m atoms:
!    First order:
!    ------------
       visited_staggering(j,l) = .true. ; visited_staggering(l,j) = .true.
       staggering_angle(1:9)   = 999.99 ; n_staggering_angles = 0
       staggering_angle_2(1:9) = 999.99 ; n_staggering_angles_2 = 0
! <
       do_i_search_staggering: do i = 1, CIF%n_atoms                       ! O,  i
        if(CIF%atom(i)%element == 8.and.ConnectedAtoms(i,j).and.i/=k) then
         do_m_search_staggering: do m = 1, CIF%n_atoms                     ! O,  m
          if(CIF%atom(m)%element == 8.and.ConnectedAtoms(m,l).and.m/=i.and.m/=k) then
! {{ First order Stagering angle
!    ---------------------------
           n_staggering_angles = n_staggering_angles + 1
!
           call FourCoordinatesInSameSpace(CIF%rv,&
              Vector2Array(CIF%atom(i)%uCoordinate),Vector2Array(CIF%atom(j)%uCoordinate),&
              Vector2Array(CIF%atom(l)%uCoordinate),Vector2Array(CIF%atom(m)%uCoordinate),&
              r1,r2,r3,r4 )
!
           v1=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r1))
           v2=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r2))
           v3=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r3))
           v4=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r4))
!
           angle = dihedral_angle4vectors_ijkl(v1, v2, v3, v4)/degtorad
           staggering_angle(n_staggering_angles) = angle
! }}
! {{ Second order:
!     ---------------
          do_o_search_staggering: do o = 1, CIF%n_atoms                       ! 
           if(CIF%atom(o)%element == 14 .and.ConnectedAtoms(o,i).and.o/=j.and.o/=l) then
            do_n_search_staggering: do n = 1, CIF%n_atoms                     ! O,  m
             if(CIF%atom(n)%element == 14.and.ConnectedAtoms(m,n).and.n/=l.and.n/=j.and.n/=o) then
             n_staggering_angles_2 = n_staggering_angles_2 + 1
             call FourCoordinatesInSameSpace(CIF%rv,&
              Vector2Array(CIF%atom(o)%uCoordinate),&
              Vector2Array(CIF%atom(j)%uCoordinate),&
              Vector2Array(CIF%atom(l)%uCoordinate),&
              Vector2Array(CIF%atom(n)%uCoordinate),&
              r1,r2,r3,r4 )
             v1=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r1))
             v2=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r2))
             v3=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r3))
             v4=Array2Vector(Crystal2BoxCoordinates(CIF%rv,r4))
             angle = dihedral_angle4vectors_ijkl(v1, v2, v3, v4)/degtorad
             staggering_angle_2(n_staggering_angles_2) = angle
             end if
            end do do_n_search_staggering
           end if
           end do do_o_search_staggering
!   }}
          end if
         end do do_m_search_staggering
        end if
       end do do_i_search_staggering
!
       write(staggering_unit,*)j,l,minval(staggering_angle),minval(staggering_angle_2)
! }}
      end if
      end if
     end do do_l_search_staggering
    end if
    end do do_k_search_staggering
    ! Write NMR estimation from literature:
    theta  = CIF%atom(j)%Sum_SiOSi_angle/4.0           ! [deg.]
    rho    = CIF%atom(j)%Sum_rho/4.0 
    dd     = CIF%atom(j)%Sum_SiSi_d/4.0                ! [A]
    r_rho = CIF%atom(j)%Sum_SiOSi_dot_SiO_product/4.0
    write(NMR_SiSi_unit,*) j, &
     -25.44-0.5793*theta,&  ! Thomas et al., Chem. Phys. Lett., 102, 1983, 158-162, 10.1016/0009-2614(83)87384-9
     247.6-116.7*dd !,     &  ! Fyfe et al., Nature, 326, 281-283, 1987, 10.1038/326281a0
    !   2.19-247.05*rho,  &  ! Engelhardt and Radeglia, Chem. Phys. Lett., 108, 1984, 271-274, 10.1016/0009-2614(84)87063-3
    !  48.54-216.96*r_rho    ! Davis et al., J. Phys. Chem., 100 (1996), 5039-5049, 10.1021/jp9530055
   end if
  end do do_j_search_staggering
  !
  !
! Close units:
  close(unit=sio_unit, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file sio_unit"
  close(unit=NMR_SiSi_unit, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file NMR_unit"
  close(unit=sisi_unit, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file sisi_unit unit "
  close(unit=bends_unit, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file bends_unit unit "
  close(unit=sisisi_unit, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file sisisi unit "
  close(unit=ooo_unit, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file ooo unit "
  close(unit=siosi_unit, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file siosi unit "
  close(unit=osio_unit, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file osio unit "
  close(staggering_unit, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file Staggering unit "
  close(q_unit, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file Q unit "
! Histograms:
  label="Histo: Si-O"
  call HistogramsByFile(sio_unit,sio_file_name,label,Descriptor(1))
  label="Histo: Si-Si"
  call HistogramsByFile(sisi_unit,sisi_filename,label,Descriptor(2))
  label="Histo: O-Si-O   "
  call HistogramsByFile(osio_unit,osio_file_name,label,Descriptor(3))
  label="Histo: Si-O-Si  "
  call HistogramsByFile(siosi_unit,siosi_file_name,label,Descriptor(4))
  label="Histo: Si-Si-Si "
  call HistogramsByFile(sisisi_unit,sisisi_file_name,label,Descriptor(5))
  label="Histo: O-O-O    "
  call HistogramsByFile(ooo_unit,ooo_file_name,label,Descriptor(6))
  label="Histo: Q        "
  call HistogramsByFile(q_unit,q_filename,label,Descriptor(7))
  label="Histo: NMR 29Si "
  call HistogramsByFile(NMR_SiSi_unit,NMR_SiSi_filename,label,Descriptor(8))
!
  return
 end subroutine TopologicalAnalysis
!
 subroutine HistogramsByFile(descriptor_unit,descriptor_file,label,D)
  use histograms
  implicit none
  integer, intent(in)                  :: descriptor_unit
  type(CollectiveVariable),intent(out) :: D
  character(len=100),intent(in)        :: descriptor_file
  integer                              :: ios, i, j, n_datos, n_boxes = 128
  character(len=100)                   :: line
  character(len=20)                    :: label
  real                                 :: ave,adev,sdev,var,skew,curt,suma,max_,min_
  real, allocatable                    :: data_(:)
  open(unit=descriptor_unit, file=descriptor_file, iostat=ios)
  if ( ios /= 0 ) stop "Error opening file"
  n_datos=0
  do
   read(descriptor_unit,'(a)',iostat=ios) line
   if (ios/=0) exit
   n_datos=n_datos+1
  end do
  rewind(descriptor_unit)
  allocate(data_(1:n_datos))
  do i=1,n_datos
   read(descriptor_unit,*) data_(i)
  end do
  call histogram(data_,n_datos,n_boxes,ave,adev,sdev,var,skew,curt,suma,max_,min_)
  max_ = maxval(data_)
  min_ = minval(data_)
  deallocate(data_)
  close(descriptor_unit)
  D%k1   = ave ; D%k2   = sdev ; D%k3   = skew ; D%k4   = curt
  D%min_ = min_ ; D%max_ = max_
  write(6,'(a,1x,6(e20.10,1x))')label,ave,max_,min_,sdev,skew,curt
  label(1:20)=" "
  return
 end subroutine
!
 pure real function norm(x,xmin,xmax)
  implicit none
  real,intent(in) :: x,xmin,xmax
  norm = (x - xmin)/abs(xmax - xmin)
  return
 end function  norm
!
 real function LEnergy(FD,D)
  implicit none
  real,intent(in)                     :: FD
  type(CollectiveVariable),intent(in) :: D(1:7)
  real,parameter                      :: Quartz = -128.703350846667

  LEnergy = 0.0
! Variables:
! label="Histo: Si-O"            1 ! TO
! label="Histo: Si-Si"           2 ! TT
! label="Histo: O-Si-O   "       3 ! OTO
! label="Histo: Si-O-Si  "       4 ! TOT
! label="Histo: Si-Si-Si "       5 ! TTT
! label="Histo: O-O-O    "       6 ! OOO
! label="Histo: Q        "       7 ! Q

  LEnergy = 4.9554 &
  -2.0333 * norm(D(1)%min_,1.3282,1.6548)          & ! TO
  +0.3717 * norm(D(1)%k1,1.5742,1.6660)            &
  +0.6529 * norm(D(1)%max_,1.5744,1.8713)          &
  +0.5799 * (norm(D(1)%min_,1.3282,1.6548))**2     &
  -0.8941 * (norm(D(1)%k1,1.5742,1.6660))**2       &
  -2.0339 * (norm(D(1)%max_,1.5744,1.8713))**2     &
  +0.4708 * (norm(D(1)%min_,1.3282,1.6548))**3     &
  +1.5727 * (norm(D(1)%k1,1.5742,1.6660))**3       &
  +1.8244 * (norm(D(1)%max_,1.5744,1.8713))**3     &
  -3.5413 * norm(D(2)%min_,2.2673,3.2072)          & ! TT
  -3.5034 * norm(D(2)%k1,2.8329,3.2351)            &
  -0.4083 * norm(D(2)%max_,2.9105,3.6102)          &
  +3.5786 * (norm(D(2)%min_,2.2673,3.2072))**2     &
  +2.1687 * (norm(D(2)%k1,2.8329,3.2351))**2       &
  +0.7027 * (norm(D(2)%max_,2.9105,3.6102))**2     &
  -1.1826 * (norm(D(2)%min_,2.2673,3.2072))**3     &
  -0.1351 * (norm(D(2)%k1,2.8329,3.2351))**3       &
  -0.1132 * (norm(D(2)%max_,2.9105,3.6102))**3     &
  -0.1377 * norm(D(3)%min_,63.5395,109.4616)       &  ! OTO
  +0.4955 * norm(D(3)%k1,105.9392,109.9492)        &
  -0.4629 * norm(D(3)%max_,109.4808,163.6060)      &
  +0.5487 * (norm(D(3)%min_,63.5395,109.4616))**2  &
  -5.2814 * (norm(D(3)%k1,105.9392,109.9492))**2   &
  +1.5215 * (norm(D(3)%max_,109.4808,163.6060))**2 &
  -0.2603 * (norm(D(3)%min_,63.5395,109.4616))**3  &
  +3.1999 * (norm(D(3)%k1,105.9392,109.9492))**3   &
  -0.8001 * (norm(D(3)%max_,109.4808,163.6060))**3 &
  +0.3608 * norm(D(4)%min_,90.9591,170.6157)       & ! TOT
  +2.0877 * norm(D(4)%k1,123.9087,175.7053)        &
  +0.4396 * norm(D(4)%max_,124.6023,180.0000)      &
  -0.0377 * (norm(D(4)%min_,90.9591,170.6157))**2  &
  -1.9925 * (norm(D(4)%k1,123.9087,175.7053))**2   &
  -0.6033 * (norm(D(4)%max_,124.6023,180.0000))**2 &
  -0.1869 * (norm(D(4)%min_,90.9591,170.6157))**3  &
  +0.7249 * (norm(D(4)%k1,123.9087,175.7053))**3   &
  +0.2432 * (norm(D(4)%max_,124.6023,180.0000))**3 &
  -0.5293 * norm(D(5)%min_,48.74126,108.6046)      & ! TTT
  +0.3632 * norm(D(5)%k1,101.5302,112.7831)        &
  +0.0459 * norm(D(5)%max_,109.9063,180.0000)      &
  +1.1257 * (norm(D(5)%min_,48.74126,108.6046))**2 &
  -0.7074 * (norm(D(5)%k1,101.5302,112.7831))**2   &
  -0.0536 * (norm(D(5)%max_,109.9063,180.0000))**2 &
  -0.7438 * (norm(D(5)%min_,48.74126,108.6046))**3 &
  +0.4935 * (norm(D(5)%k1,101.5302,112.7831))**3   &
  +0.0475 * (norm(D(5)%max_,109.9063,180.0000))**3 &
  -0.3741 * norm(FD,7.7461,37.9034)                &  ! FD
  -0.7575 * norm(D(7)%k1,0.2100,1.0000)               ! Q
!
  LEnergy = LEnergy + Quartz
! Data for normalization:
! -----------------------
!Descriptor type & AVE min  &  AVE max &  Min min & Min max  & Max min  & Max max
!----------------------------------------------------------------------------------
! FD             & 7.7461   &  37.9034 &  -       & -        & -        & -
! Q              & 0.2100   &  1.0000  &  -       & -        & -        & -
! TO             & 1.5742   &  1.6660  &  1.3282  & 1.6548   & 1.5744   & 1.8713
! OTO            & 105.9392 & 109.9492 &  63.5395 & 109.4616 & 109.4808 & 163.6060
! TOT            & 123.9087 & 175.7053 &  90.9591 & 170.6157 & 124.6023 & 180.0000
! TT             & 2.8329   &   3.2351 &   2.2673 &   3.2072 &   2.9105 &   3.6102
! TTT            & 101.5302 & 112.7831 & 48.74126 & 108.6046 & 109.9063 & 180.0000
  return
 end function LEnergy
!
 subroutine CheckAtom(label,m,s,z,zlabel)
  implicit none
  character(len=4),intent(in)  :: label
  real,intent(out)             :: m, s
  integer,intent(out)          :: z
  character(len=2),intent(out) :: zlabel
  select case(label)
   case('S   ','S0  ':'S999')
    z=8
    m=32.06   ! conventional weight
    s=1.05    ! covalent radious
    zlabel=' O'
   case('O   ','O0  ':'O999')
    z=8
    m=15.999  ! conventional weight
    s=0.66    ! covalent radious
    zlabel=' O'
   case('C   ','C0  ':'C999')
    z=6
    m=12.00  ! conventional weight
    s=0.66    ! covalent radious
    zlabel=' C'
   case('N   ','N0  ':'N999')
    z=7
    m=14.00  ! conventional weight
    s=0.66    ! covalent radious
    zlabel=' N'
   case('H   ','H0  ':'H999')
    z = 1
    m = 1.00784
    s = 0.320
    zlabel=' H'
   case('F   ','F0  ':'F999')
    z = 9
    m = 1.00784
    s = 0.0
    zlabel=' F'
   case('Ge  ','Ge0 ':'Ge99')
    z = 14
    m = 72.630
    s = 1.22  ! covalent radious
    zlabel='Si'
   case('Si  ','Si0 ':'Si99')
    z = 14
    m = 28.0855
    s = 1.11  ! covalent radious
    zlabel='Si'
   case default
    write(6,'(a1,a4,a1)')"'",label,"'"
    stop 'atom unknowed'
  end select
 end subroutine CheckAtom
!
end module GetStructures
!
program ZeoAnalyser
 use, intrinsic                               :: iso_fortran_env
 use GetStructures                            
 use GeneralVariables                         
 use functions
 implicit none                                
 integer                                      :: i
 integer                                      :: num_args
 character(len=100),dimension(:), allocatable :: args
! read arguments in line:
 num_args = command_argument_count()
 allocate(args(num_args))
 do i = 1, num_args
  call get_command_argument(i,args(i))
 end do
 write(6,'(a,1x,i2)')'Arguments:',command_argument_count()
 write(6,'(a)')( args(i),i=1,num_args)
 if(num_args==0) then
  call print_help()
  stop
 end if
 do i=1,num_args
  select case(args(i))
   case ('-nl','--no-create-list')
    generate_list_file = .false.
    opti_flag = .false.
    shellonly_flag = .false.
    restart_flag = .false.
   case ('-r','--restart')
    preoptimization_flag = .true.
   case ('-s','--shellonly')
    shellonly_flag = .true.
   case default
    generate_list_file = .true.
    opti_flag = .false.
    shellonly_flag = .false.
    restart_flag = .false.
  end select
 end do
!
 if (generate_list_file) call GenerateCIFFileList()
 call ReadListOfCIFFiles(n_files)
 allocate( CIFFiles(1:n_files) )
 write(6,'(a,1x,i6)')'Detected CIF Files:', n_files
 write(6,'(a)')'======= Starting Analysis ======='
 call AnalysisCIFFileList(n_files,CIFFiles)
 deallocate(CIFFiles)
 stop "Program finished"
end program ZeoAnalyser
