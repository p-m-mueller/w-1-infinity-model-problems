program main

  use iso_fortran_env
  use quadrature

  implicit none

  double precision, parameter :: pi = 2.d0*asin(1.d0)

  integer, parameter :: nx = 101
  integer, parameter :: ny = 101
  integer, parameter :: nv = nx * ny
  integer, parameter :: ne = (nx - 1) * (ny - 1)

  integer, parameter :: nq_1d = 2
  double precision :: xq_1d(nq_1d), wq_1d(nq_1d)
  integer, parameter :: nq = nq_1d*nq_1d
  double precision :: xq(2,nq), wq(nq)
  
  integer :: element(4,ne)
  double precision :: Jac(2,2,nq,ne), invJac(2,2,nq,ne), detJac(nq,ne)
  double precision :: area(ne)

  integer, parameter :: ne_bnd = 2 * (nx - 1 + ny - 1)

  integer, parameter :: nq_bnd = nq_1d
  integer :: element_bnd(2,ne_bnd)
  double precision :: l_bnd(nq_bnd,ne_bnd)
  double precision :: xq_bnd(nq_bnd), wq_bnd(nq_bnd)

  double precision :: coord(2,nv)
  double precision :: dx, dy

  double precision :: u(nv), old_u(nv) ! solution
  double precision :: grad_u(2,ne), grad_u_(2,ne), norm_grad_u(ne)
  double precision :: norm_q(ne), norm_q1(ne), q(2,ne), lambda(2,ne), lambda_(2,ne)
  double precision :: max_norm_grad_u

  double precision, parameter :: tau = 1.d0

  double precision :: f(nv) ! volume forcing term
  double precision :: g(nv) ! boundary data
  integer :: bctype(nv)
  integer, parameter :: INTERNAL = 0, ESSENTIAL = 1, NATURAL = 2

  double precision :: grad_phi_i(2), grad_phi_j(2)
  double precision :: phi_i, phi_j

  double precision :: Me(4,4,ne) ! element mass matrix
  double precision :: Ke(4,4,ne) ! element stiffnes matrix
  double precision :: fe(4,ne)  ! element right hand side
  double precision :: bc(nv) ! right hand side modification for boundary conditions
  double precision :: pc(nv) ! preconditioner for linear solver
  double precision :: r(nv), resPoisson, resPoisson0, resADMM, resADMM0 ! residual
  
  integer :: i, j, k, iv, jv, ie, iq
  integer :: iterPoisson, iterADMM

  write(*,'(a,i0)') 'Number of elements: ', ne
  write(*,'(a,i0)') 'Number of vertices: ', nv
  write(*,'(a,i0)') 'Number of unknowns: ', nv

  write(*,'(a)') 'Compute quadrature points and weights for boundary'
  call points(nq_bnd,xq_bnd)
  call weights(nq_bnd,xq_bnd,wq_bnd)
  
  write(*,'(a)') 'Compute quadrature points and weights for domain'
  call points(nq_1d,xq_1d)
  call weights(nq_1d,xq_1d,wq_1d)

  do j = 1, nq_1d
    do i = 1, nq_1d
      k = (j-1)*nq_1d+i
      xq(1,k) = xq_1d(j)
      xq(2,k) = xq_1d(i)
      wq(k) = wq_1d(i)*wq_1d(j)
    end do
  end do

  dx = 1.d0 / (nx - 1)
  dy = 1.d0 / (ny - 1)

  ! coordinaes
  do j = 1, ny
    do i = 1, nx
      iv = lex_map(nx, ny, i, j)
      coord(1,iv) = dx * (i-1)
      coord(2,iv) = dy * (j-1)
    end do
  end do

  ! volume forcing term
  do j = 1, ny
    do i = 1, nx
      iv = lex_map(nx, ny, i, j)
      f(iv) = 0.d0
      !f(iv) = 20.d0*sin(3.d0*pi*coord(1,iv))*sin(2.d0*pi*coord(2,iv))
    end do
  end do

  ! boundary type and values
  g = 0.d0
  bctype = INTERNAL
  j = ny
  do i = 1, nx
    iv = lex_map(nx, ny, i, j)
    bctype(iv) = NATURAL 
    g(iv) = 1.d0
  end do
  i = nx
  do j = 1, ny
    iv = lex_map(nx, ny, i, j)
    bctype(iv) = ESSENTIAL 
    g(iv) = 0.d0
  end do
  i = 1
  do j = 1, ny
    iv = lex_map(nx, ny, i, j)
    bctype(iv) = ESSENTIAL 
    g(iv) = 0.d0
  end do
  j = 1
  do i = 1, nx
    iv = lex_map(nx, ny, i, j)
    bctype(iv) = ESSENTIAL
    g(iv) = 0.d0
  end do

  ! element connectivity
  ie = 1
  do j = 1, ny - 1
    do i = 1, nx - 1
      element(1,ie) = lex_map(nx, ny, i, j)
      element(2,ie) = lex_map(nx, ny, i+1, j)
      element(3,ie) = lex_map(nx, ny, i+1, j+1)
      element(4,ie) = lex_map(nx, ny, i, j+1)

      ie = ie + 1
    end do
  end do

  ! boundary element connectivity
  ie = 1
  j = 1
  do i = 1, nx - 1
    element_bnd(1,ie) = lex_map(nx, ny, i, j)
    element_bnd(2,ie) = lex_map(nx, ny, i+1, j)
    
    ie = ie + 1 
  end do
  i = nx
  do j = 1, ny - 1
    element_bnd(1,ie) = lex_map(nx, ny, i, j)
    element_bnd(2,ie) = lex_map(nx, ny, i, j+1)
    
    ie = ie + 1
  end do
  j = ny
  do i = nx - 1, 1, -1
    element_bnd(1,ie) = lex_map(nx, ny, i+1, j)
    element_bnd(2,ie) = lex_map(nx, ny, i, j)

    ie = ie + 1 
  end do
  i = 1
  do j = ny - 1, 1, -1
    element_bnd(1,ie) = lex_map(nx, ny, i, j+1)
    element_bnd(2,ie) = lex_map(nx, ny, i, j)

    ie = ie + 1
  end do

  ! finite element values
  do ie = 1, ne
    area(ie) = 0.d0
    do iq = 1, nq
      Jac(:,:,iq,ie) = 0.d0
      do iv = 1, 4
        j = element(iv,ie)
        grad_phi_i = grad_lin_shape_fun(iv, xq(:,iq))
        Jac(:,1,iq,ie) = Jac(:,1,iq,ie) + grad_phi_i * coord(1,j)
        Jac(:,2,iq,ie) = Jac(:,2,iq,ie) + grad_phi_i * coord(2,j)
      end do

      detJac(iq,ie) = Jac(1,1,iq,ie) * Jac(2,2,iq,ie) - Jac(1,2,iq,ie) * Jac(2,1,iq,ie)

      invJac(1,1,iq,ie) =  Jac(2,2,iq,ie)
      invJac(1,2,iq,ie) = -Jac(1,2,iq,ie)
      invJac(2,1,iq,ie) = -Jac(2,1,iq,ie)
      invJac(2,2,iq,ie) =  Jac(1,1,iq,ie)
      invJac(:,:,iq,ie) = invJac(:,:,iq,ie) / detJac(iq,ie)
      
      area(ie) = area(ie) + detJac(iq,ie) * wq(iq)
    end do
  end do

  ! finite element values for boundary elements
  do ie = 1, ne_bnd
    l_bnd(:,ie) = norm2(coord(:,element_bnd(2,ie)) - coord(:,element_bnd(1,ie)))
  end do

  ! Compute elemet mass and stiffness matrix
  pc = 0.d0
  do ie = 1, ne
  
    Me(:,:,ie) = 0.d0
    Ke(:,:,ie) = 0.d0  
    fe(:,ie) = 0.d0
    do iq = 1, nq
      do jv = 1, 4

        phi_j = lin_shape_fun(jv,xq(:,iq))
        grad_phi_j = matmul(transpose(invJac(:,:,iq,ie)), grad_lin_shape_fun(jv,xq(:,iq)))

        do iv = 1, 4
          phi_i = lin_shape_fun(iv,xq(:,iq))
          Me(iv,jv,ie) = Me(iv,jv,ie) + phi_i * phi_j * detJac(iq,ie) * wq(iq)

          grad_phi_i = matmul(transpose(invJac(:,:,iq,ie)), grad_lin_shape_fun(iv,xq(:,iq)))
          Ke(iv,jv,ie) = Ke(iv,jv,ie) + dot_product(grad_phi_i, grad_phi_j) * detJac(iq,ie) * wq(iq)
        end do 
      end do 
    end do

    do iv = 1, 4
      pc(element(iv,ie)) = pc(element(iv,ie)) + Ke(iv,iv,ie) 
    end do

  end do

  ! Set initial values
  u = g
  grad_u = 0.d0; grad_u_ = 0.d0
  lambda = 0.d0; lambda_ = 0.d0
  
  open(unit=99,file='res.dat')

  iterADMM = 0
  admmLoop: do
    
    ! 1. Compute q: -tau*(Du - q, p) + -(lambda, p) = 0  <=> tau*(q,p) = (lambda + tau*Du,p) = 0 for all p in Q
    q = grad_u - lambda / tau

    ! 2. Project q
    !$omp parallel do
    do ie = 1, ne
      norm_q(ie) = maxval(abs(q(:,ie)))
      if (norm_q(ie) > 1.d0) q(:,ie) = q(:,ie) / norm_q(ie)
      norm_q1(ie) = maxval(abs(q(:,ie)))
    end do
    !$omp end parallel do

    ! 3. Solve Poisson problem for u: F'(u)v + G'(Du)v + tau*(Du - q,Dv) + (lambda, Dv) = 0 for all v in U
    iterPoisson = 0
    poissonLoop: do 

      fe = 0.d0
      !$omp parallel do reduction(+:fe) private(iq, iv, phi_i, grad_phi_i)
      do ie = 1, ne
        do iq = 1, nq
          do iv = 1, 4
            phi_i = lin_shape_fun(iv,xq(:,iq))
            grad_phi_i = matmul(transpose(invJac(:,:,iq,ie)), grad_lin_shape_fun(iv,xq(:,iq)))
            fe(iv,ie) = fe(iv,ie) + (f(element(iv,ie)) * phi_i&
                                     + dot_product(tau*q(:,ie) + lambda(:,ie), grad_phi_i))* detJac(iq,ie) * wq(iq)
          end do
        end do 
      end do
      !$omp end parallel do

      r = 0.d0
      !$omp parallel do reduction(+:r)
      do ie = 1, ne
        r(element(:,ie)) = r(element(:,ie)) + fe(:,ie) - tau*matmul(Ke(:,:,ie), u(element(:,ie)))
      end do
      !$omp end parallel do
 
      bc = 0.d0
      !$omp parallel do private(iv, jv)
      do ie = 1, ne_bnd
        do iv = 1, 2
          jv = element_bnd(iv,ie)
          if (bctype(jv) == ESSENTIAL) then
            bc(jv) = -r(jv) + (g(jv) - u(jv))
            pc(jv) = 1.d0
          end if
        end do
      end do
      !$omp end parallel do
  
      !$omp parallel do private(iv, jv, phi_i)
      do ie = 1, ne_bnd
        do iv = 1, 2
          jv = element_bnd(iv,ie)
          if (bctype(jv) == NATURAL) then
            do iq = 1, nq_bnd
              phi_i = lin_shape_fun_1d(iv,xq_bnd(iq))
              bc(jv) = bc(jv) + (phi_i * g(jv)) * l_bnd(iq,ie) * wq_bnd(iq)
            end do
          end if
        end do
      end do
      !$omp end parallel do

      !$omp parallel
      !$omp do
      do iv = 1, nv 
        r(iv) = r(iv) + bc(iv)
      end do
      !$omp end do
  
      !$omp do
      do iv = 1, nv 
        u(iv) = u(iv) + r(iv) / pc(iv) * 1.d0
      end do
      !$omp end do
      !$omp end parallel
  
      resPoisson = sqrt(L2InnerProduct(r,r))
      if (iterPoisson == 0) resPoisson0 = resPoisson
  
      if (resPoisson <= resPoisson0 * 1.d-6) exit poissonLoop
      if (iterPoisson >= 100000) exit poissonLoop
      iterPoisson = iterPoisson + 1

    end do poissonLoop
    write(*,'(a,i0,a,e14.7)') 'Poisson iterations ', iterPoisson, ' iterations with residual ', resPoisson / resPoisson0
    
    ! Compute Du
    grad_u = 0.d0
    max_norm_grad_u = 0.d0
    !$omp parallel do private(iv, iq, grad_phi_i) reduction(+:grad_u) reduction(max:max_norm_grad_u)
    do ie = 1, ne
      do iv = 1, 4
        do iq = 1, nq
          grad_phi_i = matmul(transpose(invJac(:,:,iq,ie)), grad_lin_shape_fun(iv,xq(:,iq)))
          grad_u(:,ie) = grad_u(:,ie) + grad_phi_i * u(element(iv,ie)) * detJac(iq,ie) * wq(iq)
        end do
      end do
      grad_u(:,ie) = grad_u(:,ie) / area(ie)

      norm_grad_u(ie) = maxval(abs(grad_u(:,ie)))
      if (max_norm_grad_u < norm_grad_u(ie)) max_norm_grad_u = norm_grad_u(ie)
    end do
    !$omp end parallel do
    
    ! 4. Update multiplier
    !$omp parallel do
    do ie = 1, ne
      lambda(:,ie) = lambda_(:,ie) + tau*(q(:,ie) - grad_u(:,ie))
    end do
    !$omp end parallel do

    resADMM = ADMMresidual(grad_u, grad_u_, lambda, lambda_)
    if (iterADMM == 0) resADMM0 = resADMM

    write(*,'(a,i0,a,e14.7)') 'ADMM iteration ', iterADMM, ' residual ', resADMM / resADMM0
    write(*,'(a,f18.15)') 'Max. norm grad u: ', max_norm_grad_u
    write(99,'(i10,1x,e23.15)') iterADMM, resADMM / resADMM0
    flush(99)

    if (resADMM <= resADMM0 * 1.d-6) exit admmLoop
    if (iterADMM >= 4000) exit admmLoop
    iterADMM = iterADMM + 1

    !$omp parallel do
    do ie = 1, ne
      grad_u_(:,ie) = grad_u(:,ie)
      lambda_(:,ie) = lambda(:,ie)
    end do
    !$omp end parallel do
  
    call write_vtk(iterADMM) 

  end do admmLoop

  close(99)

  call write_vtk(iterADMM) 

  contains

    function lex_map(nx, ny, i, j) result(idx)
      integer, intent(in) :: nx, ny, i, j
      integer :: idx
      idx = (j - 1)*nx + i
    end function lex_map

    function lin_shape_fun_1d(i,x) result(phi)
      integer, intent(in) :: i
      double precision, intent(in) :: x
      double precision :: phi

      if (i == 1) then 
        phi = 0.5d0 * (1.d0 - x)
      else if (i == 2) then
        phi = 0.5d0 * (1.d0 + x)
      else 
        phi = 0.d0
      end if

    end function lin_shape_fun_1d
    
    function grad_lin_shape_fun_1d(i,x) result(grad_phi)
      integer, intent(in) :: i
      double precision, intent(in) :: x
      double precision :: grad_phi

      if (i == 1) then 
        grad_phi = -0.5d0
      else if (i == 2) then
        grad_phi = 0.5d0
      else 
        grad_phi = 0.d0
      end if

    end function grad_lin_shape_fun_1d

    function lin_shape_fun(i,x) result(phi)
      integer, intent(in) :: i
      double precision, intent(in) :: x(2)
      double precision :: phi

      if (i == 1) then
        phi = 0.25d0 * (1.d0 - x(1)) * (1.d0 - x(2))
      else if (i == 2) then
        phi = 0.25d0 * (1.d0 + x(1)) * (1.d0 - x(2))
      else if (i == 3) then 
        phi = 0.25d0 * (1.d0 + x(1)) * (1.d0 + x(2))
      else if (i == 4) then 
        phi = 0.25d0 * (1.d0 - x(1)) * (1.d0 + x(2))
      else 
        phi = 0.d0
      end if
    end function lin_shape_fun

    function grad_lin_shape_fun(i,x) result(grad_phi)
      integer, intent(in) :: i
      double precision, intent(in) :: x(2)
      double precision :: grad_phi(2)

      if (i == 1) then 
        grad_phi(1) = -0.25d0 * (1.d0 - x(2))
        grad_phi(2) = -0.25d0 * (1.d0 - x(1))
      else if (i == 2) then 
        grad_phi(1) =  0.25d0 * (1.d0 - x(2))
        grad_phi(2) = -0.25d0 * (1.d0 + x(1))
      else if (i == 3) then 
        grad_phi(1) =  0.25d0 * (1.d0 + x(2))
        grad_phi(2) =  0.25d0 * (1.d0 + x(1))
      else if (i == 4) then 
        grad_phi(1) = -0.25d0 * (1.d0 + x(2))
        grad_phi(2) =  0.25d0 * (1.d0 - x(1))
      else 
        grad_phi = 0.d0
      end if


    end function grad_lin_shape_fun

    function L2InnerProduct(u,v) result(a)
      double precision, intent(in) :: u(:), v(:)
      double precision :: a, phi_i, phi_j
      integer :: ie, iq, iv, jv

      a = 0.d0
      !$omp parallel do reduction(+:a)
      do ie = 1, ne
        a = a + dot_product(u(element(:,ie)), matmul(Me(:,:,ie), v(element(:,ie))))
      end do
      !$omp end parallel do
    end function L2InnerProduct


    function ADMMresidual(grad_u, grad_u_, lambda, lambda_) result(res)
      double precision, intent(in) :: grad_u(:,:), grad_u_(:,:), lambda(:,:), lambda_(:,:)
      double precision :: res, ru, rl, du(2), dl(2)
      integer :: ie, iq, iv
      ru = 0.d0
      rl = 0.d0
      !$omp parallel do reduction(+:ru, rl)
      do ie = 1, ne
        du = grad_u(:,ie) - grad_u_(:,ie)
        dl = lambda(:,ie) - lambda_(:,ie)
        do iq = 1, nq
          ru = ru + dot_product(du, du) * detJac(iq,ie) * wq(iq)
          rl = rl + dot_product(dl, dl) * detJac(iq,ie) * wq(iq)
        end do
      end do
      !$omp end parallel do
      res = sqrt(ru + rl)
    end function ADMMresidual

    subroutine write_vtk(id)

      integer, intent(in) :: id 
      integer :: iv, ie, file_unit
      character(len=128) :: filename

      write(filename,'(a,i04.4,a)') 'out/out_',id,'.vtu'
      open(newunit=file_unit, file=trim(filename))

      write(file_unit,'(A)') '<?xml version="1.0"?>'
      write(file_unit,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
      write(file_unit,'(A)') ' <UnstructuredGrid>'
      write(file_unit,'(A,I0,A,I0,A)') '  <Piece NumberOfPoints="', nv, '" NumberOfCells="', ne, '">'
      write(file_unit,'(A)') '  <Points>'
      write(file_unit,'(A)') '   <DataArray type="Float64" NumberOfComponents="3" Format="ascii">'
      do iv = 1, nv
        write(file_unit,'(4X,3(E23.16,1X))') coord(:,iv), 0.d0
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '  </Points>'
      write(file_unit,'(A)') '  <Cells>'
      write(file_unit,'(A)') '   <DataArray type="Int32" Name="connectivity" Format="ascii">'
      do ie = 1, ne
        write(file_unit,'(4X,4(I0,1X))') element(1:4,ie)-1
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '   <DataArray type="Int32" Name="offsets" Format="ascii">'
      do ie = 1, ne
        write(file_unit,'(I0)') ie*4
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '   <DataArray type="Int32" Name="types" Format="ascii">'
      do ie = 1, ne
        write(file_unit,'(I0)') 9
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '  </Cells>'

      write(file_unit,'(A)') '  <PointData Scalars="Solution">'
      write(file_unit,'(A)') '   <DataArray Name="u" type="Float64" Format="ascii">'
      do iv = 1, nv
        write(file_unit,'(4X,E23.16)') u(iv)
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '   <DataArray Name="f" type="Float64" Format="ascii">'
      do iv = 1, nv
        write(file_unit,'(4X,E23.16)') f(iv)
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '  </PointData>'

      write(file_unit,'(A)') '  <CellData Scalars="Norm" Vectors="Gradient">'
      write(file_unit,'(A)') '   <DataArray Name="norm_q" type="Float64" Format="ascii">'
      do ie = 1, ne
        write(file_unit,'(4X,E23.16)') norm_q(ie)
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '   <DataArray Name="norm_grad_u" type="Float64" Format="ascii">'
      do ie = 1, ne
        write(file_unit,'(4X,E23.16)') norm_grad_u(ie)
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '   <DataArray Name="norm_q1" type="Float64" Format="ascii">'
      do ie = 1, ne
        write(file_unit,'(4X,E23.16)') norm_q1(ie)
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '   <DataArray Name="grad_u" NumberOfComponents="2" type="Float64" Format="ascii">'
      do ie = 1, ne
        write(file_unit,'(4X,E23.16,1X,E23.16)') grad_u(:,ie)
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '   <DataArray Name="lambda" NumberOfComponents="2" type="Float64" Format="ascii">'
      do ie = 1, ne
        write(file_unit,'(4X,E23.16,1X,E23.16)') lambda(:,ie)
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '  </CellData>'
      write(file_unit,'(A)') '  </Piece>'
      write(file_unit,'(A)') ' </UnstructuredGrid>'
      write(file_unit,'(A)') '</VTKFile>'
      close(file_unit)

    end subroutine write_vtk

end program main
