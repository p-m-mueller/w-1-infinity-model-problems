program main

  use iso_fortran_env
  use quadrature

  implicit none

  double precision, parameter :: pi = 2.d0*asin(1.d0)
  double precision, parameter :: eps = 1.d-8

  integer, parameter :: nx = 33
  integer, parameter :: ny = 33
  integer, parameter :: nv = nx * ny
  integer, parameter :: ne = (nx - 1) * (ny - 1)
  double precision, parameter :: p_min = 2.d0
  double precision, parameter :: p_max = 8.0d0
  double precision, parameter :: p_inc = (p_max - p_min) / 8.d0
  integer :: cnt

  integer, parameter :: nq_1d = 5
  double precision :: xq_1d(nq_1d), wq_1d(nq_1d)
  integer, parameter :: nq = nq_1d*nq_1d
  double precision :: xq(2,nq), wq(nq)
  
  integer :: element(4,ne)
  double precision :: Jac(2,2,nq,ne), invJac(2,2,nq,ne), detJac(nq,ne), area(ne)

  integer, parameter :: ne_bnd = 2 * (nx - 1 + ny - 1)

  integer, parameter :: nq_bnd = nq_1d
  integer :: element_bnd(2,ne_bnd)
  double precision :: l_bnd(nq_bnd,ne_bnd)
  double precision :: xq_bnd(nq_bnd), wq_bnd(nq_bnd)

  double precision :: coord(2,nv)
  double precision :: dx, dy

  double precision :: u(nv), old_u(nv) ! solution
  double precision :: delta_u(nv)
  double precision :: grad_u(2,ne), norm_grad_u(ne)
  
  double precision :: Fu, old_Fu 
  double precision :: alpha = 1.d0

  double precision :: f(nv)
  double precision :: g(nv) ! boundary data
  integer :: bctype(nv)
  integer, parameter :: INTERNAL = 0, ESSENTIAL = 1, NATURAL = 2

  double precision :: grad_phi_i(2), grad_phi_j(2)
  double precision :: phi_i, phi_j

  double precision :: grad_u_q(2,nq), u_q(nq), eta_q(nq), f_q(nq)
  double precision :: p

  double precision :: Ke(4,4,ne) ! element stiffnes matrix
  double precision :: fe(4,ne)  ! element right hand side
  double precision :: bc(nv) ! right hand side modification for boundary conditions
  double precision :: pc(nv) ! preconditioner for linear solver
  double precision :: r(nv)
  double precision :: resL1, resL, resNL1, resNL ! residual
  
  integer :: i, j, k, iv, jv, ie, iq
  integer :: iterL, iterNL
  
  write(*,'(a,i0)') 'Number of elements: ', ne
  write(*,'(a,i0)') 'Number of vertices: ', nv
  write(*,'(a,i0)') 'Number of unknowns: ', nv


  WRITE(*,'(a)') 'Compute quadrature points and weights for boundary'
  call points(nq_bnd,xq_bnd)
  call weights(nq_bnd,xq_bnd,wq_bnd)
  
  WRITE(*,'(a)') 'Compute quadrature points and weights for domain'
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

  u = 0.d0
  do ie = 1, ne_bnd
    do iv = 1, 2
      jv = element_bnd(iv,ie)
      if (bctype(jv) == ESSENTIAL) u(jv) = g(jv)
    end do
  end do

  Fu = 0.d0
  old_Fu = Fu 

  open(unit=99,file='res.dat')
 
  p = p_min
  cnt = 1
  pLoop: do

  write(*,'("# ",A,F6.3)') 'p: ', p
  
  nonlinearSolverLoop: do iterNL = 1, 100

    ! assemble operator, right hand side and preconditioner
    pc = 0.d0
    do ie = 1, ne
  
        Ke(:,:,ie) = 0.d0
        fe(:,ie) = 0.d0
        do iq = 1, nq
          grad_u_q(:,iq) = 0.d0
          u_q(iq) = 0.d0
          f_q(iq) = 0.d0
          do iv = 1, 4
            grad_phi_i = matmul(transpose(invJac(:,:,iq,ie)), grad_lin_shape_fun(iv, xq(:,iq)))
            phi_i = lin_shape_fun(iv,xq(:,iq))

            grad_u_q(:,iq) = grad_u_q(:,iq) + grad_phi_i * u(element(iv,ie))
            u_q(iq) = u_q(iq) + phi_i * u(element(iv,ie))
            f_q(iq) = f_q(iq) + phi_i * f(element(iv,ie))
          end do
          eta_q(iq) = dot_product(grad_u_q(:,iq), grad_u_q(:,iq))
          if (p < 4.d0) eta_q(iq) = eta_q(iq) + 1.e-6
        end do

        do iq = 1, nq
          do jv = 1, 4
            
            phi_j = lin_shape_fun(jv,xq(:,iq))
            grad_phi_j = matmul(transpose(invJac(:,:,iq,ie)), grad_lin_shape_fun(jv,xq(:,iq)))

            do iv = 1, 4
              grad_phi_i = matmul(transpose(invJac(:,:,iq,ie)), grad_lin_shape_fun(iv,xq(:,iq)))
              phi_i = lin_shape_fun(iv,xq(:,iq))

              if (abs(p - 2.d0) < eps) then 
                Ke(iv,jv,ie) = Ke(iv,jv,ie) + dot_product(grad_phi_i, grad_phi_j) * detJac(iq,ie) * wq(iq)
              else
                Ke(iv,jv,ie) = Ke(iv,jv,ie)&
                               + ( (p-2.d0) * eta_q(iq)**(0.5d0*(p-4.d0)) * dot_product(grad_u_q(:,iq), grad_phi_i)&
                               * dot_product(grad_u_q(:,iq), grad_phi_j) + eta_q(iq)**(0.5d0*(p-2.d0))&
                               * dot_product(grad_phi_i, grad_phi_j) ) * detJac(iq,ie) * wq(iq)
              end if
            end do
          end do 
        end do

        do iq = 1, nq
          do iv = 1, 4
            phi_i = lin_shape_fun(iv,xq(:,iq))
            grad_phi_i = matmul(transpose(invJac(:,:,iq,ie)), grad_lin_shape_fun(iv,xq(:,iq)))

            if (abs(p - 2.d0) < eps) then 
              fe(iv,ie) = fe(iv,ie) - (dot_product(grad_u_q(:,iq), grad_phi_i) - f_q(iq) * phi_i) * detJac(iq,ie) * wq(iq)
            else
              fe(iv,ie) = fe(iv,ie) - ( eta_q(iq)**(0.5d0*(p-2.d0)) * dot_product(grad_u_q(:,iq), grad_phi_i)&
                                    - f_q(iq) * phi_i) * detJac(iq,ie) * wq(iq)
            end if
          end do
        end do

        do iv = 1, 4
          pc(element(iv,ie)) = pc(element(iv,ie)) + Ke(iv, iv,ie)
        end do
    end do

    delta_u = 0.d0
    linearSolverLoop: do iterL = 1, 100000 
    
      r = 0.d0
      !$omp parallel do reduction(+:r)
      do ie = 1, ne
        r(element(:,ie)) = r(element(:,ie)) + fe(:,ie) - matmul(Ke(:,:,ie), delta_u(element(:,ie)))
      end do
      !$omp end parallel do
 
      bc = 0.d0
      !$omp parallel do private(iv, jv)
      do ie = 1, ne_bnd
        do iv = 1, 2
          jv = element_bnd(iv,ie)
          if (bctype(jv) == ESSENTIAL) then 
            bc(jv) = -r(jv) + (0.d0 - delta_u(jv))
            pc(jv) = 1.d0
          end if
        end do
      end do
      !$omp end parallel do

      !$omp parallel do reduction(+:bc) private(iv, jv)
      do ie = 1, ne_bnd
        do iv = 1, 2
          jv = element_bnd(iv,ie)
          if (bctype(jv) == NATURAL) then 
            do iq = 1, nq_bnd
              phi_i = lin_shape_fun_1d(iv, xq_bnd(iq))
              bc(jv) = bc(jv) + (phi_i * g(jv)) * l_bnd(iq,ie) * wq_bnd(iq)
            end do
          end if
        end do
      end do
      !$omp end parallel do
 
      !$omp parallel do reduction(+:r)
      do iv = 1, nv 
        r(iv) = r(iv) + bc(iv)
      end do
      !$omp end parallel do
   
      !$omp parallel do reduction(+:delta_u)
      do iv = 1, nv 
        delta_u(iv) = delta_u(iv) + r(iv) / pc(iv) * 0.8d0
      end do
      !$omp end parallel do
    
      ! compute L2-norm of the residual
      resL = 0.d0; resNL = 0.d0
      !$omp parallel do reduction(+:resL, resNL) private(iq, iv)
      do ie = 1, ne
        do iq = 1, nq
          do iv = 1, 4
            phi_i = lin_shape_fun(iv, xq(:,iq))
            resL = resL + (r(element(iv,ie)) * phi_i)**2 * detJac(iq,ie) * wq(iq)
            resNL = resNL + (alpha*delta_u(element(iv,ie)) * phi_i)**2 * detJac(iq,ie) * wq(iq)
          end do
        end do
      end do
      !$omp end parallel do
      resL = sqrt(resL); resNL = sqrt(resNL)
      if (iterL == 1) resL1 = resL

      if (resL <= 1.e-7 * resL1 .or. resL > 1.e6) exit linearSolverLoop
    
    end do linearSolverLoop 

    if (iterNL == 1) resNL1 = resNL

    write(*,'(2(I10,1X),3(E16.8,1X))') iterNL, iterL, resL / resL1, resNL / resNL1, alpha
    write(99,'(I10,2(1X,E23.15),1X,F7.5)') iterNL, resNL / resNL1, alpha, p
    flush(99)

    if (resNL < 1.e-6 * resNL1 .or. resNL > 1.e2) exit nonlinearSolverLoop

    alpha = 1.d0
    lineSearchLoop: do k = 1, 100
      !$omp parallel do
      do iv = 1, nv
        u(iv) = old_u(iv) + delta_u(iv) * alpha
      end do
      !$omp end parallel do
    
      Fu = 0.d0
      !$omp parallel do reduction(+:Fu) private(iq, iv, grad_u_q, u_q, f_q, grad_phi_i, phi_i, eta_q)
      do ie = 1, ne
        do iq = 1, nq
          grad_u_q(:,iq) = 0.d0
          u_q(iq) = 0.d0
          f_q(iq) = 0.d0
          do iv = 1, 4
            grad_phi_i = matmul(transpose(invJac(:,:,iq,ie)), grad_lin_shape_fun(iv, xq(:,iq)))
            phi_i = lin_shape_fun(iv, xq(:,iq))

            grad_u_q(:,iq) = grad_u_q(:,iq) + grad_phi_i * u(element(iv,ie)) 
            u_q(iq) = u_q(iq) + phi_i * u(element(iv,ie))
            f_q(iq) = f_q(iq) + phi_i * f(element(iv,ie))
          end do
          eta_q(iq) = dot_product(grad_u_q(:,iq), grad_u_q(:,iq))
          if (p < 4.d0) eta_q(iq) = eta_q(iq) + 1.d-6
        end do

        do iq = 1, nq
          Fu = Fu + (1.d0/p*eta_q(iq)**(0.5d0*p) - f_q(iq) * u_q(iq)) * detJac(iq,ie) * wq(iq)
        end do
      end do
      !$omp end parallel do

      !$omp parallel do reduction(+:Fu) private(iv, jv, iq, phi_i)
      do ie = 1, ne_bnd
        do iv = 1, 2
          jv = element_bnd(iv,ie)
          if (bctype(jv) == NATURAL) then 
            do iq = 1, nq_bnd
              phi_i = lin_shape_fun_1d(iv, xq_bnd(iq))
              Fu = Fu - (u(jv) * phi_i * g(jv)) * l_bnd(iq,ie) * wq_bnd(iq)
            end do
          end if
        end do
      end do
      !$omp end parallel do

      if (Fu < old_Fu) exit lineSearchLoop

      alpha = 0.5d0 * alpha 
      !write(*,'("#",I5,A,E16.8)') k, ' Line search faile. New step size: ', alpha

      if (alpha < 1.e-8) exit pLoop 
      
    end do lineSearchLoop

    grad_u = 0.d0
    
    call write_vtk(cnt) 

    !$omp parallel do
    do iv = 1, nv
      old_u(iv) = u(iv)
    end do
    !$omp end parallel do
    old_Fu = Fu

  end do nonlinearSolverLoop

    write(99,*)
  
    call write_vtk(cnt) 
    cnt = cnt+1

    p = p + p_inc
    if (p > p_max) exit pLoop

  end do pLoop

  close(99)

!  call write_vtk(cnt) 

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


    subroutine write_vtk(id)
      integer, intent(in) :: id
      integer :: iv, ie, file_unit
      character(len=128) :: filename

      write(filename,'(A,I04.4,a)') 'out/out_', id, '.vtu'
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
      write(file_unit,'(A)') '  <PointData Scalars="Solution" Vectors="Gradient">'
      write(file_unit,'(A)') '   <DataArray Name="u" type="Float64" Format="ascii">'
      do iv = 1, nv
        write(file_unit,'(4X,E23.16)') u(iv)
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '   <DataArray Name="delta_u" type="Float64" Format="ascii">'
      do iv = 1, nv
        write(file_unit,'(4X,E23.16)') delta_u(iv)
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '  </PointData>'
      write(file_unit,'(A)') '  <CellData Scalars="Norm" Vectors="Gradient">'
      write(file_unit,'(A)') '   <DataArray Name="grad_u" NumberOfComponents="2" type="Float64" Format="ascii">'
      do ie = 1, ne
        write(file_unit,'(4X,E23.16,1X,E23.16)') grad_u(:,ie)
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '   <DataArray Name="norm_grad_u" type="Float64" Format="ascii">'
      do ie = 1, ne
        write(file_unit,'(4X,E23.16)') norm_grad_u(ie)
      end do
      write(file_unit,'(A)') '   </DataArray>'
      write(file_unit,'(A)') '  </CellData>'
      write(file_unit,'(A)') '  </Piece>u'
      write(file_unit,'(A)') ' </UnstructuredGrid>'
      write(file_unit,'(A)') '</VTKFile>'
      close(file_unit)

    end subroutine write_vtk

end program main
