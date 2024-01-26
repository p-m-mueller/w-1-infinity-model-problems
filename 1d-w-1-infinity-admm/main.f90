program main
    implicit none
    
    double precision, parameter :: pi = 2.d0*acos(0.d0)

    integer, parameter :: nv = 33
    integer, parameter :: ne = nv-1
    integer, parameter :: maxIter = 20000
    double precision, parameter :: tol = 1.d-10
    double precision, parameter :: tau = 1.d0
    double precision, parameter :: xL = 0.d0
    double precision, parameter :: xR = 2.d0*pi
    double precision, parameter :: uL = 0.d0
    double precision, parameter :: uR = 0.d0
    double precision, allocatable, dimension(:) :: x, f, g, bc
    double precision, allocatable, dimension(:) :: u, u_
    double precision, allocatable, dimension(:) :: q, abs_q, du, du_
    double precision, allocatable, dimension(:) :: lambda, lambda_
    integer, allocatable, dimension(:) :: bctype
    integer, parameter :: INTERNAL = 0, ESSENTIAL = 1, NATURAL = 2
    double precision, parameter, dimension(2) :: xe_ref = [-1.d0, 1.d0]
    integer, parameter :: nq = 8
    double precision, parameter :: xq(nq) = [0.96028986, 0.79666648, 0.52553241, 0.18343464, -0.18343464, -0.52553241, -0.79666648,&
                                             -0.96028986]
    double precision, parameter :: wq(nq) = [0.10122854, 0.22238103, 0.31370665, 0.36268378,  0.36268378,  0.31370665,  0.22238103,&
                                              0.10122854]
    double precision :: dx, res, res0, obj, obj0
    integer :: j
    integer :: res_unit, obj_unit

    call init()

    call system("mkdir -p out")
    call system("rm -rf out/steps")
    call system("mkdir -p out/steps")
    open(newunit=res_unit, file='out/res.dat')
    open(newunit=obj_unit, file='out/obj.dat')
    j = 0
    admm : do
        u_ = u
        call solve_u()

        du_ = du
        call compute_du()
        
        call update_q()

        call project_q()

        lambda_ = lambda
        call update_lambda()

        call write_step(j+1)

        res = admmResidual()
        obj = objective()
        if (j == 0) then 
            res0 = res
            obj0 = obj 
        end if

        write(*,'(i10,2(1x,e14.7))') j, res/res0, obj
        write(res_unit,'(i10,1x,e23.16)') j, res/res0
        write(obj_unit,'(i10,1x,e23.16)') j, obj

        if (res <= tol * res0) exit admm

        if (j >= maxIter) exit admm

        j = j + 1

    end do admm
    close(res_unit)
    close(obj_unit) 

    call write_final()

    call cleanup()
   
    contains 

    subroutine cleanup()

        deallocate(x, f, g, bc, u, u_,&
                   q, abs_q, du, du_, lambda, lambda_,&
                   bctype)
        
    end subroutine cleanup

    subroutine init()
        integer :: iv, ie, element(2)

        allocate(x(nv), f(nv), g(nv), bc(nv), u(nv), u_(nv),&
                 q(ne), abs_q(ne), du(ne), du_(ne), lambda(ne), lambda_(ne),&
                 bctype(nv))

        dx = (xR - xL) / ne

        u = 0.d0
        u_ = 0.d0
        q = 0.d0
        du = 0.d0
        lambda = 0.d0
        lambda_ = 0.d0
        
        do iv = 1, nv
            x(iv) = xL + (iv-1)*dx
            if (x(iv) <= pi) then
                f(iv) = 0.1d0
            else
                f(iv) = -0.1d0
            end if
           g(iv) = -0.1d0*abs(x(iv)-pi)
        end do
        
        bctype = INTERNAL
        bctype(1) = ESSENTIAL
        bctype(nv) = NATURAL

        bc = 0.d0 
        bc(1) = uL
        bc(nv) = uR

    end subroutine init

    subroutine update_q() 
        
        q = du - lambda / tau

    end subroutine update_q 

    subroutine project_q()
        integer :: ie

        do ie = 1, ne 
            abs_q(ie) = abs(q(ie))
            if (abs_q(ie) > 1.d0) q(ie) = q(ie) / abs_q(ie)
        end do

    end subroutine project_q

    subroutine solve_u()
        integer :: ie, iq, iv, jv
        integer :: element(2)
        double precision :: Ae(2,2), be(2)
        double precision, allocatable :: A(:,:), B(:)
        integer, allocatable :: ipiv(:)
        integer :: idof, ndof, info
        double precision :: dphi_i, dphi_j, phi_i, phi_j

        ndof = nv+1
        allocate(A(ndof,ndof), B(ndof), ipiv(ndof))

        A = 0.d0
        B = 0.d0
        ! Assemble system
        do ie = 1, ne
            
            element = [ie, ie+1]

            Ae = 0.d0
            be = 0.d0
            do iq = 1, nq
                do iv = 1, 2

                    phi_i = phi(iv,xq(iq))
                    dphi_i = dphi(iv,xq(iq)) / dx

                    do jv = 1, 2

                        phi_j = phi(jv,xq(iq))
                        dphi_j = dphi(jv,xq(iq)) / dx

                        Ae(iv,jv) = Ae(iv,jv) - tau * dphi_i * dphi_j * dx * wq(iq)
                        be(iv) = be(iv) + (f(element(jv)) * phi_j * phi_i + g(element(jv)) * phi_j * dphi_i )* dx * wq(iq)
                    end do
                    be(iv) = be(iv) - (tau * q(ie) + lambda(ie)) * dphi_i * dx * wq(iq)
                end do
            end do

            A(element,element) = A(element,element) + Ae
            B(element) = B(element) + be

            do iv = 1, 2
                jv = element(iv)
                do iq = 1, nq 
                    phi_i = phi(iv,xq(iq))
                    A(jv,nv+1) = A(jv,nv+1) + phi_i * dx * wq(iq)
                end do
            end do
        end do
        A(nv+1,:) = A(:,nv+1)

        if (bctype(1) == NATURAL) B(1) = B(1) + bc(1)
        if (bctype(nv) == NATURAL) B(nv) = B(nv) + bc(nv)

        do ie = 1, ne
            element = [ie, ie+1]
            do iv = 1, 2
                if (bctype(element(iv)) == ESSENTIAL) then
                    idof = element(iv)
                    A(idof,:) = 0.d0
                    A(:,idof) = 0.d0
                    A(idof,idof) = 1.d0
                    B(idof) = bc(idof)
                end if
            end do
        end do

        call dgesv(ndof, 1, A, ndof, ipiv, B, ndof, info)
        if (info /= 0) then 
            if (info < 0) then 
                write(*,'(a,i0)') 'Error: LAPACK routine DGESV: Illigal input for argument ', -info
            else 
                write(*,'(a,i0)') 'Error: LAPACK routine DGESV: Matrix close to singular at entry ', info 
            end if
        else 
            u = B(1:nv)
        end if

        deallocate(A,b)

    end subroutine solve_u 

    subroutine compute_du()
        integer :: ie, iq, iv, element(2)
        double precision :: dphi_i, len

        du = 0.d0 
        do ie = 1, ne 
            element = [ie, ie+1]
            len = 0.d0
            do iq = 1, nq
                do iv = 1, 2 
                    dphi_i = dphi(iv,xq(iq)) / dx
                    du(ie) = du(ie) + u(element(iv)) * dphi_i * dx * wq(iq)
                end do
                len = len + dx * wq(iq)
            end do
            du(ie) = du(ie) / len
        end do
    end subroutine compute_du

    subroutine update_lambda()

        lambda = lambda_ + tau * (q - du)

    end subroutine update_lambda

    double precision function admmResidual()
        integer :: ie, iv, iq, element(2)

        admmResidual = 0.d0 
        do ie = 1, ne
            element = [ie, ie+1] 
            do iv = 1, 2
                do iq = 1, nq
                    admmResidual = admmResidual + ((lambda(ie) - lambda_(ie))**2 + (du(ie) - du_(ie))**2) * dx * wq(iq)
                end do
            end do
        end do
        admmResidual = sqrt(admmResidual)
    end function admmResidual


    double precision function objective() 
        integer :: ie, iv, jv, iq, element(2)
        double precision :: phi_i, phi_j, dphi_i, dphi_j, du_q(nq), u_q(nq), f_q(nq), g_q(nq)

        objective = 0.d0 
        do ie = 1, ne 
            element = [ie, ie+1]
            do iq = 1, nq
                u_q(iq) = 0.d0
                du_q(iq) = 0.d0 
                f_q(iq) = 0.d0
                g_q(iq) = 0.d0
                do iv = 1, 2
                  phi_i = phi(iv,xq(iq))
                  dphi_i = dphi(iv,xq(iq))
                  u_q(iq) = u_q(iq) + u(element(iv)) * phi_i
                  du_q(iq) = du_q(iq) + u(element(iv)) * dphi_i
                  f_q(iq) = f_q(iq) + f(element(iv)) * phi_i
                  g_q(iq) = g_q(iq) + g(iv) * phi_i
                end do
            end do

            do iq = 1, nq
                do iv = 1, 2
                    phi_i = phi(iv,xq(iq))
                    dphi_i = dphi(iv,xq(iq)) / dx
                    objective = objective + (f_q(iq)*u_q(iq) + g_q(iq)*du_q(iq)) * dx * wq(iq)
                end do
            end do

        end do
    end function objective


    subroutine write_final()
        integer :: file_unit, iv, ie, element(2)

        open(newunit=file_unit, file='out/point_data.dat')
        write(file_unit,'(4(1x,a23))') 'x', 'u', 'f', 'g'
        do iv = 1, nv
            write(file_unit,'(4(1x,e23.16))') x(iv), u(iv), f(iv), g(iv)
        end do
        close(file_unit)

        open(newunit=file_unit, file='out/element_data.dat')
        write(file_unit,'(4(1x,a23))') 'x', 'du', 'q'
        do ie = 1, ne
            element = [ie, ie+1]
            do iv = 1, 2
                write(file_unit,'(4(1x,e23.16))') x(element(iv)), du(ie), q(ie)
            end do
        end do
        close(file_unit)

    end subroutine write_final

    subroutine write_step(j)
        integer, intent(in) :: j
        integer :: file_unit, iv, ie, element(2)
        character(len=128) :: filename

        write(filename,'(a,i0.4,a)') 'out/steps/point_data_',j,'.dat'
        open(newunit=file_unit, file=trim(filename))
        write(file_unit,'(4(1x,a23))') 'x', 'u', 'f', 'g'
        do iv = 1, nv
            write(file_unit,'(4(1x,e23.16))') x(iv), u(iv), f(iv), g(iv)
        end do
        close(file_unit)

        write(filename,'(a,i0.4,a)') 'out/steps/element_data_',j,'.dat'
        open(newunit=file_unit, file=trim(filename))
        write(file_unit,'(4(1x,a23))') 'x', 'du', 'q', 'lambda'
        do ie = 1, ne
            element = [ie, ie+1]
            write(file_unit,'(4(1x,e23.16))') 0.5d0*sum(x(element)), du(ie), q(ie), lambda(ie)
        end do
        close(file_unit)

    end subroutine write_step

    

    double precision function phi(i,x)
        integer, intent(in) :: i 
        double precision, intent(in) :: x

        if (i == 1) then 
            phi = 0.5d0 * (1.d0 - x)
        else if (i == 2) then 
            phi = 0.5d0 * (1.d0 + x)
        else 
            phi = 0.d0
        end if

    end function phi
    
    double precision function dphi(i,x)
        integer, intent(in) :: i 
        double precision, intent(in) :: x

        if (i == 1) then 
            dphi = -0.5d0 
        else if (i == 2) then 
            dphi = 0.5d0
        else 
            dphi = 0.d0
        end if

    end function dphi

end program main
