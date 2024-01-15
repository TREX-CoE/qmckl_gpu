module qmckl_gpu_f
    use, intrinsic :: iso_c_binding
    integer, parameter :: qmckl_context_device = c_int64_t
    integer, parameter :: qmckl_exit_code_device = c_int32_t
    integer(qmckl_exit_code_device), parameter :: QMCKL_SUCCESS_DEVICE  = 0
    interface

    !!!!!!!!!!!
        ! CONTEXT
    !!!!!!!!!!!

        integer(qmckl_context_device) function qmckl_context_touch_device(context) &
            bind(C, name="qmckl_context_touch_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
        end function qmckl_context_touch_device

        integer(qmckl_context_device) function qmckl_context_create_device(device_id) &
            bind(C, name="qmckl_context_create_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(c_int32_t), intent(in), value :: device_id
        end function qmckl_context_create_device

        integer(qmckl_exit_code_device) function qmckl_context_destroy_device(context) &
            bind(C, name="qmckl_context_destroy_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
        end function qmckl_context_destroy_device

    !!!!!!!!!!!
        ! MEMORY
    !!!!!!!!!!!

        type(c_ptr) function qmckl_malloc_device(context, size) &
            bind(C, name="qmckl_malloc_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: size
        end function qmckl_malloc_device

        type(c_ptr) function qmckl_malloc_host(context, size) &
            bind(C, name="qmckl_malloc_host")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: size
        end function qmckl_malloc_host


        integer(qmckl_exit_code_device) function qmckl_free_device(context, ptr) &
            bind(C, name="qmckl_free_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            type(c_ptr), intent(in), value :: ptr
        end function qmckl_free_device

        integer(qmckl_exit_code_device) function qmckl_free_host(context, ptr) &
            bind(C, name="qmckl_free_host")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            type(c_ptr), intent(in), value :: ptr
        end function qmckl_free_host


        integer(qmckl_exit_code_device) function qmckl_memcpy_H2D(context, dest, src, size) &
            bind(C, name="qmckl_memcpy_H2D")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            type(c_ptr), intent(in), value :: dest
            type(c_ptr), intent(in), value :: src
            integer(c_int64_t), intent(in), value :: size
        end function qmckl_memcpy_H2D

        integer(qmckl_exit_code_device) function qmckl_memcpy_H2D_double(context, dest, src, size) &
            bind(C, name="qmckl_memcpy_H2D_double")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            type(c_ptr), intent(in), value :: dest
            double precision, intent(in) :: src(*)
            integer(c_int64_t), intent(in), value :: size
        end function qmckl_memcpy_H2D_double

        integer(qmckl_exit_code_device) function qmckl_memcpy_H2D_int32(context, dest, src, size) &
            bind(C, name="qmckl_memcpy_H2D_int32")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            type(c_ptr), intent(in), value :: dest
            integer(c_int32_t), intent(in) :: src(*)
            integer(c_int64_t), intent(in), value :: size
        end function qmckl_memcpy_H2D_int32

        integer(qmckl_exit_code_device) function qmckl_memcpy_H2D_int64(context, dest, src, size) &
            bind(C, name="qmckl_memcpy_H2D_int64")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            type(c_ptr), intent(in), value :: dest
            integer(c_int64_t), intent(in) :: src(*)
            integer(c_int64_t), intent(in), value :: size
        end function qmckl_memcpy_H2D_int64

        integer(qmckl_exit_code_device) function qmckl_memcpy_D2H(context, dest, src, size) &
            bind(C, name="qmckl_memcpy_D2H")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            type(c_ptr), intent(in), value :: dest
            type(c_ptr), intent(in), value :: src
            integer(c_int64_t), intent(in), value :: size
        end function qmckl_memcpy_D2H

        integer(qmckl_exit_code_device) function qmckl_memcpy_D2H_double(context, dest, src, size) &
            bind(C, name="qmckl_memcpy_D2H_double")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            double precision, intent(in) :: dest(*)
            type(c_ptr), intent(in), value :: src
            integer(c_int64_t), intent(in), value :: size
        end function qmckl_memcpy_D2H_double

        integer(qmckl_exit_code_device) function qmckl_memcpy_D2H_int32(context, dest, src, size) &
            bind(C, name="qmckl_memcpy_D2H_int32")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            integer(c_int32_t), intent(in) :: dest(*)
            type(c_ptr), intent(in), value :: src
            integer(c_int64_t), intent(in), value :: size
        end function qmckl_memcpy_D2H_int32

        integer(qmckl_exit_code_device) function qmckl_memcpy_D2H_int64(context, dest, src, size) &
            bind(C, name="qmckl_memcpy_D2H_int64")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            integer(c_int64_t), intent(in) :: dest(*)
            type(c_ptr), intent(in), value :: src
            integer(c_int64_t), intent(in), value :: size
        end function qmckl_memcpy_D2H_int64

        integer(qmckl_exit_code_device) function qmckl_memcpy_D2D(context, dest, src, size) &
            bind(C, name="qmckl_memcpy_D2D")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            type(c_ptr), intent(in), value :: dest
            type(c_ptr), intent(in), value :: src
            integer(c_size_t), value :: size
        end function qmckl_memcpy_D2D

    !!!!!!!!!!!
        ! TREXIO
    !!!!!!!!!!!

        integer(qmckl_exit_code_device) function qmckl_trexio_read_device &
            (context, file_name, size_max) &
            bind(C, name="qmckl_trexio_read_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character(c_char), intent(in), value :: file_name
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_trexio_read_device

    !!!!!!!!!!!
        ! NUCLEUS
    !!!!!!!!!!!

        ! Setters

        integer(qmckl_exit_code_device) function qmckl_set_nucleus_num_device(context, num) &
            bind(C, name="qmckl_set_nucleus_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: num
        end function qmckl_set_nucleus_num_device

        integer(qmckl_exit_code_device) function qmckl_set_nucleus_coord_device(context, transp, nucl_coord, size_max) &
            bind(C, name="qmckl_set_nucleus_coord_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character(c_char), intent(in), value :: transp
            real(c_double), intent(in), value :: nucl_coord
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_set_nucleus_coord_device

        integer(qmckl_exit_code_device) function qmckl_set_nucleus_charge_device(context, nucl_charge, nucl_num) &
            bind(C, name="qmckl_set_nucleus_charge_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: nucl_charge
            integer(c_int64_t), intent(in), value :: nucl_num
        end function qmckl_set_nucleus_charge_device

        ! Getters

        integer(qmckl_exit_code_device) function qmckl_get_nucleus_num_device(context, num) &
            bind(C, name="qmckl_get_nucleus_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(out), value :: num ! Elements of type int64_t
        end function qmckl_get_nucleus_num_device

        integer(qmckl_exit_code_device) function qmckl_get_nucleus_coord_device(context, transp, nucl_coord, size_max) &
            bind(C, name="qmckl_get_nucleus_coord_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character(c_char), intent(in), value :: transp
            real(c_double), intent(out), value :: nucl_coord ! Elements of type double*
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_nucleus_coord_device

        integer(qmckl_exit_code_device) function qmckl_get_nucleus_charge_device(context, nucl_charge, nucl_num) &
            bind(C, name="qmckl_get_nucleus_charge_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            type(c_ptr), intent(out), value :: nucl_charge ! Elements of type double
            integer(c_int64_t), intent(in), value :: nucl_num
        end function qmckl_get_nucleus_charge_device

    !!!!!!!!!!!
        ! POINT
    !!!!!!!!!!!

        integer(qmckl_exit_code_device) function qmckl_set_point_device(context, transp, num, coord, size_max) &
            bind(C, name="qmckl_set_point_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character(c_char), intent(in), value :: transp
            integer(c_int64_t), intent(in), value :: num
            real(c_double), intent(in), value :: coord ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_set_point_device

        integer(qmckl_exit_code_device) function qmckl_set_point_device_from_host(context, transp, num, coord, size_max) &
            bind(C, name="qmckl_set_point_device_from_host")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character(c_char), intent(in), value :: transp
            integer(c_int64_t), intent(in), value :: num
            double precision, intent(in) :: coord(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_set_point_device_from_host

        integer(qmckl_exit_code_device) function qmckl_get_point_device(context, transp, coord, size_max) &
            bind(C, name="qmckl_get_point_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character(c_char), intent(in), value :: transp
            real(c_double), intent(out), value :: coord ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_point_device

        integer(qmckl_exit_code_device) function qmckl_get_point_device_to_host(context, transp, coord, size_max) &
            bind(C, name="qmckl_get_point_device_to_host")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character(c_char), intent(in), value :: transp
            double precision, intent(out) :: coord(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_point_device_to_host

    !!!!!!!!!!!
        ! ELECTRON
    !!!!!!!!!!!

        ! Setters
        integer(qmckl_exit_code_device) function qmckl_set_electron_num_device(context, elec_up_num, elec_dn_num) &
            bind(C, name="qmckl_set_electron_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: elec_up_num
            integer(c_int64_t), intent(in), value :: elec_dn_num
        end function qmckl_set_electron_num_device

        integer(qmckl_exit_code_device) function qmckl_set_electron_coord_device(context, transp, walk_num, coord, size_max) &
            bind(C, name="qmckl_set_electron_coord_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character, intent(in), value :: transp
            integer(c_int64_t), intent(in), value :: walk_num
            real(c_double), intent(in), value :: coord
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_set_electron_coord_device

        integer(qmckl_exit_code_device) function qmckl_set_electron_coord_device_from_host(context, transp, walk_num, coord, size_max) &
            bind(C, name="qmckl_set_electron_coord_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character, intent(in), value :: transp
            integer(c_int64_t), intent(in), value :: walk_num
            double precision, intent(in) :: coord(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_set_electron_coord_device_from_host


        ! Getters

        integer(qmckl_exit_code_device) function qmckl_get_electron_coord_device(context, transp, coord, size_max) &
            bind(C, name="qmckl_get_electron_coord_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character, intent(in), value :: transp
            real(c_double), intent(out), value :: coord ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_electron_coord_device

    !!!!!!!!!!!
        ! AO
    !!!!!!!!!!!

        ! Basis setters

        integer(qmckl_exit_code_device) function qmckl_set_ao_basis_ao_num_device(context, num) &
            bind(C, name="qmckl_set_ao_basis_ao_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: num
        end function qmckl_set_ao_basis_ao_num_device

        integer(qmckl_exit_code_device) function qmckl_set_ao_basis_type_device(context, type) &
            bind(C, name="qmckl_set_ao_basis_type_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character(c_signed_char), intent(in), value :: type
        end function qmckl_set_ao_basis_type_device

        integer(qmckl_exit_code_device) function qmckl_set_ao_basis_shell_num_device(context, shell_num) &
            bind(C, name="qmckl_set_ao_basis_shell_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: shell_num
        end function qmckl_set_ao_basis_shell_num_device

        integer(qmckl_exit_code_device) function qmckl_set_ao_basis_prim_num_device(context, prim_num) &
            bind(C, name="qmckl_set_ao_basis_prim_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: prim_num
        end function qmckl_set_ao_basis_prim_num_device

        integer(qmckl_exit_code_device) function qmckl_set_ao_basis_nucleus_index_device(context, nucleus_index, size_max) &
            bind(C, name="qmckl_set_ao_basis_nucleus_index_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: nucleus_index ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_set_ao_basis_nucleus_index_device

        integer(qmckl_exit_code_device) function qmckl_set_ao_basis_nucleus_shell_num_device(context, nucleus_shell_num, nucl_num) &
            bind(C, name="qmckl_set_ao_basis_nucleus_shell_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: nucleus_shell_num ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: nucl_num
        end function qmckl_set_ao_basis_nucleus_shell_num_device

        integer(qmckl_exit_code_device) function qmckl_set_ao_basis_shell_ang_mom_device(context, shell_ang_mom, shell_num) &
            bind(C, name="qmckl_set_ao_basis_shell_ang_mom_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: shell_ang_mom ! Elements of type int32_t
            integer(c_int64_t), intent(in), value :: shell_num
        end function qmckl_set_ao_basis_shell_ang_mom_device

        integer(qmckl_exit_code_device) function qmckl_set_ao_basis_shell_factor_device(context, shell_factor, shell_num) &
            bind(C, name="qmckl_set_ao_basis_shell_factor_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: shell_factor ! Elements of type double
            integer(c_int64_t), intent(in), value :: shell_num
        end function qmckl_set_ao_basis_shell_factor_device

        integer(qmckl_exit_code_device) function qmckl_set_ao_basis_shell_prim_num_device(context, shell_prim_num, shell_num) &
            bind(C, name="qmckl_set_ao_basis_shell_prim_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: shell_prim_num ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: shell_num
        end function qmckl_set_ao_basis_shell_prim_num_device

        integer(qmckl_exit_code_device) function qmckl_set_ao_basis_shell_prim_index_device(context, shell_prim_index, shell_num) &
            bind(C, name="qmckl_set_ao_basis_shell_prim_index_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: shell_prim_index ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: shell_num
        end function qmckl_set_ao_basis_shell_prim_index_device

        integer(qmckl_exit_code_device) function qmckl_set_ao_basis_exponent_device(context, exponent, prim_num) &
            bind(C, name="qmckl_set_ao_basis_exponent_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: exponent ! Elements of type double
            integer(c_int64_t), intent(in), value :: prim_num
        end function qmckl_set_ao_basis_exponent_device

        integer(qmckl_exit_code_device) function qmckl_set_ao_basis_coefficient_device(context, coefficient, prim_num) &
            bind(C, name="qmckl_set_ao_basis_coefficient_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(double), intent(in), value :: coefficient ! Elements of type double
            integer(c_int64_t), intent(in), value :: prim_num
        end function qmckl_set_ao_basis_coefficient_device

        integer(qmckl_exit_code_device) function qmckl_set_ao_basis_prim_factor_device(context, prim_factor, prim_num) &
            bind(C, name="qmckl_set_ao_basis_prim_factor_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: prim_factor ! Elements of type double
            integer(c_int64_t), intent(in), value :: prim_num
        end function qmckl_set_ao_basis_prim_factor_device

        integer(qmckl_exit_code_device) function qmckl_set_ao_basis_ao_factor_device(context, ao_factor, ao_num) &
            bind(C, name="qmckl_set_ao_basis_ao_factor_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: ao_factor ! Elements of type double
            integer(c_int64_t), intent(in), value :: ao_num
        end function qmckl_set_ao_basis_ao_factor_device

        ! Getters (calling compute)

        integer(qmckl_exit_code_device) function qmckl_get_ao_basis_ao_num_device(context, ao_num) &
            bind(C, name="qmckl_get_ao_basis_ao_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: ao_num
        end function qmckl_get_ao_basis_ao_num_device


        integer(qmckl_exit_code_device) function qmckl_get_ao_basis_ao_vgl_device(context, ao_vgl, size_max) &
            bind(C, name="qmckl_get_ao_basis_ao_vgl_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: ao_vgl
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_ao_basis_ao_vgl_device



        integer(qmckl_exit_code_device) function qmckl_get_ao_basis_ao_value_device(context, ao_value, size_max) &
            bind(C, name="qmckl_get_ao_basis_ao_value_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value         :: ao_value
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_ao_basis_ao_value_device

        integer(qmckl_exit_code_device) function qmckl_get_ao_basis_ao_value_inplace_device(context, ao_value, size_max) &
            bind(C, name="qmckl_get_ao_basis_ao_value_inplace_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value         :: ao_value
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_ao_basis_ao_value_inplace_device


        integer(qmckl_exit_code_device) function qmckl_get_ao_basis_ao_vgl_inplace_device(context, ao_value, size_max) &
            bind(C, name="qmckl_get_ao_basis_ao_vgl_inplace_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value         :: ao_value
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_ao_basis_ao_vgl_inplace_device



    !!!!!!!!!!!
        ! MO
    !!!!!!!!!!!

        ! Basis setters

        integer(qmckl_exit_code_device) function qmckl_set_mo_basis_mo_num_device(context, num) &
            bind(C, name="qmckl_set_mo_basis_mo_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: num
        end function qmckl_set_mo_basis_mo_num_device

        integer(qmckl_exit_code_device) function qmckl_set_mo_basis_coefficient_device(context, mo_coefficient) &
            bind(C, name="qmckl_set_mo_basis_coefficient_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: mo_coefficient ! Elements of type double
        end function qmckl_set_mo_basis_coefficient_device

        ! Basis getters

        integer(qmckl_exit_code_device) function qmckl_get_mo_basis_mo_num_device(context, mo_num) &
            bind(C, name="qmckl_get_mo_basis_mo_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: mo_num
        end function qmckl_get_mo_basis_mo_num_device

        ! Getters (triggering computes)

        integer(qmckl_exit_code_device) function qmckl_get_mo_basis_mo_vgl_device(context, mo_vgl, size_max) &
            bind(C, name="qmckl_get_mo_basis_mo_vgl_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: mo_vgl
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_mo_basis_mo_vgl_device

        integer(qmckl_exit_code_device) function qmckl_get_mo_basis_mo_value_device(context, mo_value, size_max) &
            bind(C, name="qmckl_get_mo_basis_mo_value_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: mo_value
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_mo_basis_mo_value_device

        integer(qmckl_exit_code_device) function qmckl_get_mo_basis_mo_vgl_inplace_device(context, mo_vgl, size_max) &
            bind(C, name="qmckl_get_mo_basis_mo_vgl_inplace_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: mo_vgl
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_mo_basis_mo_vgl_inplace_device

        integer(qmckl_exit_code_device) function qmckl_get_mo_basis_mo_value_inplace_device(context, mo_value, size_max) &
            bind(C, name="qmckl_get_mo_basis_mo_value_inplace_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: mo_value
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_mo_basis_mo_value_inplace_device

        integer(qmckl_exit_code_device) function qmckl_mo_basis_select_mo_device(context, keep, size_max) &
            bind(C, name="qmckl_mo_basis_select_mo_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int32_t), intent(in), value :: keep
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_mo_basis_select_mo_device

    !!!!!!!!!!!
        ! MO
    !!!!!!!!!!!

        ! Setters

        integer(qmckl_exit_code_device) function qmckl_set_jastrow_rescale_factor_ee_device(context, kappa_ee) &
            bind(C, name="qmckl_set_jastrow_rescale_factor_ee_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: kappa_ee
        end function qmckl_set_jastrow_rescale_factor_ee_device

        integer(qmckl_exit_code_device) function qmckl_set_jastrow_rescale_factor_en_device(context, kappa_en, size_max) &
            bind(C, name="qmckl_set_jastrow_rescale_factor_en_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: kappa_en ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_set_jastrow_rescale_factor_en_device

        integer(qmckl_exit_code_device) function qmckl_set_jastrow_aord_num_device(context, aord_num) &
            bind(C, name="qmckl_set_jastrow_aord_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: aord_num
        end function qmckl_set_jastrow_aord_num_device

        integer(qmckl_exit_code_device) function qmckl_set_jastrow_bord_num_device(context, bord_num) &
            bind(C, name="qmckl_set_jastrow_bord_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: bord_num
        end function qmckl_set_jastrow_bord_num_device

        integer(qmckl_exit_code_device) function qmckl_set_jastrow_cord_num_device(context, cord_num) &
            bind(C, name="qmckl_set_jastrow_cord_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: cord_num
        end function qmckl_set_jastrow_cord_num_device

        integer(qmckl_exit_code_device) function qmckl_set_jastrow_type_nucl_num_device(context, type_nucl_num) &
            bind(C, name="qmckl_set_jastrow_type_nucl_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: type_nucl_num
        end function qmckl_set_jastrow_type_nucl_num_device

        integer(qmckl_exit_code_device) function qmckl_set_jastrow_type_nucl_vector_device(context, type_nucl_vector, nucl_num) &
            bind(C, name="qmckl_set_jastrow_type_nucl_vector_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integet(c_int64_t), intent(in), value :: type_nucl_vector ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: nucl_num
        end function qmckl_set_jastrow_type_nucl_vector_device

        integer(qmckl_exit_code_device) function qmckl_set_jastrow_a_vector_device(context, a_vector, size_max) &
            bind(C, name="qmckl_set_jastrow_a_vector_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: a_vector ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_set_jastrow_a_vector_device

        integer(qmckl_exit_code_device) function qmckl_set_jastrow_b_vector_device(context, b_vector, size_max) &
            bind(C, name="qmckl_set_jastrow_b_vector_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: b_vector ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_set_jastrow_b_vector_device

        integer(qmckl_exit_code_device) function qmckl_set_jastrow_c_vector_device(context, c_vector, size_max) &
            bind(C, name="qmckl_set_jastrow_c_vector_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: c_vector ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_set_jastrow_c_vector_device

        ! Getters (basic)

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_aord_num_device(context, aord_num) &
            bind(C, name="qmckl_get_jastrow_aord_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: aord_num ! Type int64_t
        end function qmckl_get_jastrow_aord_num_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_bord_num_device(context, bord_num) &
            bind(C, name="qmckl_get_jastrow_bord_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: bord_num ! Type int64_t
        end function qmckl_get_jastrow_bord_num_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_cord_num_device(context, cord_num) &
            bind(C, name="qmckl_get_jastrow_cord_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: cord_num ! Type int64_t
        end function qmckl_get_jastrow_cord_num_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_type_nucl_num_device(context, type_nucl_num) &
            bind(C, name="qmckl_get_jastrow_type_nucl_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: type_nucl_num ! Type int64_t
        end function qmckl_get_jastrow_type_nucl_num_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_type_nucl_vector_device(context, type_nucl_vector, size_max) &
            bind(C, name="qmckl_get_jastrow_type_nucl_vector_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: type_nucl_vector ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_type_nucl_vector_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_a_vector_device(context, a_vector, size_max) &
            bind(C, name="qmckl_get_jastrow_a_vector_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: a_vector ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_a_vector_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_b_vector_device(context, b_vector, size_max) &
            bind(C, name="qmckl_get_jastrow_b_vector_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: b_vector ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_b_vector_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_c_vector_device(context, c_vector, size_max) &
            bind(C, name="qmckl_get_jastrow_c_vector_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: c_vector ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_c_vector_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_rescale_factor_ee_device(context, rescale_factor_ee) &
            bind(C, name="qmckl_get_jastrow_rescale_factor_ee_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(out), value :: rescale_factor_ee ! Elements of type double
        end function qmckl_get_jastrow_rescale_factor_ee_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_rescale_factor_en_device(context, rescale_factor_en, size_max) &
            bind(C, name="qmckl_get_jastrow_rescale_factor_en_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: rescale_factor_en ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_rescale_factor_en_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_dim_c_vector_device(context, rescale_factor_en) &
            bind(C, name="qmckl_get_jastrow_dim_c_vector_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: rescale_factor_en ! Elements of type int64_t
        end function qmckl_get_jastrow_dim_c_vector_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_asymp_jasa_device(context, asymp_jasa, size_max) &
            bind(C, name="qmckl_get_jastrow_asymp_jasa_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: asymp_jasa ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_asymp_jasa_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_asymp_jasb_device(context, asymp_jasb, size_max) &
            bind(C, name="qmckl_get_jastrow_asymp_jasb_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: asymp_jasb ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_asymp_jasb_device

        ! Getters (for compute)

        ! Total Jastrow
        integer(qmckl_exit_code_device) function qmckl_get_jastrow_value_device(context, jastrow_value, size_max) &
            bind(C, name="qmckl_get_jastrow_value_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: jastrow_value ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_value_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_gl_device(context, jastrow_value, size_max) &
            bind(C, name="qmckl_get_jastrow_gl_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: jastrow_value ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_gl_device

        ! Electron/nucleus component
        integer(qmckl_exit_code_device) function qmckl_get_jastrow_factor_ee_device(context, factor_ee, size_max) &
            bind(C, name="qmckl_get_jastrow_factor_ee_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: factor_ee ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_factor_ee_device

        ! Electron/nucleus component
        integer(qmckl_exit_code_device) function qmckl_get_jastrow_factor_en_device(context, factor_en, size_max) &
            bind(C, name="qmckl_get_jastrow_factor_en_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: factor_en ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_factor_en_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_factor_en_deriv_e_device(context, factor_en_deriv_e, size_max) &
            bind(C, name="qmckl_get_jastrow_factor_en_deriv_e_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: factor_en_deriv_e ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_factor_en_deriv_e_device

        ! Electron/electron/nucleus component
        integer(qmckl_exit_code_device) function qmckl_get_jastrow_factor_een_device(context, factor_een, size_max) &
            bind(C, name="qmckl_get_jastrow_factor_een_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: factor_een ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_factor_een_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_factor_een_deriv_e_device &
        (context, factor_een_deriv_e, size_max) &
            bind(C, name="qmckl_get_jastrow_factor_een_deriv_e_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: factor_een_deriv_e ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_factor_een_deriv_e_device

        ! Distances
        integer(qmckl_exit_code_device) function qmckl_get_jastrow_ee_distance_rescaled_device(context, distance_rescaled) &
            bind(C, name="qmckl_get_jastrow_ee_distance_rescaled_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: distance_rescaled ! Elements of type double
        end function qmckl_get_jastrow_ee_distance_rescaled_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_ee_distance_rescaled_deriv_e_device &
            (context, distance_rescaled_deriv_e) &
            bind(C, name="qmckl_get_jastrow_ee_distance_rescaled_deriv_e_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: distance_rescaled_deriv_e ! Elements of type double
        end function qmckl_get_jastrow_ee_distance_rescaled_deriv_e_device

        integer(qmckl_exit_code_device) function qmckl_get_electron_en_distance_rescaled_device(context, distance_rescaled) &
            bind(C, name="qmckl_get_electron_en_distance_rescaled_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: distance_rescaled ! Elements of type double
        end function qmckl_get_electron_en_distance_rescaled_device

        integer(qmckl_exit_code_device) function qmckl_get_electron_en_distance_rescaled_deriv_e_device &
            (context, distance_rescaled_deriv_e) &
            bind(C, name="qmckl_get_electron_en_distance_rescaled_deriv_e_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: distance_rescaled_deriv_e ! Elements of type double
        end function qmckl_get_electron_en_distance_rescaled_deriv_e_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_een_rescaled_e_device(context, distance_rescaled, size_max) &
            bind(C, name="qmckl_get_jastrow_een_rescaled_e_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: distance_rescaled ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_een_rescaled_e_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_een_rescaled_e_deriv_e_device &
        (context, distance_rescaled, size_max) &
            bind(C, name="qmckl_get_jastrow_een_rescaled_e_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: distance_rescaled ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_een_rescaled_e_deriv_e_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_een_rescaled_n_device(context, distance_rescaled, size_max) &
            bind(C, name="qmckl_get_jastrow_een_rescaled_n_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: distance_rescaled ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_een_rescaled_n_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_een_rescaled_n_deriv_e_device &
        (context, distance_rescaled, size_max) &
            bind(C, name="qmckl_get_jastrow_een_rescaled_n_deriv_e_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: distance_rescaled ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_een_rescaled_n_deriv_e_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_tmp_c_device(context, tmp_c) &
            bind(C, name="qmckl_get_jastrow_tmp_c_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: tmp_c ! Elements of type double
        end function qmckl_get_jastrow_tmp_c_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_dtmp_c_device(context, dtmp_c) &
            bind(C, name="qmckl_get_jastrow_dtmp_c_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: dtmp_c ! Elements of type double
        end function qmckl_get_jastrow_dtmp_c_device

        integer(qmckl_exit_code_device) function qmckl_get_jastrow_factor_ee_deriv_e_device(context, factor_ee_deriv_e, size_max) &
            bind(C, name="qmckl_get_jastrow_factor_ee_deriv_e_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(in), value :: factor_ee_deriv_e ! Elements of type double
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_jastrow_factor_ee_deriv_e_device

    end interface
end module qmckl_gpu_f
