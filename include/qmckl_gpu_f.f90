module qmckl_gpu
    use, intrinsic :: iso_c_binding
    integer, parameter :: qmckl_context_device = c_int64_t
    integer, parameter :: qmckl_exit_code = c_int32_t
    interface

    !!!!!!!!!!!
        ! CONTEXT
    !!!!!!!!!!!

        integer(qmckl_context_device) function qmckl_context_touch_device(context) &
            bind(C, name="qmckl_context_touch_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
        end function qmckl_context_touch_device

        integer(qmckl_context_device) function qmckl_context_create_device(device_id) &
            bind(C, name="qmckl_context_create_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(c_int32_t), value :: device_id
        end function qmckl_context_create_device

        integer(qmckl_exit_code) function qmckl_context_destroy_device(context) &
            bind(C, name="qmckl_context_destroy_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
        end function qmckl_context_destroy_device

    !!!!!!!!!!!
        ! MEMORY
    !!!!!!!!!!!

        type(c_ptr) function qmckl_malloc_device(context, size) &
            bind(C, name="qmckl_malloc_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            integer(c_size_t), value :: size
        end function qmckl_malloc_device

        integer(qmckl_exit_code) function qmckl_free_device(context, ptr) &
            bind(C, name="qmckl_free_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            type(c_ptr), value :: ptr
        end function qmckl_free_device

        integer(qmckl_exit_code) function qmckl_memcpy_H2D(context, dest, src, size) &
            bind(C, name="qmckl_memcpy_H2D")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            type(c_ptr), value :: dest
            type(c_ptr), value :: src
            integer(c_size_t), value :: size
        end function qmckl_memcpy_H2D

        integer(qmckl_exit_code) function qmckl_memcpy_D2H(context, dest, src, size) &
            bind(C, name="qmckl_memcpy_D2H")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            type(c_ptr), value :: dest
            type(c_ptr), value :: src
            integer(c_size_t), value :: size
        end function qmckl_memcpy_D2H

        integer(qmckl_exit_code) function qmckl_memcpy_D2D(context, dest, src, size) &
            bind(C, name="qmckl_memcpy_D2D")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            type(c_ptr), value :: dest
            type(c_ptr), value :: src
            integer(c_size_t), value :: size
        end function qmckl_memcpy_D2D

    !!!!!!!!!!!
        ! TREXIO
    !!!!!!!!!!!

        integer(qmckl_exit_code) function qmckl_trexio_read_device &
            (context, file_name, size_max) &
            bind(C, name="qmckl_trexio_read_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: size_max
            character(c_char), intent(in)          :: file_name(size_max)
        end function qmckl_trexio_read_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_ao_num_device(context, num) &
            bind(C, name="qmckl_set_ao_basis_ao_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: num
        end function qmckl_set_ao_basis_ao_num_device

        integer(qmckl_exit_code) function qmckl_set_electron_coord_device(context, transp, walk_num, coord, size_max) &
            bind(C, name="qmckl_set_electron_coord_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character, intent(in), value :: transp
            integer(c_int64_t), intent(in), value :: walk_num
            real(c_double), intent(in)          :: coord(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_set_electron_coord_device

        integer(qmckl_exit_code) function qmckl_set_point_device(context, transp, num, coord, size_max) &
            bind(C, name="qmckl_set_point_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character(c_char), intent(in), value :: transp
            integer(c_int64_t), intent(in), value :: num
            real(c_double), intent(in)          :: coord(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_set_point_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_type_device(context, type) &
            bind(C, name="qmckl_set_ao_basis_shell_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character(c_signed_char), intent(in), value :: type
        end function qmckl_set_ao_basis_type_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_shell_num_device(context, shell_num) &
            bind(C, name="qmckl_set_ao_basis_shell_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: shell_num
        end function qmckl_set_ao_basis_shell_num_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_prim_num_device(context, prim_num) &
            bind(C, name="qmckl_set_ao_basis_prim_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: prim_num
        end function qmckl_set_ao_basis_prim_num_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_nucleus_index_device(context, nucleus_index, size_max) &
            bind(C, name="qmckl_set_ao_basis_nucleus_index_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            type(c_ptr), intent(in), value :: nucleus_index ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_set_ao_basis_nucleus_index_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_nucleus_shell_num_device(context, nucleus_shell_num, nucl_num) &
            bind(C, name="qmckl_set_ao_basis_nucleus_shell_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            type(c_ptr), intent(in), value :: nucleus_shell_num ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: nucl_num
        end function qmckl_set_ao_basis_nucleus_shell_num_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_shell_ang_mom_device(context, shell_ang_mom, shell_num) &
            bind(C, name="qmckl_set_ao_basis_shell_ang_mom_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            type(c_ptr), intent(in), value :: shell_ang_mom ! Elements of type int32_t
            integer(c_int64_t), intent(in), value :: shell_num
        end function qmckl_set_ao_basis_shell_ang_mom_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_shell_factor_device(context, shell_factor, shell_num) &
            bind(C, name="qmckl_set_ao_basis_shell_factor_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            type(c_ptr), intent(in), value :: shell_factor ! Elements of type double
            integer(c_int64_t), intent(in), value :: shell_num
        end function qmckl_set_ao_basis_shell_factor_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_shell_prim_num_device(context, shell_prim_num, shell_num) &
            bind(C, name="qmckl_set_ao_basis_shell_prim_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            type(c_ptr), intent(in), value :: shell_prim_num ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: shell_num
        end function qmckl_set_ao_basis_shell_prim_num_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_shell_prim_index_device(context, shell_prim_index, shell_num) &
            bind(C, name="qmckl_set_ao_basis_shell_prim_index_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            type(c_ptr), intent(in), value :: shell_prim_index ! Elements of type int64_t
            integer(c_int64_t), intent(in), value :: shell_num
        end function qmckl_set_ao_basis_shell_prim_index_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_exponent_device(context, exponent, prim_num) &
            bind(C, name="qmckl_set_ao_basis_exponent_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            type(c_ptr), intent(in), value :: exponent ! Elements of type double
            integer(c_int64_t), intent(in), value :: prim_num
        end function qmckl_set_ao_basis_exponent_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_coefficient_device(context, coefficient, prim_num) &
            bind(C, name="qmckl_set_ao_basis_coefficient_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            type(c_ptr), intent(in), value :: coefficient ! Elements of type double
            integer(c_int64_t), intent(in), value :: prim_num
        end function qmckl_set_ao_basis_coefficient_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_prim_factor_device(context, prim_factor, prim_num) &
            bind(C, name="qmckl_set_ao_basis_prim_factor_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            type(c_ptr), intent(in), value :: prim_factor ! Elements of type double
            integer(c_int64_t), intent(in), value :: prim_num
        end function qmckl_set_ao_basis_prim_factor_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_ao_factor_device(context, ao_factor, ao_num) &
            bind(C, name="qmckl_set_ao_basis_ao_factor_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            type(c_ptr), intent(in), value :: ao_factor ! Elements of type double
            integer(c_int64_t), intent(in), value :: ao_num
        end function qmckl_set_ao_basis_ao_factor_device

    !!!!!!!!!!!
        ! AO
    !!!!!!!!!!!

        integer(qmckl_exit_code) function qmckl_get_ao_basis_ao_vgl_device(context, ao_vgl, size_max) &
            bind(C, name="qmckl_get_ao_basis_ao_vgl_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(out)         :: ao_vgl(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_ao_basis_ao_vgl_device

        integer(qmckl_exit_code) function qmckl_get_ao_basis_ao_value_device(context, ao_value, size_max) &
            bind(C, name="qmckl_get_ao_basis_ao_value_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(out)         :: ao_value(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_ao_basis_ao_value_device

    !!!!!!!!!!!
        ! MO
    !!!!!!!!!!!

        integer(qmckl_exit_code) function qmckl_get_mo_basis_mo_vgl_device(context, mo_vgl, size_max) &
            bind(C, name="qmckl_get_mo_basis_mo_vgl_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(out)         :: mo_vgl(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_mo_basis_mo_vgl_device

        integer(qmckl_exit_code) function qmckl_get_mo_basis_mo_value_device(context, mo_value, size_max) &
            bind(C, name="qmckl_get_mo_basis_mo_value_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(out)         :: mo_value(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_mo_basis_mo_value_device

        integer(qmckl_exit_code) function qmckl_get_mo_basis_mo_vgl_inplace_device(context, mo_vgl, size_max) &
            bind(C, name="qmckl_get_mo_basis_mo_vgl_inplace_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(out)         :: mo_vgl(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_mo_basis_mo_vgl_inplace_device

        integer(qmckl_exit_code) function qmckl_get_mo_basis_mo_value_inplace_device(context, mo_value, size_max) &
            bind(C, name="qmckl_get_mo_basis_mo_value_inplace_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            real(c_double), intent(out)         :: mo_value(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_mo_basis_mo_value_inplace_device

        integer(qmckl_exit_code) function qmckl_get_mo_basis_mo_num_device(context, mo_num) &
            bind(C, name="qmckl_get_mo_basis_mo_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: mo_num
        end function qmckl_get_mo_basis_mo_num_device

        integer(qmckl_exit_code) function qmckl_mo_basis_select_mo_device(context, keep, size_max) &
            bind(C, name="qmckl_get_mo_basis_mo_num_device")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int32_t), intent(out)         :: keep(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_mo_basis_select_mo_device


    !!!!!!!!!!!
        ! ELECTRON
    !!!!!!!!!!!

    end interface
end module qmckl_gpu
