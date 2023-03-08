module qmckl_gpu
    integer, parameter :: qmckl_context_device = c_int64_t
    integer, parameter :: qmckl_exit_code = c_int32_t
    interface

    !!!!!!!!!!!
        ! CONTEXT
    !!!!!!!!!!!

        integer(qmckl_context_device) function qmckl_context_touch_device(context) &
            bind(C, name="qmckl_context_touch_device_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
        end function qmckl_context_touch_device

        integer(qmckl_context_device) function qmckl_context_create_device(device_id) &
            bind(C, name="qmckl_context_create_device_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(c_int32_t), value :: device_id
        end function qmckl_context_create_device

        integer(qmckl_exit_code) function qmckl_context_destroy_device(context) &
            bind(C, name="qmckl_context_destroy_device_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
        end function qmckl_context_destroy_device

    !!!!!!!!!!!
        ! MEMORY
    !!!!!!!!!!!

        integer(qmckl_exit_code) function qmckl_malloc_device(context, size) &
            bind(C, name="qmckl_malloc_device_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            integer(c_size_t), value :: size
        end function qmckl_malloc_device

        integer(qmckl_exit_code) function qmckl_free_device(context, ptr) &
            bind(C, name="qmckl_free_device_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            type(c_ptr), value :: ptr
        end function qmckl_free_device

        integer(qmckl_exit_code) function qmckl_memcpy_H2D(context, ptr) &
            bind(C, name="qmckl_memcpy_H2D_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            type(c_ptr), value :: dest
            type(c_ptr), value :: src
            integer(c_size_t), value :: size
        end function qmckl_memcpy_H2D

        integer(qmckl_exit_code) function qmckl_memcpy_D2H(context, ptr) &
            bind(C, name="qmckl_memcpy_D2H_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), value :: context
            type(c_ptr), value :: dest
            type(c_ptr), value :: src
            integer(c_size_t), value :: size
        end function qmckl_memcpy_D2H

        integer(qmckl_exit_code) function qmckl_memcpy_D2D(context, ptr) &
            bind(C, name="qmckl_memcpy_D2D_f")
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
            bind(C, name="qmckl_trexio_read_device_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: size_max
            character(c_char), intent(in)          :: file_name(size_max)
        end function qmckl_trexio_read_device

        integer(qmckl_exit_code) function qmckl_set_ao_basis_ao_num_device(context, num) &
            bind(C, name="qmckl_set_ao_basis_ao_num_device_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: num
        end function qmckl_set_ao_basis_ao_num_device

        integer(qmckl_exit_code) function qmckl_set_electron_coord_device(context, transp, walk_num, coord, size_max) &
            bind(C, name="qmckl_set_electron_coord_device_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character, intent(in), value :: transp
            integer(c_int64_t), intent(in), value :: walk_num
            double precision, intent(in)          :: coord(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function

        integer(qmckl_exit_code) function qmckl_set_point_device(context, transp, num, coord, size_max) &
            bind(C, name="qmckl_set_point_device_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character(c_char), intent(in), value :: transp
            integer(c_int64_t), intent(in), value :: num
            real(c_double), intent(in)          :: coord(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function

    !!!!!!!!!!!
        ! AO
    !!!!!!!!!!!

        integer(qmckl_exit_code) function qmckl_get_ao_basis_ao_vgl_device(context, ao_vgl, size_max) &
            bind(C, name="qmckl_get_ao_basis_ao_vgl_device_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            double precision, intent(out)         :: ao_vgl(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_ao_basis_ao_vgl_device

        integer(qmckl_exit_code) function qmckl_get_ao_basis_ao_value_device(context, ao_value, size_max) &
            bind(C, name="qmckl_get_ao_basis_ao_value_device_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            double precision, intent(out)         :: ao_value(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_ao_basis_ao_value_device


    !!!!!!!!!!!
        ! MO
    !!!!!!!!!!!

        integer(qmckl_exit_code) function qmckl_get_mo_basis_mo_vgl_device(context, mo_vgl, size_max) &
            bind(C, name="qmckl_get_mo_basis_mo_vgl_device_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            double precision, intent(out)         :: mo_vgl(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_mo_basis_mo_vgl_device

        integer(qmckl_exit_code) function qmckl_get_mo_basis_mo_value_device(context, mo_value, size_max) &
            bind(C, name="qmckl_get_mo_basis_mo_value_device_f")
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            double precision, intent(out)         :: mo_vgl(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_mo_basis_mo_value_device


    !!!!!!!!!!!
        ! ELECTRON
    !!!!!!!!!!!

    end interface
end module qmckl_gpu
