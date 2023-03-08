module qmckl_gpu
    use, intrinsic :: iso_c_binding
    use qmckl
    integer, parameter :: qmckl_context_device = c_int64_t

    interface
        integer(qmckl_context_device) function qmckl_context_create_omp_device(device_id) bind(C)
            use, intrinsic :: iso_c_binding
            import
            integer(c_int32_t), value :: device_id
        end function qmckl_context_create_omp_device
    end interface

    interface
        integer(c_int32_t) function qmckl_trexio_read_omp_device &
            (context, file_name, size_max) &
            bind(C)
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: size_max
            character(c_char), intent(in)          :: file_name(size_max)

        end function qmckl_trexio_read_omp_device
    end interface

    interface
        integer(c_int32_t) function qmckl_set_ao_basis_ao_num_omp_device(context, &
                                                                         num) bind(C)
            use, intrinsic :: iso_c_binding
            import
            implicit none
            integer(qmckl_context_device), intent(in), value :: context
            integer(c_int64_t), intent(in), value :: num
        end function qmckl_set_ao_basis_ao_num_omp_device
    end interface

    interface
        integer(c_int32_t) function qmckl_get_ao_basis_ao_vgl_omp_device(context, &
                                                                         ao_vgl, size_max) bind(C)
            use, intrinsic :: iso_c_binding
            import
            implicit none
            integer(qmckl_context_device), intent(in), value :: context
            double precision, intent(out)         :: ao_vgl(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function qmckl_get_ao_basis_ao_vgl_omp_device
    end interface

    interface
        integer(c_int32_t) function qmckl_set_electron_coord_omp_device(context, transp, walk_num, coord, size_max) bind(C)
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character, intent(in), value :: transp
            integer(c_int64_t), intent(in), value :: walk_num
            double precision, intent(in)          :: coord(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function
    end interface

    interface
        integer(c_int32_t) function qmckl_set_point_omp_device(context, &
                                                               transp, num, coord, size_max) bind(C)
            use, intrinsic :: iso_c_binding
            import
            implicit none

            integer(qmckl_context_device), intent(in), value :: context
            character(c_char), intent(in), value :: transp
            integer(c_int64_t), intent(in), value :: num
            real(c_double), intent(in)          :: coord(*)
            integer(c_int64_t), intent(in), value :: size_max
        end function
    end interface

end module qmckl_gpu

