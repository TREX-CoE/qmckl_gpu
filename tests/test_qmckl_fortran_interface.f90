include '../include/qmckl_gpu_f.f90'

program qmckl_test_fortran_interface

    use, intrinsic :: iso_c_binding
    use :: qmckl_gpu

    !!!
    ! Interface to read chbrclf.h header
    !!!

    interface
        integer(kind=c_int64_t) function get_elec_up_num() &
            bind(C, name="get_elec_up_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_elec_up_num

        integer(kind=c_int64_t) function get_elec_dn_num() &
            bind(C, name="get_elec_dn_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_elec_dn_num

        integer(kind=c_int64_t) function get_nucl_num() &
            bind(C, name="get_nucl_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_nucl_num

        integer(kind=c_int64_t) function get_walk_num() &
            bind(C, name="get_walk_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_walk_num

        integer(kind=c_int64_t) function get_elec_num() &
            bind(C, name="get_elec_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_elec_num

        integer(kind=c_int64_t) function get_shell_num() &
            bind(C, name="get_shell_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_shell_num

        integer(kind=c_int64_t) function get_prim_num() &
            bind(C, name="get_prim_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_prim_num

        integer(kind=c_int64_t) function get_ao_num() &
            bind(C, name="get_ao_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_ao_num

        integer(kind=c_int64_t) function get_mo_num() &
            bind(C, name="get_mo_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_mo_num

        type(c_ptr) function get_elec_coord() result(elec_coord) &
            bind(C, name="get_elec_coord")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_elec_coord

        type(c_ptr) function get_nucl_charge() result(nucl_charge) &
            bind(C, name="get_nucl_charge")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_nucl_charge

        type(c_ptr) function get_nucl_coord() result(nucl_coord) &
            bind(C, name="get_nucl_coord")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_nucl_coord

        type(c_ptr) function get_nucleus_index() result(nucleus_index) &
            bind(C, name="get_nucleus_index")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_nucleus_index

        type(c_ptr) function get_nucleus_shell_num() result(nucleus_shell_num) &
            bind(C, name="get_nucleus_shell_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_nucleus_shell_num

        type(c_ptr) function get_shell_ang_mom() result(shell_ang_mom) &
            bind(C, name="get_shell_ang_mom")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_shell_ang_mom


        type(c_ptr) function get_shell_prim_num() result(shell_prim_num) &
            bind(C, name="get_shell_prim_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_shell_prim_num

        type(c_ptr) function get_shell_prim_index() result(shell_prim_index) &
            bind(C, name="get_shell_prim_index")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_shell_prim_index

        type(c_ptr) function get_shell_factor() result(shell_factor) &
            bind(C, name="get_shell_factor")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_shell_factor

        type(c_ptr) function get_exponent() result(exponent) &
            bind(C, name="get_exponent")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_exponent

        type(c_ptr) function get_coefficient() result(coefficient) &
            bind(C, name="get_coefficient")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_coefficient

        type(c_ptr) function get_prim_factor() result(prim_factor) &
            bind(C, name="get_prim_factor")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_prim_factor

        type(c_ptr) function get_ao_factor() result(ao_factor) &
            bind(C, name="get_ao_factor")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_ao_factor

    end interface

    !!!
    ! Declarations
    !!!

    integer(qmckl_context_device) context
    integer(qmckl_exit_code) rc

    character(c_signed_char) :: typ

    integer(c_int64_t) elec_up_num
    integer(c_int64_t) elec_dn_num
    integer(c_int64_t) nucl_num
    integer(c_int64_t) walk_num
    integer(c_int64_t) point_num
    integer(c_int64_t) shell_num
    integer(c_int64_t) ao_num
    integer(c_int64_t) mo_num
    integer(c_int64_t) prim_num
    integer(c_int64_t) elec_num

    type(c_ptr) elec_coord
    type(c_ptr) nucl_charge
    type(c_ptr) nucl_coord
    type(c_ptr) elec_coord_d
    type(c_ptr) nucl_charge_d
    type(c_ptr) nucl_coord_d

    type(c_ptr) nucleus_index
    type(c_ptr) nucleus_shell_num
    type(c_ptr) shell_ang_mom
    type(c_ptr) shell_prim_num
    type(c_ptr) shell_prim_index
    type(c_ptr) shell_factor
    type(c_ptr) exponent
    type(c_ptr) coefficient
    type(c_ptr) prim_factor
    type(c_ptr) ao_factor
    type(c_ptr) nucleus_index_d
    type(c_ptr) nucleus_shell_num_d
    type(c_ptr) shell_ang_mom_d
    type(c_ptr) shell_prim_num_d
    type(c_ptr) shell_prim_index_d
    type(c_ptr) shell_factor_d
    type(c_ptr) exponent_d
    type(c_ptr) coefficient_d
    type(c_ptr) prim_factor_d
    type(c_ptr) ao_factor_d

    type(c_ptr) ao_vgl
    type(c_ptr) mo_vgl
    type(c_ptr) ao_vgl_d
    type(c_ptr) mo_vgl_d

    context = qmckl_context_create_device(0)


    !!!
    ! Read CPU Fortran arrays from the .h file
    !!!

    walk_num = get_walk_num();
    elec_num = get_elec_num();
    shell_num = get_shell_num();
    ao_num = get_ao_num();
    prim_num = get_prim_num();

    elec_up_num = get_elec_up_num();
    elec_dn_num = get_elec_dn_num();
    elec_coord = get_elec_coord();
    nucl_num = get_nucl_num();
    nucl_charge = get_nucl_charge();
    nucl_coord = get_nucl_coord();
    point_num = elec_num * walk_num;
    mo_num = get_mo_num();

    nucleus_index = get_nucleus_index();
    nucleus_shell_num = get_nucleus_shell_num();
    shell_ang_mom = get_shell_ang_mom();
    shell_prim_num = get_shell_prim_num();
    shell_prim_index = get_shell_prim_index();
    shell_factor = get_shell_factor();
    exponent = get_exponent();
    coefficient = get_coefficient();
    prim_factor = get_prim_factor();
    ao_factor = get_ao_factor();


    !!!
    ! Allocate GPU arrays and copy the CPU arrays onto them
    !!!
    elec_coord_d = qmckl_malloc_device(context, point_num * 3 * c_sizeof(c_double) * 2);
    rc = qmckl_memcpy_H2D(context, elec_coord_d, elec_coord, point_num * 3 * c_sizeof(c_double) * 2);

    nucl_coord_d = qmckl_malloc_device(context, 3 * nucl_num * c_sizeof(c_double) * 2);
    rc = qmckl_memcpy_H2D(context, nucl_coord_d, nucl_coord, 3 * nucl_num * c_sizeof(c_double) * 2);

    nucl_charge_d = qmckl_malloc_device(context, nucl_num * c_sizeof(c_double) * 2);
    rc = qmckl_memcpy_H2D(context, nucl_charge_d, nucl_charge, nucl_num * c_sizeof(c_double) * 2);


    nucleus_index_d = qmckl_malloc_device(context, nucl_num * c_sizeof(c_int64_t) * 2);
    rc = qmckl_memcpy_H2D(context, nucleus_index_d, nucleus_index, nucl_num * c_sizeof(c_int64_t) * 2);

    nucleus_shell_num_d = qmckl_malloc_device(context, nucl_num * c_sizeof(c_int64_t) * 2);
    rc = qmckl_memcpy_H2D(context, nucleus_shell_num_d, nucleus_shell_num, nucl_num * c_sizeof(c_int64_t) * 2);

    shell_ang_mom_d = qmckl_malloc_device(context, shell_num * c_sizeof(c_int32_t));
    rc = qmckl_memcpy_H2D(context, shell_ang_mom_d, shell_ang_mom, shell_num * c_sizeof(c_int32_t));

    shell_prim_num_d = qmckl_malloc_device(context, shell_num * c_sizeof(c_int64_t) * 2);
    rc = qmckl_memcpy_H2D(context, shell_prim_num_d, shell_prim_num, shell_num * c_sizeof(c_int64_t) * 2);

    shell_prim_index_d = qmckl_malloc_device(context, shell_num * c_sizeof(c_int64_t) * 2);
    rc = qmckl_memcpy_H2D(context, shell_prim_index_d, shell_prim_index, shell_num * c_sizeof(c_int64_t) * 2);

    shell_factor_d = qmckl_malloc_device(context, shell_num * c_sizeof(c_double) * 2);
    rc = qmckl_memcpy_H2D(context, shell_factor_d, shell_factor, shell_num * c_sizeof(c_double) * 2);

    exponent_d = qmckl_malloc_device(context, prim_num * c_sizeof(c_double) * 2);
    rc = qmckl_memcpy_H2D(context, exponent_d, exponent, prim_num * c_sizeof(c_double) * 2);

    coefficient_d = qmckl_malloc_device(context, prim_num * c_sizeof(c_double) * 2);
    rc = qmckl_memcpy_H2D(context, coefficient_d, coefficient, prim_num * c_sizeof(c_double) * 2);

    prim_factor_d = qmckl_malloc_device(context, prim_num * c_sizeof(c_double) * 2);
    rc = qmckl_memcpy_H2D(context, prim_factor_d, prim_factor, prim_num * c_sizeof(c_double) * 2);

    ao_factor_d = qmckl_malloc_device(context, ao_num * c_sizeof(c_double) * 2);
    rc = qmckl_memcpy_H2D(context, ao_factor_d, ao_factor, ao_num * c_sizeof(c_double) * 2);

    !!!
    ! Set context values
    !!!

    ! electron / nucleus
    rc = qmckl_set_electron_num_device(context, elec_up_num, elec_dn_num);
    rc = qmckl_set_point_device(context, 'N', point_num, elec_coord_d, point_num * 3);
    rc = qmckl_set_nucleus_num_device(context, nucl_num);
    rc = qmckl_set_nucleus_coord_device(context, 'T', nucl_coord_d, nucl_num * 3);
    rc = qmckl_set_nucleus_charge_device(context, nucl_charge_d, nucl_num);

    ! ao_basis
    rc = qmckl_set_ao_basis_type_device(context, 'G');
    rc = qmckl_set_ao_basis_shell_num_device(context, shell_num);
    rc = qmckl_set_ao_basis_prim_num_device(context, prim_num);
    rc = qmckl_set_ao_basis_nucleus_index_device(context, nucleus_index_d, nucl_num);
    rc = qmckl_set_ao_basis_nucleus_shell_num_device(context, nucleus_shell_num_d, nucl_num);
    rc = qmckl_set_ao_basis_shell_ang_mom_device(context, shell_ang_mom_d, shell_num);
    rc = qmckl_set_ao_basis_shell_factor_device(context, shell_factor_d, shell_num);
    rc = qmckl_set_ao_basis_shell_prim_num_device(context, shell_prim_num_d, shell_num);
    rc = qmckl_set_ao_basis_shell_prim_index_device(context, shell_prim_index_d, shell_num);
    rc = qmckl_set_ao_basis_exponent_device(context, exponent_d, prim_num);
    rc = qmckl_set_ao_basis_coefficient_device(context, coefficient_d, prim_num);
    rc = qmckl_set_ao_basis_prim_factor_device(context, prim_factor_d, prim_num);

    rc = qmckl_set_ao_basis_ao_num_device(context, ao_num);
    rc = qmckl_set_ao_basis_ao_factor_device(context, ao_factor_d, ao_num);

    ! mo_basis
    rc = qmckl_set_mo_basis_mo_num_device(context, mo_num);
    rc = qmckl_set_mo_basis_coefficient_device(context, mo_coefficient_d);


    !!!
    ! AO computations
    !!!

    ao_vgl_d = qmckl_malloc_device(context, 5 * point_num * ao_num * c_sizeof(c_double) * 2);
    rc = qmckl_get_ao_basis_ao_vgl_device(context, ao_vgl_d, 5 * point_num * ao_num);


    !!!
    ! MO computations
    !!!

    mo_vgl_d = qmckl_malloc_device(context, 5 * point_num * mo_num * c_sizeof(c_double) * 2);
    rc = qmckl_get_mo_basis_mo_vgl_device(context, ao_vgl_d, 5 * point_num * mo_num);

end program
