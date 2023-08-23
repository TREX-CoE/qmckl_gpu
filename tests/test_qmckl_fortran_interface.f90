include '../include/qmckl_gpu_f.f90'

program qmckl_test_fortran_interface
    use, intrinsic :: iso_c_binding
    use :: qmckl_gpu_f
    implicit none

    !!!
    ! Interface to read chbrclf.h and n2.H header
    !!!

    interface

       !
       ! chbrclf
       !

        integer(kind=c_int64_t) function get_chbrclf_elec_up_num() &
            bind(C, name="get_chbrclf_elec_up_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_elec_up_num

        integer(kind=c_int64_t) function get_chbrclf_elec_dn_num() &
            bind(C, name="get_chbrclf_elec_dn_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_elec_dn_num

        integer(kind=c_int64_t) function get_chbrclf_nucl_num() &
            bind(C, name="get_chbrclf_nucl_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_nucl_num

        integer(kind=c_int64_t) function get_chbrclf_walk_num() &
            bind(C, name="get_chbrclf_walk_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_walk_num

        integer(kind=c_int64_t) function get_chbrclf_elec_num() &
            bind(C, name="get_chbrclf_elec_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_elec_num

        integer(kind=c_int64_t) function get_chbrclf_shell_num() &
            bind(C, name="get_chbrclf_shell_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_shell_num

        integer(kind=c_int64_t) function get_chbrclf_prim_num() &
            bind(C, name="get_chbrclf_prim_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_prim_num

        integer(kind=c_int64_t) function get_chbrclf_ao_num() &
            bind(C, name="get_chbrclf_ao_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_ao_num

        integer(kind=c_int64_t) function get_chbrclf_mo_num() &
            bind(C, name="get_chbrclf_mo_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_mo_num

        type(c_ptr) function get_chbrclf_elec_coord() result(elec_coord) &
            bind(C, name="get_chbrclf_elec_coord")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_elec_coord

        type(c_ptr) function get_chbrclf_nucl_charge() result(nucl_charge) &
            bind(C, name="get_chbrclf_nucl_charge")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_nucl_charge

        type(c_ptr) function get_chbrclf_nucl_coord() result(nucl_coord) &
            bind(C, name="get_chbrclf_nucl_coord")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_nucl_coord

        type(c_ptr) function get_chbrclf_nucleus_index() result(nucleus_index) &
            bind(C, name="get_chbrclf_nucleus_index")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_nucleus_index

        type(c_ptr) function get_chbrclf_nucleus_shell_num() result(nucleus_shell_num) &
            bind(C, name="get_chbrclf_nucleus_shell_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_nucleus_shell_num

        type(c_ptr) function get_chbrclf_shell_ang_mom() result(shell_ang_mom) &
            bind(C, name="get_chbrclf_shell_ang_mom")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_shell_ang_mom

        type(c_ptr) function get_chbrclf_shell_prim_num() result(shell_prim_num) &
            bind(C, name="get_chbrclf_shell_prim_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_shell_prim_num

        type(c_ptr) function get_chbrclf_shell_prim_index() result(shell_prim_index) &
            bind(C, name="get_chbrclf_shell_prim_index")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_shell_prim_index

        type(c_ptr) function get_chbrclf_shell_factor() result(shell_factor) &
            bind(C, name="get_chbrclf_shell_factor")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_shell_factor

        type(c_ptr) function get_chbrclf_exponent() result(exponent) &
            bind(C, name="get_chbrclf_exponent")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_exponent

        type(c_ptr) function get_chbrclf_coefficient() result(coefficient) &
            bind(C, name="get_chbrclf_coefficient")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_coefficient

        type(c_ptr) function get_chbrclf_prim_factor() result(prim_factor) &
            bind(C, name="get_chbrclf_prim_factor")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_prim_factor

        type(c_ptr) function get_chbrclf_ao_factor() result(ao_factor) &
            bind(C, name="get_chbrclf_ao_factor")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_ao_factor

        type(c_ptr) function get_chbrclf_mo_coef() result(mo_coef) &
            bind(C, name="get_chbrclf_mo_coef")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_chbrclf_mo_coef


       !
       ! n2
       !

        integer(kind=c_int64_t) function get_n2_nucl_num() &
            bind(C, name="get_n2_nucl_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_nucl_num

        type(c_ptr) function get_n2_charge() result(nucl_charge) &
            bind(C, name="get_n2_charge")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_charge

        integer(kind=c_int64_t) function get_n2_elec_up_num() &
            bind(C, name="get_n2_elec_up_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_elec_up_num

        integer(kind=c_int64_t) function get_n2_elec_dn_num() &
            bind(C, name="get_n2_elec_dn_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_elec_dn_num

        integer(kind=c_int64_t) function get_n2_walk_num() &
            bind(C, name="get_n2_walk_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_walk_num

        integer(kind=c_int64_t) function get_n2_elec_num() &
            bind(C, name="get_n2_elec_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_elec_num

        type(c_ptr) function get_n2_elec_coord() result(elec_coord) &
            bind(C, name="get_n2_elec_coord")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_elec_coord

        integer(kind=c_int64_t) function get_n2_type_nucl_num() &
            bind(C, name="get_n2_type_nucl_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_type_nucl_num

        type(c_ptr) function get_n2_aord_vector() result(aord_vector) &
            bind(C, name="get_n2_aord_vector")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_aord_vector

        type(c_ptr) function get_n2_bord_vector() result(bord_vector) &
            bind(C, name="get_n2_aord_vector")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_bord_vector

        type(c_ptr) function get_n2_cord_vector() result(cord_vector) &
            bind(C, name="get_n2_cord_vector")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_cord_vector

        type(c_ptr) function get_n2_cord_vector_full() result(cord_vector_full) &
            bind(C, name="get_n2_cord_vector_full")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_cord_vector_full

        type(c_ptr) function get_n2_lkpm_of_cindex() result(lkpm_of_cindex) &
            bind(C, name="get_n2_lkpm_of_cindex")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_lkpm_of_cindex

        type(c_ptr) function get_n2_rescale_factor_en() result(rescale_factor_en) &
            bind(C, name="get_n2_rescale_factor_en")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_rescale_factor_en
        real(kind=c_double) function get_n2_rescale_factor_ee() result(rescale_factor_ee) &
            bind(C, name="get_n2_rescale_factor_ee")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_rescale_factor_ee

        integer(kind=c_int64_t) function get_n2_aord_num() &
            bind(C, name="get_n2_aord_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_aord_num
        integer(kind=c_int64_t) function get_n2_bord_num() &
            bind(C, name="get_n2_bord_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_bord_num
        integer(kind=c_int64_t) function get_n2_cord_num() &
            bind(C, name="get_n2_cord_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_cord_num
        integer(kind=c_int64_t) function get_n2_dim_cord_vec() &
            bind(C, name="get_n2_dim_cord_vec")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_dim_cord_vec

        type(c_ptr) function get_n2_type_nucl_vector() result(type_nucl_vector) &
            bind(C, name="get_n2_type_nucl_vector")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_type_nucl_vector

        type(c_ptr) function get_n2_nucl_coord() result(nucl_coord) &
            bind(C, name="get_n2_nucl_coord")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_n2_nucl_coord

    end interface

    !!!
    ! Declarations
    !!!

    integer(qmckl_context_device) context
    integer(qmckl_exit_code_device) rc

    integer i
    integer j
    integer k

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
    type(c_ptr) mo_coef
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
    type(c_ptr) mo_coef_d

    type(c_ptr) ao_vgl_h
    type(c_ptr) mo_vgl_h
    type(c_ptr) ao_vgl_d
    type(c_ptr) mo_vgl_d
    real(8), pointer :: ao_vgl(:,:,:)
    real(8), pointer :: mo_vgl(:,:,:)

    real(8) ao_ref
    real(8) mo_ref

    ! n2
    type(c_ptr) factor_ee_d
    type(c_ptr) factor_en_d
    type(c_ptr) factor_een_d
    type(c_ptr) total_j_d

    type(c_ptr) factor_ee_h
    type(c_ptr) factor_en_h
    type(c_ptr) factor_een_h
    type(c_ptr) total_j_h

    real(8), pointer :: factor_ee(:)
    real(8), pointer :: factor_en(:)
    real(8), pointer :: factor_een(:)
    real(8), pointer :: total_j(:)

    integer(c_int64_t) type_nucl_num
    type(c_ptr) type_nucl_vector
    type(c_ptr) type_nucl_vector_d

    type(c_ptr) rescale_factor_en
    type(c_ptr) rescale_factor_en_d
    real(c_double) rescale_factor_ee


    integer(c_int64_t) aord_num
    integer(c_int64_t) bord_num
    integer(c_int64_t) cord_num
    integer(c_int64_t) dim_cord_vec
    type(c_ptr) a_vector
    type(c_ptr) b_vector
    type(c_ptr) c_vector
    type(c_ptr) a_vector_d
    type(c_ptr) b_vector_d
    type(c_ptr) c_vector_d

    context = qmckl_context_create_device(0)

    !!!
    ! Read CPU Fortran arrays from the chbrclf.h file
    !!!

    walk_num = get_chbrclf_walk_num();
    elec_num = get_chbrclf_elec_num();
    shell_num = get_chbrclf_shell_num();
    ao_num = get_chbrclf_ao_num();
    prim_num = get_chbrclf_prim_num();
    elec_up_num = get_chbrclf_elec_up_num();
    elec_dn_num = get_chbrclf_elec_dn_num();
    elec_coord = get_chbrclf_elec_coord();
    nucl_num = get_chbrclf_nucl_num();
    nucl_charge = get_chbrclf_nucl_charge();
    nucl_coord = get_chbrclf_nucl_coord();
    point_num = elec_num;
    mo_num = get_chbrclf_mo_num();
    nucleus_index = get_chbrclf_nucleus_index();
    nucleus_shell_num = get_chbrclf_nucleus_shell_num();
    shell_ang_mom = get_chbrclf_shell_ang_mom();
    shell_prim_num = get_chbrclf_shell_prim_num();
    shell_prim_index = get_chbrclf_shell_prim_index();
    shell_factor = get_chbrclf_shell_factor();
    exponent = get_chbrclf_exponent();
    coefficient = get_chbrclf_coefficient();
    prim_factor = get_chbrclf_prim_factor();
    ao_factor = get_chbrclf_ao_factor();
    mo_coef = get_chbrclf_mo_coef();

    !!!
    ! Allocate GPU arrays and copy the CPU arrays onto them
    !!!
    elec_coord_d = qmckl_malloc_device(context, point_num*3*c_sizeof(c_double)*2); 
    rc = qmckl_memcpy_H2D(context, elec_coord_d, elec_coord, point_num*3*c_sizeof(c_double)*2); 
    nucl_coord_d = qmckl_malloc_device(context, 3*nucl_num*c_sizeof(c_double)*2); 
    rc = qmckl_memcpy_H2D(context, nucl_coord_d, nucl_coord, 3*nucl_num*c_sizeof(c_double)*2); 
    nucl_charge_d = qmckl_malloc_device(context, nucl_num*c_sizeof(c_double)*2); 
    rc = qmckl_memcpy_H2D(context, nucl_charge_d, nucl_charge, nucl_num*c_sizeof(c_double)*2); 
    nucleus_index_d = qmckl_malloc_device(context, nucl_num*c_sizeof(c_int64_t)*2); 
    rc = qmckl_memcpy_H2D(context, nucleus_index_d, nucleus_index, nucl_num*c_sizeof(c_int64_t)*2); 
    nucleus_shell_num_d = qmckl_malloc_device(context, nucl_num*c_sizeof(c_int64_t)*2); 
    rc = qmckl_memcpy_H2D(context, nucleus_shell_num_d, nucleus_shell_num, nucl_num*c_sizeof(c_int64_t)*2); 
    shell_ang_mom_d = qmckl_malloc_device(context, shell_num*c_sizeof(c_int32_t)); 
    rc = qmckl_memcpy_H2D(context, shell_ang_mom_d, shell_ang_mom, shell_num*c_sizeof(c_int32_t)); 
    shell_prim_num_d = qmckl_malloc_device(context, shell_num*c_sizeof(c_int64_t)*2); 
    rc = qmckl_memcpy_H2D(context, shell_prim_num_d, shell_prim_num, shell_num*c_sizeof(c_int64_t)*2); 
    shell_prim_index_d = qmckl_malloc_device(context, shell_num*c_sizeof(c_int64_t)*2); 
    rc = qmckl_memcpy_H2D(context, shell_prim_index_d, shell_prim_index, shell_num*c_sizeof(c_int64_t)*2); 
    shell_factor_d = qmckl_malloc_device(context, shell_num*c_sizeof(c_double)*2); 
    rc = qmckl_memcpy_H2D(context, shell_factor_d, shell_factor, shell_num*c_sizeof(c_double)*2); 
    exponent_d = qmckl_malloc_device(context, prim_num*c_sizeof(c_double)*2); 
    rc = qmckl_memcpy_H2D(context, exponent_d, exponent, prim_num*c_sizeof(c_double)*2); 
    coefficient_d = qmckl_malloc_device(context, prim_num*c_sizeof(c_double)*2); 
    rc = qmckl_memcpy_H2D(context, coefficient_d, coefficient, prim_num*c_sizeof(c_double)*2); 
    prim_factor_d = qmckl_malloc_device(context, prim_num*c_sizeof(c_double)*2); 
    rc = qmckl_memcpy_H2D(context, prim_factor_d, prim_factor, prim_num*c_sizeof(c_double)*2); 
    ao_factor_d = qmckl_malloc_device(context, ao_num*c_sizeof(c_double)*2); 
    rc = qmckl_memcpy_H2D(context, ao_factor_d, ao_factor, ao_num*c_sizeof(c_double)*2);
    mo_coef_d = qmckl_malloc_device(context, mo_num*ao_num*c_sizeof(c_double)*2);
    rc = qmckl_memcpy_H2D(context, mo_coef_d, mo_coef, mo_num*ao_num*c_sizeof(c_double)*2);


    !!!
    ! Set context values
    !!!

    ! electron / nucleus
    rc = qmckl_set_electron_num_device(context, elec_up_num, elec_dn_num);
    rc = qmckl_set_point_device(context, 'N', point_num, elec_coord_d, point_num*3);
    rc = qmckl_set_nucleus_num_device(context, nucl_num);
    rc = qmckl_set_nucleus_coord_device(context, 'T', nucl_coord_d, nucl_num*3);
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
    rc = qmckl_set_mo_basis_coefficient_device(context, mo_coef_d);


    !!!
    ! AO computations
    !!!

    ao_vgl_d = qmckl_malloc_device(context, 5*point_num*ao_num*c_sizeof(c_double)*2);
    rc = qmckl_get_ao_basis_ao_vgl_device(context, ao_vgl_d, 5*point_num*ao_num);

    ! Copy values back to CPU in Fortran ptr
    ao_vgl_h = qmckl_malloc_host(context, 5*point_num*ao_num*c_sizeof(c_double)*2);
    rc = qmckl_memcpy_D2H(context, ao_vgl_h, ao_vgl_d, 5*point_num*ao_num*c_sizeof(c_double)*2);

    call c_f_pointer(ao_vgl_h, ao_vgl, [int(point_num, kind(4)), 5, int(ao_num, kind(4))]);

    ! Compare to reference
    open (unit=1, file='tests/ao_reference.txt', status='old', action='read')
    do k = 1, ao_num
        do j = 1, 5
            do i = 1, point_num
                read(1, *), ao_ref
                if (abs(ao_vgl(i,j,k) - ao_ref) > 1e-12) then
                   print *, "Error at (i,j,k)=", i, j, k
                   print *, "ao_vgl =", ao_vgl(i,j,k)
                   print *, "ao_ref =", ao_ref
                   call exit(1)
                end if
            end do
        end do
    end do

    !!!
    ! MO computations
    !!!

    mo_vgl_d = qmckl_malloc_device(context, 5*point_num*mo_num*c_sizeof(c_double)*2);
    rc = qmckl_get_mo_basis_mo_vgl_device(context, mo_vgl_d, 5*point_num*mo_num);

    ! Copy values back to CPU in Fortran ptr
    mo_vgl_h = qmckl_malloc_host(context, 5*point_num*mo_num*c_sizeof(c_double)*2);
    rc = qmckl_memcpy_D2H(context, mo_vgl_h, mo_vgl_d, 5*point_num*mo_num*c_sizeof(c_double)*2);

    call c_f_pointer(mo_vgl_h, mo_vgl, [int(point_num, kind(4)), 5, int(mo_num, kind(4))]);

    ! Compare to reference
    open (unit=2, file='tests/mo_reference.txt', status='old', action='read')
    do k = 1, mo_num
        do j = 1, 5
            do i = 1, point_num
                read(2, *), mo_ref
                if (abs(mo_vgl(i,j,k) - mo_ref) > 1e-12) then
                   print *, "Error at (i,j,k)=", i, j, k
                   print *, "mo_vgl =", mo_vgl(i,j,k)
                   print *, "mo_ref =", mo_ref
                   call exit(1)
                end if
            end do
        end do
    end do


    !!!
    ! Jastrow computations
    !!!

    ! Destroy old chbrclf context (this will free all GPU arrays)
    ! rc = qmckl_context_destroy_device(context);

    ! Reinitialize context with the n2 dataset

    ! Get CPU values
    walk_num = get_n2_walk_num();
    elec_num = get_n2_elec_num();
    elec_up_num = get_n2_elec_up_num();
    elec_dn_num = get_n2_elec_dn_num();
    nucl_num = get_n2_nucl_num();
    type_nucl_num = get_n2_type_nucl_num();
    rescale_factor_ee = 0.6;

    rescale_factor_en = get_n2_rescale_factor_en();
    rescale_factor_ee = get_n2_rescale_factor_ee();

    aord_num = get_n2_aord_num();
    bord_num = get_n2_bord_num();
    cord_num = get_n2_cord_num();
    dim_cord_vec = get_n2_dim_cord_vec();
    a_vector = get_n2_aord_vector();
    b_vector = get_n2_bord_vector();
    c_vector = get_n2_cord_vector();

    elec_coord = get_n2_elec_coord();
    nucl_charge = get_n2_charge();
    nucl_coord = get_n2_nucl_coord();
    type_nucl_vector = get_n2_type_nucl_vector();

    ! Move CPU arrays to GPU
    elec_coord_d = qmckl_malloc_device(context, walk_num*3*elec_num*c_sizeof(c_double)*2);
    rc = qmckl_memcpy_H2D(context, elec_coord_d, elec_coord, walk_num*3*elec_num*c_sizeof(c_double)*2);
    nucl_coord_d = qmckl_malloc_device(context, 3*nucl_num*c_sizeof(c_double)*2);
    rc = qmckl_memcpy_H2D(context, nucl_coord_d, nucl_coord, 3*nucl_num*c_sizeof(c_double)*2);
    nucl_charge_d = qmckl_malloc_device(context, nucl_num*c_sizeof(c_double)*2);
    rc = qmckl_memcpy_H2D(context, nucl_charge_d, nucl_charge, nucl_num*c_sizeof(c_double)*2);
    rescale_factor_en_d = qmckl_malloc_device(context, 2*c_sizeof(c_double)*2);
    rc = qmckl_memcpy_H2D(context, rescale_factor_en_d, rescale_factor_en, 2*c_sizeof(c_double)*2);

    type_nucl_vector_d = qmckl_malloc_device(context, nucl_num*c_sizeof(c_double)*2);
    rc = qmckl_memcpy_H2D(context, type_nucl_vector_d, type_nucl_vector, 2*c_sizeof(c_double)*2);
    a_vector_d = qmckl_malloc_device(context, (aord_num+1)*type_nucl_num*c_sizeof(c_double)*2);
    rc = qmckl_memcpy_H2D(context, a_vector_d, a_vector, 2*c_sizeof(c_double)*2);
    b_vector_d = qmckl_malloc_device(context, (bord_num+1)*c_sizeof(c_double)*2);
    rc = qmckl_memcpy_H2D(context, b_vector_d, b_vector, 2*c_sizeof(c_double)*2);
    c_vector_d = qmckl_malloc_device(context, dim_cord_vec*type_nucl_num*c_sizeof(c_double)*2);
    rc = qmckl_memcpy_H2D(context, c_vector_d, c_vector, 2*c_sizeof(c_double)*2);

    ! Set general Jastrow data

    rc = qmckl_set_electron_num_device(context, elec_up_num, elec_dn_num);

    rc = qmckl_set_electron_coord_device(context, 'N', walk_num, elec_coord_d, walk_num*3*elec_num);
    rc = qmckl_set_nucleus_num_device(context, nucl_num);
    rc = qmckl_set_nucleus_coord_device(context, 'T', nucl_coord_d, 3*nucl_num);

    ! Initialize Jastrow specific data

    rc = qmckl_set_jastrow_aord_num_device(context, aord_num);
    rc = qmckl_set_jastrow_bord_num_device(context, bord_num);
    rc = qmckl_set_jastrow_cord_num_device(context, cord_num);
    rc = qmckl_set_jastrow_type_nucl_num_device(context, type_nucl_num);
    rc = qmckl_set_jastrow_type_nucl_vector_device(context, type_nucl_vector_d, nucl_num);
    rc = qmckl_set_jastrow_a_vector_device(context, a_vector_d, (aord_num+1)*type_nucl_num);
    rc = qmckl_set_jastrow_b_vector_device(context, b_vector_d, (bord_num+1));
    rc = qmckl_set_jastrow_c_vector_device(context, c_vector_d, dim_cord_vec*type_nucl_num);

    rc = qmckl_set_jastrow_rescale_factor_en_device(context, rescale_factor_en_d, type_nucl_num);
    rc = qmckl_set_jastrow_rescale_factor_ee_device(context, rescale_factor_ee);

    ! TODO Get factors & total Jastrow
    factor_ee_d = qmckl_malloc_device(context, walk_num*c_sizeof(c_double)*2);
    factor_en_d = qmckl_malloc_device(context, walk_num*c_sizeof(c_double)*2);
    factor_een_d = qmckl_malloc_device(context, walk_num*c_sizeof(c_double)*2);
    total_j_d = qmckl_malloc_device(context, walk_num*c_sizeof(c_double)*2);

    factor_ee_h = qmckl_malloc_host(context, walk_num*c_sizeof(c_double)*2);
    factor_en_h = qmckl_malloc_host(context, walk_num*c_sizeof(c_double)*2);
    factor_een_h = qmckl_malloc_host(context, walk_num*c_sizeof(c_double)*2);
    total_j_h = qmckl_malloc_host(context, walk_num*c_sizeof(c_double)*2);

    rc = qmckl_get_jastrow_factor_ee_device(context, factor_ee_d, walk_num);
    rc = qmckl_get_jastrow_factor_en_device(context, factor_en_d, walk_num);
    rc = qmckl_get_jastrow_factor_een_device(context, factor_een_d, walk_num);
    rc = qmckl_get_jastrow_value_device(context, total_j_d, walk_num);

    rc = qmckl_memcpy_D2H(context, factor_ee_h, factor_ee_d, walk_num*c_sizeof(c_double)*2);
    rc = qmckl_memcpy_D2H(context, factor_en_h, factor_en_d, walk_num*c_sizeof(c_double)*2);
    rc = qmckl_memcpy_D2H(context, factor_een_h, factor_een_d, walk_num*c_sizeof(c_double)*2);
    rc = qmckl_memcpy_D2H(context, total_j_h, total_j_d, walk_num*c_sizeof(c_double)*2);


    call c_f_pointer(factor_ee_h, factor_ee, [int(walk_num, kind(4))]);
    call c_f_pointer(factor_en_h, factor_en, [int(walk_num, kind(4))]);
    call c_f_pointer(factor_een_h, factor_een, [int(walk_num, kind(4))]);
    call c_f_pointer(total_j_h, total_j, [int(walk_num, kind(4))]);

    ! Check values
    do i = 1, walk_num
        if (abs(total_j(i) - factor_ee(i) - factor_en(i) - factor_een(i)) > 1e-12) then
                   print *, "Error at (i)=", i
                   print *, "factor_ee =", factor_ee(i)
                   print *, "factor_en =", factor_en(i)
                   print *, "factor_een =", factor_een(i)
                   print *, "total_j ", total_j(i)
                   call exit(1)
        end if
    end do

end program
