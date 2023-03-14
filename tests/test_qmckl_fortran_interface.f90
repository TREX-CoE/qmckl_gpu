include '../include/qmckl_gpu_f.f90'

program qmckl_test_fortran_interface


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

        integer(kind=c_int64_t) function get_ao_num() &
            bind(C, name="get_ao_num")
            use, intrinsic :: iso_c_binding
            import
            implicit none

        end function get_ao_num

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

    ! For now, this test does nothing and  only ensures that the interface compiles
    ! We might want to add simple calculations in the future
    print *, "Interface compiles"


end program
