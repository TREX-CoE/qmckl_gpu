!
!    ------------------------------------------
!     QMCkl - Quantum Monte Carlo kernel library
!     ------------------------------------------
!    
!     Documentation : https://trex-coe.github.io/qmckl
!     Issues        : https://github.com/trex-coe/qmckl/issues
!    
!     BSD 3-Clause License
!     
!     Copyright (c) 2020, TREX Center of Excellence
!     All rights reserved.
!     
!     Redistribution and use in source and binary forms, with or without
!     modification, are permitted provided that the following conditions are met:
!     
!     1. Redistributions of source code must retain the above copyright notice, this
!        list of conditions and the following disclaimer.
!     
!     2. Redistributions in binary form must reproduce the above copyright notice,
!        this list of conditions and the following disclaimer in the documentation
!        and/or other materials provided with the distribution.
!     
!     3. Neither the name of the copyright holder nor the names of its
!        contributors may be used to endorse or promote products derived from
!        this software without specific prior written permission.
!     
!     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!     AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!     IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!     FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!     DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!     OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!     
!     
!    
!    
!
module qmckl_gpu
  use, intrinsic :: iso_c_binding
integer  , parameter :: qmckl_context = c_int64_t
integer  , parameter :: qmckl_context_device = c_int64_t
integer*8, parameter :: QMCKL_NULL_CONTEXT = 0
integer  , parameter :: qmckl_exit_code = c_int32_t

integer(qmckl_exit_code), parameter :: QMCKL_SUCCESS                  = 0
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_1            = 1
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_2            = 2
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_3            = 3
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_4            = 4
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_5            = 5
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_6            = 6
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_7            = 7
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_8            = 8
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_9            = 9
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_10           = 10
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_11           = 11
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_12           = 12
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_13           = 13
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_14           = 14
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_15           = 15
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_16           = 16
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_17           = 17
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_18           = 18
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_19           = 19
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_20           = 20
integer(qmckl_exit_code), parameter :: QMCKL_FAILURE                  = 101
integer(qmckl_exit_code), parameter :: QMCKL_ERRNO                    = 102
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_CONTEXT          = 103
integer(qmckl_exit_code), parameter :: QMCKL_ALLOCATION_FAILED        = 104
integer(qmckl_exit_code), parameter :: QMCKL_DEALLOCATION_FAILED      = 105
integer(qmckl_exit_code), parameter :: QMCKL_NOT_PROVIDED             = 106
integer(qmckl_exit_code), parameter :: QMCKL_OUT_OF_BOUNDS            = 107
integer(qmckl_exit_code), parameter :: QMCKL_ALREADY_SET              = 108
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_EXIT_CODE        = 109
! Fortran interface


interface
   integer (qmckl_context_device) function qmckl_context_create_omp_device(device_id) bind(C)
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

    integer  (qmckl_context_device) , intent(in)  , value :: context
    integer  (c_int64_t) , intent(in)  , value :: size_max
    character(c_char   ) , intent(in)          :: file_name(size_max)

  end function qmckl_trexio_read_omp_device
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_ao_num_omp_device(context, &
       num) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (qmckl_context_device) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: num
  end function qmckl_set_ao_basis_ao_num
end interface

interface
   integer(c_int32_t) function qmckl_get_ao_basis_ao_vgl_omp_device (context, &
        ao_vgl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context_device) , intent(in)  , value :: context
     double precision,     intent(out)         :: ao_vgl(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_ao_basis_ao_vgl
end interface

interface
  integer(c_int32_t) function qmckl_set_electron_coord_omp_device(context, transp, walk_num, coord, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context_device) , intent(in)  , value :: context
    character           , intent(in)  , value :: transp
    integer (c_int64_t) , intent(in)  , value :: walk_num
    double precision    , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_set_point_omp_device(context, &
       transp, num, coord, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context_device) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    integer (c_int64_t) , intent(in)  , value :: num
    real    (c_double ) , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

end module qmckl_gpu

