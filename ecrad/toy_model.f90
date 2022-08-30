!
! (C) Copyright 1996- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

program my_program
use inferof
use iso_c_binding, only : c_double, c_int, c_float, c_char, c_null_char, c_ptr
implicit none

integer, parameter :: n_inference_reps = 10
integer :: i

! Command line arguments
character(1024) :: model_path
character(1024) :: model_type
character(1024) :: yaml_config

! model of infero model
type(infero_model) :: model

! input and output tensors
real(c_float) :: it2f(1,60,1,1)
real(c_float) :: it3f(1,1,1)
real(c_float) :: ot2f(1,60,1) = 0
type(infero_tensor_set) :: iset
type(infero_tensor_set) :: oset

call random_number(it2f)
call random_number(it3f)

call infero_check(iset%initialise())
call infero_check(iset%push_tensor(it2f, "serving_default_input_3d"))
call infero_check(iset%push_tensor(it3f, "serving_default_input_2d"))

call infero_check(oset%initialise())
call infero_check(oset%push_tensor(ot2f, "StatefulPartitionedCall"))

call infero_check(iset%print())
call infero_check(oset%print())

! Get CL arguments
CALL get_command_argument(1, model_path)
CALL get_command_argument(2, model_type)

! 0) init infero
call infero_check(infero_initialise())


! YAML config string
yaml_config = "---"//NEW_LINE('A') &
  //"  path: "//TRIM(model_path)//NEW_LINE('A') &
  //"  type: "//TRIM(model_type)//c_null_char

! get a infero model
  call infero_check(model%initialise_from_yaml_string(yaml_config))

! run inference
do i=1,n_inference_reps
  call infero_check(model%infer(iset, oset ))
end do

! free the model
  call infero_check(model%free())

! finalise infero library
  call infero_check(infero_finalise())

end program

