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
use iso_c_binding, only : c_float, c_null_char
implicit none

integer, parameter :: n_inference_reps = 1
integer :: i

! Command line arguments
character(1024) :: model_path
character(1024) :: model_type
character(1024) :: yaml_config

! model of infero model
type(infero_model) :: model

! input and output tensors
real(c_float) :: input_3d(1,60,4,6)
real(c_float) :: input_2d(1,4,8)
real(c_float) :: pred_swflx_dn(1,60,2) = 0
type(infero_tensor_set) :: iset
type(infero_tensor_set) :: oset

! netcdf
character(1024) :: nc_name,varname,icon_grid
character(50) :: dim_name(6),grid_dim_name(14)
integer :: dim_len(6),grid_dim_len(14)

! indices from notebook
integer :: x_cell_index=1234,x_time_index=5

real(c_float), allocatable :: from_netcdf_3d(:,:,:,:), from_netcdf_2d(:,:,:),neighbor_cell_index(:,:)

nc_name='/scratch/snx3000/juckerj/WG_1_tasks/ecrad/regress-radiation/&
      &online-datasets/datasets/APE_JABW_R02B05_G_ECRAD_icongrid_atm_3d_DOM01_ml_0024_lonlat.nc'

icon_grid='/scratch/snx3000/juckerj/WG_1_tasks/ecrad/grids/icon_grid_0008_R02B05_G.nc'

call icon_grid_dims(icon_grid,grid_dim_name,grid_dim_len)
ALLOCATE(neighbor_cell_index(grid_dim_len(1),grid_dim_len(5)))

varname="neighbor_cell_index"
call readgrid_1d(icon_grid,varname,neighbor_cell_index(:,:),grid_dim_len(1),grid_dim_len(5))

call griddims(nc_name,dim_name,dim_len)

ALLOCATE(from_netcdf_3d(dim_len(1),dim_len(3),dim_len(6),6))
ALLOCATE(from_netcdf_2d(dim_len(1),dim_len(6),8))

! 2D FIELDS
varname="pres_sfc"
call readgrid_1d(nc_name,varname,from_netcdf_2d(:,:,1),dim_len(1),dim_len(6))

varname="cosmu0"
call readgrid_1d(nc_name,varname,from_netcdf_2d(:,:,2),dim_len(1),dim_len(6))

varname="qv_s"
call readgrid_1d(nc_name,varname,from_netcdf_2d(:,:,3),dim_len(1),dim_len(6))

varname="albvisdir"
call readgrid_1d(nc_name,varname,from_netcdf_2d(:,:,4),dim_len(1),dim_len(6))

varname="albnirdir"
call readgrid_1d(nc_name,varname,from_netcdf_2d(:,:,5),dim_len(1),dim_len(6))

varname="tsfctrad"
call readgrid_1d(nc_name,varname,from_netcdf_2d(:,:,6),dim_len(1),dim_len(6))

varname="albvisdif"
call readgrid_1d(nc_name,varname,from_netcdf_2d(:,:,7),dim_len(1),dim_len(6))

varname="albnirdif"
call readgrid_1d(nc_name,varname,from_netcdf_2d(:,:,8),dim_len(1),dim_len(6))

! 3D FIELDS
varname="clc"
call readgrid_2d(nc_name,varname,from_netcdf_3d(:,:,:,1),dim_len(1),dim_len(3),dim_len(6))

varname="temp"
call readgrid_2d(nc_name,varname,from_netcdf_3d(:,:,:,2),dim_len(1),dim_len(3),dim_len(6))

varname="pres"
call readgrid_2d(nc_name,varname,from_netcdf_3d(:,:,:,3),dim_len(1),dim_len(3),dim_len(6))

varname="qc"
call readgrid_2d(nc_name,varname,from_netcdf_3d(:,:,:,4),dim_len(1),dim_len(3),dim_len(6))

varname="qi"
call readgrid_2d(nc_name,varname,from_netcdf_3d(:,:,:,5),dim_len(1),dim_len(3),dim_len(6))

varname="qv"
call readgrid_2d(nc_name,varname,from_netcdf_3d(:,:,:,6),dim_len(1),dim_len(3),dim_len(6))

! check values from notebook
!print*, cosmu0(x_cell_index+1,x_time_index+1)
!print*, MAXVAL(temp(x_cell_index+1,:,x_time_index+1))

input_2d(1,1,:) = from_netcdf_2d(x_cell_index+1,x_time_index+1,:)
input_2d(1,2,:) = from_netcdf_2d(INT(neighbor_cell_index(x_cell_index+1,1)),x_time_index+1,:)
input_2d(1,3,:) = from_netcdf_2d(INT(neighbor_cell_index(x_cell_index+1,2)),x_time_index+1,:)
input_2d(1,4,:) = from_netcdf_2d(INT(neighbor_cell_index(x_cell_index+1,3)),x_time_index+1,:)

input_3d(1,:,1,:) = from_netcdf_3d(x_cell_index+1,:,x_time_index+1,:)
input_3d(1,:,2,:) = from_netcdf_3d(x_cell_index+1,:,x_time_index+1,:)
input_3d(1,:,3,:) = from_netcdf_3d(x_cell_index+1,:,x_time_index+1,:)
input_3d(1,:,4,:) = from_netcdf_3d(x_cell_index+1,:,x_time_index+1,:)

call infero_check(iset%initialise())
call infero_check(iset%push_tensor(input_3d, "serving_default_input_3d"))
call infero_check(iset%push_tensor(input_2d, "serving_default_input_2d"))

call infero_check(oset%initialise())
call infero_check(oset%push_tensor(pred_swflx_dn, "StatefulPartitionedCall"))

!call infero_check(iset%print())
!call infero_check(oset%print())

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

  ! check values from notebook
  print*, MAXVAL(pred_swflx_dn)

end program

SUBROUTINE icon_grid_dims(infile,dim_name,dim_len)
  use netcdf

  INTEGER(KIND=4), INTENT(OUT) :: dim_len(14)
  CHARACTER(LEN=50), INTENT(OUT) :: dim_name(14)
  INTEGER(KIND=4) :: ncid,i
  CHARACTER(LEN=1024), INTENT(IN) :: infile
  !Open netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_open(TRIM(infile), nf90_nowrite, ncid))
  !Inquire about the dimensions
  !:-------:-------:-------:-------:-------:-------:-------:
  !CALL check(nf90_inquire_dimension(ncid,1,dim_name(1),dim_len(1)))
  !CALL check(nf90_inquire_dimension(ncid,2,dim_name(2),dim_len(2)))
  !CALL check(nf90_inquire_dimension(ncid,3,dim_name(3),dim_len(3)))
  !CALL check(nf90_inquire_dimension(ncid,4,dim_name(4),dim_len(4)))
  !CALL check(nf90_inquire_dimension(ncid,5,dim_name(5),dim_len(5)))
  !CALL check(nf90_inquire_dimension(ncid,6,dim_name(6),dim_len(6)))
  !CALL check(nf90_inquire_dimension(ncid,7,dim_name(7),dim_len(7)))
  !CALL check(nf90_inquire_dimension(ncid,8,dim_name(8),dim_len(8)))

  print *, 'Dimension of NetCDF:'
  DO i=1,14
    CALL check(nf90_inquire_dimension(ncid,i,dim_name(i),dim_len(i)))
    print *, '  ',TRIM(dim_name(i)), ':', dim_len(i)
  ENDDO
  !Close netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_close(ncid))
END SUBROUTINE icon_grid_dims

SUBROUTINE griddims(infile,dim_name,dim_len)
  use netcdf

  INTEGER(KIND=4), INTENT(OUT) :: dim_len(6)
  CHARACTER(LEN=50), INTENT(OUT) :: dim_name(6)
  INTEGER(KIND=4) :: ncid,i
  CHARACTER(LEN=1024), INTENT(IN) :: infile
  !Open netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_open(TRIM(infile), nf90_nowrite, ncid))
  !Inquire about the dimensions
  !:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_inquire_dimension(ncid,1,dim_name(1),dim_len(1)))
  CALL check(nf90_inquire_dimension(ncid,2,dim_name(2),dim_len(2)))
  CALL check(nf90_inquire_dimension(ncid,3,dim_name(3),dim_len(3)))
  CALL check(nf90_inquire_dimension(ncid,4,dim_name(4),dim_len(4)))
  CALL check(nf90_inquire_dimension(ncid,5,dim_name(5),dim_len(5)))
  CALL check(nf90_inquire_dimension(ncid,6,dim_name(6),dim_len(6)))

  print *, 'Dimension of NetCDF:'
  DO i=1,6
    print *, '  ',TRIM(dim_name(i)), ':', dim_len(i)
  ENDDO
  !Close netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_close(ncid))
END SUBROUTINE griddims

SUBROUTINE readgrid_2d(infile,varname,idata,nx,ny,nz)
  use iso_c_binding, only : c_float
  use netcdf

  REAL(c_float), DIMENSION(nx,ny,nz), INTENT(OUT) :: idata
  INTEGER(KIND=4), INTENT(IN) :: nx,ny,nz
  INTEGER(KIND=4), DIMENSION(3) :: dimids
  INTEGER(KIND=4) :: ncid, ndims, varid
  CHARACTER(LEN=1024), INTENT(IN) :: infile,varname
  !Open netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_open(TRIM(infile), nf90_nowrite, ncid))

  !Get the values of the coordinates and put them in xpos & ypos
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  !CALL check(nf90_inquire_variable(ncid,1,vname,xtype,ndims,dimids))
  CALL check(nf90_inq_varid(ncid,TRIM(varname),varid))
  CALL check(nf90_inquire_variable(ncid=ncid,varid=varid,ndims=ndims))
  CALL check(nf90_inquire_variable(ncid=ncid,varid=varid,ndims=ndims,dimids=dimids))

  print *, TRIM(varname)
  PRINT *, '   ndims:',ndims
  PRINT *, '   dimids:',dimids

  CALL check(nf90_get_var(ncid,varid,idata))
  !Close netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_close(ncid))

END SUBROUTINE readgrid_2d

SUBROUTINE readgrid_1d(infile,varname,idata,nx,ny)
  use iso_c_binding, only : c_float
  use netcdf

  REAL(c_float), DIMENSION(nx,ny), INTENT(OUT) :: idata
  INTEGER(KIND=4), INTENT(IN) :: nx,ny
  INTEGER(KIND=4), DIMENSION(2) :: dimids
  INTEGER(KIND=4) :: ncid, ndims, varid
  CHARACTER(LEN=1024), INTENT(IN) :: infile,varname
  !Open netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_open(TRIM(infile), nf90_nowrite, ncid))
  CALL check(nf90_inq_varid(ncid,TRIM(varname),varid))
  CALL check(nf90_inquire_variable(ncid=ncid,varid=varid,ndims=ndims,dimids=dimids))

  print *, TRIM(varname)
  PRINT *, '   ndims:',ndims
  PRINT *, '   dimids:',dimids

  CALL check(nf90_get_var(ncid,varid,idata))
  !Close netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_close(ncid))

END SUBROUTINE readgrid_1d

SUBROUTINE check(istatus)
  use netcdf

  INTEGER, INTENT (IN) :: istatus
  IF (istatus /= nf90_noerr) THEN
  write(*,*) TRIM((nf90_strerror(istatus)))
  END IF
END SUBROUTINE check
!:=========================================================================

!END MODULE read_file
