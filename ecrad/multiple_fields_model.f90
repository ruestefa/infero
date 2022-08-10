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

  ! Command line arguments
  character(1024) :: model_path_sw,model_path_lw
  character(1024) :: model_type
  character(1024) :: yaml_config

  ! model of infero model
  type(infero_model) :: model_sw,model_lw

  integer,parameter :: batch_size = 1000

  ! input and output tensors
  real(c_float) :: input_3d(batch_size,60,4,6)
  real(c_float) :: input_2d(batch_size,4,8)
  real(c_float) :: pred_swflx(batch_size,60,2), pred_lwflx(batch_size,60,2)
  type(infero_tensor_set) :: iset
  type(infero_tensor_set) :: oset_lw, oset_sw

  ! netcdf
  character(1024) :: nc_name,varname,icon_grid
  character(50) :: dim_name(6),grid_dim_name(14)
  character(19) :: timestamp(4)
  integer :: dim_len(6),grid_dim_len(14)

  ! indices from notebook
  integer :: step, idx,n_step, s_idx,e_idx, step_idx(4)

  real(c_float), allocatable :: from_netcdf_3d(:,:,:,:), from_netcdf_2d(:,:,:), neighbor_cell_index(:,:)
  real(c_float), allocatable :: swflx(:,:,:,:), lwflx(:,:,:,:)

  ! Get CL arguments
  CALL get_command_argument(1, model_path_lw)
  CALL get_command_argument(2, model_path_sw)
  CALL get_command_argument(3, model_type)

  ! 0) init infero
  call infero_check(infero_initialise())

  ! YAML config string -> LW
  yaml_config = "---"//NEW_LINE('A') &
    //"  path: "//TRIM(model_path_lw)//NEW_LINE('A') &
    //"  type: "//TRIM(model_type)//c_null_char

  ! get a infero model
  call infero_check(model_lw%initialise_from_yaml_string(yaml_config))

  ! YAML config string -> SW
  yaml_config = "---"//NEW_LINE('A') &
    //"  path: "//TRIM(model_path_sw)//NEW_LINE('A') &
    //"  type: "//TRIM(model_type)//c_null_char

  ! get a infero model
  call infero_check(model_sw%initialise_from_yaml_string(yaml_config))

  ! init tensor sets
  call infero_check(iset%initialise())
  call infero_check(oset_lw%initialise())
  call infero_check(oset_sw%initialise())

  ! input is identical for both models
  call infero_check(iset%push_tensor(input_3d, "serving_default_input_3d"))
  call infero_check(iset%push_tensor(input_2d, "serving_default_input_2d"))

  ! LW
  call infero_check(oset_lw%push_tensor(pred_lwflx, "StatefulPartitionedCall"))

  ! SW
  call infero_check(oset_sw%push_tensor(pred_swflx, "StatefulPartitionedCall"))

  n_step = 4
  idx=1

  timestamp(1) = "2000-12-14 00:00:00"
  timestamp(2) = "2000-12-14 06:00:00"
  timestamp(3) = "2000-12-14 12:00:00"
  timestamp(4) = "2000-12-14 18:00:00"

  ALLOCATE(swflx(batch_size,60,2,n_step))
  ALLOCATE(lwflx(batch_size,60,2,n_step))

  
  ! strange behaviour with init to 0 -> all values stay 0 even after infer
  !pred_swflx(1,:,:) = 0
  !pred_lwflx(1,:,:) = 0

  nc_name='/scratch/snx3000/juckerj/WG_1_tasks/ecrad/regress-radiation/&
        &online-datasets/datasets/APE_JABW_R02B05_G_ECRAD_icongrid_atm_3d_DOM01_ml_0024_lonlat.nc'

  icon_grid='/scratch/snx3000/juckerj/WG_1_tasks/ecrad/grids/icon_grid_0008_R02B05_G.nc'

  ! icon grid
  call get_nc_dims(icon_grid,grid_dim_name,grid_dim_len,14)

  ALLOCATE(neighbor_cell_index(grid_dim_len(1),grid_dim_len(5)))

  varname="neighbor_cell_index"
  call readgrid_1d(icon_grid,varname,neighbor_cell_index(:,:),grid_dim_len(1),grid_dim_len(5))

  ! data
  call get_nc_dims(nc_name,dim_name,dim_len,6)

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


  s_idx = idx
  e_idx = idx + batch_size
  step_idx(1) = 1
  step_idx(2) = 3
  step_idx(3) = 5
  step_idx(4) = 7
  DO step=1,n_step
    input_2d(:,1,:) = from_netcdf_2d(s_idx:e_idx,step_idx(step),:)
    input_2d(:,2,:) = from_netcdf_2d(INT(neighbor_cell_index(s_idx:e_idx,1)),step_idx(step),:)
    input_2d(:,3,:) = from_netcdf_2d(INT(neighbor_cell_index(s_idx:e_idx,2)),step_idx(step),:)
    input_2d(:,4,:) = from_netcdf_2d(INT(neighbor_cell_index(s_idx:e_idx,3)),step_idx(step),:)

    input_3d(:,:,1,:) = from_netcdf_3d(s_idx:e_idx,:,step,:)
    input_3d(:,:,2,:) = from_netcdf_3d(INT(neighbor_cell_index(s_idx:e_idx,1)),:,step_idx(step),:)
    input_3d(:,:,3,:) = from_netcdf_3d(INT(neighbor_cell_index(s_idx:e_idx,2)),:,step_idx(step),:)
    input_3d(:,:,4,:) = from_netcdf_3d(INT(neighbor_cell_index(s_idx:e_idx,3)),:,step_idx(step),:)


    ! apply model
    call infero_check(model_lw%infer(iset, oset_lw ))
    call infero_check(model_sw%infer(iset, oset_sw ))

    swflx(:,:,:,step) = pred_swflx
    lwflx(:,:,:,step) = pred_lwflx
  ENDDO

! free the model
  call infero_check(model_lw%free())
  call infero_check(model_sw%free())

! finalise infero library
  call infero_check(infero_finalise())


  print*, 'Statistics:'
  DO step=1,n_step
    print*, timestamp(step)
    print*, 'SW (all levels): ',MAXVAL(swflx(:,:,:,step)), MINVAL(swflx(:,:,:,step))
    print*, 'LW (all levels): ',MAXVAL(lwflx(:,:,:,step)), MINVAL(lwflx(:,:,:,step))
    print*, 'SW (sfc): ',MAXVAL(swflx(:,1,:,step)), MINVAL(swflx(:,1,:,step))
    print*, 'LW (sfc): ',MAXVAL(lwflx(:,1,:,step)), MINVAL(lwflx(:,1,:,step))
    print*, 'SW (toa): ',MAXVAL(swflx(:,60,:,step)), MINVAL(swflx(:,60,:,step))
    print*, 'LW (toa): ',MAXVAL(lwflx(:,60,:,step)), MINVAL(lwflx(:,60,:,step))
    print*, ''
  ENDDO


end program

SUBROUTINE get_nc_dims(infile,dim_name,dim_len,nr_dims)
  use netcdf

  INTEGER(KIND=4), INTENT(OUT) :: dim_len(nr_dims)
  CHARACTER(LEN=50), INTENT(OUT) :: dim_name(nr_dims)
  INTEGER(KIND=4) :: ncid,i
  CHARACTER(LEN=1024), INTENT(IN) :: infile
  !Open netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_open(TRIM(infile), nf90_nowrite, ncid))
  !Inquire about the dimensions
  !:-------:-------:-------:-------:-------:-------:-------:

  print *, 'Dimension of NetCDF: ',TRIM(infile)
  DO i=1,nr_dims
    CALL check(nf90_inquire_dimension(ncid,i,dim_name(i),dim_len(i)))
    print *, '  ',TRIM(dim_name(i)), ':', dim_len(i)
  ENDDO
  !Close netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_close(ncid))
END SUBROUTINE get_nc_dims

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
