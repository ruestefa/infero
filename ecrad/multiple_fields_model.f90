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

  !integer,parameter :: batch_size = 100
  integer,parameter :: batch_size = 1000
  !integer,parameter :: batch_size = 81919
  REAL (c_float), PARAMETER ::  pi = 3.14159265358979323846264338327950288
  REAL (c_float), PARAMETER ::  rad2deg   = 180.0/pi

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
  integer :: step, idx,n_step, s_idx,e_idx, step_idx(4),i,k,id_d, counter(4)

  real(c_float), allocatable :: from_netcdf_3d(:,:,:,:), from_netcdf_2d(:,:,:), neighbor_cell_index(:,:), diff_to_ref(:,:,:,:)
  real(c_float), allocatable :: swflx(:,:,:,:), lwflx(:,:,:,:)
  real(c_float) :: percentage, lat_min,lat_max,lon_min,lon_max
  real(c_float), allocatable :: clat(:), clon(:),lon(:),lat(:)

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
  ALLOCATE(clat(grid_dim_len(1)))
  ALLOCATE(clon(grid_dim_len(1)))
  ALLOCATE(lat(batch_size + 1))
  ALLOCATE(lon(batch_size + 1))

  varname="neighbor_cell_index"
  call readgrid_1d(icon_grid,varname,neighbor_cell_index(:,:),grid_dim_len(1),grid_dim_len(5))

  varname="clat"
  call readgrid_0d(icon_grid,varname,clat,grid_dim_len(1))
  varname="clon"
  call readgrid_0d(icon_grid,varname,clon,grid_dim_len(1))

  ! compute domain extent
  DO i=1,batch_size+1
    lon(i) = rad2deg * clon(i)
    lat(i) = rad2deg * clat(i)
  ENDDO

  lon_max = MAXVAL(lon(:))
  lon_min = MINVAL(lon(:))
  lat_max = MAXVAL(lat(:))
  lat_min = MINVAL(lat(:))



  write(*,'(A,I7,A)') 'Domain extent for batch_size', batch_size,':'
  write(*,'(A,F8.2,A,F8.2)') '  latmax:', lat_max, ' latmin', lat_min
  write(*,'(A,F8.2,A,F8.2)') '  lonmax:', lon_max, ' lonmin', lon_min

  ! data
  call get_nc_dims(nc_name,dim_name,dim_len,6)

  ALLOCATE(from_netcdf_3d(dim_len(1),dim_len(3),dim_len(6),10))
  ALLOCATE(diff_to_ref(dim_len(1),dim_len(3),dim_len(6),4))
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

  varname="lwflx_up"
  call readgrid_2d(nc_name,varname,from_netcdf_3d(:,:,:,7),dim_len(1),dim_len(3),dim_len(6))

  varname="lwflx_dn"
  call readgrid_2d(nc_name,varname,from_netcdf_3d(:,:,:,8),dim_len(1),dim_len(3),dim_len(6))

  varname="swflx_up"
  call readgrid_2d(nc_name,varname,from_netcdf_3d(:,:,:,9),dim_len(1),dim_len(3),dim_len(6))

  varname="swflx_dn"
  call readgrid_2d(nc_name,varname,from_netcdf_3d(:,:,:,10),dim_len(1),dim_len(3),dim_len(6))


  s_idx = idx
  e_idx = idx + batch_size

  ! cover complete diurnal solar cycle
  step_idx(1) = 1
  step_idx(2) = 3
  step_idx(3) = 5
  step_idx(4) = 7

  DO step=1,n_step

    input_2d(:,1,:) = from_netcdf_2d(s_idx:e_idx,step_idx(step),:)

    ! assign neigbours 2D
    input_2d(:,2,:) = from_netcdf_2d(INT(neighbor_cell_index(s_idx:e_idx,1)),step_idx(step),:)
    input_2d(:,3,:) = from_netcdf_2d(INT(neighbor_cell_index(s_idx:e_idx,2)),step_idx(step),:)
    input_2d(:,4,:) = from_netcdf_2d(INT(neighbor_cell_index(s_idx:e_idx,3)),step_idx(step),:)

    ! assign neigbours 3D
    input_3d(:,:,1,:) = from_netcdf_3d(s_idx:e_idx,:,step_idx(step),:)
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

  counter(:) = 0
  DO step=1,n_step
    DO i=1,batch_size
      DO k=1,60
        DO id_d=1,2
          IF(swflx(i,k,id_d,step) < 0.0) THEN
            counter(step) = counter(step) + 1
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  write(*,'(A,I7,A)') 'Domain extent for batch_size', batch_size,':'
  write(*,'(A,F8.2,A,F8.2)') '  latmax:', lat_max, ' latmin', lat_min
  write(*,'(A,F8.2,A,F8.2)') '  lonmax:', lon_max, ' lonmin', lon_min

  write(*,'(I7,A)') batch_size*60*2,' datapoints scanned for each step of SW-flux'
  DO step=1,n_step
    percentage = MAX(0.0,100.0 * REAL(counter(step))/REAL((batch_size*60*2)))
    write(*,'(F6.1,A,A)') percentage, '% values below 0.0 for ',timestamp(step)
  ENDDO

  DO step=1,n_step
    ! lwflux_up
    diff_to_ref(:,:,1,step) = ABS(lwflx(:,:,1,step) - from_netcdf_3d(s_idx:e_idx,:,step_idx(step),7))
    ! lwflux_dn
    diff_to_ref(:,:,2,step) = ABS(lwflx(:,:,2,step) - from_netcdf_3d(s_idx:e_idx,:,step_idx(step),8))
    ! swflux_up
    diff_to_ref(:,:,3,step) = ABS(swflx(:,:,1,step) - from_netcdf_3d(s_idx:e_idx,:,step_idx(step),9))
    ! swflux_dn
    diff_to_ref(:,:,4,step) = ABS(swflx(:,:,2,step) - from_netcdf_3d(s_idx:e_idx,:,step_idx(step),10))
  ENDDO


  write(*,'(A)') ''
  write(*,'(A)') 'Statistics'
  DO step=1,n_step
    write(*,'(A)') ''
    write(*,'(A,A)') '  ',timestamp(step)
    write(*,'(A)') '    Upwards:'
    write(*,'(A,F8.2,F8.2)') '     SW (all levels): ',MAXVAL(swflx(:,:,1,step)), MINVAL(swflx(:,:,1,step))
    write(*,'(A,F8.2,F8.2)') '     LW (all levels): ',MAXVAL(lwflx(:,:,1,step)), MINVAL(lwflx(:,:,1,step))
    write(*,'(A,F8.2,F8.2)') '     SW (sfc): ',MAXVAL(swflx(:,1,1,step)), MINVAL(swflx(:,1,1,step))
    write(*,'(A,F8.2,F8.2)') '     LW (sfc): ',MAXVAL(lwflx(:,1,1,step)), MINVAL(lwflx(:,1,1,step))
    write(*,'(A,F8.2,F8.2)') '     SW (toa): ',MAXVAL(swflx(:,60,1,step)), MINVAL(swflx(:,60,1,step))
    write(*,'(A,F8.2,F8.2)') '     LW (toa): ',MAXVAL(lwflx(:,60,1,step)), MINVAL(lwflx(:,60,1,step))
    write(*,'(A)') ''
    write(*,'(A)') '    Downwards:'
    write(*,'(A,F8.2,F8.2)') '     SW (all levels): ',MAXVAL(swflx(:,:,2,step)), MINVAL(swflx(:,:,2,step))
    write(*,'(A,F8.2,F8.2)') '     LW (all levels): ',MAXVAL(lwflx(:,:,2,step)), MINVAL(lwflx(:,:,2,step))
    write(*,'(A,F8.2,F8.2)') '     SW (sfc): ',MAXVAL(swflx(:,1,2,step)), MINVAL(swflx(:,1,2,step))
    write(*,'(A,F8.2,F8.2)') '     LW (sfc): ',MAXVAL(lwflx(:,1,2,step)), MINVAL(lwflx(:,1,2,step))
    write(*,'(A,F8.2,F8.2)') '     SW (toa): ',MAXVAL(swflx(:,60,2,step)), MINVAL(swflx(:,60,2,step))
    write(*,'(A,F8.2,F8.2)') '     LW (toa): ',MAXVAL(lwflx(:,60,2,step)), MINVAL(lwflx(:,60,2,step))
    write(*,'(A)') ''
    write(*,'(A)') '    Absolute difference to reference'
    write(*,'(A)') '      Upwards:'
    write(*,'(A,F8.2,F8.2)') '        SW (all levels): ',MAXVAL(diff_to_ref(:,:,3,step)), MINVAL(diff_to_ref(:,:,3,step))
    write(*,'(A,F8.2,F8.2)') '        LW (all levels): ',MAXVAL(diff_to_ref(:,:,1,step)), MINVAL(diff_to_ref(:,:,1,step))
    write(*,'(A,F8.2,F8.2)') '        SW (sfc): ',MAXVAL(diff_to_ref(:,1,3,step)), MINVAL(diff_to_ref(:,1,3,step))
    write(*,'(A,F8.2,F8.2)') '        LW (sfc): ',MAXVAL(diff_to_ref(:,1,1,step)), MINVAL(diff_to_ref(:,1,1,step))
    write(*,'(A,F8.2,F8.2)') '        SW (toa): ',MAXVAL(diff_to_ref(:,60,3,step)), MINVAL(diff_to_ref(:,60,3,step))
    write(*,'(A,F8.2,F8.2)') '        LW (toa): ',MAXVAL(diff_to_ref(:,60,1,step)), MINVAL(diff_to_ref(:,60,1,step))
    write(*,'(A)') '      Downwards:'
    write(*,'(A,F8.2,F8.2)') '        SW (all levels): ',MAXVAL(diff_to_ref(:,:,4,step)), MINVAL(diff_to_ref(:,:,4,step))
    write(*,'(A,F8.2,F8.2)') '        LW (all levels): ',MAXVAL(diff_to_ref(:,:,2,step)), MINVAL(diff_to_ref(:,:,2,step))
    write(*,'(A,F8.2,F8.2)') '        SW (sfc): ',MAXVAL(diff_to_ref(:,1,4,step)), MINVAL(diff_to_ref(:,1,4,step))
    write(*,'(A,F8.2,F8.2)') '        LW (sfc): ',MAXVAL(diff_to_ref(:,1,2,step)), MINVAL(diff_to_ref(:,1,2,step))
    write(*,'(A,F8.2,F8.2)') '        SW (toa): ',MAXVAL(diff_to_ref(:,60,4,step)), MINVAL(diff_to_ref(:,60,4,step))
    write(*,'(A,F8.2,F8.2)') '        LW (toa): ',MAXVAL(diff_to_ref(:,60,2,step)), MINVAL(diff_to_ref(:,60,2,step))
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

  write(*,'(A,A)') 'Dimension of NetCDF: ',TRIM(infile)
  DO i=1,nr_dims
    CALL check(nf90_inquire_dimension(ncid,i,dim_name(i),dim_len(i)))
    write(*,'(A,A,A,I7)') '  ',TRIM(dim_name(i)), ':', dim_len(i)
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

  write(*,'(A)') TRIM(varname)
  write(*,'(A,I2)') '   ndims:',ndims
  write(*,*) '   dimids:',dimids

  CALL check(nf90_get_var(ncid,varid,idata))
  !Close netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_close(ncid))

END SUBROUTINE readgrid_2d

SUBROUTINE readgrid_0d(infile,varname,idata,nx)
  use iso_c_binding, only : c_float
  use netcdf

  REAL(c_float), DIMENSION(nx), INTENT(OUT) :: idata
  INTEGER(KIND=4), INTENT(IN) :: nx
  INTEGER(KIND=4), DIMENSION(1) :: dimids
  INTEGER(KIND=4) :: ncid, ndims, varid
  CHARACTER(LEN=1024), INTENT(IN) :: infile,varname
  !Open netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_open(TRIM(infile), nf90_nowrite, ncid))
  CALL check(nf90_inq_varid(ncid,TRIM(varname),varid))
  CALL check(nf90_inquire_variable(ncid=ncid,varid=varid,ndims=ndims,dimids=dimids))

  write(*,'(A)') TRIM(varname)
  write(*,'(A,I2)') '   ndims:',ndims
  write(*,*) '   dimids:',dimids

  CALL check(nf90_get_var(ncid,varid,idata))
  !Close netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_close(ncid))

END SUBROUTINE readgrid_0d

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

  write(*,'(A)') TRIM(varname)
  write(*,'(A,I2)') '   ndims:',ndims
  write(*,*) '   dimids:',dimids

  CALL check(nf90_get_var(ncid,varid,idata))
  !Close netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  CALL check(nf90_close(ncid))

END SUBROUTINE readgrid_1d

SUBROUTINE check(istatus)
  use netcdf

  INTEGER, INTENT (IN) :: istatus
  IF (istatus /= nf90_noerr) THEN
  write(*,'(A)') TRIM((nf90_strerror(istatus)))
  END IF
END SUBROUTINE check
!:=========================================================================

!END MODULE read_file
