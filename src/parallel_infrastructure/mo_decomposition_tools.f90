! Set of methods for grid decomposition.
! It runs only on a single process
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

!-------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------
MODULE mo_decomposition_tools
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: i8, wp
  USE mo_exception,          ONLY: finish
  USE mo_io_units,           ONLY: find_next_free_unit
  USE mo_math_utilities
  USE mo_util_sort,          ONLY: quicksort
  USE mo_impl_constants,     ONLY: success
  USE mo_dist_dir,           ONLY: t_dist_dir
  USE ppm_extents,           ONLY: extent, extent_start, extent_size

  IMPLICIT NONE

  PUBLIC :: t_glb2loc_index_lookup

  PUBLIC :: divide_cells_by_location
  PUBLIC :: t_cell_info, sort_cell_info_by_cell_number
  PUBLIC :: read_ascii_decomposition
  PUBLIC :: t_grid_domain_decomp_info
  PUBLIC :: get_local_index
  PUBLIC :: get_valid_local_index
  PUBLIC :: generate_glb_to_loc_lookup
  PUBLIC :: init_glb2loc_index_lookup
  PUBLIC :: set_inner_glb_index
  PUBLIC :: set_outer_glb_index
  PUBLIC :: deallocate_glb2loc_index_lookup
  PUBLIC :: partidx_of_elem_uniform_deco, uniform_partition, &
    &       uniform_partition_start

  PRIVATE


  !------------------------------
  TYPE t_cell_info
    INTEGER :: lat !< latitude coordinate, scaled to integer value range
    INTEGER :: lon !< longitude coordinate, scaled to integer value range
    INTEGER :: cell_number !< cell number (for back-sorting at the end)
    INTEGER :: owner !< will be set to the owner
  END TYPE t_cell_info
  !------------------------------

  TYPE t_glb2loc_index_lookup
    INTEGER, ALLOCATABLE :: inner_glb_index(:)
    INTEGER, ALLOCATABLE :: inner_glb_index_to_loc(:)
    INTEGER, ALLOCATABLE :: outer_glb_index(:)
    INTEGER, ALLOCATABLE :: outer_glb_index_to_loc(:)
    INTEGER :: global_size
  END TYPE  t_glb2loc_index_lookup

  TYPE t_grid_domain_decomp_info

    ! Owner mask:
    ! For cells this is the same as decomp_domain(:,:)==0
    ! index1=nproma, index2=1,nblks_c
    ! For edges, this can not be derived from decomp_domain:
    ! edges at the border are assigned to the PE with the bigger number
    ! index1=nproma, index2=1,nblks_e
    ! For verts, this can not be derived from decomp_domain:
    ! verts at the border are assigned the PE with the bigger number
    ! index1=nproma, index2=1,nblks_v
    LOGICAL, ALLOCATABLE :: owner_mask(:,:)

    ! The following is only used internally for the coupler
    ! and the icon_comm_lib
    INTEGER, ALLOCATABLE :: owner_local(:)

    ! The following is only used internally for domain decomposition
    INTEGER, ALLOCATABLE :: glb_index(:)
    TYPE (t_glb2loc_index_lookup) :: glb2loc_index

    ! Distributed directory containing owner information
    TYPE(t_dist_dir) :: owner_dist_dir

    ! Domain decomposition flag:
    ! decomp_domain==0: inner domain, decomp_domain>0: boundary, decomp_domain<0: undefined
    ! For cells:
    ! 0=owned, 1=shared edge with owned, 2=shared vertex with owned
    ! index1=nproma, index2=1,nblks_c
    ! For edges:
    ! 0=owned, 1=on owned cell=in domain, 2=exactly one shared vertex with owned cells
    ! index1=nproma, index2=1,nblks_e
    ! For verts:
    ! 0=owned, 1=on owned cell=in domain, 2=on level 1 cells
    ! index1=nproma, index2=1,nblks_v
    INTEGER, POINTER :: decomp_domain(:,:)

    INTEGER, POINTER :: halo_level(:,:)! just points to the decomp_domain as a more accurate name
  END TYPE

  INTERFACE get_valid_local_index
    MODULE PROCEDURE get_valid_local_index_prev
    MODULE PROCEDURE get_valid_local_index_next
  END INTERFACE

  PUBLIC :: OPERATOR(==), OPERATOR(/=)

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE t_cell_info_eq
  END INTERFACE OPERATOR(==)

  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE t_cell_info_ne
  END INTERFACE OPERATOR(/=)

  !------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_ascii_decomposition(ascii_file_name, cell_owner, no_of_cells)
    INTEGER, POINTER :: cell_owner(:)
    CHARACTER(LEN=*) :: ascii_file_name
    INTEGER, INTENT(in) :: no_of_cells

    INTEGER :: file_id, return_status, cell_no
    CHARACTER(*), PARAMETER :: method_name = "read_ascii_decomposition"

    WRITE(0,*) "Read decomposition from file: ", TRIM(ascii_file_name)
    file_id = find_next_free_unit(10,99)

    OPEN(file_id, FILE=TRIM(ascii_file_name), STATUS='OLD', IOSTAT=return_status)
    IF(return_status /= 0) CALL finish(method_name,&
      & 'Unable to open input file: '//TRIM(ascii_file_name))

    DO cell_no = 1, no_of_cells
      READ(file_id, *, IOSTAT=return_status) cell_owner(cell_no)
      IF(return_status /= 0) CALL finish(method_name,'Error reading: '//TRIM(ascii_file_name))
    ENDDO
    CLOSE(file_id)

  END SUBROUTINE read_ascii_decomposition
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Actually divides geometrically by location on cpu_a .. cpu_b
  !!
  RECURSIVE SUBROUTINE divide_cells_by_location(n_cells,cell_desc,cpu_a,cpu_b)

    INTEGER, INTENT(in) :: n_cells, cpu_a, cpu_b

    TYPE(t_cell_info), INTENT(inout) :: cell_desc(n_cells)

    INTEGER :: cpu_m, n_cells_m, i, j
    INTEGER(i8) :: lat_sum
    INTEGER :: max_lat_i, max_lon_i, min_lat_i, min_lon_i
    REAL(wp) :: max_lat, max_lon, min_lat, min_lon, avglat, scalexp
    !-----------------------------------------------------------------------

    ! If there is only 1 CPU for distribution, we are done

    IF(cpu_a==cpu_b .OR. n_cells == 0) THEN
      cell_desc(:)%owner = cpu_a
      RETURN
    ENDIF

    ! Get geometric extensions and total number of points of all patches
    min_lat_i = HUGE(min_lat_i)
    max_lat_i = -HUGE(max_lat_i)
    min_lon_i = HUGE(min_lon_i)
    max_lon_i = -HUGE(max_lon_i)
    lat_sum = 0_i8
    DO i = 1, n_cells
      min_lat_i = MIN(min_lat_i, cell_desc(i)%lat)
      max_lat_i = MAX(max_lat_i, cell_desc(i)%lat)
      min_lon_i = MIN(min_lon_i, cell_desc(i)%lon)
      max_lon_i = MAX(max_lon_i, cell_desc(i)%lon)
      lat_sum = lat_sum + INT(cell_desc(i)%lat, i8)
    END DO
    ! min_lat_i = MINVAL(cell_desc(:)%lat)
    ! min_lon_i = MINVAL(cell_desc(:)%lon)
    ! max_lat_i = MAXVAL(cell_desc(:)%lat)
    ! max_lon_i = MAXVAL(cell_desc(:)%lon)

    min_lon = flp_lon(min_lon_i)
    min_lat = flp_lat(min_lat_i)
    max_lon = flp_lon(max_lon_i)
    max_lat = flp_lat(max_lat_i)

    ! average latitude in patch
    avglat  = flp_lat(INT(lat_sum/INT(n_cells,i8)))

    ! account somehow for convergence of meridians - this formula is just empiric
    scalexp = 1._wp - MAX(0._wp,ABS(max_lat)-1._wp,ABS(min_lat)-1._wp)
    min_lon = min_lon*(COS(avglat))**scalexp
    max_lon = max_lon*(COS(avglat))**scalexp

    ! Get dimension with biggest distance from min to max
    ! and sort cells in this dimension

    IF(max_lat-min_lat >= max_lon-min_lon) THEN
      CALL sort_cell_info_by_lat(cell_desc, n_cells)
    ELSE
      CALL sort_cell_info_by_lon(cell_desc, n_cells)
    ENDIF

    ! CPU number where to split CPU set

    cpu_m = (cpu_a+cpu_b-1)/2

    ! If the number of CPUs is not even, we have to split the set
    ! of cells accordingly into to differently sized halfes
    ! in order to end with an equal number of points on every CPU.
    ! Note that DOUBLE arithmetic is used since the integer size
    ! may be exceeded in this calculation!

    n_cells_m = INT(DBLE(n_cells)*DBLE((cpu_b-cpu_a+1)/2)/DBLE(cpu_b-cpu_a+1))

    ! find cells with same lat/lon coordinate to left and right ...
    i = n_cells_m
    j = n_cells_m
    IF (max_lat-min_lat >= max_lon-min_lon) THEN
      DO WHILE (i > 1)
        IF (cell_desc(i - 1)%lat /= cell_desc(n_cells_m)%lat) EXIT
        i = i - 1
      END DO
      DO WHILE (j < n_cells)
        IF (cell_desc(j + 1)%lat /= cell_desc(n_cells_m)%lat) EXIT
        j = j + 1
      END DO
    ELSE
      DO WHILE (i > 1)
        IF (cell_desc(i - 1)%lon /= cell_desc(n_cells_m)%lon) EXIT
        i = i - 1
      END DO
      DO WHILE (j < n_cells)
        IF (cell_desc(j + 1)%lon /= cell_desc(n_cells_m)%lon) EXIT
        j = j + 1
      END DO
    END IF
    ! ... and sort those by cell number to make assignment to partition
    ! repeatable in parallel version
    CALL sort_cell_info_by_cell_number(cell_desc(i:j), j - i + 1)

    ! If there are only two CPUs, we are done

    IF(cpu_b == cpu_a+1) THEN
      cell_desc(1:n_cells_m)%owner = cpu_a
      cell_desc(n_cells_m+1:n_cells)%owner = cpu_b
      RETURN
    ENDIF

    ! Further divide both halves recursively

    CALL divide_cells_by_location(n_cells_m,cell_desc(1:n_cells_m),&
      cpu_a,cpu_m)
    CALL divide_cells_by_location(n_cells-n_cells_m,cell_desc(n_cells_m+1:n_cells),&
      cpu_m+1,cpu_b)

  END SUBROUTINE divide_cells_by_location

  !-------------------------------------------------------------------------
  !> generates the information required to get the local index for a given
  !  global index. The decomp_info%glb_index array needs to be set in order
  !  for this routine to work.
  !  The routine can be used reconstruct the respective information.
  !  Originally this information is generated by
  !  mo_setup_subdivion->divide_patch. The data generated by this routine might
  !  slightly differ from the original data...
  SUBROUTINE generate_glb_to_loc_lookup(n, n_inner, decomp_info)

    INTEGER, INTENT(IN)    :: n        ! Number of local points
    INTEGER, INTENT(IN)    :: n_inner  ! Number of inner global points
    TYPE(t_grid_domain_decomp_info), INTENT(INOUT) :: decomp_info

    INTEGER :: i

    CALL set_inner_glb_index(decomp_info%glb2loc_index, &
      &                      decomp_info%glb_index(1:n_inner), &
      &                      (/(i, i = 1, n_inner)/))
    CALL set_outer_glb_index(decomp_info%glb2loc_index, &
      &                      decomp_info%glb_index(n_inner+1:), &
      &                      (/(i, i = n_inner+1, n)/))

  END SUBROUTINE generate_glb_to_loc_lookup

  !-------------------------------------------------------------------------

  PURE FUNCTION binary_search(array, key)

    INTEGER, INTENT(IN) :: array(:), key
    INTEGER :: binary_search

    INTEGER :: lb, ub, middle

    !$ACC ROUTINE SEQ

    lb = 1
    ub = SIZE(array)
    middle = ub / 2

    IF (ub == 0) THEN
      binary_search = 0
      RETURN
    END IF

    DO WHILE (ub >= lb)

      middle = (ub + lb) / 2;

      IF (array(middle) < key) THEN
        lb = middle + 1
      ELSE IF (array(middle) > key) THEN
        ub = middle - 1
      ELSE
        EXIT
      END IF
    END DO

    IF (array(middle) == key) THEN
      binary_search = middle
    ELSE IF (array(middle) > key) THEN
      binary_search = -middle + 1
    ELSE
      binary_search = -middle
    END IF
  END FUNCTION binary_search

  !-------------------------------------------------------------------------
  ! quicksort implementations to sort an array of type(t_cell_info)
  ! according to the lat, lon and cell_number fields

  ! returns the local index for a given global index
  ! in case the global index is not available locally, -1 is returned
  ! in case the global index is invalid, 0 is returned
  ELEMENTAL FUNCTION get_local_index(glb2loc_index, glb_index)

    TYPE(t_glb2loc_index_lookup), INTENT(in) :: glb2loc_index
    INTEGER, INTENT(in) :: glb_index
    INTEGER :: get_local_index

    INTEGER :: temp

    !$ACC ROUTINE SEQ

    IF (glb_index > glb2loc_index%global_size .OR. glb_index < 1) THEN
      get_local_index = 0
    ELSE
      ! find in outer indices
      temp = binary_search(glb2loc_index%outer_glb_index(:), glb_index)
      IF (temp > 0) temp = glb2loc_index%outer_glb_index_to_loc(temp)
      ! find in inner indices
      IF (temp <= 0) THEN
        temp = binary_search(glb2loc_index%inner_glb_index(:), glb_index)
        IF (temp > 0) temp = glb2loc_index%inner_glb_index_to_loc(temp)
      END IF

      get_local_index = MERGE(temp, -1, temp > 0)
    END IF

  END FUNCTION get_local_index

  ! returns the local index for a given global index
  ! in case the global index is in the valid range but locally not
  ! available, the index of the next smaller index is return, in
  ! case there is no smaller index or the global index is
  ! invalid 0 is returned
  ELEMENTAL FUNCTION get_valid_local_index_prev(glb2loc_index, glb_index)

    TYPE(t_glb2loc_index_lookup), INTENT(in) :: glb2loc_index
    INTEGER, INTENT(in) :: glb_index
    INTEGER :: get_valid_local_index_prev

    INTEGER :: temp

    !$ACC ROUTINE SEQ

    IF (glb_index > glb2loc_index%global_size .OR. glb_index < 1) THEN

      get_valid_local_index_prev = 0
    ELSE

      ! find in outer indices
      temp = binary_search(glb2loc_index%outer_glb_index(:), glb_index)
      IF (temp > 0) THEN
        get_valid_local_index_prev = glb2loc_index%outer_glb_index_to_loc(temp)
      ELSE
        ! find in inner indices
        temp = binary_search(glb2loc_index%inner_glb_index(:), glb_index)
        IF (temp > 0) THEN
          get_valid_local_index_prev = glb2loc_index%inner_glb_index_to_loc(temp)
        ELSE
          get_valid_local_index_prev = -temp
        END IF
      END IF
    END IF
  END FUNCTION get_valid_local_index_prev

  ! returns the local index for a given global index
  ! in case the global index is in the valid range but locally not
  ! available, the index of the next greater index is return, in
  ! case there is no greater index the maximum local index + 1, in
  ! case the global index is invalid 0 is returned
  ELEMENTAL FUNCTION get_valid_local_index_next(glb2loc_index, glb_index, use_next)

    TYPE(t_glb2loc_index_lookup), INTENT(in) :: glb2loc_index
    INTEGER, INTENT(in) :: glb_index
    INTEGER :: get_valid_local_index_next
    LOGICAL, INTENT(in) :: use_next
    INTEGER :: temp


    !$ACC ROUTINE SEQ

    IF (.NOT. use_next) THEN
      get_valid_local_index_next = get_valid_local_index_prev(glb2loc_index, glb_index)
      RETURN
    ENDIF

    IF (glb_index > glb2loc_index%global_size .OR. &
      & glb_index < 1) THEN

      get_valid_local_index_next = 0
    ELSE

      ! find in outer indices
      temp = binary_search(glb2loc_index%outer_glb_index(:), glb_index)
      IF (temp > 0) THEN
        get_valid_local_index_next = glb2loc_index%outer_glb_index_to_loc(temp)
      ELSE
        ! find in inner indices
        temp = binary_search(glb2loc_index%inner_glb_index(:), glb_index)
        IF (temp > 0) THEN
          get_valid_local_index_next = glb2loc_index%inner_glb_index_to_loc(temp)
        ELSE
          get_valid_local_index_next = - temp + 1
        END IF
      END IF
    END IF
  END FUNCTION get_valid_local_index_next

  !-------------------------------------------------------------------------

  SUBROUTINE init_glb2loc_index_lookup(glb2loc, global_size)

    TYPE (t_glb2loc_index_lookup), INTENT(OUT) :: glb2loc
    INTEGER, INTENT(IN) :: global_size

    INTEGER :: ist

    ALLOCATE(glb2loc%inner_glb_index(0), &
      &      glb2loc%inner_glb_index_to_loc(0), &
      &      glb2loc%outer_glb_index(0), &
      &      glb2loc%outer_glb_index_to_loc(0), stat=ist)
    IF(ist/=success) &
      CALL finish  ("init_glb2loc_index_lookup", &
        &           'allocate in init_glb2loc_index_lookup failed')

    glb2loc%global_size = global_size

  END SUBROUTINE init_glb2loc_index_lookup

  !-------------------------------------------------------------------------

  SUBROUTINE set_inner_glb_index(glb2loc, glb_index, loc_index)

    TYPE (t_glb2loc_index_lookup), INTENT(INOUT) :: glb2loc
    INTEGER, INTENT(IN), CONTIGUOUS :: glb_index(:), loc_index(:)
    INTEGER :: ist, num_indices

    num_indices = SIZE(glb_index(:))

    IF (ANY(glb_index(:) < 1 .OR. glb_index(:) > glb2loc%global_size)) &
      CALL finish("set_inner_glb_index", "invalid global index")

    DEALLOCATE(glb2loc%inner_glb_index, &
      &        glb2loc%inner_glb_index_to_loc, stat=ist)
    IF (ist /= success) &
      CALL finish("set_inner_glb_index", "deallocate failed")

    ALLOCATE(glb2loc%inner_glb_index(num_indices), &
      &      glb2loc%inner_glb_index_to_loc(num_indices), stat=ist)
    IF (ist /= success) &
      CALL finish("set_inner_glb_index", "allocate failed")

    glb2loc%inner_glb_index(:) = glb_index(:)
    glb2loc%inner_glb_index_to_loc(:) = loc_index(1:num_indices)

    CALL quicksort(glb2loc%inner_glb_index(:), &
      &            glb2loc%inner_glb_index_to_loc(:))
  END SUBROUTINE set_inner_glb_index

  !-------------------------------------------------------------------------

  SUBROUTINE set_outer_glb_index(glb2loc, glb_index, loc_index)

    TYPE (t_glb2loc_index_lookup), INTENT(INOUT) :: glb2loc
    INTEGER, INTENT(IN), CONTIGUOUS :: glb_index(:), loc_index(:)
    INTEGER :: ist, num_indices

    num_indices = SIZE(glb_index(:))

    IF (ANY(glb_index(:) < 1 .OR. glb_index(:) > glb2loc%global_size)) &
      CALL finish("set_outer_glb_index", "invalid global index")

    DEALLOCATE(glb2loc%outer_glb_index, &
      &        glb2loc%outer_glb_index_to_loc, stat=ist)
    IF (ist /= success) &
      CALL finish("set_outer_glb_index", "deallocate failed")

    ALLOCATE(glb2loc%outer_glb_index(num_indices), &
      &      glb2loc%outer_glb_index_to_loc(num_indices), stat=ist)
    IF (ist /= success) &
      CALL finish("set_outer_glb_index", "allocate failed")

    glb2loc%outer_glb_index(:) = glb_index(:)
    glb2loc%outer_glb_index_to_loc(:) = loc_index(1:num_indices)

    CALL quicksort(glb2loc%outer_glb_index(:), &
      &            glb2loc%outer_glb_index_to_loc(:))
  END SUBROUTINE set_outer_glb_index

  !-------------------------------------------------------------------------

  SUBROUTINE deallocate_glb2loc_index_lookup(glb2loc)

    TYPE (t_glb2loc_index_lookup), INTENT(INOUT) :: glb2loc

    INTEGER :: ist

    DEALLOCATE(glb2loc%inner_glb_index, &
      &        glb2loc%inner_glb_index_to_loc, &
      &        glb2loc%outer_glb_index, &
      &        glb2loc%outer_glb_index_to_loc, stat=ist)
    IF(ist/=success) &
      CALL finish  ('deallocate_glb2loc_index_lookup', 'deallocate failed')

  END SUBROUTINE deallocate_glb2loc_index_lookup

  !-------------------------------------------------------------------------

  FUNCTION med3_cell_info_lat(a, i, j, k) RESULT(median)
    TYPE(t_cell_info), INTENT(in) :: a(0:)
    INTEGER, INTENT(in) :: i, j, k
    INTEGER :: median
    IF (a(i)%lat < a(j)%lat) THEN
      IF (a(j)%lat < a(k)%lat) THEN
        median = j
      ELSE
        IF (a(i)%lat < a(k)%lat) THEN
          median = k
        ELSE
          median = i
        END IF
      END IF
    ELSE
      IF (a(j)%lat > a(k)%lat) THEN
        median = j
      ELSE
        IF (a(i)%lat < a(k)%lat) THEN
          median = i
        ELSE
          median = k
        END IF
      END IF
    END IF
  END FUNCTION med3_cell_info_lat

  SUBROUTINE cell_info_vec_swap(a, b, n)
    INTEGER, INTENT(in) :: n
    TYPE(t_cell_info), INTENT(inout) :: a(1:n), b(1:n)
    TYPE(t_cell_info) :: t
    INTEGER :: i
    DO i = 1, n
      t = a(i)
      a(i) = b(i)
      b(i) = t
    END DO
  END SUBROUTINE cell_info_vec_swap



  RECURSIVE SUBROUTINE sort_cell_info_by_lat(a, n)
    INTEGER, INTENT(in) :: n
    TYPE(t_cell_info), INTENT(inout) :: a(0:n - 1)
    TYPE(t_cell_info) :: temp

    LOGICAL :: swap_cnt

    INTEGER :: d, pa, pb, pc, pd, pl, pm, pn, pdiff

    IF (n < 7) THEN
      DO pm = 1, n - 1
        pl = pm
        DO WHILE (pl > 0)
          IF (a(pl - 1)%lat <= a(pl)%lat) EXIT
          temp = a(pl) ; a(pl) = a( pl - 1); a( pl - 1) = temp
          pl = pl - 1
        END DO
      END DO
      RETURN
    END IF
    pm = n / 2
    IF (n > 7) THEN
      pl = 0
      pn = n - 1
      IF (n > 40) THEN
        d = n / 8
        pl = med3_cell_info_lat(a, pl, pl + d, pl + 2 * d)
        pm = med3_cell_info_lat(a, pm - d, pm, pm + d)
        pn = med3_cell_info_lat(a, pn - 2 * d, pn - d, pn)
      END IF
      pm = med3_cell_info_lat(a, pl, pm, pn)
    END IF
    temp = a(0) ; a(0) = a( pm); a( pm) = temp
    pb = 1
    pa = pb
    pd = n - 1
    pc = pd
    swap_cnt = .FALSE.
    DO WHILE (.TRUE.)
      DO WHILE (pb <= pc)
        IF (a(pb)%lat > a(0)%lat) EXIT
        IF (a(pb)%lat == a(0)%lat) THEN
          swap_cnt = .TRUE.
          temp = a(pa) ; a(pa) = a( pb); a( pb) = temp
          pa = pa + 1
        END IF
        pb = pb + 1
      END DO
      DO WHILE (pb <= pc)
        IF (a(pc)%lat < a(0)%lat) EXIT
        IF (a(pc)%lat == a(0)%lat) THEN
          swap_cnt = .TRUE.
          temp = a(pc) ; a(pc) = a( pd); a( pd) = temp
          pd = pd - 1
        END IF
        pc = pc - 1
      END DO
      IF (pb > pc) EXIT
      temp = a(pb) ; a(pb) = a( pc); a( pc) = temp
      swap_cnt = .TRUE.
      pb = pb + 1
      pc = pc - 1
    END DO
    IF (.NOT. swap_cnt) THEN  ! Switch to insertion sort
      DO pm = 1, n - 1
        pl = pm
        DO WHILE(pl > 0)
          IF (a(pl - 1)%lat <= a(pl)%lat) EXIT
          temp = a(pl) ; a(pl) = a( pl - 1); a( pl - 1) = temp
          pl = pl - 1
        END DO
      END DO
      RETURN
    END IF
    pn =  n
    pdiff = MIN(pa, pb - pa)
    IF (pdiff > 0) &
         CALL cell_info_vec_swap(a(0:pdiff - 1), a(pb - pdiff:pb-1), pdiff)
    pdiff = MIN(pd - pc, pn - pd - 1)
    if (pdiff > 0) &
         CALL cell_info_vec_swap(a(pb:pb + pdiff - 1), a(pn - pdiff:pn - 1), &
         pdiff)
    pdiff = pb - pa
    IF (pdiff > 1) &
         CALL sort_cell_info_by_lat(a, pdiff)
    pdiff = pd - pc
    ! hope the compiler can tail-recurse and save stack space
    IF (pdiff > 1) &
         CALL sort_cell_info_by_lat(a(pn - pdiff:), pdiff)
  END SUBROUTINE sort_cell_info_by_lat

  RECURSIVE SUBROUTINE sort_cell_info_by_lon(a, n)
    INTEGER, INTENT(in) :: n
    TYPE(t_cell_info), INTENT(inout) :: a(0:n - 1)
    TYPE(t_cell_info) :: temp

    LOGICAL :: swap_cnt

    INTEGER :: d, pa, pb, pc, pd, pl, pm, pn, pdiff

    IF (n < 7) THEN
      DO pm = 1, n - 1
        pl = pm
        DO WHILE (pl > 0)
          IF (a(pl - 1)%lon <= a(pl)%lon) EXIT
          temp = a(pl) ; a(pl) = a( pl - 1); a( pl - 1) = temp
          pl = pl - 1
        END DO
      END DO
      RETURN
    END IF
    pm = n / 2
    IF (n > 7) THEN
      pl = 0
      pn = n - 1
      IF (n > 40) THEN
        d = n / 8
        pl = med3_cell_info_lon(a, pl, pl + d, pl + 2 * d)
        pm = med3_cell_info_lon(a, pm - d, pm, pm + d)
        pn = med3_cell_info_lon(a, pn - 2 * d, pn - d, pn)
      END IF
      pm = med3_cell_info_lon(a, pl, pm, pn)
    END IF
    temp = a(0) ; a(0) = a( pm); a( pm) = temp
    pb = 1
    pa = pb
    pd = n - 1
    pc = pd
    swap_cnt = .FALSE.
    DO WHILE (.TRUE.)
      DO WHILE (pb <= pc)
        IF (a(pb)%lon > a(0)%lon) EXIT
        IF (a(pb)%lon == a(0)%lon) THEN
          swap_cnt = .TRUE.
          temp = a(pa) ; a(pa) = a( pb); a( pb) = temp
          pa = pa + 1
        END IF
        pb = pb + 1
      END DO
      DO WHILE (pb <= pc)
        IF (a(pc)%lon < a(0)%lon) EXIT
        IF (a(pc)%lon == a(0)%lon) THEN
          swap_cnt = .TRUE.
          temp = a(pc) ; a(pc) = a( pd); a( pd) = temp
          pd = pd - 1
        END IF
        pc = pc - 1
      END DO
      IF (pb > pc) EXIT
      temp = a(pb) ; a(pb) = a( pc); a( pc) = temp
      swap_cnt = .TRUE.
      pb = pb + 1
      pc = pc - 1
    END DO
    IF (.NOT. swap_cnt) THEN  ! Switch to insertion sort
      DO pm = 1, n - 1
        pl = pm
        DO WHILE(pl > 0)
          IF (a(pl - 1)%lon <= a(pl)%lon) EXIT
          temp = a(pl) ; a(pl) = a( pl - 1); a( pl - 1) = temp
          pl = pl - 1
        END DO
      END DO
      RETURN
    END IF
    pn =  n
    pdiff = MIN(pa, pb - pa)
    IF (pdiff > 0) &
         CALL cell_info_vec_swap(a(0:pdiff - 1), a(pb - pdiff:pb-1), pdiff)
    pdiff = MIN(pd - pc, pn - pd - 1)
    if (pdiff > 0) &
         CALL cell_info_vec_swap(a(pb:pb + pdiff - 1), a(pn - pdiff:pn - 1), &
         pdiff)
    pdiff = pb - pa
    IF (pdiff > 1) &
         CALL sort_cell_info_by_lon(a, pdiff)
    pdiff = pd - pc
    ! hope the compiler can tail-recurse and save stack space
    IF (pdiff > 1) &
         CALL sort_cell_info_by_lon(a(pn - pdiff:), pdiff)
  END SUBROUTINE sort_cell_info_by_lon

  FUNCTION med3_cell_info_lon(a, i, j, k) RESULT(median)
    TYPE(t_cell_info), INTENT(in) :: a(0:)
    INTEGER, INTENT(in) :: i, j, k
    INTEGER :: median
    IF (a(i)%lon < a(j)%lon) THEN
      IF (a(j)%lon < a(k)%lon) THEN
        median = j
      ELSE
        IF (a(i)%lon < a(k)%lon) THEN
          median = k
        ELSE
          median = i
        END IF
      END IF
    ELSE
      IF (a(j)%lon > a(k)%lon) THEN
        median = j
      ELSE
        IF (a(i)%lon < a(k)%lon) THEN
          median = i
        ELSE
          median = k
        END IF
      END IF
    END IF
  END FUNCTION med3_cell_info_lon

  RECURSIVE SUBROUTINE sort_cell_info_by_cell_number(a, n)
    INTEGER, INTENT(in) :: n
    TYPE(t_cell_info), INTENT(inout) :: a(0:n - 1)
    TYPE(t_cell_info) :: temp

    LOGICAL :: swap_cnt

    INTEGER :: d, pa, pb, pc, pd, pl, pm, pn, pdiff

    IF (n < 7) THEN
      DO pm = 1, n - 1
        pl = pm
        DO WHILE (pl > 0)
          IF (a(pl - 1)%cell_number <= a(pl)%cell_number) EXIT
          temp = a(pl) ; a(pl) = a( pl - 1); a( pl - 1) = temp
          pl = pl - 1
        END DO
      END DO
      RETURN
    END IF
    pm = n / 2
    IF (n > 7) THEN
      pl = 0
      pn = n - 1
      IF (n > 40) THEN
        d = n / 8
        pl = med3_cell_info_cell_number(a, pl, pl + d, pl + 2 * d)
        pm = med3_cell_info_cell_number(a, pm - d, pm, pm + d)
        pn = med3_cell_info_cell_number(a, pn - 2 * d, pn - d, pn)
      END IF
      pm = med3_cell_info_cell_number(a, pl, pm, pn)
    END IF
    temp = a(0) ; a(0) = a( pm); a( pm) = temp
    pb = 1
    pa = pb
    pd = n - 1
    pc = pd
    swap_cnt = .FALSE.
    DO WHILE (.TRUE.)
      DO WHILE (pb <= pc)
        IF (a(pb)%cell_number > a(0)%cell_number) EXIT
        IF (a(pb)%cell_number == a(0)%cell_number) THEN
          swap_cnt = .TRUE.
          temp = a(pa) ; a(pa) = a( pb); a( pb) = temp
          pa = pa + 1
        END IF
        pb = pb + 1
      END DO
      DO WHILE (pb <= pc)
        IF (a(pc)%cell_number < a(0)%cell_number) EXIT
        IF (a(pc)%cell_number == a(0)%cell_number) THEN
          swap_cnt = .TRUE.
          temp = a(pc) ; a(pc) = a( pd); a( pd) = temp
          pd = pd - 1
        END IF
        pc = pc - 1
      END DO
      IF (pb > pc) EXIT
      temp = a(pb) ; a(pb) = a( pc); a( pc) = temp
      swap_cnt = .TRUE.
      pb = pb + 1
      pc = pc - 1
    END DO
    IF (.NOT. swap_cnt) THEN  ! Switch to insertion sort
      DO pm = 1, n - 1
        pl = pm
        DO WHILE(pl > 0)
          IF (a(pl - 1)%cell_number <= a(pl)%cell_number) EXIT
          temp = a(pl) ; a(pl) = a( pl - 1); a( pl - 1) = temp
          pl = pl - 1
        END DO
      END DO
      RETURN
    END IF
    pn =  n
    pdiff = MIN(pa, pb - pa)
    IF (pdiff > 0) &
         CALL cell_info_vec_swap(a(0:pdiff - 1), a(pb - pdiff:pb-1), pdiff)
    pdiff = MIN(pd - pc, pn - pd - 1)
    IF (pdiff > 0) &
         CALL cell_info_vec_swap(a(pb:pb + pdiff - 1), a(pn - pdiff:pn - 1), &
         pdiff)
    pdiff = pb - pa
    IF (pdiff > 1) &
         CALL sort_cell_info_by_cell_number(a, pdiff)
    pdiff = pd - pc
    ! hope the compiler can tail-recurse and save stack space
    IF (pdiff > 1) &
         CALL sort_cell_info_by_cell_number(a(pn - pdiff:), pdiff)
  END SUBROUTINE sort_cell_info_by_cell_number

  FUNCTION med3_cell_info_cell_number(a, i, j, k) RESULT(median)
    TYPE(t_cell_info), INTENT(in) :: a(0:)
    INTEGER, INTENT(in) :: i, j, k
    INTEGER :: median
    IF (a(i)%cell_number < a(j)%cell_number) THEN
      IF (a(j)%cell_number < a(k)%cell_number) THEN
        median = j
      ELSE
        IF (a(i)%cell_number < a(k)%cell_number) THEN
          median = k
        ELSE
          median = i
        END IF
      END IF
    ELSE
      IF (a(j)%cell_number > a(k)%cell_number) THEN
        median = j
      ELSE
        IF (a(i)%cell_number < a(k)%cell_number) THEN
          median = i
        ELSE
          median = k
        END IF
      END IF
    END IF
  END FUNCTION med3_cell_info_cell_number

  FUNCTION t_cell_info_eq(a, b) RESULT(p)
    LOGICAL :: p
    TYPE(t_cell_info), INTENT(in) :: a, b
    p = a%lat == b%lat .AND. a%lon == b%lon &
         .AND. a%cell_number == b%cell_number &
         .AND. a%owner == b%owner
  END FUNCTION t_cell_info_eq

  FUNCTION t_cell_info_ne(a, b) RESULT(p)
    LOGICAL :: p
    TYPE(t_cell_info), INTENT(in) :: a, b
    p = a%lat /= b%lat .OR. a%lon /= b%lon &
         .OR. a%cell_number /= b%cell_number &
         .OR. a%owner /= b%owner
  END FUNCTION t_cell_info_ne

  !> compute start integer of uniform interval partition
  ELEMENTAL FUNCTION uniform_partition_start(set_interval, nparts, &
    &                                        part_idx) RESULT(start)
    INTEGER, INTENT(in) :: nparts
    TYPE(extent), INTENT(in) :: set_interval
    INTEGER, INTENT(in) :: part_idx
    INTEGER :: start, part_offset

    part_offset = INT((INT(extent_size(set_interval), i8) &
         &             * INT(part_idx - 1, i8)) / INT(nparts, i8))
    start = extent_start(set_interval) + part_offset

  END FUNCTION uniform_partition_start

  !> compute nth part of integer set interval
  !!
  !! The interval is divided into roughly same sized sub-intervals
  !! forming a uniform partition.
  !! @param set_interval global domain
  !! @param nparts number of parts to decompose into
  !! @param part_idx number of sub-interval to compute
  !! <tt>SIZE(uniform_partition(i)) == SIZE(uniform_partition(nparts - i + 1))</tt>
  !! @return part range corresponding to part_idx
  ELEMENTAL FUNCTION uniform_partition(set_interval, nparts, part_idx) &
    & RESULT(interval)
    INTEGER, INTENT(in) :: nparts
    TYPE(extent), INTENT(in) :: set_interval
    INTEGER, INTENT(in) :: part_idx
    TYPE(extent) :: interval
    INTEGER :: start_part, start_next_part

    start_part = uniform_partition_start(set_interval, nparts, part_idx)
    start_next_part = uniform_partition_start(set_interval, nparts, &
         part_idx + 1)
    interval = extent(start_part, start_next_part - start_part)
  END FUNCTION uniform_partition

  ELEMENTAL FUNCTION partidx_of_elem_uniform_deco(set_interval, nparts, &
    &                                             elem_idx) RESULT(part_idx)
    TYPE(extent), INTENT(in) :: set_interval
    INTEGER, INTENT(in) :: nparts, elem_idx
    INTEGER :: part_idx

    part_idx = INT((INT(elem_idx - extent_start(set_interval), i8) &
                    * INT(nparts, i8) + INT(nparts, i8) - 1) &
                   / INT(extent_size(set_interval), i8)) + 1

  END FUNCTION partidx_of_elem_uniform_deco

END MODULE mo_decomposition_tools

