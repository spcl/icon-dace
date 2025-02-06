!> Data structures and procedures resembling an ICON variable list.
!
!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  Please see the file LICENSE in the root of the source tree for this code.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE vars

  USE mo_utilities,         ONLY : wp, mpi_rank, finish
  USE comin_host_interface, ONLY : t_comin_var_ptr, t_comin_var_descriptor,         &
    &                              comin_var_list_append,                           &
    &                              comin_var_get_descr_list_head,                   &
    &                              comin_request_get_list_head,                     &
    &                              comin_var_update, COMIN_ZAXIS_2D,                &
    &                              COMIN_ZAXIS_3D, t_comin_var_descr_list_item,     &
    &                              t_var_request_list_item,                         &
    &                              comin_metadata_set, comin_metadata_get_or,       &
    &                              COMIN_HGRID_UNSTRUCTURED_CELL,                   &
    &                              COMIN_HGRID_UNSTRUCTURED_EDGE,                   &
    &                              COMIN_METADATA_TYPEID_UNDEFINED,                 &
    &                              COMIN_METADATA_TYPEID_INTEGER,                   &
    &                              COMIN_METADATA_TYPEID_REAL,                      &
    &                              COMIN_METADATA_TYPEID_CHARACTER,                 &
    &                              COMIN_METADATA_TYPEID_LOGICAL
  USE iso_c_binding, ONLY: C_LOC, C_PTR, C_NULL_PTR

  IMPLICIT NONE
  PRIVATE

  ! ICON-specific: separator for varname and time level
  CHARACTER(LEN=3), PARAMETER :: TIMELEVEL_SUFFIX = '.TL'

  PUBLIC :: wp
  PUBLIC :: get_var_name, get_var_timelevel, find_list_element
  PUBLIC :: add_var, add_ref
  PUBLIC :: t_state, t_var, t_var_item
  PUBLIC :: deallocate_vars
  PUBLIC :: varlist
  PUBLIC :: get_ntracer_comin
  PUBLIC :: expose_variables
  PUBLIC :: append_comin_variables
  PUBLIC :: update_exposed_variables

  TYPE :: t_metadata
    LOGICAL            :: restart = .FALSE.
    LOGICAL            :: lcontainer = .FALSE.

    ! index for containers (tracer fields)
    ! ... if this variable is a container, then "ncontained = SIZE(dataarray,var_ref_pos)"
    ! ... if this variable is a reference, then "dataarray => container%dataarray(:,:,:,:,ncontained)"
    ! ... else: "ncontained = 0"
    INTEGER               :: ncontained = 0

    !> vertical and horizontal limiters
    INTEGER              :: itype_vlimit = -1
    INTEGER              :: itype_hlimit = -1

    ! tracer variables correspond to 4D slices of the 5D array, where
    ! the position of the slicing dimension, and the slice index are
    ! stored in the variable's meta-data:
    !        - `info%ncontained`: slice index
    !        - `info%var_ref_pos`: slice position
    INTEGER               :: var_ref_pos = -1
    CHARACTER(LEN=32)     :: units
  END TYPE t_metadata

  TYPE :: t_var
    CHARACTER(LEN=256)    :: name = ''
    INTEGER               :: jg
    TYPE(t_metadata)      :: metadata
    REAL(wp), POINTER     :: dataarray(:,:,:,:,:) => NULL()
  END TYPE t_var

  TYPE :: t_var_item
    TYPE(t_var) :: this
    TYPE(t_var_item), POINTER :: next => NULL()
  END TYPE t_var_item

  TYPE t_state
    REAL(wp), POINTER, DIMENSION(:,:,:,:,:) :: u,v,pres,temp,rho,pres_sfc,vn
  END TYPE t_state

  TYPE(t_var_item), POINTER :: varlist => NULL()

CONTAINS

  !> @return plain variable name (i.e. without TIMELEVEL_SUFFIX)
  !  ... this resembles how it is done in ICON
  CHARACTER(LEN=256) FUNCTION get_var_name(name)
    CHARACTER(LEN=256), INTENT(IN) :: name
    INTEGER :: idx

    idx = INDEX(name,TIMELEVEL_SUFFIX)
    IF (idx .EQ. 0) THEN
      get_var_name = name
    ELSE
      get_var_name = name(1:idx-1)
    END IF
  END FUNCTION get_var_name

  !> @return time level (extracted from time level suffix) or "-1"
  !  ... this resembles how it is done in ICON
  INTEGER FUNCTION get_var_timelevel(vname) RESULT(tl)
    CHARACTER(*), INTENT(IN) :: vname

    tl = INDEX(vname,TIMELEVEL_SUFFIX)
    IF (tl .EQ. 0) THEN
      tl = -1
    ELSE
      tl = ICHAR(vname(tl+3:tl+3)) - ICHAR('0')
    END IF
  END FUNCTION get_var_timelevel

  !> @return ICON variable by list search
  !  ... this resembles how it is done in ICON
  FUNCTION find_list_element(vname, tl)  RESULT(element)
    TYPE(t_var_item), POINTER :: element
    CHARACTER(*), INTENT(IN) :: vname  !< variable name (without time level suffix)
    INTEGER,      INTENT(IN) :: tl     !< time level

    CHARACTER(LEN=256) :: query_name
    TYPE(t_var_item), POINTER :: vn

    element => NULL()
    WRITE(query_name, '(a,a,i0)') TRIM(vname),TIMELEVEL_SUFFIX,tl
    vn => varlist
    DO WHILE (ASSOCIATED(vn))
      IF (TRIM(vn%this%name) == TRIM(query_name)) THEN
        element => vn; RETURN
      END IF
      vn => vn%next
    END DO
  END FUNCTION find_list_element

  ! The following code simulates ICON's "add_var"
  !
  SUBROUTINE add_var(name, tl, lcontainer, dimshape, units, vc)
    CHARACTER(LEN=*),          INTENT(IN)  :: name
    INTEGER,                   INTENT(IN)  :: tl    !< time level
    LOGICAL,                   INTENT(IN)  :: lcontainer
    INTEGER,                   INTENT(IN)  :: dimshape(5)
    CHARACTER(LEN=*),          INTENT(IN)  :: units
    TYPE(t_var_item), POINTER, INTENT(OUT), OPTIONAL :: vc
    !
    TYPE(t_var_item), POINTER :: vc_new, vn, element

    element => find_list_element(TRIM(name), tl)
    IF (ASSOCIATED(element)) THEN
      CALL finish("vars::add_var", "Duplicate variable exists!")
    END IF

    ALLOCATE(vc_new)
    WRITE(vc_new%this%name, '(2a,i0)') TRIM(name),TIMELEVEL_SUFFIX,tl
    ALLOCATE (vc_new%this%dataarray(dimshape(1), dimshape(2), dimshape(3), &
      &                             dimshape(4), dimshape(5)))
    vc_new%this%metadata%lcontainer = lcontainer
    vc_new%this%metadata%ncontained = 0
    vc_new%this%metadata%units = units
    vc_new%next => NULL()
    IF (.NOT. ASSOCIATED(varlist)) THEN
      varlist => vc_new
    ELSE
      vn => varlist
      DO WHILE (ASSOCIATED(vn%next))
        vn => vn%next
      END DO
      vn%next => vc_new
    END IF
    IF (PRESENT(vc))  vc => vc_new

  END SUBROUTINE add_var

  ! The following code simulates ICON's "add_ref"
  !
  SUBROUTINE add_ref(targetname, refname, tl, hl, vl, vc, is_3d_var, units, opt_var_ref_pos)
    CHARACTER(LEN=*),          INTENT(IN)  :: targetname, refname
    INTEGER,                   INTENT(IN)  :: tl    !< time level
    INTEGER,         OPTIONAL, INTENT(IN)  :: hl
    INTEGER,         OPTIONAL, INTENT(IN)  :: vl
    LOGICAL,                   INTENT(IN)  :: is_3d_var
    CHARACTER(LEN=*),          INTENT(IN)  :: units
    TYPE(t_var_item), POINTER, INTENT(OUT), OPTIONAL :: vc
    INTEGER,                   INTENT(IN),  OPTIONAL :: opt_var_ref_pos
    !
    TYPE(t_var_item), POINTER :: vc_new, vn, tgt_element, element
    INTEGER :: icontained, ndims

    element => find_list_element(TRIM(refname), tl)
    IF (ASSOCIATED(element)) THEN
      CALL finish("vars::add_ref", "Duplicate variable exists!")
    END IF

    tgt_element => find_list_element(TRIM(targetname), tl)
    IF (.NOT. ASSOCIATED(tgt_element)) THEN
      CALL finish("vars::add_ref", "Target variable '"//TRIM(targetname)//"' not found!")
    END IF
    IF (.NOT. tgt_element%this%metadata%lcontainer) THEN
      CALL finish("vars::add_ref", "Target variable is not a container!")
    END IF
    icontained = tgt_element%this%metadata%ncontained + 1
    tgt_element%this%metadata%ncontained = icontained

    ALLOCATE(vc_new)
    WRITE(vc_new%this%name, '(2a,i0)') TRIM(refname),TIMELEVEL_SUFFIX,tl
    vc_new%this%metadata%lcontainer = .FALSE.
    vc_new%this%metadata%ncontained = icontained
    vc_new%this%metadata%units = units
    IF (PRESENT(hl)) vc_new%this%metadata%itype_hlimit = hl
    IF (PRESENT(vl)) vc_new%this%metadata%itype_vlimit = vl
    vc_new%this%dataarray => tgt_element%this%dataarray

    IF (.NOT. ASSOCIATED(varlist)) THEN
      varlist => vc_new
      vn => varlist
    ELSE
      vn => varlist
      DO WHILE (ASSOCIATED(vn%next))
        vn => vn%next
      END DO
      vn%next => vc_new
    END IF
    IF (PRESENT(vc))  vc => vc_new

    IF (is_3d_var) THEN
      ndims = 3
    ELSE
      ndims = 2
    END IF

    ! `info%var_ref_pos`: slice position
    IF (.NOT. PRESENT(opt_var_ref_pos)) THEN
      vc_new%this%metadata%var_ref_pos      = ndims + 1
      IF (tgt_element%this%metadata%var_ref_pos == -1) THEN
        tgt_element%this%metadata%var_ref_pos = ndims + 1
      END IF
    ELSE
      vc_new%this%metadata%var_ref_pos      = opt_var_ref_pos
      IF (tgt_element%this%metadata%var_ref_pos == -1) THEN
        tgt_element%this%metadata%var_ref_pos = opt_var_ref_pos
      END IF
    END IF
  END SUBROUTINE add_ref

  SUBROUTINE deallocate_vars()
    TYPE(t_var_item), POINTER :: vn

    vn => varlist
    DO WHILE (ASSOCIATED(vn))
      varlist => vn
      vn => vn%next
      ! clean up only true variables (not references):
      IF (varlist%this%metadata%lcontainer .OR. &
        & (varlist%this%metadata%ncontained==0)) THEN
        DEALLOCATE(varlist%this%dataarray)
      END IF
      DEALLOCATE(varlist)
    END DO
  END SUBROUTINE deallocate_vars

  !> Expose mock-up-ICON's variables to the ComIn infrastructure.
  SUBROUTINE expose_variables(nlev)
    INTEGER, INTENT(IN) :: nlev
    !
    LOGICAL :: is_3d_field, tracer
    INTEGER :: zaxis_id, pos_jcjkjb(3)
    TYPE(t_var_item), POINTER :: vn
    TYPE(t_comin_var_descriptor) :: var_descriptor
    TYPE(t_comin_var_ptr),  POINTER :: comin_var_ptr
    REAL(wp), POINTER :: dataptr(:,:,:,:,:)

    IF (mpi_rank == 0)  WRITE (0,*) "     Associate ICON varlist and ComIn"
    vn => varlist
    DO WHILE (ASSOCIATED(vn))
      ! for time-level dependent variables: register only once
      IF (get_var_timelevel(vn%this%name) == 1) THEN

        ! Provide data pointer, together with index positions (for
        ! correct interpretation).
        !
        ! for this mockup, we assume that this is an atmospheric (3D)
        ! field if the second array dimension equals `nlev`:
        is_3d_field = SIZE(vn%this%dataarray,2) == nlev
        !
        ! In ICON, arrays have the implicit ordering: `jc`, `jk`, `jb`.
        !
        ! - in the case of 2D arrays, `jk` is omitted
        ! - the dimensions start with position 1 (i.e. usually:
        !   dimension position 1 = `jc`, position 2 = `jk`, position 3 = `jb`)
        ! - tracer variables correspond to 4D slices of the 5D array,
        !   where the position of the slicing dimension, and the slice
        !   index are stored in the variable's meta-data (`info%var_ref_pos`)

        IF (.NOT. is_3d_field) THEN
          pos_jcjkjb =  [1, -1, 2]
        ELSE
          pos_jcjkjb =  [1, 2, 3]

          IF (vn%this%metadata%ncontained > 0) THEN
            IF (vn%this%metadata%var_ref_pos <= pos_jcjkjb(1)) THEN
              pos_jcjkjb(1:3) = pos_jcjkjb(1:3) + 1
            ELSEIF (vn%this%metadata%var_ref_pos <= pos_jcjkjb(2)) THEN
              pos_jcjkjb(2:3) = pos_jcjkjb(2:3) + 1
            ELSEIF (vn%this%metadata%var_ref_pos <= pos_jcjkjb(3)) THEN
              pos_jcjkjb(3:3) = pos_jcjkjb(3:3) + 1
            END IF
          END IF
        END IF

        ! for references into container variables, the array slice
        ! information is stored implicity, by the `ncontained` index
        ! (similar to ICON).
        var_descriptor = t_comin_var_descriptor(name=get_var_name(vn%this%name), id=1)

        ALLOCATE(comin_var_ptr)
        dataptr => get_var_ptr(vn%this)
        comin_var_ptr = t_comin_var_ptr(ptr=dataptr,                &
             &                      pos_jc = pos_jcjkjb(1),         &
             &                      pos_jk = pos_jcjkjb(2),         &
             &                      pos_jb = pos_jcjkjb(3),         &
             &                      pos_jn = vn%this%metadata%var_ref_pos,    &
             &                      lcontainer = vn%this%metadata%lcontainer, &
             &                      ncontained=vn%this%metadata%ncontained,   &
             &                      descriptor = var_descriptor)

        IF (is_3d_field) THEN
          zaxis_id = COMIN_ZAXIS_3D
        ELSE
          zaxis_id = COMIN_ZAXIS_2D
        END IF

        tracer = (.NOT. vn%this%metadata%lcontainer) .AND. (vn%this%metadata%ncontained > 0)

        CALL comin_var_list_append( &
          &       p          = comin_var_ptr)
        CALL comin_metadata_set( var_descriptor, "restart", vn%this%metadata%restart)
        CALL comin_metadata_set( var_descriptor, "tracer", tracer)
        CALL comin_metadata_set( var_descriptor, "zaxis_id", zaxis_id)
        CALL comin_metadata_set( var_descriptor, "units", vn%this%metadata%units)
        CALL comin_metadata_set( var_descriptor, "multi_timelevel", .FALSE.)

        IF (TRIM(get_var_name(vn%this%name)) == "vn") THEN
          CALL comin_metadata_set( var_descriptor, "hgrid_id", COMIN_HGRID_UNSTRUCTURED_EDGE)
        ELSE
          CALL comin_metadata_set( var_descriptor, "hgrid_id", COMIN_HGRID_UNSTRUCTURED_CELL )
        END IF

      END IF
      vn => vn%next
    END DO
  END SUBROUTINE expose_variables

  !> loop over the list of requested variables and count the number of
  !  additional tracer variables:
  INTEGER FUNCTION get_ntracer_comin() RESULT(ntracer)
    !
    TYPE(t_var_request_list_item), POINTER :: ptr
    LOGICAL :: tracer

    ntracer = 0
    ptr => comin_request_get_list_head()
    DO WHILE (ASSOCIATED(ptr))
      ! Note: only single-time-level variables are added
      ASSOCIATE (var_list_request_element => ptr%item_value)
        CALL comin_metadata_get_or(var_list_request_element%metadata, "tracer", tracer, .FALSE.)
        IF (tracer) THEN
          ntracer = ntracer + 1
        END IF
      END ASSOCIATE
      ptr => ptr%next()
    END DO
  END FUNCTION get_ntracer_comin

  !> loop over the total list of additional requested variables and
  !  perform `add_var` / `add_ref` operations needed.
  SUBROUTINE append_comin_variables(nproma, nlev, nblks_c)
    INTEGER, INTENT(IN) :: nproma, nlev, nblks_c
    !
    INTEGER :: dimshape(5)
    TYPE(t_var_request_list_item), POINTER :: ptr
    LOGICAL :: tracer
    INTEGER :: zaxis_id
    CHARACTER(LEN=:), ALLOCATABLE :: units

    ptr => comin_request_get_list_head()
    DO WHILE (ASSOCIATED(ptr))
      ! Note: only single-time-level variables are added
      ASSOCIATE (var_list_request_element => ptr%item_value)
        CALL comin_metadata_get_or(var_list_request_element%metadata, "zaxis_id", zaxis_id, COMIN_ZAXIS_3D)
        CALL comin_metadata_get_or(var_list_request_element%metadata, "tracer", tracer, .FALSE.)
        CALL comin_metadata_get_or(var_list_request_element%metadata, "units", units, "")
        IF (tracer) THEN
          CALL add_ref("tracer", var_list_request_element%descriptor%name, tl=1,  &
               &    is_3d_var = (zaxis_id == COMIN_ZAXIS_3D), &
               &    units=units)
          ! ComIn is informed about the tracer index through expose_variables
        ELSE
          IF (zaxis_id == COMIN_ZAXIS_3D) THEN
            dimshape=[nproma,nlev,nblks_c,1,1]
            CALL add_var(var_list_request_element%descriptor%name, tl=1, &
              &          lcontainer=.FALSE.,                             &
              &          dimshape=dimshape,                              &
              &          units=units)
          ELSE
            dimshape=[nproma,nblks_c,1,1,1]
            CALL add_var(var_list_request_element%descriptor%name, tl=1, &
              &          lcontainer=.FALSE.,                             &
              &          dimshape=dimshape,                              &
              &          units=units)
          END IF
        END IF
        IF (mpi_rank == 0)  WRITE (0,*) "        add variable '", TRIM(var_list_request_element%descriptor%name), &
          &         "', requested by third party plugins."
      END ASSOCIATE
      ptr => ptr%next()
    END DO
  END SUBROUTINE append_comin_variables

  !> Update of exposed variables (ComIn)
  SUBROUTINE update_exposed_variables(nnow)
    INTEGER, INTENT(IN) :: nnow
    !
    TYPE(t_comin_var_descriptor) :: var_descr
    TYPE(t_comin_var_descr_list_item), POINTER :: ptr
    TYPE(t_var_item), POINTER :: vn

    ptr => comin_var_get_descr_list_head()
    DO WHILE (ASSOCIATED(ptr))
      ASSOCIATE (descriptor => ptr%item_value)
        vn => find_list_element(descriptor%name, nnow)
        IF (ASSOCIATED(vn)) THEN
          var_descr%name = TRIM(ADJUSTL(descriptor%name))
          var_descr%id   = 1
          CALL comin_var_update(var_descr, get_var_ptr(vn%this), C_NULL_PTR)
        ELSE
          IF (mpi_rank == 0)  WRITE (0,*) "WARNING: Variable not found: name, nnow = ", &
               &               TRIM(descriptor%name), nnow
        END IF
      END ASSOCIATE
      ptr => ptr%next()
    END DO
  END SUBROUTINE update_exposed_variables

  FUNCTION get_var_ptr(v) RESULT(dataptr)
    TYPE(t_var), INTENT(IN)          :: v
    REAL(wp),    POINTER             :: dataptr(:,:,:,:,:)
    !
    INTEGER :: ncontained, slice_pos = -1
    ncontained = v%metadata%ncontained
    IF (.NOT. v%metadata%lcontainer) THEN
      slice_pos = v%metadata%var_ref_pos
    END IF
    SELECT CASE (slice_pos)
    CASE (-1)
      dataptr => v%dataarray
    CASE (1)
      dataptr => v%dataarray(ncontained:ncontained,:,:,:,:)
    CASE (2)
      dataptr => v%dataarray(:,ncontained:ncontained,:,:,:)
    CASE (3)
      dataptr => v%dataarray(:,:,ncontained:ncontained,:,:)
    CASE (4)
      dataptr => v%dataarray(:,:,:,ncontained:ncontained,:)
    CASE (5)
      dataptr => v%dataarray(:,:,:,:,ncontained:ncontained)
    CASE DEFAULT
      CALL finish("vars::get_var_ptr", "invaild value for var_ref_pos.")
    END SELECT
  END FUNCTION get_var_ptr

END MODULE vars
