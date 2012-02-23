module plot

  use constants
  use error,           only: fatal_error
  use geometry,        only: find_cell, distance_to_boundary, cross_surface, &
                             cross_lattice, cell_contains
  use datatypes,       only: dict_create, dict_add_key, dict_get_key,         &
                              dict_has_key, dict_keys
  use geometry_header, only: Universe, BASE_UNIVERSE
  use global
  use particle_header, only: LocalCoord, deallocate_coord
  use source,          only: initialize_particle

  implicit none

contains

!===============================================================================
! RUN_PLOT
!===============================================================================

  subroutine run_plot()

    !TODO: deal with cell universes

    call find_cell_pointclouds(0)
    !call find_cell_pointclouds(1)
    !call find_cell_pointclouds(2)

    call dump_point_cloud()
    !call dump_cell_quadrics()



  end subroutine run_plot



!===============================================================================
! FIND_CELL_POINTCLOUDS
!===============================================================================
  subroutine find_cell_pointclouds(basis)

    integer, intent(in) :: basis ! sweep direction 0: xy, 1: xz, 2: yz

    integer :: i,n             ! indices
    real(8) :: direction(3)    ! sweep direction
    integer :: outerI, innerI  ! index of xyz for outer and inner sweeps
    integer :: surface_crossed ! surface which particle is on
    integer :: lattice_crossed ! is surface crossing in lattice?
    integer :: last_cell       ! most recent cell particle was in
    integer :: enter_surface   ! entrance surface
    real(8) :: xyz(3)          ! starting coordinates
    real(8) :: last_track_coord    ! bounding x coordinate
    real(8) :: last_inner_coord    ! bounding y coordinate
    real(8) :: last_outer_coord! bounding z coordinate
    real(8) :: d               ! distance to boundary
    real(8) :: distance        ! distance particle travels
    logical :: found_cell      ! found cell which particle is in?
    type(Cell),       pointer :: c    => null()
    type(Universe),   pointer :: univ => null()
    type(LocalCoord), pointer :: coord => null()

    if (basis == 0 ) then
      direction = (/ 1, 0, 0 /)
      innerI = 2 ! y
      outerI = 3 ! z
      ! Determine bounding x and y coordinates for plot
      xyz(1) = plot_origin(1) - plot_width(1) / 2.0
      xyz(2) = plot_origin(2) + plot_width(2) / 2.0
      xyz(3) = plot_origin(3) + plot_width(3) / 2.0
      last_track_coord = plot_origin(1) + plot_width(1) / 2.0
      last_inner_coord = plot_origin(2) - plot_width(2) / 2.0
      last_outer_coord = plot_origin(3) - plot_width(3) / 2.0
    else if (basis == 1) then
      direction = (/ 0, 1, 0 /)
      innerI = 1 ! x
      outerI = 3 ! z
      ! Determine bounding x and y coordinates for plot
      xyz(1) = plot_origin(1) + plot_width(1) / 2.0
      xyz(2) = plot_origin(2) - plot_width(2) / 2.0
      xyz(3) = plot_origin(3) + plot_width(3) / 2.0
      last_inner_coord = plot_origin(1) - plot_width(1) / 2.0
      last_track_coord = plot_origin(2) + plot_width(2) / 2.0
      last_outer_coord = plot_origin(3) - plot_width(3) / 2.0
    else if (basis == 2) then
      direction = (/ 0, 0, 1 /)
      innerI = 2 ! y 
      outerI = 1 ! x
      ! Determine bounding x and y coordinates for plot
      xyz(1) = plot_origin(1) + plot_width(1) / 2.0
      xyz(2) = plot_origin(2) - plot_width(2) / 2.0
      xyz(3) = plot_origin(3) + plot_width(3) / 2.0
      last_outer_coord = plot_origin(1) - plot_width(1) / 2.0
      last_inner_coord = plot_origin(3) - plot_width(3) / 2.0
      last_track_coord = plot_origin(2) + plot_width(2) / 2.0
    else
      return
    end if

    ! allocate and initialize particle
    allocate(p)

    do while(xyz(outerI) > last_outer_coord)
      ! loop over horizontal rays
      do while(xyz(innerI) > last_inner_coord)

         ! initialize the particle and set starting coordinate and direction
         call initialize_particle()

         p % coord % xyz = xyz
         p % coord % uvw = direction

         ! Find cell that particle is currently in
         call find_cell(found_cell)

         if (found_cell) then
            c => cells(p % coord0 % cell)
            call add_pointcloud_point(c)
         end if

         ! =======================================================================
         ! MOVE PARTICLE FORWARD TO NEXT CELL

         if (.not. found_cell) then
            ! Clear any coordinates beyond first level
            call deallocate_coord(p % coord0 % next)
            p % coord => p % coord0

            distance = INFINITY
            univ => universes(BASE_UNIVERSE)
            do i = 1, univ % n_cells
               p % coord0 % xyz = xyz
               p % coord0 % cell = univ % cells(i)

               call distance_to_boundary(d, surface_crossed, lattice_crossed)
               if (d < distance) then
                  ! Check to make sure particle is actually going into this cell
                  ! by moving it slightly forward and seeing if the cell contains
                  ! that coordinate

                  p % coord0 % xyz = p % coord0 % xyz + (d + TINY_BIT) * p % coord0 % uvw

                  c => cells(p % coord0 % cell)
                  if (.not. cell_contains(c)) cycle

                  last_cell = p % coord0 % cell
                  ! Set new distance and retain pointer to this cell
                  distance = d
                  enter_surface = surface_crossed
               end if
            end do

            ! No cell was found on this horizontal ray
            if (distance == INFINITY) then
               p % coord0 % xyz(1) = last_track_coord

               ! Move to next ray start position
               xyz(innerI) = xyz(innerI) - plot_aspect

               cycle
            end if

            ! Write coordinate where next cell begins
            p % coord0 % xyz = xyz + distance * p % coord0 % uvw
            c => cells(last_cell)
            call add_pointcloud_point(c)

            ! Process surface crossing for next cell
            p % coord0 % cell = NONE
            p % surface = -enter_surface
            call cross_surface(NONE)
         end if ! end if (.not. found_cell)

         ! =======================================================================
         ! MOVE PARTICLE ACROSS HORIZONTAL TRACK

         do while (p % alive)
            ! save particle's current cell
            last_cell = p % coord % cell

            c => cells(last_cell)
            call add_pointcloud_point(c)

            ! Calculate distance to next boundary
            call distance_to_boundary(distance, surface_crossed, lattice_crossed)

            ! Advance particle
            coord => p % coord0
            do while (associated(coord))
               coord % xyz = coord % xyz + distance * coord % uvw
               coord => coord % next
            end do

            ! If next boundary crossing is out of range of the plot, only include
            ! the visible portion and move to next horizontal ray
            if (p % coord0 % xyz(1) >= last_track_coord) then
              p % alive = .false.
              p % coord0 % xyz(1) = last_track_coord

              ! If there is no cell beyond this boundary, mark it as cell 0
              if (distance == INFINITY) then
                p % coord % cell = 0
              else
                c => cells(last_cell)
                call add_pointcloud_point(c)
              end if

              cycle
            end if ! end if (p % coord0 % xyz(1) >= last_track_coord)

            c => cells(last_cell)
            call add_pointcloud_point(c)

            p % coord % cell = 0
            if (lattice_crossed /= NONE) then
               p % surface = NONE
               call cross_lattice(lattice_crossed)
            else
               p % surface = surface_crossed
               call cross_surface(last_cell)

               if (surfaces(abs(surface_crossed)) % bc /= BC_TRANSMIT) then
                  c => cells(last_cell)
                  call add_pointcloud_point(c)
                  exit
               end if
            end if

         end do ! (end do while particle is alive)

         ! Move y-coordinate to next ray start position
         xyz(innerI) = xyz(innerI) - plot_aspect
      end do

      if (basis == 0) then
        ! Move z-coordinate to next position and reset x and y
        xyz(1) = plot_origin(1) - plot_width(1) / 2.0
        xyz(2) = plot_origin(2) + plot_width(2) / 2.0
        xyz(3) = xyz(3) - plot_aspect
      else if (basis == 1) then
        ! Move z-coordinate to next position and reset x and y
        xyz(1) = plot_origin(1) + plot_width(1) / 2.0
        xyz(2) = plot_origin(2) - plot_width(2) / 2.0
        xyz(3) = xyz(3) - plot_aspect     
      else if (basis == 2) then
        ! Move x-coordinate to next position and reset y and z
        xyz(2) = plot_origin(2) + plot_width(2) / 2.0
        xyz(3) = plot_origin(3) - plot_width(3) / 2.0
        xyz(1) = xyz(1) - plot_aspect    
      end if

    end do


  end subroutine find_cell_pointclouds

!===============================================================================
! add_pointcloud_point
!===============================================================================
  subroutine add_pointcloud_point(c)

    type(Cell),       pointer :: c

    integer :: i

    call dict_add_key(plot_cells_dict,p % coord % cell,p % coord % cell)
    if (.not. allocated(c%pointcloud)) then
      allocate(c%pointcloud(100000))
      c%n_points = 1
    else
      c%n_points = c%n_points + 1
    endif
    if (c%n_points >=100000 ) return
    c%pointcloud(c%n_points)%xyz = p%coord0%xyz

    call update_cell_limits(c)

  end subroutine add_pointcloud_point


!===============================================================================
! UPDATE_CELL_LIMITS
!===============================================================================
  subroutine update_cell_limits(c)

    type(Cell),       pointer :: c

    integer :: i

    call dict_add_key(plot_cells_dict,p % coord % cell,p % coord % cell)
    if (.not. allocated(c%limits)) then
      allocate(c%limits(6))
      do i=1,3
        c%limits(i) = p%coord0%xyz(i)
        c%limits(i+3) = p%coord0%xyz(i)
      end do
    else
      do i=1,3
        if (p % coord0 % xyz(i) < c % limits(i)) then
          c % limits(i) = p % coord0 % xyz(i)
        end if
        if (p % coord0 % xyz(i) > c % limits(i+3)) then
          c % limits(i+3) = p % coord0 % xyz(i)
        end if
      end do
    endif

  end subroutine update_cell_limits


!===============================================================================
! DUMP_POINT_CLOUD
!===============================================================================
  subroutine dump_point_cloud()

    integer                     :: n, t
    character(MAX_LINE_LEN)     :: path_plot ! unit for binary plot file

    ! Open plot file for binary writing
    path_plot = trim(path_input) // "plot.out"
    open(UNIT=UNIT_PLOT, FILE=path_plot, STATUS="replace", ACCESS="stream")

    do n=1,n_cells
      if (dict_has_key(plot_cells_dict,n)) then
        write(UNIT=UNIT_PLOT) cells(n)%id
        write(UNIT=UNIT_PLOT) cells(n)%n_points
        do t=1,cells(n)%n_points
          write(UNIT=UNIT_PLOT) cells(n)%pointcloud(t)
        end do
      end if
    end do

    ! Close plot file
    close(UNIT=UNIT_PLOT)

  end subroutine dump_point_cloud

!===============================================================================
! DUMP_CELL_QUADRICS
!===============================================================================
  subroutine dump_cell_quadrics()
    integer :: n,s,i             ! indices
    type(Cell),       pointer :: c    => null()
    type(Surface), pointer      :: surf => null()
    
    character(MAX_LINE_LEN)     :: path_plot ! unit for binary plot file

    ! Open plot file for binary writing
    path_plot = trim(path_input) // "plot.out"
    open(UNIT=UNIT_PLOT, FILE=path_plot, STATUS="replace", ACCESS="stream")

    do n=1,n_cells
      if (dict_has_key(plot_cells_dict,n)) then
          c => cells(n)
          write(UNIT=UNIT_PLOT)c%id
          write(UNIT=UNIT_PLOT)c%limits
          write(UNIT=UNIT_PLOT)c%n_surfaces
          do s=1,c%n_surfaces
            if (c%surfaces(s) >= OP_DIFFERENCE) cycle ! TODO: skipped operators for now -- we *will* need this information
            surf => surfaces(abs(c%surfaces(s)))
            write(UNIT=UNIT_PLOT)surf%type
            write(UNIT=UNIT_PLOT)-1*sign(1,c%surfaces(s))
            write(UNIT=UNIT_PLOT)size(surf%coeffs)
            write(UNIT=UNIT_PLOT)surf%coeffs
          end do
      end if
    end do

    ! Close plot file
    close(UNIT=UNIT_PLOT)

  end subroutine dump_cell_quadrics


end module plot
