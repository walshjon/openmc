module plot

  use constants
  use error,           only: fatal_error
  use geometry,        only: find_cell, distance_to_boundary, cross_surface, &
                             cross_lattice, cell_contains, on_cell_surf
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
    call find_cell_pointclouds(1)
    call find_cell_pointclouds(2)

    call dump_point_cloud()
    !call dump_cell_quadrics()

  end subroutine run_plot



!===============================================================================
! FIND_CELL_POINTCLOUDS
!===============================================================================
  subroutine find_cell_pointclouds(basis)

    integer, intent(in) :: basis ! sweep direction 0: xy, 1: xz, 2: yz

    real(8) :: direction(3)    ! sweep direction
    integer :: outerI, innerI, trackI  ! index of xyz for outer, inner, and track sweeps
    integer :: enter_surface   ! entrance surface
    integer :: surface_crossed ! surface which particle is on
    integer :: lattice_crossed ! is surface crossing in lattice?
    integer :: last_cell       ! most recent cell particle was in
    real(8) :: xyz(3)          ! starting coordinates
    real(8) :: last_track_coord    ! bounding x coordinate
    real(8) :: last_inner_coord    ! bounding y coordinate
    real(8) :: last_outer_coord! bounding z coordinate
    real(8) :: distance        ! distance particle travels
    logical :: found_cell      ! found cell which particle is in?
    type(Cell),       pointer :: c    => null()
    type(LocalCoord), pointer :: coord => null()



    if (basis == 0 ) then
      direction = (/ 1, 0, 0 /)
      trackI = 1 ! x
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
      trackI = 2 ! y
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
      trackI = 3 ! z
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

    do while(xyz(outerI) >= last_outer_coord)
      ! loop over horizontal rays
      do while(xyz(innerI) >= last_inner_coord)

         write(*,*)"starting ray at ",xyz

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

            call track_through_void(distance, enter_surface, xyz)

            ! No cell was found on this horizontal ray
            if (distance == INFINITY) then

               p % coord0 % xyz(trackI) = last_track_coord

               ! Move to next ray start position
               xyz(innerI) = xyz(innerI) - plot_aspect

               cycle
            end if

            ! Advance particle to surface
            p % coord0 % xyz = xyz + distance * p % coord0 % uvw

            ! Process surface crossing for next cell
            p % coord0 % cell = NONE
            p % surface = -enter_surface
            call cross_surface(NONE)

         end if

         ! =======================================================================
         ! MOVE PARTICLE ACROSS HORIZONTAL TRACK

         do while (p % alive)
            write(*,*)"at:",p % coord0 % xyz
            write(*,*)"in:",p % coord % cell

            if (p % coord % cell == 0) then
              message = "Voids inside the plotting area, or bad boundary conditions."
              call fatal_error()

              !call track_through_void(distance, enter_surface, p % coord % xyz)

              ! No cell was found on this horizontal ray
              !if (distance == INFINITY) then
              !   exit
              !end if

              ! Write coordinate where next cell begins
              !write(*,*)"moving: ",p % coord % xyz, distance
              !p % coord0 % xyz = p % coord % xyz + distance * p % coord % uvw
              !write(*,*)"   about to add ",p % coord0 % xyz
              !c => cells(p % coord0 % cell)
              !call add_pointcloud_point(c)

              ! Process surface crossing for next cell
              !p % coord0 % cell = NONE
              !p % surface = -enter_surface
              !call cross_surface(NONE)
              !write(*,*)"here"
              !cycle

            else

              c => cells(p % coord % cell)
              call add_pointcloud_point(c)

              ! save particle's current cell
              last_cell = p % coord % cell
              call distance_to_boundary(distance, surface_crossed, lattice_crossed)
            end if

            ! Advance particle
            coord => p % coord0
            do while (associated(coord))
               coord % xyz = coord % xyz + distance * coord % uvw
               coord => coord % next
            end do

            ! If next boundary crossing is out of range of the plot, only include
            ! the visible portion and move to next horizontal ray
            if (p % coord0 % xyz(trackI) >= last_track_coord) then
              p % alive = .false.
              p % coord0 % xyz(trackI) = last_track_coord

              !c => cells(last_cell)
              !call add_pointcloud_point(c)

              exit

            end if

            c => cells(p % coord % cell)
            call add_pointcloud_point(c)

            p % coord % cell = 0
            if (lattice_crossed /= NONE) then
               p % surface = NONE
               call cross_lattice(lattice_crossed)
            else
               p % surface = surface_crossed
               call cross_surface(last_cell)
                write(*,*)"SURFACE",surface_crossed
               if (surfaces(abs(surface_crossed)) % bc /= BC_TRANSMIT) then
                  exit
                  !c => cells(last_cell)
                  !call add_pointcloud_point(c)

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
! track_through_void
!===============================================================================
  subroutine track_through_void(distance, enter_surface, xyz)

    real(8), intent(inout) :: distance        ! distance particle travels
    integer, intent(inout) :: enter_surface   ! entrance surface
    real(8), intent(in)    :: xyz(3)          ! starting coordinates

    integer :: i               ! indices
    integer :: surface_crossed ! surface which particle is on
    integer :: lattice_crossed ! is surface crossing in lattice?
    integer :: last_cell       ! most recent cell particle was in
    real(8) :: d               ! distance to boundary

    type(Cell),       pointer :: c    => null()
    type(Universe),   pointer :: univ => null()

    write(*,*)"searching from void p:",xyz
    write(*,*)"along ",p % coord0 % uvw

    ! Clear any coordinates beyond first level
    call deallocate_coord(p % coord0 % next)
    p % coord => p % coord0

    last_cell = 0
    enter_surface = NONE
    distance = INFINITY
    univ => universes(BASE_UNIVERSE)
    do i = 1, univ % n_cells
       p % coord0 % xyz = xyz
       p % coord0 % cell = univ % cells(i)

       call distance_to_boundary(d, surface_crossed, lattice_crossed)
       if (d < distance) then
          ! Check to make sure particle is actually going into this cell
          c => cells(p % coord0 % cell)
          p % coord0 % xyz = p % coord0 % xyz + d * p % coord0 % uvw

          if (.not. on_cell_surf(c,p % coord0 % xyz)) cycle

          last_cell = p % coord0 % cell
          ! Set new distance and retain pointer to this cell
          distance = d
          enter_surface = surface_crossed
       end if
    end do

    write(*,*)"found dist",distance

    p % coord0 % cell = last_cell
    p % coord0 % xyz = xyz
    p % coord % xyz = xyz



  end subroutine track_through_void

!===============================================================================
! add_pointcloud_point
!===============================================================================
  subroutine add_pointcloud_point(c)

    type(Cell),       pointer :: c

    call dict_add_key(plot_cells_dict,c%id,c%id)
    if (.not. allocated(c%pointcloud)) then
      allocate(c%pointcloud(100000))
    end if
    c%n_points = c%n_points + 1

    if (c%n_points > 100000 ) return

    write(*,*)c%id
    write(*,*)"adding ",c%n_points,p%coord0%xyz
    write(*,*)""

    c%pointcloud(c%n_points)%xyz = p%coord0%xyz

    call update_cell_limits(c)

  end subroutine add_pointcloud_point


!===============================================================================
! UPDATE_CELL_LIMITS
!===============================================================================
  subroutine update_cell_limits(c)

    type(Cell),       pointer :: c

    integer :: i

    call dict_add_key(plot_cells_dict,c%id,c%id)
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
      write(*,*)"checking cell",cells(n)%id
      if (dict_has_key(plot_cells_dict,cells(n)%id)) then
        write(*,*)"    dumping cell ",cells(n)%id
        write(*,*)"        n_points",cells(n)%n_points
        write(UNIT=UNIT_PLOT) cells(n)%id
        write(UNIT=UNIT_PLOT) cells(n)%n_points
        do t=1,cells(n)%n_points
          !write(*,*)"            p:",cells(n)%pointcloud(t)
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
    integer :: n,s            ! indices
    type(Cell),       pointer :: c    => null()
    type(Surface), pointer      :: surf => null()
    
    character(MAX_LINE_LEN)     :: path_plot ! unit for binary plot file

    ! Open plot file for binary writing
    path_plot = trim(path_input) // "plot.out"
    open(UNIT=UNIT_PLOT, FILE=path_plot, STATUS="replace", ACCESS="stream")

    do n=1,n_cells
      c => cells(n)
      if (dict_has_key(plot_cells_dict,c%id)) then
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
