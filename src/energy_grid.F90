module energy_grid

  use ace_header,       only: Nuclide
  use constants,        only: MAX_LINE_LEN
  use global
  use list_header,      only: ListReal
  use output,           only: write_message
  use search,           only: binary_search

  implicit none

  ! number of unionized global energy grid hash bins
  integer :: n_hash_bins

  ! indices on unionized global grid for hash bin boundaries
  integer, allocatable :: hash_indices(:)

  ! maximum energy on unionized global energy grid
  real(8) :: E_grid_max

  ! global hash table bin width in lethargy
  real(8) :: du_hash

  ! global hash table bin width in energy
  real(8) :: dE_hash

  implicit none

contains

!===============================================================================
! HASH_TABLE creates an energy hash table to accelerate searches on cross
! section energy grids (globally unionized or nuclide-specific)
!===============================================================================

  subroutine hash_table()

    type(Nuclide), pointer :: nuc => null()
    real(8), allocatable   :: e_loc(:)  ! energy bounds on the hash bins
    real(8)                :: E_min_loc ! minimum energy of the grid
    integer                :: i ! hash bin loop index
    integer                :: j ! nuclide loop index
    integer                :: k ! material loop index

    select case(grid_method)
    case(GRID_UNION)
      allocate(e_loc(n_hash_bins + 1))
      allocate(hash_indices(n_hash_bins + 1))
      E_min_loc  = e_grid(1)
      E_grid_max = e_grid(n_grid)
      du_hash = log(E_grid_max / E_min_loc) / dble(n_hash_bins)
      dE_hash = (E_grid_max - E_min_loc) / dble(n_hash_bins)
      select case(hash_spacing)
      case(LETHARGY)
        e_loc = E_grid_max / exp((/(dble(i) * du_hash, &
          & i = n_hash_bins, 0, -1)/))
      case(ENERGY)
        e_loc = E_min_loc + (/(dble(i) * dE_hash, i = 0, n_hash_bins)/)
      end select
      hash_indices(1) = 1
      hash_indices(n_hash_bins + 1) = n_grid - 1
      do i = 2, n_hash_bins
        hash_indices(i) = binary_search(e_grid, n_grid, e_loc(i))
      end do
      deallocate(e_loc)

    case(GRID_NUCLIDE)
      do j = 1, n_nuclides_total
        nuc => nuclides(j)
        nuc % n_hash_bins = n_hash_bins
        allocate(e_loc(nuc % n_hash_bins + 1))
        allocate(nuc % hash_indices(nuc % n_hash_bins + 1))
        E_min_loc = nuc % energy(1)
        nuc % E_grid_max = nuc % energy(nuc % n_grid)
        nuc % du_hash = log(nuc % E_grid_max / E_min_loc) &
          & / dble(nuc % n_hash_bins)
        nuc % dE_hash = (nuc % E_grid_max - E_min_loc) &
          & / dble(nuc % n_hash_bins)
        select case(hash_spacing)
        case(LETHARGY)
          e_loc = nuc % E_grid_max / exp((/(dble(i) * nuc % du_hash, &
            & i = nuc % n_hash_bins, 0, -1)/))
        case(ENERGY)
          e_loc = E_min_loc + (/(dble(i) * nuc % dE_hash, &
            & i = 0, nuc % n_hash_bins)/)
        end select
        nuc % hash_indices(1) = 1
        nuc % hash_indices(nuc % n_hash_bins + 1) = nuc % n_grid - 1
        do i = 2, nuc % n_hash_bins
          nuc % hash_indices(i) = &
            & binary_search(nuc % energy, nuc % n_grid, e_loc(i))
        end do
        deallocate(e_loc)
      end do
    end select
  end subroutine hash_table

!===============================================================================
! UNIONIZED_GRID creates a single unionized energy grid combined from each
! nuclide of each material. Right now, the grid for each nuclide is added into a
! linked list one at a time with an effective insertion sort. Could be done with
! a hash for all energy points and then a quicksort at the end (what hash
! function to use?)
!===============================================================================

  subroutine unionized_grid()

    integer :: i ! index in nuclides array
    type(ListReal), pointer :: list => null()
    type(Nuclide),  pointer :: nuc => null()

    call write_message("Creating unionized energy grid...", 5)

    ! Add grid points for each nuclide in the problem
    do i = 1, n_nuclides_total
      nuc => nuclides(i)
      call add_grid_points(list, nuc % energy)
    end do

    ! Set size of unionized energy grid 
    n_grid = list % size() 

    ! create allocated array from linked list
    allocate(e_grid(n_grid))
    do i = 1, n_grid
      e_grid(i) = list % get_item(i)
    end do

    ! delete linked list and dictionary
    call list % clear()
    deallocate(list)

    ! Set pointers to unionized energy grid for each nuclide
    call grid_pointers()

  end subroutine unionized_grid

!===============================================================================
! ADD_GRID_POINTS adds energy points from the 'energy' array into a linked list
! of points already stored from previous arrays.
!===============================================================================

  subroutine add_grid_points(list, energy)

    type(ListReal), pointer :: list
    real(8), intent(in) :: energy(:)

    integer :: i       ! index in energy array
    integer :: n       ! size of energy array
    integer :: current ! current index 
    real(8) :: E       ! actual energy value

    i = 1
    n = size(energy)

    ! If the original list is empty, we need to allocate the first element and
    ! store first energy point
    if (.not. associated(list)) then
      allocate(list)
      do i = 1, n
        call list % append(energy(i))
      end do
      return
    end if

    ! Set current index to beginning of the list 
    current = 1

    do while (i <= n)
      E = energy(i)

      ! If we've reached the end of the grid energy list, add the remaining
      ! energy points to the end
      if (current > list % size()) then
        ! Finish remaining energies
        do while (i <= n)
          call list % append(energy(i))
          i = i + 1
        end do
        exit
      end if

      if (E < list % get_item(current)) then

        ! Insert new energy in this position
        call list % insert(current, E)

        ! Advance index in linked list and in new energy grid
        i = i + 1
        current = current + 1

      elseif (E == list % get_item(current)) then
        ! Found the exact same energy, no need to store duplicates so just
        ! skip and move to next index
        i = i + 1
        current = current + 1
      else
        current = current + 1
      end if

    end do

  end subroutine add_grid_points

!===============================================================================
! GRID_POINTERS creates an array of pointers (ints) for each nuclide to link
! each point on the nuclide energy grid to one on the unionized energy grid
!===============================================================================

  subroutine grid_pointers()

    integer :: i            ! loop index for nuclides
    integer :: j            ! loop index for nuclide energy grid
    integer :: index_e      ! index on union energy grid
    real(8) :: union_energy ! energy on union grid
    real(8) :: energy       ! energy on nuclide grid
    type(Nuclide), pointer :: nuc => null()

    do i = 1, n_nuclides_total
      nuc => nuclides(i)
      allocate(nuc % grid_index(n_grid))

      index_e = 1
      energy = nuc % energy(index_e)

      do j = 1, n_grid
        union_energy = e_grid(j)
        if (union_energy >= energy .and. index_e < nuc % n_grid) then
          index_e = index_e + 1
          energy = nuc % energy(index_e)
        end if
        nuc % grid_index(j) = index_e - 1
      end do
    end do

  end subroutine grid_pointers

end module energy_grid
