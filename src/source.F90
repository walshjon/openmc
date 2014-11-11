module source

  use bank_header,      only: Bank
  use constants
  use error,            only: fatal_error
  use geometry,         only: find_cell
  use geometry_header,  only: BASE_UNIVERSE
  use global
  use math,             only: maxwell_spectrum, watt_spectrum
  use output,           only: write_message
  use output_interface, only: BinaryOutput
  use particle_header,  only: Particle
  use random_lcg,       only: prn, set_particle_seed, prn_set_stream
  use string,           only: to_str

#ifdef MPI
  use mpi
#endif

  implicit none

contains

!===============================================================================
! INITIALIZE_SOURCE initializes particles in the source bank
!===============================================================================

  subroutine initialize_source()

    character(MAX_FILE_LEN) :: filename
    integer(8) :: i          ! loop index over bank sites
    integer(8) :: id         ! particle id
    integer(4) :: itmp       ! temporary integer
    type(Bank), pointer :: src => null() ! source bank site
    type(BinaryOutput) :: sp ! statepoint/source binary file

    call write_message("Initializing source particles...", 6)

    if (path_source /= '') then
      ! Read the source from a binary file instead of sampling from some
      ! assumed source distribution

      call write_message('Reading source file from ' // trim(path_source) &
           &// '...', 6)

      ! Open the binary file
      call sp % file_open(path_source, 'r', serial = .false.)

      ! Read the file type
      call sp % read_data(itmp, "filetype")

      ! Check to make sure this is a source file
      if (itmp /= FILETYPE_SOURCE) then
        call fatal_error("Specified starting source file not a source file &
             &type.")
      end if

      ! Read in the source bank
      call sp % read_source_bank()

      ! Close file
      call sp % file_close()

    else
      ! Generation source sites from specified distribution in user input
      do i = 1, work
        ! Get pointer to source bank site
        src => source_bank(i)

        ! initialize random number seed
        id = work_index(rank) + i
        call set_particle_seed(id)

        ! sample external source distribution
        call sample_external_source(src)
      end do
    end if

    ! Write out initial source
    if (write_initial_source) then
      call write_message('Writing out initial source guess...', 1)
#ifdef HDF5
      filename = trim(path_output) // 'initial_source.h5'
#else
      filename = trim(path_output) // 'initial_source.binary'
#endif
      call sp % file_create(filename, serial = .false.)
      call sp % write_source_bank()
      call sp % file_close()
    end if

  end subroutine initialize_source

!===============================================================================
! SAMPLE_EXTERNAL_SOURCE
!===============================================================================

  subroutine sample_external_source(site)

    type(Bank), pointer :: site ! source site

    integer :: i          ! dummy loop index
    integer :: j          ! Cartesian axis index
    real(8) :: r          ! sampled coordinates
    real(8) :: phi        ! azimuthal angle
    real(8) :: mu         ! cosine of polar angle
    real(8) :: p_min      ! minimum coordinates of source
    real(8) :: p_max      ! maximum coordinates of source
    real(8) :: a          ! Arbitrary parameter 'a'
    real(8) :: b          ! Arbitrary parameter 'b'
    logical :: found      ! Does the source particle exist within geometry?
    type(Particle) :: p   ! Temporary particle for using find_cell
    integer, save :: num_resamples = 0 ! Number of resamples encountered

    ! Set weight to one by default
    site % wgt = ONE

    ! Set the random number generator to the source stream.
    call prn_set_stream(STREAM_SOURCE)

    ! Sample position
    select case (external_source % type_space)
    case (SRC_SPACE_BOX)
      ! Set particle defaults
      call p % initialize()
      ! Repeat sampling source location until a good site has been found
      found = .false.
      do while (.not.found)

        ! loop over axes
        do j = 1, 3
          p_min = external_source % params_space(j)
          p_max = external_source % params_space(3+j)
          r = prn()
          select case(external_source % src_dist_xyz(j))
          case(SRC_DIST_UNIFORM)
            site % xyz(j) = p_min + r * (p_max - p_min)
          case(SRC_DIST_COSINE)
            site % xyz(j) = p_min &
              & + (TWO * asin(TWO * r - ONE) + PI) / (TWO * PI) * (p_max - p_min)
          end select
        end do

        ! Fill p with needed data
        p % coord0 % xyz = site % xyz
        p % coord0 % uvw = [ ONE, ZERO, ZERO ]

        ! Now search to see if location exists in geometry
        call find_cell(p, found)
        if (.not. found) then
          num_resamples = num_resamples + 1
          if (num_resamples == MAX_EXTSRC_RESAMPLES) then
            call fatal_error("Maximum number of external source spatial &
                 &resamples reached!")
          end if
        end if
      end do
      call p % clear()

    case (SRC_SPACE_FISSION)
      ! Repeat sampling source location until a good site has been found
      found = .false.
      do while (.not.found)
        ! Set particle defaults
        call p % initialize()

        ! loop over axes
        do j = 1, 3
          p_min = external_source % params_space(j)
          p_max = external_source % params_space(3+j)
          r = prn()
          select case(external_source % src_dist_xyz(j))
          case(SRC_DIST_UNIFORM)
            site % xyz(j) = p_min + r * (p_max - p_min)
          case(SRC_DIST_COSINE)
            site % xyz(j) = p_min &
              & + (TWO * asin(TWO * r - ONE) + PI) / (TWO * PI) * (p_max - p_min)
          end select
        end do

        ! Fill p with needed data
        p % coord0 % xyz = site % xyz
        p % coord0 % uvw = [ ONE, ZERO, ZERO ]

        ! Now search to see if location exists in geometry
        call find_cell(p, found)
        if (.not. found) then
          num_resamples = num_resamples + 1
          if (num_resamples == MAX_EXTSRC_RESAMPLES) then
            call fatal_error("Maximum number of external source spatial &
                 &resamples reached!")
          end if
          cycle
        end if
        if (p % material == MATERIAL_VOID) then
          found = .false.
          cycle
        end if
        if (.not. materials(p % material) % fissionable) found = .false.
      end do
      call p % clear()

    case (SRC_SPACE_POINT)
      ! Point source
      site % xyz = external_source % params_space

    end select

    ! Sample angle
    select case (external_source % type_angle)
    case (SRC_ANGLE_ISOTROPIC)
      ! Sample isotropic distribution
      phi = TWO*PI*prn()
      mu = TWO*prn() - ONE
      site % uvw(1) = mu
      site % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
      site % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

    case (SRC_ANGLE_MONO)
      ! Monodirectional source
      site % uvw = external_source % params_angle

    case default
      call fatal_error("No angle distribution specified for external source!")
    end select

    ! Sample energy distribution
    select case (external_source % type_energy)
    case (SRC_ENERGY_MONO)
      ! Monoenergtic source
      site % E = external_source % params_energy(1)

    case (SRC_ENERGY_MAXWELL)
      a = external_source % params_energy(1)
      do
        ! Sample Maxwellian fission spectrum
        site % E = maxwell_spectrum(a)

        ! resample if energy is >= 20 MeV
        if (site % E < 20) exit
      end do

    case (SRC_ENERGY_WATT)
      a = external_source % params_energy(1)
      b = external_source % params_energy(2)
      do
        ! Sample Watt fission spectrum
        site % E = watt_spectrum(a, b)

        ! resample if energy is >= 20 MeV
        if (site % E < 20) exit
      end do

    case default
      call fatal_error("No energy distribution specified for external source!")
    end select

    ! Set the random number generator back to the tracking stream.
    call prn_set_stream(STREAM_TRACKING)

  end subroutine sample_external_source

!===============================================================================
! GET_SOURCE_PARTICLE returns the next source particle
!===============================================================================

  subroutine get_source_particle(p, index_source)

    type(Particle), intent(inout) :: p
    integer(8),     intent(in)    :: index_source

    integer(8) :: particle_seed  ! unique index for particle
    integer :: i
    type(Bank), pointer, save :: src => null()
!$omp threadprivate(src)

    ! set defaults
    call p % initialize()

    ! Copy attributes from source to particle
    src => source_bank(index_source)
    call copy_source_attributes(p, src)

    ! set identifier for particle
    p % id = work_index(rank) + index_source

    ! set random number seed
    particle_seed = (overall_gen - 1)*n_particles + p % id
    call set_particle_seed(particle_seed)

    ! set particle trace
    trace = .false.
    if (current_batch == trace_batch .and. current_gen == trace_gen .and. &
         p % id == trace_particle) trace = .true.

    ! Set particle track.
    p % write_track = .false.
    if (write_all_tracks) then
      p % write_track = .true.
    else if (allocated(track_identifiers)) then
      do i=1, size(track_identifiers(1,:))
        if (current_batch == track_identifiers(1,i) .and. &
             &current_gen == track_identifiers(2,i) .and. &
             &p % id == track_identifiers(3,i)) then
          p % write_track = .true.
          exit
        end if
      end do
    end if

  end subroutine get_source_particle

!===============================================================================
! COPY_SOURCE_ATTRIBUTES
!===============================================================================

  subroutine copy_source_attributes(p, src)

    type(Particle), intent(inout) :: p
    type(Bank),     pointer       :: src

    ! copy attributes from source bank site
    p % wgt         = src % wgt
    p % last_wgt    = src % wgt
    p % coord % xyz = src % xyz
    p % coord % uvw = src % uvw
    p % last_xyz    = src % xyz
    p % last_uvw    = src % uvw
    p % E           = src % E
    p % last_E      = src % E

  end subroutine copy_source_attributes

!===============================================================================
! SRC_DIST_FUNCTION translates the string from the user input to the
! corresponding spatial source distribution
!===============================================================================

  function src_dist_function(inp) result(dist_func)

    character(MAX_LINE_LEN) :: inp
    integer :: dist_func

    select case(trim(adjustl(inp)))
    case('uniform')
      dist_func = SRC_DIST_UNIFORM

    case('cosine')
      dist_func = SRC_DIST_COSINE

    case default
      call fatal_error('Unknown spatial source distribution function')

    end select

  end function src_dist_function

end module source
