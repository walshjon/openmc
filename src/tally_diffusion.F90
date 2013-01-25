module tally_diffusion

  use constants

  ! module options
  private
  public :: create_diffusion_tally, calculate_diffusion

  ! number of energy groups
  integer, parameter :: N_GRPS = 70

  ! set up vectors for tallies
  real(8), allocatable :: flux(:,:,:,:)
  real(8), allocatable :: total(:,:,:,:)
  real(8), allocatable :: totalH1(:,:,:,:)
  real(8), allocatable :: p1scatt(:,:,:,:)
  real(8), allocatable :: p1scattH1(:,:,:,:)
  real(8), allocatable :: diffcoef(:,:,:,:)
  real(8), allocatable :: trans(:,:,:,:)
  real(8), allocatable :: transH1(:,:,:,:)

  ! energy grid for correction
  real(8), parameter :: EGRID(N_GRPS+1) = (/ 0.0_8,    &
                                   5.0000e-09_8, &
                                   1.0000e-08_8, &
                                   1.5000e-08_8, &
                                   2.0000e-08_8, &
                                   2.5000e-08_8, &
                                   3.0000e-08_8, &
                                   3.5000e-08_8, &
                                   4.2000e-08_8, &
                                   5.0000e-08_8, &
                                   5.8000e-08_8, &
                                   6.7000e-08_8, &
                                   8.0000e-08_8, &
                                   1.0000e-07_8, &
                                   1.4000e-07_8, &
                                   1.8000e-07_8, &
                                   2.2000e-07_8, &
                                   2.5000e-07_8, &
                                   2.8000e-07_8, &
                                   3.0000e-07_8, &
                                   3.2000e-07_8, &
                                   3.5000e-07_8, &
                                   4.0000e-07_8, &
                                   5.0000e-07_8, &
                                   6.2500e-07_8, &
                                   7.8000e-07_8, &
                                   8.5000e-07_8, &
                                   9.1000e-07_8, &
                                   9.5000e-07_8, &
                                   9.7200e-07_8, &
                                   9.9600e-07_8, &
                                   1.0200e-06_8, &
                                   1.0450e-06_8, &
                                   1.0710e-06_8, &
                                   1.0970e-06_8, &
                                   1.1230e-06_8, &
                                   1.1500e-06_8, &
                                   1.3000e-06_8, &
                                   1.5000e-06_8, &
                                   1.8550e-06_8, &
                                   2.1000e-06_8, &
                                   2.6000e-06_8, &
                                   3.3000e-06_8, &
                                   4.0000e-06_8, &
                                   9.8770e-06_8, &
                                   1.5968e-05_8, &
                                   2.7700e-05_8, &
                                   4.8052e-05_8, &
                                   7.5501e-05_8, &
                                   1.4873e-04_8, &
                                   3.6726e-04_8, &
                                   9.0690e-04_8, &
                                   1.4251e-03_8, &
                                   2.2395e-03_8, &
                                   3.5191e-03_8, &
                                   5.5300e-03_8, &
                                   9.1180e-03_8, &
                                   1.5030e-02_8, &
                                   2.4780e-02_8, &
                                   4.0850e-02_8, &
                                   6.7340e-02_8, &
                                   1.1100e-01_8, &
                                   1.8300e-01_8, &
                                   3.0250e-01_8, &
                                   5.0000e-01_8, &
                                   8.2100e-01_8, &
                                   1.3530e+00_8, &
                                   2.2310e+00_8, &
                                   3.6790e+00_8, &
                                   6.0655e+00_8, &
                                   1.0000e+01_8 /)

  ! correction curve 
  real(8), parameter :: inscatt(N_GRPS) = (/               &
                                         9.5199e-01_8, &
                                         1.0351e+00_8, &
                                         9.3946e-01_8, &
                                         8.7042e-01_8, &
                                         8.3910e-01_8, &
                                         8.3244e-01_8, &
                                         7.8438e-01_8, &
                                         7.7327e-01_8, &
                                         7.3919e-01_8, &
                                         7.1333e-01_8, &
                                         7.1264e-01_8, &
                                         6.8051e-01_8, &
                                         6.7389e-01_8, &
                                         6.4561e-01_8, &
                                         5.8493e-01_8, &
                                         5.0988e-01_8, &
                                         4.8314e-01_8, &
                                         4.0685e-01_8, &
                                         4.3215e-01_8, &
                                         4.0903e-01_8, &
                                         3.8992e-01_8, &
                                         3.6148e-01_8, &
                                         3.6706e-01_8, &
                                         3.5664e-01_8, &
                                         3.5371e-01_8, &
                                         3.3869e-01_8, &
                                         3.4995e-01_8, &
                                         3.6857e-01_8, &
                                         3.4995e-01_8, &
                                         3.2715e-01_8, &
                                         3.2658e-01_8, &
                                         3.5098e-01_8, &
                                         3.2105e-01_8, &
                                         3.2322e-01_8, &
                                         3.8905e-01_8, &
                                         3.5944e-01_8, &
                                         3.4043e-01_8, &
                                         3.3408e-01_8, &
                                         3.3611e-01_8, &
                                         3.3832e-01_8, &
                                         3.2368e-01_8, &
                                         3.2260e-01_8, &
                                         3.2856e-01_8, &
                                         3.2824e-01_8, &
                                         3.2864e-01_8, &
                                         3.2942e-01_8, &
                                         3.2221e-01_8, &
                                         3.1897e-01_8, &
                                         3.1546e-01_8, &
                                         3.1089e-01_8, &
                                         2.9724e-01_8, &
                                         2.8777e-01_8, &
                                         2.7847e-01_8, &
                                         2.6265e-01_8, &
                                         2.6034e-01_8, &
                                         2.4364e-01_8, &
                                         2.3206e-01_8, &
                                         2.2144e-01_8, &
                                         2.1931e-01_8, &
                                         2.1471e-01_8, &
                                         2.2009e-01_8, &
                                         2.3376e-01_8, &
                                         2.6024e-01_8, &
                                         3.0118e-01_8, &
                                         3.6275e-01_8, &
                                         4.4399e-01_8, &
                                         5.4645e-01_8, &
                                         6.5621e-01_8, &
                                         7.5675e-01_8, &
                                         8.3578e-01_8 /)

contains

!===============================================================================
! CREATE_DIFFUSION_TALLY
!===============================================================================

  subroutine create_diffusion_tally(t)

    use datatypes,     only: dict_has_key, dict_get_key
    use error,         only: fatal_error
    use global,        only: default_xs, nuclide_dict, message, &
                             difcof_mesh, mesh_dict, meshes, &
                             n_user_analog_tallies, n_user_tallies, &
                             analog_tallies
    use mesh_header,   only: StructuredMesh
    use tally_header,  only: TallyObject, TallyScore, TallyFilter

    ! arguments
    type(TallyObject) :: t

    ! local variables
    integer :: n_filters
    integer :: i_mesh 
    character(MAX_WORD_LEN) :: hydrogen
    type(StructuredMesh), pointer :: m => null()
    type(TallyFilter) :: filters(N_FILTER_TYPES)

    ! check if hydogen is included in problem
    hydrogen = "H-1" // "." // default_xs
    if (.not. dict_has_key(nuclide_dict, hydrogen)) then
      message = "Could not find nuclide " // trim(hydrogen) // &
                " needed for diffusion coefficients."
      call fatal_error()
    end if

    ! check if mesh exists
    if (.not. dict_has_key(mesh_dict, difcof_mesh)) then
      message = "Mesh for diffusion coefficient does not exist."
      call fatal_error()
    else
      i_mesh = dict_get_key(mesh_dict, difcof_mesh)
      m => meshes(i_mesh)
    end if

    ! initialize filters
    n_filters = 0

    ! set tally type to volume
    t % type = TALLY_VOLUME

    ! set estimator to analog
    t % estimator = ESTIMATOR_ANALOG

    ! set tally id
    t % id = n_user_tallies + 4 

    ! set tally label
    t % label = 'DIFFUSION COEFFICIENT'

    ! set mesh filter bins
    n_filters = n_filters + 1
    filters(n_filters) % type = FILTER_MESH
    filters(n_filters) % n_bins = product(m % dimension)
    allocate(filters(n_filters) % int_bins(1))
    filters(n_filters) % int_bins(1) = i_mesh
    t % find_filter(FILTER_MESH) = n_filters

    ! set energy in filter bins
    n_filters = n_filters + 1
    filters(n_filters) % type = FILTER_ENERGYIN
    filters(n_filters) % n_bins = N_GRPS
    allocate(filters(n_filters) % real_bins(N_GRPS+1))
    filters(n_filters) % real_bins = EGRID
    t % find_filter(FILTER_ENERGYIN) = n_filters

    ! allocate and set filters
    t % n_filters = n_filters
    allocate(t % filters(n_filters))
    t % filters = filters(1:n_filters)

    ! set nuclide data
    allocate(t % nuclide_bins(2)) ! 1 for hydrogen and 2 for total
    t % n_nuclide_bins = 2
    t % nuclide_bins(1) = dict_get_key(nuclide_dict, hydrogen)
    t % nuclide_bins(2) = -1

    ! set scores
    allocate(t % score_bins(3))
    t % n_score_bins = 3
    t % score_bins(1) = SCORE_FLUX
    t % score_bins(2) = SCORE_TOTAL
    t % score_bins(3) = SCORE_SCATTER_1

  end subroutine create_diffusion_tally

!===============================================================================
! CALCULATE_DIFFUSION
!===============================================================================

  subroutine calculate_diffusion(diff_out,ng,nx,ny,nz)

    use datatypes,     only: dict_get_key
    use error,         only: fatal_error
    use global,        only: meshes, tallies, difcof_mesh, n_user_tallies, &
                             mesh_dict, message
    use mesh,          only: mesh_indices_to_bin
    use mesh_header,   only: StructuredMesh
    use tally_header,  only: TallyObject

    ! arguments
    integer :: ng                     ! number of energy groups
    integer :: nx                     ! number of mesh cells in z direction
    integer :: ny                     ! number of mesh cells in y direction
    integer :: nz                     ! number of mesh cells in x direction
    real(8) :: diff_out(2,nx,ny,nz)   ! output diffusion coefficients

    ! local variables
    integer :: g                    ! iteration counter for groups
    integer :: i                    ! iteration counter for x
    integer :: j                    ! iteration counter for y
    integer :: k                    ! iteration counter for z
    integer :: i_mesh               ! mesh index
    integer :: ijk(3)               ! indices for mesh cells
    integer :: filter_index         ! index to pull from tally object
    integer :: score_index          ! index to pull from tally object
    integer :: bins(N_FILTER_TYPES) ! bins for filters
    integer :: i_score              ! score number
    integer :: i_nuclide            ! nuclide number
    integer :: i_filter_mesh        ! index for mesh filter
    integer :: i_filter_ein         ! index for incoming energy filter
    integer :: i_filter_nuclide     ! index for nuclide filter
    type(StructuredMesh), pointer :: m => null()
    type(TallyObject), pointer :: t => null()

    ! set pointers
    t => tallies(n_user_tallies + 4)
    i_mesh = t % filters(t % find_filter(FILTER_MESH)) % int_bins(1)
    m => meshes(i_mesh)

    ! set up filters
    i_filter_mesh = t % find_filter(FILTER_MESH)
    i_filter_ein  = t % find_filter(FILTER_ENERGYIN)

    ! allocate variables
    allocate(flux(N_GRPS,nx,ny,nz))
    allocate(total(N_GRPS,nx,ny,nz))
    allocate(totalH1(N_GRPS,nx,ny,nz))
    allocate(p1scatt(N_GRPS,nx,ny,nz))
    allocate(p1scattH1(N_GRPS,nx,ny,nz))
    allocate(diffcoef(N_GRPS,nx,ny,nz))
    allocate(trans(N_GRPS,nx,ny,nz))
    allocate(transH1(N_GRPS,nx,ny,nz))

    ! begin loop around space and energy
    ZLOOP: do k = 1,nz

      YLOOP: do j = 1,ny

        XLOOP: do i = 1,nx

          GLOOP: do g = 1,N_GRPS

            ! reset all bins to 1
            t % matching_bins = 1

            ! set ijk as mesh indices
            ijk = (/i,j,k/)
            t % matching_bins(i_filter_mesh) = mesh_indices_to_bin(m,ijk)

            ! apply energy in filter
            t % matching_bins(i_filter_ein) = g

            ! calculate filter index from bins
            filter_index = sum((t % matching_bins - 1) * t%stride) + 1

            ! calculate score index for H-1 total
            i_nuclide = 1
            i_score = 2
            score_index = (i_nuclide - 1)*t % n_score_bins + i_score 
            totalH1(g,i,j,k) = t % scores(score_index,filter_index) % sum 

            ! calculate score index for H-1 p1scatt
            i_nuclide = 1
            i_score = 3
            score_index = (i_nuclide - 1)*t % n_score_bins + i_score
            p1scattH1(g,i,j,k) = t % scores(score_index,filter_index) % sum

            ! calculate score index for all nuclides flux
            i_nuclide = 2
            i_score = 1
            score_index = (i_nuclide - 1)*t % n_score_bins + i_score
            flux(g,i,j,k) = t % scores(score_index,filter_index) % sum

            ! calculate score index for all nuclides total 
            i_nuclide = 2
            i_score = 2
            score_index = (i_nuclide - 1)*t % n_score_bins + i_score
            total(g,i,j,k) = t % scores(score_index,filter_index) % sum

            ! calculate score index for all nuclides flux
            i_nuclide = 2
            i_score = 3
            score_index = (i_nuclide - 1)*t % n_score_bins + i_score
            p1scatt(g,i,j,k) = t % scores(score_index,filter_index) % sum 

          end do GLOOP

          ! if flux is zero skip
          if (abs(minval(flux(:,i,j,k)) - ZERO) < TINY_BIT) then
            diff_out(:,i,j,k) = ZERO_FLUX 
            cycle
          end if

          ! compute transport reaction rate
          trans(:,i,j,k) = total(:,i,j,k) - p1scatt(:,i,j,k)
          transH1(:,i,j,k) = totalH1(:,i,j,k) - p1scattH1(:,i,j,k)

          ! subtract out hydrogen from transport reaction rate
          trans(:,i,j,k) = trans(:,i,j,k) - transH1(:,i,j,k)

          ! apply in-scatter correction curve to hydrogen total
          transH1(:,i,j,k) = totalH1(:,i,j,k)*inscatt

          ! add in corrected hydrogen to transport reaction rate
          trans(:,i,j,k) = trans(:,i,j,k) + transH1(:,i,j,k)

          ! compute transport cross section
          trans(:,i,j,k) = trans(:,i,j,k) / flux(:,i,j,k)

          ! compute diffusion coefficients
          diffcoef(:,i,j,k) = 1/(3*trans(:,i,j,k))

          ! collapse diffusion coefficient
          if (ng == 2) then
            diff_out(2,i,j,k) = sum(diffcoef(1:24,i,j,k)*flux(1:24,i,j,k)) / &
                                sum(flux(1:24,i,j,k))
            diff_out(1,i,j,k) = sum(diffcoef(25:70,i,j,k)*flux(25:70,i,j,k)) / &
                                sum(flux(25:70,i,j,k))
          else if (ng == 1) then
            diff_out(1,i,j,k) = sum(diffcoef(1:70,i,j,k)*flux(1:70,i,j,k)) / &
                                sum(flux(1:70,i,j,k))
          else
            message = 'Only 1 or 2 group calculation allowed for diff. coeffs'
            call fatal_error()
          end if

        end do XLOOP

      end do YLOOP

    end do ZLOOP

    ! deallocate variables
    deallocate(flux)
    deallocate(total)
    deallocate(totalH1)
    deallocate(p1scatt)
    deallocate(p1scattH1)
    deallocate(diffcoef)
    deallocate(trans)
    deallocate(transH1)

  end subroutine calculate_diffusion

end module tally_diffusion
