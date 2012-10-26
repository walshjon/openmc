module tally_diffusion

  use constants

  ! module options
  private
  public :: create_diffusion_tally

  ! number of energy groups
  integer, parameter :: N_GRPS = 70

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
                             n_user_analog_tallies
    use mesh_header,   only: StructuredMesh
    use tally_header,  only: TallyObject

    ! arguments
    type(TallyObject) :: t

    ! local variables
    integer :: n_filters
    integer :: filters(N_FILTER_TYPES)
    integer :: i_mesh 
    character(MAX_WORD_LEN) :: hydrogen
    type(StructuredMesh), pointer :: m => null()

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

    ! allocate arrays for number of bins and stride in scores array
    allocate(t % n_filter_bins(N_FILTER_TYPES))
    allocate(t % stride(N_FILTER_TYPES))

    ! initialize number of bins and stride
    t % n_filter_bins = 0
    t % stride = 0

    ! initialize filters
    n_filters = 0
    filters = 0

    ! set tally type to volume
    t % type = TALLY_VOLUME

    ! set estimator to analog
    t % estimator = ESTIMATOR_ANALOG

    ! set tally id
    t % id = 900

    ! set tally label
    t % label = 'DIFFUSION COEFFICIENT'

    ! set mesh filter bins
    t % n_filter_bins(FILTER_MESH) = t % n_filter_bins(FILTER_MESH) + &
                                     product(m % dimension)
    n_filters = n_filters + 1
    filters(n_filters) = FILTER_MESH

    ! set energy in filter bins
    allocate(t % energy_in(N_GRPS+1))
    t % energy_in = EGRID
    t % n_filter_bins(FILTER_ENERGYIN) = N_GRPS
    n_filters = n_filters + 1
    filters(n_filters) = FILTER_ENERGYIN

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

    ! increment number of analog tallies
    n_user_analog_tallies = n_user_analog_tallies + 1

  end subroutine create_diffusion_tally

end module tally_diffusion
