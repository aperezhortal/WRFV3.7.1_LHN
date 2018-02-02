!WRF:MEDIATION_LAYER:PHYSICS
!
! Copyright (C) 2018 Andres Perez Hortal, McGill University, Montreal, Canada
!
!> This module implements the Latent Heat Nudging scheme for radar data assimilation.

#define NEARBY_PROFILE_FOUND 0
#define NO_PROFILE_FOUND 1

MODULE module_latent_heat_nudging

   USE module_model_constants, ONLY:SVP1_Pa, SVP2, SVP3, SVPT0, EP_2
   ! Saturation vapor pressure constants for Bolton, 1980 formula
   ! ew[kPa] = 0.6112 exp[ 17.67 * T / (T+243.5) ] ; temperature T in Celsius
   ! ew[kPa] = SVP1 exp[ SVP2 * ( T - SVPT0) / (T - SVP3) ] ; temperature T in Kelvin
   ! SVP1_Pa constant for saturation vapor pressure calculation (Pa)
   ! SVP2 constant for saturation vapor pressure calculation (K)
   ! SVP3 constant for saturation vapor pressure calculation(K)
   ! SVPT0 constant for saturation vapor pressure calculation(K)
   ! EP_2=R_d/R_v constant for specific humidity calculation

   IMPLICIT NONE

   REAL, PARAMETER :: minimum_rainrate = 0.1 ! Minimum rainrate value. Lower values are considered zero [mm/h]
   REAL, PARAMETER :: rainrate_tolerance = 0.01 ! Rainrate tolerance when comparing wo precipitation values [mm/h]
   REAL, PARAMETER :: grid_points_tolerance = 0.1 ! Minimum tolerance when comparing grid points

CONTAINS

   !> Latent Heat Nudging A latent heat nudging scheme for the assimilation of precipitation data.
   !!
   !!  "It is based on latent heat nudging, in which the model profiles of latent heating are scaled by the ratio of observed
   !!and mo del precipitation rates. This causes the model to adjust so that the diagnosed precipitation rate
   !! agrees more closely with observations."
   !! Jones, C. D., & Macpherson, B. (1997). A latent heat nudging scheme for the assimilation of precipitation data
   !! into an operational mesoscale model. Meteorological Applications, 4(3), 269-277.
   SUBROUTINE latent_heat_nudging_JM97(current_time_in_seconds, dt, dx, &
                                       microphysics_latent_heating, &
                                       cumulus_latent_heating, &
                                       raincv, rainshv, rainncv, &
                                       temperature, &
                                       qvapor, pressure, &
                                       exner_function, &
                                       rainrate_observations, &
                                       lhn_start, lhn_end, &
                                       latent_heat_tendency, qvapor_adjustment_tend, &
                                       minimum_profile_scaling, maximum_profile_scaling, &
                                       search_range, &
                                       lhn_smooth_filter_sigma, &
                                       lhn_limit, &
                                       ids, ide, jds, jde, kds, kde, &
                                       ims, ime, jms, jme, kms, kme, &
                                       num_tiles, i_start, i_end, j_start, j_end, &
                                       k_start, k_end)

      IMPLICIT NONE

      !------------------------------ Subroutine Parameters Declarations --------------------------

      INTEGER, INTENT(IN) :: current_time_in_seconds ! Simulation time in seconds from the beggiging [s].
      REAL, INTENT(IN) :: dt ! Timestep [s].
      REAL, INTENT(IN) :: dx ! Horizontal resolution [m].

      !> Minimum factor to scales the model latent heating profile.
      REAL, INTENT(IN) :: minimum_profile_scaling

      !> Minimum factor to scales the model latent heating profile.
      REAL, INTENT(IN) :: maximum_profile_scaling

      !> Maximum search range in meters where a nearby profile is searched [grid points].
      ! If the range is equal to zero, no profiles are search.
      INTEGER, INTENT(IN) :: search_range

      !> Latent heat nudging absolute limit for pot. temperature [K/s]
      REAL, INTENT(IN) :: lhn_limit

      !> Latent heating tendency (Potential temperature) due to the microphysics scheme [K/s].
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: microphysics_latent_heating

      !> Latent heating tendency (Potential temperature) due to the cumulus scheme [K/s].
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: cumulus_latent_heating

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: raincv ! Time-step cumulus precipitation [mm].
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: rainshv ! Time-step shallow cumulus precipitation [mm].
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: rainncv ! Time-step non convective precipitation [mm].

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: qvapor !Vapor mixing ration [kg/kg].
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: pressure ! Total presure [Pa].
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: temperature ! Temperature [K].
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: exner_function ! Exner function.

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: rainrate_observations ! Rainrate observations [mm/h].

      INTEGER, INTENT(IN) :: lhn_start ! LHN start time from the beggining of the simulation [s].
      INTEGER, INTENT(IN) :: lhn_end ! LHN end time from the beggining of the simulation [s].

      !> Standard deviation of the gaussian filter in the horizontal dimension [grid points].
      ! If sigma = 0 no filter is applied.
      REAL, INTENT(IN) :: lhn_smooth_filter_sigma

      !> Latent heat Nudging Tendency for the Potential temperature [K/s].
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(OUT) :: latent_heat_tendency

      !> Vapor mixing ratio tendency due to the Latent heat Nudging [Kg/Kg/s].
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(OUT) :: qvapor_adjustment_tend

      INTEGER, INTENT(IN) :: ids ! Start index for i in domain.
      INTEGER, INTENT(IN) :: ide ! End index for i in domain (staggered grid).
      INTEGER, INTENT(IN) :: jds ! Start index for j in domain.
      INTEGER, INTENT(IN) :: jde ! End index for j in domain (staggered grid).
      INTEGER, INTENT(IN) :: kds ! Start index for k in domain.
      INTEGER, INTENT(IN) :: kde ! End index for  in domain (staggered grid).
      INTEGER, INTENT(IN) :: ims ! Start index for i in memory.
      INTEGER, INTENT(IN) :: ime ! End index for i in memory.
      INTEGER, INTENT(IN) :: jms ! Start index for j in memory.
      INTEGER, INTENT(IN) :: jme ! End index for j in memory.
      INTEGER, INTENT(IN) :: kms ! Start index for k in memory.
      INTEGER, INTENT(IN) :: kme ! End index for k in memory.
      INTEGER, INTENT(IN) :: num_tiles ! Number of tiles.
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: i_start ! Start index for i in each tile.
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: i_end ! End index for i in each tile.
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: j_start ! Start index for j in each tile.
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: j_end ! End index for j in each tile.
      INTEGER, INTENT(IN) :: k_start ! Memory patch start in k.
      INTEGER, INTENT(IN) :: k_end ! Memory patch end in k.

      ! -------------------End Subroutine Parameters Declarations ---------------------------------

      ! -------------------- Local variables Declarations -------------------------------------------
      REAL, DIMENSION(ims:ime, jms:jme) :: model_rainrate
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme) :: model_latent_heating !(K/s)
      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      INTEGER :: its
      INTEGER :: ite
      INTEGER :: jts
      INTEGER :: jte
      INTEGER :: kts
      INTEGER :: kte_unstag
      INTEGER :: tile
      REAL :: latent_heat_factor
      LOGICAL :: apply_latent_heat_nudging
      LOGICAL :: noModelPrecipitation
      LOGICAL :: noObservedPrecipitation
      INTEGER :: max_number_of_grid_points
      INTEGER :: profile_found
      INTEGER :: i_profile
      INTEGER :: j_profile
      REAL :: min_latent_heat_factor

      ! -------------------- End of Local variables Declarations ------------------------------------

      apply_latent_heat_nudging = (lhn_start <= current_time_in_seconds) .AND. (current_time_in_seconds <= lhn_end)

      min_latent_heat_factor = minimum_profile_scaling - 1.0

      kts = k_start

      kte_unstag = min(k_end, kde - 1) ! correspond to the end of the domain in the unstaggered Z dimension!

      ! Initialize the LH and vapor adjustment tendency to zero
      latent_heat_tendency(:, :, :) = 0.

      qvapor_adjustment_tend(:, :, :) = 0.

      ! Just in case that a rounding or trucation error is present
      IF (search_range > grid_points_tolerance ) THEN ! Search range bigger than 1, otherwise is considered as zero
         max_number_of_grid_points = search_range**2
      ELSE
         max_number_of_grid_points = 0
      ENDIF

      IF (apply_latent_heat_nudging) THEN

         model_latent_heating(:, :, :) = 0.

         write (*, *) "Applying latent heat nudging. t=", current_time_in_seconds, "seconds"

         ! First compute the model rainrate and model latent heating

!!$OMP PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!!$OMP PRIVATE ( tile ,its,ite,jts,jte, i,j,k)
         DO tile = 1, num_tiles
            its = i_start(tile)
            ite = i_end(tile)
            jts = j_start(tile)
            jte = j_end(tile)

            DO i = its, ite
            DO j = jts, jte

               model_rainrate(i, j) = (RAINCV(i, j) + RAINSHV(i, j) + RAINNCV(i, j))*3600./dt ! in mm/h

               IF (model_rainrate(i, j) < minimum_rainrate) THEN
                  model_rainrate(i, j) = 0.
               ENDIF

               DO k = kts, kte_unstag
                  model_latent_heating(i, k, j) = microphysics_latent_heating(i, k, j) + &
                                                  cumulus_latent_heating(i, k, j)

               ENDDO
            ENDDO
            ENDDO
         ENDDO
!!$OMP END PARALLEL DO

         ! Check which model grid points meet the desired criteria

!!$OMP PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!!$OMP PRIVATE ( tile ,its,ite,jts,jte, i,j,k,i_profile, j_profile,profile_found)
         DO tile = 1, num_tiles
            its = i_start(tile)
            ite = i_end(tile)
            jts = j_start(tile)
            jte = j_end(tile)

            DO i = its, ite
            DO j = jts, jte

               !---------------------------------------------------------------------------------------------
               ! Check if the target grid point meet the criteria

               noModelPrecipitation = model_rainrate(i, j) < minimum_rainrate
               noObservedPrecipitation = rainrate_observations(i, j) < minimum_rainrate

               ! If We have no precipitaiton in the model and the observations, or if the observed and the model
               ! values are similar, leave LHN as zero.
               ! That is, don't do anything. LHN was already initialized to zero.

               IF (noObservedPrecipitation) THEN
                  ! No observed precipitation treatment
                  IF (.NOT. noModelPrecipitation) THEN
                     DO k = kts, kte_unstag
                        latent_heat_tendency(i, k, j) = model_latent_heating(i, k, j)*min_latent_heat_factor
                     ENDDO

                  ENDIF

               ELSE
                  ! Some observed precipitation treatment

                  IF ((.NOT. noModelPrecipitation) .AND. & ! Avoid divisions close to zero
                      (rainrate_observations(i, j) < model_rainrate(i, j)*maximum_profile_scaling) .AND. &
                      (rainrate_observations(i, j) > model_rainrate(i, j)*minimum_profile_scaling)) THEN

                     latent_heat_factor = (rainrate_observations(i, j)/model_rainrate(i, j)) - 1.0

                     ! Scale the profile according to the criteria that the profiles matches
                     DO k = kts, kte_unstag
                        latent_heat_tendency(i, k, j) = model_latent_heating(i, k, j)*latent_heat_factor
                     ENDDO

                  ELSE
                     ! IF criteria is not met, perform an spiral search in the nearby profiles

                     CALL spiral_profile_search(i, j, max_number_of_grid_points, &
                                                rainrate_observations, model_rainrate, &
                                                minimum_profile_scaling, maximum_profile_scaling, &
                                                ids, ide, jds, jde, &
                                                ims, ime, jms, jme, &
                                                i_profile, j_profile, profile_found)

                     IF (profile_found == NEARBY_PROFILE_FOUND) THEN

                        latent_heat_factor = (rainrate_observations(i, j)/model_rainrate(i_profile, j_profile))

                        DO k = kts, kte_unstag
                           latent_heat_tendency(i, k, j) = (model_latent_heating(i_profile, k, j_profile)*latent_heat_factor - &
                                                            model_latent_heating(i, k, j))

                        ENDDO

                     ELSE
                        !profile_found == NO_PROFILE_FOUND !TODO: Use climatological profile!

                        ! No profile that match the criteria found, search unsuccessful
                        IF (rainrate_observations(i, j) > model_rainrate(i, j)) THEN
                           latent_heat_factor = maximum_profile_scaling - 1.
                        ELSE
                           latent_heat_factor = minimum_profile_scaling - 1.
                        ENDIF

                        DO k = kts, kte_unstag
                           latent_heat_tendency(i, k, j) = model_latent_heating(i, k, j)*latent_heat_factor
                        ENDDO

                     ENDIF

                  ENDIF
               ENDIF

               !---------------------------------------------------------------------------------------------

               !---------------------------------------------------------------------------------------------

            ENDDO
            ENDDO

         ENDDO
!!$OMP END PARALLEL DO

         ! Apply the Gaussian smoothing
         IF (lhn_smooth_filter_sigma > 0.0001) THEN
            CALL recursive_filter_2d(latent_heat_tendency, &
                                     lhn_smooth_filter_sigma, &
                                     ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme)
         ENDIF

         ! Keep the LHN between bounds.The lhn_limit is for temperature increments,
         ! not for potential temperature. Hence the latent heat tendency is multiplied
         ! by the exner function to convert into a temperature tendency

!!$OMP PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!!$OMP PRIVATE ( tile ,its,ite,jts,jte, i,j,k)
         DO tile = 1, num_tiles
            its = i_start(tile)
            ite = i_end(tile)
            jts = j_start(tile)
            jte = j_end(tile)

            DO i = its, ite
            DO j = jts, jte
!               DO k = kts, kte_unstag
!                  IF (latent_heat_tendency(i, k, j)*exner_function(i, k, j) > lhn_limit) THEN
!                     latent_heat_tendency(i, k, j) = lhn_limit/exner_function(i, k, j)
!                  ENDIF
!
!                  IF (latent_heat_tendency(i, k, j)*exner_function(i, k, j) < -lhn_limit) THEN
!                     latent_heat_tendency(i, k, j) = -lhn_limit/exner_function(i, k, j)
!                  ENDIF
!               ENDDO

               DO k = kts, kte_unstag
                  IF (latent_heat_tendency(i, k, j) > lhn_limit) THEN
                     latent_heat_tendency(i, k, j) = lhn_limit
                  ENDIF

                  IF (latent_heat_tendency(i, k, j) < -lhn_limit) THEN
                     latent_heat_tendency(i, k, j) = -lhn_limit
                  ENDIF
               ENDDO

            ENDDO
            ENDDO
         ENDDO
!!$OMP END PARALLEL DO

         ! Adjust humidity to keep the Relative humidity constant
         CALL adjust_humidity(temperature, pressure, qvapor, exner_function, dt, &
                              model_latent_heating, latent_heat_tendency, &
                              qvapor_adjustment_tend, &
                              ids, ide - 1, jds, jde - 1, kds, kde - 1, &
                              ims, ime, jms, jme, kms, kme, &
                              num_tiles, i_start, i_end, j_start, j_end, &
                              kts, kte_unstag)

         write (*, *)
         write (*, *) "minimum_profile_scaling", minimum_profile_scaling
         write (*, *) "maximum_profile_scaling", maximum_profile_scaling
         write (*, *)
         write (*, *) "Min latent_heat_tendency:", MINVAL(latent_heat_tendency(ids:ide, kds:kde, jds:jde))
         write (*, *) "Min cumulus_latent_heating CU:", MINVAL(cumulus_latent_heating(ids:ide, kds:kde, jds:jde))
         write (*, *) "Min microphysics_latent_heating CU:", MINVAL(microphysics_latent_heating(ids:ide, kds:kde, jds:jde))
         write (*, *)
         write (*, *) "Max latent_heat_tendency:", MAXVAL(latent_heat_tendency(ids:ide, kds:kde, jds:jde))
         write (*, *) "Max cumulus_latent_heating CU:", MAXVAL(cumulus_latent_heating(ids:ide, kds:kde, jds:jde))
         write (*, *) "Max microphysics_latent_heating CU:", MAXVAL(microphysics_latent_heating(ids:ide, kds:kde, jds:jde))
         write (*, *)
         write (*, *) "Min model_rainrate:", MINVAL(model_rainrate(ids:ide - 1, jds:jde - 1))
         write (*, *) "Min rainrate_observations:", MINVAL(rainrate_observations(ids:ide - 1, jds:jde - 1))
         write (*, *)
         write (*, *) "Max model_rainrate:", MAXVAL(model_rainrate(ids:ide - 1, jds:jde - 1))
         write (*, *) "Max rainrate_observations:", MAXVAL(rainrate_observations(ids:ide - 1, jds:jde - 1))
         write (*, *)
         write (*, *) "Min qvapor_adjustment_tend:", MINVAL(qvapor_adjustment_tend(ids:ide, kds:kde, jds:jde))*1000
         write (*, *) "Max qvapor_adjustment_tend:", MAXVAL(qvapor_adjustment_tend(ids:ide, kds:kde, jds:jde))*1000

      ENDIF

   END SUBROUTINE latent_heat_nudging_JM97

!> Algorithm for iterating over an outward spiral on a discrete 2D grid from the origin.
! The observed precipitation is assument greater than zero!
! Adapted from:
! https://stackoverflow.com/questions/3706219/algorithm-for-iterating-over-an-outward-spiral-on-a-discrete-2d-grid-from-the-or
   SUBROUTINE spiral_profile_search(i_center, j_center, max_number_of_grid_points, &
                                    rainrate_observations, model_rainrate, &
                                    minimum_profile_scaling, maximum_profile_scaling, &
                                    ids, ide, jds, jde, &
                                    ims, ime, jms, jme, &
                                    i_profile, j_profile, profile_found)
      IMPLICIT NONE

      !------------------------------ Subroutine Parameters Declarations --------------------------

      INTEGER, INTENT(IN) :: i_center ! i index corresponding to the center of the spiral.
      INTEGER, INTENT(IN) :: j_center ! j index corresponding to the center of the spiral.
      INTEGER, INTENT(IN) :: max_number_of_grid_points ! Maximum number of points to do the spiral search.

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: model_rainrate ! Model Rainrate [mm/h].
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: rainrate_observations ! Rainrate observations [mm/h].

      !> Minimum factor to scales the model latent heating profile.
      REAL, INTENT(IN) :: minimum_profile_scaling

      !> Minimum factor to scales the model latent heating profile.
      REAL, INTENT(IN) :: maximum_profile_scaling

      INTEGER, INTENT(IN) :: ids ! Start index for i in domain.
      INTEGER, INTENT(IN) :: ide ! End index for i in domain (staggered grid).
      INTEGER, INTENT(IN) :: jds ! Start index for j in domain.
      INTEGER, INTENT(IN) :: jde ! End index for j in domain (staggered grid).
      INTEGER, INTENT(IN) :: ims ! Start index for i in memory.
      INTEGER, INTENT(IN) :: ime ! End index for i in memory.
      INTEGER, INTENT(IN) :: jms ! Start index for j in memory.
      INTEGER, INTENT(IN) :: jme ! End index for j in memory.

      INTEGER, INTENT(OUT) :: i_profile ! i index of the profile found (if any).
      INTEGER, INTENT(OUT) :: j_profile ! j index of the profile found (if any).
      INTEGER, INTENT(OUT) :: profile_found ! Flag indicating the result of the profile search.

      ! -------------------End Subroutine Parameters Declarations ---------------------------------

      ! -------------------- Local variables Declarations -------------------------------------------
      INTEGER :: di
      INTEGER :: dj
      INTEGER :: buffer
      INTEGER :: number_of_steps
      INTEGER :: steps
      INTEGER :: segment_points_visited
      INTEGER :: segment_length

      ! -------------------- End of Local variables Declarations ------------------------------------

      IF (max_number_of_grid_points < grid_points_tolerance) THEN
         profile_found = NO_PROFILE_FOUND
         RETURN
      ENDIF

      di = 1
      dj = 0
      number_of_steps = 1

      i_profile = i_center
      j_profile = j_center

      profile_found = NO_PROFILE_FOUND

      ! Search for nearby profiles

      segment_points_visited = 0
      segment_length = 1

      DO steps = 0, max_number_of_grid_points

         ! make a step, add 'direction' vector (di, dj) to current position (i, j)
         i_profile = i_profile + di
         j_profile = j_profile + dj

         ! Ignore points outside domain
         if ((i_profile < ids) .OR. &
             (i_profile > ide) .OR. &
             (j_profile < jds) .OR. &
             (j_profile > jde)) THEN

            CYCLE
         ENDIF

         segment_points_visited = segment_points_visited + 1

         if ((model_rainrate(i_profile, j_profile) > minimum_rainrate) .AND. & ! Avoid divisions close to zero
             (rainrate_observations(i_center, j_center) < model_rainrate(i_profile, j_profile)*maximum_profile_scaling) .AND. &
             (rainrate_observations(i_center, j_center) > model_rainrate(i_profile, j_profile)*minimum_profile_scaling)) THEN

            profile_found = NEARBY_PROFILE_FOUND

            RETURN

         ENDIF

         IF (segment_points_visited == segment_length) THEN
            ! Done with current segment
            segment_points_visited = 0

            ! 'rotate' directions
            buffer = di
            di = -dj
            dj = buffer

            ! Increase segment length if necessary
            IF (dj == 0) THEN
               segment_length = segment_length + 1
            ENDIF

         ENDIF

      ENDDO

   END SUBROUTINE spiral_profile_search

!> Gaussian Fiilter: 3rd order recursive filer to approximate a 2D convolution with a Gaussian Kernel of a given
! standard deviation.
   SUBROUTINE recursive_filter_2d(input_3d, &
                                  sigma_in_gridpoints, &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme)

      IMPLICIT NONE

      !------------------------------ Subroutine Parameters Declarations --------------------------

      !> Input field to convolve with a Gaussian Kernel
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: input_3d

      !> Standard deviation of the gaussian filter in the horizontal dimension [grid points].
      ! If sigma = 0 no filter is applied.
      REAL, INTENT(IN) :: sigma_in_gridpoints

      INTEGER, INTENT(IN) :: ids ! Start index for i in domain.
      INTEGER, INTENT(IN) :: ide ! End index for i in domain (staggered grid).
      INTEGER, INTENT(IN) :: jds ! Start index for j in domain.
      INTEGER, INTENT(IN) :: jde ! End index for j in domain (staggered grid).
      INTEGER, INTENT(IN) :: kds ! Start index for k in domain.
      INTEGER, INTENT(IN) :: kde ! End index for  in domain (staggered grid).
      INTEGER, INTENT(IN) :: ims ! Start index for i in memory.
      INTEGER, INTENT(IN) :: ime ! End index for i in memory.
      INTEGER, INTENT(IN) :: jms ! Start index for j in memory.
      INTEGER, INTENT(IN) :: jme ! End index for j in memory.
      INTEGER, INTENT(IN) :: kms ! Start index for k in memory.
      INTEGER, INTENT(IN) :: kme ! End index for k in memory.
      ! -------------------End Subroutine Parameters Declarations ---------------------------------

      ! -------------------- Local variables Declarations -------------------------------------------
      INTEGER :: i
      INTEGER :: j
      INTEGER :: k

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme) :: forward_filter !, backward_filter
      REAL, DIMENSION(1:3) :: alpha
      REAL :: beta
      REAL :: ai
      ! -------------------- End of Local variables Declarations ------------------------------------

      ai = 3.738128+5.788982*sigma_in_gridpoints + 3.382473*sigma_in_gridpoints**2 + sigma_in_gridpoints**3

      alpha(1) = (5.788982*sigma_in_gridpoints + 6.764946*(sigma_in_gridpoints**2) + 3*(sigma_in_gridpoints**3))/ai
      alpha(2) = -(3.382473*(sigma_in_gridpoints**2) + 3*(sigma_in_gridpoints**3))/ai
      alpha(3) = (sigma_in_gridpoints**3)/ai

      beta = 3.738128/ai

!!$OMP PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!!$OMP PRIVATE ( k ,j, i)
      DO k = kds, kde

         ! Apply gaussian filter in the "x" direction for each j
         DO j = jds, jde

            ! Apply forward filter

            forward_filter(ids, k, j) = beta*input_3d(ids, k, j)

            forward_filter(ids + 1, k, j) = beta*input_3d(ids + 1, k, j) + alpha(1)*forward_filter(ids, k, j)

            forward_filter(ids + 2, k, j) = beta*input_3d(ids + 2, k, j) + &
                                            alpha(1)*forward_filter(ids + 1, k, j) + &
                                            alpha(2)*forward_filter(ids, k, j)

            DO i = ids + 3, ide

               forward_filter(i, k, j) = beta*input_3d(i, k, j) + &
                                         alpha(1)*forward_filter(i - 1, k, j) + &
                                         alpha(2)*forward_filter(i - 2, k, j) + &
                                         alpha(3)*forward_filter(i - 3, k, j)
            ENDDO

            ! Apply backward filter "x" direction

            input_3d(ide, k, j) = beta*forward_filter(ide, k, j)

            input_3d(ide - 1, k, j) = beta*forward_filter(ide - 1, k, j) + alpha(1)*input_3d(ide, k, j)

            input_3d(ide - 2, k, j) = beta*forward_filter(ide - 2, k, j) + &
                                      alpha(1)*input_3d(ide - 1, k, j) + &
                                      alpha(2)*input_3d(ide, k, j)

            DO i = ide - 3, ids, -1
               input_3d(i, k, j) = beta*forward_filter(i, k, j) + &
                                   alpha(1)*input_3d(i + 1, k, j) + &
                                   alpha(2)*input_3d(i + 2, k, j) + &
                                   alpha(3)*input_3d(i + 3, k, j)
            ENDDO

         ENDDO

         ! Apply gaussian filter in the "y" direction for each i
         DO i = ids, ide

            ! Apply forward filter

            forward_filter(i, k, jds) = beta*input_3d(i, k, jds)

            forward_filter(i, k, jds + 1) = beta*input_3d(i, k, jds + 1) + alpha(1)*forward_filter(i, k, jds)

            forward_filter(i, k, jds + 1) = beta*input_3d(i, k, jds + 2) + &
                                            alpha(1)*forward_filter(i, k, jds + 1) + &
                                            alpha(2)*forward_filter(i, k, jds)

            DO j = jds + 3, jde

               forward_filter(i, k, j) = beta*input_3d(i, k, j) + &
                                         alpha(1)*forward_filter(i, k, j - 1) + &
                                         alpha(2)*forward_filter(i, k, j - 2) + &
                                         alpha(3)*forward_filter(i, k, j - 3)
            ENDDO

            ! Apply backward filter "x" direction

            input_3d(i, k, jde) = beta*forward_filter(i, k, jde)

            input_3d(i, k, jde - 1) = beta*forward_filter(i, k, jde - 1) + alpha(1)*input_3d(i, k, jde)

            input_3d(i, k, jde - 2) = beta*forward_filter(i, k, jde - 2) + &
                                      alpha(1)*input_3d(i, k, jde - 1) + &
                                      alpha(2)*input_3d(i, k, jde)

            DO j = jde - 3, jds, -1
               input_3d(i, k, j) = beta*forward_filter(i, k, j) + &
                                   alpha(1)*input_3d(i, k, j + 1) + &
                                   alpha(2)*input_3d(i, k, j + 2) + &
                                   alpha(3)*input_3d(i, k, j + 3)
            ENDDO

         ENDDO

      ENDDO
!!$OMP END PARALLEL DO

   END SUBROUTINE recursive_filter_2d

   SUBROUTINE adjust_humidity(temperature, pressure, qvapor, exner_function, dt, &
                              model_latent_heating, latent_heat_tendency, &
                              qvapor_adjustment_tend, &
                              ids, ide, jds, jde, kds, kde, &
                              ims, ime, jms, jme, kms, kme, &
                              num_tiles, i_start, i_end, j_start, j_end, &
                              kts, kte_unstag)

      IMPLICIT NONE
      !------------------------------ Subroutine Parameters Declarations --------------------------

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: qvapor !Vapor mixing ration [kg/kg].
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: pressure ! Total presure [Pa].
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: temperature ! Temperature [K].
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: exner_function ! Exner function.

      !> Total model latent heating (potential temperature) due to Cumulus and Microphysics parametrizations [K/s].
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: model_latent_heating

      !> Latent heat Nudging Tendency for the Potential temperature [K/s].
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(OUT) :: latent_heat_tendency

      !> Vapor mixing ratio tendency due to the Latent heat Nudging [Kg/Kg/s].
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(OUT) :: qvapor_adjustment_tend

      REAL, INTENT(IN) :: dt ! Timestep [s].

      INTEGER, INTENT(IN) :: num_tiles ! Number of tiles.
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: i_start ! Start index for i in each tile.
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: i_end ! End index for i in each tile.
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: j_start ! Start index for j in each tile.
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: j_end ! End index for j in each tile.
      INTEGER, INTENT(IN) :: kts ! Memory patch start in z direction.
      INTEGER, INTENT(IN) :: kte_unstag ! Memory patch end in z direction (unstaggered).
      INTEGER, INTENT(IN) :: ids ! Start index for i in domain.
      INTEGER, INTENT(IN) :: ide ! End index for i in domain (staggered grid).
      INTEGER, INTENT(IN) :: jds ! Start index for j in domain.
      INTEGER, INTENT(IN) :: jde ! End index for j in domain (staggered grid).
      INTEGER, INTENT(IN) :: kds ! Start index for k in domain.
      INTEGER, INTENT(IN) :: kde ! End index for  in domain (staggered grid).
      INTEGER, INTENT(IN) :: ims ! Start index for i in memory.
      INTEGER, INTENT(IN) :: ime ! End index for i in memory.
      INTEGER, INTENT(IN) :: jms ! Start index for j in memory.
      INTEGER, INTENT(IN) :: jme ! End index for j in memory.
      INTEGER, INTENT(IN) :: kms ! Start index for k in memory.
      INTEGER, INTENT(IN) :: kme ! End index for k in memory.

      !------------------------------ Local Parameters Declarations --------------------------

      REAL :: newTemperature !  Expected temperature at the end of the timestep due to model latent heating.
      REAL :: newNudgedTemperature ! Expected temperature at the end of the timestep after nudging.
      REAL :: es ! Saturation vapor pressure
      REAL :: lhn_qvapor_saturation
      REAL :: model_qvapor_saturation
      REAL, PARAMETER :: no_heating_tolerance = 1e-10 ![k/s]

      INTEGER :: i, j, k, tile
      INTEGER :: its
      INTEGER :: ite
      INTEGER :: jts
      INTEGER :: jte

!!$OMP PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!!$OMP PRIVATE ( tile,its,ite,jts,jte,i,j,k)&
!!$OMP PRIVATE ( newTemperature,newNudgedTemperature,es,model_qvapor_saturation,lhn_qvapor_saturation)
      DO tile = 1, num_tiles
         its = i_start(tile)
         ite = i_end(tile)
         jts = j_start(tile)
         jte = j_end(tile)

         DO i = its, ite
         DO j = jts, jte

            DO k = kts, kte_unstag

               IF (ABS(latent_heat_tendency(i, k, j)) > no_heating_tolerance) THEN

                  newTemperature = temperature(i, k, j) + model_latent_heating(i, k, j)*exner_function(i, k, j)*dt

                  newNudgedTemperature = temperature(i, k, j) + &
                                         (model_latent_heating(i, k, j) + &
                                          latent_heat_tendency(i, k, j))*exner_function(i, k, j)*dt

                  ! Bolton, 1980 formula
                  es = SVP1_Pa*EXP(SVP2*(newTemperature - SVPT0)/(newTemperature - SVP3))
                  model_qvapor_saturation = EP_2*es/(pressure(i, k, j) - es)

                  es = SVP1_Pa*EXP(SVP2*(newNudgedTemperature - SVPT0)/(newNudgedTemperature - SVP3))
                  lhn_qvapor_saturation = EP_2*es/(pressure(i, k, j) - es)

                  qvapor_adjustment_tend(i, k, j) = qvapor(i, k, j)*(lhn_qvapor_saturation/model_qvapor_saturation - 1.0)/dt
               ELSE
                  qvapor_adjustment_tend(i, k, j) = 0.
               ENDIF

            ENDDO
         ENDDO
         ENDDO
      ENDDO
!!$OMP END PARALLEL DO

   END SUBROUTINE adjust_humidity

END MODULE module_latent_heat_nudging
