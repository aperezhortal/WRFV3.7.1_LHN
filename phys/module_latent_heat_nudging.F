!WRF:MEDIATION_LAYER:PHYSICS
!

#define PROFILE_CRITERIA_MEET 0
#define NEARBY_PROFILE_FOUND 1
#define NO_PROFILE_FOUND 2

MODULE module_latent_heat_nudging

   IMPLICIT NONE

CONTAINS

   ! Jones, C. D., & Macpherson, B. (1997). A latent heat nudging scheme for the assimilation of precipitation data
   ! into an operational mesoscale model. Meteorological Applications, 4(3), 269-277.
   SUBROUTINE latent_heat_nudging_1(current_time_in_seconds, dt, dx, &
                                    microphysics_latent_heating, &
                                    cumulus_latent_heating, &
                                    model_precip_accum, precip_accum_obs, &
                                    lhn_start, lhn_end, latent_heat_tendency, &
                                    minimum_profile_scaling, maximum_profile_scaling, &
                                    search_range, &
                                    ids, ide, jds, jde, kds, kde, &
                                    ims, ime, jms, jme, kms, kme, &
                                    num_tiles, i_start, i_end, j_start, j_end, &
                                    k_start, k_end)

      IMPLICIT NONE

      !------------------------------ Subroutine Parameters Declarations --------------------------

      INTEGER, INTENT(IN) :: current_time_in_seconds
      REAL, INTENT(IN) :: dt ! in seconds
      REAL, INTENT(IN) :: dx

      REAL, INTENT(IN) :: minimum_profile_scaling
      REAL, INTENT(IN) :: maximum_profile_scaling
      REAL, INTENT(IN) :: search_range ![meters]

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: microphysics_latent_heating !(K/s)
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: cumulus_latent_heating !(K/s)

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: precip_accum_obs ! mm/hr

      INTEGER, INTENT(IN) :: lhn_start ! Latent Heat Nudging start time in seconds from simulation start time
      INTEGER, INTENT(IN) :: lhn_end ! Latent Heat Nudging end time in seconds from simulation start time

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(OUT) :: latent_heat_tendency

      INTEGER, INTENT(IN) :: ids
      INTEGER, INTENT(IN) :: ide
      INTEGER, INTENT(IN) :: jds
      INTEGER, INTENT(IN) :: jde
      INTEGER, INTENT(IN) :: kds
      INTEGER, INTENT(IN) :: kde
      INTEGER, INTENT(IN) :: ims
      INTEGER, INTENT(IN) :: ime
      INTEGER, INTENT(IN) :: jms
      INTEGER, INTENT(IN) :: jme
      INTEGER, INTENT(IN) :: kms
      INTEGER, INTENT(IN) :: kme
      INTEGER, INTENT(IN) :: num_tiles
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: i_start
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: i_end
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: j_start
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: j_end
      INTEGER, INTENT(IN) :: k_start
      INTEGER, INTENT(IN) :: k_end
      ! -------------------End Subroutine Parameters Declarations ---------------------------------

      ! -------------------- Local variables Declarations -------------------------------------------
      REAL, DIMENSION(ims:ime, jms:jme) :: model_precip_accum
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme) :: model_latent_heating !(K/s)
      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      INTEGER :: its
      INTEGER :: ite
      INTEGER :: jts
      INTEGER :: jte
      INTEGER :: kts
      INTEGER :: kte
      INTEGER :: current_profile_num
      INTEGER :: ii
      INTEGER :: jj
      INTEGER :: tile
      REAL :: latent_heat_factor
      LOGICAL :: apply_latent_heat_nudging
      INTEGER :: max_number_of_grid_points
      INTEGER :: profile_found
      REAL :: minimum_distance
      INTEGER :: i_profile
      INTEGER :: j_profile

      ! -------------------- End of Local variables Declarations ------------------------------------

      minimum_distance = search_range/dx ! minimum distance in grind points length

      apply_latent_heat_nudging = (lhn_start <= current_time_in_seconds) .AND. (current_time_in_seconds <= lhn_end)

      kts = k_start

      kte = min(k_end, kde - 1)

      max_number_of_grid_points = (search_range/dx)**2

      IF (apply_latent_heat_nudging) THEN
         model_latent_heating(:, :, :) = 0.
         write (*, *) "Applying latent heat nudging. t=", current_time_in_seconds, "seconds"

         ! First compute the model rainrate and check which model grid points meet the desired criteria (accepted profiles)


!$OMP PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!$OMP PRIVATE ( tile ,its,ite,jts,jte, i,j,i_profile, j_profile,profile_found)
         DO tile = 1, num_tiles
            its = i_start(tile)
            ite = i_end(tile)
            jts = j_start(tile)
            jte = j_end(tile)

            DO i = its, ite
            DO j = jts, jte

               profile_found = NO_PROFILE_FOUND

               CALL spiral_profile_search(i, j, max_number_of_grid_points, &
                                          precip_accum_obs, model_precip_accum, &
                                          minimum_profile_scaling, maximum_profile_scaling, &
                                          ids, ide, jds, jde, &
                                          ims, ime, jms, jme, &
                                          i_profile, j_profile, profile_found)

               IF (profile_found == PROFILE_CRITERIA_MEET) THEN

                  latent_heat_factor = (precip_accum_obs(i, j)/model_precip_accum(i, j)) - 1.0

                  DO k = kts, kte
                     latent_heat_tendency(i, k, j) = (microphysics_latent_heating(i, k, j) + &
                                                      cumulus_latent_heating(i, k, j))*latent_heat_factor
                  ENDDO

               ELSE IF (profile_found == NEARBY_PROFILE_FOUND) THEN

                  latent_heat_factor = (precip_accum_obs(i, j)/model_precip_accum(i_profile, j_profile))

                  DO k = kts, kte
                     latent_heat_tendency(i, k, j) = ((microphysics_latent_heating(i_profile, k, j_profile) + &
                                                       cumulus_latent_heating(i_profile, k, j_profile))*latent_heat_factor - &
                                                      (microphysics_latent_heating(i, k, j) + &
                                                       cumulus_latent_heating(i, k, j)))

                     latent_heat_tendency(i, k, j) = 0.

                  ENDDO

               ELSE

                  IF (abs(precip_accum_obs(i, j) - model_precip_accum(i, j)) < 0.001) THEN
                     DO k = kts, kte
                        latent_heat_tendency(i, k, j) = 0.

                     ENDDO
                  ELSE

                     ! No profile that match the criteria found, search unsuccessful
                     IF (precip_accum_obs(i, j) > model_precip_accum(i, j)) THEN
                        latent_heat_factor = maximum_profile_scaling - 1.
                     ELSE
                        latent_heat_factor = minimum_profile_scaling - 1.
                     ENDIF

                     DO k = kts, kte
                        latent_heat_tendency(i, k, j) = (microphysics_latent_heating(i, k, j) + &
                                                         cumulus_latent_heating(i, k, j))*latent_heat_factor
                     ENDDO
                  ENDIF

               ENDIF

            ENDDO
            ENDDO

         ENDDO
!$OMP END PARALLEL DO

         write (*, *)
         write (*, *) "minimum_profile_scaling", minimum_profile_scaling
         write (*, *) "maximum_profile_scaling", maximum_profile_scaling
         write (*, *)
         write (*, *) "Min latent_heat_tendency:", MINVAL(latent_heat_tendency)
         write (*, *) "Min cumulus_latent_heating CU:", MINVAL(cumulus_latent_heating)
         write (*, *) "Min microphysics_latent_heating CU:", MINVAL(microphysics_latent_heating)
         write (*, *)
         write (*, *) "Max latent_heat_tendency:", MAXVAL(latent_heat_tendency)
         write (*, *) "Max cumulus_latent_heating CU:", MAXVAL(cumulus_latent_heating)
         write (*, *) "Max microphysics_latent_heating CU:", MAXVAL(microphysics_latent_heating)
         write (*, *)
         write (*, *) "Min model_precip_accum:", MINVAL(model_precip_accum)
         write (*, *) "Min precip_accum_obs:", MINVAL(precip_accum_obs)
         write (*, *)
         write (*, *) "Max model_precip_accum:", MAXVAL(model_precip_accum)
         write (*, *) "Max precip_accum_obs:", MAXVAL(precip_accum_obs)

      ELSE ! Do not apply latent heat nudging

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!$OMP PRIVATE ( tile ,its,ite,jts,jte, i,j,k)
         DO tile = 1, num_tiles
            its = i_start(tile)
            ite = i_end(tile)
            jts = j_start(tile)
            jte = j_end(tile)

            DO i = its, ite
            DO k = kts, kte
            DO j = jts, jte
               latent_heat_tendency(i, k, j) = 0.
            ENDDO
            ENDDO
            ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

   END SUBROUTINE latent_heat_nudging_1

! https://stackoverflow.com/questions/3706219/algorithm-for-iterating-over-an-outward-spiral-on-a-discrete-2d-grid-from-the-or
   SUBROUTINE spiral_profile_search(i_center, j_center, max_number_of_grid_points, &
                                    precip_accum_obs, model_precip_accum, &
                                    minimum_profile_scaling, maximum_profile_scaling, &
                                    ids, ide, jds, jde, &
                                    ims, ime, jms, jme, &
                                    i_profile, j_profile, profile_found)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i_center
      INTEGER, INTENT(IN) :: j_center
      INTEGER, INTENT(IN) :: max_number_of_grid_points

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: model_precip_accum
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN):: precip_accum_obs !(K/s)
      REAL, INTENT(IN) :: minimum_profile_scaling
      REAL, INTENT(IN) :: maximum_profile_scaling
      INTEGER, INTENT(IN) :: ids
      INTEGER, INTENT(IN) :: ide
      INTEGER, INTENT(IN) :: jds
      INTEGER, INTENT(IN) :: jde
      INTEGER, INTENT(IN) :: ims
      INTEGER, INTENT(IN) :: ime
      INTEGER, INTENT(IN) :: jms
      INTEGER, INTENT(IN) :: jme

      INTEGER, INTENT(OUT) :: i_profile
      INTEGER, INTENT(OUT) :: j_profile
      INTEGER, INTENT(OUT) :: profile_found

      ! Local Declarations
      INTEGER :: di
      INTEGER :: dj
      INTEGER :: buffer
      INTEGER :: number_of_steps
      INTEGER :: steps
      INTEGER :: segment_points_visited
      INTEGER :: segment_length

      di = 1
      dj = 0
      number_of_steps = 1

      i_profile = i_center
      j_profile = j_center

      profile_found = NO_PROFILE_FOUND

      ! Check if the target grid point meet the criteria
      if ((model_precip_accum(i_center, j_center) > 0.001) .AND. & ! Avoid divisions close to zero
          (precip_accum_obs(i_center, j_center) < model_precip_accum(i_center, j_center)*maximum_profile_scaling) .AND. &
          (precip_accum_obs(i_center, j_center) > model_precip_accum(i_center, j_center)*minimum_profile_scaling)) THEN

         profile_found = PROFILE_CRITERIA_MEET
         RETURN
      ENDIF

      ! Otherwise:  search for nearby profiles

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

         if ((model_precip_accum(i_profile, j_profile) > 0.001) .AND. & ! Avoid divisions close to zero
             (precip_accum_obs(i_center, j_center) < model_precip_accum(i_profile, j_profile)*maximum_profile_scaling) .AND. &
             (precip_accum_obs(i_center, j_center) > model_precip_accum(i_profile, j_profile)*minimum_profile_scaling)) THEN

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

   ! Update the model rain accumulation over an interval
   SUBROUTINE update_accumulation(RAINCV,RAINSHV,RAINNCV, &
                                  model_previous_timestep_precip_accum, &
                                  model_precip_accum, &
                                  accum_time_levels, &
                                  ids, ide, jds, jde, &
                                  ims, ime, jms, jme, &
                                  num_tiles, i_start, i_end, j_start, j_end)
      IMPLICIT NONE

      !------------------------------ Subroutine Parameters Declarations --------------------------

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: RAINCV
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: RAINSHV
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: RAINNCV

      REAL, DIMENSION(1:accum_time_levels, ims:ime, jms:jme), INTENT(INOUT) :: model_previous_timestep_precip_accum

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(OUT) :: model_precip_accum

      INTEGER, INTENT(IN) :: accum_time_levels
      INTEGER, INTENT(IN) :: ids
      INTEGER, INTENT(IN) :: ide
      INTEGER, INTENT(IN) :: jds
      INTEGER, INTENT(IN) :: jde
      INTEGER, INTENT(IN) :: ims
      INTEGER, INTENT(IN) :: ime
      INTEGER, INTENT(IN) :: jms
      INTEGER, INTENT(IN) :: jme

      INTEGER, INTENT(IN) :: num_tiles
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: i_start
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: i_end
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: j_start
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: j_end

      !--------------------------------------------------------------------------------------------
      ! Local declarations
      REAL, DIMENSION(ims:ime, jms:jme) :: model_timestep_precip_accum
      INTEGER :: t
      INTEGER :: i
      INTEGER :: j
      INTEGER :: its
      INTEGER :: ite
      INTEGER :: jts
      INTEGER :: jte
      INTEGER :: tile


!$OMP PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!$OMP PRIVATE ( tile ,its,ite,jts,jte, i,j)
         DO tile = 1, num_tiles
            its = i_start(tile)
            ite = i_end(tile)
            jts = j_start(tile)
            jte = j_end(tile)

            DO i = its, ite
            DO j = jts, jte
               model_timestep_precip_accum(i, j) = RAINCV(i, j) + RAINSHV(i, j) + RAINNCV(i, j)
            ENDDO
            ENDDO
         ENDDO
!$OMP END PARALLEL DO

      IF (accum_time_levels == 1) THEN

!$OMP     PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!$OMP PRIVATE ( tile ,its,ite,jts,jte, i,j)
         DO tile = 1, num_tiles
            its = i_start(tile)
            ite = i_end(tile)
            jts = j_start(tile)
            jte = j_end(tile)

            DO i = its, ite
            DO j = jts, jte

               model_previous_timestep_precip_accum(1, i, j) = model_timestep_precip_accum(i, j)
               model_precip_accum(i, j) = model_previous_timestep_precip_accum(1, i, j)
            ENDDO
            ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSE !accum_time_levels > 1

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!$OMP PRIVATE ( tile ,its,ite,jts,jte, i,j)
         DO tile = 1, num_tiles
            its = i_start(tile)
            ite = i_end(tile)
            jts = j_start(tile)
            jte = j_end(tile)

            DO i = its, ite
            DO j = jts, jte
               ! Update previous rainrate
               model_precip_accum(i, j) = 0.

               DO t = 1, accum_time_levels - 1
                  model_previous_timestep_precip_accum(t, i, j) = model_previous_timestep_precip_accum(t + 1, i, j)

                  model_precip_accum(i, j) = model_precip_accum(i, j) + model_previous_timestep_precip_accum(t, i, j)
               ENDDO

               model_previous_timestep_precip_accum(accum_time_levels, i, j) = model_timestep_precip_accum(i, j)

               model_precip_accum(i, j) = model_precip_accum(i, j) + model_timestep_precip_accum(i, j)
            ENDDO
            ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

   END SUBROUTINE update_accumulation

END MODULE module_latent_heat_nudging
