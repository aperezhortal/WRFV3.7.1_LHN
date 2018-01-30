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
   SUBROUTINE latent_heat_nudging_JM97(current_time_in_seconds, dt, dx, &
                                       microphysics_latent_heating, &
                                       cumulus_latent_heating, &
                                       raincv, rainshv, rainncv, rainrate_observations, &
                                       lhn_start, lhn_end, latent_heat_tendency, &
                                       minimum_profile_scaling, maximum_profile_scaling, &
                                       search_range, &
                                       lhn_smooth_filter_sigma, &
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

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: raincv !time-step cumulus precipitation in mm
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: rainshv !time-step shallow cumulus precipitation in mm
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: rainncv !time-step non convective precipitation in mm
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: rainrate_observations ! mm/hr

      INTEGER, INTENT(IN) :: lhn_start ! LHN start time in seconds from simulation start time
      INTEGER, INTENT(IN) :: lhn_end ! LHN end time in seconds from simulation start time
      REAL, INTENT(IN) :: lhn_smooth_filter_sigma ! LHN horizontal recursive filter standard deviation in grid points

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
      INTEGER :: kte
      INTEGER :: current_profile_num
      INTEGER :: ii
      INTEGER :: jj
      INTEGER :: tile
      REAL :: latent_heat_factor
      LOGICAL :: apply_latent_heat_nudging
      INTEGER :: max_number_of_grid_points
      INTEGER :: profile_found
      INTEGER :: i_profile
      INTEGER :: j_profile

      ! -------------------- End of Local variables Declarations ------------------------------------



      apply_latent_heat_nudging = (lhn_start <= current_time_in_seconds) .AND. (current_time_in_seconds <= lhn_end)

      kts = k_start

      kte = min(k_end, kde - 1)

      ! Just in case that a rounding or trucation error is present
      IF (search_range > 1.) THEN ! Search range bigger tha 1, otherwise is considered as zero
         max_number_of_grid_points = (search_range/dx)**2
      ELSE
         max_number_of_grid_points = 0.
      ENDIF

      IF (apply_latent_heat_nudging) THEN

         model_latent_heating(:, :, :) = 0.

         write (*, *) "Applying latent heat nudging. t=", current_time_in_seconds, "seconds"

         ! First compute the model rainrate

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!$OMP PRIVATE ( tile ,its,ite,jts,jte, i,j)
         DO tile = 1, num_tiles
            its = i_start(tile)
            ite = i_end(tile)
            jts = j_start(tile)
            jte = j_end(tile)

            DO i = its, ite
            DO j = jts, jte
               model_rainrate(i, j) = (RAINCV(i, j) + RAINSHV(i, j) + RAINNCV(i, j))*3600./dt ! in mm/hr
            ENDDO
            ENDDO
         ENDDO
!$OMP END PARALLEL DO


         ! Check which model grid points meet the desired criteria
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!$OMP PRIVATE ( tile ,its,ite,jts,jte, i,j,i_profile, j_profile,profile_found)
         DO tile = 1, num_tiles
            its = i_start(tile)
            ite = i_end(tile)
            jts = j_start(tile)
            jte = j_end(tile)

            DO i = its, ite
            DO j = jts, jte


               !---------------------------------------------------------------------------------------------
               ! Check if the target grid point meet the criteria
               IF ((model_rainrate(i, j) > 0.001) .AND. & ! Avoid divisions close to zero
                   (rainrate_observations(i, j) < model_rainrate(i, j)*maximum_profile_scaling) .AND. &
                   (rainrate_observations(i, j) > model_rainrate(i, j)*minimum_profile_scaling)) THEN

                  profile_found = PROFILE_CRITERIA_MEET
!                  write(*,*) "Found!"

               ELSE IF (max_number_of_grid_points > .0) THEN
                  ! IF criteria is not met, perform an spiral search in the nearby profiles

                  CALL spiral_profile_search(i, j, max_number_of_grid_points, &
                                             rainrate_observations, model_rainrate, &
                                             minimum_profile_scaling, maximum_profile_scaling, &
                                             ids, ide, jds, jde, &
                                             ims, ime, jms, jme, &
                                             i_profile, j_profile, profile_found)

               ELSE

                  profile_found = NO_PROFILE_FOUND

               ENDIF
               !---------------------------------------------------------------------------------------------



               !---------------------------------------------------------------------------------------------
               ! Scale the profile according to the criteria that the profiles matches
               IF (profile_found == PROFILE_CRITERIA_MEET) THEN

                  latent_heat_factor = (rainrate_observations(i, j)/model_rainrate(i, j)) - 1.0

                  DO k = kts, kte
                     latent_heat_tendency(i, k, j) = (microphysics_latent_heating(i, k, j) + &
                                                      cumulus_latent_heating(i, k, j))*latent_heat_factor
                  ENDDO

               ELSE IF (profile_found == NEARBY_PROFILE_FOUND) THEN
                  !write(*,*) 'NEARBY_PROFILE_FOUND','i=',i,'j=',j,'i_profile',i_profile,'j_profile',j_profile
                  latent_heat_factor = (rainrate_observations(i, j)/model_rainrate(i_profile, j_profile))

                  DO k = kts, kte
                     latent_heat_tendency(i, k, j) = ((microphysics_latent_heating(i_profile, k, j_profile) + &
                                                       cumulus_latent_heating(i_profile, k, j_profile))*latent_heat_factor - &
                                                      (microphysics_latent_heating(i, k, j) + &
                                                       cumulus_latent_heating(i, k, j)))

                  ENDDO

               ELSE !profile_found == NO_PROFILE_FOUND

                  IF (abs(rainrate_observations(i, j) - model_rainrate(i, j)) < 0.001) THEN
                     DO k = kts, kte
                        latent_heat_tendency(i, k, j) = 0.
                     ENDDO
                  ELSE

                     ! No profile that match the criteria found, search unsuccessful
                     IF (rainrate_observations(i, j) > model_rainrate(i, j)) THEN
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
               !---------------------------------------------------------------------------------------------

            ENDDO
            ENDDO

         ENDDO
!$OMP END PARALLEL DO

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
         write (*, *) "Min model_rainrate:", MINVAL(model_rainrate(ids:ide, jds:jde))
         write (*, *) "Min rainrate_observations:", MINVAL(rainrate_observations(ids:ide, jds:jde))
         write (*, *)
         write (*, *) "Max model_rainrate:", MAXVAL(model_rainrate(ids:ide, jds:jde))
         write (*, *) "Max rainrate_observations:", MAXVAL(rainrate_observations(ids:ide, jds:jde))

         IF (lhn_smooth_filter_sigma > 0.0001) THEN
            CALL recursive_filter_2d(latent_heat_tendency, &
                                     lhn_smooth_filter_sigma, &
                                     ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme)
         ENDIF
      ELSE ! Do not apply latent heat nudging

         latent_heat_tendency(:, :, :) = 0.

      ENDIF

   END SUBROUTINE latent_heat_nudging_JM97





! https://stackoverflow.com/questions/3706219/algorithm-for-iterating-over-an-outward-spiral-on-a-discrete-2d-grid-from-the-or
   SUBROUTINE spiral_profile_search(i_center, j_center, max_number_of_grid_points, &
                                    rainrate_observations, model_rainrate, &
                                    minimum_profile_scaling, maximum_profile_scaling, &
                                    ids, ide, jds, jde, &
                                    ims, ime, jms, jme, &
                                    i_profile, j_profile, profile_found)

      INTEGER, INTENT(IN) :: i_center
      INTEGER, INTENT(IN) :: j_center
      INTEGER, INTENT(IN) :: max_number_of_grid_points

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: model_rainrate
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN):: rainrate_observations !(K/s)
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

         if ((model_rainrate(i_profile, j_profile) > 0.001) .AND. & ! Avoid divisions close to zero
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




   SUBROUTINE recursive_filter_2d(input_3d, &
                                  sigma_in_gridpoints, &
                                  ids, ide, jds, jde, kds, kde, &
                                  ims, ime, jms, jme, kms, kme)

      IMPLICIT NONE

      !------------------------------ Subroutine Parameters Declarations --------------------------
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT) :: input_3d
      REAL, INTENT(IN) :: sigma_in_gridpoints

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
      ! -------------------End Subroutine Parameters Declarations ---------------------------------

      ! -------------------- Local variables Declarations -------------------------------------------
      INTEGER :: i
      INTEGER :: j
      INTEGER :: k

      REAL, DIMENSION(ims:ime, kms:kme, jms:jme) :: forward_filter
      REAL, DIMENSION(ims:ime, kms:kme, jms:jme) :: backward_filter
      REAL :: E
      REAL :: a_coef
      REAL :: b_coef
      ! -------------------- End of Local variables Declarations ------------------------------------
      write(*,*) "Smoothing"
      E = 1.0/(sigma_in_gridpoints*sigma_in_gridpoints)
      a_coef = 1.0+E - sqrt(E*(E + 2.0))
      b_coef = 1.0-a_coef

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!$OMP PRIVATE ( k )
      DO k = kds, kde
         DO j = jds, jde
            ! Apply forward filter "x" direction first
            forward_filter(ids, k, j) = b_coef*input_3d(ids, k, j)
            DO i = ids + 1, ide
               forward_filter(i, k, j) = b_coef*input_3d(i, k, j) + a_coef*forward_filter(i - 1, k, j)
            ENDDO

            ! Apply backward filter "x" direction
            backward_filter(ide, k, j) = b_coef*input_3d(ide, k, j)
            DO i = ide - 1, ids, -1
               backward_filter(i, k, j) = b_coef*forward_filter(i, k, j) + a_coef*backward_filter(i + 1, k, j)
            ENDDO

         ENDDO

         DO i = ids, ide
            ! Apply forward filter y direction
            forward_filter(i, k, jds) = b_coef*input_3d(i, k, jds)
            DO j = jds + 1, jde
               forward_filter(i, k, j) = b_coef*input_3d(i, k, j) + a_coef*forward_filter(i, k, j - 1)
            ENDDO

            ! Apply backward filter "y" direction
            backward_filter(i, k, jde) = b_coef*forward_filter(i, k, jde)
            DO j = jde - 1, jds, -1
               backward_filter(i, k, j) = b_coef*forward_filter(i, k, j) + a_coef*backward_filter(i, k, j + 1)
            ENDDO

         ENDDO

      ENDDO
!$OMP END PARALLEL DO

      input_3d = backward_filter

   END SUBROUTINE recursive_filter_2d

END MODULE module_latent_heat_nudging
