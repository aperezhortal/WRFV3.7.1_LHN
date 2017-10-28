!WRF:MEDIATION_LAYER:PHYSICS
!

MODULE module_latent_heat_nudging

   USE kdtree2_module, ONLY:kdtree2_create, kdtree2_n_nearest, kdtree2_result, kdtree2_destroy, kdtree2

   IMPLICIT NONE

CONTAINS

   ! Jones, C. D., & Macpherson, B. (1997). A latent heat nudging scheme for the assimilation of precipitation data
   ! into an operational mesoscale model. Meteorological Applications, 4(3), 269-277.
   SUBROUTINE latent_heat_nudging_1(current_time_in_seconds, dt, dx, &
                                    microphysics_latent_heating, &
                                    cumulus_latent_heating, &
                                    raincv, rainshv, rainncv, rainrate_observations, &
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

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: raincv !time-step cumulus precipitation in mm
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: rainshv !time-step shallow cumulus precipitation in mm
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: rainncv !time-step non convective precipitation in mm
      REAL, DIMENSION(ims:ime, jms:jme), INTENT(IN) :: rainrate_observations ! mm/hr

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
      LOGICAL, DIMENSION(ims:ime, jms:jme) :: criteria_meet
      LOGICAL :: apply_latent_heat_nudging
      LOGICAL :: profile_seach_flag
      INTEGER :: number_of_accepted_profiles
      INTEGER :: number_of_accepted_profiles_in_tile
      REAL :: minimum_distance

      REAL, ALLOCATABLE :: accepted_profiles_position(:, :)
      REAL :: grid_point_position(2)
      type(kdtree2), pointer :: my_kdtree
      TYPE(kdtree2_result) :: nearest_profile(1)

      ! -------------------- End of Local variables Declarations ------------------------------------

      minimum_distance = search_range/dx ! minimum distance in grind points length

      apply_latent_heat_nudging = (lhn_start <= current_time_in_seconds) .AND. (current_time_in_seconds <= lhn_end)

      kts = k_start

      kte = min(k_end, kde - 1)

      model_rainrate(:,:)=0.
      model_latent_heating(:,:,:)=0.
      !latent_heat_tendency(:,:,:)=0.
      number_of_accepted_profiles = 0


      IF (apply_latent_heat_nudging) THEN

         write(*,*) "Applying latent heat nudging. t=",current_time_in_seconds , "seconds"

         ! First compute the model rainrate and check which model grid points meet the desired criteria (accepted profiles)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!$OMP PRIVATE ( tile ,its,ite,jts,jte, i,j,number_of_accepted_profiles_in_tile) &
!$OMP REDUCTION(+:number_of_accepted_profiles)
         DO tile = 1, num_tiles
            its = i_start(tile)
            ite = i_end(tile)
            jts = j_start(tile)
            jte = j_end(tile)
            number_of_accepted_profiles_in_tile = 0

            DO i = its, ite
            DO j = jts, jte

               model_rainrate(i, j) = (RAINCV(i, j) + RAINSHV(i, j) + RAINNCV(i, j))*3600./dt ! in mm/hr

               if ( model_rainrate(i, j) > 0.01 )  then ! Avoid divisions close to zero

                 criteria_meet(i, j) = (rainrate_observations(i, j) < model_rainrate(i, j)*maximum_profile_scaling) .AND. &
                                       (rainrate_observations(i, j) > model_rainrate(i, j)*minimum_profile_scaling)

               else
                  criteria_meet(i, j) = .False.
               endif


               IF (criteria_meet(i, j)) THEN
                  number_of_accepted_profiles_in_tile = number_of_accepted_profiles_in_tile + 1
               ENDIF

            ENDDO
            ENDDO
            number_of_accepted_profiles = number_of_accepted_profiles + number_of_accepted_profiles_in_tile

         ENDDO
!$OMP END PARALLEL DO


         ! If there are only a few profiles that meet the criteria do not perform the profile search
         IF (number_of_accepted_profiles < 4) THEN
            profile_seach_flag = .FALSE.
         ELSE
            profile_seach_flag = .TRUE.
           ! Second: Create the kdtree with the accepted profiles
             ALLOCATE (accepted_profiles_position(2, number_of_accepted_profiles))
             current_profile_num = 1

             DO i = ids, ide
             DO j = jds, jde

                IF (criteria_meet(i, j)) THEN
                   accepted_profiles_position(1, current_profile_num) = real(i)
                   accepted_profiles_position(2, current_profile_num) = real(j)
                   current_profile_num = current_profile_num + 1
                ENDIF

             ENDDO
             ENDDO

             my_kdtree => kdtree2_create(accepted_profiles_position, 2, sort=.TRUE., rearrange=.TRUE.)

         ENDIF


!$OMP PARALLEL DO SCHEDULE(DYNAMIC) default(shared) &
!$OMP PRIVATE ( tile ,its,ite,jts,jte, i,j,k,ii,jj,latent_heat_factor,nearest_profile,grid_point_position)
         DO tile = 1, num_tiles
            its = i_start(tile)
            ite = i_end(tile)
            jts = j_start(tile)
            jte = j_end(tile)
            ! Third: Compute the latent heat tendencies

            DO i = its, ite
            DO j = jts, jte

               IF (criteria_meet(i, j)) THEN
                  ! If criteria is meet, apply the correction factor

                  !write(*,*) "rainrate_observations(i, j)",rainrate_observations(i, j)
                  !write(*,*) "model_rainrate(i, j)",model_rainrate(i, j)


                  latent_heat_factor = (rainrate_observations(i, j)/model_rainrate(i, j)) - 1.0


                  DO k = kts, kte

                     latent_heat_tendency(i, k, j) = (microphysics_latent_heating(i, k, j) + &
                                                      cumulus_latent_heating(i, k, j))*latent_heat_factor
                  ENDDO

               ELSE
                  ! If Criteria is not meet, perform a search in nearby points for a profile that meet the criteria

                  IF (profile_seach_flag) THEN
                  !IF (.FALSE.) THEN

                     grid_point_position(1) = i
                     grid_point_position(2) = j

                     ! If criteria is not meet, search profiles nearby that meet the .criteria
!$OMP CRITICAL
                    ! This should be executed by one single thread at a time, otherwise, an error occurs who knows why...
                    CALL kdtree2_n_nearest(tp=my_kdtree, qv=grid_point_position, nn=1, results=nearest_profile)
!$OMP END CRITICAL
                     !write(*,*) "i",i,"j",j
                     !write(*,*) "nearest_profile(1)%idx",nearest_profile(1)%idx
                     !write(*,*) "nearest_profile(1)%dis",nearest_profile(1)%dis


                     ! Found profile position indexes = (ii,jj)
                     ii = nint(accepted_profiles_position(1, nearest_profile(1)%idx))
                     jj = nint(accepted_profiles_position(2, nearest_profile(1)%idx))

                     IF (sqrt(nearest_profile(1)%dis) < minimum_distance) THEN

                        latent_heat_factor = (rainrate_observations(i, j)/model_rainrate(ii, jj))

                        DO k = kts, kte
                           latent_heat_tendency(i, k, j) = ((microphysics_latent_heating(ii, k, jj) + &
                                                             cumulus_latent_heating(ii, k, jj))*latent_heat_factor - &
                                                            (microphysics_latent_heating(i, k, j) + &
                                                             cumulus_latent_heating(i, k, j)))

                        ENDDO

                     ELSE
                        ! The profile is far away, search unsuccessful
                        IF (rainrate_observations(i, j) > model_rainrate(i, j)) THEN
                           latent_heat_factor = maximum_profile_scaling - 1.
                        ELSE
                           latent_heat_factor = minimum_profile_scaling - 1.
                        ENDIF

                        IF (abs(rainrate_observations(i, j)-model_rainrate(i, j)) > 0.01) THEN
                            latent_heat_factor = 0
                        ELSE
                            write(*,*)
                            write(*,*) "rainrate_observations(i, j)", rainrate_observations(i, j)
                            write(*,*) "model_rainrate(i, j)", model_rainrate(i, j)

                        ENDIF

                        DO k = kts, kte
                           latent_heat_tendency(i, k, j) = (microphysics_latent_heating(i, k, j) + &
                                                            cumulus_latent_heating(i, k, j))*latent_heat_factor

                        ENDDO

                     ENDIF
                  ELSE
                     DO k = kts, kte
                        latent_heat_tendency(i, k, j) = 0.
                     ENDDO
                  ENDIF

               ENDIF

            ENDDO
            ENDDO

         ENDDO
!$OMP END PARALLEL DO

        IF (profile_seach_flag) THEN
         DEALLOCATE (accepted_profiles_position)
         CALL kdtree2_destroy(my_kdtree)
        ENDIF

   write(*,*)
   write(*,*) "minimum_profile_scaling",minimum_profile_scaling
   write(*,*) "maximum_profile_scaling",maximum_profile_scaling
   write(*,*)
   write(*,*) "Min latent_heat_tendency:",MINVAL(latent_heat_tendency)
   write(*,*) "Min cumulus_latent_heating CU:",MINVAL(cumulus_latent_heating)
   write(*,*) "Min microphysics_latent_heating CU:",MINVAL(microphysics_latent_heating)
   write(*,*)
   write(*,*) "Max latent_heat_tendency:",MAXVAL(latent_heat_tendency)
   write(*,*) "Max cumulus_latent_heating CU:",MAXVAL(cumulus_latent_heating)
   write(*,*) "Max microphysics_latent_heating CU:",MAXVAL(microphysics_latent_heating)
   write(*,*)
   write(*,*) "Min model_rainrate:",MINVAL(model_rainrate)
   write(*,*) "Min rainrate_observations:",MINVAL(rainrate_observations)
   write(*,*)
   write(*,*) "Max model_rainrate:",MAXVAL(model_rainrate)
   write(*,*) "Max rainrate_observations:",MAXVAL(rainrate_observations)


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

END MODULE module_latent_heat_nudging
