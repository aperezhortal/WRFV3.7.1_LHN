
!
MODULE module_analysis_nudging

CONTAINS

!
!-------------------------------------------------------------------
!
   SUBROUTINE nudge_towards_analysis(itimestep, dt, xtime, &
                                     analysis_time, relaxation_time,nudging_gain, &
                                     nudge_u, nudge_v, nudge_w, nudge_q, nudge_t, nudge_ph, nudge_mu, &
                                     u_analysis, v_analysis, t_analysis, q_analysis, w_analysis, &
                                     ph_analysis, mu_analysis, &
                                     nudging_mask, &
                                     u_staggered, v_staggered, w_staggered, qvapor, temperature, &
                                     ph_staggered, mu, &
                                     rundgdten, rvndgdten, rwndgdten, rqvndgdten, rthndgdten, &
                                     rphndgdten, rmundgdten, &
                                     nca,cu_reinit, cu_reinit_dt, cu_reinit_duration,cu_reinit_counter, &
                                     ids, ide, jds, jde, kds, kde, &
                                     ims, ime, jms, jme, kms, kme, &
                                     i_start, i_end, j_start, j_end, kts, kte, num_tiles)
!-------------------------------------------------------------------
      IMPLICIT NONE
      !-------------------------------------------------------------------

      !-- xtime time in minutes from the start of the simulation
      !-- u3d 3d u-velocity staggered on u points
      !-- v3d 3d v-velocity staggered on v points
      !-- th3d 3d potential temperature (k)
      !-- t3d temperature (k)
      !-- qv3d 3d water vapor mixing ratio (kg/kg)
      !-- p3d 3d pressure (pa)
      !-- pi3d 3d exner function (dimensionless)
      !-- rundgdten staggered u tendency due to
      ! fdda grid nudging (m/s/s)
      !-- rvndgdten staggered v tendency due to
      ! fdda grid nudging (m/s/s)
      !-- rthndgdten theta tendency due to
      ! fdda grid nudging (K/s)
      !-- rqvndgdten qv tendency due to
      ! fdda grid nudging (kg/kg/s)
      !-- rmundgdten mu tendency due to
      ! fdda grid nudging (Pa/s)
      !-- ids start index for i in domain
      !-- ide end index for i in domain
      !-- jds start index for j in domain
      !-- jde end index for j in domain
      !-- kds start index for k in domain
      !-- kde end index for k in domain
      !-- ims start index for i in memory
      !-- ime end index for i in memory
      !-- jms start index for j in memory
      !-- jme end index for j in memory
      !-- kms start index for k in memory
      !-- kme end index for k in memory
      !-- its start index for i in tile
      !-- ite end index for i in tile
      !-- jts start index for j in tile
      !-- jte end index for j in tile
      !-- kts start index for k in tile
      !-- kte end index for k in tile
      !-------------------------------------------------------------------
      !

      INTEGER, INTENT(IN) :: itimestep
      INTEGER, INTENT(IN) :: analysis_time
      INTEGER, INTENT(IN) :: num_tiles
      REAL, INTENT(IN) :: dt
      REAL, INTENT(IN) :: xtime
      REAL, INTENT(IN) :: relaxation_time
      REAL, INTENT(IN) :: nudging_gain

      LOGICAL, INTENT(IN) :: nudge_u
      LOGICAL, INTENT(IN) :: nudge_v
      LOGICAL, INTENT(IN) :: nudge_w
      LOGICAL, INTENT(IN) :: nudge_q
      LOGICAL, INTENT(IN) :: nudge_t
      LOGICAL, INTENT(IN) :: nudge_ph
      LOGICAL, INTENT(IN) :: nudge_mu

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
      INTEGER, INTENT(IN) :: kts
      INTEGER, INTENT(IN) :: kte

      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: i_start
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: i_end
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: j_start
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: j_end

      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: u_staggered
      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: v_staggered
      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: w_staggered
      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: qvapor
      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: temperature
      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: ph_staggered
      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: u_analysis
      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: v_analysis
      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: t_analysis
      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: q_analysis
      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: w_analysis
      REAL, INTENT(IN), DIMENSION(ims:ime, kms:kme, jms:jme) :: ph_analysis
      REAL, INTENT(INOUT), DIMENSION( ims:ime , jms:jme )    :: nca

      LOGICAL, INTENT(IN) :: cu_reinit
      REAL, INTENT(IN) :: cu_reinit_dt
      REAL, INTENT(IN) :: cu_reinit_duration
      REAL, INTENT(INOUT) :: cu_reinit_counter

      REAL, INTENT(IN), DIMENSION(ims:ime, jms:jme) :: mu
      REAL, INTENT(INOUT), DIMENSION(ims:ime, jms:jme) :: mu_analysis
      REAL, INTENT(INOUT), DIMENSION(ims:ime, jms:jme) :: nudging_mask

      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: rundgdten
      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: rvndgdten
      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: rwndgdten
      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: rqvndgdten
      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: rthndgdten
      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: rphndgdten
      REAL, INTENT(INOUT), DIMENSION(ims:ime, jms:jme) :: rmundgdten

      ! Local Variables
      INTEGER :: i
      INTEGER :: j
      INTEGER :: k
      INTEGER :: ij
      REAL :: nudgingstarttimeinseconds
      REAL :: nudgingendtimeinseconds
      REAL :: currenttimeinseconds
      REAL :: nudging_factor
      REAL :: cu_reinit_end

      CHARACTER(LEN=256):: message

      !DO k = k_start, k_end
      ! z_at_w(k) = (grid%phb(ix,k,iy)+grid%ph_2(ix,k,iy))/g
      !END DO

      !DO k = k_start, k_end-1
      ! z(k) = 0.5*(z_at_w(k) + z_at_w(k+1))
      !END DO


      nudgingstarttimeinseconds = analysis_time - relaxation_time

      nudgingendtimeinseconds = analysis_time

      currenttimeinseconds = xtime*60.

      IF (currenttimeinseconds .GE. nudgingstarttimeinseconds .AND. &
          currenttimeinseconds .LT. nudgingendtimeinseconds) THEN

          WRITE (*,*) 'Apllying nudging. currentTimeInSeconds :', currenttimeinseconds

          nudging_factor = nudging_gain / relaxation_time

!$OMP PARALLEL DO &
!$OMP PRIVATE ( ij,i,j,k )
         DO ij = 1, num_tiles

            IF (nudge_u) THEN
               DO j = j_start(ij), j_end(ij)
               DO k = kts, kte
               DO i = i_start(ij), i_end(ij) + 1 ! U staggered in X -> Need to add 1 to i_end
                  if (nudging_mask(i, j) > 0.) then
                     rundgdten(i, k, j) = (U_ANALYSIS(i, k, j) - u_staggered(i, k, j))*nudging_factor
                  else
                     rundgdten(i, k, j) = 0.
                  endif
               ENDDO
               ENDDO
               ENDDO
            ENDIF

            IF (nudge_v) THEN
               DO j = j_start(ij), j_end(ij) + 1 ! V staggered in Y -> Need to add 1 to j_end
               DO k = kts, kte
               DO i = i_start(ij), i_end(ij)

                  if (nudging_mask(i, j) > 0.) then
                     rvndgdten(i, k, j) = (V_ANALYSIS(i, k, j) - v_staggered(i, k, j))*nudging_factor
                  else
                     rvndgdten(i, k, j) = 0.
                  endif

               ENDDO
               ENDDO
               ENDDO
            ENDIF

            IF (nudge_w) THEN
               DO j = j_start(ij), j_end(ij)
               DO k = kts, kte + 1 ! W staggered in Z -> Need to add 1 to kte
               DO i = i_start(ij), i_end(ij)
                  if (nudging_mask(i, j) > 0.) then
                     rwndgdten(i, k, j) = (W_ANALYSIS(i, k, j) - W_staggered(i, k, j))*nudging_factor
                  else
                     rwndgdten(i, k, j) = 0.
                  endif
               ENDDO
               ENDDO
               ENDDO
            ENDIF

            IF (nudge_ph) THEN
               DO j = j_start(ij), j_end(ij)
               DO k = kts, kte + 1 ! W staggered in Z -> Need to add 1 to kte
               DO i = i_start(ij), i_end(ij)
                  if (nudging_mask(i, j) > 0.) then
                     rphndgdten(i, k, j) = (PH_ANALYSIS(i, k, j) - ph_staggered(i, k, j))*nudging_factor
                  else
                     rphndgdten(i, k, j) = 0.
                  endif
               ENDDO
               ENDDO
               ENDDO
            ENDIF

            IF (nudge_t) THEN

               DO j = j_start(ij), j_end(ij)
               DO k = kts, kte
               DO i = i_start(ij), i_end(ij)
                  if (nudging_mask(i, j) > 0.) then
                     rthndgdten(i, k, j) = (T_ANALYSIS(i, k, j) + 300.-temperature(i, k, j))*nudging_factor
                  else
                     rthndgdten(i, k, j) = 0.
                  endif
               ENDDO
               ENDDO
               ENDDO
            ENDIF

            IF (nudge_q) THEN
               DO j = j_start(ij), j_end(ij)
               DO k = kts, kte
               DO i = i_start(ij), i_end(ij)
                  if (nudging_mask(i, j) > 0.) then
                     rqvndgdten(i, k, j) = (q_ANALYSIS(i, k, j) - qvapor(i, k, j))*nudging_factor
                  else
                     rqvndgdten(i, k, j) = 0.
                  endif
               ENDDO
               ENDDO
               ENDDO
            ENDIF

            IF (nudge_mu) THEN
               DO j = j_start(ij), j_end(ij)
               DO i = i_start(ij), i_end(ij)
                  if (nudging_mask(i, j) > 0.) then
                     rmundgdten(i, j) = (MU_ANALYSIS(i, j) - mu(i, j))*nudging_factor
                  else
                     rmundgdten(i, j) = 0.
                  endif
               ENDDO
               ENDDO
            ENDIF

         ENDDO
!$OMP END PARALLEL DO

      ELSE

         CALL initialize_nudging_tendencies(rundgdten, rvndgdten, rwndgdten, rqvndgdten, rthndgdten, &
                                            rphndgdten, rmundgdten, &
                                            ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            i_start, i_end, j_start, j_end, kts, kte, num_tiles)

      ENDIF


      IF (cu_reinit) THEN
          IF (currenttimeinseconds .LT. analysis_time) THEN
                cu_reinit_counter = 0.
          ENDIF

          ! Reinit Cumulus parametrization periodically
          IF ( ABS(currenttimeinseconds - analysis_time) .lt. dt ) THEN
              cu_reinit_counter = cu_reinit_dt
              !Initialize at the analysis time (end of nudging period)
          ENDIF

          cu_reinit_end = analysis_time + cu_reinit_duration
          IF ( (currenttimeinseconds .GE. analysis_time) .AND. &
               (currenttimeinseconds .LE. cu_reinit_end) ) THEN

                IF ( cu_reinit_counter .ge. cu_reinit_dt) THEN
                    WRITE (6, *) 'Reinit CU counter'
                    cu_reinit_counter = 0.
                    nca(:,:)=-100.
                ENDIF
          ENDIF

          cu_reinit_counter = cu_reinit_counter + dt
      ENDIF


      IF (.TRUE.) THEN
         WRITE (6, '(A,F15.2)') ' currentTimeInSeconds :', currenttimeinseconds
         WRITE (6, '(A,F15.2)') ' nudgingStartTimeInSeconds :', nudgingstarttimeinseconds
         WRITE (6, '(A,F15.2)') ' nudgingEndTimeInSeconds :', nudgingendtimeinseconds
         WRITE (6, '(A,I15)') ' analysis_time :', analysis_time
         WRITE (6, '(A,F15.2)') ' itimestep*dt :', itimestep*dt
         WRITE (6, '(A,F15.2)') ' xtime :', xtime
         WRITE (6, '(A,F15.2)') ' nudgingRelaxation :', relaxation_time
         WRITE (6, '(A,F15.2)') ' dt :', dt
         WRITE (6, '(A,F15.2)') ' cu_reinit_counter :', cu_reinit_counter
      ENDIF

   END SUBROUTINE nudge_towards_analysis

   SUBROUTINE initialize_nudging_tendencies(rundgdten, rvndgdten, rwndgdten, rqvndgdten, rthndgdten, &
                                            rphndgdten, rmundgdten, &
                                            ids, ide, jds, jde, kds, kde, &
                                            ims, ime, jms, jme, kms, kme, &
                                            i_start, i_end, j_start, j_end, kts, kte, num_tiles)

      IMPLICIT NONE

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
      INTEGER, INTENT(IN) :: kts
      INTEGER, INTENT(IN) :: kte
      INTEGER, INTENT(IN) :: num_tiles

      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: i_start
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: i_end
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: j_start
      INTEGER, DIMENSION(num_tiles), INTENT(IN) :: j_end

      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: rundgdten
      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: rvndgdten
      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: rwndgdten
      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: rqvndgdten
      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: rthndgdten
      REAL, INTENT(INOUT), DIMENSION(ims:ime, kms:kme, jms:jme) :: rphndgdten
      REAL, INTENT(INOUT), DIMENSION(ims:ime, jms:jme) :: rmundgdten

      ! Local Variables
      INTEGER:: i
      INTEGER:: j
      INTEGER:: k
      INTEGER:: ij

!$OMP PARALLEL DO &
!$OMP PRIVATE ( ij,i,j,k )
      DO ij = 1, num_tiles

         DO j = j_start(ij), j_end(ij)
         DO k = kts, kte
         DO i = i_start(ij), i_end(ij) + 1 ! U staggered in X -> Need to add 1 to i_end

            rundgdten(i, k, j) = 0.

         ENDDO
         ENDDO
         ENDDO

         DO j = j_start(ij), j_end(ij) + 1 ! V staggered in Y -> Need to add 1 to j_end
         DO k = kts, kte
         DO i = i_start(ij), i_end(ij)

            rvndgdten(i, k, j) = 0.

         ENDDO
         ENDDO
         ENDDO

         DO j = j_start(ij), j_end(ij)
         DO k = kts, kte + 1 ! W staggered in Z -> Need to add 1 to kte
         DO i = i_start(ij), i_end(ij)

            rwndgdten(i, k, j) = 0.
            rphndgdten(i, k, j) = 0.

         ENDDO
         ENDDO
         ENDDO

         DO j = j_start(ij), j_end(ij)
         DO k = kts, kte
         DO i = i_start(ij), i_end(ij)

            rthndgdten(i, k, j) = 0.

            rqvndgdten(i, k, j) = 0.

         ENDDO
         ENDDO
         ENDDO

         DO j = j_start(ij), j_end(ij)
         DO i = i_start(ij), i_end(ij)

            rmundgdten(i, j) = 0.

         ENDDO
         ENDDO

      ENDDO
!$OMP END PARALLEL DO

   END SUBROUTINE initialize_nudging_tendencies
!-------------------------------------------------------------------

!-------------------------------------------------------------------
END MODULE module_analysis_nudging
