      MODULE ez_driver ! This is the procedural interface to the stellar evolution code.
      IMPLICIT NONE
      
      CONTAINS
      
      SUBROUTINE EZ_Start ( metals_dir )
      ! The metalicity is determined by metals_dir, the path to one of the "metals" subdirectories.
      ! You need to call this routine at the start before doing anything else.
      USE star_constants, ONLY: strlen
      USE ez_setup, ONLY: SETUP
      USE ez_data, ONLY: EZ_DATA_DIR
      CHARACTER (LEN=strlen) :: metals_dir
      EZ_DATA_DIR = trim(metals_dir)
      CALL SETUP
      END SUBROUTINE EZ_Start
      

      LOGICAL FUNCTION EZ_ZAMS ( initial_Mass, Initialize_Params )
      ! EZ_ZAMS initializes a ZAMS model of the requested mass.  (Remember to call EZ_Start first!)
      ! Then it sets up default values in star_data and calls Initialize_Params to let it modify those defaults
      ! before finishing the initialization of the model.
      USE ez_setup, ONLY: Read_ZAMS
      USE ez_data, ONLY: init_M, init_Z
      USE ez_shell_data, ONLY: CZS
      DOUBLE PRECISION, INTENT(IN) :: initial_Mass
      INTERFACE
         SUBROUTINE Initialize_Params
         END SUBROUTINE Initialize_Params
      END INTERFACE
      EZ_ZAMS = Read_ZAMS ( initial_Mass, Initialize_Params )
      init_M = initial_Mass
      init_Z = CZS
      END FUNCTION EZ_ZAMS
      

      SUBROUTINE EZ_Evolve ( Check_Model )
      ! After you call EZ_ZAMS to set up a model, the next thing is to call EZ_Evolve to evolve it.
      ! It calls your Check_Model routine at the end of each timestep to inspect the star and decide what to do.
      ! Anything in star_data is fair game for Check_Model to consult to make its decisions.
      ! EZ_ZAMS keeps going until either Check_Model says to quit, or it is forced to quit by something such as
      ! an attempt to reduce the timestep to below the dynamical timescale.
      USE ez_cycle, ONLY: Evolve
      INTERFACE
         INTEGER FUNCTION Check_Model() ! return BACKUP, TERMINATE, or KEEP_GOING
         END FUNCTION Check_Model
      END INTERFACE
      INTEGER :: JO
      CALL Evolve ( JO, Check_Model )
      END SUBROUTINE EZ_Evolve
      

      LOGICAL FUNCTION EZ_Save ( filename )
      ! Saves the system state to the specified file for future restoring.  Saves both the star and the parameters.
      ! You can do this at any stage of evolution, either from your top-level exec or from your Check_Model routine.
      ! The save files are binary and are NOT meant to be used for anything but save and restore.
      ! They are valid only for a specific version of the system and are not intended for long term storage of data.
      ! That said, they are very useful in searching to meet some goal since you can save an evolved base state
      ! and then try different parameter values by restoring the saved state instead of recomputing it. 
      USE ez_data, ONLY: SysOut
      USE star_constants, ONLY: strlen
      CHARACTER (LEN=strlen), INTENT(IN) :: filename
      INTEGER :: IOS
      CALL SysOut ( filename, IOS )
      EZ_Save = IOS .EQ. 0
      END FUNCTION EZ_Save
      

      LOGICAL FUNCTION EZ_Restore ( filename )
      ! Restores the system state from the specified file so you can go back to try again with different settings.
      ! NOTE 1: You must call EZ_Start at startup even if you are going to call EZ_Restore next.
      ! NOTE 2: If you restore a file saved by a different version of the system, the code will get very confused.
      USE ez_data, ONLY: SysIn
      USE star_constants, ONLY: strlen
      CHARACTER (LEN=strlen), INTENT(IN) :: filename
      INTEGER :: IOS
      CALL SysIn ( filename, IOS )
      EZ_Restore = IOS .EQ. 0
      END FUNCTION EZ_Restore
      

      LOGICAL FUNCTION EZ_Trade_X_for_Y ( target_Y, Check_Model )
      ! EZ_Trade_X_for_Y trades hydrogen for helium to reach the requested helium mass fraction Y.
      ! It only does small steps so that the model will have a chance to adjust gradually to the change.
      ! Your Check_Model can follow the changes as they happen.
      ! This routine should typically be called after EZ_ZAMS and before EZ_Evolve.
      USE ez_setup, ONLY: Trade_X_for_Y
      DOUBLE PRECISION, INTENT(IN) :: target_Y
      INTERFACE
         INTEGER FUNCTION Check_Model() ! return BACKUP, TERMINATE, or KEEP_GOING
         END FUNCTION Check_Model
      END INTERFACE
      EZ_Trade_X_for_Y = Trade_X_for_Y ( target_Y, Check_Model )
      END FUNCTION EZ_Trade_X_for_Y
      

      SUBROUTINE EZ_Extras
      ! EZ_Extras fills in the data in the star_extras module.
      ! The info in star_extras is more expensive to calculate than that in star_data, 
      ! so it is only created on demand.
      USE ez_report, ONLY: Get_Extra_Model_Data
      CALL Get_Extra_Model_Data
      END SUBROUTINE EZ_Extras
      

      SUBROUTINE EZ_Shell_Extras ( shell )
      ! Get the array SX information for the specified shell only.
      ! NB: gets the info that depends only on the one shell, but not that which depends on other shells.
      ! for example, this gets the local burn rates at the shell for the various nuclear reactions, 
      !    but it doesn not get the integrated power generated interior to the shell.
      ! If you need the integrated values, you'll have to use EZ_Extras and calculate everything.
      USE ez_report, ONLY: Get_Extra_Shell_Data
      INTEGER, INTENT(IN) :: shell
      CALL Get_Extra_Shell_Data ( shell )
      END SUBROUTINE EZ_Shell_Extras      
      

      LOGICAL FUNCTION EZ_Post_He_Flash ( total_Mass, core_Mass, age_Offset, Check_Model )
      ! This function is intended for post-Helium flash evolution.
      ! EZ_Post_He_Flash initializes a model of the requested total mass and helium core mass.
      ! The new model will also have the same ZAMS metallicity, but the detailed composition profile will of course be different.
      ! For example, the C-N-O equilibrium abundances will certainly be different in the new model after the flash settles.
      ! The age for the resulting model will be the current age plus the age offset argument.  (e.g,. I use 1 or 2% of star_Age for this)
      ! Your Check_Model can follow the changes as the system builds the new model for you.
      ! Here's how I use this routine.
      !    For masses in the helium flash range, my Check_Model notices when helium burning reaches break even with neutrino losses.
      !    Rather than letting the flash run its course, which will cause the system to terminate, I call this function to switch models.
      !    The necessary information for this call is easily available in star_data.  I use star_Mass, mass_He_Core, and star_Age.
      USE ez_flash, ONLY: Get_Post_He_Flash
      DOUBLE PRECISION, INTENT(IN) :: total_Mass, core_Mass, age_Offset
      INTERFACE
         INTEGER FUNCTION Check_Model() ! return BACKUP, TERMINATE, or KEEP_GOING
         END FUNCTION Check_Model
      END INTERFACE
      EZ_Post_He_Flash = Get_Post_He_Flash ( total_Mass, core_Mass, age_Offset, Check_Model )
      END FUNCTION EZ_Post_He_Flash
      

      SUBROUTINE EZ_Reduce_DTY ( next_DTY )
      USE ez_cycle_data, ONLY: REQUESTED_DTY
      DOUBLE PRECISION, INTENT(IN) :: next_DTY ! next time step (in years)
      ! This will apply to the next timestep only.  Then system automatically goes back to normal.
      ! If the requested timestep is larger than what the system would have done, it is ignored.
      ! In other words, this con only reduce the time step, not increase it.
      ! The system also limits how small you can go to twice the current dynamical timescale.
      ! Why have such a control?  Sometimes you need a small time step so can hit a a desired mass, for example.
      ! By controlling both the mass loss rate AND the time step, you can make sure you hit the target.
      REQUESTED_DTY = next_DTY
      END SUBROUTINE EZ_Reduce_DTY
      

      SUBROUTINE Null_Initial_Params
      END SUBROUTINE Null_Initial_Params
      

      END MODULE ez_driver
