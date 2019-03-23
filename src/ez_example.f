      MODULE ez_example
      USE ez_do_one
      IMPLICIT NONE
      
      ! This does standard evolution of a single star and produces the standard log files as output.
      ! Edit Do_Example below to specify the mass and metallicity.
      ! The output log files will be written to the current (run) directory (see 'demos/DATA_README' for info on the logs).
      ! Edit Initial_Params if you need to change the default control settings for the run.
      
      CONTAINS
      
      SUBROUTINE Do_Example
      DOUBLE PRECISION :: mass
      metals_dir = '../metals/z02'  ! metallicity -- use z0001, z0003, z001, z004, z01, z02, or z03
      mass = 1.9D0  ! initial mass -- anything from 0.3 to 100
      CALL Init_Do_One_Utils(mass)
      CALL Get_Star_Info(mass, Dummy_Build_Filename, Initial_Params, Dummy_Before_Evolve, Do_One_Check_Model)
      END SUBROUTINE Do_Example
      
      SUBROUTINE Initial_Params ! see star_controls.f for details of what can be set in this routine
      wind_Eta = 1D0  ! set Reimers' wind eta
      END SUBROUTINE Initial_Params
      
      END MODULE ez_example
