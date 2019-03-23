      subroutine dsLEinit
      implicit none
      include 'dssusy.h'
      include 'dsversion.h'  ! set DarkSUSY version
      include 'dssubversion.h'
      include 'dsdir.h'      ! set DarkSUSY root directory

      integer i

c...Startup

      write(*,*)
      write(*,*) '    DarkSUSY found at '//dsinstall
      write(*,*)
     &  '    *********************************************************'
      write(*,*) 
     &  '    *** Welcome to DarkSUSY version                       ***'
      write(*,*) '    *** ',dsversion,'***'
      write(*,*) '    *** ',dssubversion,'***'
      write(*,*) 
     &  '    *********************************************************'
      write(*,*)
      write(*,*) '    Initializing DarkSUSYLE capture routines...'      
      
      call dscapsetup

      write(*,*) '    done.'
      write(*,*)

      return
      end


