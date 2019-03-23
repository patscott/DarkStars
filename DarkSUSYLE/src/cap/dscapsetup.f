      subroutine dscapsetup

***********************************************************************
*** Setup common blocks needed for general capture rate calculations
*** Author: Joakim Edsjo
*** Date: 2007-06-26
*** Modified: 2007-07-11 Pat Scott (pat@physto.se)
***   Accurate atomic masses and nuclear mass formula from AME2003
***   (Audi et al 2003, NucPhysA 729:129 & 729:337)
***********************************************************************

      implicit none
      include 'dscapstar.h'
      include 'dsparam.h'

      real*8 binding, Z
      integer i
      
      logical called
      data called/.false./
      save called

      binding(Z) = 1.44381d-8*Z**2.39d0+1.55468d-15*Z**5.35d0

      if (called) return

      staraa(1)        = 1.0d0            ! H1; mass number A
      starza(1)        = 1.0d0            ! proton number Z
      starma_atomic(1) = 1.00782503207d0  ! atomic mass [a.m.u.]
      
      staraa(2)        = 4.0d0            ! He4
      starza(2)        = 2.0d0   
      starma_atomic(2) = 4.00260325415d0        
      
      staraa(3)        = 12.0d0           ! C12
      starza(3)        = 6.0d0  
      starma_atomic(3) = 12.0d0       
      
      staraa(4)        = 14.0d0           ! N14
      starza(4)        = 7.0d0  
      starma_atomic(4) = 14.00307400478d0       
      
      staraa(5)        = 16.0d0           ! O16
      starza(5)        = 8.0d0  
      starma_atomic(5) = 15.99491461956d0       
      
      staraa(6)        = 20.0d0           ! Ne20
      starza(6)        = 10.0d0  
      starma_atomic(6) = 19.99244017542d0       
      
      staraa(7)        = 24.0d0           ! Mg24
      starza(7)        = 12.0d0  
      starma_atomic(7) = 23.985041700d0       
      
      staraa(8)        = 23.0d0           ! Na23
      starza(8)        = 11.0d0  
      starma_atomic(8) = 22.98976928087d0       
      
      staraa(9)        = 27.0d0           ! Al27
      starza(9)        = 13.0d0  
      starma_atomic(9) = 26.98153863d0       
      
      staraa(10)       = 28.0d0           ! Si28
      starza(10)       = 14.0d0  
      starma_atomic(10)= 27.97692653246d0       
      
      staraa(11)       = 32.0d0           ! S32
      starza(11)       = 16.0d0  
      starma_atomic(11)= 31.97207100d0       
      
      staraa(12)       = 40.0d0           ! Ar40
      starza(12)       = 18.0d0  
      starma_atomic(12)= 39.96238312251d0       
      
      staraa(13)       = 40.0d0           ! Ca40
      starza(13)       = 20.0d0  
      starma_atomic(13)= 39.96259098d0       
      
      staraa(14)       = 56.0d0           ! Fe56
      starza(14)       = 26.0d0  
      starma_atomic(14)= 55.9349375d0       
      
      staraa(15)       = 58.0d0           ! Ni58
      starza(15)       = 28.0d0  
      starma_atomic(15)= 57.9353429d0    
      
      staraa(16)       = 60.0d0           ! Ni60
      starza(16)       = 28.0d0  
      starma_atomic(16)= 59.9307864d0
      
      
      if (n_species.gt.16) then

        staraa(17)       = 3.0d0            ! He3
        starza(17)       = 2.0d0  
        starma_atomic(17)= 3.01602931914d0       
        
        staraa(18)       = 13.0d0           ! C13
        starza(18)       = 6.0d0  
        starma_atomic(18)= 13.00335483778d0       
        
        staraa(19)       = 18.0d0           ! O18
        starza(19)       = 8.0d0  
        starma_atomic(19)= 17.9991610d0       
        
        staraa(20)       = 208.0d0           ! Pb208
        starza(20)       = 82.0d0  
        starma_atomic(20)= 207.9766521d0
        
        staraa(21)       = 207.0d0           ! Pb207
        starza(21)       = 82.0d0  
        starma_atomic(21)= 206.9758969d0       
        
        staraa(22)       = 206.0d0           ! Pb206
        starza(22)       = 82.0d0  
        starma_atomic(22)= 205.9744653d0       
        
      endif

      do i=1,n_species
        starma_atomic(i) = starma_atomic(i)*amu2GeV ! convert atomic masses to GeV
        starma(i) = starma_atomic(i) -starza(i)*m_e +binding(starza(i)) ! nuclear masses [GeV]
      end do

      called=.true.
      
      return

      end      

