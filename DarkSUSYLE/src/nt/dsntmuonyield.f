*****************************************************************************
*** function dsntmuonyield gives the total yield of muons above threshold for
*** a given neutralino mass or the differential muon yield for a given
*** energy and a given angle. put yieldk=3 for integrated yield above
*** given thresholds and put yieldk=103 for differntial yield.
*** the annihilation branching ratios and
*** higgs parameters are extracted from susy.h and by calling dsandwdcosnn
*** wh='su' corresponds to annihilation in the sun and wh='ea' corresponds
*** to annihilation in the earth. if istat=1 upon return,
*** some inaccesible parts the differential muon spectra has been wanted,
*** and the returned yield should then be treated as a lower bound.
*** if istat=2 energetically forbidden annihilation channels have been
*** wanted. if istat=3 both of these things has happened.
*** units: 1.0e-30 m**-2 (annihilation)**-1 for integrated yield.
***        1.0e-30 m**-2 gev**-1 (degree)^-1 (annihilation)**-1 for
***        differential yield.
*** author: joakim edsjo  edsjo@physto.se
*** date: 96-03-19
*** modified: 96-09-03 to include new index order
*** modified: 97-12-03 to include new muyield routines (v3.21)
*****************************************************************************

      real*8 function dsntmuonyield(emuth0,thmax,wh,yieldk,istat)
      implicit none
      include 'dssusy.h'
      include 'dsmucom.h'
      include 'dsidtag.h'
      include 'dsprep.h'

c------------------------ functions ------------------------------------

      real*8 dsmuyield

c------------------------ variables ------------------------------------

      real*8 mneu,emuth0,yield,thmax
      integer ch,istat,yieldk
      character*2 wh

c----------------------------------------------- set-up common variables
      mneu=mass(kn(1))

      if (.not.dsprepcalled) then
        write(*,*) 'error in dsntmuonyield: dsprep must be called',
     &    ' before any rate calculation'
        write(*,*) 'begins for every new model. program stopping.'
      endif

c...loop through different channels and calculate the yield above threshold
c...for each channel.
c      call wirate(6,6,1)

      yield=0.0d0
      mfistat=0
      do 100 ch=1,14
        if (abr(ch).gt.0.0d0) then
          yield=yield+abr(ch)*dsmuyield(mneu,emuth0,
     &      thmax,ch,wh,yieldk,istat)
          mfistat=or(mfistat,istat)
        endif
  100 continue
      dsntmuonyield=yield
      istat=mfistat

      end


















