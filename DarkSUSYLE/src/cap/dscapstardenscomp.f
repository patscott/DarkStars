***********************************************************************
*** dscapstardenscomp calculates the number density of a given element 
*** as a function of radius in a given star
*** Input: r = radius (meters)
***        itype = element type (as defined in dscapsetup)           
*** Output: number density (number/cm^3)
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: June 26, 2007
*** Modified: 2007-07-05 Pat Scott (pat@physto.se); comments fixed 07-24
***           ;added espsilon 2008-02-19
***********************************************************************

      real*8 function dscapstardenscomp(r_in,itype)
      implicit none
      
      include 'dscapstar.h'
      real*8 epsilon
	  real*8 r,r_in,rpl
      real*8 starmfr_temp(meshpoints),starmfrd2_temp(meshpoints)
      integer itype,j
      !real*8 dsntsundenscomp  ! Sun-specific density tables    
      parameter(epsilon=1d-10)

      if (itype.lt.1.or.itype.gt.n_species) then
        write(*,*) 
     &    'WARNING in dscapstardenscomp: illegal element type: ',
     &    itype
        dscapstardenscomp=0.0d0
        return
      endif

      r = r_in

      if (r.gt.r_star) then
        if (r.gt.r_star*(1.d0+epsilon)) then
          dscapstardenscomp=0.0d0
          return
        else
	      r = r_star
        endif
      endif

      if (r.lt.0.0d0) then
        dscapstardenscomp=starmfr(itype,1)
        return
      endif


      if (interpmode) then

        do j=1,meshpoints
          starmfr_temp(j) = starmfr(itype,j)
          starmfrd2_temp(j) = starmfrd2(itype,j)
        enddo
      
        call dssplint(starr,starmfr_temp,starmfrd2_temp,meshpoints,
     &   r/r_star,dscapstardenscomp)
        
        if (dscapstardenscomp.lt.0.0d0) dscapstardenscomp = 0.0d0
        
      else

c...Interpolate in table

        call dshunt(starr,meshpoints,r/r_star,j)
        if (j.lt.meshpoints) goto 20

        dscapstardenscomp=starmfr(itype,meshpoints)
        return

 20     rpl=(r-starr(j)*r_star)/(starr(j+1)*r_star-starr(j)*r_star)

        dscapstardenscomp=starmfr(itype,j)*(1.0d0-rpl)+
     &   starmfr(itype,j+1)*rpl

        !dscapstardenscomp=dsntsundenscomp(r,itype) ! Sun-specific density tables

      endif
      
      return
      
      end
