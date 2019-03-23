      module ez_magnitude
      USE ez_magnitude_data
      IMPLICIT none
      private
      public :: Get_Magnitudes

      type :: fe_h_list ! sorted in decreasing order of fe_h
         real :: fe_h
         type (fe_h_list), pointer :: nxt
         real :: mags(n_mags)
      end type

      type :: lgt_list ! sorted in decreasing order of lgt
         real :: lgt
         type (lgt_list), pointer :: nxt
         type (fe_h_list), pointer :: zlist
      end type

      type :: lgg_list ! sorted in decreasing order of lgg
         real :: lgg
         type (lgg_list), pointer :: nxt
         type (lgt_list), pointer :: tlist
      end type
      
      type (lgg_list), pointer :: ghead
      

      contains


      subroutine Get_Magnitudes(log_Teff, log_L, mass, Fe_H_in, results)
         double precision, intent(IN) :: log_Teff ! log10 of surface temp
         double precision, intent(IN) :: log_L ! log10 of luminosity in solar units
         double precision, intent(IN) :: mass ! mass in solar units
         double precision, intent(IN) :: Fe_H_in ! [Fe/H]
         double precision, intent(OUT) :: results(n_mags)
         
         double precision, parameter :: Zsol = 0.02d0, mag_bol_sun = 4.746d0
         double precision :: mag_bol
         
         real :: lgg, fe_h, lgt, results1(n_mags), results2(n_mags), alfa, beta
         type (lgg_list), pointer :: glist, gnxt
         
         lgg = log10(mass) + 4.0D0*log_Teff - log_L - 10.6071D0
         lgt = log_Teff
         fe_h = Fe_H_in
         
         mag_bol = mag_bol_sun - 2.5d0 * log_L
         
         if (.not. associated(ghead)) call Read_Magnitudes_Data
         
         if (lgg >= ghead % lgg) then ! use the largest lgg
            call get_tlist_results(ghead % tlist, lgt, fe_h, results1)
            results = results1
         else
            
            glist => ghead
            do while (associated(glist % nxt))
               gnxt => glist % nxt
               if (lgg == gnxt % lgg) then ! use gnxt
                  call get_tlist_results(gnxt % tlist, lgt, fe_h, results1)
                  results = results1
                  exit
               end if
               if (lgg >= gnxt % lgg) then ! interpolate between glist and gnxt
                  call get_tlist_results(glist % tlist, lgt, fe_h, results1)
                  call get_tlist_results(gnxt % tlist, lgt, fe_h, results2)
                  alfa = (lgg - gnxt % lgg) / (glist % lgg - gnxt % lgg)
                  beta = 1 - alfa
                  results = alfa * results1 + beta * results2
                  exit
               end if
               glist => gnxt
            end do
      
            if (.not. (associated(glist % nxt))) then
               ! use the smallest lgg
               call get_tlist_results(glist % tlist, lgt, fe_h, results1)
               results = results1
            end if
         
         end if
         
         results(bol) = mag_bol
         
      end subroutine Get_Magnitudes
      
      
      subroutine get_tlist_results(tlist, lgt, fe_h, results)
         type (lgt_list), pointer :: tlist
         real, intent(IN) :: lgt, fe_h
         real, intent(OUT) :: results(n_mags)
         
         type (lgt_list), pointer :: tnxt
         real :: results1(n_mags), results2(n_mags), alfa, beta
         
         if (.not. associated(tlist)) stop 'bad tlist for get_tlist_results'
         
         if (lgt >= tlist % lgt) then ! use the largest lgt
            call get_zlist_results(tlist % zlist, fe_h, results1)
            results = results1
            return
         end if

         do while (associated(tlist % nxt))
            tnxt => tlist % nxt
            if (lgt == tnxt % lgt) then ! use tnxt
               call get_zlist_results(tnxt % zlist, fe_h, results1)
               results = results1
               return
            end if
            if (lgt >= tnxt % lgt) then ! interpolate between tlist and tnxt
               call get_zlist_results(tlist % zlist, fe_h, results1)
               call get_zlist_results(tnxt % zlist, fe_h, results2)
               alfa = (lgt - tnxt % lgt) / (tlist % lgt - tnxt % lgt)
               beta = 1 - alfa
               results = alfa * results1 + beta * results2
               return
            end if
            tlist => tnxt
         end do
      
         ! use the smallest lgt
         call get_zlist_results(tlist % zlist, fe_h, results1)
         results = results1
      
      end subroutine get_tlist_results
      
      
      subroutine get_zlist_results(zlist, fe_h, results)
         type (fe_h_list), pointer :: zlist
         real, intent(IN) :: fe_h
         real, intent(OUT) :: results(n_mags)
         
         type (fe_h_list), pointer :: znxt
         real :: alfa, beta
         
         if (.not. associated(zlist)) stop 'bad zlist for get_zlist_results'
         
         if (fe_h >= zlist % fe_h) then ! use the largest fe_h
            results = zlist % mags
            return
         end if

         do while (associated(zlist % nxt))
            znxt => zlist % nxt
            if (fe_h == znxt % fe_h) then ! use znxt
               results = znxt % mags
               return
            end if
            if (fe_h >= znxt % fe_h) then ! interpolate between zlist and znxt
               alfa = (fe_h - znxt % fe_h) / (zlist % fe_h - znxt % fe_h)
               beta = 1 - alfa
               results = alfa * zlist % mags + beta * znxt % mags
               return
            end if
            zlist => znxt
         end do
      
         ! use the smallest fe_h
         results = zlist % mags
      
      end subroutine get_zlist_results
      
      
      subroutine Read_Magnitudes_Data
         use ez_data
         ! read file and build lists
         integer :: ios, istat, cnt
         character (len=256) :: fname
         real :: lgt, lgg, fe, mags(n_mags), fe_h
         type (lgg_list), pointer :: glist
         type (lgt_list), pointer :: tlist
         type (fe_h_list), pointer :: zlist
         integer :: num_entries, num_made

         IF (LEN(TRIM(EZ_DATA_DIR)) .GT. 0) THEN
            fname = TRIM(EZ_DATA_DIR) // '/../lcb98cor.dat'
         ELSE
            fname = 'lcb98cor.dat'
         END IF
         OPEN(UNIT=IO_UBV, FILE=TRIM(fname), ACTION='READ', STATUS='OLD', IOSTAT=ios)
         if (ios /= 0) stop 'failed to open lcb98cor.dat'
         
         num_entries = 0
         cnt = 0
         do while (.true.)
         
            read(IO_UBV,fmt=*,iostat=ios) lgt, lgg, fe_h, mags
            if (ios /= 0) then
               close(IO_UBV)
               exit
            end if
            cnt = cnt + 1
            lgt = log10(lgt)           
            
            call get_glist(ghead, lgg, glist)
            call get_tlist(glist % tlist, lgt, tlist)
            call get_zlist(tlist % zlist, fe_h, zlist, num_entries)
            
            zlist % mags = mags
            
         end do
         
         num_made = 0
         
         glist => ghead
         lgg = 1d30
         do while (associated(glist))
            if (glist % lgg >= lgg) stop 'bad glist order'
            lgg = glist % lgg
            tlist => glist % tlist
            lgt = 1d30
            do while (associated(tlist))
               if (tlist % lgt >= lgt) stop 'bad tlist order'
               lgt = tlist % lgt
               zlist => tlist % zlist
               fe_h = 1d30
               do while (associated(zlist))
                  if (zlist % fe_h >= fe_h) stop 'bad zlist order'
                  fe_h = zlist % fe_h
                  num_made = num_made + 1
                  zlist => zlist % nxt
               end do
               tlist => tlist % nxt
            end do
            glist => glist % nxt
         end do
         
         if (num_made /= num_entries) stop 'error in Read_Magnitudes_Data'
         
      end subroutine Read_Magnitudes_Data


      subroutine get_glist(head, lgg, glist)
         type (lgg_list), pointer :: head
         real, intent(IN) :: lgg
         type (lgg_list), pointer :: glist
         
         type (lgg_list), pointer :: g1, g2
         
         if (.not. associated(head)) then ! first time
            call alloc_glist
            head => glist
            return
         end if
         
         if (head % lgg == lgg) then ! matches head of list
            glist => head
            return
         end if
         
         if (head % lgg < lgg) then ! becomes new head of list
            call alloc_glist
            glist % nxt => head
            head => glist
            return
         end if
         
         ! check list
         g1 => head
         do while (associated(g1 % nxt))
            g2 => g1 % nxt
            if (g2 % lgg == lgg) then
               glist => g2; return
            end if
            if (g2 % lgg < lgg) then ! insert new one before g2 
               call alloc_glist
               glist % nxt => g2
               g1 % nxt => glist
               return
            end if
            g1 => g2
         end do
         ! add to end of list after g1
         call alloc_glist
         g1 % nxt => glist
         
         contains
         
         subroutine alloc_glist
            integer :: istat
            allocate(glist,stat=istat)
            if (istat /= 0) stop 'allocate failed in get_glist'
            nullify(glist % tlist)
            nullify(glist % nxt)
            glist % lgg = lgg       
         end subroutine alloc_glist
         
      end subroutine get_glist
                     

      subroutine get_tlist(head, lgt, tlist)
         type (lgt_list), pointer :: head
         real, intent(IN) :: lgt
         type (lgt_list), pointer :: tlist
         
         type (lgt_list), pointer :: t1, t2
         
         if (.not. associated(head)) then ! first time
            call alloc_tlist
            head => tlist
            return
         end if
         
         if (head % lgt == lgt) then ! matches head of list
            tlist => head
            return
         end if
         
         if (head % lgt < lgt) then ! becomes new head of list
            call alloc_tlist
            tlist % nxt => head
            head => tlist
            return
         end if
         
         ! check list
         t1 => head
         do while (associated(t1 % nxt))
            t2 => t1 % nxt
            if (t2 % lgt == lgt) then
               tlist => t2; return
            end if
            if (t2 % lgt < lgt) then ! insert new one before t2 
               call alloc_tlist
               tlist % nxt => t2
               t1 % nxt => tlist
               return
            end if
            t1 => t2
         end do
         ! add to end of list after t1
         call alloc_tlist
         t1 % nxt => tlist
         
         contains
         
         subroutine alloc_tlist
            integer :: istat
            allocate(tlist,stat=istat)
            if (istat /= 0) stop 'allocate failed in get_tlist'
            nullify(tlist % zlist)
            nullify(tlist % nxt)
            tlist % lgt = lgt       
         end subroutine alloc_tlist
         
      end subroutine get_tlist
                     

      subroutine get_zlist(head, fe_h, zlist, num_created)
         type (fe_h_list), pointer :: head
         real, intent(IN) :: fe_h
         type (fe_h_list), pointer :: zlist
         integer, intent(OUT) :: num_created
         
         type (fe_h_list), pointer :: z1, z2
         
         if (.not. associated(head)) then ! first time
            call alloc_zlist
            head => zlist
            return
         end if
         
         if (head % fe_h == fe_h) then ! matches head of list
            zlist => head
            return
         end if
         
         if (head % fe_h < fe_h) then ! becomes new head of list
            call alloc_zlist
            zlist % nxt => head
            head => zlist
            return
         end if
         
         ! check list
         z1 => head
         do while (associated(z1 % nxt))
            z2 => z1 % nxt
            if (z2 % fe_h == fe_h) then
               zlist => z2; return
            end if
            if (z2 % fe_h < fe_h) then ! insert new one before z2 
               call alloc_zlist
               zlist % nxt => z2
               z1 % nxt => zlist
               return
            end if
            z1 => z2
         end do
         ! add to end of list after z1
         call alloc_zlist
         z1 % nxt => zlist
         
         contains
         
         subroutine alloc_zlist
            integer :: istat
            allocate(zlist,stat=istat)
            if (istat /= 0) stop 'allocate failed in get_zlist'
            nullify(zlist % nxt)
            zlist % fe_h = fe_h     
            num_created = num_created + 1 
         end subroutine alloc_zlist
         
      end subroutine get_zlist
                     

      end module ez_magnitude

