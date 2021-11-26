      implicit none
c--
c-- nmaxpass is the maximum nr of passes that can be stored over the to be analysed period
c-- nmaxpointing is the maximum nr of timesteps that can be stored per single pass
c-- nsatmax is the maximum nr of satellites that can be analysed in one single run
c
      integer*4 nmaxpass, nmaxpointing, nsatmax, nmaxcross
      parameter (nmaxpass=1000, nmaxpointing=40000, nsatmax=50, nmaxcross=50)
c
      real*8    sat_az_pass(nmaxpointing,nmaxpass)    , sat_el_pass(nmaxpointing,nmaxpass)    , time_pass(nmaxpointing,nmaxpass)
      real*8    sat_az_track_la(nmaxpointing,nmaxpass), sat_el_track_la(nmaxpointing,nmaxpass), time_track_la(nmaxpointing,nmaxpass)
      real*8    sat_pass_start(nmaxpass), sat_pass_end(nmaxpass), sat_el_max_pass(nmaxpass)
      integer*4 npointspass(nmaxpass), n_track(nmaxpass), n_track_la(nmaxpass), i_sat_prio(nsatmax), i_start(nsatmax), i_end(nsatmax)
      integer*4 nrpasses, key(nmaxpass), ok(nmaxpass), sat_pass_prio(nmaxpass)
      integer*1 ovlflag(nmaxpointing), sat_co_pass(nmaxpointing,nmaxpass), sat_co_track_la(nmaxpointing,nmaxpass)
      character*20 c_sat_pass(nmaxpass), c_sat_use(nsatmax), c_sat_type(nsatmax)
      logical   swns(nmaxpass)
c-- Diagnostic variables - can be deleted if plots no longer needed
      integer*4 i1, i2, i3, i4, i5, i6
      real*8    angle1, angle2, azpark, elpark, aztst, eltst
c
      real*8    doy_tle, time, ro(3), vo(3), long, lat, theta0g, time_start, sat_az, sat_el, rsc, sat_el_max, time_max, el_track_min
      real*8    torad, todeg, track_angle, x1, y1, z1, x2, y2, z2, time_save, secsubsteps, time_offset, rsc1, rsc2, rsc3, rsc4, elpriothresh
      real*8    azrate, elrate, azoff, eloff, azrange(2), elrange(2), sunaz, sunel, sunaz1, sunel1, sunaz2, sunel2, anglesun, sunrate
      real*8    suntim1, suntim2, sunrepoint, sunahead, elpower, azpower, delta_t_la, azhor, elhor, azpoint, elpoint, latsav
      real*8    sat_co_az, sat_co_el, rsc5, rsc6
c-- Raster variables
      real*8    grid_az, grid_el, grid_az_range, grid_el_range, grid_az_step, grid_el_step, grid_delta_t, az_cmd, el_cmd, t_move, az1, el1, azavg, elavg
      integer*4 ncrossaz, ncrossel
      real*4    azresp(nmaxcross), elresp(nmaxcross), azx(nmaxcross), elx(nmaxcross), sumy, sumxy, p1, p2, p3
      logical   swgrid, swcross
c--
      integer*4 dattim(8), luntle, iyr_tle, doy, i, j, k, l, ihr, imn, i_year, i_month, i_day, i_pass, i_track, lunsat, nsatuse, ios, isat
      integer*4 i_arg1, i_arg, ihr1, ihr2, ihr3, ihr4, imn1, imn2, imn3, imn4, ihr_ahead, lunhor, ih1, ih2, ih3, ih4, i_mask1, i_mask2
      integer*4 ihr5, imn5, ihr6, imn6, i_co1, i_co2
      integer*2 npos
      logical foundtle, swdebug, swlive, swallowflip, swnoovl, swboth, swoffset, swazrange, swelrange, swtiltcor, swsun, swcalib, swazelalt, swhorcmd
      logical swnoload, swfound, swpoint, swdoublec, swjumpin, swhmask, check_horizon, swcovis
      character*1    csq, cdq
      character*3    c_month(12), cloc
      character*5    c_zone
      character*8    c_date, c_datvis
      character*10   c_time, c_horcmd, c_datehor1, c_datehor2, c_time1, c_time2, c_time3, c_time4
      character*17   c_datehormatch
      character*20   c_sat
      character*50   c_file, c_command, cfmt, clong, clat
      character*69   tlestring
      character*250  argstring, tledir, command, chortxt
      character*2000 horizons_command
c
      data c_month/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/
c
c-- Run the program as follows:
c-- ./tletrack.exe \"AQUA\" resource beam=1.8 etc - TBD
c
c--   First argument satellite name as in the TLE between double quotes (in UNIX command line speak) or an indirect to a file using @
c--   Examples:
c--   ./tletrack.exe \"NOAA\ 19\" weather beam=1.8
c--   ./tletrack.exe @./examples/tlesats.txt beam=1.8
c--
c-- Where tlesats.txt contained:
c--
c-- "FENGYUN 3B"   weather   9
c-- "FENGYUN 3C"   weather   9
c-- "METOP-A"      weather   5
c-- "METOP-B"      weather   8
c-- "METOP-C"      weather   8
c-- "NOAA 15"      weather   4
c-- "NOAA 18"      weather   7
c-- "NOAA 19"      weather   7
c--
c-- with the first " left aligned in column 1 !
c--
c--   Second argument is the type  of TLE : resource or weather IF the first argument is not an indirect !
c--                   NOT required if first argument is reference to indirect list (using @)
c--
c--   Further arguments:  
c--
c--     beam=1.8  as beamwidth of dish to be used
c--     debug     show the commands that would be issued for reactive and pro-active (look ahead)
c--     diag      show the az el values at each time step for reactive and pro-active (look ahead) - L....o....n....g
c--     live      actually command the hardware
c--     flip      allow dish flip mode if crossing the N vector
c--     noovl     removes Az, El commands overlapping in time - hope this cures the MD-02 M2 pulse timeout problem .......
c--     both      do the calculations with and without look-ahead - only useful with debug or diag
c
c-- The TLE environment variable defines where to find the TLE files
c-- Test tiltcor
c      do i = 1, 36
c        eltst = 30.0
c        aztst = dble(i-1)* 10.0D0
c        call tilt_correct_cos(aztst, eltst)
c        write (*,*) aztst, eltst
c      enddo
c      read (*,*)
c-- Test tiltcor_rot - the complete version
c      aztst = 115.0D0
c      eltst = 30.0D0
c      call tilt_correct_rot(aztst, eltst)
c      write (*,*) aztst, eltst
c      read (*,*)      
c      aztst = 192.3D0
c      eltst = 29.8D0
c      call tilt_correct_rot(aztst, eltst)
c      write (*,*) aztst, eltst
c      read (*,*)      
c-- WSL clock drift test
c      call date_and_time(c_date, c_time, c_zone, dattim)
c      i = dattim(7)
c      do while (.true.)
c        call date_and_time(c_date, c_time, c_zone, dattim)
c        if (dattim(7).ne.i) then
c          write (*,*) (dattim(k), k=5,8)
c          i = dattim(7)
c        endif
c      enddo
c
      torad           = 2.0D0 * dasin(1.0D0) / 180.0D0
      todeg           = 180.0D0 / (2.0D0 * dasin(1.0D0))
      el_track_min    = 3.0D0
      track_angle     = 1.8D0
c-- Use 0.05 sec time steps (required for X-band)
      secsubsteps     = 20.0D0
      swdebug         = .false.
      swlive          = .false.
      swallowflip     = .false.
      swnoovl         = .false.
c-- Power setting for M1 (Az) and M2 (El) on SPID MD-0n controller.
      azpower         = 100.0D0 / 100.0D0
      elpower         =  70.0D0 / 100.0D0
      azrate          = 4.5D0 * azpower
      elrate          = 5.2D0 * elpower
c-- ELRate @ 70% max power is 3.6 deg / sec - a tiny bit lower then expected - probably not significant.
      swboth          = .false.
      swoffset        = .false.
      azoff           = 0.0D0
      eloff           = 0.0D0
      swgrid          = .false.
      swcross         = .false.
      ncrossaz        = 0
      ncrossel        = 0
      swazrange       = .false.
      swelrange       = .false.
      azrange(1)      =   0.0D0
      azrange(2)      = 360.0D0
      elrange(1)      =   0.0D0
      elrange(2)      =  90.0D0
      swtiltcor       = .false.
      azpark          = 209.0D0
      elpark          = 0.0D0
      swsun           = .false.
      swcalib         = .false.
      swazelalt       = .false.
      swhorcmd        = .false.
      swnoload        = .false.
      swpoint         = .false.
      swdoublec       = .false.
      swjumpin        = .false.
      swhmask         = .false.
      swcovis         = .false.
c-- If two satellites passes have a conflict and their priorities differ by 1 and their max elevations are above this limit : the max elevation wins
      elpriothresh    = 30.0D0
      delta_t_la      = 0.0D0
      ihr_ahead       = 24
c
      do i = 1, nmaxpass
        sat_pass_start(i) = 0.0D0
        sat_pass_end(i)   = 0.0D0
        ok(i)             = 1
        swns(i)           = .false.
      enddo
      call getarg(1, argstring)
      if (index(argstring(1:lnblnk(argstring)),'-help').ne.0) then
        call tletrack_help()
        stop '** Help end **'
      endif
      if (argstring(1:1).eq.'@') then
        if (argstring(2:4).eq.'sun') then
          swsun = .true.
        else if (argstring(2:8).eq.'horcmd=') then
          swhorcmd = .true.
          do i = 1, len(c_horcmd)
            c_horcmd(i:i) = ' '
          enddo
          c_horcmd   = argstring(9:lnblnk(argstring))
        else
          call get_lun(lunsat)
          open (unit=lunsat,file=argstring(2:lnblnk(argstring)),form='formatted')
          nsatuse = 0
          do while (nsatuse.lt.nsatmax)
            read (lunsat,'(a)',iostat=ios) argstring
            if (ios.ne.0) goto 2
            nsatuse = nsatuse + 1
            i1 = 1
            do while (argstring(i1:i1).ne.'"')
              i1 = i1 + 1
            enddo
            i2 = i1 + 1
            do while (argstring(i2:i2).ne.'"')
              i2 = i2 + 1
            enddo
            c_sat_use(nsatuse) = argstring(i1:i2)
            i1 = i2 + 1
            do while (argstring(i1:i1).eq.' ')
              i1 = i1 + 1
            enddo
            i2 = i1
            do while (argstring(i2:i2).ne.' ')
              i2 = i2 + 1
            enddo
            c_sat_type(nsatuse) = argstring(i1:i2-1)
            i1 = i2 + 1
            do while (argstring(i1:i1).eq.' ')
              i1 = i1 + 1
            enddo
            i2 = i1
            do while (argstring(i2:i2).ne.' ')
              i2 = i2 + 1
            enddo
            call string_to_i4(argstring(i1:i2-1)//' ', npos, i_sat_prio(nsatuse))
          enddo
 2        close (unit=lunsat)
          call free_lun(lunsat)
        endif
        i_arg1 = 2
      else
        nsatuse = 1
        c_sat_use(nsatuse)  = argstring(1:lnblnk(argstring))
        call getarg(2, argstring)
        c_sat_type(nsatuse) = argstring(1:lnblnk(argstring))
        i_sat_prio(nsatuse) = 9
        i_arg1 = 3
      endif
      do k = i_arg1, 30
        call getarg(k, argstring)
        if (index(argstring(1:lnblnk(argstring)),'debug').ne.0)  swdebug      = .true.
        if (index(argstring(1:lnblnk(argstring)),'beam=').ne.0) then
          call string_to_r8(argstring(index(argstring,'beam=')+5:lnblnk(argstring)), npos, track_angle)
          write (*,'(''Modified track angle                     : '',F10.3)') track_angle
        endif
        if (index(argstring(1:lnblnk(argstring)),'live').ne.0)   swlive       = .true.
        if (index(argstring(1:lnblnk(argstring)),'flip').ne.0)   swallowflip  = .true.
        if (index(argstring(1:lnblnk(argstring)),'noovl').ne.0)  swnoovl      = .true.
        if (index(argstring(1:lnblnk(argstring)),'both').ne.0)   swboth       = .true.
        if (index(argstring(1:lnblnk(argstring)),'azoff=').ne.0) then
          call string_to_r8(argstring(index(argstring,'azoff=')+6:lnblnk(argstring)), npos, azoff)
          write (*,'(''Azimuth Offset                           : '',F10.3)') azoff
          swoffset = .true.
        endif
        if (index(argstring(1:lnblnk(argstring)),'eloff=').ne.0) then
          call string_to_r8(argstring(index(argstring,'eloff=')+6:lnblnk(argstring)), npos, eloff)
          write (*,'(''Elevation Offset                         : '',F10.3)') eloff
          swoffset = .true.
        endif
        if (index(argstring(1:lnblnk(argstring)),'raster=').ne.0) then
          if (swcross) stop '** raster= and cross= cannot be combined **'
          i = index(argstring,'=')
          j = index(argstring,',')
          call string_to_r8(argstring(i+1:j-1)//' ',npos, grid_az)
          i = j
          if (index(argstring(j+1:),',').eq.0) stop '** Improper raster format **'
          j = index(argstring(j+1:),',') + j
          call string_to_r8(argstring(i+1:j-1)//' ',npos, grid_el)
          i = j
          if (index(argstring(j+1:),',').eq.0) stop '** Improper raster format **'
          j = index(argstring(j+1:),',') + j
          call string_to_r8(argstring(i+1:j-1)//' ',npos, grid_az_range)
          i = j
          if (index(argstring(j+1:),',').eq.0) stop '** Improper raster format **'
          j = index(argstring(j+1:),',') + j
          call string_to_r8(argstring(i+1:j-1)//' ',npos, grid_el_range)
          i = j
          if (index(argstring(j+1:),',').eq.0) stop '** Improper raster format **'
          j = index(argstring(j+1:),',') + j
          call string_to_r8(argstring(i+1:j-1)//' ',npos, grid_az_step)
          i = j
          if (index(argstring(j+1:),',').eq.0) stop '** Improper raster format **'
          j = index(argstring(j+1:),',') + j
          call string_to_r8(argstring(i+1:j-1)//' ',npos, grid_el_step)
          i = j
          call string_to_r8(argstring(i+1:lnblnk(argstring))//' ',npos, grid_delta_t)
          write (*,'(7F10.2)') grid_az, grid_az_range, grid_az_step, grid_el, grid_el_range, grid_el_step, grid_delta_t
          swgrid = .true.
        endif
        if (index(argstring(1:lnblnk(argstring)),'cross=').ne.0) then
          if (swgrid) stop '** raster= and cross= cannot be combined **'
          i = index(argstring,'=')
          j = index(argstring,',')
          call string_to_r8(argstring(i+1:j-1)//' ',npos, grid_az)
          i = j
          if (index(argstring(j+1:),',').eq.0) stop '** Improper raster format **'
          j = index(argstring(j+1:),',') + j
          call string_to_r8(argstring(i+1:j-1)//' ',npos, grid_el)
          i = j
          if (index(argstring(j+1:),',').eq.0) stop '** Improper raster format **'
          j = index(argstring(j+1:),',') + j
          call string_to_r8(argstring(i+1:j-1)//' ',npos, grid_az_range)
          i = j
          if (index(argstring(j+1:),',').eq.0) stop '** Improper raster format **'
          j = index(argstring(j+1:),',') + j
          call string_to_r8(argstring(i+1:j-1)//' ',npos, grid_el_range)
          i = j
          if (index(argstring(j+1:),',').eq.0) stop '** Improper raster format **'
          j = index(argstring(j+1:),',') + j
          call string_to_r8(argstring(i+1:j-1)//' ',npos, grid_az_step)
          i = j
          call string_to_r8(argstring(i+1:lnblnk(argstring))//' ',npos, grid_el_step)
          write (*,'(6F10.2)') grid_az, grid_az_range, grid_az_step, grid_el, grid_el_range, grid_el_step
          swcross = .true.
          if (int(grid_az_range/grid_az_step)+1.gt.nmaxcross) stop '** Requested nr of Az cross steps exceeds maximum **'
          if (int(grid_el_range/grid_el_step)+1.gt.nmaxcross) stop '** Requested nr of Az cross steps exceeds maximum **'
        endif
        if (index(argstring(1:lnblnk(argstring)),'azrange=').ne.0) then
          i = index(argstring,'=')
          j = index(argstring,'-')
          if (j.eq.0) then
            j = lnblnk(argstring)
            call string_to_r8(argstring(i+1:j)//' '  ,npos, azrange(1))
          else
            i = index(argstring,'=')
            call string_to_r8(argstring(i+1:j-1)//' ',npos, azrange(1))
            i = j
            j = lnblnk(argstring)
            call string_to_r8(argstring(i+1:j)//' '  ,npos, azrange(2))
          endif
          write (*,'(''Lower Azimuth limit                      : '',F10.3)') azrange(1)
          write (*,'(''Upper Azimuth limit                      : '',F10.3)') azrange(2)
          swazrange = .true.
          call SPID_MD02_AzLim(azrange)
        endif
        if (index(argstring(1:lnblnk(argstring)),'elrange=').ne.0) then
          i = index(argstring,'=')
          j = index(argstring,'-')
          if (j.eq.0) then
            j = lnblnk(argstring)
            call string_to_r8(argstring(i+1:j)//' '  ,npos, elrange(1))
          else
            i = index(argstring,'=')
            call string_to_r8(argstring(i+1:j-1)//' ',npos, elrange(1))
            i = j
            j = lnblnk(argstring)
            call string_to_r8(argstring(i+1:j)//' '  ,npos, elrange(2))
          endif
          write (*,'(''Lower Elevation limit                    : '',F10.3)') elrange(1)
          write (*,'(''Upper Elevation limit                    : '',F10.3)') elrange(2)
          swelrange = .true.
          call SPID_MD02_ElLim(elrange)
        endif
        if (index(argstring(1:lnblnk(argstring)),'tiltcor').ne.0) then   
          swtiltcor = .true.
          call SPID_MD02_tiltcor()
        endif
        if (index(argstring(1:lnblnk(argstring)),'azpark=').ne.0) then
          call string_to_r8(argstring(index(argstring,'azpark=')+7:lnblnk(argstring)), npos, azpark)
          write (*,'(''Modified Azimuth Park Position           : '',F10.3)') azpark
        endif
        if (index(argstring(1:lnblnk(argstring)),'elpark=').ne.0) then
          call string_to_r8(argstring(index(argstring,'elpark=')+7:lnblnk(argstring)), npos, elpark)
          write (*,'(''Modified Elevation Park Position         : '',F10.3)') elpark
        endif
        if (index(argstring(1:lnblnk(argstring)),'calib').ne.0)  swcalib = .true.
        if (index(argstring(1:lnblnk(argstring)),'deltat=').ne.0) then
          call string_to_r8(argstring(index(argstring,'deltat=')+7:lnblnk(argstring)), npos, delta_t_la)
          write (*,'(''Delta Time for commands                  : '',F10.3)') delta_t_la
        endif
        if (index(argstring(1:lnblnk(argstring)),'azelalt').ne.0) swazelalt = .true.
        if (index(argstring(1:lnblnk(argstring)),'eltrackmin=').ne.0) then
          call string_to_r8(argstring(index(argstring,'eltrackmin=')+11:lnblnk(argstring)), npos, el_track_min)
          write (*,'(''Max Elevation of track has to be .gt.    : '',F10.3)') el_track_min
        endif
        if (index(argstring(1:lnblnk(argstring)),'predict=').ne.0) then
          call string_to_i4(argstring(index(argstring,'predict=')+8:lnblnk(argstring))//' ', npos, ihr_ahead)
          write (*,'(''Modified predict window (hr)             : '',I10)') ihr_ahead
        endif
        if (index(argstring(1:lnblnk(argstring)),'noload').ne.0) swnoload = .true.
        if (index(argstring(1:lnblnk(argstring)),'point=').ne.0) then
          i = index(argstring,'=')
          j = index(argstring,',')
          if (j.eq.0) then
            stop '** Format error in point= command, example: point=192.3,29.7 **'
          else
            swpoint = .true.
            i = index(argstring,'=')
            call string_to_r8(argstring(i+1:j-1)//' ',npos, azpoint)
            i = j
            j = lnblnk(argstring)
            call string_to_r8(argstring(i+1:j)//' '  ,npos, elpoint)
          endif
        endif
        if (index(argstring(1:lnblnk(argstring)),'doublec').ne.0) swdoublec = .true.
        if (index(argstring(1:lnblnk(argstring)),'jumpin').ne.0) swjumpin = .true.
        if (index(argstring(1:lnblnk(argstring)),'hmask').ne.0) swhmask = .true.
        if (index(argstring(1:lnblnk(argstring)),'covis=').ne.0) then
          swcovis = .true.
          cloc = argstring(7:9)
          if (cloc.ne.'SVB'.and.cloc.ne.'KIR') stop '** Unknown co-visibility location **'
        endif
      enddo
c-- Just execute a simple point and stare and exit
      if (swpoint) then
        write (*,'(''Az pointing to                           : '',F10.3)') azpoint
        write (*,'(''El pointing to                           : '',F10.3)') elpoint
        if (swlive) call SPID_MD02(azpoint, elpoint, swdebug, .false.)
        if (swlive.and.swdoublec) then
          write (c_command,'(''sleep '',F8.1)') 1.0
          call system(c_command)
          call SPID_MD02(azpoint, elpoint, swdebug, .false.)
        endif
        stop '** Pointing completed **'
      endif
c-- Check if dish alive and aligned - using visual inspection from my study :-)
c-- Remove this as it is not really used
      if (swlive) then
c        call SPID_MD02(250.1D0, 0.0D0, swdebug, .false.)
c        write (*,'(''If dish moved and aligned press enter :'')')
c        read (*,*)
        call SPID_MD02(azpark + azoff, elpark + eloff, swdebug, .false.)
c-- Wait for dish to return to the park position
        t_move = abs(250.1D0 - azpark)/azrate
        write (c_command,'(''sleep '',F8.1)') t_move
        call system(c_command)
        if (swcalib) call calibrate(azpark, elpark, azoff, eloff)
      endif
c-- Retrieve JPL Horizons Az/El if requested
      if (swhorcmd) then
        call date_and_time(c_date, c_time, c_zone, dattim)
        call datetodoy(dattim(1), dattim(2), dattim(3), doy)
        call doytodate(doy + 2, dattim(1), c_datvis)
        c_datehor1 = c_date(1:4)//'-'//c_date(5:6)//'-'//c_date(7:8)
        c_datehor2 = c_datvis(1:4)//'-'//c_datvis(5:6)//'-'//c_datvis(7:8) 
        csq = char(ichar("'"))
        cdq = char(ichar('"'))
        if (.not.swnoload) then
          horizons_command = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='//csq//c_horcmd(1:lnblnk(c_horcmd))//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'CENTER='//csq//'coord'//csq//'&'
          call getenv ('LONG',clong)
          call getenv ('LAT',clat)
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'SITE_COORD='//csq//clong(1:lnblnk(clong))//','//clat(1:lnblnk(clat))//',0.00'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'MAKE_EPHEM='//csq//'YES'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'TABLE_TYPE='//csq//'OBS'//csq//'&'
c-- Remember, half the step size is used some 76 lines down for an offset calc !!!
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'STEP_SIZE='//csq//'1%20m'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'CAL_FORMAT='//csq//'CAL'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'TIME_DIGITS='//csq//'MINUTES'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'ANG_FORMAT='//csq//'HMS'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'OUT_UNITS='//csq//'KM-S'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'RANGE_UNITS='//csq//'KM'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'APPARENT='//csq//'AIRLESS'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'SOLAR_ELONG='//csq//'0,180'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'SUPPRESS_RANGE_RATE='//csq//'NO'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'SKIP_DAYLT='//csq//'NO'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'REF_SYSTEM='//csq//'J2000'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'CSV_FORMAT='//csq//'YES'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'OBJ_DATA='//csq//'YES'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'QUANTITIES='//csq//'2,3,4,20'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'START_TIME='//csq//c_datehor1//'%2000:00'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'STOP_TIME='//csq//c_datehor2//'%2023:59'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'REF_PLANE='//csq//'ECLIPTIC'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'COORD_TYPE='//csq//'GEODETIC'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'AIRMASS='//csq//'38.0'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'EXTRA_PREC='//csq//'NO'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'ELM_LABELS='//csq//'YES'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'TP_TYPE='//csq//'ABSOLUTE'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'R_T_S_ONLY='//csq//'NO'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'CA_TABLE_TYPE='//csq//'STANDARD'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'TCA3SG_LIMIT='//csq//'14400'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'CALIM_SB='//csq//'0.1'//csq//'&'
          write (horizons_command(lnblnk(horizons_command)+1:),'(a)') 'CALIM_PL='//csq//'.1,.1,.1,.1,1.0,1.0,1.0,1.0,.1,.003'//csq
          write (*,'(a)') '........ Querying JPL Horizons - Pls wait'
          call system('curl -s -k '//cdq//horizons_command(1:lnblnk(horizons_command))//cdq//' -o tempeph.txt')
          write (*,'(a)') 'Finished Querying JPL Horizons'
        endif
        call get_lun(lunhor)
        open (unit=lunhor, file='tempeph.txt',form='formatted',access='sequential')
        read (lunhor,'(a)') chortxt
        read (lunhor,'(a)') chortxt
        write (*,'(a)') chortxt(1:lnblnk(chortxt))
        do while (index(chortxt,'$$SOE').eq.0)
          read (lunhor,'(a)') chortxt
        enddo
c-- Read the actual data lines and crudely :-( isolate Az/El
        call date_and_time(c_date, c_time, c_zone, dattim)
        write (c_datehormatch,'(a4,''-'',a3,''-'',a2,1x,I2.2,'':'',I2.2)') c_date(1:4), c_month(dattim(2)), c_date(7:8), dattim(5) - (dattim(4)/60), dattim(6)
        read (lunhor,'(a)') chortxt
        swfound = .false.
        do while (index(chortxt,'$$EOE').eq.0.and.(.not.swfound))
          ih1 = index(chortxt,',')
          ih2 = ih1 + 1
          do i = 1,6
            k = index(chortxt(ih2:),',')
            ih2 = ih2 + k
          enddo
          ih3 = ih2
          k = index(chortxt(ih3:),',')
          ih3 = ih3 + k
          ih4 = ih3
          k = index(chortxt(ih4:),',')
          ih4 = ih4 + k
          if (index(chortxt(1:ih1-1),c_datehormatch).ne.0) then
            swfound = .true.
            call string_to_r8(chortxt(ih2:ih3-2), npos, azhor)
            call string_to_r8(chortxt(ih3:ih4-2), npos, elhor)
            if (azrange(1).le.azhor.and.azhor.le.azrange(2).and.elrange(1).le.elhor.and.elhor.le.elrange(2)) then
              if (elhor.ge.el_track_min) write (*,'('' Time (UT), Az, El                            : '',a,2F10.3)') chortxt(1:ih1-1), azhor, elhor
              if (swlive) then
                call SPID_MD02(azhor + azoff, elhor + eloff, swdebug, .false.)
              endif
            endif
          endif
          read (lunhor,'(a)') chortxt
        enddo
c-- First entry found. Now loop over entries
        do while (index(chortxt,'$$EOE').eq.0)
 11       call date_and_time(c_date, c_time, c_zone, dattim)
c-- Modify the system time half a STEP_SIZE into the future to execute the commands 30 seconds early - achieves look ahead - fails accross a leap second :-))
c-- A fool proof method would be to convert time to Julian date add/subtract the offset and convert back again. This would be robust against leap seconds
          call offset_time(dattim, 30)
          write (c_datehormatch,'(a4,''-'',a3,''-'',a2,1x,I2.2,'':'',I2.2)') c_date(1:4), c_month(dattim(2)), c_date(7:8), dattim(5) - (dattim(4)/60), dattim(6)
          ih1 = index(chortxt,',')
          ih2 = ih1 + 1
          do i = 1,6
            k = index(chortxt(ih2:),',')
            ih2 = ih2 + k
          enddo
          ih3 = ih2
          k = index(chortxt(ih3:),',')
          ih3 = ih3 + k
          ih4 = ih3
          k = index(chortxt(ih4:),',')
          ih4 = ih4 + k
          if (index(chortxt(1:ih1-1),c_datehormatch).eq.0) then
            call system('sleep 1.0')
            goto 11
          endif
c-- Now decode Az/El and check constraints and if OK launch command
          call string_to_r8(chortxt(ih2:ih3-2), npos, azhor)
          call string_to_r8(chortxt(ih3:ih4-2), npos, elhor)
          if (azrange(1).le.azhor.and.azhor.le.azrange(2).and.elrange(1).le.elhor.and.elhor.le.elrange(2)) then
            if (elhor.ge.el_track_min) write (*,'('' Time (UT), Az, El                            : '',a,2F10.3)') chortxt(1:ih1-1), azhor, elhor
            if (swlive) then
              call SPID_MD02(azhor + azoff, elhor + eloff, swdebug, .false.)
              if (swdoublec) then
                write (c_command,'(''sleep '',F8.1)') 1.0
                call system(c_command)
                call SPID_MD02(azhor + azoff, elhor + eloff, swdebug, .false.)
              endif
            endif
          endif
          read (lunhor,'(a)') chortxt
        enddo
        close (unit=lunhor)
        call free_lun(lunhor)
        read (*,*)
      endif
c-- Track the sun if requested
      if (swsun) then
        call sunazel_init()
        call sunazel(dattim, sunaz1, sunel1, 0.0D0)
        suntim1 = float(dattim(5)) * 3600.0 + float(dattim(6)) * 60.0 + float(dattim(7)) + (float(dattim(8))/1000.0)
        call system('sleep 20.0')
        call sunazel(dattim, sunaz2, sunel2, 0.0D0)
        suntim2 = float(dattim(5)) * 3600.0 + float(dattim(6)) * 60.0 + float(dattim(7)) + (float(dattim(8))/1000.0)
        call sunangle(sunaz1, sunel1, sunaz2, sunel2, anglesun)
        sunrate    = anglesun / (suntim2 - suntim1)
        sunrepoint = track_angle / sunrate
        sunahead   = sunrepoint / 2.0
        write (*,'(''Solar repointing every                   : '',F10.3,'' seconds'')') sunrepoint
        write (*,'(''Solar look ahead time                    : '',F10.3,'' seconds'')') sunahead
        do while (.true.)
          call sunazel(dattim, sunaz, sunel, sunahead)
          write (*,'(2I4,3F10.3)') dattim(5), dattim(6), float(dattim(7)) + (float(dattim(8))/1000.0), sunaz, sunel
          if (swoffset) then
            sunaz = sunaz + azoff
            sunel = sunel + eloff
          endif
          if (swlive) then
            call SPID_MD02(sunaz, sunel, swdebug, .false.)
          endif
c-- You could include the rate calibration here again and execute a sleep of sunrepoint - (suntim2-suntim1)
c          call sunazel(dattim, sunaz1, sunel1, 0.0D0)
c          suntim1 = float(dattim(5)) * 3600.0 + float(dattim(6)) * 60.0 + float(dattim(7)) + (float(dattim(8))/1000.0)
c          call system('sleep 20.0')
c          call sunazel(dattim, sunaz2, sunel2, 0.0D0)
c          suntim2 = float(dattim(5)) * 3600.0 + float(dattim(6)) * 60.0 + float(dattim(7)) + (float(dattim(8))/1000.0)
c          call sunangle(sunaz1, sunel1, sunaz2, sunel2, anglesun)
c          sunrate    = anglesun / (suntim2 - suntim1)
c          sunrepoint = track_angle / sunrate
c          sunahead   = sunrepoint / 2.0
c-- End of rate calibration
c          if (sunrepoint - (suntim2 - suntim1).gt.0.0D0) write (c_command,'(''sleep '',F8.1)') sunrepoint - (suntim2 - suntim1)
          write (c_command,'(''sleep '',F8.1)') sunrepoint
          call system(c_command)
        enddo
      endif
c-- Implement the grid option
      if (swgrid) then
        az1 = azpark
        el1 = elpark
        k   = nint((grid_az_range/2.0D0)/grid_az_step)
        l   = nint((grid_el_range/2.0D0)/grid_el_step)
        do i = -l, l
          el_cmd = grid_el + dble(i) * grid_el_step
          do j = -k, k
            az_cmd = grid_az + dble(j) * grid_az_step
c-- estimate a 0.3 sec commanding/response overhead for the SPID_MD02
            t_move = max(abs(az_cmd - az1)/azrate,abs(el_cmd - el1)/elrate) + grid_delta_t + 0.3
            write (c_command,'(''sleep '',F8.1)') t_move
            call system(c_command)
            az1 = az_cmd
            el1 = el_cmd
c-- Retrieve time of command
            call date_and_time(c_date, c_time, c_zone, dattim)
c
            write (*,'(''Az : '',F10.2,'' El : '',F10.2,3x,I4.4,1x,I2.2,1x,I2.2,3x,I2.2,1x,I2.2,1x,I2.2)') az_cmd, el_cmd, (dattim(k),k=1,3), (dattim(k),k=5,7)
            if (swlive) call SPID_MD02(az_cmd, el_cmd, swdebug, .false.)
          enddo
        enddo
c-- Send to rest position
        write (c_command,'(''sleep '',F8.1)') grid_delta_t + 1.0
        call system(c_command)
        call SPID_MD02(azpark, elpark, swdebug, .false.)
        if (swcalib) call calibrate(azpark, elpark, azoff, eloff)
        stop '** End Grid **'
      endif
c-- Implement the cross option
      if (swcross) then
        k   = nint((grid_az_range/2.0D0)/grid_az_step)
        l   = nint((grid_el_range/2.0D0)/grid_el_step)
        el_cmd = grid_el
        do i = -k, k
          ncrossaz = ncrossaz + 1
          az_cmd = grid_az + dble(i) * grid_az_step
          write (*,'(''Az : '',F10.2,'' El : '',F10.2)') az_cmd, el_cmd
          azx(ncrossaz) = sngl(az_cmd)
          if (swlive) call SPID_MD02(az_cmd, el_cmd, swdebug, .false.)
          write (*,'(''Enter response : '')')
          read  (*,*) azresp(ncrossaz)
        enddo
        az_cmd = grid_az
        do i = -l, l
          ncrossel = ncrossel + 1
          el_cmd   = grid_el + dble(i) * grid_el_step
          write (*,'(''Az : '',F10.2,'' El : '',F10.2)') az_cmd, el_cmd
          elx(ncrossel) = sngl(el_cmd)
          if (swlive) call SPID_MD02(az_cmd, el_cmd, swdebug, .false.)
          write (*,'(''Enter response : '')')
          read  (*,*) elresp(ncrossel)
        enddo
c-- Calculate and send to peak position
        write (c_command,'(''sleep '',F8.1)') 3.0
        call system(c_command)
c
        call estimate_peak(azx, azresp, ncrossaz, p1, p2, p3)
        az_cmd = p3
        call estimate_peak(elx, elresp, ncrossel, p1, p2, p3)
        el_cmd = p3
c
        write (*,'(''Azimuth   - Response'')')
        do i = 1, ncrossaz
          write (*,'(F10.2,3x,F10.2)') azx(i), azresp(i)
        enddo
        write (*,'(''Elevation - Response'')')
        do i = 1, ncrossel
          write (*,'(F10.2,3x,F10.2)') elx(i), elresp(i)
        enddo
        call SPID_MD02(az_cmd, el_cmd, swdebug, .false.)
        write (*,'(''Peak system response at : '')')
        write (*,'(''Az : '',F10.2,'' El : '',F10.2)') az_cmd, el_cmd
        stop '** End Cross **'
      endif
c
      write (*,'(''Satellites in use (name, type, priority) : '')')
      do i = 1,nsatuse
        write (*,'(a,a,i5)') c_sat_use(i), c_sat_type(i), i_sat_prio(i)
      enddo
c-- Loop over the sats requested
      call date_and_time(c_date, c_time, c_zone, dattim)
      write (c_file(1:4),'(I4.4)') dattim(1)
      write (c_file(5:6),'(I2.2)') dattim(2)
      write (c_file(7:8),'(I2.2)') dattim(3)
      call datetodoy(dattim(1), dattim(2), dattim(3), doy)
      time_start = dble(dattim(5) * 3600 + dattim(6) * 60 + dattim(7)) - dattim(4) * 60.0D0
      c_file(9:9) = '-'
c-- Start the loop over the satellites
      nrpasses = 0
      call get_lun(luntle)
      do isat = 1, nsatuse
        call datetodoy(dattim(1), dattim(2), dattim(3), doy)
        i_start(isat)         = nrpasses + 1
        c_sat = c_sat_use(isat)(1:lnblnk(c_sat_use(isat)))
        write (c_file(10:),'(a,a)') c_sat_type(isat)(1:lnblnk(c_sat_type(isat))), '.tle'
        call getenv('TLE',tledir)
        foundtle = .false.
        call exist(tledir(1:lnblnk(tledir))//c_file, foundtle)
        if (foundtle) then
          call system('cp '//tledir(1:lnblnk(tledir))//c_file//' ./weather.txt')
        else
          stop 'TLE file not found'
        endif
        call system('cat weather.txt | grep -A 2 -e '//c_sat(1:lnblnk(c_sat))//' - > hrpt.tmp')
        call system('tail -n +2 hrpt.tmp > hrpt.tle')
        call system('rm hrpt.tmp')
        open (unit=luntle,file='hrpt.tle',form='formatted')
        read (luntle,'(a)') tlestring
        rewind (luntle)
        read (tlestring(19:20),'(i2)') iyr_tle
        iyr_tle = iyr_tle + 2000
        read (tlestring(21:32),*) doy_tle
        close (unit=luntle)
        if (swdebug) write (*,'(''TLE Reference date : '',I5,2x,F8.2)') iyr_tle, doy_tle
        call sgp4_init(ro, vo, long, lat)
        time        = time_start
        time_offset = 0.0D0
        i_pass   = 0
        do i = 1, ihr_ahead * 3600 * nint(secsubsteps)
          time = time + (1.0D0 / secsubsteps)
          if (time.gt.86400.0D0) then
            doy = doy + 1
            time        = time        - 86400.0D0
            time_offset = time_offset + 86400.0D0
          endif
          call run_TLE(time, doy, ro, vo, long, lat, theta0g, doy_tle, iyr_tle, dattim(1))
          call azel(long, lat, ro, sat_az, sat_el, swazelalt, 'LOC')
          if (sat_el.gt.0.0D0) then
            if (i_pass.eq.0) then
              i_pass                   = 1
              nrpasses                 = nrpasses + 1
              npointspass(nrpasses)    = 0
              sat_pass_start(nrpasses) = time  + time_offset + dattim(4) * 60.0D0
              c_sat_pass(nrpasses)     = c_sat_use(isat)
              sat_pass_prio(nrpasses)  = i_sat_prio(isat) 
            endif
            if (npointspass(nrpasses).eq.nmaxpointing) then
              if (i_pass.eq.1) then
                if (sat_pass_end(nrpasses).eq.0.0D0) sat_pass_end(nrpasses) = time  + time_offset + dattim(4) * 60.0D0
              endif
              i_pass = 0
              goto 6
            endif
            npointspass(nrpasses) = npointspass(nrpasses) + 1
c- Crude determination whether a pass is S => N or N => S
            if (npointspass(nrpasses).eq.1) latsav = lat
            if (npointspass(nrpasses).eq.nint(secsubsteps).and.lat.lt.latsav) swns(nrpasses) = .true.
c
            call doytodate(doy, dattim(1), c_datvis)
            ihr = int(time/3600.0D0)
            imn = (int(time - dble(ihr)*3600.0D0)/60.0D0)
            rsc = time - dble(ihr)*3600.0D0 - imn * 60.0D0
            sat_az_pass(npointspass(nrpasses), nrpasses) = sat_az
            sat_el_pass(npointspass(nrpasses), nrpasses) = sat_el
            time_pass  (npointspass(nrpasses), nrpasses) = time + dattim(4) * 60.0D0
            if (swdebug.and.(time-dble(int(time))).lt.0.001D0) write (*,'(2F10.3,1x,a8,1x,2I5,F10.3, I5, I6, 2F10.3)') sat_az, sat_el, c_datvis, ihr, imn, rsc, doy, npointspass(nrpasses), long, lat
            if (swcovis) then
              call azel(long, lat, ro, sat_co_az, sat_co_el, swazelalt, cloc)
              if (sat_co_el.gt.3.0D0) sat_co_pass(npointspass(nrpasses), nrpasses) = 1
              if (sat_co_el.le.3.0D0) sat_co_pass(npointspass(nrpasses), nrpasses) = 0
            endif
          else
            if (i_pass.eq.1) then
              if (sat_pass_end(nrpasses).eq.0.0D0) sat_pass_end(nrpasses) = time  + time_offset + dattim(4) * 60.0D0
c-- Try to eliminate the 0 second passes that creeped in after the change to the Az/El calculation. (Do not understand where they come from)
              if (nint(dble(npointspass(nrpasses))/secsubsteps).eq.0) nrpasses = nrpasses - 1
            endif
            i_pass = 0
          endif
 6        continue
        enddo
        i_end(isat) = nrpasses
        call sgp4_end()
        write (*,*)
        write (*,'(''Pass # Secs   Max El   Time of max El        SAT ID       #LA Commands'')')
        do i = i_start(isat), i_end(isat)
          n_track(i)    = 0
          n_track_la(i) = 0
          sat_el_max    = 0.0
          do j = 1, npointspass(i)
            if (sat_el_pass(j,i).gt.sat_el_max) then
              sat_el_max = sat_el_pass(j,i)
              time_max   = time_pass(j,i)
              sat_el_max_pass(i) = sat_el_max
            endif
          enddo
c-- Flag the passes to be ignored for planning/conflict resolution below
          if (sat_el_max.lt.el_track_min) ok(i) = 0
c
          ihr = int(time_max/3600.0D0)
          imn = (int(time_max - dble(ihr)*3600.0D0)/60.0D0)
          rsc = time_max - dble(ihr)*3600.0D0 - imn * 60.0D0
c-- Pro-active tracking - look ahead (arrays with _la at end of name)
          i_track = 0
          j       = 0
          do while (j.lt.npointspass(i))
            j = j + 1
            if (sat_el_pass(j,i).gt.el_track_min.and.i_track.eq.0) then
              i_track                          = 1
              n_track_la(i)                    = 1
              sat_az_track_la(n_track_la(i),i) = sat_az_pass(j,i)
              sat_el_track_la(n_track_la(i),i) = sat_el_pass(j,i)
              time_track_la(n_track_la(i),i)   = time_pass(j,i)
              if (swcovis) then
                sat_co_track_la(n_track_la(i),i) = sat_co_pass(j,i)
              endif
            endif
            if (sat_el_pass(j,i).lt.el_track_min.and.i_track.eq.1) then
              i_track                          = 0
              n_track_la(i)                    = n_track_la(i) + 1
              sat_az_track_la(n_track_la(i),i) = sat_az_pass(j,i)
              sat_el_track_la(n_track_la(i),i) = sat_el_pass(j,i)
              time_track_la(n_track_la(i),i)   = time_pass(j,i)
              if (swcovis) then
                sat_co_track_la(n_track_la(i),i) = sat_co_pass(j,i)
              endif
            endif
            if (i_track.eq.1) then
              x1 = cos(torad*sat_az_track_la(n_track_la(i),i)) * cos(torad*sat_el_track_la(n_track_la(i),i))
              y1 = sin(torad*sat_az_track_la(n_track_la(i),i)) * cos(torad*sat_el_track_la(n_track_la(i),i))
              z1 =                                               sin(torad*sat_el_track_la(n_track_la(i),i))
              x2 = cos(torad*sat_az_pass(j,i)) * cos(torad*sat_el_pass(j,i))
              y2 = sin(torad*sat_az_pass(j,i)) * cos(torad*sat_el_pass(j,i))
              z2 =                               sin(torad*sat_el_pass(j,i))
c-- Find the half angle point
              if (todeg * acos(min(max((x1*x2 + y1*y2 + z1*z2),-1.0D0),1.0D0)).ge.track_angle/2.0D0) then
                time_save = time_pass(j,i)
                do while (j.lt.npointspass(i))
                  j = j + 1
                  x1 = cos(torad*sat_az_track_la(n_track_la(i),i)) * cos(torad*sat_el_track_la(n_track_la(i),i))
                  y1 = sin(torad*sat_az_track_la(n_track_la(i),i)) * cos(torad*sat_el_track_la(n_track_la(i),i))
                  z1 =                                               sin(torad*sat_el_track_la(n_track_la(i),i))
                  x2 = cos(torad*sat_az_pass(j,i)) * cos(torad*sat_el_pass(j,i))
                  y2 = sin(torad*sat_az_pass(j,i)) * cos(torad*sat_el_pass(j,i))
                  z2 =                               sin(torad*sat_el_pass(j,i))
c-- Find the full angle point and store with the time of the half-angle point - so point earlier to a position of the satellite in the future
                  if (todeg * acos(min(max((x1*x2 + y1*y2 + z1*z2),-1.0D0),1.0D0)).ge.track_angle) then
                    n_track_la(i)                    = n_track_la(i) + 1
                    sat_az_track_la(n_track_la(i),i) = sat_az_pass(j,i)
                    sat_el_track_la(n_track_la(i),i) = sat_el_pass(j,i)
                    time_track_la(n_track_la(i),i)   = time_save
                    if (swcovis) then
                      sat_co_track_la(n_track_la(i),i) = sat_co_pass(j,i)
                    endif
                    goto 1
                  endif
                enddo
 1              continue
              endif
            endif
          enddo
          if (swdebug.and.n_track_la(i).gt.0) write (*,'(''Pro-active tracking - look ahead'')')
          if (swdebug.and.n_track_la(i).gt.0) write (*,'(''Pointing #  Az          El     Time'')')
          do j = 1, n_track_la(i)
            ihr = int(time_track_la(j,i)/3600.0D0)
            imn = (int(time_track_la(j,i) - dble(ihr)*3600.0D0)/60.0D0)
            rsc = time_track_la(j,i) - dble(ihr)*3600.0D0 - imn * 60.0D0
            if (swdebug) write (*,'(I5,3x,2F10.3,2I5,F10.3)') j, sat_az_track_la(j,i), sat_el_track_la(j,i), ihr, imn, rsc
          enddo
          ihr = int(time_max/3600.0D0)
          imn = (int(time_max - dble(ihr)*3600.0D0)/60.0D0)
          rsc = time_max - dble(ihr)*3600.0D0 - imn * 60.0D0
c-- Print and also try to eliminate the 0 second passes that creeped in after the change to the Az/El calculation. Do not understand where they come from)
          if (sat_el_max.gt.el_track_min) then
            k = lnblnk(c_sat_use(isat))
            if (swdebug) then
              write (cfmt,'(a,I2.2,a)') '(2i5,F10.3,2I5,F10.3,5x,A,',14-k,'x,I5,3F10.2)'
              write (*,cfmt) i, nint(dble(npointspass(i))/secsubsteps), sat_el_max, ihr, imn, rsc, c_sat_use(isat)(2:lnblnk(c_sat_use(isat))-1), n_track_la(i), sat_pass_start(i), sat_pass_end(i), sat_pass_end(i) - sat_pass_start(i)
            else
              write (cfmt,'(a,I2.2,a)') '(2i5,F10.3,2I5,F10.3,5x,A,',14-k,'x,I5)'
              write (*,cfmt) i, nint(dble(npointspass(i))/secsubsteps), sat_el_max, ihr, imn, rsc, c_sat_use(isat)(2:lnblnk(c_sat_use(isat))-1), n_track_la(i)
            endif
          endif          
        enddo
      enddo
      call free_lun(luntle)
c-- Sort on the start times of the passes and determine if there is overlap between passes
c-- If so apply priorities to eliminate conflicts and re-print a list of passes. 
      do i = 1, nrpasses
        key(i) = i
      enddo
      call qsort_r8(sat_pass_start, key, nrpasses)
      i = 0
      do while (i.lt.nrpasses)
        i = i + 1
        if (ok(key(i)).eq.0) goto 4
        j = i
        do while (j.lt.nrpasses)
          j = j + 1
          if (ok(key(i)).eq.0) goto 4
          if (ok(key(j)).eq.0) goto 5
          if (sat_pass_start(j).gt.sat_pass_end(key(i))) goto 4
c-- Identify conflicts - remember time_pass_start is time sorted !
          if (sat_pass_start(i).le.sat_pass_start(j).and.sat_pass_start(j).le.sat_pass_end(key(i))) then
c
c-- Conflict resolution logic:
c-- Conflict free                                                                                      prio flag 1
c-- If equal priority ==> highest elevation wins                                                       prio flag 2
c-- If priorities differ by 1 point, but elevation of both > 30 degrees ==> highest elevation wins     prio flag 3
c-- If not the previous bullet  ==> the highest priority wins                                          prio flag 4
c
            if (swdebug) write (*,'(''Conflict between passes : '',6I5)')  i, j, key(i), key(j), sat_pass_prio(key(i)), sat_pass_prio(key(j))
            if (sat_pass_prio(key(i)).eq.sat_pass_prio(key(j))) then
              if (sat_el_max_pass(key(i)).gt.sat_el_max_pass(key(j))) then
                ok(key(j)) = 0
                ok(key(i)) = 2
              else
                ok(key(i)) = 0
                ok(key(j)) = 2
              endif
            else
              if (abs(sat_pass_prio(key(i))-sat_pass_prio(key(j))).eq.1.and.sat_el_max_pass(key(i)).gt.elpriothresh.and.sat_el_max_pass(key(j)).gt.elpriothresh) then
                if (sat_el_max_pass(key(i)).gt.sat_el_max_pass(key(j))) then
                  ok(key(j)) = 0
                  ok(key(i)) = 3
                else
                  ok(key(i)) = 0
                  ok(key(j)) = 3
                endif
              else
                if (sat_pass_prio(key(i)).gt.sat_pass_prio(key(j))) then
                  ok(key(j)) = 0
                  ok(key(i)) = 4
                else
                  ok(key(i)) = 0
                  ok(key(j)) = 4
                endif
              endif
            endif
          endif
 5        continue
        enddo
 4      continue
      enddo
      write (*,*)
      do i = 1, nrpasses
        if (ok(key(i)).ge.1) then
          ihr1 = int(mod(sat_pass_start(i),86400.0D0) / 3600.0D0)
          imn1 =    (mod(sat_pass_start(i),86400.0D0) - dble(ihr1)*3600.0D0)/60.0D0
          rsc1 =     mod(sat_pass_start(i),86400.0D0) - dble(ihr1)*3600.0D0 - imn1 * 60.0D0
          ihr2 = int(mod(sat_pass_end(key(i)),86400.0D0) / 3600.0D0)
          imn2 =    (mod(sat_pass_end(key(i)),86400.0D0) - dble(ihr2)*3600.0D0)/60.0D0
          rsc2 =     mod(sat_pass_end(key(i)),86400.0D0) - dble(ihr2)*3600.0D0 - imn2 * 60.0D0
          if (swhmask) then
            i_mask1 = -1
            i_mask2 = -1
            do j = 1, n_track_la(key(i))
              if (check_horizon(sat_az_track_la(j, key(i)), sat_el_track_la(j, key(i))).and.i_mask1.lt.0) i_mask1 = j
            enddo
            do j = n_track_la(key(i)), 1, -1
              if (check_horizon(sat_az_track_la(j, key(i)), sat_el_track_la(j, key(i))).and.i_mask2.lt.0) i_mask2 = j
            enddo
c-- If co-visibility requested, find the visibility window
            if (swcovis) then
              i_co1 = -1
              i_co2 = -1
              do j = 1, n_track_la(key(i))
                if (sat_co_track_la(j, key(i)).eq.1.and.i_co1.lt.0) i_co1 = j
              enddo
              do j = n_track_la(key(i)), 1, -1
                if (sat_co_track_la(j, key(i)).eq.1.and.i_co2.lt.0) i_co2 = j
              enddo
              if (i_co2.le.i_mask1.or.i_co1.gt.i_mask2) then
                i_co1 = -1
                i_co2 = -1
              else
                i_co1 = max(i_mask1, i_co1)
                i_co2 = min(i_mask2, i_co2)
              endif
            endif
c--
            if (i_mask1.eq.-1.or.i_mask2.eq.-1) then
              write (*,'(I4,F10.3,2(2I5,F10.3,F10.3,2x),3x,a,'' No El commands above Horizon mask - skipping pass'')') ok(key(i)), sat_el_max_pass(key(i)), ihr1, imn1, rsc1, sat_az_track_la(1,key(i)), ihr2, imn2, rsc2, sat_az_track_la(n_track_la(key(i)),key(i)), c_sat_pass(key(i))
            else
              ihr3 = int(mod(time_track_la(max(i_mask1-1,1)                 , key(i)),86400.0D0) / 3600.0D0)
              imn3 =    (mod(time_track_la(max(i_mask1-1,1)                 , key(i)),86400.0D0) - dble(ihr3)*3600.0D0)/60.0D0
              rsc3 =     mod(time_track_la(max(i_mask1-1,1)                 , key(i)),86400.0D0) - dble(ihr3)*3600.0D0 - imn3 * 60.0D0
              ihr4 = int(mod(time_track_la(min(i_mask2+1,n_track_la(key(i))), key(i)),86400.0D0) / 3600.0D0)
              imn4 =    (mod(time_track_la(min(i_mask2+1,n_track_la(key(i))), key(i)),86400.0D0) - dble(ihr4)*3600.0D0)/60.0D0
              rsc4 =     mod(time_track_la(min(i_mask2+1,n_track_la(key(i))), key(i)),86400.0D0) - dble(ihr4)*3600.0D0 - imn4 * 60.0D0
              write (c_time1,'(I2.2,'':'',I2.2,'':'',F4.1)') ihr3, imn3, rsc3
              write (c_time2,'(I2.2,'':'',I2.2,'':'',F4.1)') ihr4, imn4, rsc4
              if (swcovis.and.(.not.(i_co1.eq.-1.or.i_co2.eq.-1))) then
                ihr5 = int(mod(time_track_la(max(i_co1-1,1)                 , key(i)),86400.0D0) / 3600.0D0)
                imn5 =    (mod(time_track_la(max(i_co1-1,1)                 , key(i)),86400.0D0) - dble(ihr5)*3600.0D0)/60.0D0
                rsc5 =     mod(time_track_la(max(i_co1-1,1)                 , key(i)),86400.0D0) - dble(ihr5)*3600.0D0 - imn5 * 60.0D0
                ihr6 = int(mod(time_track_la(min(i_co2+1,n_track_la(key(i))), key(i)),86400.0D0) / 3600.0D0)
                imn6 =    (mod(time_track_la(min(i_co2+1,n_track_la(key(i))), key(i)),86400.0D0) - dble(ihr6)*3600.0D0)/60.0D0
                rsc6 =     mod(time_track_la(min(i_co2+1,n_track_la(key(i))), key(i)),86400.0D0) - dble(ihr6)*3600.0D0 - imn6 * 60.0D0
                write (c_time3,'(I2.2,'':'',I2.2,'':'',F4.1)') ihr5, imn5, rsc5
                write (c_time4,'(I2.2,'':'',I2.2,'':'',F4.1)') ihr6, imn6, rsc6
                write (*,'(I4,F10.3,2(2I5,F10.3,F10.3,2x),3x,a,'' Visible from : '',a,'' To : '',a,'' Co-Visibility from : '',a,'' To : '',a)') 
     *            ok(key(i)), sat_el_max_pass(key(i)), ihr1, imn1, rsc1, sat_az_track_la(1,key(i)), ihr2, imn2, rsc2, sat_az_track_la(n_track_la(key(i)),key(i)), c_sat_pass(key(i)), c_time1, c_time2, c_time3, c_time4
              else
                write (*,'(I4,F10.3,2(2I5,F10.3,F10.3,2x),3x,a,'' Visible from : '',a,'' To : '',a)') 
     *            ok(key(i)), sat_el_max_pass(key(i)), ihr1, imn1, rsc1, sat_az_track_la(1,key(i)), ihr2, imn2, rsc2, sat_az_track_la(n_track_la(key(i)),key(i)), c_sat_pass(key(i)), c_time1, c_time2
              endif
              call azel_plot(sat_az_track_la(1,key(i)), sat_el_track_la(1,key(i)), i_mask1, i_mask2, n_track_la(key(i)), ihr3, imn3, rsc3, ihr4, imn4, rsc4, swhmask, swlive, c_sat_pass(key(i)), sat_pass_start(i), swns(key(i)), sat_co_track_la(1,key(i)), swcovis, cloc)
            endif
          else
            if (swcovis) then
              i_co1 = -1
              i_co2 = -1
              do j = 1, n_track_la(key(i))
                if (sat_co_track_la(j, key(i)).eq.1.and.i_co1.lt.0) i_co1 = j
              enddo
              do j = n_track_la(key(i)), 1, -1
                if (sat_co_track_la(j, key(i)).eq.1.and.i_co2.lt.0) i_co2 = j
              enddo
            endif
            if (swcovis.and.(.not.(i_co1.eq.-1.or.i_co2.eq.-1))) then
              ihr5 = int(mod(time_track_la(max(i_co1-1,1)                 , key(i)),86400.0D0) / 3600.0D0)
              imn5 =    (mod(time_track_la(max(i_co1-1,1)                 , key(i)),86400.0D0) - dble(ihr5)*3600.0D0)/60.0D0
              rsc5 =     mod(time_track_la(max(i_co1-1,1)                 , key(i)),86400.0D0) - dble(ihr5)*3600.0D0 - imn5 * 60.0D0
              ihr6 = int(mod(time_track_la(min(i_co2+1,n_track_la(key(i))), key(i)),86400.0D0) / 3600.0D0)
              imn6 =    (mod(time_track_la(min(i_co2+1,n_track_la(key(i))), key(i)),86400.0D0) - dble(ihr6)*3600.0D0)/60.0D0
              rsc6 =     mod(time_track_la(min(i_co2+1,n_track_la(key(i))), key(i)),86400.0D0) - dble(ihr6)*3600.0D0 - imn6 * 60.0D0
              write (c_time3,'(I2.2,'':'',I2.2,'':'',F4.1)') ihr5, imn5, rsc5
              write (c_time4,'(I2.2,'':'',I2.2,'':'',F4.1)') ihr6, imn6, rsc6
              write (*,'(I4,F10.3,2(2I5,F10.3,F10.3,2x),3x,a,'' Co-Visibility from : '',a,'' To : '',a)') ok(key(i)), sat_el_max_pass(key(i)), ihr1, imn1, rsc1, sat_az_track_la(1,key(i)), ihr2, imn2, rsc2, sat_az_track_la(n_track_la(key(i)),key(i)), c_sat_pass(key(i)), c_time3, c_time4
            else
              write (*,'(I4,F10.3,2(2I5,F10.3,F10.3,2x),3x,a)') ok(key(i)), sat_el_max_pass(key(i)), ihr1, imn1, rsc1, sat_az_track_la(1,key(i)), ihr2, imn2, rsc2, sat_az_track_la(n_track_la(key(i)),key(i)), c_sat_pass(key(i))
            endif
c-- Try and plot the azel profiles
            call azel_plot(sat_az_track_la(1,key(i)), sat_el_track_la(1,key(i)), 1, n_track_la(key(i)), n_track_la(key(i)), ihr1, imn1, rsc1, ihr2, imn2, rsc2, swhmask, swlive, c_sat_pass(key(i)), sat_pass_start(i), swns(key(i)), sat_co_track_la(1,key(i)), swcovis, cloc)
          endif
        endif
c-- Try and plot the azel profiles
      enddo
c-- The tracks have been sorted and conflicts settled - now implement the commanding as an option ?
      write (*,*)
      do i = 1, nrpasses
        if (ok(key(i)).ge.1) then
          if (swlive) then
            if (swnoovl) then
              do j = 1, nmaxpointing
                ovlflag(j) = 0
              enddo
            endif
            ihr1 = int(mod(sat_pass_start(i),86400.0D0) / 3600.0D0)
            imn1 =    (mod(sat_pass_start(i),86400.0D0) - dble(ihr1)*3600.0D0)/60.0D0
            rsc1 =     mod(sat_pass_start(i),86400.0D0) - dble(ihr1)*3600.0D0 - imn1 * 60.0D0
            ihr2 = int(mod(sat_pass_end(key(i)),86400.0D0) / 3600.0D0)
            imn2 =    (mod(sat_pass_end(key(i)),86400.0D0) - dble(ihr2)*3600.0D0)/60.0D0
            rsc2 =     mod(sat_pass_end(key(i)),86400.0D0) - dble(ihr2)*3600.0D0 - imn2 * 60.0D0
            write (*,'(I4,F10.2,2(2I5,F10.3,2x),3x,A)') ok(key(i)), sat_el_max_pass(key(i)), ihr1, imn1, rsc1, ihr2, imn2, rsc2, c_sat_pass(key(i))
            call azel_command(sat_pass_start(i), sat_pass_end(key(i)), time_track_la(1,key(i)), sat_az_track_la(1,key(i)), sat_el_track_la(1,key(i)), n_track_la(key(i)), swdebug, swallowflip, swnoovl, ovlflag, azrate, elrate, azoff, eloff, swoffset, azpark, elpark, delta_t_la, swjumpin, swhmask)
            if (swcalib) call calibrate(azpark, elpark, azoff, eloff)
          endif
        endif
      enddo
      stop
      end

      subroutine run_tle(linetime,linedoy, ro, vo, long, lat, theta0g, doy_tle, iyr_tle, iyr_now)
      implicit none
      integer*4 linedoy, iyr_tle, iyr_now, diy
      real*8 linetime, ro(3), vo(3), long, lat, doy_tle
c
      real*8    Tmfe, theta0g
      save
      diy = 0
      if (iyr_now.ne.iyr_tle) then
        diy = 365
        if (mod(iyr_tle,4).eq.0.and.iyr_tle.ne.2000) diy = diy + 1
      endif
      Tmfe = ((linetime/86400.0D0) + dble(linedoy + diy) - doy_tle) * 1440.0D0
      call sgp4_run(Tmfe, ro, vo, long, lat, theta0g)
      return
      end

c-- For some reason inquire does not work so I implement here a bit of a shaky replacement
      subroutine exist(file, doesexist)
      implicit none
      character*(*) file
      logical       doesexist
c
      integer*4 status
      status    = -1
      doesexist = .false.
      call system ('ls '//file(1:lnblnk(file))//' > test.out', status)
      if (status.eq.0) doesexist = .true.
      call system('rm test.out ')
      return
      end

      subroutine azel(long, lat, ro, az, el, swazelalt, cloc)
      implicit none
      real*8      long, lat, ro(3), az, el
      logical     swazelalt
      character*3 cloc
c
      real*8 xyz(3), lla(3), xyzsat(3), xyzobs(3), a(3), b(3), c1000, torad, todeg, a_par_b(3), a_ortho_b(3), const
      real*8 mylong, mylat
      integer*2 npos
      character*50 instring
c-- Test variables - declare separately to be able to remove quickly
      real*8 t_e, t_n, t_u, t_az, t_el, t_re
c
      c1000           = 1000.0D0
      torad           = 2.0D0 * dasin(1.0D0) / 180.0D0
      todeg           = 180.0D0 / (2.0D0 * dasin(1.0D0))
      if (cloc.eq.'LOC') then
        call getenv('LONG',instring)
        call string_to_r8(instring(1:lnblnk(instring)), npos, mylong)
        call getenv('LAT',instring)
        call string_to_r8(instring(1:lnblnk(instring)), npos, mylat)
      else
c- Co-analyse az/el for Svallbard
        if (cloc.eq.'SVB') then
          mylong = 15.392766
          mylat  = 78.230457
        endif
c- Co-analyse az/el for Kiruna
        if (cloc.eq.'KIR') then
          mylong = 20.966881
          mylat  = 67.858429
        endif
      endif
c
      xyz(1) = ro(1) * c1000
      xyz(2) = ro(2) * c1000
      xyz(3) = ro(3) * c1000
      call xyz2lla(xyz, lla)
      lla(1) = long * torad
      lla(2) = lat  * torad
      call lla2xyz(lla, xyzsat)
      lla(1) =  mylong * torad
      lla(2) =  mylat  * torad
      lla(3) =  0.0D0
      call lla2xyz(lla, xyzobs)
c
      b(1) = xyzobs(1)
      b(2) = xyzobs(2)
      b(3) = xyzobs(3)
      a(1) = xyzsat(1) - b(1)
      a(2) = xyzsat(2) - b(2)
      a(3) = xyzsat(3) - b(3)
c-- Code to test a different approach to getting Az, El from Soler et al - see d:\Tools\Documents\
      if (.not.swazelalt) then
        t_e = -1.0D0 * a(1) * sin(mylong * torad)                       + a(2) * cos(mylong * torad)
        t_n = -1.0D0 * a(1) * sin(mylat  * torad) * cos(mylong * torad) - a(2) * sin(mylat  * torad) * sin(mylong * torad) + a(3) * cos(mylat * torad)
        t_u =          a(1) * cos(mylat  * torad) * cos(mylong * torad) + a(2) * cos(mylat  * torad) * sin(mylong * torad) + a(3) * sin(mylat * torad)
        t_az = todeg * atan(t_e / t_n)
        t_el = todeg * atan(t_u / sqrt(t_e * t_e + t_n * t_n))
c-- Atmospheric refraction correction - for 1010 hPa, 10 deg C and yellow light - is this valid for L- and X-band ????
c-- Only execute the refraction correction for t_el > 0 !!
        if (t_el.gt.0.0D0) then
          t_re = 1.02D0 /tan(torad * (t_el + (10.3D0/(t_el + 5.11D0))))
          t_re = t_re / 60.00D0
          t_el = t_el + t_re
        endif
        if (t_n.le.0.0D0) t_az = 180.0D0 + t_az
        if (t_n.gt.0.0D0) t_az = 360.0D0 + t_az
        t_az = mod(t_az, 360.0D0)
        az = t_az
        el = t_el
c--
      else
c-- My own approach to calculating Az/El - gives a slightly different result than the Soler et al - which is correct ?
        el = 90.0D0 - todeg * dacos( (a(1)*b(1)+a(2)*b(2)+a(3)*b(3)) / (dsqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3)) * dsqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3)) ) )
c-- Assuming the refraction correction is correct - also for X-band - you could also add it here ?
        const = (a(1) * b(1) + a(2) * b(2) + a(3) * b(3)) / (b(1) * b(1) + b(2) * b(2) + b(3) * b(3))
        a_par_b(1)   = const * b(1)
        a_par_b(2)   = const * b(2)
        a_par_b(3)   = const * b(3)
        a_ortho_b(1) = a(1) - a_par_b(1)
        a_ortho_b(2) = a(2) - a_par_b(2)
        a_ortho_b(3) = a(3) - a_par_b(3)
        a(1) = 0.0D0 - b(1)
        a(2) = 0.0D0 - b(2)
        a(3) = ( (b(1)*b(1) + b(2)*b(2) + b(3)*b(3)) / b(3)) - b(3)
        b(1) = a_ortho_b(1)
        b(2) = a_ortho_b(2)
        b(3) = a_ortho_b(3)
        az = todeg * dacos( (a(1)*b(1)+a(2)*b(2)+a(3)*b(3)) / (dsqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3)) * dsqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3)) ) )
        if (long.lt.lla(1)*todeg) az = 360.0D0 - az
      endif
c
      return
      end

      subroutine azel_plot(az, el, i1, i2, n, ihr1, imn1, rsc1, ihr2, imn2, rsc2, swhmask, swlive, csat, tstart, swns, co, swcovis, cloc)
      implicit none
      logical   swhmask, swlive, swns, swcovis
      integer*4 i1, i2, n, ihr1, imn1, ihr2, imn2
      real*8    az(n), el(n), rsc1, rsc2, tstart
      integer*1 co(n)
      character*3 cloc
      character*(*) csat
c
      logical        check_horizon
      integer*1      R, G, B
      integer*4      dattim(8), i, doy
      character*1    csq, cdq
      character*5    c_zone
      character*6    c_time1, c_time2, cns
      character*8    c_date, c_time3, c_time4
      character*10   c_time, c_date1
      character*250  tledir, cstring
      character*1000 command
c
      call date_and_time(c_date, c_time, c_zone, dattim)
      if (tstart.gt.86400.0D0) then
        call datetodoy(dattim(1), dattim(2), dattim(3), doy)
        doy = doy + int(tstart/86400.0D0)
        call doytodate(doy, dattim(1), c_date)
      endif
      c_date1 = c_date(1:4)//'-'//c_date(5:6)//'-'//c_date(7:8)
      write (c_time1(1:2),'(I2.2)') ihr1
      write (c_time1(3:4),'(I2.2)') imn1
      write (c_time1(5:6),'(I2.2)') nint(rsc1)
      write (c_time2(1:2),'(I2.2)') ihr2
      write (c_time2(3:4),'(I2.2)') imn2
      write (c_time2(5:6),'(I2.2)') nint(rsc2)
      c_date1 = c_date(1:4)//'-'//c_date(5:6)//'-'//c_date(7:8)
      write (c_time3,'(I2.2,'':'',I2.2,'':'',I2.2)') ihr1, imn1, nint(rsc1)
      write (c_time4,'(I2.2,'':'',I2.2,'':'',I2.2)') ihr2, imn2, nint(rsc2)
      call getenv('TLE',tledir)
c
      call map_azel_init()
      if (swhmask) then
        r = 'FF'x
        g = '00'x
        b = '00'x
        do i = 1, min(max(1, i1 - 1), n)
          call map_azel_place(az(i), el(i), r, g, b)
        enddo
        do i = min(max(1, i1), n), max(min(i2, n), 1)
          if (check_horizon(az(i), el(i))) then
            r = '00'x
            g = 'FF'x
            b = '00'x
            if (swcovis.and.co(i).eq.1) then
              r = 'FF'x
              g = '00'x
              b = 'FF'x
            endif
            call map_azel_place(az(i), el(i), r, g, b)
          else
            r = 'FF'x
            g = '00'x
            b = '00'x
            call map_azel_place(az(i), el(i), r, g, b)
          endif
        enddo
        r = 'FF'x
        g = 'FF'x
        b = '00'x
        do i = max(min(i2 + 1, n), 1), n
          call map_azel_place(az(i), el(i), r, g, b)
        enddo
        call map_azel_write()
        csq = char(ichar("'"))
        cdq = char(ichar('"'))
        cns = 'S => N'
        if (swns) cns = 'N => S'
        cstring = '-pointsize 20 -fill red -draw '//csq//'text 20,20 '//cdq//c_date1//' '//c_time3//' '//c_time4//cdq//csq
        call system('convert satmap-azel.png '//cstring(1:lnblnk(cstring))//' satmap-azel.png')
        cstring = '-pointsize 20 -fill red -draw '//csq//'text 20,50 '//csat(1:lnblnk(csat))//csq
        call system('convert satmap-azel.png '//cstring(1:lnblnk(cstring))//' satmap-azel.png')
        cstring = '-pointsize 20 -fill red -draw '//csq//'text 20,75 '//cdq//cns//cdq//csq
        call system('convert satmap-azel.png '//cstring(1:lnblnk(cstring))//' satmap-azel.png')
        if (swcovis) then
          cstring = '-pointsize 20 -fill red -draw '//csq//'text 20,100 '//cdq//'covis='//cloc//cdq//csq
          call system('convert satmap-azel.png '//cstring(1:lnblnk(cstring))//' satmap-azel.png')
        endif
        if (swlive) then
          do i = lnblnk(csat), 1, -1
            if (csat(i:i).eq.' ') csat(i:i) = '_'
          enddo
          write (command,'(''cp satmap-azel.png '',a,a,''_'',a,''_'',a,''_azel_live_'',a,''.png'')') tledir(1:lnblnk(tledir)), c_date, c_time1, c_time2, csat(1:lnblnk(csat))
        else
          do i = lnblnk(csat), 1, -1
            if (csat(i:i).eq.' ') csat(i:i) = '_'
          enddo
          write (command,'(''cp satmap-azel.png '',a,a,''_'',a,''_'',a,''_azel_predict_'',a,''.png'')') tledir(1:lnblnk(tledir)), c_date, c_time1, c_time2, csat(1:lnblnk(csat))
        endif
        call system(command)
      else
        do i = min(max(1, i1), n), max(min(i2, n), 1)
          r = '00'x
          g = 'FF'x
          b = '00'x
          if (swcovis.and.co(i).eq.1) then
            r = 'FF'x
            g = '00'x
            b = 'FF'x
          endif
          call map_azel_place(az(i), el(i), r, g, b)
        enddo
c-- To be able to establish the direction of the pass when not using hmask, label the first point red and the last one yellow
        r = 'FF'x
        g = '00'x
        b = '00'x
        call map_azel_place(az(1), el(1), r, g, b)
        r = 'FF'x
        g = 'FF'x
        b = '00'x
        call map_azel_place(az(n), el(n), r, g, b)
        call map_azel_write()
        csq = char(ichar("'"))
        cdq = char(ichar('"'))
        cns = 'S => N'
        if (swns) cns = 'N => S'
        cstring = '-pointsize 20 -fill red -draw '//csq//'text 70,50 '//csat(1:lnblnk(csat))//csq
        call system('convert satmap-azel.png '//cstring(1:lnblnk(cstring))//' satmap-azel.png')
        cstring = '-pointsize 20 -fill red -draw '//csq//'text 70,75 '//cdq//cns//cdq//csq
        call system('convert satmap-azel.png '//cstring(1:lnblnk(cstring))//' satmap-azel.png')
        if (swcovis) then
          cstring = '-pointsize 20 -fill red -draw '//csq//'text 70,100 '//cdq//'covis='//cloc//cdq//csq
          call system('convert satmap-azel.png '//cstring(1:lnblnk(cstring))//' satmap-azel.png')
        endif
        if (swlive) then
          write (command,'(''cp satmap-azel.png '',a,a,''_'',a,''_'',a,''_azel_live_'',a,''.png'')') tledir(1:lnblnk(tledir)), c_date, c_time1, c_time2, csat(1:lnblnk(csat))
        else
          write (command,'(''cp satmap-azel.png '',a,a,''_'',a,''_'',a,''_azel_predict_'',a,''.png'')') tledir(1:lnblnk(tledir)), c_date, c_time1, c_time2, csat(1:lnblnk(csat))
        endif
        call system(command)
      endif
c
      return
      end

      subroutine map_azel_init()
      implicit none
c--
      integer*4 nazel
      parameter (nazel=1400)
c--
      integer*1 mapsline (3,nazel), maps(3,nazel,nazel), r, g, b
      integer*4 lun, i, j, ix, iy
      real*8    az, el, torad, todeg
      character*6 cmap
      character*11 csize
      character*500 cstring
c
      save
c
      cmap    = 'satmap' 
      torad   = 2.0D0 * dasin(1.0D0) / 180.0
      todeg   = 180.0D0 / 2.0D0 / dasin(1.0D0)
      write (csize,'(I5.5,''x'',I5.5)') nazel, nazel
c
      do i = 1, nazel
        do j = 1, nazel
          maps(1,j,i) = 0
          maps(2,j,i) = 0
          maps(3,j,i) = 0
        enddo
      enddo
c-- Az-El test
      do i = 1, nazel/2 - 1
        j  = (nazel/2) - nint(float(i) / tan(asin(2.0*float(i)/float(nazel))))
        ix = j + 1
        iy = (nazel/2) - i
        ix = min(max(ix,1),nazel)
        iy = min(max(iy,1),nazel)
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
        iy = (nazel/2) + i
        iy = min(max(iy,1),nazel)
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
        j  = (nazel/2) + nint(float(i) / tan(asin(2.0*float(i)/float(nazel))))
        ix = j + 1
        iy = (nazel/2) - i
        ix = min(max(ix,1),nazel)
        iy = min(max(iy,1),nazel)
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
        iy = (nazel/2) + i
        iy = min(max(iy,1),nazel)
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
      enddo
      do i = 1, nazel
        iy = i
        ix = nazel / 2
        ix = min(max(ix,1),nazel)
        iy = min(max(iy,1),nazel)
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
      enddo
      do i = 1, nazel
        ix = i
        iy = nazel / 2
        ix = min(max(ix,1),nazel)
        iy = min(max(iy,1),nazel)
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
      enddo
      do i = 1, 720
        ix = (nazel/2) + int(dsin(dble(i/2)*torad)*dble(nazel/6))
c iy is counterintuitive as the image is vertically flipped for display !
        iy = (nazel/2) + int(dcos(dble(i/2)*torad)*dble(nazel/6))
        ix = min(max(ix,1),nazel)
        iy = min(max(iy,1),nazel)
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
      enddo
      do i = 1, 720
        ix = (nazel/2) + int(dsin(dble(i/2)*torad)*dble(nazel/3))
c iy is counterintuitive as the image is vertically flipped for display !
        iy = (nazel/2) + int(dcos(dble(i/2)*torad)*dble(nazel/3))
        ix = min(max(ix,1),nazel)
        iy = min(max(iy,1),nazel)
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
      enddo
      return
c
      entry map_azel_place(az, el, r, g, b)
      if (el.gt.0.0) then
        ix = (nazel/2) + int(dsin(az*torad)*dble(nazel/2)*(1.0D0-(el/90.0D0)))
c iy is counterintuitive as the image is vertically flipped for display !
        iy = (nazel/2) + int(dcos(az*torad)*dble(nazel/2)*(1.0D0-(el/90.0D0)))
        ix = min(max(ix,1),nazel)
        iy = min(max(iy,1),nazel)
        do i = ix - 1, ix + 1
          do j = iy - 1, iy + 1
            maps(1,i,j) = r
            maps(2,i,j) = g
            maps(3,i,j) = b
          enddo
        enddo
      endif
      return
c
      entry map_azel_write()
      open (unit=lun, file=cmap//'-azel.rgb', access='direct', form='unformatted', recl=3*nazel)
      do i = 1, nazel
        do j = 1, nazel
          mapsline(1,j) = maps(1,j,i)
          mapsline(2,j) = maps(2,j,i)
          mapsline(3,j) = maps(3,j,i)
        enddo
        write (lun,rec=i) mapsline
      enddo
      close (unit=lun)
      call free_lun(lun)
c      write (*,'(a)') ' convert -depth 8 -size '//csize//' rgb:'//cmap//'-azel.rgb -flip '//' '//cmap//'-azel.png'
      call system('convert -depth 8 -size '//csize//' rgb:'//cmap//'-azel.rgb -flip '//' '//cmap//'-azel.png')
      return
      end

      subroutine doytodate(linedoy, iyr, cdate)
      implicit none
      integer*4 linedoy, iyr
      character*8 cdate
c
      integer*4    diy, monthdays(12), i, j
      character*5  c_zone
      character*8  c_date
      character*10 c_time
c
      data monthdays/31,0,31,30,31,30,31,31,30,31,30,31/
c
      save
c
      cdate = '        '
      diy = 365
      if (mod(iyr,4).eq.0.and.iyr.ne.2000) then
        diy = diy + 1
        monthdays(2) = 29
      else
        monthdays(2) = 28
      endif
      j = 0
      do i = 1,12
        if (j + monthdays(i).ge.linedoy) then
          write (cdate,'(I4.4,I2.2,I2.2)') iyr, i, linedoy - j
          goto 1
        endif
        j = j + monthdays(i)
      enddo
 1    continue
      return
      end
      
      subroutine datetodoy(iyr, imn, ida, doy)
      implicit none
      integer*4 iyr, imn, ida, doy
c
      integer*4 monthdays(12), i
c
      data monthdays/31,0,31,30,31,30,31,31,30,31,30,31/
c
      save
c
      if (mod(iyr,4).eq.0.and.iyr.ne.2000) then
        monthdays(2) = 29
      else
        monthdays(2) = 28
      endif
      doy = 0
      do i = 1, imn - 1
        doy = doy + monthdays(i)
      enddo
      doy = doy + ida
      return
      end

      subroutine offset_time(dattim, isec_off_in)
      integer*4 dattim(8), isec_off_in
c-- Offset the system time array by a fixed amount in the range between -86399 and +86399 seconds
      integer*4 ihr_off, imn_off, monthdays(12), isec_off
c
      data monthdays/31,0,31,30,31,30,31,31,30,31,30,31/
c
      save
c
      isec_off = isec_off_in
      if (isec_off.lt.0) goto 2
      if (isec_off.ge.3600) then
        ihr_off  = isec_off / 3600
        isec_off = mod(isec_off, 3600)
        if (ihr_off.ge.24) stop 'subroutine offset_time: offsets >= 86400 seconds not supported'
        dattim(5) = dattim(5) + ihr_off
        if (dattim(5).ge.24) then
          dattim(5) = mod(dattim(5), 24)
          if (mod(dattim(1),4).eq.0.and.dattim(1).ne.2000) then
            monthdays(2) = 29
          else
            monthdays(2) = 28
          endif
          dattim(3) = dattim(3) + 1
          if (dattim(3).gt.monthdays(dattim(2))) then
            dattim(3) = mod(dattim(3),monthdays(dattim(2)))
            dattim(2) = dattim(2) + 1
            if (dattim(2).gt.12) then
              dattim(2) = mod(dattim(2),12)
              dattim(1) = dattim(1) + 1
              if (mod(dattim(1),4).eq.0.and.dattim(1).ne.2000) then
                monthdays(2) = 29
              else
                monthdays(2) = 28
              endif
            endif
          endif
        endif
      endif
      if (isec_off.ge.60) then
        imn_off = isec_off / 60
        isec_off = mod(isec_off, 60)
        dattim(6) = dattim(6) + imn_off
        if (dattim(6).ge.60) then
          dattim(6) = mod(dattim(6), 60)
          dattim(5) = dattim(5) + 1
          if (dattim(5).ge.24) then
            dattim(5) = mod(dattim(5), 24)
            if (mod(dattim(1),4).eq.0.and.dattim(1).ne.2000) then
              monthdays(2) = 29
            else
              monthdays(2) = 28
            endif
            dattim(3) = dattim(3) + 1
            if (dattim(3).gt.monthdays(dattim(2))) then
              dattim(3) = mod(dattim(3),monthdays(dattim(2)))
              dattim(2) = dattim(2) + 1
              if (dattim(2).gt.12) then
                dattim(2) = mod(dattim(2),12)
                dattim(1) = dattim(1) + 1
                if (mod(dattim(1),4).eq.0.and.dattim(1).ne.2000) then
                  monthdays(2) = 29
                else
                  monthdays(2) = 28
                endif
              endif
            endif
          endif
        endif
      endif
      dattim(7) = dattim(7) + isec_off
      if (dattim(7).ge.60) then
        dattim(7) = mod(dattim(7), 60)
        dattim(6) = dattim(6) + 1
        if (dattim(6).ge.60) then
          dattim(6) = mod(dattim(6), 60)
          dattim(5) = dattim(5) + 1
          if (dattim(5).ge.24) then
            dattim(5) = mod(dattim(5), 24)
            if (mod(dattim(1),4).eq.0.and.dattim(1).ne.2000) then
              monthdays(2) = 29
            else
              monthdays(2) = 28
            endif
            dattim(3) = dattim(3) + 1
            if (dattim(3).gt.monthdays(dattim(2))) then
              dattim(3) = mod(dattim(3),monthdays(dattim(2)))
              dattim(2) = dattim(2) + 1
              if (dattim(2).gt.12) then
                dattim(2) = mod(dattim(2),12)
                dattim(1) = dattim(1) + 1
                if (mod(dattim(1),4).eq.0.and.dattim(1).ne.2000) then
                  monthdays(2) = 29
                else
                  monthdays(2) = 28
                endif
              endif
            endif
          endif
        endif
      endif
      goto 3
c-- In a cowards way, I deal with negative offsets in a separate way
 2    isec_off = abs(isec_off)
      if (isec_off.ge.3600) then
        ihr_off  = isec_off / 3600
        isec_off = mod(isec_off, 3600)
        if (ihr_off.ge.24) stop 'subroutine offset_time: offsets >= 86400 seconds not supported'
        dattim(5) = dattim(5) - ihr_off
        if (dattim(5).lt.0) then
          dattim(5) = dattim(5) + 24
          if (mod(dattim(1),4).eq.0.and.dattim(1).ne.2000) then
            monthdays(2) = 29
          else
            monthdays(2) = 28
          endif
          dattim(3) = dattim(3) - 1
          if (dattim(3).le.0) then
            dattim(2) = dattim(2) - 1
            if (dattim(2).le.0) then
              dattim(2) = dattim(2) + 12
              dattim(3) = dattim(3) + monthdays(dattim(2))
              dattim(1) = dattim(1) - 1
              if (mod(dattim(1),4).eq.0.and.dattim(1).ne.2000) then
                monthdays(2) = 29
              else
                monthdays(2) = 28
              endif
            else
              dattim(3) = dattim(3) + monthdays(dattim(2))
            endif
          endif
        endif
      endif
      if (isec_off.ge.60) then
        imn_off = isec_off / 60
        isec_off = mod(isec_off, 60)
        dattim(6) = dattim(6) - imn_off
        if (dattim(6).lt.0) then
          dattim(6) = dattim(6) + 60
          dattim(5) = dattim(5) - 1
          if (dattim(5).lt.0) then
            dattim(5) = dattim(5) + 24
            if (mod(dattim(1),4).eq.0.and.dattim(1).ne.2000) then
              monthdays(2) = 29
            else
              monthdays(2) = 28
            endif
            dattim(3) = dattim(3) - 1
            if (dattim(3).le.0) then
              dattim(2) = dattim(2) - 1
              if (dattim(2).le.0) then
                dattim(2) = dattim(2) + 12
                dattim(3) = dattim(3) + monthdays(dattim(2))
                dattim(1) = dattim(1) - 1
                if (mod(dattim(1),4).eq.0.and.dattim(1).ne.2000) then
                  monthdays(2) = 29
                else
                  monthdays(2) = 28
                endif
              else
                dattim(3) = dattim(3) + monthdays(dattim(2))
              endif
            endif
          endif          
        endif
      endif
      dattim(7) = dattim(7) - isec_off
      if (dattim(7).lt.0) then
        dattim(7) = dattim(7) + 60
        dattim(6) = dattim(6) - 1
        if (dattim(6).lt.0) then
          dattim(6) = dattim(6) + 60
          dattim(5) = dattim(5) - 1
          if (dattim(5).lt.0) then
            dattim(5) = dattim(5) + 24
            if (mod(dattim(1),4).eq.0.and.dattim(1).ne.2000) then
              monthdays(2) = 29
            else
              monthdays(2) = 28
            endif
            dattim(3) = dattim(3) - 1
            if (dattim(3).le.0) then
              dattim(2) = dattim(2) - 1
              if (dattim(2).le.0) then
                dattim(2) = dattim(2) + 12
                dattim(3) = dattim(3) + monthdays(dattim(2))
                dattim(1) = dattim(1) - 1
                if (mod(dattim(1),4).eq.0.and.dattim(1).ne.2000) then
                  monthdays(2) = 29
                else
                  monthdays(2) = 28
                endif
              else
                dattim(3) = dattim(3) + monthdays(dattim(2))
              endif
            endif
          endif          
        endif
      endif
 3    return
      end

      subroutine get_lun(ilun)
      implicit none
      integer*4 ilun,lun(100), i
      byte      occupado(100)
      data lun/  1,  2,  3,  4,  7,  8,  9, 10, 11, 12,
     *          13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
     *          23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
     *          33, 34, 35, 36, 37, 38, 39, 40, 41, 42,
     *          43, 44, 45, 46, 47, 48, 49, 50, 51, 52,
     *          53, 54, 55, 56, 57, 58, 59, 60, 61, 62,
     *          63, 64, 65, 66, 67, 68, 69, 70, 71, 72,
     *          73, 74, 75, 76, 77, 78, 79, 80, 81, 82,
     *          83, 84, 85, 86, 87, 88, 89, 90, 91, 92,
     *          93, 94, 95, 96, 97, 98, 99,100,101,102/
      data occupado/0,0,0,0,0,0,0,0,0,0,
     *              0,0,0,0,0,0,0,0,0,0,
     *              0,0,0,0,0,0,0,0,0,0,
     *              0,0,0,0,0,0,0,0,0,0,
     *              0,0,0,0,0,0,0,0,0,0,
     *              0,0,0,0,0,0,0,0,0,0,
     *              0,0,0,0,0,0,0,0,0,0,
     *              0,0,0,0,0,0,0,0,0,0,
     *              0,0,0,0,0,0,0,0,0,0,
     *              0,0,0,0,0,0,0,0,0,0/
c
      save
c
      ilun = -1
      i    = 1
      do while (i.le.100.and.occupado(i).ne.0)
        i = i + 1
      enddo
      if (i.le.100) then
        ilun        = lun(i)
	occupado(i) = 1
      endif
      return
c
      entry free_lun(ilun)
      do i = 1,100
        if (lun(i).eq.ilun) then
          ilun        = -1
          occupado(i) = 0
        endif
      enddo
      return
      end

c-- Code borrowed from : https://ea4eoz.blogspot.com/2015/11/simple-wgs-84-ecef-conversion-functions.html?m=1
      subroutine lla2xyz(lla,xyz)
      implicit none
      real*8 lla(3), xyz(3)
c
      real*8 a, e, asq, esq, lon, lat, alt, N
      a = 6378137.0D0
      e = 8.1819190842622D-2
      asq = a * a
      esq = e * e
      lon = lla(1)
      lat = lla(2)
      alt = lla(3)
      N   = a / dsqrt(1.0D0 - esq * sin(lat) * sin(lat))
      xyz(1) = (N + alt) * dcos(lat) * dcos(lon)
      xyz(2) = (N + alt) * dcos(lat) * dsin(lon)
      xyz(3) = ((1.0D0 - esq) * N + alt) * dsin(lat)
c
      return
      end

      subroutine xyz2lla(xyz,lla)
      implicit none
      real*8 xyz(3), lla(3)
c
      real*8 a, e, asq, esq, x, y, z, b, bsq, ep, p, th, lon, lat, alt, N, g(3), gm, am
      a   = 6378137.0D0
      e   = 8.1819190842622D-2
      asq = a * a
      esq = e * e
      x   = xyz(1)
      y   = xyz(2)
      z   = xyz(3)
      b   = dsqrt(asq * (1.0D0 - esq))
      bsq = b * b
      ep  = dsqrt(( asq - bsq) / bsq)
      p   = dsqrt(x * x + y * y)
      th  = datan2(a * z, b * p)
      lon = datan2(y, x)
      lat = datan2((z + ep * ep * b * dsin(th) * dsin(th) * dsin(th)),(p - esq * a * dcos(th) * dcos(th) * dcos(th)))
      N   = a / (dsqrt(1.0D0 - esq * dsin(lat) * dsin(lat)))
      lla(1) = lon
      lla(2) = lat
      lla(3) = 0.0D0
      call lla2xyz(lla, g)
      gm  = dsqrt(g(1) * g(1) + g(2) * g(2) + g(3) * g(3))
      am  = dsqrt(x * x + y * y + z * z)
      alt = am - gm
c-- Second iteration suggested on web page
      lla(1) = lon
      lla(2) = lat
      lla(3) = alt
      call lla2xyz(lla, g)
      gm  = dsqrt(g(1) * g(1) + g(2) * g(2) + g(3) * g(3))
      alt = alt + (am - gm)
      lla(1) = lon
      lla(2) = lat
      lla(3) = alt
c
      return
      end

      subroutine string_to_i4(string,npos,value)
      integer*2 npos,ilen
      integer*4 value
      character*(*) string
      character*10  cfmt
      ilen = len(string)
      do 1 j = 1,ilen
        if (string(j:j).ne.' ')then
          do 2 i = j,ilen
            if (string(i:i).eq.' ') then
              npos = i + 1
              write (cfmt,'(''(I'',i1,'')'')') i-j
              goto 3
            endif
 2        continue
        endif
 1    continue
 3    read (string(j:i-1),cfmt) value
      return
      end

      subroutine string_to_r8(string,npos,value)
      logical*1 swexp
      integer*2 npos,ilen,i_mantisse,i_main
      real*8    value
      character*1 cf
      character*10 cfmt
      character*(*) string
      ilen = len(string)
      i1 = 1
      i2 = ilen
      do 1 i = 1,ilen
        if (string(i:i).ne.' ') then
          i1 = i
          do 2 j = i,ilen
            if (string(j:j).eq.' ') then
              i2 = j - 1
              goto 3
            endif
 2       continue
         i2 = ilen
         goto 3
       endif
 1    continue
 3    npos = i2 + 1
      i_exp = 0
      swexp = .false.
      cf    = 'F'
      do 4 i = i2,i1,-1
        if (string(i:i).eq.'e'.or.string(i:i).eq.'E') then
          swexp = .true.
          i_exp = i
          cf    = 'E'
          goto 5
        else if (string(i:i).eq.'d'.or.string(i:i).eq.'D') then
          swexp = .true.
          i_exp = i
          cf    = 'D'
          goto 5
        endif
 4    continue
 5    if (swexp) then
        i_mantisse = 0
        i_main     = 0
        do 8 i = i_exp-1,i1,-1
          if (string(i:i).eq.'.') then
            i_mantisse = i_exp - 1 - i 
            i_main     = i2 - i1 + 1
            goto 7
          endif
 8      continue
        i_mantisse = 0
        i_main     = i2 - i1 + 1
      else
        i_mantisse = 0
        i_main     = 0
        do 6 i = i2,i1,-1
          if (string(i:i).eq.'.') then
            i_mantisse = i2 - i 
            i_main     = i2 - i1 + 1
            goto 7
          endif
 6      continue
        i_mantisse = 0
        i_main     = i2 - i1 + 1
      endif
 7    if (i_main.ge.10) then
        if (i_mantisse.ge.10) then
          write (cfmt,'(''('',A1,I2,''.'',I2,'')'')') 
     *           cf,i_main,i_mantisse
        else
          write (cfmt,'(''('',A1,I2,''.'',I1,'')'')') 
     *           cf,i_main,i_mantisse
        endif
      else
        if(i_mantisse.ge.10) then
          write (cfmt,'(''('',A1,I1,''.'',I2,'')'')') 
     *           cf,i_main,i_mantisse
        else
          write (cfmt,'(''('',A1,I1,''.'',I1,'')'')') 
     *           cf,i_main,i_mantisse
        endif
      endif
      read (string(i1:i2),cfmt) value
      return
      end

      subroutine qsort_r8(sort,key,nsort)
      integer   nsort
      real*8    sort(nsort)
      integer   key(nsort)
c
      integer   st(2,17),l,r,s,dbuf,lmed
      integer   i, j
      real*8    rbuf,dmed
      s       = 1
      st(1,s) = 1
      st(2,s) = nsort
 1    l = st(1,s)
      r = st(2,s)
      s = s - 1
      if (l.ge.r) goto 12
 2    i    = l
      j    = r
      lmed = (l+r) / 2
      dmed = sort(lmed)
      if (i.gt.j) goto 15
 5    if (sort(i).ge.dmed) goto 13
 3    i = i + 1
      if (sort(i).lt.dmed) goto 3
 13   if (sort(j).le.dmed) goto 14
 4    j = j - 1
      if (sort(j).gt.dmed) goto 4
 14   if (i.le.j) then
         rbuf    = sort(i)
         sort(i) = sort(j)
         sort(j) = rbuf
         dbuf    = key(i)
         key(i)  = key(j)
         key(j)  = dbuf
         i       = i + 1
         j       = j - 1
      endif
      if (i.le.j) goto 5
 15   if((j-l) .lt. (r-i)) then
         if(i.lt.r) then
            s       = s + 1
            st(1,s) = i
            st(2,s) = r
         endif
         r = j
      else
         if (l.lt.j) then
            s       = s + 1
            st(1,s) = l
            st(2,s) = j
         endif
         l = i
      endif
      if (l.lt.r) goto 2
 12   if (s.ne.0) goto 1
      return
      end

      subroutine azel_command(pass_start, pass_end, time_track, sat_az_track, sat_el_track, n_in, swdebug, swallowflip, swnoovl, ovlflag, azrate, elrate, azoff, eloff, swoffset, azpark, elpark, delta_t_la, swjumpin, swhmask)
      implicit none
      integer*4 n_in
      real*8    pass_start, pass_end, time_track(n_in), sat_az_track(n_in), sat_el_track(n_in), azrate, elrate, azoff, eloff, azpark, elpark, delta_t_la
      integer*1 ovlflag(n_in)
      logical   swdebug, swallowflip, swnoovl, swoffset, swjumpin, swhmask
c
      logical   swflip, check_horizon, visible
      real*8    rsc, time_now, time_sleep, time_start, az1, az2, el1, el2, t_end_move
      integer*4 n, i, j, ihr, imn, dattim(8), doy, i_check, i_mask1, i_mask2, i_start
      character*5  c_zone
      character*7  c_visstr
      character*8  c_date
      character*10 c_time
      character*50 c_command
c
      save
c
      n = n_in
c-- If Horizon mask is to be applied check the range here
      if (swhmask) then
        i_mask1 = -1
        i_mask2 = -1
        do j = 1, n
          if (check_horizon(sat_az_track(j), sat_el_track(j)).and.i_mask1.lt.0) i_mask1 = j
        enddo
        do j = n, 1, -1
          if (check_horizon(sat_az_track(j), sat_el_track(j)).and.i_mask2.lt.0) i_mask2 = j
        enddo
        if (i_mask1.eq.-1.or.i_mask2.eq.-1) then
          write (*,'('' No El commands above Horizon mask - skipping pass '')')
          return
        endif
        write (*,'('' Commands : '',I5, '' through : '',I5,'' delimit the window with El above Horizon mask (there may be gaps in between)'')') i_mask1, i_mask2
      else
        i_mask1 = 1
        i_mask2 = n
      endif
c-- Apply offset if needed
      if (swoffset) then
        do j = 1, n
          sat_az_track(j) = sat_az_track(j) + azoff
          sat_el_track(j) = sat_el_track(j) + eloff
        enddo
      endif
c-- Add delta_t_la (=0.0D0 unless set on command line with deltat=)
      do j = 1, n
        time_track(j) = time_track(j) - delta_t_la
      enddo
c-- Check if Flip mode required
      swflip = .false.
      do j = 1, n-1
        if (abs(sat_az_track(j+1) - sat_az_track(j)).gt.250.0D0.and.swallowflip) then
          swflip = .true.
          write (*,'(''Using FLIP mode !'')')
        endif
      enddo
c-- Check overlapping commands if requested - the algorithm works but is a bit flaky. After a large Az turn eg you skip say 6 commands and then the first 
c-- OK command is again a fairly large Az-El rotation and this currently is not checked. This would require a second loop where i+1 is replaced by the first 
c-- index after i with ovlflag=0 - I think I have fixed the problem by a double loop and a collapse of the to be ignored commands.
      if (swnoovl) then
        i_check = 0
 3      i_check = i_check + 1
        do i = 1, n-2
          if (ovlflag(i).eq.0) then
            az1 = sat_az_track(i)
            az2 = sat_az_track(i+1)
            el1 = sat_el_track(i)
            el2 = sat_el_track(i+1)
            az1 = min(max(az1,0.0D0),360.0D0)
            if (swflip) az1 = mod((az1 + 180.0D0),360.0D0)
            az2 = min(max(az2,0.0D0),360.0D0)
            if (swflip) az2 = mod((az2 + 180.0D0),360.0D0)
            el1 = min(max(el1,0.0D0), 90.0D0)
            if (swflip) el1 = 180.0D0 - el1
            el2 = min(max(el2,0.0D0), 90.0D0)
            if (swflip) el2 = 180.0D0 - el2
            t_end_move = time_track(i) + max(abs(az2-az1)/azrate,abs(el2-el1)/elrate)
            do j = i+2, n
              if (time_track(j).lt.t_end_move) ovlflag(j) = 1
            enddo
          endif
        enddo
c-- Collapse the commands to only leave those with ovlflag=0
        j = 0
        do i = 1,n
          if (ovlflag(i).eq.0) then
            j = j + 1
            if (i.ne.j) then
              time_track(j)   = time_track(i)
              sat_az_track(j) = sat_az_track(i)
              sat_el_track(j) = sat_el_track(i)
              ovlflag(j)      = ovlflag(i)
            endif
          endif
        enddo
        if (swdebug) write (*,'('' Initial # of commands : '',I5,'' # remaining after overlap removal : '',I5)') n, j
        n = j
c-- Do the time check loop twice to try and fix the problem spotted in the comment above
        if (i_check.eq.1) goto 3
c        do j = 1,n
c          ihr = int(time_track(j)/3600.0D0)
c          imn = (int(time_track(j) - dble(ihr)*3600.0D0)/60.0D0)
c          rsc = time_track(j) - dble(ihr)*3600.0D0 - imn * 60.0D0
c          write (*,'(I5,3x,2F10.3,2I7,F10.3,I5)') j, sat_az_track(j), sat_el_track(j), ihr, imn, rsc, ovlflag(j)
c        enddo
c        read (*,*)
c
c-- If there is a swnoovl and an swhmask we need to check the hmask again, after the collapse !
        if (swhmask) then
          i_mask1 = -1
          i_mask2 = -1
          do j = 1, n
            if (swoffset) then
              if (check_horizon(sat_az_track(j) - azoff, sat_el_track(j) - eloff).and.i_mask1.lt.0) i_mask1 = j
            else
              if (check_horizon(sat_az_track(j), sat_el_track(j)).and.i_mask1.lt.0) i_mask1 = j
            endif
          enddo
          do j = n, 1, -1
            if (swoffset) then
              if (check_horizon(sat_az_track(j) - azoff, sat_el_track(j) - eloff).and.i_mask2.lt.0) i_mask2 = j
            else
              if (check_horizon(sat_az_track(j), sat_el_track(j)).and.i_mask2.lt.0) i_mask2 = j
            endif
          enddo
          if (i_mask1.eq.-1.or.i_mask2.eq.-1) then
            write (*,'('' After overlap removal, there are no El commands above Horizon mask - skipping pass '')')
            return
          endif
          write (*,'('' Commands : '',I5, '' through : '',I5,'' remain after overlap removal and delimit the window with El above Horizon mask (there may be gaps in between)'')') i_mask1, i_mask2
        else
          i_mask1 = 1
          i_mask2 = n
        endif
      endif
      ihr = int(time_track(max(i_mask1-1,1))/3600.0D0)
      imn = (int(time_track(max(i_mask1-1,1)) - dble(ihr)*3600.0D0)/60.0D0)
      rsc = time_track(max(i_mask1-1,1)) - dble(ihr)*3600.0D0 - imn * 60.0D0
      write (*,'(''Pass from : '',2I5,F10.3)') ihr, imn, rsc
      ihr = int(time_track(min(i_mask2+1,n))/3600.0D0)
      imn = (int(time_track(min(i_mask2+1,n)) - dble(ihr)*3600.0D0)/60.0D0)
      rsc = time_track(min(i_mask2+1,n)) - dble(ihr)*3600.0D0 - imn * 60.0D0
      write (*,'(''       to : '',2I5,F10.3)') ihr, imn, rsc
c--
      call date_and_time(c_date, c_time, c_zone, dattim)
c      write (*,*) dattim(1), dattim(2), dattim(3), dattim(4), dattim(5), dattim(6), dattim(7), dattim(8)
      time_now = dble(dattim(5) * 3600.0 + dattim(6) * 60.0 + dattim(7) + dattim(8)/1000.0)
      time_sleep = pass_start - time_now - 60.0D0
c-- Skip ongoing passes ?
      if (time_sleep.lt.0.0.and.(.not.swjumpin)) return
      if (time_sleep.gt.0.0) then
        write (c_command,'(''sleep '',F8.1)') time_sleep
        call system(c_command)
      endif
      if (swhmask) then
        j = max(i_mask1 - 1, 1)
      else
        j = 1
      endif
      ihr = int(time_track(j)/3600.0D0)
      imn = (int(time_track(j) - dble(ihr)*3600.0D0)/60.0D0)
      rsc = time_track(j) - dble(ihr)*3600.0D0 - imn * 60.0D0
      if (swoffset) then
        visible = check_horizon(sat_az_track(j) - azoff, sat_el_track(j) - eloff)
        c_visstr = '       '
        if (visible) c_visstr = 'Visible'
        write (*,'(I5,3x,'' Commanded : '',2F10.3,2I5,F10.3,''   Actual satellite : '',2F10.3,3x,A7)') j, sat_az_track(j), sat_el_track(j), ihr, imn, rsc, sat_az_track(j) - azoff, sat_el_track(j) - eloff, c_visstr
      else
        visible = check_horizon(sat_az_track(j), sat_el_track(j))
        c_visstr = '       '
        if (visible) c_visstr = 'Visible'
        write (*,'(I5,3x,'' Commanded : '',2F10.3,2I5,F10.3,3x,A7)') j, sat_az_track(j), sat_el_track(j), ihr, imn, rsc, c_visstr
      endif
      call SPID_MD02(sat_az_track(j), sat_el_track(j), swdebug, swflip)
c      i_start = j + 1
      i_start = 2
      do j = i_start, n
        call date_and_time(c_date, c_time, c_zone, dattim)
c        write (*,*) dattim(1), dattim(2), dattim(3), dattim(4), dattim(5), dattim(6), dattim(7), dattim(8)
        time_now = dble(dattim(5) * 3600.0 + dattim(6) * 60.0 + dattim(7) + dattim(8)/1000.0)
        time_sleep = time_track(j) - time_now - 0.2
        if (time_sleep.lt.0.0) goto 2
        write (c_command,'(''sleep '',F8.1)') time_sleep
c        write (*,*) c_command(1:lnblnk(c_command))
        call system(c_command)
c-- Track real time
 1      call date_and_time(c_date, c_time, c_zone, dattim)
        time_now = dble(dattim(5) * 3600.0 + dattim(6) * 60.0 + dattim(7) + dattim(8)/1000.0)
        time_sleep = time_track(j) - time_now
        if (time_sleep.lt.0.) then
c          write (*,*) dattim(1), dattim(2), dattim(3), dattim(4), dattim(5), dattim(6), dattim(7), dattim(8)
          ihr = int(time_track(j)/3600.0D0)
          imn = (int(time_track(j) - dble(ihr)*3600.0D0)/60.0D0)
          rsc = time_track(j) - dble(ihr)*3600.0D0 - imn * 60.0D0
          if (swoffset) then
            visible = check_horizon(sat_az_track(j) - azoff, sat_el_track(j) - eloff)
            c_visstr = '       '
            if (visible) c_visstr = 'Visible'
            write (*,'(I5,3x,'' Commanded : '',2F10.3,2I5,F10.3,''   Actual satellite : '',2F10.3,3x,A7)') j, sat_az_track(j), sat_el_track(j), ihr, imn, rsc, sat_az_track(j) - azoff, sat_el_track(j) - eloff, c_visstr
          else
            visible = check_horizon(sat_az_track(j), sat_el_track(j))
            c_visstr = '       '
            if (visible) c_visstr = 'Visible'
            write (*,'(I5,3x,'' Commanded : '',2F10.3,2I5,F10.3,3x,A7)') j, sat_az_track(j), sat_el_track(j), ihr, imn, rsc, c_visstr
          endif
          if (i_mask1.le.j.and.j.le.i_mask2) then
            call SPID_MD02(sat_az_track(j), sat_el_track(j), swdebug, swflip)
          endif
          goto 2
        endif
        goto 1
 2      continue
      enddo
c-- Estimate time needed for the last move to complete and sleep that long, to prevent a M2 pulse timeout
c-- Could also decide to sleep until pass_end - but that will fail again if end elevation is set to zero in the main program
c      az1 = sat_az_track(n-1)
c      az2 = sat_az_track(n)
c      el1 = sat_el_track(n-1)
c      el2 = sat_el_track(n)
c      az1 = min(max(az1,0.0D0),360.0D0)
c      if (swflip) az1 = mod((az1 + 180.0D0),360.0D0)
c      az2 = min(max(az2,0.0D0),360.0D0)
c      if (swflip) az2 = mod((az2 + 180.0D0),360.0D0)
c      el1 = min(max(el1,0.0D0), 90.0D0)
c      if (swflip) el1 = 180.0D0 - el1
c      el2 = min(max(el2,0.0D0), 90.0D0)
c      if (swflip) el2 = 180.0D0 - el2
c      write (c_command,'(''sleep '',F8.1)') max(abs(az2-az1)/azrate,abs(el2-el1)/elrate) + 0.2
      call date_and_time(c_date, c_time, c_zone, dattim)
      time_now = dble(dattim(5) * 3600.0 + dattim(6) * 60.0 + dattim(7) + dattim(8)/1000.0)
      time_sleep = pass_end - time_now
      write (c_command,'(''sleep '',F8.1)') time_sleep
      call system(c_command)
c-- Send to rest position
      call SPID_MD02(azpark, elpark, swdebug, .false.)
      call system('sleep 0.2')
c-- Remove offset if used
      if (swoffset) then
        do j = 1, n
          sat_az_track(j) = sat_az_track(j) - azoff
          sat_el_track(j) = sat_el_track(j) - eloff
        enddo
      endif
      return
      end

      subroutine SPID_MD02(Az, El, swdebug, swflip)
      implicit none
      real*8  Az, El
      logical swdebug, swflip
c
      integer*1 devbuf(13), statbuf(13), unkbuf(13), inbuf1(13), inbuf2(13), dumbuf(13)
      integer*4 lun, init, i, irec, ios, nstepsperdeg
      real*8 Azuse, Eluse, Azcorr, Elcorr, spd
      character*4 cangle
      data init/0/
      data statbuf/'57'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'1F'x,'20'x/
      data unkbuf /'57'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'3F'x,'20'x/
c
      real*8 azrange(2), elrange(2), azlow, azhig, ellow, elhig
      logical swazlim, swellim, swok, swtiltcor
c Reminder - nstepsperdeg can be 5 or 10 - if 10 is used you HAVE to disable the newline to carriage return-newline conversion for the COM port using stty -F /dev/ttyS4 -onlcr !
      data swazlim/.false./, swellim/.false./, swtiltcor/.false./, nstepsperdeg/10/
c
      real*8    azstor, elstor
      integer*4 dattim(8)
      character*5   c_zone
      character*8   c_date
      character*10  c_time
c
      save
c
      if (init.eq.0) then
        init = 1
        irec = 1
        call get_lun(lun)
        open (unit=lun,file='/dev/ttyS4',access='direct',recl=13)
      endif
      spd   = dble(nstepsperdeg)
      Azuse = min(max(Az,0.0D0),360.0D0)
      Eluse = min(max(El,0.0D0), 90.0D0)
      if (swflip) then
        Azuse = mod((Azuse + 180.0D0),360.0D0)
        Eluse = 180.0D0 - Eluse
      endif
c-- Check Az/El range if requested
      swok = .true.
      if (swazlim) then
        if (azuse.lt.azlow.or.azuse.gt.azhig) swok = .false.
      endif
      if (swellim) then
        if (eluse.lt.ellow.or.eluse.gt.elhig) swok = .false.
      endif
c-- Trying to find the pointing error
      azstor = Azuse
      elstor = Eluse
c-- Correct for system tilt - if requested - remember : If in flip mode I suspect this fails. Azuse is OK, but the sign of the tiltcor correction should be opposite ?
      if (swtiltcor) then
        Azcorr = Azuse
        Elcorr = Eluse
c        call tilt_correct(Azcorr, Elcorr)
c        call tilt_correct_cos(Azcorr, Elcorr)
        call tilt_correct_rot(Azcorr, Elcorr)
        Eluse = Elcorr
        Azuse = Azcorr
      endif
c-- Trying to find the pointing error
      if (swok) then
        call date_and_time(c_date, c_time, c_zone, dattim)
        if (swdebug) write (*,'(4I5,3(A,2F10.3),2I10)') (dattim(i), i=5,8), ' In : ', Az, El,' Used: ', azstor, elstor,' Tilt corrected : ', Azuse, Eluse, nint(spd * (Azuse + 360.0D0)), nint(spd * (Eluse + 360.0D0))
      endif
      if (.not.swok) goto 1
c      write (*,*) '==> OK'
c--
      write (cangle,'(I4.4)') nint(spd * (Azuse + 360.0D0))
      devbuf( 1) = '57'x
      devbuf( 2) = ichar(cangle(1:1))
      devbuf( 3) = ichar(cangle(2:2))
      devbuf( 4) = ichar(cangle(3:3))
      devbuf( 5) = ichar(cangle(4:4))
      if (nstepsperdeg.eq.5)  devbuf( 6) = '05'x
      if (nstepsperdeg.eq.10) devbuf( 6) = '0A'x
      write (cangle,'(I4.4)') nint(spd * (Eluse + 360.0D0))
      devbuf( 7) = ichar(cangle(1:1))
      devbuf( 8) = ichar(cangle(2:2))
      devbuf( 9) = ichar(cangle(3:3))
      devbuf(10) = ichar(cangle(4:4))
      if (nstepsperdeg.eq.5)  devbuf(11) = '05'x
      if (nstepsperdeg.eq.10) devbuf(11) = '0A'x
      devbuf(12) = '2F'x
      devbuf(13) = '20'x
      if (swdebug) write (*,'(4(Z2.2,1x),2x,4(Z2.2,1x))') (devbuf(i), i = 2,5), (devbuf(i), i = 7,10)
      write (lun,rec=irec) devbuf
      irec = irec + 1
      write (lun,rec=irec) statbuf
      irec = irec + 1
      write (lun,rec=irec) unkbuf
      irec = irec + 1
 1    return
c
      entry SPID_MD02_AzLim(azrange)
      azlow = azrange(1)
      azhig = azrange(2)
      swazlim = .true.
      return
c
      entry SPID_MD02_ElLim(elrange)
      ellow = elrange(1)
      elhig = elrange(2)
      swellim = .true.
      return
c
      entry SPID_MD02_tiltcor()
      swtiltcor = .true.
      return
c
      end
      
      subroutine tilt_correct(Azcorr, Elcorr)
      implicit none
      real*8 Azcorr, Elcorr
c
      integer*4 i
      real*4 azarr(13), elarr(13), frac, fcor
      data azarr/ 0.0,   30.0,  60.0,  90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 360.0/
      data elarr/0.66, -0.020, -0.58, -0.96, -1.03, -0.76, -0.20,  0.44,  1.01,  1.41,  1.44,  1.17,  0.65/

c
      do i = 1,12
        if (azarr(i).le.Azcorr.and.Azcorr.le.azarr(i+1)) then
          frac = (sngl(Azcorr) - azarr(i)) / (azarr(i+1) - azarr(i))
          fcor = elarr(i) + frac * (elarr(i+1) - elarr(i))
          Elcorr = Elcorr - dble(fcor)
        endif
      enddo
c
      return
      end

      subroutine tilt_correct_cos(Azcorr, Elcorr)
      implicit none
      real*8 Azcorr, Elcorr
c
      real*8 torad
c-- This is an anlytical model of the El variation over 360 deg Az (measured in 10 deg steps) of my set up. Better then 0.1 deg accurate.
      torad           = 2.0D0 * dasin(1.0D0) / 180.0D0
      Elcorr = Elcorr - (cos((azcorr - 295.0D0) * torad) * 1.17D0 + 0.175D0)
      return
      end

      subroutine tilt_correct_rot(Azcorr, Elcorr)
      implicit none
      real*8 Azcorr, Elcorr
c
      real*4 Azcorr_r4, Elcorr_r4
      real*4 nullvec(3), point_calc(3), point_rot(3), nullaz, nullel, torad, todeg, azrot, azoff, eloff, azcal, elcal
c
c
      Azcorr_r4  = sngl(Azcorr)
      Elcorr_r4  = sngl(Elcorr)
      torad      = 2.0 * asin(1.0) / 180.0
      todeg      = 180.0 / (2.0 * asin(1.0))
c-- This is the Azimuth in between the two extremes of the measured inclinometer curve for El=0
      nullaz     = 203.0
      nullel     = 0.0
      nullvec(1) = cos(nullel * torad) * cos (nullaz * torad)
      nullvec(2) = cos(nullel * torad) * sin (nullaz * torad)
      nullvec(3) = sin(nullel * torad)
c-- Assume SYRACUSE 3B is used for calibration - fix later
      point_calc(1) = cos(29.8 * torad) * cos (192.3 * torad)
      point_calc(2) = cos(29.8 * torad) * sin (192.3 * torad)
      point_calc(3) = sin(29.8 * torad)
      call rotate_vector_a_around_b(point_calc, nullvec, 1.25, point_rot)
      azrot = atan2(point_rot(2), point_rot(1)) * todeg
      if (azrot.lt.0) azrot = azrot + 360.0
      elcal = asin(point_rot(3)) * todeg - 29.8
      azcal = azrot - 192.3
c
      point_calc(1) = cos(Elcorr_r4 * torad) * cos (Azcorr_r4 * torad)
      point_calc(2) = cos(Elcorr_r4 * torad) * sin (Azcorr_r4 * torad)
      point_calc(3) = sin(Elcorr_r4 * torad)
      call rotate_vector_a_around_b(point_calc, nullvec, 1.25, point_rot)
      azrot = atan2(point_rot(2), point_rot(1)) * todeg
      if (azrot.lt.0) azrot = azrot + 360.0
c      write (*,*) - ((azrot - Azcorr_r4) - azcal), - ((asin(point_rot(3)) * todeg - Elcorr_r4) - elcal)
      Azcorr_r4 = Azcorr_r4 - ((azrot - Azcorr_r4) - azcal)
      Elcorr_r4 = Elcorr_r4 - ((asin(point_rot(3)) * todeg - Elcorr_r4) - elcal)
c
      Azcorr = dble(Azcorr_r4)
      Elcorr = dble(Elcorr_r4)
c
      return
      end

      subroutine rotate_vector_a_around_b(a, b, theta, c)
      implicit none
      real*4 a(3), b(3), c(3), theta, const, x1, x2
c
      real*4 a_par_b(3), a_ortho_b(3), a_ortho_b_theta(3), w(3)
      real*4 torad, angle, rnorm_ortho
c
      torad = 2.0 * asin(1.0) / 180.0
      angle = theta * torad
c
      const = (a(1) * b(1) + a(2) * b(2) + a(3) * b(3)) / (b(1) * b(1) + b(2) * b(2) + b(3) * b(3))
      a_par_b(1)   = const * b(1)
      a_par_b(2)   = const * b(2)
      a_par_b(3)   = const * b(3)
      a_ortho_b(1) = a(1) - a_par_b(1)
      a_ortho_b(2) = a(2) - a_par_b(2)
      a_ortho_b(3) = a(3) - a_par_b(3)
      w(1)         = b(2) * a_ortho_b(3) - b(3) * a_ortho_b(2)
      w(2)         = b(3) * a_ortho_b(1) - b(1) * a_ortho_b(3)
      w(3)         = b(1) * a_ortho_b(2) - b(2) * a_ortho_b(1)
      rnorm_ortho  = sqrt(a_ortho_b(1) * a_ortho_b(1) + a_ortho_b(2) * a_ortho_b(2) + a_ortho_b(3) * a_ortho_b(3))
      x1           = cos(angle) / rnorm_ortho
      x2           = sin(angle) / sqrt(w(1) * w(1) + w(2) * w(2) + w(3) * w(3))
      a_ortho_b_theta(1) = rnorm_ortho * (x1 * a_ortho_b(1) + x2 * w(1))
      a_ortho_b_theta(2) = rnorm_ortho * (x1 * a_ortho_b(2) + x2 * w(2))
      a_ortho_b_theta(3) = rnorm_ortho * (x1 * a_ortho_b(3) + x2 * w(3))
c
      c(1)         = a_ortho_b_theta(1) + a_par_b(1)
      c(2)         = a_ortho_b_theta(2) + a_par_b(2)
      c(3)         = a_ortho_b_theta(3) + a_par_b(3)
c
      return
      end

      subroutine sunazel_init()
      implicit none
      integer*4 linedoy, dattim_out(8)
      real*8    solaz, solel, solaz1, solel1, solaz2, solel2, anglesun, delta_t
c
      integer*2 npos
      integer*4 dattim(8), monthdays(12), diy, doy, i, timezone
      real*4 days_in_year, pi, gamma, eqtime, decl, time_offset, tst, ha, long, lat, ha_rise, ha_set, sunrise, sunset, solzenang, rsc, tstval
      real*8 torad, todeg, mylong, mylat, x1, y1, z1, x2, y2, z2
      character*5  c_zone
      character*8  c_date
      character*10 c_time
      character*50 instring
c
      data monthdays/31,0,31,30,31,30,31,31,30,31,30,31/
c
      save
c
      pi              = 2.0 * asin(1.0)
      torad           = 2.0D0 * dasin(1.0D0) / 180.0D0
      todeg           = 180.0D0 / (2.0D0 * dasin(1.0D0))
      call getenv('LONG',instring)
      call string_to_r8(instring(1:lnblnk(instring)), npos, mylong)
      call getenv('LAT',instring)
      call string_to_r8(instring(1:lnblnk(instring)), npos, mylat)

c
      call date_and_time(c_date, c_time, c_zone, dattim)
      timezone = nint(float(dattim(4))/60.0)
c
      diy = 365
      if (mod(dattim(1),4).eq.0.and.dattim(1).ne.2000) diy = diy + 1
      days_in_year    = float(diy)
      return
c
      entry sunazel(dattim_out, solaz, solel, delta_t)
      long       = mylong * torad
      lat        = mylat  * torad
      call date_and_time(c_date, c_time, c_zone, dattim)
      do i = 1,8
        dattim_out(i) = dattim(i)
      enddo
      call datetodoy(dattim(1), dattim(2), dattim(3), doy)
      rsc        = delta_t + dble(dattim(7)) + dble(dattim(8)) / 1000.0D0
c
      gamma      = ((2.0 * pi) / days_in_year) * (float(doy) - 1.0 + ((float(dattim(5))-12.0)/24.0) + float(dattim(6))/1440.)
c
      eqtime     = 229.18*(0.000075 + 0.001868*cos(gamma) - 0.032077*sin(gamma) - 0.014615*cos(2.0*gamma) - 0.040849*sin(2.0*gamma))
      decl       = 0.006918 - 0.399912*cos(gamma) + 0.070257*sin(gamma) - 0.006758*cos(2.0*gamma) + 0.000907*sin(2.0*gamma) - 0.002697*cos(3.0*gamma) + 0.00148*sin(3.0*gamma)
      tstval     = float(dattim(5)) * 60.0 + float(dattim(6)) + (rsc / 60.0)
c
      time_offset= eqtime + 4.0 * mylong - 60.0 * float(timezone)
      tst        = tstval + time_offset
      ha         = ((tst / 4.0) - 180.0) * torad
      solzenang  = acos(sin(lat) * sin(decl) + cos(lat) * cos(decl) * cos(ha))
      if (ha.ge.0.0D0) then
        solaz      = dble(360.0 - acos(-((sin(lat) * cos(solzenang) - sin(decl))/(cos(lat) * sin(solzenang)))) * todeg)
      else
        solaz      = dble(acos(-((sin(lat) * cos(solzenang) - sin(decl))/(cos(lat) * sin(solzenang)))) * todeg)
      endif
      solzenang  = solzenang * todeg
      solel      = dble(90.0 - solzenang)
c
      return
c
      entry sunangle(solaz1, solel1, solaz2, solel2, anglesun)
      x1 = cos(solaz1 * torad) * cos(solel1 * torad)
      y1 = sin(solaz1 * torad) * cos(solel1 * torad)
      z1 =                       sin(solel1 * torad)
      x2 = cos(solaz2 * torad) * cos(solel2 * torad)
      y2 = sin(solaz2 * torad) * cos(solel2 * torad)
      z2 =                       sin(solel2 * torad)
      anglesun = todeg * acos((x1*x2 + y1*y2 + z1*z2) / (sqrt(x1*x1 + y1*y1 + z1*z1) * sqrt(x2*x2 + y2*y2 + z2*z2)))
      end
      
      subroutine calibrate(azpark_in, elpark_in, azoff_in, eloff_in)
      implicit none
      real*8 azpark_in, elpark_in, azoff_in, eloff_in
c
      real*8 azstep, elstep, azpark, elpark, azoff, eloff
      character*1 cin
c
      azstep = 0.25D0
      elstep = 0.25D0
      azpark = azpark_in
      elpark = elpark_in
      azoff  = azoff_in
      eloff  = eloff_in
      write (*,'(''Commands possible u,U=up, d,D=down, r,R=right, l,L=left, z,Z=reset, q,Q=quit and return'')')
      write (*,'(''Commands in capitals use 0.5 deg steps, otherwise 0.25 deg'')')
      do while (.true.)
        read(*,'(a1)') cin
        if (cin(1:1).eq.'q'.or.cin(1:1).eq.'Q') goto 1
        if (cin(1:1).eq.'z'.or.cin(1:1).eq.'Z') then
          azoff  = azoff_in
          eloff  = eloff_in
        endif
        if (cin(1:1).eq.'u') then
          eloff = eloff + elstep
        endif
        if (cin(1:1).eq.'d') then
          eloff = eloff - elstep
        endif
        if (cin(1:1).eq.'l') then
          azoff = azoff - azstep
        endif
        if (cin(1:1).eq.'r') then
          azoff = azoff + azstep
        endif
        if (cin(1:1).eq.'U') then
          eloff = eloff + 2.0D0 * elstep
        endif
        if (cin(1:1).eq.'D') then
          eloff = eloff - 2.0D0 * elstep
        endif
        if (cin(1:1).eq.'L') then
          azoff = azoff - 2.0D0 * azstep
        endif
        if (cin(1:1).eq.'R') then
          azoff = azoff + 2.0D0 * azstep
        endif
c-- As there is a read in between commands I ignore here the time it takes to move and don't do a sleep
        call SPID_MD02(azpark + azoff, elpark + eloff, .false., .false.)
      enddo
 1    azoff_in  = azoff
      eloff_in  = eloff
      write (*,'(''Re-calibrated Azimuth   Offset           : '',F10.3)') azoff_in
      write (*,'(''Re-calibrated Elevation Offset           : '',F10.3)') eloff_in
      return
      end

      logical function check_horizon(az, el)
      implicit none
      real*8 az, el
c
      real*4 elmask(38),azvals(38), elcheck
      integer*4 i
c
c-- Lazy approach at the start and end point the 355 deg (start) and 5 deg (end) values are repeated, to aid interpolating.
c-- Original data
c      data elmask/40.9  ,39.2  ,30.2  ,19.5  , 9.7  ,16.4  ,17.8  ,20.4  ,27.6  ,30.2  ,23.1  ,
c     *                   19.5,   4.4,   3.0,   5.3,   7.0,   5.3,   8.8,   8.8,   9.7,  10.6  ,
c     *                   11.5,  10.6,   9.7,  11.1,  11.1,  18.2,  19.1,  10.6,  13.3,  15.1  ,
c     *                   16.4,  18.6,  30.2,  37.4,  39.2,  40.9,                               39.2/
c      data azvals/  -5.0,   5.0,  15.0,  25.0,  35.0,  45.0,  55.0,  65.0,  75.0,  85.0,  95.0,
c     *                    105.0, 115.0, 125.0, 135.0, 145.0, 155.0, 165.0, 175.0, 185.0, 195.0,
c     *                    205.0, 215.0, 225.0, 235.0, 245.0, 255.0, 265.0, 275.0, 285.0, 295.0,
c     *                    305.0, 315.0, 325.0, 335.0, 345.0, 355.0,                            365.0/
c
c-- Modded data (based on observations in X-band)
      data elmask/40.9  ,39.2  ,30.7  ,19.5  , 9.7  ,16.4  ,17.8  ,20.4  ,27.6  ,30.2  ,23.1  ,
     *                   19.5,   4.4,   3.0,   5.3,   6.5,   5.3,   8.8,   8.8,   9.7,  10.6  ,
     *                   11.0,  10.6,   9.7,  11.1,  11.1,  18.2,  19.1,  10.6,  13.3,  15.1  ,
     *                   16.4,  18.1,  30.2,  37.4,  39.2,  40.9,                               39.2/
      data azvals/  -5.0,   5.0,  15.0,  25.0,  35.0,  45.0,  55.0,  65.0,  75.0,  85.0,  95.0,
     *                    105.0, 115.0, 125.0, 135.0, 145.0, 155.0, 165.0, 175.0, 185.0, 195.0,
     *                    205.0, 215.0, 225.0, 235.0, 245.0, 255.0, 265.0, 275.0, 285.0, 295.0,
     *                    305.0, 315.0, 325.0, 335.0, 345.0, 355.0,                            365.0/
c
      check_horizon = .false.
      if (az.le.azvals(1).or.az.gt.azvals(38)) return
      do i = 1, 37
        if (azvals(i).lt.az.and.az.le.azvals(i+1)) goto 1
      enddo
c-- linear interpolation
 1    elcheck = ((az - azvals(i)) / (azvals(i+1) - azvals(i))) * (elmask(i+1) - elmask(i)) + elmask(i)
      if (el.gt.dble(elcheck)) check_horizon = .true.
c
      return
      end

      subroutine tletrack_help()
c
      write (*,'(''-- Run the program as follows: '')')
      write (*,'('' '')')
      write (*,'(''./tletrack.exe \"AQUA\" resource beam=1.8 etc - see below for switches '')')
      write (*,'('' '')')
      write (*,'(''   First argument satellite name as in the TLE between double quotes (in UNIX command line speak) or an indirect to a file using @ '')')
      write (*,'(''   Examples: '')')
      write (*,'(''   ./tletrack.exe \"NOAA\ 19\" weather beam=1.8 '')')
      write (*,'(''   ./tletrack.exe @./examples/tlesats.txt beam=1.8 '')')
      write (*,'('' '')')
      write (*,'('' Where tlesats.txt contained: '')')
      write (*,'('' '')')
      write (*,'('' "FENGYUN 3B"   weather   9'')')
      write (*,'('' "FENGYUN 3C"   weather   9'')')
      write (*,'('' "METOP-A"      weather   5'')')
      write (*,'('' "METOP-B"      weather   8'')')
      write (*,'('' "METOP-C"      weather   8'')')
      write (*,'('' "NOAA 15"      weather   4'')')
      write (*,'('' "NOAA 18"      weather   7'')')
      write (*,'('' "NOAA 19"      weather   7'')')
      write (*,'('' '')')
      write (*,'('' with the first " left aligned in column 1 ! '')')
      write (*,'('' In this file, like on the command line, the second argument is the type of TLE - resource or weather'')')
      write (*,'('' In this file the 3rd argument is the priority for a certain satellite. Used in de-conflicting'')')
      write (*,'('' '')')
      write (*,'('' A special use case is :'')')
      write (*,'(''   ./tletrack.exe @sun beam=1.0'')')
      write (*,'('' which enables look ahead tracking of the sun - JPL method below - using @horcmd=10 - is much better !'')')
      write (*,'('' '')')
      write (*,'('' A second special use case is :'')')
      write (*,'(''   ./tletrack.exe @horcmd=499 beam=1.0 tiltcor'')')
      write (*,'('' which enables look ahead tracking of the JPL HORIZONS objects (here Mars ; which has ID # 499)'')')
      write (*,'('' NOTE 1 : Also see noload command below'')')
      write (*,'('' NOTE 2 : At the moment it is a primitive implementation that points to the object position every minute with a 30 second lookahead'')')
      write (*,'('' NOTE 3 : Time accuracy  is below 1 second'')')
      write (*,'('' '')')
      write (*,'(''   Second argument is the type  of TLE : resource or weather IF the first argument is not an indirect ! '')')
      write (*,'(''                   NOT required if first argument is reference to indirect list (using @) '')')
      write (*,'('' '')')
      write (*,'(''   Further arguments:   '')')
      write (*,'('' '')')
      write (*,'(''     beam=1.8        beamwidth of dish to be used (in deg)'')')
      write (*,'(''     debug           show the commands that would be issued for reactive and pro-active (look ahead) '')')
      write (*,'(''     live            actually command the hardware '')')
      write (*,'(''     flip            allow dish flip mode if crossing the N vector '')')
      write (*,'(''     noovl           removes Az, El commands overlapping in time - hope this cures the MD-02 M2 pulse timeout problem ....... '')')
      write (*,'(''     both            do the calculations with and without look-ahead - only useful with debug or diag '')')
      write (*,'('' '')')
      write (*,'(''     raster=200.5,28.3,3.0,2.0,0.4,0.3,10.0 parameters are Az, El, Az range, El range, Az step, El step, delta-T (sec) - perform a raster scan'')')
      write (*,'('' '')')
      write (*,'('' '')')
      write (*,'(''     cross=200.5,28.3,3.0,2.0,0.4,0.3 parameters are Az, El, Az range, El range, Az step, El step'')')
      write (*,'(''                     This is to perform an Az/El antenna profile centroiding and asks for a signal strength after every pointing'')')
      write (*,'(''                     then calculates the maximum Az/El response point and moves the antenna there'')')
      write (*,'('' '')')
      write (*,'(''     NOTE : The program terminates after a raster or a cross'')')
      write (*,'('' '')')
      write (*,'(''     azrange=130-290  only allow commands with an Azimuth   in this range - used when running tests at the dish'')')
      write (*,'(''     elrange=0-45     only allow commands with an Elevation in this range - used when running tests at the dish'')')
      write (*,'(''     eltrackmin=8.0   limits the minimum elevation where tracking starts to a user defined value (here 8.0 deg - default is 3.0 deg)'')')
      write (*,'(''     tiltcor          enable the tilt correction algorithm - rotation of the dish-object vector around a user defined vector by a user defined amount'')')
      write (*,'(''                      this allows to compensate for not-perfect verticality of the center pole'')')
      write (*,'(''     azpark=180.0     Override the default park position of Az=209.0'')')
      write (*,'(''     elpark=45.0      Override the default park position of El=0.0'')')
      write (*,'(''     azoff=-1.0       Define an Az offset - found from pre-calibration - to add to the calculated values'')')
      write (*,'(''     eloff=0.9        Define an El offset - found from pre-calibration - to add to the calculated values'')')
      write (*,'(''     calib            inserts an Az/El calibration routine - assumes you defined azpark, elpark, azoff and eloff !!'')')
      write (*,'(''     azelalt          Enables my own Az/El calculation - which is slightly different - but likely incorrect.'')')
      write (*,'(''     deltat=0.9       Defines an offset time which is added (!) to all calculatedcommanding times'')')
      write (*,'(''     predict=8        Limits the prediction window to the user defined value (here 8) instead of 24. Ignored for using JPL HORIZONS'')')
      write (*,'(''     noload           Do not reload the HORIZONS file - re-uses the existing tempeph.txt - is faster and avoids unnecessarily loading the JPL Horizons system'')')
      write (*,'(''     point=100.0,50.0 Points the dish to the defined Az, El (if combined with the live option) - ALL OTHER OPTIONS ARE IGNORED !!!'')')
      write (*,'(''                      A valid command line could be :'')')
      write (*,'(''                      ./eltrack.exe \"AQUA\" resource point=200.0,20.0 live'')')
      write (*,'(''                      This would point the dish to that position and then stop'')')
      write (*,'(''     doublec          In certain conditions (point, @horcmd) execute the pointing command twice. This allows the MD-02 to converge better'')')
      write (*,'(''     jumpin           When asking for a single satellite to be tracked, this tries to jump in live if the sat is above the horizon'')')
      write (*,'(''     hmask            Use the user defined (in the code) horizon mask to determine the start and end of a pass - my location derived from L-band'')')
      write (*,'(''                      This switch is also useful if running like : ./tletrack.exe @./examples/tlesats.txt beam=0.6 predict=9 hmask'')')
      write (*,'(''                      This generates the visible passes in the next 9 hours - with azel plots in TLE directory: green - visible, red - invisible at start of pass'')')
      write (*,'(''                      yellow - invisible at end of pass, purple - when green, visibility for co-visibility station requested'')')
      write (*,'(''     covis=xxx        Analyse the co-visibility of other groundstations - used to catch satellites when dumping to commercial groundstations.'')')
      write (*,'(''                      Currently xxx=KIR (Kiruna) and SVB (Svallbard) are supported'')')
      write (*,'('' '')')
      write (*,'('' The TLE environment variable defines where to find the TLE files and where the azel plots will be deposited '')')
      write (*,'('' The LONG and LAT environment variables define the geographic location of the dish (also used for JPL HORIZONS)'')')
c
      return
      end

      subroutine estimate_peak(x, resp, n, p1, p2, p3)
      integer*4 n
      real*4 x(n), resp(n), p1, p2, p3
c
      real*4 sum, sumy, sumxy, yval, x_a, y_a
c
      p2 = (resp(n) - resp(1)) / (x(n) - x(1))
      p1 = resp(1) - p2 * x(1)
      sum   = 0.
      sumy  = 0.
      sumxy = 0.
      do i = 1, n
        yval = resp(i) - (p1 + p2 * x(i))
        sumy  = sumy  + yval
        sumxy = sumxy + x(i) * yval
      enddo
      p3 = 0.
      if (sumy.ne.0.) p3 = sumxy / sumy
      do i = 1, n
        yval = resp(i) - (p1 + p2 * x(i))
        sum  = sum + yval * ((x(i)-p3)**2.)
      enddo
      p2 = 0.
      if (sumy.gt.1.) p2 = sqrt(sum/(sumy-1.))
      p1 = sumy
      return
      end

