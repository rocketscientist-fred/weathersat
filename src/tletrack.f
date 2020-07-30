      implicit none
c--
c-- nmaxpass is the maximum nr of passes that can be stored over the to be analysed period
c-- nmaxpointing is the maximum nr of timesteps that can be stored per single pass
c-- nsatmax is the maximum nr of satellites that can be analysed in one single run
c
      integer*4 nmaxpass, nmaxpointing, nsatmax
      parameter (nmaxpass=1000, nmaxpointing=40000, nsatmax=50)
c
      real*8    sat_az_pass(nmaxpointing,nmaxpass)    , sat_el_pass(nmaxpointing,nmaxpass)    , time_pass(nmaxpointing,nmaxpass)
      real*8    sat_az_track(nmaxpointing,nmaxpass)   , sat_el_track(nmaxpointing,nmaxpass)   , time_track(nmaxpointing,nmaxpass)
      real*8    sat_az_track_la(nmaxpointing,nmaxpass), sat_el_track_la(nmaxpointing,nmaxpass), time_track_la(nmaxpointing,nmaxpass)
      real*8    sat_pass_start(nmaxpass), sat_pass_end(nmaxpass), sat_el_max_pass(nmaxpass)
      integer*4 npointspass(nmaxpass), n_track(nmaxpass), n_track_la(nmaxpass), i_sat_prio(nsatmax), i_start(nsatmax), i_end(nsatmax)
      integer*4 nrpasses, key(nmaxpass), ok(nmaxpass), sat_pass_prio(nmaxpass)
      integer*1 ovlflag(nmaxpointing)
      character*20 c_sat_pass(nmaxpass), c_sat_use(nsatmax), c_sat_type(nsatmax)
c-- Diagnostic variables - can be deleted if plots no longer needed
      integer*4 i1, i2, i3, i4, i5, i6
      real*8    angle1, angle2, azpark, elpark
c
      real*8    doy_tle, time, ro(3), vo(3), long, lat, theta0g, time_start, sat_az, sat_el, rsc, sat_el_max, time_max, el_track_min
      real*8    torad, todeg, track_angle, x1, y1, z1, x2, y2, z2, time_save, secsubsteps, time_offset, rsc1, rsc2, elpriothresh
      real*8    azrate, elrate, azoff, eloff, azrange(2), elrange(2), sunaz, sunel, sunaz1, sunel1, sunaz2, sunel2, anglesun, sunrate
      real*8    suntim1, suntim2, sunrepoint, sunahead
c-- Raster variables
      real*8  grid_az, grid_el, grid_az_range, grid_el_range, grid_az_step, grid_el_step, grid_delta_t, az_cmd, el_cmd, t_move, az1, el1
      logical swgrid
c--
      integer*4 dattim(8), luntle, iyr_tle, doy, i, j, k, l, ihr, imn, i_year, i_month, i_day, i_pass, i_track, lunsat, nsatuse, ios, isat
      integer*4 i_arg1, i_arg, ihr1, ihr2, imn1, imn2
      integer*2 npos
      logical foundtle, swdebug, swdiag, swlive, swallowflip, swnoovl, swboth, swoffset, swazrange, swelrange, swtiltcor, swsun, swcalib
      character*5   c_zone
      character*8   c_date, c_datvis
      character*10  c_time
      character*20  c_sat
      character*50  c_file, c_command
      character*69  tlestring
      character*250 argstring, tledir, command
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
c
      torad           = 2.0D0 * dasin(1.0D0) / 180.0D0
      todeg           = 180.0D0 / (2.0D0 * dasin(1.0D0))
      el_track_min    = 3.0D0
      track_angle     = 1.8D0
c-- Use 0.05 sec time steps (required for X-band)
      secsubsteps     = 20.0D0
      swdebug         = .false.
      swdiag          = .false.
      swlive          = .false.
      swallowflip     = .false.
      swnoovl         = .false.
      azrate          = 4.5D0
      elrate          = 5.2D0
      swboth          = .false.
      swoffset        = .false.
      azoff           = 0.0D0
      eloff           = 0.0D0
      swgrid          = .false.
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
c-- If two satellites passes have a conflict and their priorities differ by 1 and their max elevations are above this limit : the max elevation wins
      elpriothresh    = 30.0D0
c
      do i = 1, nmaxpass
        sat_pass_start(i) = 0.0D0
        sat_pass_end(i)   = 0.0D0
        ok(i)             = 1
      enddo
      call getarg(1, argstring)
      if (index(argstring(1:lnblnk(argstring)),'-help').ne.0) then
        call tletrack_help()
        stop '** Help end **'
      endif
      if (argstring(1:1).eq.'@') then
        if (argstring(2:4).eq.'sun') then
          swsun = .true.
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
      do k = i_arg1, 20
        call getarg(k, argstring)
        if (index(argstring(1:lnblnk(argstring)),'debug').ne.0)  swdebug      = .true.
        if (index(argstring(1:lnblnk(argstring)),'diag').ne.0)   swdiag       = .true.
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
      enddo
c-- Check if dish alive and aligned - using visual inspection from my study :-)
      if (swlive) then
        call SPID_MD02(250.1D0, 0.0D0, swdebug, .false.)
        write (*,'(''If dish moved and aligned press enter :'')')
        read (*,*)
        call SPID_MD02(azpark, elpark, swdebug, .false.)
c-- Wait for dish to return to the park position
        t_move = abs(250.1D0 - azpark)/azrate
        write (c_command,'(''sleep '',F8.1)') t_move
        call system(c_command)
        if (swcalib) call calibrate(azpark, elpark, azoff, eloff)
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
c        call map_azel_init()
        call sgp4_init(ro, vo, long, lat)
        time        = time_start
        time_offset = 0.0D0
        i_pass   = 0
        do i = 1, 86400 * nint(secsubsteps)
          time = time + (1.0D0 / secsubsteps)
          if (time.gt.86400.0D0) then
            doy = doy + 1
            time        = time        - 86400.0D0
            time_offset = time_offset + 86400.0D0
          endif
          call run_TLE(time, doy, ro, vo, long, lat, theta0g, doy_tle, iyr_tle, dattim(1))
          call azel(long, lat, ro, sat_az, sat_el)
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
            call doytodate(doy, dattim(1), c_datvis)
c            call map_azel_place(sat_az, sat_el)
            ihr = int(time/3600.0D0)
            imn = (int(time - dble(ihr)*3600.0D0)/60.0D0)
            rsc = time - dble(ihr)*3600.0D0 - imn * 60.0D0
            sat_az_pass(npointspass(nrpasses), nrpasses) = sat_az
            sat_el_pass(npointspass(nrpasses), nrpasses) = sat_el
            time_pass  (npointspass(nrpasses), nrpasses) = time + dattim(4) * 60.0D0
c          write (*,'(2F10.3,1x,a8,1x,2I5,F10.3, I5)') sat_az, sat_el, c_datvis, ihr, imn, rsc, doy
          else
            if (i_pass.eq.1) then
              if (sat_pass_end(nrpasses).eq.0.0D0) sat_pass_end(nrpasses) = time  + time_offset + dattim(4) * 60.0D0
            endif
            i_pass = 0
          endif
 6        continue
        enddo
        i_end(isat) = nrpasses
        call sgp4_end()
c        call map_azel_write(command)
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
c-- Reactive tracking - the usual way
          if (.not.swboth) goto 10
          i_track = 0
          do j = 1, npointspass(i)
            if (sat_el_pass(j,i).gt.el_track_min.and.i_track.eq.0) then
              i_track                    = 1
              n_track(i)                 = 1
              sat_az_track(n_track(i),i) = sat_az_pass(j,i)
              sat_el_track(n_track(i),i) = sat_el_pass(j,i)
              time_track(n_track(i),i)   = time_pass(j,i)
            endif
            if (sat_el_pass(j,i).lt.el_track_min.and.i_track.eq.1) then
              i_track                    = 0
              n_track(i)                 = n_track(i) + 1
              sat_az_track(n_track(i),i) = sat_az_pass(j,i)
              sat_el_track(n_track(i),i) = sat_el_pass(j,i)
              time_track(n_track(i),i)   = time_pass(j,i)
            endif
            if (i_track.eq.1) then
              x1 = cos(torad*sat_az_track(n_track(i),i)) * cos(torad*sat_el_track(n_track(i),i))
              y1 = sin(torad*sat_az_track(n_track(i),i)) * cos(torad*sat_el_track(n_track(i),i))
              z1 =                                         sin(torad*sat_el_track(n_track(i),i))
              x2 = cos(torad*sat_az_pass(j,i)) * cos(torad*sat_el_pass(j,i))
              y2 = sin(torad*sat_az_pass(j,i)) * cos(torad*sat_el_pass(j,i))
              z2 =                               sin(torad*sat_el_pass(j,i))
              if (todeg * acos(min(max((x1*x2 + y1*y2 + z1*z2),-1.0D0),1.0D0)).ge.track_angle) then
                n_track(i)                 = n_track(i) + 1
                sat_az_track(n_track(i),i) = sat_az_pass(j,i)
                sat_el_track(n_track(i),i) = sat_el_pass(j,i)
                time_track(n_track(i),i)   = time_pass(j,i)
              endif
            endif
          enddo
          if (swdebug.and.n_track(i).gt.0) write (*,'(''Reactive tracking'')')
          if (swdebug.and.n_track(i).gt.0) write (*,'(''Pointing #  Az          El     Time'')')
          do j = 1, n_track(i)
            ihr = int(time_track(j,i)/3600.0D0)
            imn = (int(time_track(j,i) - dble(ihr)*3600.0D0)/60.0D0)
            rsc = time_track(j,i) - dble(ihr)*3600.0D0 - imn * 60.0D0
            if (swdebug) write (*,'(I5,3x,2F10.3,2I5,F10.3)') j, sat_az_track(j,i), sat_el_track(j,i), ihr, imn, rsc
          enddo
 10       continue
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
            endif
            if (sat_el_pass(j,i).lt.el_track_min.and.i_track.eq.1) then
              i_track                          = 0
              n_track_la(i)                    = n_track_la(i) + 1
              sat_az_track_la(n_track_la(i),i) = sat_az_pass(j,i)
              sat_el_track_la(n_track_la(i),i) = sat_el_pass(j,i)
              time_track_la(n_track_la(i),i)   = time_pass(j,i)
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
          if (sat_el_max.gt.el_track_min) then
            write (*,'(2i5,F10.3,2I5,F10.3,5x,A,I5)') i, nint(dble(npointspass(i))/secsubsteps), sat_el_max, ihr, imn, rsc, c_sat_use(isat)(2:lnblnk(c_sat_use(isat))-1), n_track_la(i)
          endif          
        enddo
      enddo
      call free_lun(luntle)
c
c-- Diagnostics section - plot difference angles dish <-> sat for reactive and pro-active (look ahead) pointing
c
      if (.not.swdiag) goto 3
      do i = 1, nrpasses
        write (*,'(''Diagnostic data for pass '',I3)') i
        write (*,*)
        write (*,'(''                                                  Delta-   Delta-'')')
        write (*,'(''Serial #  Az         El     Time                  angle    angle LA'')')
        if (n_track(i).gt.0) then
          do j = 1, npointspass(i)
            i1 = 1
            i2 = n_track(i)
            do while (i2-i1.gt.1)
              i3 = (i2 + i1) / 2
              if (time_pass(j,i).gt.time_track(i3,i)) then
                i1 = i3
              else
                i2 = i3
              endif
            enddo
c            write (*,*) i, j, i1, i2, i3, time_pass(j,i), time_track(i1,i), time_track(i2,i)
c            read (*,*)
            x1 = cos(torad*sat_az_track(i1,i)) * cos(torad*sat_el_track(i1,i))
            y1 = sin(torad*sat_az_track(i1,i)) * cos(torad*sat_el_track(i1,i))
            z1 =                                 sin(torad*sat_el_track(i1,i))
            x2 = cos(torad*sat_az_pass(j,i)) * cos(torad*sat_el_pass(j,i))
            y2 = sin(torad*sat_az_pass(j,i)) * cos(torad*sat_el_pass(j,i))
            z2 =                               sin(torad*sat_el_pass(j,i))
            angle1 = todeg * acos(min(max((x1*x2 + y1*y2 + z1*z2),-1.0D0),1.0D0))
            i4 = 1
            i5 = n_track_la(i)
            do while (i5-i4.gt.1)
              i6 = (i5 + i4) / 2
              if (time_pass(j,i).gt.time_track_la(i6,i)) then
                i4 = i6
              else
                i5 = i6
              endif
            enddo
            x1 = cos(torad*sat_az_track_la(i4,i)) * cos(torad*sat_el_track_la(i4,i))
            y1 = sin(torad*sat_az_track_la(i4,i)) * cos(torad*sat_el_track_la(i4,i))
            z1 =                                    sin(torad*sat_el_track_la(i4,i))
            x2 = cos(torad*sat_az_pass(j,i)) * cos(torad*sat_el_pass(j,i))
            y2 = sin(torad*sat_az_pass(j,i)) * cos(torad*sat_el_pass(j,i))
            z2 =                               sin(torad*sat_el_pass(j,i))
            angle2 = todeg * acos(min(max((x1*x2 + y1*y2 + z1*z2),-1.0D0),1.0D0))
            ihr = int(time_pass(j,i)/3600.0D0)
            imn = (int(time_pass(j,i) - dble(ihr)*3600.0D0)/60.0D0)
            rsc = time_pass(j,i) - dble(ihr)*3600.0D0 - imn * 60.0D0
            write (*,'(I5,2F10.3,2I5,3F10.3)') j, sat_az_pass(j,i), sat_el_pass(j,i), ihr, imn, rsc, angle1, angle2
          enddo
        endif
      enddo
 3    continue
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
            if (swdebug) write (*,'(''Conflict between passes : '',4I5)')  i, j, key(i), key(j)
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
c      do i = 1, nrpasses
c        if (sat_el_max_pass(key(i)).gt.el_track_min) write (*,*) key(i), sat_pass_start(i), sat_pass_end(key(i)), c_sat_pass(key(i))
c      enddo
      write (*,*)
      do i = 1, nrpasses
        if (ok(key(i)).ge.1) then
          ihr1 = int(mod(sat_pass_start(i),86400.0D0) / 3600.0D0)
          imn1 =    (mod(sat_pass_start(i),86400.0D0) - dble(ihr1)*3600.0D0)/60.0D0
          rsc1 =     mod(sat_pass_start(i),86400.0D0) - dble(ihr1)*3600.0D0 - imn1 * 60.0D0
          ihr2 = int(mod(sat_pass_end(key(i)),86400.0D0) / 3600.0D0)
          imn2 =    (mod(sat_pass_end(key(i)),86400.0D0) - dble(ihr2)*3600.0D0)/60.0D0
          rsc2 =     mod(sat_pass_end(key(i)),86400.0D0) - dble(ihr2)*3600.0D0 - imn2 * 60.0D0
          write (*,'(I4,F10.2,2(2I5,F10.3,2x),3x,A)') ok(key(i)), sat_el_max_pass(key(i)), ihr1, imn1, rsc1, ihr2, imn2, rsc2, c_sat_pass(key(i))
        endif
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
            call azel_command(sat_pass_start(i), sat_pass_end(key(i)), time_track_la(1,key(i)), sat_az_track_la(1,key(i)), sat_el_track_la(1,key(i)), n_track_la(key(i)), swdebug, swallowflip, swnoovl, ovlflag, azrate, elrate, azoff, eloff, swoffset, azpark, elpark)
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

      subroutine azel(long, lat, ro, az, el)
      implicit none
      real*8 long, lat, ro(3), az, el
c
      real*8 xyz(3), lla(3), xyzsat(3), xyzobs(3), a(3), b(3), c1000, torad, todeg, a_par_b(3), a_ortho_b(3), const
      real*8 mylong, mylat
      integer*2 npos
      character*50 instring
c
      c1000           = 1000.0D0
      torad           = 2.0D0 * dasin(1.0D0) / 180.0D0
      todeg           = 180.0D0 / (2.0D0 * dasin(1.0D0))
      call getenv('LONG',instring)
      call string_to_r8(instring(1:lnblnk(instring)), npos, mylong)
      call getenv('LAT',instring)
      call string_to_r8(instring(1:lnblnk(instring)), npos, mylat)
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
      el = 90.0D0 - todeg * dacos( (a(1)*b(1)+a(2)*b(2)+a(3)*b(3)) / (dsqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3)) * dsqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3)) ) )
c
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
c
      return
      end

      subroutine map_azel_init()
      implicit none
c--
      integer*4 nazel
      parameter (nazel=1400)
c--
      integer*1 mapsline (3,nazel), maps(3,nazel,nazel)
      integer*4 lun, i, j, ix, iy
      real*8    az, el, torad, todeg
      character*6 cmap
      character*11 csize
      character*500 cstring
      character*(*) command
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
      entry map_azel_place(az, el)
      if (el.gt.0.0) then
        ix = (nazel/2) + int(dsin(az*torad)*dble(nazel/2)*(1.0D0-(el/90.0D0)))
c iy is counterintuitive as the image is vertically flipped for display !
        iy = (nazel/2) + int(dcos(az*torad)*dble(nazel/2)*(1.0D0-(el/90.0D0)))
        ix = min(max(ix,1),nazel)
        iy = min(max(iy,1),nazel)
        maps(1,ix,iy) = -1
        maps(2,ix,iy) =  0
        maps(3,ix,iy) = -1
      endif
      return
c
      entry map_azel_write(command)
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

      subroutine azel_command(pass_start, pass_end, time_track, sat_az_track, sat_el_track, n_in, swdebug, swallowflip, swnoovl, ovlflag, azrate, elrate, azoff, eloff, swoffset, azpark, elpark)
      implicit none
      integer*4 n_in
      real*8    pass_start, pass_end, time_track(n_in), sat_az_track(n_in), sat_el_track(n_in), azrate, elrate, azoff, eloff, azpark, elpark
      integer*1 ovlflag(n_in)
      logical   swdebug, swallowflip, swnoovl, swoffset
c
      logical   swflip
      real*8    rsc, time_now, time_sleep, time_start, az1, az2, el1, el2, t_end_move
      integer*4 n, i, j, ihr, imn, dattim(8), doy, i_check
      character*5  c_zone
      character*8  c_date
      character*10 c_time
      character*50 c_command
c
      n = n_in
c-- Apply offset if needed
      if (swoffset) then
        do j = 1, n
          sat_az_track(j) = sat_az_track(j) + azoff
          sat_el_track(j) = sat_el_track(j) + eloff
        enddo
      endif
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
      endif
c--
      call date_and_time(c_date, c_time, c_zone, dattim)
c      write (*,*) dattim(1), dattim(2), dattim(3), dattim(4), dattim(5), dattim(6), dattim(7), dattim(8)
      time_now = dble(dattim(5) * 3600.0 + dattim(6) * 60.0 + dattim(7) + dattim(8)/1000.0)
      time_sleep = pass_start - time_now - 60.0D0
c-- Skip ongoing passes ?
      if (time_sleep.lt.0.0) return
      write (c_command,'(''sleep '',F8.1)') time_sleep
c      write (*,*) c_command(1:lnblnk(c_command))
      call system(c_command)
      j = 1
      ihr = int(time_track(j)/3600.0D0)
      imn = (int(time_track(j) - dble(ihr)*3600.0D0)/60.0D0)
      rsc = time_track(j) - dble(ihr)*3600.0D0 - imn * 60.0D0
      write (*,'(I5,3x,2F10.3,2I5,F10.3)') j, sat_az_track(j), sat_el_track(j), ihr, imn, rsc
      call SPID_MD02(sat_az_track(j), sat_el_track(j), swdebug, swflip)
      do j = 2,n
        call date_and_time(c_date, c_time, c_zone, dattim)
c        write (*,*) dattim(1), dattim(2), dattim(3), dattim(4), dattim(5), dattim(6), dattim(7), dattim(8)
        time_now = dble(dattim(5) * 3600.0 + dattim(6) * 60.0 + dattim(7) + dattim(8)/1000.0)
        time_sleep = time_track(j) - time_now - 0.2
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
          write (*,'(I5,3x,2F10.3,2I5,F10.3)') j, sat_az_track(j), sat_el_track(j), ihr, imn, rsc
          call SPID_MD02(sat_az_track(j), sat_el_track(j), swdebug, swflip)
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
      integer*4 lun, init, i, irec, ios
      real*8 Azuse, Eluse, Azcorr, Elcorr
      character*4 cangle
      data init/0/
      data statbuf/'57'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'1F'x,'20'x/
      data unkbuf /'57'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'00'x,'3F'x,'20'x/
c
      real*8 azrange(2), elrange(2), azlow, azhig, ellow, elhig
      logical swazlim, swellim, swok, swtiltcor
      data swazlim/.false./, swellim/.false./, swtiltcor/.false./
c
      save
c
      if (init.eq.0) then
        init = 1
        irec = 1
        call get_lun(lun)
        open (unit=lun,file='/dev/ttyS4',access='direct',recl=13)
      endif
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
c-- Correct for system tilt - if requested
      if (swtiltcor) then
        Azcorr = Azuse
        Elcorr = Eluse
        call tilt_correct(Azcorr, Elcorr)
        Eluse = Elcorr
      endif
      if (.not.swok) goto 1
      write (*,*) '==> OK'
c--
      write (cangle,'(I4.4)') nint(5.0D0 * (Azuse + 360.0D0))
      devbuf( 1) = '57'x
      devbuf( 2) = ichar(cangle(1:1))
      devbuf( 3) = ichar(cangle(2:2))
      devbuf( 4) = ichar(cangle(3:3))
      devbuf( 5) = ichar(cangle(4:4))
      devbuf( 6) = '05'x
      write (cangle,'(I4.4)') nint(5.0D0 * (Eluse + 360.0D0))
      devbuf( 7) = ichar(cangle(1:1))
      devbuf( 8) = ichar(cangle(2:2))
      devbuf( 9) = ichar(cangle(3:3))
      devbuf(10) = ichar(cangle(4:4))
      devbuf(11) = '05'x
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
      write (*,'('' which enables look ahead tracking of the sun'')')
      write (*,'('' '')')
      write (*,'(''   Second argument is the type  of TLE : resource or weather IF the first argument is not an indirect ! '')')
      write (*,'(''                   NOT required if first argument is reference to indirect list (using @) '')')
      write (*,'('' '')')
      write (*,'(''   Further arguments:   '')')
      write (*,'('' '')')
      write (*,'(''     beam=1.8  as beamwidth of dish to be used '')')
      write (*,'(''     debug     show the commands that would be issued for reactive and pro-active (look ahead) '')')
      write (*,'(''     diag      show the az el values at each time step for reactive and pro-active (look ahead) - L....o....n....g '')')
      write (*,'(''     live      actually command the hardware '')')
      write (*,'(''     flip      allow dish flip mode if crossing the N vector '')')
      write (*,'(''     noovl     removes Az, El commands overlapping in time - hope this cures the MD-02 M2 pulse timeout problem ....... '')')
      write (*,'(''     both      do the calculations with and without look-ahead - only useful with debug or diag '')')
      write (*,'('' '')')
      write (*,'(''     raster=200.5,28.3,3.0,2.0,0.4,0.3,10.0 parameters are Az, El, Az range, El range, Az step, El step, delta-T (sec) - perform a raster scan'')')
      write (*,'('' '')')
      write (*,'(''     azrange=130-290 only allow commands with an Azimuth   in this range - used when running tests at the dish'')')
      write (*,'(''     elrange=0-45    only allow commands with an Elevation in this range - used when running tests at the dish'')')
      write (*,'(''     tiltcor         enable the tilt correction measured with the level inclinometer at 30 deg intervals'')')
      write (*,'(''     azpark=180.0    Override the default park position of Az=209.0'')')
      write (*,'(''     elpark=45.0     Override the default park position of El=0.0'')')
      write (*,'(''     tiltcor   enable the tilt correction measured with the level inclinometer at 30 deg intervals'')')
      write (*,'(''     azoff=-1.0      Define an Az offset - found from pre-calibration - to add to the calculated values'')')
      write (*,'(''     eloff=0.9       Define an El offset - found from pre-calibration - to add to the calculated values'')')
      write (*,'(''     calib     inserts an Az/El calibration routine - assumes you defined azpark, elpark, azoff and eloff !!'')')
      write (*,'('' '')')
      write (*,'('' The TLE environment variable defines where to find the TLE files '')')
c
      return
      end
