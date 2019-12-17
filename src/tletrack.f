      implicit none
c--
c-- nmaxpass is the maximum nr of passes that can be stored over the to be analysed period
c-- nmaxpointing is the maximum nr of timesteps that can be stored per single pass
c-- nsatmax is the maximum nr of satellites that can be analysed in one single run
c
      integer*4 nmaxpass, nmaxpointing, nsatmax
      parameter (nmaxpass=1000, nmaxpointing=10000, nsatmax=50)
c
      real*8    sat_az_pass(nmaxpointing,nmaxpass)    , sat_el_pass(nmaxpointing,nmaxpass)    , time_pass(nmaxpointing,nmaxpass)
      real*8    sat_az_track(nmaxpointing,nmaxpass)   , sat_el_track(nmaxpointing,nmaxpass)   , time_track(nmaxpointing,nmaxpass)
      real*8    sat_az_track_la(nmaxpointing,nmaxpass), sat_el_track_la(nmaxpointing,nmaxpass), time_track_la(nmaxpointing,nmaxpass)
      real*8    sat_pass_start(nmaxpass), sat_pass_end(nmaxpass), sat_el_max_pass(nmaxpass)
      integer*4 npointspass(nmaxpass), n_track(nmaxpass), n_track_la(nmaxpass), i_sat_prio(nsatmax), i_start(nsatmax), i_end(nsatmax)
      integer*4 nrpasses, key(nmaxpass), ok(nmaxpass), sat_pass_prio(nmaxpass)
      character*20 c_sat_pass(nmaxpass), c_sat_use(nsatmax), c_sat_type(nsatmax)
c-- Diagnostic variables - can be deleted if plots no longer needed
      integer*4 i1, i2, i3, i4, i5, i6
      real*8    angle1, angle2
c
      real*8    doy_tle, time, ro(3), vo(3), long, lat, theta0g, time_start, sat_az, sat_el, rsc, sat_el_max, time_max, el_track_min
      real*8    torad, todeg, track_angle, x1, y1, z1, x2, y2, z2, time_save, secsubsteps, time_offset, rsc1, rsc2, elpriothresh
      integer*4 dattim(8), luntle, iyr_tle, doy, i, j, ihr, imn, i_year, i_month, i_day, i_pass, i_track, lunsat, nsatuse, ios, isat
      integer*4 i_arg1, i_arg, ihr1, ihr2, imn1, imn2
      integer*2 npos
      logical foundtle, swdebug, swdiag
      character*5   c_zone
      character*8   c_date, c_datvis
      character*10  c_time
      character*20  c_sat
      character*50  c_file
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
c--
c--   Further arguments:  
c--
c--     beam=1.8  as beamwidth of dish to be used
c
c-- The TLE environment variable defines where to find the TLE files
c
      torad           = 2.0D0 * dasin(1.0D0) / 180.0D0
      todeg           = 180.0D0 / (2.0D0 * dasin(1.0D0))
      el_track_min    = 3.0D0
      track_angle     = 1.8D0
      secsubsteps     = 5.0D0
      swdebug         = .false.
      swdiag          = .false.
c-- If two satellites passes have a conflict and their priorities differ by 1 and their max elevations are above this limit : the max elevation wins
      elpriothresh    = 30.0D0
c
      do i = 1, nmaxpass
        sat_pass_start(i) = 0.0D0
        sat_pass_end(i)   = 0.0D0
        ok(i)             = 1
      enddo
c
      call getarg(1, argstring)
      if (argstring(1:1).eq.'@') then
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
 2      close (unit=lunsat)
        call free_lun(lunsat)
        i_arg1 = 2
      else
        nsatuse = 1
        c_sat_use(nsatuse)  = argstring(1:lnblnk(argstring))
        call getarg(2, argstring)
        c_sat_type(nsatuse) = argstring(1:lnblnk(argstring))
        i_sat_prio(nsatuse) = 9
        i_arg1 = 3
      endif
      do i = i_arg1, 20
        call getarg(i, argstring)
        if (index(argstring(1:lnblnk(argstring)),'debug').ne.0) swdebug = .true.
        if (index(argstring(1:lnblnk(argstring)),'diag').ne.0)  swdiag  = .true.
        if (index(argstring(1:lnblnk(argstring)),'beam=').ne.0) then
          call string_to_r8(argstring(index(argstring,'beam=')+5:lnblnk(argstring)), npos, track_angle)
          write (*,'(''Modified track angle                     : '',F10.3)') track_angle
        endif
      enddo
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
        call map_azel_init()
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
          call run_TLE(time, doy, ro, vo, long, lat, theta0g, doy_tle)
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
            npointspass(nrpasses) = npointspass(nrpasses) + 1
            call doytodate(doy, dattim(1), c_datvis)
            call map_azel_place(sat_az, sat_el)
            ihr = int(time/3600.0D0)
            imn = (int(time - dble(ihr)*3600.0D0)/60.0D0)
            rsc = time - dble(ihr)*3600.0D0 - imn * 60.0D0
            sat_az_pass(npointspass(nrpasses), nrpasses) = sat_az
            sat_el_pass(npointspass(nrpasses), nrpasses) = sat_el
            time_pass  (npointspass(nrpasses), nrpasses) = time + dattim(4) * 60.0D0
c          write (*,'(2F10.3,1x,a8,1x,2I5,F10.3, I5)') sat_az, sat_el, c_datvis, ihr, imn, rsc, doy
          else
            if (i_pass.eq.1) then
              if (sat_pass_end(nrpasses).eq.0.0D0)sat_pass_end(nrpasses) = time  + time_offset + dattim(4) * 60.0D0
            endif
            i_pass = 0
          endif
        enddo
        i_end(isat) = nrpasses
        call sgp4_end()
        call map_azel_write(command)
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
c            write (*,'(''Conflict between passes : '',4I5)')  i, j, key(i), key(j)
            if (sat_pass_prio(key(i)).eq.sat_pass_prio(key(j))) then
              if (sat_el_max_pass(key(i)).gt.sat_el_max_pass(key(j))) then
c              write (*,*) i, key(i), ' Wins '
                ok(key(j)) = 0
                ok(key(i)) = 2
              else
c              write (*,*) j, key(j), ' Wins '
                ok(key(i)) = 0
                ok(key(j)) = 2
              endif
            else
              if (abs(sat_pass_prio(key(i))-sat_pass_prio(key(j))).eq.1.and.sat_el_max_pass(key(i)).gt.elpriothresh.and.sat_el_max_pass(key(j)).gt.elpriothresh) then
                if (sat_el_max_pass(key(i)).gt.sat_el_max_pass(key(j))) then
c              write (*,*) i, key(i), ' Wins '
                  ok(key(j)) = 0
                  ok(key(i)) = 3
                else
c              write (*,*) j, key(j), ' Wins '
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
      stop
      end

      subroutine run_tle(linetime,linedoy, ro, vo, long, lat, theta0g, doy_tle)
      implicit none
      integer*4 linedoy
      real*8 linetime, ro(3), vo(3), long, lat, doy_tle
c
      real*8    Tmfe, theta0g
      save
      Tmfe = ((linetime/86400.0D0) + dble(linedoy) - doy_tle) * 1440.0
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
