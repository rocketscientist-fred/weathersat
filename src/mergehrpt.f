c
c-- Dr Fred A.Jansen - ESA/ESTEC - written for personal use and in private time, using no ESA resources
c-- Use at your own risk !!!
c--
c-- HRPT reading and merging s/w for NOAA, Metop and FengYun .raw16, .hpt and .C10 images
c
c-- Major overhaul ongoing - 30-Dec-2018
c-- In 'current' version I step through the images line by line 3 times to build the miscelaneous histograms 
c-- and combination images.
c
c-- Program usage:
c
c-- ./hrpt.exe filename xyz
c-- filename is one of the FY, Metop or NOAA images in full - can have a directory in front of it
c--
c-- North or South is determined from the TLE calculation
c--
c-- xyz is a 3 digit code for generating the merged ch11 images with x for R, y for G and z for B - e.g. 221 is @petermeteor 's 
c--   standard (and is the default, except for FY where default is 621) - default is invoked e.g. when giving xxx
c-- Parameters 3-20 are optional and have no specific order:
c--   debug  ==> Nr of diagnostic print outs - use >& hrpt.log to route to file
c--   project ==> project the satellite image on the zoomed swath area :options: project=4 ==> ch 4 as channel in B/W image. project=221 use ch2 for R, ch 2 for G and ch 1 for B in projected image. Default w. project is project=221 - requires gapcor
c--               Only in case of FY satellites 0 => ch 10 ! 
c--   border  ==> overplot country borders on projected satellite image - requires project (and project requires gapcor)
c--   gamma=x.y ==> gamma correction for corrected image. Example : gamma=1.0 - default gamma=1.4
c--   theta=x.y ==> Override the thetascan value below for testing purposes (x.y is free format)
c--   linedelay=x.y ==> # of lines of delay introduced on the linetime (can be negative) (x.y is free format but only int(x.y) is internally used)
c--   bestres ==> Places the best resolution pixels in the merged images - use with caution if gaps present
c--
c-- Auxiliary files needed:
c
c--   worldmaphr.rgb - the MODIS 5400 by 2700 pixel map (credit: Retro Stöckli, NASA Earth Observatory) from https://visibleearth.nasa.gov/view.php?id=73580
c--   worldmapuhr.rgb - the MODIS 21600 by 10800 pixel map converted to polar stereographic (credit: Reto Stöckli, NASA Earth Observatory) from https://visibleearth.nasa.gov/view.php?id=73580
c--   color1.tab - color file for the plot software to generate the histograms
c--   color2.tab - color file for the plot software to generate the histograms
c--   color3.tab - color file for the plot software to generate the histograms
c--   color4.tab - color file for the plot software to generate the histograms
c--   weather.txt - a TLE file from CELESTRAK for a date close to when the image file was recorded - need to add a script to download these files
c--
c-- Output files produced:
c
c--   To be described
c
c-- New developments:
c
c
c-- And now for the code:
c-- Top section - until get_lun and free_lun - all written fully by me using wikipedia, NOAA and Eumetsat public documentation
c-- Plot code - some 50% written in the early 80's by Hans Deutekom at the Cosmic Ray Working Group in Leiden - heavily modified by me for postscript
c--             the rest all written by me (pltextbl, plraster, plhist, plpois etc etc)
c-- SGP4 code from the 70's written by David Vallado (see T.S Kelso's celectrak.com) in ancient style fortran and where needed adapted by me.
c-- xyz2lla and lla2xyz taken from https://ea4eoz.blogspot.com/2015/11/simple-wgs-84-ecef-conversion-functions.html
c--             with an implemented improvement from the footnote
c
      implicit none
c
      integer*4 nmaxlin, nmaxpix, nmaxch, nstereo, nmergemax
      parameter (nmaxlin=6000, nmaxpix=2048, nmaxch=10, nstereo=20000000, nmergemax=3)
c
      character*250 filedatain(nmergemax)
      integer*4 nmerge
c
      integer*1 linevalid(nmaxlin), satch_valid(nmaxlin)
      integer*2 satch_image(nmaxpix,nmaxch,nmaxlin,nmergemax), stereo_image(3,nstereo)
      integer*4 satch_linedoy(nmaxlin), linedoy(nmaxlin)
      real*8    satch_linetime(nmaxlin), linetime(nmaxlin), longlat(3,2049,nmaxlin,nmergemax), satch_az(nmaxlin), satch_el(nmaxlin)
c
      integer*1 inbuf(27728)
      integer*2 channel_image(2048,10), imgbuf(2048), rgbbuf(3,2048), rgbbuf_correct(3,5000), pixcol(3), npos
      integer*4 lun, irec, i, j, k, l, lunimg, ios, k1, k2, j1, j2, i_merge, n_merge, ilim_low(nmaxch,nmergemax), ilim_high(nmaxch,nmergemax)
      integer*4 visch(4), nvisch, jch1(4), jch2(4), xchoff(10), nsat_pix
      integer*4 nchan, nrecl, doy, i_arg, lunborder, nxs, nys, ilinedelay, luninput, ilinedelay_in(nmergemax), doyfile, doy_off
      integer*4 nblank, nlines_gapcor(nmergemax), k3, nwrite, i_off, irp, igp, ibp, nvalid, npix1, npix2, perc1, perc2, hist1(1024), hist2(1024)
      real*4    sat_avg_alt, sat_fov, lat_start, lat_end
      real*4    dr, dg, db, dy, di, dq, sum1, sum2
      real*8    timestamp, timstart, timprev, ro(3), vo(3), long, lat, long_stereo, lat_stereo
      real*8    theta0g, thetascan, theta_in(nmergemax), w1, w2, w3, rlinedelay
      logical   swnorth, swfy, swnoaa, swmetop, swgapcor, swswath, swdebug, swproject, swborder, swgamma, swtheta, swlinedelay, swbestres
      character*4 cfy, cmetop
      character*6 cnoaa
      character*9 cgamma
      character*20  ctheta, clinedelay
      character*250 filedata, argstring, imgname, command, outstring
      call get_command(command)
      write (*,'(a)') command(1:lnblnk(command))
      if (index(command,'--help').ne.0) then
        write (*,'(a)') 'Program usage:'
        write (*,'(a)') './hrpt.exe filename dir '
        write (*,'(a)') '           filename is one of the FY, Metop or NOAA images in full - can have a directory in front of it, should at least start with ./ or the like '
        write (*,'(a)') ''
        write (*,'(a)') 'North or South is determined from the TLE calculation'
        write (*,'(a)') ''
        write (*,'(a)') '           xyz is a 3 digit code for generating the merged ch11 images with x for R, y for G and z for B - e.g. 221 is @petermeteor  '
        write (*,'(a)') '             standard (and is the default, except for FY where default is 621) - default is invoked e.g. when giving xxx'
        write (*,'(a)') ''
        write (*,'(a)') 'Parameters 2-20 are optional and have no specific order:'
        write (*,'(a)') ''
        write (*,'(a)') '  debug  ==> Nr of diagnostic print outs - use >& hrpt.log to route to file'
        write (*,'(a)') ''
        write (*,'(a)') '  project ==> project the satellite image on the zoomed swath area :options: project=4 ==> ch 4 as channel in B/W image. project=221 use ch2 for R, ch 2 for G and ch 1 for B in projected image. '
        write (*,'(a)') '              Only in case of FY satellites 0 => ch 10 ! Default w. project is project=221 - requires gapcor'        
        write (*,'(a)') ''
        write (*,'(a)') '  border  ==> overplot country borders on projected satellite image - requires project (and project requires gapcor)'
        write (*,'(a)') ''
        write (*,'(a)') '  gamma=x.y ==> gamma correction for corrected image. Example : gamma=1.0 - default gamma=1.4'
        write (*,'(a)') ''
        write (*,'(a)') '  theta=x.y ==> Override the thetascan value below for testing purposes (x.y is free format)'
        write (*,'(a)') ''
        write (*,'(a)') '  linedelay=x.y ==> # of lines of delay introduced on the linetime (can be negative) (x.y is free format but only int(x.y) is internally used)'
        write (*,'(a)') ''
        write (*,'(a)') '  bestres ==> Places the best resolution pixels in the merged images - use with caution if gaps present'
        write (*,'(a)') ''
        write (*,'(a)') 'Auxiliary files needed:'
        write (*,'(a)') ''
        write (*,'(a)') '  worldmaphr.rgb - the MODIS 5400 by 2700 pixel map (credit: Retro Stöckli, NASA Earth Observatory) from https://visibleearth.nasa.gov/view.php?id=73580'
        write (*,'(a)') '  worldmapuhr.rgb - the MODIS 21600 by 10800 pixel map converted to polar stereographic (credit: Reto Stöckli, NASA Earth Observatory) from https://visibleearth.nasa.gov/view.php?id=73580'
        write (*,'(a)') '  color1.tab - color file for the plot software to generate the histograms'
        write (*,'(a)') '  color2.tab - color file for the plot software to generate the histograms'
        write (*,'(a)') '  color3.tab - color file for the plot software to generate the histograms'
        write (*,'(a)') '  color4.tab - color file for the plot software to generate the histograms'
        write (*,'(a)') '  weather.txt - a TLE file from CELESTRAK for a date close to when the image file was recorded - currently in QNAP crontab as : '
        write (*,'(a)') ''
        write (*,'(a)') '#!/bin/bash'
        write (*,'(a)') 'mytlefile=/share/Public/TLE/`(date +%Y%m%d)`-weather.tle'
        write (*,'(a)') 'wget --secure-protocol=auto https://www.celestrak.com/NORAD/elements/weather.txt'
        write (*,'(a)') 'mv weather.txt $mytlefile'
        write (*,'(a)') ''
        write (*,'(a)') 'Output files produced:'
        write (*,'(a)') ''
        write (*,'(a)') '  Per channel B/W images : for NOAA/Metop ch1.png through ch5.png 16 bit images using the full 10 bit dynamic range'
        write (*,'(a)') '                           for FengYun ch1.png through chA.png 16 bit images using the full 10 bit dynamic range'
        write (*,'(a)') '  Two channel colour images using the GOES LUT (might suppress this in future)'
        write (*,'(a)') '  For FengYun only "true" colour images based on the four visible channel scans translated center wavelength to RGB'
        write (*,'(a)') '    and merged into one image - need to determine the weight'
        write (*,'(a)') '  RGB image ch11.png using for RGB the channels from the xyz command line parameter'
        write (*,'(a)') '  RGB image ch11_corrected.png using for RGB the channels from the xyz command line parameter and corrected to a fixed pixel size on the earths surface'
        write (*,'(a)') '  satmap.png image showing the image swath grayed out with a purple subsatellite line over the modis world map'
        write (*,'(a)') ''
        write (*,'(a)') 'Environment variables used :'
        write (*,'(a)') 'export TLE=/home/fjansen/psa/TLE/               ==> Points to the TLE repository'
        write (*,'(a)') 'export HRPTOUT=/home/fjansen/Win7/SDR/Output/   ==> Points to the where the plots, images and data files are to be stored'
        stop
      endif
c
      do i = 1, nmergemax
        theta_in(i)      = 0.0
        ilinedelay_in(i) = 0
      enddo
      swnorth   = .false.
c-- Arguments from 4 onwards have no fixed order
      swgapcor  = .true.
      swdebug   = .false.
      swproject = .false.
      irp       = 2
      igp       = 2
      ibp       = 1
      swborder  = .false.
      swgamma   = .false.
      swtheta   = .false.
      swlinedelay = .false.
      swbestres = .false.
      do i_arg = 2,20
        call getarg(i_arg, argstring)
        if (index (argstring,'debug').ne.0)  swdebug  = .true.
        if (index (argstring,'project').ne.0) then
          swproject = .true.
          if (index(argstring,'=').ne.0) then
            i = index(argstring,'=')
            if (lnblnk(argstring)-i.eq.1) then
              if (ichar('0').le.ichar(argstring(i+1:i+1)).and.ichar(argstring(i+1:i+1)).le.ichar('9')) then
                read (argstring(i+1:i+1),'(i1)') irp
                igp = irp
                ibp = irp
              endif
            else if (lnblnk(argstring)-i.eq.3) then
              if (ichar('0').le.ichar(argstring(i+1:i+1)).and.ichar(argstring(i+1:i+1)).le.ichar('9').and.
     *            ichar('0').le.ichar(argstring(i+2:i+2)).and.ichar(argstring(i+2:i+2)).le.ichar('9').and.
     *            ichar('0').le.ichar(argstring(i+3:i+3)).and.ichar(argstring(i+3:i+3)).le.ichar('9')) then
                read (argstring(i+1:i+1),'(i1)') irp
                read (argstring(i+2:i+2),'(i1)') igp
                read (argstring(i+3:i+3),'(i1)') ibp
              endif
            endif
          else
            irp = 2
            igp = 2
            ibp = 1
          endif
        endif
        if (index (argstring,'border').ne.0)  swborder  = .true.
        if (index (argstring,'gamma=').ne.0) then
          swgamma  = .true.
          cgamma   = argstring(1:9)
          i = index(cgamma,'=')
          cgamma(i:i) = ' '
        endif
        if (index (argstring,'bestres').ne.0)  swbestres  = .true.
      enddo
c-- Borders require project
      if (swborder) swproject = .true.
c-- Polar stereograhic projection requires gapcor
      if (swproject) swgapcor = .true.
c
c-- Loop over the input file(s)
c
      call get_lun(luninput)
      call getarg(1, argstring)
      open (unit=luninput, file=argstring(1:lnblnk(argstring)), form='formatted')
      n_merge = 0
      ios     = 0
      do while (ios.eq.0) 
        read (luninput,'(a)',iostat=ios) argstring
        if (ios.ne.0) goto 102
        i = 1
        do while (argstring(i:i).eq.' '.and.i.le.len(argstring))
          i = i + 1
        enddo
        k1 = i
        do while (argstring(i:i).ne.' '.and.i.le.len(argstring))
          i = i + 1
        enddo
        k2 = i - 1
        n_merge = n_merge + 1
        do i = 1, len(filedatain(n_merge))
          filedatain(n_merge)(i:i) = ' '
        enddo
        filedatain(n_merge) = argstring(k1:k2)
c
        if (index(argstring,'theta=').ne.0) then
          k1 = index(argstring,'theta=') + 6
          i  = k1
          do while (argstring(i:i).ne.' '.and.i.le.len(argstring))
            i = i + 1
          enddo
          k2 = i - 1
          swtheta  = .true.
          call string_to_r8(argstring(k1:k2)//' ',npos, theta_in(n_merge))
        endif
c
        if (index(argstring,'linedelay=').ne.0) then
          k1 = index(argstring,'linedelay=') + 10
          i  = k1
          do while (argstring(i:i).ne.' '.and.i.le.len(argstring))
            i = i + 1
          enddo
          k2 = i - 1
          swlinedelay = .true.
          call string_to_r8(argstring(k1:k2)//' ',npos, rlinedelay)
          ilinedelay_in(n_merge) = int(rlinedelay)
        endif
      enddo
 102  close (unit=luninput)
      call free_lun(luninput)
c
      do i = 1, len(outstring)
        outstring(i:i) = ' '
      enddo
      call getenv('HRPTOUT',outstring)
c
c-- Start the loop over the input files !
c-- Been lazy, don't want to indent all the lines by two so it's a poor man's loop
c
      i_merge = 0
 100  i_merge = i_merge + 1
      if (i_merge.gt.n_merge) goto 101
      j = 0
      i = lnblnk(filedatain(i_merge))
      do while (i.ge.1.and.filedatain(i_merge)(i:i).ne.'/')
        i = i - 1
      enddo
      j = i + 1
      filedata = outstring(1:lnblnk(outstring))//filedatain(i_merge)(j:lnblnk(filedatain(i_merge)))
c
c-- Set the mission dependent variables
c
      swfy    = .false.
      swnoaa  = .false.
      swmetop = .false.
      if (index(filedatain(i_merge),'FY').ne.0) then
        nchan = 10
        swfy  = .true.
        cfy   = '.C10'
        if (irp.eq.0) irp = 10
        if (igp.eq.0) igp = 10
        if (ibp.eq.0) ibp = 10
        nrecl = 27728
c-- Initialise the projection correction arrays
c-- Remember that the IFOV (= pixel size) is BIGGER than the sampling interval !
        sat_avg_alt  = 833.0
        sat_fov      = 55.37 * 2.0
        nsat_pix     = 2048
        thetascan    = 0.00D0
c-- Define the FY3B and FY3C channels to be used for trucolor
        nvisch   = 4
        visch(1) = 1
        visch(2) = 7
        visch(3) = 8
        visch(4) = 9
c-- Define the FY3B and FY3C x offsets ! Used "Razende Bol" in Anders 16/Oct/2018 FY3B 1408 image to align
c-- IR channels 3, 4, 5 and 10 not yet calibrated - poor spatial resolution
        xchoff(1) =  0
        xchoff(2) =  0
        xchoff(3) =  0
        xchoff(4) =  0
        xchoff(5) =  0
c-- Ch 6 also appears to have (as one of the few) a +1 Y-offset
        xchoff(6) = +2
        xchoff(7) = -4
        xchoff(8) = +3
        xchoff(9) = -2
        xchoff(10)=  0
      endif
      if (index(filedatain(i_merge),'NOAA').ne.0) then
        nchan  = 5
        swnoaa = .true.
        cnoaa  = '.raw16'
        nrecl  = 22180
c-- Initialise the projection correction arrays
c-- Remember that the IFOV (= pixel size) is BIGGER than the sampling interval !
        sat_avg_alt  = 833.0
        sat_fov      = 55.37 * 2.0
        nsat_pix     = 2048
        nvisch       = 1
        thetascan    = 0.00D0
c-- Define the NOAA x offsets ! 
        xchoff(1) =  0
        xchoff(2) =  0
        xchoff(3) =  0
        xchoff(4) =  0
        xchoff(5) =  0
      endif
      if (index(filedatain(i_merge),'Metop').ne.0.or.index(filedatain(i_merge),'M02').ne.0) then
        nchan   = 5
        swmetop = .true.
        cmetop  = '.hpt'
        nrecl   = 13864
c-- Initialise the projection correction arrays
c-- Remember that the IFOV (= pixel size) is BIGGER than the sampling interval !
        sat_avg_alt  = 833.0
        sat_fov      = 55.37 * 2.0
        nsat_pix     = 2048
        nvisch       = 1
        thetascan    = -2.70D0
c        thetascan    = 0.0D0
c-- Define the Metop x offsets ! 
        xchoff(1) =  0
        xchoff(2) =  0
        xchoff(3) =  0
        xchoff(4) =  0
        xchoff(5) =  0
      endif
c-- Theta override
      ilinedelay = 0
      if (swlinedelay) ilinedelay = ilinedelay_in(i_merge)
c
      if (i_merge.eq.1) then
        do i = 1, nmaxlin
          satch_valid(i) = 0
          do l = 1, nmergemax
            do j = 1,nmaxch
              do k = 1,nmaxpix
                satch_image(k, j, i, l) = 0
              enddo
            enddo
            do j = 1, nmaxpix + 1
              longlat(1,j,i,l) = 0.0D0
              longlat(2,j,i,l) = 0.0D0
              longlat(3,j,i,l) = 0.0D0
            enddo
          enddo
        enddo
      endif
c
      call get_lun(lun)
      open (unit=lun,file=filedatain(i_merge)(1:lnblnk(filedatain(i_merge))),form='unformatted',access='direct',recl=nrecl)
c
c-- First loop through the image to retrieve the timestamp per scanline (needed to reject lines later)
c
      irec = 1
      do while (.true.)
        ios = 0
        if (swfy)    call get_scanline(lun, irec, ios, channel_image, 2048, nchan, 'FY', inbuf, nrecl, timestamp, doy, swdebug)
        if (swnoaa)  call get_scanline(lun, irec, ios, channel_image, 2048, nchan, 'NO', inbuf, nrecl, timestamp, doy, swdebug)
        if (swmetop) call get_scanline(lun, irec, ios, channel_image, 2048, nchan, 'ME', inbuf, nrecl, timestamp, doy, swdebug)
        if (ios.ne.0) goto 2
        linetime(irec) = timestamp
        linedoy(irec)  = doy
        irec = irec + 1
      enddo
 2    irec = irec - 1
      call check_linetimes(linetime, linedoy, linevalid, irec, timstart, swdebug)
c-- Allow for the TLE being off in time just a bit. Specified in # of lines.
      if (swlinedelay) then
        do i = 1, irec
          if (linevalid(i).eq.1) then
            linetime(i) = linetime(i) + dble(ilinedelay) * 0.1667D0
          endif
        enddo
      endif
c
c-- Init the TLE code, initialise OLD SGP4 code and initialise the worldmap code
c
      i = 1
      do while (linevalid(i).ne.1.and.i.le.irec)
        i = i + 1
      enddo
      if (swfy)    call init_TLE('FY', filedatain(i_merge)(1:lnblnk(filedatain(i_merge))), linedoy(i), doyfile)
      if (swnoaa)  call init_TLE('NO', filedatain(i_merge)(1:lnblnk(filedatain(i_merge))), linedoy(i), doyfile)
      if (swmetop) call init_TLE('ME', filedatain(i_merge)(1:lnblnk(filedatain(i_merge))), linedoy(i), doyfile)
c-- I have seen the in file doy on NOAA being wrong - correct with date in filename
      if (swnoaa) then
        if (doyfile.ne.linedoy(i)) then
          write (*,'(''ERROR doy in data  : '',I5,'' is NOT doy in filename : '',I5,'' using the latter ! '')') linedoy(i), doyfile
          doy_off = doyfile - linedoy(i)
          do j = 1, irec
            if (linevalid(i).eq.1) linedoy(j) = linedoy(j) + doy_off
          enddo
          call init_TLE('NO', filedatain(i_merge)(1:lnblnk(filedatain(i_merge))), linedoy(i), doyfile)
        endif
      endif
      call sgp4_init(ro, vo, long, lat)
      if (i_merge.eq.1) then
        call map_init()
        call map_azel_init()
        call map_stereo_init()
      endif
c-- Map out the subsatellite point and then the swath
c      swswath = .true.
      lat_start = 0.0
      lat_end   = 0.0
      do i = 1, irec
        if (linevalid(i).eq.1) then
          call run_TLE(linetime(i),linedoy(i), ro, vo, long, lat, theta0g)
          call map_mark(long, lat, 0)
          call map_stereo_mark(long, lat, 0)
          if (lat_start.eq.0.0) lat_start = lat
          lat_end = lat
        endif
      enddo
      if (lat_start.lt.lat_end) swnorth = .true.
      if (lat_start.gt.lat_end) swnorth = .false.
      write (*,'(''Latitude start     : '',F8.2,'' Latitude end           : '',F8.2,'' swnorth : '',L1)') lat_start, lat_end, swnorth
      if (swnorth) thetascan = - thetascan
c-- Theta override
      if (swtheta) thetascan = theta_in(i_merge)
c-- Second loop through the image(s) - write the channel images with or without gapcor and fill the satch_image (and related) taking gaps into account
      timprev                = 0.0D0
      nlines_gapcor(i_merge) = 0      
      irec                   = 1
      do while (.true.)
        ios = 0
        if (swfy)    call get_scanline(lun, irec, ios, channel_image, 2048, nchan, 'FY', inbuf, nrecl, timestamp, doy, swdebug)
        if (swnoaa)  call get_scanline(lun, irec, ios, channel_image, 2048, nchan, 'NO', inbuf, nrecl, timestamp, doy, swdebug)
        if (swmetop) call get_scanline(lun, irec, ios, channel_image, 2048, nchan, 'ME', inbuf, nrecl, timestamp, doy, swdebug)
        if (ios.ne.0) goto 1
c-- store timref of previous valid irec and calc delta-lines from that. If 1 just write. If more add blank lines !!
        if (swgapcor) then
          if (linevalid(irec).eq.1) then
            if (timprev.eq.0.0D0) then
              timprev = linetime(irec)
            else
              nblank = nint((linetime(irec)-timprev)/0.1667D0) - 1
              if (nblank.gt.0) then
                do k3 = 1, nblank
                  nlines_gapcor(i_merge) = nlines_gapcor(i_merge) + 1
                  do j = 1,nchan
                    do i = 1,2048
                      imgbuf(i) = 0
                      satch_image(i,j,nlines_gapcor(i_merge),i_merge) = 0
                    enddo
                  enddo
                  satch_linetime(nlines_gapcor(i_merge)) = timprev + dble(k3) * 0.1667D0
                  satch_linedoy(nlines_gapcor(i_merge))  = linedoy(irec)
                  swswath                       = .false.
                  call run_TLE(satch_linetime(nlines_gapcor(i_merge)),satch_linedoy(nlines_gapcor(i_merge)), ro, vo, long, lat, theta0g)
                  call scan_track(ro, vo, long, lat, satch_linetime(nlines_gapcor(i_merge)),satch_linedoy(nlines_gapcor(i_merge)), longlat(1,1,nlines_gapcor(i_merge),i_merge), swswath, theta0g, thetascan, swnorth, swtheta)
                enddo
              endif
              timprev = linetime(irec)
            endif
            nlines_gapcor(i_merge) = nlines_gapcor(i_merge) + 1
            do j = 1,nchan
              do i = 1,2048
                if (swnorth)      l = 2049 - i
                if (.not.swnorth) l = i
                if (0.le.channel_image(i,j).and.channel_image(i,j).le.1023) then
                  imgbuf(l) = channel_image(i,j)
                else
                  imgbuf(l) = 0
                endif
                satch_image(l,j,nlines_gapcor(i_merge),i_merge) = imgbuf(l)
              enddo
            enddo
            satch_linetime(nlines_gapcor(i_merge)) = linetime(irec)
            satch_linedoy(nlines_gapcor(i_merge))  = linedoy(irec)
            satch_valid(nlines_gapcor(i_merge))    = 1
            swswath                       = .true.
            call run_TLE(satch_linetime(nlines_gapcor(i_merge)),satch_linedoy(nlines_gapcor(i_merge)), ro, vo, long, lat, theta0g)
            call scan_track(ro, vo, long, lat, satch_linetime(nlines_gapcor(i_merge)),satch_linedoy(nlines_gapcor(i_merge)), longlat(1,1,nlines_gapcor(i_merge),i_merge), swswath, theta0g, thetascan, swnorth, swtheta)
            call azel(long, lat, ro, satch_az(nlines_gapcor(i_merge)), satch_el(nlines_gapcor(i_merge)))
            call map_azel_place(satch_az(nlines_gapcor(i_merge)), satch_el(nlines_gapcor(i_merge)))
            call map_stereo_place_azel(satch_az(nlines_gapcor(i_merge)), satch_el(nlines_gapcor(i_merge)))
          endif
        endif
        irec = irec + 1
      enddo
 1    continue
      if (swdebug) then
        write (*,*) '** Azimuth & Elevation for observer **'
        do i = 1, nlines_gapcor(i_merge)
          if (satch_valid(i).eq.1) write (*,*) i, satch_linetime(i), longlat(1,1025,i,i_merge), longlat(2,1025,i,i_merge), satch_az(i), satch_el(i)
        enddo
      endif
c-- Release the logical Unit nrs = lun
      close (unit=lun)
      call free_lun(lun)
      call sgp4_close()
      goto 100
 101  continue
c
c-- Find the zoom parameters for the polar stereographic projection and execute the projection
c-- First zero out the image area
c
      do i = 1, nstereo
        stereo_image(1,i) = 0
        stereo_image(2,i) = 0
        stereo_image(3,i) = 0
      enddo
      if (swproject) then
        call map_stereo_zoom_init(2048,irec,nxs,nys)
        if (nxs*nys.gt.nstereo) stop ' ** Enlarge stereo Image in source code (nstereo) ** '
      endif
c--
      if (swproject) then
        i_off = 1
        do l = 1, n_merge
          call get_hist_limits(satch_image(1,1,1,l), 2048, 10, 8000, ilim_low(1,l), ilim_high(1,l), nchan, nlines_gapcor(l))
        enddo
        do l = 1, n_merge
          do i = 1, nlines_gapcor(l)
            if (satch_valid(i).eq.1) then
              do k = 1,2048
                if (ilim_low(irp,l).le.satch_image(k,irp,i,l).and.satch_image(k,irp,i,l).le.ilim_high(irp,l).and.
     *              ilim_low(igp,l).le.satch_image(k,igp,i,l).and.satch_image(k,igp,i,l).le.ilim_high(igp,l).and.
     *              ilim_low(ibp,l).le.satch_image(k,ibp,i,l).and.satch_image(k,ibp,i,l).le.ilim_high(ibp,l))     then
                  if (l.eq.1) then
                    pixcol(1) = satch_image(k,irp,i,l)
                    pixcol(2) = satch_image(k,igp,i,l)
                    pixcol(3) = satch_image(k,ibp,i,l)
                  else
                    pixcol(1) = ilim_low(irp,1) + int(float(satch_image(k,irp,i,l) - ilim_low(irp,l)) * float(ilim_high(irp,1) - ilim_low(irp,1)) / float(ilim_high(irp,l) - ilim_low(irp,l)))
                    pixcol(2) = ilim_low(igp,1) + int(float(satch_image(k,igp,i,l) - ilim_low(igp,l)) * float(ilim_high(igp,1) - ilim_low(igp,1)) / float(ilim_high(igp,l) - ilim_low(igp,l)))
                    pixcol(3) = ilim_low(ibp,1) + int(float(satch_image(k,ibp,i,l) - ilim_low(ibp,l)) * float(ilim_high(ibp,1) - ilim_low(ibp,1)) / float(ilim_high(ibp,l) - ilim_low(ibp,l)))
                  endif
                  do j = 1,8
                    long_stereo = longlat(1,k,i,l) + dble(j-1) * (longlat(1,k+1,i,l)-longlat(1,k,i,l)) / 8.0D0
                    lat_stereo  = longlat(2,k,i,l) + dble(j-1) * (longlat(2,k+1,i,l)-longlat(2,k,i,l)) / 8.0D0
                    if (swbestres) then
                      call map_stereo_zoom_place(stereo_image,3,nxs,nys,long_stereo, lat_stereo, pixcol, max(abs(k-1024),1))
                    else
                      call map_stereo_zoom_place(stereo_image,3,nxs,nys,long_stereo, lat_stereo, pixcol, 1)
                    endif
                  enddo
c-- Separately place the original value
                  long_stereo = (longlat(1,k,i,l)+longlat(1,k+1,i,l)) / 2.0D0
                  lat_stereo  = (longlat(2,k,i,l)+longlat(2,k+1,i,l)) / 2.0D0
                  if (swbestres) then
                    call map_stereo_zoom_place(stereo_image,3,nxs,nys,long_stereo, lat_stereo, pixcol, max(abs(k-1024),1))
                  else
                    call map_stereo_zoom_place(stereo_image,3,nxs,nys,long_stereo, lat_stereo, pixcol, 1)
                  endif
                else if (satch_image(k,irp,i,l).gt.ilim_high(irp,l).or.satch_image(k,igp,i,l).gt.ilim_high(igp,l).or.satch_image(k,ibp,i,l).gt.ilim_high(ibp,l)) then
                  if (satch_image(k,irp,i,l).gt.ilim_high(irp,l)) pixcol(1) = ilim_high(irp,l)
                  if (satch_image(k,igp,i,l).gt.ilim_high(igp,l)) pixcol(2) = ilim_high(igp,l)
                  if (satch_image(k,ibp,i,l).gt.ilim_high(ibp,l)) pixcol(3) = ilim_high(ibp,l)
                  do j = 1,8
                    long_stereo = longlat(1,k,i,l) + dble(j-1) * (longlat(1,k+1,i,l)-longlat(1,k,i,l)) / 8.0D0
                    lat_stereo  = longlat(2,k,i,l) + dble(j-1) * (longlat(2,k+1,i,l)-longlat(2,k,i,l)) / 8.0D0
                    if (swbestres) then
                      call map_stereo_zoom_place(stereo_image,3,nxs,nys,long_stereo, lat_stereo, pixcol, max(abs(k-1024),1))
                    else
                      call map_stereo_zoom_place(stereo_image,3,nxs,nys,long_stereo, lat_stereo, pixcol, 1)
                    endif
                  enddo
c-- Separately place the original value
                  long_stereo = (longlat(1,k,i,l)+longlat(1,k+1,i,l)) / 2.0D0
                  lat_stereo  = (longlat(2,k,i,l)+longlat(2,k+1,i,l)) / 2.0D0
                  if (swbestres) then
                    call map_stereo_zoom_place(stereo_image,3,nxs,nys,long_stereo, lat_stereo, pixcol, max(abs(k-1024),1))
                  else
                    call map_stereo_zoom_place(stereo_image,3,nxs,nys,long_stereo, lat_stereo, pixcol, 1)
                  endif
                endif
              enddo
            endif
          enddo
        enddo
        if (swborder) then
          call get_borders(stereo_image,3,nxs,nys)
          call get_lun(lunborder)
          if (swfy)    k = index(filedata,cfy)
          if (swnoaa)  k = index(filedata,cnoaa)
          if (swmetop) k = index(filedata,cmetop)
          write (imgname,'(a,a)') filedata(1:k-1),'_border.rgb'
          open (unit=lunborder,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=3*nxs*2)
          close (unit=lunborder,status='delete')
          open (unit=lunborder,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=3*nxs*2)
          do i = 1, nys
            call write_buffer(lunborder,i,stereo_image(1,i_off),3*nxs)
            i_off = i_off + nxs
          enddo
          close (unit=lunborder)
          call free_lun(lunborder)
c-- This projection always requires a flip !
          if (swgamma) then
            if (swfy)    call create_png(filedata,cfy   ,'.rgb','_border.rgb',0, '-'//cgamma, .true., nxs, -nys)
            if (swnoaa)  call create_png(filedata,cnoaa ,'.rgb','_border.rgb',0, '-'//cgamma, .true., nxs, -nys)
            if (swmetop) call create_png(filedata,cmetop,'.rgb','_border.rgb',0, '-'//cgamma, .true., nxs, -nys)
          else
            if (swfy)    call create_png(filedata,cfy   ,'.rgb','_border.rgb',0,' ', .true., nxs, -nys)
            if (swnoaa)  call create_png(filedata,cnoaa ,'.rgb','_border.rgb',0,' ', .true., nxs, -nys)
            if (swmetop) call create_png(filedata,cmetop,'.rgb','_border.rgb',0,' ', .true., nxs, -nys)
          endif
        else
          call get_lun(lunborder)
          if (swfy)    k = index(filedata,cfy)
          if (swnoaa)  k = index(filedata,cnoaa)
          if (swmetop) k = index(filedata,cmetop)
          write (imgname,'(a,a)') filedata(1:k-1),'_project.rgb'
          open (unit=lunborder,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=3*nxs*2)
          close (unit=lunborder,status='delete')
          open (unit=lunborder,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=3*nxs*2)
          do i = 1, nys
            call write_buffer(lunborder,i,stereo_image(1,i_off),3*nxs)
            i_off = i_off + nxs
          enddo
          close (unit=lunborder)
          call free_lun(lunborder)
c-- This projection always requires a flip !
          if (swgamma) then
            if (swfy)    call create_png(filedata,cfy   ,'.rgb','_project.rgb',0, '-'//cgamma, .true., nxs, -nys)
            if (swnoaa)  call create_png(filedata,cnoaa ,'.rgb','_project.rgb',0, '-'//cgamma, .true., nxs, -nys)
            if (swmetop) call create_png(filedata,cmetop,'.rgb','_project.rgb',0, '-'//cgamma, .true., nxs, -nys)
          else
            if (swfy)    call create_png(filedata,cfy   ,'.rgb','_project.rgb',0,' ', .true., nxs, -nys)
            if (swnoaa)  call create_png(filedata,cnoaa ,'.rgb','_project.rgb',0,' ', .true., nxs, -nys)
            if (swmetop) call create_png(filedata,cmetop,'.rgb','_project.rgb',0,' ', .true., nxs, -nys)
          endif
        endif
      endif
c
c-- End of the projection/border part
c
c
      call map_write()
      if (swgapcor) call map_stereo_write(command)
      if (swgapcor) call map_azel_write(command)
      if (swfy)    k = index(filedata,cfy)
      if (swnoaa)  k = index(filedata,cnoaa)
      if (swmetop) k = index(filedata,cmetop)
      write (imgname,'(a,a)') filedata(1:k-1),'_satmap.png'
      call system('mv satmap.png '//imgname(1:lnblnk(imgname)))
      write (*,'(a)') ' mv satmap.png '//imgname(1:lnblnk(imgname))
      if (swgapcor) then
        write (imgname,'(a,a)') filedata(1:k-1),'_satmap-stereo.png'
        write (*,'(a)') ' mv satmap-stereo.png '//imgname(1:lnblnk(imgname))
        call system('mv satmap-stereo.png '//imgname(1:lnblnk(imgname)))
        write (imgname,'(a,a)') filedata(1:k-1),'_satmap-stereo-zoom.png'
        write (*,'(a)') ' mv satmap-stereo-zoom.png '//imgname(1:lnblnk(imgname))
        call system('mv satmap-stereo-zoom.png '//imgname(1:lnblnk(imgname)))
        write (imgname,'(a,a)') filedata(1:k-1),'_azel.png'
        write (*,'(a)') ' mv satmap-azel.png '//imgname(1:lnblnk(imgname))
        call system('mv satmap-azel.png '//imgname(1:lnblnk(imgname)))
      endif      
c
      nvalid = 0
      do i = 1, irec
        if (linevalid(i).eq.1) nvalid = nvalid + 1
      enddo
      write (*,'(''# of records : '',I8,'' of which : '',I8,'' in satellite data file and : '',I8,'' with a valid (in sequence) timestamp '')') nwrite, irec, nvalid
c
      stop
      end

      subroutine get_scanline(lun, irec, ios, channel_image, npix, nchan, cmis, buf, nbuf, timestamp, doy, swdebug)
      implicit none
      real*8    timestamp
      integer*4 lun, irec, ios, npix, nchan, nbuf, doy
      integer*2 channel_image(npix, 10)
      integer*1 buf(nbuf)
      logical   swdebug
      character*2 cmis
c
      integer*4 i, j, k, l, k1, k2, ihr, imn, init
      real*4    rsc
      real*8    offset
c
      data init/0/
c
      save
c
      if (cmis.eq.'FY') then
        read (lun,rec=irec,iostat=ios) buf
        if (ios.ne.0) then
          init = 1
          goto 1
        endif
        k = 2001
c-- Approach 1 : try what David Taylor said - continuous 10 bit packing, i.e. 4 pixels in 5 bytes
c-- Looked wrong at first, but the start value of k is critical - 2001 appears to work for FY-3
        j = 1
        i = 1
        do while (.true.)
          do l = 1,2
            k1 = buf(k)
            if (k1.lt.0) k1 = k1 + 256
            k1 = k1 * 4
            k2 = ishft(iand(buf(k+1),'C0'x),-6)
            channel_image(i,j) = k1 + k2
            j = j + 1
            if (j.eq.nchan+1) then
              if (i.eq.npix) goto 1
              j = 1
              i = i + 1
            endif
            k1 = iand(buf(k+1),'3F'x)
            k1 = k1 * 16
            k2 = ishft(iand(buf(k+2),'F0'x),-4)
            channel_image(i,j) = k1 + k2
            j = j + 1
            if (j.eq.nchan+1) then
              if (i.eq.npix) goto 1
              j = 1
              i = i + 1
            endif
            k1 = iand(buf(k+2),'0F'x)
            k1 = k1 * 64
            k2 = ishft(iand(buf(k+3),'FC'x),-2)
            channel_image(i,j) = k1 + k2
            j = j + 1
            if (j.eq.nchan+1) then
              if (i.eq.npix) goto 1
              j = 1
              i = i + 1
            endif
            k1 = iand(buf(k+3),'03'x)
            k1 = k1 * 256
            k2 = buf(k+4)
            if (k2.lt.0) k2 = k2 + 256
            channel_image(i,j) = k1 + k2
            j = j + 1
            if (j.eq.nchan+1) then
              if (i.eq.npix) goto 1
              j = 1
              i = i + 1
            endif
            k = k + 5
          enddo
        enddo
 1      continue
        call HRPT_time_FyMe(buf(12),4,timestamp)
        if (timestamp.lt.43200.0D0) offset =  43200.0D0
        if (timestamp.ge.43200.0D0) offset = -43200.0D0
        timestamp = timestamp + offset
        doy = buf(11)
        if (doy.le.0) doy = doy + 256
        doy = doy * 2 + ishft(iand(buf(12),'80'x),-7)
        ihr = int (timestamp / 3600.0D0)
        imn = int((timestamp - dble(ihr)*3600.0D0)/60.0D0)
        rsc = timestamp - dble(ihr) * 3600.0D0 - dble(imn) * 60.0D0
        if (init.eq.0.and.swdebug) write (*,'(I4,2X,24(Z2.2,2x),F20.4,3I5,F8.3)') irec, (buf(j), buf(j+1), j = 1, 23, 2), timestamp, doy, ihr, imn, rsc
      endif
      if (cmis.eq.'NO') then
        read (lun,rec=irec,iostat=ios) buf
        if (ios.ne.0) then
          init = 1
          goto 2
        endif
c        k = 751 - Is the old Integer*2 offset 
        k = 1502
        do i = 1,2048
          do j = 1,nchan
            channel_image(i,j) = buf(k) * 256
            if (buf(k-1).lt.0) channel_image(i,j) = channel_image(i,j) + buf(k-1) + 256
            if (buf(k-1).ge.0) channel_image(i,j) = channel_image(i,j) + buf(k-1)
            k = k + 2
          enddo
        enddo
        call HRPT_time_NOAA(buf(19),6,timestamp)
        doy = iand(buf(18),'03'x) * 128
        doy = doy + ishft(iand(buf(17),'FD'x),-1)
        ihr = int (timestamp / 3600.0D0)
        imn = int((timestamp - dble(ihr)*3600.0D0)/60.0D0)
        rsc = timestamp - dble(ihr) * 3600.0D0 - dble(imn) * 60.0D0
        if (init.eq.0.and.swdebug) write (*,'(I4,2X,24(Z2.2,2x),F20.4,3I5,F8.3)') irec, (buf(j), buf(j+1), j = 1, 23, 2), timestamp, doy, ihr, imn, rsc
 2      continue
      endif
      if (cmis.eq.'ME') then
        read (lun,rec=irec,iostat=ios) buf
        if (ios.ne.0) then
          init = 1
          goto 3
        endif
        k = 951
c-- Approach 1 : try what David Taylor said - continuous 10 bit packing, i.e. 4 pixels in 5 bytes
c-- Looked wrong at first, but the start value of k is critical - 941 appears to work
        j = 1
        i = 1
        do while (.true.)
          k1 = buf(k)
          if (k1.lt.0) k1 = k1 + 256
          k1 = k1 * 4
          k2 = ishft(iand(buf(k+1),'C0'x),-6)
          channel_image(i,j) = k1 + k2
          j = j + 1
          if (j.eq.nchan+1) then
            if (i.eq.npix) goto 3
            j = 1
            i = i + 1
          endif
          k1 = iand(buf(k+1),'3F'x)
          k1 = k1 * 16
          k2 = ishft(iand(buf(k+2),'F0'x),-4)
          channel_image(i,j) = k1 + k2
          j = j + 1
          if (j.eq.nchan+1) then
            if (i.eq.npix) goto 3
            j = 1
            i = i + 1
          endif
          k1 = iand(buf(k+2),'0F'x)
          k1 = k1 * 64
          k2 = ishft(iand(buf(k+3),'FC'x),-2)
          channel_image(i,j) = k1 + k2
          j = j + 1
          if (j.eq.nchan+1) then
            if (i.eq.npix) goto 3
            j = 1
            i = i + 1
          endif
          k1 = iand(buf(k+3),'03'x)
          k1 = k1 * 256
          k2 = buf(k+4)
          if (k2.lt.0) k2 = k2 + 256
          channel_image(i,j) = k1 + k2
          j = j + 1
          if (j.eq.nchan+1) then
            if (i.eq.npix) goto 3
            j = 1
            i = i + 1
          endif
          k = k + 5
        enddo
 3      continue
        call HRPT_time_FyMe(buf(12),4,timestamp)
        doy = buf(11)
        if (doy.le.0) doy = doy + 256
        doy = doy * 2 + ishft(iand(buf(12),'80'x),-7)
        ihr = int (timestamp / 3600.0D0)
        imn = int((timestamp - dble(ihr)*3600.0D0)/60.0D0)
        rsc = timestamp - dble(ihr) * 3600.0D0 - dble(imn) * 60.0D0
        if (init.eq.0.and.swdebug) write (*,'(I4,2X,24(Z2.2,2x),F20.4,3I5,F8.3)') irec, (buf(j), buf(j+1), j = 1, 23, 2), timestamp, doy, ihr, imn, rsc
      endif
c
      return
      end

      subroutine HRPT_time_FyMe(buf,n,timestamp)
      implicit none
      integer*4 n
      integer*1 buf(n)
      real*8    timestamp
c
      integer*4 i, b
      integer*8 multiplier, time_msec
c
      multiplier = 1
      time_msec  = 0
      do i = n, 1, -1
       if (i.eq.1) then
         b = iand(buf(i),'07'x)
         time_msec  = time_msec + b * multiplier
       else
         b = buf(i)
         if (b.lt.0) b = b + 256
         time_msec  = time_msec + b * multiplier
         multiplier = multiplier * 256
        endif
      enddo
      timestamp = dble(time_msec) / 1000.0D0
c
      return
      end

      subroutine HRPT_time_NOAA(buf,n,timestamp)
      implicit none
      integer*4 n
      integer*1 buf(n)
      real*8    timestamp
c
      integer*4 b
      integer*8 time_msec
c
      time_msec  = 0
      b          = buf(5)
      if (b.lt.0) b = b + 256
      time_msec  = time_msec + (iand(buf(6),'03'x) * 256 + b)
      b          = buf(3)
      if (b.lt.0) b = b + 256
      time_msec  = time_msec + (iand(buf(4),'03'x) * 256 + b) * 1024
      time_msec  = time_msec + iand(buf(1),'7F'x) * 1024 * 1024
      timestamp = dble(time_msec) / 1000.0D0
c
      return
      end      

      subroutine check_linetimes(linetime, linedoy, linevalid, irec, timstart, swdebug)
c-- Determine the minimum and maximum valid time and check 'validity' ..... difficult algorithm and NOT bullet proof
      implicit none
      integer*4 irec
      integer*1 linevalid(irec)
      logical   swdebug
      integer*4 linedoy(irec)
      real*8    linetime(irec), timstart
c
      integer*4  i, j, timhist(864), i_tmax, v_tmax, i_min, i_max, ilinelast, doyhist(365), i_valid, irecval(3)
      real*8     delta_1, delta_2, timeval(3), lastvalid
c
      do i = 1,864
        timhist(i) = 0
      enddo
      do i = 1,365
        doyhist(i) = 0
      enddo
      do i = 1,irec
        linevalid(i) = -1
      enddo
c-- Build histogram of times in scanlines, find peak and search for 0 at either side with a +/- max 20 minutes delta
      do i = 1, irec
        j = int(linetime(i) / 100.0D0) + 1
        if (1.le.j.and.j.le.864) timhist(j) = timhist(j) + 1
      enddo
      i_tmax = 0
      v_tmax = 0
      do i = 1, 864
        if (timhist(i).gt.v_tmax) then
          v_tmax = timhist(i)
          i_tmax = i
        endif
      enddo
      i_min = 0
      i_max = 0
      do i = 1, 12
        if (timhist(i_tmax-i).eq.0.and.i_min.eq.0) i_min = i_tmax - i + 1
        if (timhist(i_tmax+i).eq.0.and.i_max.eq.0) i_max = i_tmax + i - 1
      enddo
      if (i_min.eq.0) i_min = i_tmax - 12 + 1
      if (i_max.eq.0) i_max = i_tmax - 12 + 1
c-- Find three scanlines in a row which are correct and separated by one scanline delta time (0.166666666 sec) - thus create an initial linevalid seed
      do i = 2, irec-1
        if (dble(i_min-1)*100.0D0.le.linetime(i-1).and.linetime(i-1).lt.dble(i_max)*100.0D0.and.
     *      dble(i_min-1)*100.0D0.le.linetime(i  ).and.linetime(i  ).lt.dble(i_max)*100.0D0.and.
     *      dble(i_min-1)*100.0D0.le.linetime(i+1).and.linetime(i+1).lt.dble(i_max)*100.0D0)     then
          delta_1 = (linetime(i)   - linetime(i-1)) / 0.1666666D0
          delta_2 = (linetime(i+1) - linetime(i))   / 0.1666666D0
          if (dabs(delta_1-1.0D0).lt.5.0D-3.and.dabs(delta_2-1.0D0).lt.5.0D-3) then
            linevalid(i-1) = 1
            linevalid(i)   = 1
            linevalid(i+1) = 1
          endif
        endif
      enddo
c-- Find the first valid line after the previous screening
      i_valid = 0
      do i = 1, irec
        if (linevalid(i).eq.1.and.i_valid.eq.0) i_valid = i
      enddo
c-- Now go through all the lines with invalid linetime - use delta linetime as the filter with an absolute 0.05 sec criterion - refine on # of lines ?
      do i = 1,irec
        if (linevalid(i).eq.-1) then
          if (i.lt.i_valid) then
            delta_1 = (linetime(i_valid) - linetime(i)) / 0.1666666D0
            if (delta_1.gt.0.and.abs(delta_1-dble(nint(delta_1))).le.5.0D-2) then
              linevalid(i) = 1
            endif
          else
            delta_1 = (linetime(i) - linetime(i_valid)) / 0.1666666D0
            if (delta_1.gt.0.and.abs(delta_1-dble(nint(delta_1))).le.5.0D-2) then
              linevalid(i) = 1
            endif
          endif
        else
          i_valid = i
        endif
      enddo
c-- Some time values outside of the peak could have made it through, refilter
      do i = 1, irec
        if (linetime(i).lt.dble(i_min-1)*100.0D0.or.linetime(i).gt.dble(i_max)*100.0D0) linevalid(i) = -1
      enddo
c-- Check for the most appearing day-of-year (with a histogram) and filter on it +/- one doy (cater for midnight passes) either side
      do i = 1, irec
        if (linedoy(i).ge.1.and.linedoy(i).le.365) then
          doyhist(linedoy(i)) = doyhist(linedoy(i)) + 1
        endif
      enddo
      i_min = 0
      do i = 1, 365
        if (doyhist(i).gt.i_min) then
          i_min = doyhist(i)
          i_max = i
        endif
      enddo
      do i = 1, irec
        if (linedoy(i).lt.i_max-1.or.linedoy(i).gt.i_max+1) linevalid(i) = -1
      enddo
c-- Silly algorithm checking for an ascending sequence of times - only works for single/isolated outliers
      j = 1
      do i = 1, irec
        if(linevalid(i).eq.1) then
          timeval(j) = linetime(i)
          irecval(j) = i
          j = j + 1
          if (j.gt.3) then
            if ((timeval(2)-timeval(1)).lt.0.0D0.or.(timeval(3)-timeval(2)).lt.0.0D0) then
              linevalid(irecval(2)) = -1
              j = 3
              timeval(2) = timeval(3)
              irecval(2) = irecval(3)
            else
              j = 3
              timeval(1) = timeval(2)
              timeval(2) = timeval(3)
              irecval(1) = irecval(2)
              irecval(2) = irecval(3)
            endif
          endif
        endif
      enddo
c-- Determine the first valid timestamp
      timstart = 0.0D0
      do i = 1, irec
        if (linevalid(i).eq.1.and.timstart.eq.0.0D0) then 
          timstart = linetime(i)
        endif
      enddo
c--
      if (swdebug) then
        write (*,'(''Valid Lines '')')
        lastvalid = 0.0D0
        do i = 1, irec
          if (linevalid(i).eq.1) then 
            write (*,*) i, linedoy(i), linetime(i), linevalid(i), linetime(i) - lastvalid
            lastvalid = linetime(i)
          endif
        enddo
      endif
      return
      end
      
      subroutine init_TLE(csat,file,linedoy,doyfile)
      implicit none
      integer*4     linedoy, doyfile
      character*2   csat
      character*(*) file
c
      logical   foundtle
      integer*4 lun, iyr_tle, i, iyr_file, imn_file, ida_file
      real*8    doy_tle, Tmfe, linetime, ro(3), vo(3), long, lat
      real*8    theta0g, thetadum
      character*8 cdate
      character*69 tlestring
      character*250 tledir
      save
c-- Auto seek TLE - tricky stuff !
      call getenv('TLE',tledir)
      if (csat.eq.'FY') then
        i = index(file,'FY3') + 5
        read (file(i:i+3),'(I4)') iyr_file
      endif
      if (csat.eq.'NO') then
        i = index(file,'HRPT_') + 5
        read (file(i  :i+3),'(I4)') iyr_file
        read (file(i+4:i+5),'(I2)') imn_file
        read (file(i+6:i+7),'(I2)') ida_file
        call datetodoy(iyr_file, imn_file, ida_file, doyfile)
      endif
      if (csat.eq.'ME') then
        if (index(file,'_M02').ne.0) then
          i = index(file,'_M02') - 15
          read (file(i:i+3),'(I4)') iyr_file
        else
          i = index(file,'MetopB')
          if (i.eq.0) i = index(file,'MetopC')
          i = i + 7
          read (file(i:i+3),'(I4)') iyr_file
        endif
      endif
      write (*,'(''TLE Seek year      : '',I5)') iyr_file
      call doytodate(linedoy, iyr_file, cdate)
      foundtle = .false.
      call exist(tledir(1:lnblnk(tledir))//cdate//'-weather.tle', foundtle)
      if (foundtle) then
        call system('cp '//tledir(1:lnblnk(tledir))//cdate//'-weather.tle'//' ./weather.txt')
        write (*,'(''Using TLE file     : '',A)') tledir(1:lnblnk(tledir))//cdate//'-weather.tle'
      else
        call doytodate(linedoy-1, iyr_file, cdate)
        foundtle = .false.
        call exist(tledir(1:lnblnk(tledir))//cdate//'-weather.tle', foundtle)
        if (foundtle) then
          call system('cp '//tledir(1:lnblnk(tledir))//cdate//'-weather.tle'//' ./weather.txt')
          write (*,'(''Using TLE file     : '',A)') tledir(1:lnblnk(tledir))//cdate//'-weather.tle'
        else
          call doytodate(linedoy+1, iyr_file, cdate)
          foundtle = .false.
          call exist(tledir(1:lnblnk(tledir))//cdate//'-weather.tle', foundtle)
          if (foundtle) then
            call system('cp '//tledir(1:lnblnk(tledir))//cdate//'-weather.tle'//' ./weather.txt')
            write (*,'(''Using TLE file     : '',A)') tledir(1:lnblnk(tledir))//cdate//'-weather.tle'
          else
            write (*,*) 'TLE not found - using current weather.txt'
          endif
        endif
      endif
c
      if (csat.eq.'FY') then
        if (index(file,'3A').ne.0) call system('cat weather.txt | grep -A 2 -e "FENGYUN 3A" - > hrpt.tmp')
        if (index(file,'3B').ne.0) call system('cat weather.txt | grep -A 2 -e "FENGYUN 3B" - > hrpt.tmp')
        if (index(file,'3C').ne.0) call system('cat weather.txt | grep -A 2 -e "FENGYUN 3C" - > hrpt.tmp')
      endif
      if (csat.eq.'NO') then
        if (index(file,'NOAA-').ne.0) then
          if (index(file,'NOAA-15').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 15" - > hrpt.tmp')
          if (index(file,'NOAA-18').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 18" - > hrpt.tmp')
          if (index(file,'NOAA-19').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 19" - > hrpt.tmp')
        else if (index(file,'NOAA_').ne.0) then
          if (index(file,'NOAA_15').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 15" - > hrpt.tmp')
          if (index(file,'NOAA_18').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 18" - > hrpt.tmp')
          if (index(file,'NOAA_19').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 19" - > hrpt.tmp')
        endif
c        if (file(index(file,'NOAA')+4:index(file,'NOAA')+4).eq.'-') then
c          if (index(file,'NOAA-15').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 15" - > hrpt.tmp')
c          if (index(file,'NOAA-18').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 18" - > hrpt.tmp')
c          if (index(file,'NOAA-19').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 19" - > hrpt.tmp')
c        else if (file(index(file,'NOAA')+4:index(file,'NOAA')+4).eq.'_') then
c          if (index(file,'NOAA_15').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 15" - > hrpt.tmp')
c          if (index(file,'NOAA_18').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 18" - > hrpt.tmp')
c          if (index(file,'NOAA_19').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 19" - > hrpt.tmp')
c        endif
      endif
      if (csat.eq.'ME') then
        if (index(file,'opA').ne.0.or.index(file,'_M02').ne.0) call system('cat weather.txt | grep -A 2 -e "METOP-A" - > hrpt.tmp')
        if (index(file,'opB').ne.0) call system('cat weather.txt | grep -A 2 -e "METOP-B" - > hrpt.tmp')
        if (index(file,'opC').ne.0) call system('cat weather.txt | grep -A 2 -e "METOP-C" - > hrpt.tmp')
      endif
      call system('tail -n +2 hrpt.tmp > hrpt.tle')
      call system('rm hrpt.tmp')
      call get_lun(lun)
      open (unit=lun,file='hrpt.tle',form='formatted')
      read (lun,'(a)') tlestring
      rewind (lun)
      read (tlestring(19:20),'(i2)') iyr_tle
      iyr_tle = iyr_tle + 2000
      read (tlestring(21:32),*) doy_tle
      close (unit=lun)
      call free_lun(lun)
      write (*,'(''TLE Reference date : '',I5,2x,F8.2)') iyr_tle, doy_tle
c
      return
c
      entry run_tle(linetime,linedoy, ro, vo, long, lat, theta0g)
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

c-- 17-Jan-2019 - modified for WGS84
      subroutine scan_track(ro, vo, long, lat, linetime, linedoy, longlat, swswath, theta0g, thetascan, swnorth, swtheta)
      implicit none
      real*8    ro(3), vo(3), long, lat, linetime, longlat(3,2049), thetascan
      integer*4 linedoy
      logical swswath, swnorth, swtheta
c
      real*8 torad, todeg, x1, y1, z1, d1, d2, d3,  c1, c2, c3, a, b, c, xs1, ys1, zs1, angle_save
      real*8 xv1, yv1, zv1, xv2, yv2, zv2, angle,  longss, latss, angle_threshold, ystep, long_save, lat_save
      real*8 long_v, lat_v, c1000, a_r, b_r, xyz(3), lla(3), theta0g, twopi, phicor, angle1, angle2
      real*8 ro_copy(3), vo_copy(3), d3br, c21, c22, c31, c32, c33, h, onepi, d3i, scandist, depsilon
      real*8 aa(3), ba(3), ca(3), theta, long_int, lat_int, linetime_ref, erad
      real*8 lastlongplus, lastlatplus, lastlongmin, lastlatmin
      integer*4 ind_pix, ind_pix_save, isgn, nstep, nstep_save, i, j, init
      data init/0/
c
      save
c
c-- The FOV of 55.37 deg is the half width of the band scanned as seen from the satellite - 0.95 mrad sample step ! (ie .lt. IFOV)
      angle_threshold = 55.37D0
      scandist        = 2.0D0 * angle_threshold / 2048.0D0
      depsilon        = 1.0D0 * scandist
      a_r             = 6378137.0D0
      b_r             = 6356752.315D0
      c1000           = 1000.0D0
      torad           = 2.0D0 * dasin(1.0D0) / 180.0D0
      todeg           = 180.0D0 / (2.0D0 * dasin(1.0D0))
      twopi           = 4.0D0 * dasin(1.0D0)
      onepi           = 2.0D0 * dasin(1.0D0)
c
      do i = 1, 2049
        do j = 1,3
          longlat(j,i) = 0.0D0
        enddo
      enddo
      longlat(1,1025) = long
      longlat(2,1025) = lat
c-- Copy ro and vo for local use
      ro_copy(1) = ro(1) * c1000
      ro_copy(2) = ro(2) * c1000
      ro_copy(3) = ro(3) * c1000
      vo_copy(1) = vo(1) * c1000 + ro_copy(1)
      vo_copy(2) = vo(2) * c1000 + ro_copy(2)
      vo_copy(3) = vo(3) * c1000 + ro_copy(3)
c-- Rotate ro back to start on the Greenwich meridian - does not matter for lat - add vo to ro (line before) and rotate back by same angle.
c-- then subtract rotated ro and you have the vo at Greenwich meridian for normal to plane (d1, d2, d3)
      call xyz2lla(ro_copy, lla)
      phicor = lla(1)
      lla(1) = 0.0D0
      call lla2xyz( lla, ro_copy)
c
      call xyz2lla(vo_copy, lla)
      lla(1) = lla(1) - phicor
      call lla2xyz( lla, vo_copy)
      vo_copy(1) = vo_copy(1) - ro_copy(1)
      vo_copy(2) = vo_copy(2) - ro_copy(2)
      vo_copy(3) = vo_copy(3) - ro_copy(3)
c
      x1     = ro_copy(1)
      y1     = ro_copy(2)
      z1     = ro_copy(3)
      d1     = vo_copy(1)
      d2     = vo_copy(2)
      d3     = vo_copy(3)
c-- Put the satellite vector on Earth's surface by setting latitude (= lla(3) ) to zero
      xyz(1) = x1
      xyz(2) = y1
      xyz(3) = z1
      call xyz2lla(xyz, lla)
      lla(3) = 0.0D0
      call lla2xyz(lla, xyz)
      x1     = xyz(1)
      y1     = xyz(2)
      z1     = xyz(3)
c-- Experimental - Rotate the vector d around the vector from nadir point on earth to satellite
      aa(1)  = d1
      aa(2)  = d2
      aa(3)  = d3
c-- ba is THE vector from Satellite to Nadir point ! this is the vector to compare with for binning the scanline
      ba(1)  = ro_copy(1) - x1
      ba(2)  = ro_copy(2) - y1
      ba(3)  = ro_copy(3) - z1
      theta  = thetascan
c
c-- Try a default rotation of the nadir plane by using the earth rotation angle per scanline
c-- This basically is earth circumference at latitude divide by 86400 and 6 to get meters per scanline
c-- Then divide by 1100 and take arc tangent to get angle
c-- This is enabled by theta=0.0 on the command line
c
      if(dabs(theta).lt.0.001D0.and.swtheta) then
        erad = dsqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2))
        theta = - datan(twopi*erad/86400.0D0/6.0D0/1100.0D0) * todeg
        if (swnorth) theta = -theta
        write (*,*) linetime_ref, linetime, theta, long, lat
      endif
c--
      call rotate_vector_a_around_b(aa, ba, theta, ca)
      d1     = ca(1)
      d2     = ca(2)
      d3     = ca(3)
c-- Calculate h with the new vector and x1, y1, z1
      h      = d1 * x1 + d2 * y1 + d3 * z1
c
      d3br   =  1.0D0 / ( d3 * d3 * b_r * b_r)
      d3i    =  1.0D0 / d3
      c1     = (1.0D0 / (a_r * a_r)) + d1 * d1 * d3br
      c21    =  2.0D0 * d1 * d2 * d3br
      c22    =  2.0D0 * d1 * h  * d3br
      c31    = (1.0D0 / (a_r * a_r)) + d2 * d2 * d3br
      c32    =  2.0d0 * d2 * h * d3br
      c33    =  h * h * d3br - 1.0D0
c-  Shift Y - this should be OK as I rotate back the true satellite position to long = 0 and later rotate the solution back
      isgn         = 1
      ystep        = 10.0D0
c-- Silly construct to restart at original stepsize after sign reversal
 2    ystep        = dsign(10.0D0, ystep)
      ys1          = y1 - ystep
      angle        = 0.0D0
      ind_pix_save = -1
      nstep        = 0
      nstep_save   = 0
      do while (angle.le.angle_threshold + depsilon)
        a = c1
        b = c21 * ys1 - c22
        c = ys1 * ( c31 * ys1 - c32)  + c33
        if (b*b-4.0D0*a*c.gt.0.0D0) then
          xs1 = (-b + dsqrt(b*b-4.0D0*a*c)) / 2.0D0 / a
          zs1 = d3i * (h - d1 * xs1 - d2 * ys1)
c-- Check the angle between the satellite vector and the surface point, the smallest is the valid one
          xv1 = ro_copy(1)
          yv1 = ro_copy(2)
          zv1 = ro_copy(3)
          angle  = (xv1 * xs1 + yv1 * ys1 + zv1 * zs1) / dsqrt(xv1*xv1+yv1*yv1+zv1*zv1) / dsqrt(xs1*xs1+ys1*ys1+zs1*zs1)
          angle1 = todeg * acos(angle)
          xs1 = (-b - dsqrt(b*b-4.0D0*a*c)) / 2.0D0 / a
          zs1 = d3i * (h - d1 * xs1 - d2 * ys1)
          angle  = (xv1 * xs1 + yv1 * ys1 + zv1 * zs1) / dsqrt(xv1*xv1+yv1*yv1+zv1*zv1) / dsqrt(xs1*xs1+ys1*ys1+zs1*zs1)
          angle2 = todeg * acos(angle)
          if (angle2.lt.angle1) then
            xs1 = (-b - dsqrt(b*b-4.0D0*a*c)) / 2.0D0 / a
            zs1 = d3i * (h - d1 * xs1 - d2 * ys1)
          else
            xs1 = (-b + dsqrt(b*b-4.0D0*a*c)) / 2.0D0 / a
            zs1 = d3i * (h - d1 * xs1 - d2 * ys1)
          endif
          xv2 = ro_copy(1) - xs1
          yv2 = ro_copy(2) - ys1
          zv2 = ro_copy(3) - zs1
          angle  = (ba(1) * xv2 + ba(2) * yv2 + ba(3) * zv2) / dsqrt(ba(1)*ba(1)+ba(2)*ba(2)+ba(3)*ba(3)) / dsqrt(xv2*xv2+yv2*yv2+zv2*zv2)
          angle  = todeg * acos(angle)
          if (angle.le.angle_threshold + depsilon) then
            xyz(1) = xs1
            xyz(2) = ys1
            xyz(3) = zs1
            call xyz2lla(xyz, lla)
            if (lla(1).lt.0.0D0) lla(1) = lla(1) + twopi
            latss  =  lla(2)                     * todeg
            longss = (lla(1) + phicor - theta0g) * todeg
            if (longss.lt.-180.0D0) longss = longss + 360.0D0
            if (longss.gt.+180.0D0) longss = longss - 360.0D0
c-- Find pixel boundaries !
            ind_pix = int(angle/scandist)
c-- Remove this if part if I find the cause for the initial big jump
            if (init.eq.0.and.ind_pix_save.eq.-1) then
              if (ystep.gt.0.0D0) write (*,'(''Temporary diagnostic in scan_track :'')')
              if (ystep.gt.0.0D0) write (*,'(''    ind_pix    angle   scandist      x1          xs1        diff         y1          ys1        diff         z1          zs1        diff        nstep'')')
              write (*,'(I10,2F10.4,9F12.2,I8)') ind_pix, angle, scandist, x1, xs1, abs(x1-xs1), y1, ys1, abs(y1-ys1), z1, zs1, abs(z1-zs1), nstep
              if (ystep.lt.0.0D0) init = 1
            endif
            if (ind_pix_save.ne.-1.and.ind_pix.ne.ind_pix_save) then
c-- Why does this happen ?? e.g. NOAA on 30-Dec-2018 (from Anders)
c              if (ind_pix-ind_pix_save.gt.2) then
c                write (*,*) linetime, ind_pix-ind_pix_save, ind_pix, ind_pix_save, longss, latss
c              endif
              call linear_interpolate(angle_save/scandist, angle/scandist, long_save, longss, long_v)
              call linear_interpolate(angle_save/scandist, angle/scandist,  lat_save,  latss,  lat_v)
              if (ind_pix.le.1024) then
                longlat(1,1025 - isgn * ind_pix) = long_v
                longlat(2,1025 - isgn * ind_pix) = lat_v
c                if (swswath) call darkness(longss, latss, linetime, linedoy, longlat(3,1025 - isgn * ind_pix))
              endif
              if (nstep-nstep_save.ge.10) then
                ystep = ystep * 1.414214
              endif
              nstep_save   = nstep
              ind_pix_save = ind_pix
              long_save    = longss
              lat_save     = latss
              angle_save   = angle
            else
              ind_pix_save = ind_pix
              long_save    = longss
              lat_save     = latss
              angle_save   = angle
            endif
c            if (swswath) then
c              call map_mark(longss, latss, 1)
c              call map_stereo_mark(longss, latss, 1)
c            endif
         endif
        else
          goto 1
        endif
        ys1   = ys1 - ystep
        nstep = nstep + 1
      enddo
 1    if (ystep.gt.0.0) then
        ystep = - ystep
        isgn  = - isgn
        goto 2
      endif
c
c-- Capture the problem that there are "empty" bins 
c-- Investigate and/or improve by interpolation - assume itś numerical problems.
      lastlongplus = longlat(1,1025)
      lastlatplus  = longlat(2,1025)
      lastlongmin  = longlat(1,1025)
      lastlatmin   = longlat(2,1025)
      do i = 1,1024
        if (longlat(1,1025-i).eq.0.0D0) then
          longlat(1,1025-i) = lastlongmin
        else
          lastlongmin  = longlat(1,1025-i)
        endif
        if (longlat(1,1025+i).eq.0.0D0) then
          longlat(1,1025+i) = lastlongplus
        else
          lastlongplus = longlat(1,1025+i)
        endif
        if (longlat(2,1025-i).eq.0.0D0) then
          longlat(2,1025-i) = lastlatmin
        else
          lastlatmin  = longlat(2,1025-i)
        endif
        if (longlat(2,1025+i).eq.0.0D0) then
          longlat(2,1025+i) = lastlatplus
        else
          lastlatplus = longlat(2,1025+i)
        endif
      enddo
      if (swswath) then
        do i = 1, 2048
c          do j = 1,5
c            long_int = longlat(1,i) + dble(j-1) * (longlat(1,i+1)-longlat(1,i)) / 5.0D0
c            lat_int  = longlat(2,i) + dble(j-1) * (longlat(2,i+1)-longlat(2,i)) / 5.0D0
c            call map_mark(long_int, lat_int, 1)
c            call map_stereo_mark(long_int, lat_int, 1)
c          enddo
          call map_mark((longlat(1,i)+longlat(1,i+1))/2.0D0, (longlat(2,i)+longlat(2,i+1))/2.0D0, 1)
          call map_stereo_mark((longlat(1,i)+longlat(1,i+1))/2.0D0, (longlat(2,i)+longlat(2,i+1))/2.0D0, 1)
        enddo
      endif
      return
      end

      subroutine rotate_vector_a_around_b(a, b, theta, c)
      implicit none
      real*8 a(3), b(3), c(3), theta, const, x1, x2
c
      real*8 a_par_b(3), a_ortho_b(3), a_ortho_b_theta(3), w(3)
      real*8 torad, angle, rnorm_ortho
c
      torad = 2.0D0 * dasin(1.0D0) / 180.0D0
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
      rnorm_ortho  = dsqrt(a_ortho_b(1) * a_ortho_b(1) + a_ortho_b(2) * a_ortho_b(2) + a_ortho_b(3) * a_ortho_b(3))
      x1           = dcos(angle) / rnorm_ortho
      x2           = dsin(angle) / dsqrt(w(1) * w(1) + w(2) * w(2) + w(3) * w(3))
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

      subroutine azel(long, lat, ro, az, el)
      implicit none
      real*8 long, lat, ro(3), az, el
c
      real*8 xyz(3), lla(3), xyzsat(3), xyzobs(3), a(3), b(3), c1000, torad, todeg, a_par_b(3), a_ortho_b(3), const
c
      c1000           = 1000.0D0
      torad           = 2.0D0 * dasin(1.0D0) / 180.0D0
      todeg           = 180.0D0 / (2.0D0 * dasin(1.0D0))
c
      xyz(1) = ro(1) * c1000
      xyz(2) = ro(2) * c1000
      xyz(3) = ro(3) * c1000
      call xyz2lla(xyz, lla)
      lla(1) = long * torad
      lla(2) = lat  * torad
      call lla2xyz(lla, xyzsat)
      lla(1) =  4.536008 * torad
      lla(2) = 52.173131 * torad
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
      
      subroutine create_png(file, ext1, ext2, ext3, n, option, swnorth, npix, nlines_in)
      implicit none
      logical   swnorth
      integer*4 n, npix, nlines_in
      character*(*) file, ext1, ext2, ext3, option
c
      integer*4 j, nlines, i
      character*1 csq, cdq
      character*500 imgname, command, cinput, cstring
c
      if (nlines_in.lt.0) then
        call get_command(cinput)
        csq = char(ichar("'"))
        cdq = char(ichar('"'))
        cstring = ' -pointsize 20 -fill red -draw '//csq//'text 50,25 '//cdq//cinput(1:lnblnk(cinput))//cdq//csq
      else
        do i = 1, len(cstring)
          cstring(i:i) = ' '
        enddo
      endif
      nlines = abs(nlines_in)
      j = index(file,ext1(1:lnblnk(ext1)))
      if (n.ne.0) then
        write (imgname,'(a,a,z1,a)') file(1:j-1),ext3(1:lnblnk(ext3)),n,ext2(1:lnblnk(ext2))
      else
        write (imgname,'(a,a)') file(1:j-1),ext3(1:lnblnk(ext3))
      endif
      j = index(imgname,ext2(1:lnblnk(ext2)))
      if (swnorth) then
        if (index(ext2,'.rgb').eq.0) then
          write (command,'(''convert -depth 16 -endian lsb -size '',I4.4,''x'',I4.4,1x,a,'' -equalize '',a,'' -flip -depth 16 '',a,'' '',a)') 
     *      npix, nlines, 'gray:'//imgname(1:lnblnk(imgname)), option(1:lnblnk(option)), cstring(1:lnblnk(cstring)), imgname(1:j-1)//'.png'
        else
          write (command,'(''convert -depth 16 -endian lsb -size '',I4.4,''x'',I4.4,1x,a,'' -equalize '',a,'' -flip -depth 16 '',a,'' '',a)') 
     *      npix, nlines, 'rgb:'//imgname(1:lnblnk(imgname)), option(1:lnblnk(option)), cstring(1:lnblnk(cstring)), imgname(1:j-1)//'.png'
        endif
      else
        if (index(ext2,'.rgb').eq.0) then
          write (command,'(''convert -depth 16 -endian lsb -size '',I4.4,''x'',I4.4,1x,a,'' -equalize '',a,'' -depth 16 '',a,'' '',a)') 
     *      npix, nlines, 'gray:'//imgname(1:lnblnk(imgname)), option(1:lnblnk(option)), cstring(1:lnblnk(cstring)), imgname(1:j-1)//'.png'
        else
          write (command,'(''convert -depth 16 -endian lsb -size '',I4.4,''x'',I4.4,1x,a,'' -equalize '',a,'' -depth 16 '',a,'' '',a)') 
     *      npix, nlines, 'rgb:'//imgname(1:lnblnk(imgname)), option(1:lnblnk(option)), cstring(1:lnblnk(cstring)), imgname(1:j-1)//'.png'
        endif
      endif
      write (*,*) command(1:lnblnk(command))
      call system(command)
c
      return
      end
      
      subroutine map_init()
      implicit none
c--
      integer*4 nx, ny
      parameter (nx=5400, ny=2700)
c--
      integer*1 mapline (3,nx), map(3,nx,ny), mask(nx, ny)
      integer*2 pix(3)
      integer*4 lun, i, j, ix, iy, itype, value
      real*8    long, lat
      character*6 cmap
      character*9 csize
c
      save
c
      cmap = 'satmap' 
      write (csize,'(I4.4,''x'',I4.4)') nx, ny
c
      call get_lun(lun)
      open (unit=lun, file='worldmaphr.rgb', access='direct', form='unformatted', recl=3*nx)
      do i = 1, ny
        read (lun,rec=i) mapline
        do j = 1, nx
          map(1,j,i) = mapline (1,j)
          map(2,j,i) = mapline (2,j)
          map(3,j,i) = mapline (3,j)
          mask(j,i)  = 0
        enddo
      enddo
      close (unit=lun)
      return
c
      entry map_mark(long, lat, itype)
      ix = (nx/2) + int((long/180.0D0) * dble(nx/2))
      iy = (ny/2) - int((lat / 90.0D0) * dble(ny/2))
      ix = min(max(ix,1),nx)
      iy = min(max(iy,1),ny)
      if (itype.eq.0) then
        map(1,ix,iy) = -1
        map(2,ix,iy) =  0
        map(3,ix,iy) = -1
        mask(ix,iy)  =  1
      else if (itype.eq.1.and.mask(ix,iy).eq.0) then
        do j = 1,3
          value = map(j,ix,iy)
          if (value.lt.0) value = value + 256
          value = int(0.7 * float(value))
          if (value.gt.127) value = value - 256
          map(j,ix,iy) = value
        enddo
        mask(ix,iy) = 1
      else if (itype.eq.2) then
        map(1,ix,iy) = -1
        map(2,ix,iy) =  127
        map(3,ix,iy) =  127
        mask(ix,iy)  =  1
      endif
      return
c
      entry map_write()
      open (unit=lun, file=cmap//'.rgb', access='direct', form='unformatted', recl=3*nx)
      do i = 1, ny
        do j = 1, nx
          mapline(1,j) = map(1,j,i)
          mapline(2,j) = map(2,j,i)
          mapline(3,j) = map(3,j,i)
        enddo
        write (lun,rec=i) mapline
      enddo
      close (unit=lun)
      call free_lun(lun)
      write (*,'(a)') ' convert -depth 8 -size '//csize//' rgb:'//cmap//'.rgb '//cmap//'.png'
      call system('convert -depth 8 -size '//csize//' rgb:'//cmap//'.rgb '//cmap//'.png')
      return
      end

      subroutine map_stereo_init()
      implicit none
c--
      integer*4 nxs, nys, nazel
      parameter (nxs=5000, nys=5000, nazel=1400)
c--
      integer*1 mapsline (3,nxs), maps(3,nxs,nys), masks(nxs, nys)
      integer*4 lun, i, j, ix, iy, itype, value, ixmin, iymin, ixmax, iymax, nmasked, nxsat, nysat, nzoom, nxz
      integer*4 ixminz, iyminz, ixmaxz, iymaxz, nx_map, ny_map, nc, nx_in, ny_in, ix1, iy1, ix2, iy2, ixdir, iydir
      integer*4 ipxoff
      integer*2 stereo_image(nc,nx_in,ny_in), pix(3), ipx(nxs,nys)
      real*4    y_delta, x_delta, longs, lats, xs, ys
      real*8    long, lat, torad, todeg, pi, pid4, long1, lat1, long2, lat2, rzoom, az, el
      character*1 csq, cdq
      character*6 cmap
      character*9 csize
      character*19 cextract
      character*500 cstring
      character*(*) command
c
      save
c
      cmap = 'satmap' 
      write (csize,'(I4.4,''x'',I4.4)') nxs, nys
c
      torad   = 2.0D0 * dasin(1.0D0) / 180.0
      todeg   = 180.0D0 / 2.0D0 / dasin(1.0D0)
      pi      = 2.0D0 * dasin(1.0D0)
      pid4    = pi / 4.0D0
c
      ixmin   = nxs + 1
      ixmax   = -1
      iymin   = nys + 1
      iymax   = -1
      nmasked = 0
c
      call get_lun(lun)
      open (unit=lun, file='polster.rgb', access='direct', form='unformatted', recl=3*nxs)
      do i = 1, nys
        read (lun,rec=i) mapsline
        do j = 1, nxs
          maps(1,j,i) = mapsline (1,j)
          maps(2,j,i) = mapsline (2,j)
          maps(3,j,i) = mapsline (3,j)
          masks(j,i)  = 0
          ipx(j,i)    = 0
        enddo
      enddo
      close (unit=lun)
c-- Az-El test
      do i = 1, nazel/2 - 1
        j  = (nazel/2) - nint(float(i) / tan(asin(2.0*float(i)/float(nazel))))
        ix = j + 1
        iy = (nazel/2) - i
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
        iy = (nazel/2) + i
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
        j  = (nazel/2) + nint(float(i) / tan(asin(2.0*float(i)/float(nazel))))
        ix = j + 1
        iy = (nazel/2) - i
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
        iy = (nazel/2) + i
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
      enddo
      do i = 1, nazel
        iy = i
        ix = nazel / 2
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
      enddo
      do i = 1, nazel
        ix = i
        iy = nazel / 2
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
      enddo
      do i = 1, 720
        ix = (nazel/2) + int(dsin(dble(i/2)*torad)*dble(nazel/6))
c iy is counterintuitive as the image is vertically flipped for display !
        iy = (nazel/2) + int(dcos(dble(i/2)*torad)*dble(nazel/6))
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
      enddo
      do i = 1, 720
        ix = (nazel/2) + int(dsin(dble(i/2)*torad)*dble(nazel/3))
c iy is counterintuitive as the image is vertically flipped for display !
        iy = (nazel/2) + int(dcos(dble(i/2)*torad)*dble(nazel/3))
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
      enddo
      return
c
      entry map_stereo_mark(long, lat, itype)
      longs = long * torad
      lats  =  lat * torad
      xs =  2.0 * float(nxs)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
      ys = -2.0 * float(nxs)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
      ix = min(max(int(xs)+(nxs/2),1),nxs)
      iy = min(max(int(ys)+(nys/2),1),nys)
      if (ix.eq.1.or.iy.eq.1.or.ix.eq.nxs.or.iy.eq.nxs) then
c        write (*,*) ix, iy, longs*todeg, lats*todeg
c        read (*,*)
        return
      endif
      if (itype.eq.0) then
        maps(1,ix,iy) = -1
        maps(2,ix,iy) =  0
        maps(3,ix,iy) = -1
        masks(ix,iy)  =  1
      else if (itype.eq.1.and.masks(ix,iy).eq.0) then
        do j = 1,3
          value = maps(j,ix,iy)
          if (value.lt.0) value = value + 256
          value = int(0.7 * float(value))
          if (value.gt.127) value = value - 256
          maps(j,ix,iy) = value
        enddo
        masks(ix,iy) = 1
c-- map out the square surrounding the greyed out area
        if (ix.lt.ixmin) ixmin   = ix
        if (iy.lt.iymin) iymin   = iy
        if (ix.gt.ixmax) ixmax   = ix
        if (iy.gt.iymax) iymax   = iy
c-- End of square mapping
      else if (itype.eq.2) then
        maps(1,ix,iy) = -1
        maps(2,ix,iy) =  127
        maps(3,ix,iy) =  127
        masks(ix,iy)  =  1
      endif
      return
c
      entry map_stereo_place(long, lat, pix)
      longs = long * torad
      lats  =  lat * torad
      xs =  2.0 * float(nxs)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
      ys = -2.0 * float(nxs)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
      ix = min(max(int(xs)+(nxs/2),1),nxs)
      iy = min(max(int(ys)+(nys/2),1),nys)
      maps(1,ix,iy) = (pix(1) / 8) 
      maps(2,ix,iy) = (pix(2) / 8)
      maps(3,ix,iy) = (pix(3) / 8)
      return
c
      entry map_stereo_place_azel(az, el)
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
      entry map_stereo_zoom_init(nxsat, nysat, nx_map, ny_map)
      do j = 1, nys
        do i = 1,nxs
          if (masks(i,j).eq.1) nmasked = nmasked + 1
        enddo
      enddo
c-- changed to the rzoom
      nzoom = nint(sqrt(float(nxsat*nysat)/float(nmasked)))
      rzoom =      sqrt(float(nxsat*nysat)/float(nmasked))
      write (*,*) ixmin, ixmax, iymin, iymax, nmasked, nxsat*nysat, int(float(nxsat*nysat)/float(nmasked)), nzoom, rzoom
c      ixminz = ((nzoom*nxs)/2) - nzoom*((nxs/2)-ixmin)
c      iyminz = ((nzoom*nys)/2) - nzoom*((nys/2)-iymin)
c      ixmaxz = ((nzoom*nxs)/2) - nzoom*((nxs/2)-ixmax)
c      iymaxz = ((nzoom*nys)/2) - nzoom*((nys/2)-iymax)
c      nxz    = nzoom * nxs
      ixminz = int((rzoom*float(nxs))/2.0) - int(rzoom*float((nxs/2)-ixmin))
      iyminz = int((rzoom*float(nys))/2.0) - int(rzoom*float((nys/2)-iymin))
      ixmaxz = int((rzoom*float(nxs))/2.0) - int(rzoom*float((nxs/2)-ixmax))
      iymaxz = int((rzoom*float(nys))/2.0) - int(rzoom*float((nys/2)-iymax))
      nxz    = int(rzoom * float(nxs))
c
      nx_map = ixmaxz - ixminz + 1
      ny_map = iymaxz - iyminz + 1
      write (*,*) ixminz, ixmaxz, iyminz, iymaxz, nxz, nx_map*ny_map
      return
c
      entry map_stereo_zoom_mark(stereo_image, nc, nx_in, ny_in,long1, lat1, long2, lat2)
      longs = long1 * torad
      lats  =  lat1 * torad
      xs =  2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
      ys = -2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
      ix1 = min(max(int(xs)+(nxz/2),1),nxz)
      iy1 = min(max(int(ys)+(nxz/2),1),nxz)
      ix1 = ix1 - ixminz + 1
      iy1 = iy1 - iyminz + 1
      longs = long2 * torad
      lats  =  lat2 * torad
      xs =  2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
      ys = -2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
      ix2 = min(max(int(xs)+(nxz/2),1),nxz)
      iy2 = min(max(int(ys)+(nxz/2),1),nxz)
      ix2 = ix2 - ixminz + 1
      iy2 = iy2 - iyminz + 1
c-- Single Point
      if (ix1.eq.ix2.and.iy1.eq.iy2) then
        ix = ix1
        iy = iy1
        if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
          stereo_image(1,ix,iy) = 1023
          stereo_image(2,ix,iy) = 1023
          stereo_image(3,ix,iy) = 0
        endif
        goto 1
      endif
c-- Vertical  Line
      if (ix1.eq.ix2) then
        ix = ix1
        iydir = 1
        if (iy2.lt.iy1) iydir = -1
        do iy = iy1, iy2, iydir
          if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
            stereo_image(1,ix,iy) = 1023
            stereo_image(2,ix,iy) = 1023
            stereo_image(3,ix,iy) = 0
          endif
        enddo
        goto 1
      endif
c-- Horizontal  Line
      if (iy1.eq.iy2) then
        iy = iy1
        ixdir = 1
        if (ix2.lt.ix1) ixdir = -1
        do ix = ix1, ix2, ixdir
          if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
            stereo_image(1,ix,iy) = 1023
            stereo_image(2,ix,iy) = 1023
            stereo_image(3,ix,iy) = 0
          endif
        enddo
        goto 1
      endif
c-- The Rest
      if (abs(ix2-ix1).gt.abs(iy2-iy1)) then
        ixdir = 1
        if (ix2.lt.ix1) ixdir = -1
        do ix = ix1, ix2, ixdir
          iy = nint(float(iy1) + ((float(ix-ix1)/float(ix2-ix1))*float(iy2-iy1)))
          if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
            stereo_image(1,ix,iy) = 1023
            stereo_image(2,ix,iy) = 1023
            stereo_image(3,ix,iy) = 0
          endif
        enddo
        goto 1
      else
        iydir = 1
        if (iy2.lt.iy1) iydir = -1
        do iy = iy1, iy2, iydir
          ix = nint(float(ix1) + ((float(iy-iy1)/float(iy2-iy1))*float(ix2-ix1)))
          if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
            stereo_image(1,ix,iy) = 1023
            stereo_image(2,ix,iy) = 1023
            stereo_image(3,ix,iy) = 0
          endif
        enddo
        goto 1
      endif
 1    continue
      return
c
      entry map_stereo_zoom_place(stereo_image,nc,nx_in,ny_in,long,lat,pix,ipxoff)
      longs = long * torad
      lats  =  lat * torad
      xs =  2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
      ys = -2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
      ix = min(max(int(xs)+(nxz/2),1),nxz)
      iy = min(max(int(ys)+(nxz/2),1),nxz)
      ix = ix - ixminz + 1
      iy = iy - iyminz + 1
      if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
        if (ipx(ix,iy).eq.0) then
          stereo_image(1,ix,iy) = pix(1)
          stereo_image(2,ix,iy) = pix(2)
          stereo_image(3,ix,iy) = pix(3)
          ipx(ix,iy)            = ipxoff
        else if (ipx(ix,iy).ne.0.and.ipxoff.le.ipx(ix,iy)) then
          stereo_image(1,ix,iy) = pix(1)
          stereo_image(2,ix,iy) = pix(2)
          stereo_image(3,ix,iy) = pix(3)
          ipx(ix,iy)            = ipxoff
        endif
      endif
      return
c
      entry map_stereo_write(command)
      open (unit=lun, file=cmap//'-stereo.rgb', access='direct', form='unformatted', recl=3*nxs)
      do i = 1, nys
        do j = 1, nxs
          mapsline(1,j) = maps(1,j,i)
          mapsline(2,j) = maps(2,j,i)
          mapsline(3,j) = maps(3,j,i)
        enddo
        write (lun,rec=i) mapsline
      enddo
      close (unit=lun)
      call free_lun(lun)
      csq = char(ichar("'"))
      cdq = char(ichar('"'))
      write (*,'(a)') ' convert -depth 8 -size '//csize//' rgb:'//cmap//'-stereo.rgb -flip '//' '//cmap//'-stereo.png'
      call system('convert -depth 8 -size '//csize//' rgb:'//cmap//'-stereo.rgb -flip '//' '//cmap//'-stereo.png')
c-- Add command and credit lines
      cstring = '-pointsize 20 -fill red -draw '//csq//'text 70,50 '//cdq//command(1:lnblnk(command))//cdq//csq
      write (*,'(a)') ' convert  '//cmap//'-stereo.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo.png'
      call system('convert  '//cmap//'-stereo.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo.png')
      cstring = '-pointsize 20 -fill red -draw '//csq//'text 70,80 '//cdq//'Image credit: Retro StÃ¶ckli, NASA Earth Observatory'//cdq//csq
      write (*,'(a)') ' convert  '//cmap//'-stereo.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo.png'
      call system('convert  '//cmap//'-stereo.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo.png')
c-- Zoom on shaded area (for twitter size image)
      write (cextract,'(I4.4,''x'',I4.4,''+'',I4.4,''+'',I4.4)') ixmax-ixmin+1, iymax-iymin+1, ixmin, iymin
      write (*,'(a)') ' convert -depth 8 -size '//csize//' -extract '//cextract//' rgb:'//cmap//'-stereo.rgb -flip '//' '//cmap//'-stereo-zoom.png'
      call system('convert -depth 8 -size '//csize//' -extract '//cextract//' rgb:'//cmap//'-stereo.rgb -flip '//' '//cmap//'-stereo-zoom.png')
c-- Add command and credit lines
      cstring = '-pointsize 10 -fill red -draw '//csq//'text 35,25 '//cdq//command(1:lnblnk(command))//cdq//csq
      write (*,'(a)') ' convert  '//cmap//'-stereo-zoom.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo-zoom.png'
      call system('convert  '//cmap//'-stereo-zoom.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo-zoom.png')
      cstring = '-pointsize 10 -fill red -draw '//csq//'text 35,40 '//cdq//'Image credit: Retro StÃ¶ckli, NASA Earth Observatory'//cdq//csq
      write (*,'(a)') ' convert  '//cmap//'-stereo-zoom.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo-zoom.png'
      call system('convert  '//cmap//'-stereo-zoom.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo-zoom.png')
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
      character*9 csize
      character*19 cextract
      character*500 cstring
      character*(*) command
c
      save
c
      cmap    = 'satmap' 
      torad   = 2.0D0 * dasin(1.0D0) / 180.0
      todeg   = 180.0D0 / 2.0D0 / dasin(1.0D0)
      write (csize,'(I4.4,''x'',I4.4)') nazel, nazel
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
      write (*,'(a)') ' convert -depth 8 -size '//csize//' rgb:'//cmap//'-azel.rgb -flip '//' '//cmap//'-azel.png'
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
      
      subroutine get_borders(stereo_image,nc,nx_in,ny_in)
      implicit none
      integer*4 nx_in, ny_in, nc
      integer*2 stereo_image(nc,nx_in,ny_in)
c
      integer*1 buf(1048576)
      integer*4 n_actual, i, j, j1, j2, k, parts(500)
c
      integer*1 i1_shape(4), i1_xmin(8), i1_ymin(8), i1_xmax(8), i1_ymax(8), i1_recnum(4), i1_length(4)
      integer*1 i1_numparts(4), i1_numpoints(4), i1_fl(4), i1_parts(4), i1_x(8), i1_y(8)
      integer*4 i4_shape, i4_recnum, i4_length, i4_numparts, i4_numpoints, i4_fl, i4_parts
      real*8    xmin, ymin, xmax, ymax, i8_x, i8_y, i8_xsave, i8_ysave
c
      equivalence (i1_shape, i4_shape), (i1_xmin,xmin), (i1_ymin, ymin), (i1_xmax, xmax), (i1_ymax, ymax), (i1_fl, i4_fl)
      equivalence (i1_recnum, i4_recnum), (i1_length, i4_length), (i1_numparts, i4_numparts), (i1_numpoints, i4_numpoints)
      equivalence (i1_parts, i4_parts), (i1_x, i8_x), (i1_y,i8_y)
c
      n_actual = -1
      call get_shp_borders('borders.shp',buf,1048576,n_actual,100)
      call move_bytes(buf(25),i1_fl(1),4)
      call move_bytes(buf(33),i1_shape(1),4)
      do while (.true.)
        call get_shp_borders('borders.shp',buf,1048576,n_actual,52)
        if (n_actual.eq.-1) goto 1
        call move_bytes_reverse(buf(1), i1_recnum(1),4)
        call move_bytes_reverse(buf(5), i1_length(1),4)
        call move_bytes        (buf(9), i1_shape(1),4)
        call move_bytes        (buf(13),i1_xmin(1),8)
        call move_bytes        (buf(21),i1_ymin(1),8)
        call move_bytes        (buf(29),i1_xmax(1),8)
        call move_bytes        (buf(37),i1_ymax(1),8)
        call move_bytes        (buf(45),i1_numparts(1),4)
        call move_bytes        (buf(49),i1_numpoints(1),4)
        call get_shp_borders('borders.shp',buf,1048576,n_actual,i4_length*2+8-52)
        if (n_actual.eq.-1) goto 1
        do i = 1, i4_numparts
          call move_bytes(buf((i-1)*4+1), i1_parts,4)
          parts(i) = i4_parts
        enddo
        do i = 1, i4_numparts
          j1 = parts(i) + 1
          if (i.lt.i4_numparts) then
            j2 = parts(i+1)
          else
            j2 = i4_numpoints
          endif
          do j = j1, j2
            call move_bytes(buf(i4_numparts*4+(j-1)*16+1), i1_x(1),8)
            call move_bytes(buf(i4_numparts*4+(j-1)*16+9), i1_y(1),8)
            if (j.gt.j1) call map_stereo_zoom_mark(stereo_image, nc, nx_in, ny_in, i8_xsave, i8_ysave, i8_x, i8_y)
            i8_xsave = i8_x
            i8_ysave = i8_y
          enddo
        enddo
 2      continue
        enddo
 1    continue
      call rewind_shp_borders()
      return
      end

      subroutine get_shp_borders(filenm,buf,n_max,n_actual,n_req)
      implicit none
c
      integer*4 n_max, n_actual, n_req
      integer*1 buf(n_max)
      character*(*) filenm
c
      integer*1 buffer(16384)
c
      integer*4 iptr, irec, lun, n_used, i, nbytes, eof
c
      save iptr, buffer, irec, lun, nbytes, eof
c
      n_used = 0
      if (n_actual.eq.-1) then
        call get_lun(lun)
        irec = 0
        open (unit=lun,file=filenm,status='old',form='unformatted',access='direct',recl=16384)
        inquire (file=filenm, size=nbytes)
        iptr = 0
        eof  = 0
        irec = 0
        call get_file_buffer_i1(lun,16384,buffer,irec,iptr)
      endif
      do i = 1, n_req
        iptr = iptr + 1
        if (iptr.gt.16384) then
          call get_file_buffer_i1(lun,16384,buffer,irec,iptr)
          if (irec*16384.gt.nbytes) then
            eof = 1
          endif
          if (irec.eq.-1) then
            close (unit=lun)
            call free_lun(lun)
            n_actual = -1
            return
          endif
        endif
        if (eof.eq.1) then
          if (iptr+(irec-1)*16384.gt.nbytes) then
            eof = 0
            close (unit=lun)
            call free_lun(lun)
            n_actual = -1
            return
          endif
        endif
        n_used = n_used + 1
        buf(n_used) = buffer(iptr)
      enddo
      n_actual = n_used
      return
c
      entry rewind_shp_borders()
      iptr = 0
      eof  = 0
      irec = 0
      if (lun.gt.0) then
        call get_file_buffer_i1(lun,16384,buffer,irec,iptr)
        close (unit=lun)
        call free_lun(lun)
      endif
      return
      end

      subroutine get_file_buffer_i1(lun,nbytes,buffer,irec,iptr)
      implicit none
      integer*4 lun, nbytes, irec, iptr, ios
      integer*1 buffer(nbytes)
c
      irec = irec + 1
      ios = 0
      read (lun,rec=irec,iostat=ios) buffer
      if (ios.ne.0) goto 1
      if (iptr.gt.nbytes) iptr = iptr - nbytes
      return
 1    irec = -1
      return
      end

      subroutine write_buffer(lun,irec,array,n)
      implicit none
      integer*4 n, irec, lun
      integer*2 array(n)
c
      write (lun,rec=irec) array
      return
      end

      subroutine move_bytes_reverse(in,out,n)
      implicit none
      integer*4 n, i
      integer*1 in(n), out(n)
      do i = 1,n
        out(n+1-i) = in(i)
      enddo
      return
      end

      subroutine move_bytes(in,out,n)
      implicit none
      integer*4 n, i
      integer*1 in(n), out(n)
      do i = 1,n
        out(i) = in(i)
      enddo
      return
      end

      subroutine linear_interpolate(x1, x2, y1, y2, yv)
      implicit none
      real*8 x1, x2, y1, y2, yv
c
      yv = y1 + ((dble(int(x2))-x1)/(x2-x1)) * (y2 - y1)
c
      return
      end

c      subroutine wgs84_to_round(long_in, lat_in, long_out, lat_out)
c      implicit none
c      real*8 long_in, lat_in, long_out, lat_out
c
c      real*8 lla(3), xyz(3), torad, todeg, twopi, theta, rearth, rwgs84
c
c      torad   = 2.0D0 * dasin(1.0D0) / 180.0D0
c      todeg   = 180.0D0 / 2.0D0 / dasin(1.0D0)
c      twopi   = 4.0D0 * dasin(1.0D0)
c      lla(1)  = long_in * torad
c      lla(2)  = lat_in * torad
c      lla(3)  = 0.0D0
c      call lla2xyz(lla, xyz)
c      rwgs84  = dsqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2) + xyz(3) * xyz(3))
c      rearth  = 6371000.0D0
c      xyz(1)  = xyz(1) * (rearth / rwgs84)
c      xyz(2)  = xyz(2) * (rearth / rwgs84)
c      xyz(3)  = xyz(3) * (rearth / rwgs84)
c      theta   = datan2(xyz(2),xyz(1))
c      if (theta.lt.0) theta = theta + twopi
c      long_out = theta * todeg
c      lat_out  = datan(xyz(3)/dsqrt(xyz(1)* xyz(1) + xyz(2) * xyz(2))) * todeg
c
c      return
c      end

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

      subroutine get_hist_limits(image, nx, nc, ny, ilim_low, ilim_high, nch, nyl)
      implicit none
      integer*4 nx, nc, ny, nch, nyl, ilim_low(nc), ilim_high(nc)
      integer*2 image(nx, nc, ny)
c
      integer*4 i, j, k, l, hist(1024), ich
      real*4 sum, sumcheck
c
      do j = 1, nch
        do i = 1, 1024
          hist(i) = 0
        enddo
        ilim_low(j)  = 0
        ilim_high(j) = 0
        do k = 1, nyl
          do i = 1, nx
            ich = image(i, j, k)
            ich = min(max(ich + 1, 1), 1024)
            if (ich.gt.1) hist(ich) = hist(ich) + 1
          enddo
        enddo
        sum = 0.
        do i = 1, 1024
          sum = sum + float(hist(i))
        enddo
        sumcheck = 0.
        do i = 1, 1024
          sumcheck = sumcheck + float(hist(i))
          if (sumcheck.gt.0.005*sum.and.ilim_low (j).eq.0) ilim_low(j)  = i
          if (sumcheck.gt.0.995*sum.and.ilim_high(j).eq.0) ilim_high(j) = i - 1
        enddo
      enddo
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

      subroutine hist_plot(hist,nx,ny,filenm)
      implicit none
      integer*4 nx, ny, hist(nx, ny)
      character*(*) filenm
c
      integer*4 lunplot, i, j
      real*4 xr(2), yr(2), rhisty(1024), rhistx(1024)
      character*79 ctxt
      call get_lun(lunplot)
      call PS_init_colourtable(1,'color')
      call set_PS_fullpage()
      call pkg_openpl(filenm(1:lnblnk(filenm))//'.ps', lunplot)
      call newpen(6)
      xr(1) = 0.
      xr(2) = 1024.
      yr(1) = hist(1,1)
      yr(2) = yr(1)
      do i = 1,nx
        do j = 1,5
          if (hist(i,j).lt.yr(1)) yr(1) = hist(i,j)
          if (hist(i,j).gt.yr(2)) yr(2) = hist(i,j)
        enddo
      enddo
      call pkg_frame(11,-6,1.,xr,yr,'Bin #','Frequency','FY Channel histogram')
      do j = 1,5
        do i = 1,nx
          rhisty(i) = float(hist(i,j))
          rhistx(i) = float(i)
        enddo
        call newpen(7-j)
        call pkg_plhist(rhistx,rhisty,nx)
      enddo
      do j = 1,5
        write (ctxt,'(''Channel '',i1)') j
        call newpen(7-j)
        call pkg_pltextbl(ctxt(1:lnblnk(ctxt)),j)
      enddo
      call pkg_clospl()
      call free_lun(lunplot)
      return
      end


      subroutine wltoRGB(wl,RGB,rnorm)
      implicit none
c
      real*4    wl, rnorm
      integer*4 RGB(3)
c
      real*4    gamma, r, g, b, sss
      integer*4 ir, ig, ib
c
      gamma =   0.80
c
c-- Remember 380 < wl < 780 nm (!!)
c
      IF ((WL.GE.380.).AND.(WL.LE.440.)) THEN 
        R = -1.*(WL-440.)/(440.-380.)
        G = 0.
        B = 1.
      ENDIF
      IF ((WL.GE.440.).AND.(WL.LE.490.)) THEN
        R = 0.
        G = (WL-440.)/(490.-440.)
        B = 1.
      ENDIF
      IF ((WL.GE.490.).AND.(WL.LE.510.)) THEN 
        R = 0.
        G = 1.
        B = -1.*(WL-510.)/(510.-490.)
      ENDIF
      IF ((WL.GE.510.).AND.(WL.LE.580.)) THEN 
        R = (WL-510.)/(580.-510.)
        G = 1.
        B = 0.
      ENDIF
      IF ((WL.GE.580.).AND.(WL.LE.645.)) THEN
        R = 1.
        G = -1.*(WL-645.)/(645.-580.)
        B = 0.
      ENDIF
      IF ((WL.GE.645.).AND.(WL.LE.780.)) THEN
        R = 1.
        G = 0.
        B = 0.
      ENDIF
c
c      LET THE INTENSITY SSS FALL OFF NEAR THE VISION LIMITS
c
      IF (WL.GT.700.) THEN
         SSS=.3+.7* (780.-WL)/(780.-700.)
      ELSE IF (WL.LT.420.) THEN
         SSS=.3+.7*(WL-380.)/(420.-380.)
      ELSE
         SSS=1.
      ENDIF
c
c      GAMMA ADJUST
c
      IR = int(rnorm*((SSS*R)**GAMMA))
      IG = int(rnorm*((SSS*G)**GAMMA))
      IB = int(rnorm*((SSS*B)**GAMMA))
c
      RGB(1) = IR
      RGB(2) = IG
      RGB(3) = IB
      return
      end



c------------- Plotsoft ---------------
c**
c**
c**	pkg_common:
c**	-----------
c**	blockdata routine to fill the common blocks used by this package
c**	with the required data.
c**	
c**
c**	/frparm/
c**	this common block defines the coordinates of special positions
c**	for the frame as well as some special sizes.
c**
c**
c**	description of the variables in /frparm/:
c**	   pwidth - width of "paper" (max. y-size)
c**	   yfrac  - current factor to scale y-axis; 0<yfrac<1
c**	   yshft  - current up-shift of y-axis
c**	   x0     - x-coordinate of left bottom corner of frame
c**	   y0     - y-coordinate of left bottom corner of frame
c**	   yaxis  - length of y-axis
c**	   grdsz  - length of a grid mark
c**	   symsz  - size of symbols to be plotted
c**	   dvx    - distance between x-axis and its value annotation
c**	   dtx    - distance between x-axis and its title
c**	   yhd1   - distance between y-axis and first text line above plot
c**	   yhd2   - distance between y-axis and second text line above plot
c**	   dvy    - distance between y-axis and its value annotation
c**	   dty    - distance between y-axis and its title
c**	   chsz   - character size (one character is a chsz by chsz matrix)
c**	   dl	  - length of line segments for dashed plots (in mm)
c**	   ndgt   - number of significant digits available for value annotation
c**
c**
c**	/pkgpar/
c**	this common defines some overall parameters
c**
c**	description of the variables in pkgpar
c**	   swopen - logical specifying whether the plotfile is open or not
c**	   swfrst - logical specifying whether the current plot is the first
c**	            plot in the file or not
c**	   swxlog - logical specifying whether the x-axis is log or not
c**	   swylog - logical specifying whether the y-axis is log or not
c**	   xmin   - minimum x-value that can be plotted
c**	   xmax   - maximum x-value that can be plotted
c**	   ymin   - minimum y-value that can be plotted
c**	   ymax   - maximum y-value that can be plotted
c**	   xfact  - length of x-axis divided by length of y-axis
c**	   idset  - pkg_preset / pkg_reset id.
c**	   redy   - saved reduction factor for y_axis (pkg_preset argument)
c**
c**
	block data pkg_common
c**
c**
	logical swopen,swfrst,swxlog,swylog
	common /pkgpar/ swopen,swfrst,swxlog,swylog,xmin,xmax,ymin,ymax,
     *	                xfact,idset,redy
	common /frparm/ pwidth,yfrac,yshft,x0,y0,yaxis,grdsz,symsz,dvx,
     *	     	        dtx,yhd1,yhd2,dvy,dty,chsz,dl,ndgt
c**
	data pwidth,yfrac,yshft, x0 , y0 ,yaxis,grdsz,symsz,dvx,dtx/
     *	     254. ,  1. ,  0. , 35., 28., 200.,  2. ,  2. , 2.,14./,
     *	     yhd1,yhd2,dvy,dty,chsz,  dl,ndgt/
     *	     20., 5. , 2.,25., 8. , 1.5,   5 /
	data swopen,swfrst,idset/
     *	    .false.,.true.,  0  /
c**
c**
	end
c**
c**	pkg_openpl -- routine to open the plotfile and to do the initialization
c**	              this routine must be called before any other routine
c**	              is called
c**	pkg_clospl -- routine to close the plotfile
c**	              this routine must be called before the program exits
c**	pkg_dfyrng -- routine to set the y-range to its maximum
c**	pkg_gdyrng -- routine to set the y-range so, that a square frame fits
c**	              on a single gould sheet
c**	pkg_preset -- routine to preset the plotting parameters for stacked
c**		      plots
c**	pkg_reset  -- routine to reset the plotting parameters for stacked plots
c**
c**	usage:
c**	   call pkg_openpl(fname,lun)
c**	   call pkg_clospl
c**	   call pkg_dfyrng
c**	   call pkg_gdyrng
c**	   call pkg_preset(fractx,fracty,rat)
c**	   call pkg_reset(fractx,fracty,rat)
c**
c**	parameters:
c**	   fname  - character variable defining the plot filename
c**	   lun    - lun to be used for accessing the plotfile
c**	   fractx - reduction factor for x-axis
c**	   fracty - reduction factor for y-axis
c**	   rat    - (modified) x/y ratio
c**
c**
	subroutine pkg_openpl(fname,lun)
c**
c**
	parameter     (nvar=14)
	character*(*) fname
	logical       swopen,swfrst,swxlog,swylog
	common /frparm/ pwidth,yfrac,yshft,parm(nvar)
	common /pkgpar/ swopen,swfrst,swxlog,swylog,xmin,xmax,ymin,ymax,
     *	                xfact,idset,redy
	logical switch
c**
c**
c**
	switch = .true.
	if(swopen) goto 5
2	continue
	if (lun.gt.0) call plots(fname,lun)
	swopen = .true.
	swfrst = .true.
c**
c**	scale frame coordinates between 0 and 1
c**
	do i=1,nvar-1
	  parm(i) = parm(i) / pwidth
	enddo
	call plot_rotate_none()
	goto 99
c******************************************************************************
	entry pkg_clospl
c**
	switch = .true.
	if(.not.swopen) goto 99
	switch = .false.
5	call plot(0.,0.,999)
c	call route(0,0,1,ier)
	swopen = .false.
	goto 10
c**
c**	rescale frame coordinates to the original values
c**
6	do i=1,nvar-1
	  parm(i) = parm(i) * pwidth
	enddo
	if(switch) goto 2
	goto 99
c******************************************************************************
	entry pkg_dfyrng
c**
10	if(idset.eq.-1) then
c**
c**	  reset for stacked plots
c**
	  fac1 = exp(redy - 1.)
	  fac2 = fac1 / redy
	  chsz = chsz / fac2
	  yhd1 = yhd1 / fac2
	  yhd2 = yhd2 / fac2
	  dtx  = dtx  / fac2
	  dty  = dty  / fac2
	  x0   = x0   / fac2
	  y0   = y0   / fac2
	  dl   = dl * redy
	  call factor(1.)
	  idset = 0
	end if
	do i=1,nvar-1
	  parm(i) = parm(i) / yfrac
	enddo
	parm(2) = parm(2) - yshft
	yshft = 0.
	yfrac = 1.
	if(idset.eq.+1) then
c**
c**	  preset for stacked plots
c**
	  fac1 = exp(redy - 1.)
	  fac2 = fac1 / redy
	  chsz = fac2 * chsz
	  yhd1 = fac2 * yhd1
	  yhd2 = fac2 * yhd2
	  dtx  = fac2 * dtx
	  dty  = fac2 * dty
	  x0   = fac2 * x0
	  y0   = fac2 * y0
	  dl   = dl / redy
	  call factor(redy)
	  idset = 0
	end if
	if(.not.swopen) goto 6
	goto 99
c******************************************************************************
	entry pkg_gdyrng
c**
c	yfrac = 208. / pwidth	
	yfrac = 0.95	
c !plot must fit on one fan-fold sheet
	yshft = 0.5*parm(2)		
c !room for additional text under plot
	parm(2) = parm(2) + yshft	
c !new y-origin (=3*y0)
	do i=1,nvar-1
c !scale frame coordinates
	  parm(i) = parm(i) * yfrac
	enddo
	goto 99
c******************************************************************************
	entry pkg_preset(fractx,fracty,rat)
c**
	redy = fracty
	rat  = rat  * fractx / fracty
	idset = +1
	goto 99
c******************************************************************************
	entry pkg_reset(fractx,fracty,rat)
c**
	redy = fracty
	rat  = rat  * fracty / fractx
	idset = -1
	goto 99
c**
c**
99	continue
	return
	end
c**
c**	pkg_pldatr -- routine to plot real data points in a frame previously
c**	              plotted by routine pkg_frame.
c**
c**	pkg_pldati -- same as pkg_pldatr, but data points are integers.
c**
c**	pkg_plstr  -- routine to plot a character string in a frame previously
c**	              plotted by pkg_frame.
c**
c**	pkg_plstra -- routine to plot a character string under a given angle 
c**                   in a frame previously plotted by pkg_frame.
c**
c**
c**	usage:
c**	   call pkg_pldatr(itp,x,y,np)
c**	   call pkg_pldati(itp,ix,iy,np)
c**	   call pkg_plstr (x,y,ctxt)
c**        call pkg_plstra(x,y,ctxt,angt,szud)
c**
c**
c**	parameters:
c**	   itp - the type of plot to be made:
c**	         if itp =-99 no boundary checking is made
c**                     < 0, a line will be plotted between data points
c**	                >=0, a symbol will be plotted at the specified
c**	                     coordinates; itp represents the number of
c**	                     the symbol as defined in the symbol table
c**	   x   - real array with x-values
c**	   ix  - integer array with x-values
c**	   y   - real array with y-values
c**	   iy  - integer array with y-values
c**	   np  - number of elements in x- and y-arrays
c**	   ctxt- character variable with string to be plotted
c**        angt- acw angle in degrees about +x for plotting string 
c**	   szud- user defined string character size in mm (0 = default)
c**
c**
c**	note:
c**	1) if the plotfile is not yet opened, or if frame has not yet
c**	   been called, the call is merely a nop.
c**
c**	2) when plotting vectors and (normally) symbols the pen is 
c**	   raised outside the edge of the (niceed) x, y limits defined 
c**	   in the pkg_frame call.
c**
c**	2) if np < 0 then dashed lines are plotted.
c**
c**
	subroutine pkg_pldatr(itp,x,y,np)
c**
	parameter (f0=0.3, f1=1.-f0)
c**
	character*(*) ctxt
	logical       swopen,swfrst,swxlog,swylog,swreal,swstr,swfull
	logical       swbnd
	dimension     x(*),ix(*),y(*),iy(*)
c**
	common /frparm/ pwidth,yfrac,yshft,x0,y0,yaxis,grdsz,symsz,dvx,
     *	                dtx,yhd1,yhd2,dvy,dty,chsz,dl,ndgt
	common /pkgpar/ swopen,swfrst,swxlog,swylog,xmin,xmax,ymin,ymax,
     *	                xfact,idset,redy
c**
c**
c**
	swreal = .true.
	goto 1
c**
c**
c***********************************************************************
c**
	entry pkg_pldati(itp,ix,iy,np)
c**
c**
	swreal = .false.
1	swstr = .false.
	ityp = itp
	npnt = abs(np)
	swfull = np.gt.0
	swbnd = itp.ne.-99
	goto 3
c**
c**
c***********************************************************************
c**
	entry pkg_plstr(x,y,ctxt)
c**
c**
	swstr = .true.
	swreal = .true.
	ang = 0.
	csize = chsz/2.
	ityp = 999
	npnt = 1
	goto 3
c**
c**
c***********************************************************************
c**
	entry pkg_plstra(x,y,ctxt,angt,szud)
c**
c**
	swstr = .true.
	swreal = .true.
	csize = chsz/2.
	if (szud.gt.0.) csize = szud/pwidth
	ityp = 999
	npnt = 1
	ang = angt
c
c***********************************************************************
c
3	if(.not.swopen) return		
c !plotfile not open
	if(swfrst) return		
c !frame was not yet called
c**
c**	calculate x scale factor
c**
	if(.not.swxlog) then
	  xf = xmax - xmin
	else
	  xf = log10(xmax/xmin)
	endif
c**
c**	calculate y scale factor
c**
	if(.not.swylog) then
	  yf = ymax - ymin
	else
	  yf = log10(ymax/ymin)
	endif
	xf = (yaxis*xfact) / xf
	yf = yaxis / yf
c**
c** set frame boundaries
c**
	xminp = pkg_xcoord (xmin, xf)
	xmaxp = pkg_xcoord (xmax, xf)
	yminp = pkg_ycoord (ymin, yf)
	ymaxp = pkg_ycoord (ymax, yf)
	if (swbnd) then
	  call bound (xminp, xmaxp, yminp, ymaxp)       
c !limit vectors
	else
	  call bound (-1.e6,  1.e6, -1.e6,  1.e6)	
c !delimit vectors
	end if
c**
c***********************************************************************
c**
c**	translate x,y values to x,y coordinates and plot
c**
	ipen = 13
c
	do i=1,npnt
c
	  if(swreal) then
	    xp = x(i)
	    yp = y(i)
	  else
	    xp = ix(i)
	    yp = iy(i)
	  endif
	  xp = pkg_xcoord(xp,xf)
	  yp = pkg_ycoord(yp,yf)
	  if(ityp.ge.0) then
	    if(swstr) then
	      call symbol(xp,yp,csize,ctxt,ang,len(ctxt))
	    else 
	      if ((xminp.le.xp.and.xp.le.xmaxp .and.
     *	           yminp.le.yp.and.yp.le.ymaxp).or.np.lt.0)
     *	             call symbol(xp,yp,symsz,char(ityp),ang,-1)
	    endif
	  else
	    if (swfull) then
	      call plot(xp,yp,ipen)	
c !full line plot
	      ipen = 12
c
	    else if (i.eq.1) then
	      call plot(xp,yp,13)		
c !dashed line plot - position pen
	      xpp = xp
	      ypp = yp
	    else if (i.gt.1) then
	      dx = xp - xpp
	      dy = yp - ypp
	      dr = sqrt (dx**2 + dy**2)
	      nd = dr/dl + f0	
c !number of line segments of length dl
	      dll = dr/(nd + f0)
	      nd = nd + 1
	      if (dr.lt.1.e-30) dr = 1.e-30
	      dxl = f0*dll*dx/dr		
c !x-line increment
	      dxg = f1*dll*dx/dr		
c !x-gap  increment
	      dyl = f0*dll*dy/dr
	      dyg = f1*dll*dy/dr
	      do id = 1, nd
	        if (id.lt.nd) then
	          xpp = xpp + dxl
	          ypp = ypp + dyl
	        else
	          xpp = xp	
c !last segment - carry through to x,yp
	          ypp = yp
	        endif
	        call plot(xpp,ypp,12)	
c !draw line segment
	        if (id.lt.nd) then
	          xpp = xpp + dxg
	          ypp = ypp + dyg
	          call plot(xpp,ypp,13)	
c !"draw" gap segment
	        endif
	      enddo
	    endif
	  endif
	enddo	
	return
	end

c**	pkg_frame -- routine to define and plot a frame for subsequent plotting
c**	             of data
c**
c**
c**	calling sequence:
c**	   call pkg_frame(intx,inty,rat,xr,yr,cxt,cyt,ctit)
c**
c**	where:
c**	   intx - if <0: number of decades (log), otherwise number of
c**	          intervals (lin) for x-axis
c**	   inty - same for y-axis
c**	   rat  - length of x-axis divided by length of y-axis
c**	          if rat=1, the square plot will be scaled so, that it fits
c**	                    on a single gould sheet
c**               if rat=0, no annotation and grid marks are plotted
c**               if rat < 0. then you can plot on the same page using
c**               pkg_preset and pkg_reset
c**	   xr   - 2 element real array specifying the minimum and maximum
c**	          x value that can be plotted in the frame
c**	   yr   - 2 element real array specifying the minimum and maximum
c**	          y value that can be plotted in the frame
c**	   cxt  - character variable with text to be written under x-axis
c**	   cyt  - character variable with text to be written along y-axis
c**	   ctit - character variable with title text to be written above frame
c**
c**
c**	note:
c**	   if the plotfile is not yet opened, the call is merely a nop
c**
c**
	subroutine pkg_frame(intx,inty,ratt,xr,yr,cxt,cyt,ctit)
c**
	parameter     (ncdat=15)
	parameter     (xycor=1.01146)
	character*(*) cxt,cyt,ctit
	character*5   czone
	character*15  cdate
	character*16  ctime
	logical       swopen,swfrst,swxlog,swylog,swxant,swyant
	dimension     xr(2),yr(2)
	integer*4     timval(8)
c**
	common /frparm/ pwidth,yfrac,yshft,x0,y0,yaxis,grdsz,symsz,dvx,
     *	                dtx,yhd1,yhd2,dvy,dty,chsz,dl,ndgt
	common /pkgpar/ swopen,swfrst,swxlog,swylog,xmin,xmax,ymin,ymax,
     *	                xfact,idset,redy
c**
	data cdate/'date: dd-mmm-yy'/
	data ctime/'time: hhmmss.sss'/
c**
c**
c**
	if(.not.swopen) return		
c !plotfile not open
	rat = abs (ratt)		
c !make sure old programs will work !*
c
	if (ratt.ge.0.) then
	  if (.not.swfrst) call plot(0.,0.,999)	
c !end the former plot
	else
	  call plot(-x0,(1.-y0),-3)	
c !shift origin for extra plots
	endif
c
	call pkg_dfyrng			
c !set default y-range
c**
c**	determine frame type and number of gridmarks
c**
	if(rat.eq.1.) call pkg_gdyrng	
c !set to fit plot on a single gould sheet
	idfr = 1
c
	if (intx.eq.0) then
	  nx = 20
	  swxant = .false.
	else
	  nx = intx
	  swxant = .true.
	endif
	ngrx = abs(nx)
	if(nx.lt.0) idfr = 3
c
	if (inty.eq.0) then
	  ny = 20
	  swyant = .false.
	else
	  ny = inty
	  swyant = .true.
	endif
	ngry = abs(ny)
	if(ny.lt.0) idfr = idfr + 1
c**
c**	define parameters for x-axis
c**
	xfact = rat*xycor
	xaxis = yaxis*xfact		
c !size of x-axis
	xmin = xr(1)
	xmax = xr(2)
	xgr = ngrx
	if(idfr.gt.2) goto 5
c**
c**	linear x-axis
c**
	swxlog = .false.
	call nice(xmin,xmax,xgr)
	tmp = (xmax-xmin) / xgr
	nx = tmp + .5
	disx = xaxis / tmp
	goto 10
c**
c**	logarithmic x-axis
c**
5	swxlog = .true.
	ipwx = alog10(xmax)
	if((10.**ipwx).ne.xmax) then
	  if(xmax.gt.1.) ipwx = ipwx + 1
	endif
	xmax = 10.**ipwx
	if(xmin.le.0.) goto 7
	i = alog10(xmin)
	if (xmin.lt.1.) i = i - 1
	nx = ipwx - i
	if(nx.gt.ngrx) goto 7
	xmin = 10.**i
	xgr = nx
	goto 8
7	xmin = 10.**(ipwx-ngrx)
	nx = ngrx
8	nx = -nx
	disx = xaxis / xgr
c**
c**	define parameters for y-axis
c**
10	ymin = yr(1)
	ymax = yr(2)
	ygr = ngry
	if(mod(idfr,2).eq.0) goto 15
c**
c**	linear y-axis
c**
	swylog = .false.
	call nice(ymin,ymax,ygr)
	tmp = (ymax-ymin) / ygr
	ny = tmp + .5
	disy = yaxis / tmp
	goto 20
c**
c**	logarithmic y-axis
c**
15	swylog = .true.
	ipwy = alog10(ymax)
	if((10.**ipwy).ne.ymax) then
	  if(ymax.gt.1.) ipwy = ipwy + 1
	endif
	ymax = 10.**ipwy
	if(ymin.le.0.) goto 17
	i = alog10(ymin)
	if(ymin.lt.1.) i = i - 1
	ny = ipwy - i
	if(ny.gt.ngry) goto 17
	ymin = 10.**i
	ygr = ny
	goto 18
17	ymin = 10.**(ipwy-ngry)
	ny = ngry
18	ny = -ny
	disy = yaxis / ygr
c**
c**	plot frame and text
c**
20	chszd2 = chsz / 2.
	call plot(x0,y0,-3)		
c !define origin
	if (.not.swxant) disx = xaxis
	if (.not.swyant) disy = yaxis
	call plfram(idfr,xaxis,yaxis,disx,disy,grdsz,0.,0.)
	if (swxant) call axtext(0.,0.,xaxis,nx,0.,cxt,chsz,dtx,ndgt,
     *	                        xmin,xgr,chszd2,dvx)
	if (swyant) call axtext(0.,0.,-yaxis,ny,90.,cyt,chsz,dty,-ndgt,
     *	                        ymin,ygr,chszd2,dvy)
c**
c**	plot heading
c**
	ntit = len(ctit)
	xpdat = xaxis - (ncdat*chszd2)
	if(ntit.ne.0) then
	  nch = ntit
	  nchmx = (xpdat/chsz) - 1
	  indx = 1
	  if(nch.le.nchmx) goto 25
	  do i=nchmx,1,-1
	    iindx = i
	    if(ctit(i:i).eq.' ') goto 23
	  enddo
23	  indx = iindx + 1
	  nch = nch - iindx
	  if(nch.gt.nchmx) nch = nchmx
	  call symbol(0.,(yaxis+yhd1),chsz,ctit,0.,(iindx-1))
25	  call symbol(0.,(yaxis+yhd2),chsz,ctit(indx:),0.,nch)
	endif
	call date_and_time(cdate(7:15),ctime(7:16),czone,timval)
c	call date(cdate(7:15))
c	call time(ctime(7:14))
	if (swxant.or.swyant) then
	  call symbol(xpdat,(yaxis+yhd1/2.),chszd2,cdate,0.,ncdat)
	  call symbol(xpdat,(yaxis+yhd2/2.),chszd2,ctime,0.,ncdat)
	endif
c**
c**
	swfrst = .false.
	return
	end

c**
c**	pkg_xcoord / pkg_ycoord -- function routines to convert a value to
c**	                           a coordinate in a predefined frame.
c**
c**	usage:
c**	   xpos = pkg_xcoord(val,fact)
c**	   ypos = pkg_ycoord(val,fact)
c**
c**	parameters:
c**	   val  - value to be converted
c**	   fact - scale factor
c**
c**
	function pkg_xcoord(val,fact)
c**
c**
	logical swopen,swfrst,swxlog,swylog,switch
	common /pkgpar/ swopen,swfrst,swxlog,swylog,xmin,xmax,ymin,ymax,
     *	                xfact,idset,redy
c**
c**
c**
c**	convert value to x-position
c**
	id = 0
	switch = swxlog
	rmin = xmin
	goto 5
c**
c**	convert value to y-position
c**
	entry pkg_ycoord(val,fact)
c**
c**
	id = 1
	switch = swylog
	rmin = ymin
c**
c**	common code
c**
5	if(.not.switch) then
	  pkg_xcoord = val - rmin
	else
	  pkg_xcoord = val / rmin
	  if(pkg_xcoord.le.0.) then
	    pkg_xcoord = 0.
	  else
	    pkg_xcoord = log10(pkg_xcoord)
	  endif
	endif
	pkg_xcoord = pkg_xcoord * fact
	if(id.ne.0) pkg_ycoord = pkg_xcoord
	return
	end
c**
c**	nice -- routine to determine for a range of values a linear scale in
c**	        such a way, that is readable for the user. a readable scale
c**	        is defined here as a scale with interval size a product of an
c**	        integer power of 10 and 1,2 or 5, and scale values integer
c**	        multiples of the interval size.
c**	        the definition of readability used for 'nice' permits scale
c**	        values such as:
c**	           -0.5,0.0,0.5,1.0,.....
c**	           1.24,1.26,1.28,.....
c**	           100.0,200.0,300.0,.....,etc.
c**	        it prohibits the following examples:
c**	           -1.0,4.0,9.0,.....
c**	           1.2,1.31,1.42,.....
c**	           0.0,4.0,8.0,12.0,.....,etc.
c**
c**	calling sequence:
c**	   call nice(rmin,rmax,stepz)
c**
c**	parameters:
c**	   rmin  - minimum value of range (will be replaced by a new
c**	           minimum value)
c**	   rmax  - maximum value of range (will be replaced by a new
c**	           maximum value)
c**	   stepz - number of approximately required grid intervals (will
c**	           be replaced by the calculated interval size)
c**
c**	notes:
c**	   1. mind that the arguments must be variables and not constants,
c**	      since they might be replaced by new values.
c**	   2. the new rmin and rmax values are defined as an integer
c**	      multiple of the interval size
c**	   3. also see: 'communications of acm', vol. 16, num. 10, algorithm
c**	      463 scale1.
c**
c**
	subroutine nice(xmin,xmax,xstep)
c**
	data s50,s10,s2/7.071068,3.162278,1.414214/
c***
	a = (xmax-xmin)/xstep
	rm1 = alog10(a)
	nal = int(rm1)
	if(a.lt.1.) nal=nal-1
	r = a/(10.**nal)
	i = 10
	if(r.lt.s50) i=5
	if(r.lt.s10) i=2
	if(r.lt.s2)  i=1
	xstep = 10.**nal  *i
	xd = xmin/xstep
	xmin = int(xd)
	if(xd.lt.0.) xmin=xmin - 1.
	xn = 0.
	if((xmin+1.-xd).lt.2.e-5) xn=1.
	xmin = xstep*(xmin+xn)
	xd = xmax/xstep
	xmax = int(xmax/xstep+1.)
	if(xd.lt.(-1.)) xmax=xmax-1.
	xn = 0.
	if((abs(xd)+1.-xmax).lt.2.e-5) xn=1.
	xmax = xstep*(xmax-xn)
	return
	end
c**
c**	plfram -- routine to plot a frame
c**
c**	calling sequence:
c**	   call plfram(ityp,xlen,ylen,disx,disy,grdsz,orx,ory)
c**
c**	parameters:
c**	   ityp - type of frame   1 = x-axis lin  - y-axis lin
c**	                          2 = x-axis lin  - y-axis log
c**	                          3 = x-axis log  - y-axis lin
c**	                          4 = x-axis log  - y-axis log
c**	   all following parameter values must be given in the
c**	   indicated scale in user's units
c**	   xlen - length of x-axis (negative if direction negative)
c**	   ylen - length of y-axis (negative if direction negative)
c**	   disx - lin : distance between thick marks along x-axis
c**	          log : decade length along the x-axis
c**	   disy - same as disx along the y-axis
c**	   grdsz- size of gridmarks
c**	   orx  - x-coord frame origin related to page origin
c**	   ory  - y-coord frame origin related to page origin
c**
c**
	subroutine plfram(ityp,xlen,ylen,disx,disy,grdsz,orx,ory)
c**
	dimension xlog(10),ylog(10)
c**
c**
c**
	indx = 1
	indy = 1
c**
c**	indx and indy are indicators whether the direction
c**	of x- or y-axis is positive or negative.
c**	they are tested on xlen and ylen.
c**
	if(xlen.le.0) indx = -1
	if(ylen.le.0) indy = -1
c**
c**	arrays xlog and ylog are filled with distances for
c**	logarithmical scale numbers 2 to 10
c**
	do i=2,10
	  xx = i
	  xx = log10(xx)
	  xlog(i) = xx * disx * indx
	  ylog(i) = xx * disy * indy
	enddo
	nx = 10
	ny = 10
c**
c**	move pen to origin of frame
c**
	call plot(orx,ory,3)
c**
c**	pen has been moved to origin
c**
	npx = (xlen+0.001)/disx * indx
	npx = npx + 1
	rxlen = npx * disx * indx
	npy = (ylen+0.001)/disy * indy
	npy = npy + 1
	rylen = npy * disy * indy
c**
c**	npx and npy are number of marks along an axis(lin) or
c**	number of decades (log) given in integer value.
c**	rxlen and rylen are places of the last mark
c**
	goto (4,4,5,5),ityp
4	xlog(2) = disx*indx
	nx = 2
5	goto (6,7,6,7),ityp
6	ylog(2) = disy*indy
	ny = 2
7	xml = -grdsz * indy
	yml = -grdsz * indx
c**
c**	xml and yml are  lengthes of marks along a log axis.
c**	first part of x-axis
c**	--------------------
c**
	yy = 0.
	fac = -1.
	id = 0
c**
c**	draw x-axis
c**
15	do j=1,npx
	  j1 = j - 1
	  if (id.ne.0) j1 = npx - j
	  xr = j1 * disx * indx
	  do i=2,nx
	    ir = i
	    if (id.ne.0) ir = 12 - i
	    xx = xlog(ir) + xr
	    if(abs(xx).le.abs(xlen)) then
	      xp = orx + xx
	      yp = ory + yy
	      do ih=1,3
	        if(ih.eq.2) call plot(xp,(yp+fac*xml),2)
	        call plot(xp,yp,2)
	      enddo
	    endif
	  enddo
	enddo
	if(id.ne.0) goto 30
	call plot((orx+xlen),ory,2)
c**
c**	first part of y-axis
c**	--------------------
c**
	xx = xlen
	fac = +1.
c**
c**	draw y-axis
c**
25	do j=1,npy
	  j1 = j - 1
	  if(id.ne.0) j1 = npy - j
	  yr = j1 * disy * indy
	  do i=2,ny
	    ir = i
	    if (id.ne.0) ir = 12 - i
	    yy = ylog(ir) + yr
	    if (abs(yy).le.abs(ylen)) then
	      xp = orx + xx
	      yp = ory + yy
	      do ih=1,3
	        if(ih.eq.2) call plot((xp+fac*yml),yp,2)
	        call plot(xp,yp,2)
	      enddo
	    endif
	  enddo
	enddo
	if (id.ne.0) goto 35
	call plot((orx+xlen),(ory+ylen),2)
c**
c**	second part of x-axis
c**	---------------------
c**
	yy = ylen
	id = 1
	goto 15
30	call plot(orx,(ory+ylen),2)
c**
c**	second part of y-axis
c**	---------------------
c**
	xx = 0.
	fac = -1.
	goto 25
35	call plot(orx,ory,2)
c**
c**
	return
	end

	subroutine pkg_frame_actual_limits(intx,inty,xr,yr)
	integer*4 intx,inty
	real*4    xr(2),yr(2)
c
	integer*4 nx,ny,ngrx,ngry,idfr,ipwx,ipwy,i
	real*4    xmin,xmax,ymin,ymax,xgr,ygr,tmp

	idfr = 1
	nx = intx
	if (intx.eq.0) nx = 20
	ngrx = abs(nx)
	if(nx.lt.0) idfr = 3
	ny = inty
	if (inty.eq.0) ny = 20
	ngry = abs(ny)
	if(ny.lt.0) idfr = idfr + 1
	xmin = xr(1)
	xmax = xr(2)
	xgr  = ngrx
	if(idfr.gt.2) goto 5
	call nice(xmin,xmax,xgr)
	tmp = (xmax-xmin) / xgr
	nx = tmp + .5
	goto 10
5	ipwx = alog10(xmax)
	if((10.**ipwx).ne.xmax) then
	  if(xmax.gt.1.) ipwx = ipwx + 1
	endif
	xmax = 10.**ipwx
	if(xmin.le.0.) goto 7
	i = alog10(xmin)
	if (xmin.lt.1.) i = i - 1
	nx = ipwx - i
	if(nx.gt.ngrx) goto 7
	xmin = 10.**i
	xgr = nx
	goto 8
7	xmin = 10.**(ipwx-ngrx)
	nx = ngrx
8	nx = -nx
10	ymin = yr(1)
	ymax = yr(2)
	ygr = ngry
	if(mod(idfr,2).eq.0) goto 15
	call nice(ymin,ymax,ygr)
	tmp = (ymax-ymin) / ygr
	ny = tmp + .5
	goto 20
15	ipwy = alog10(ymax)
	if((10.**ipwy).ne.ymax) then
	  if(ymax.gt.1.) ipwy = ipwy + 1
	endif
	ymax = 10.**ipwy
	if(ymin.le.0.) goto 17
	i = alog10(ymin)
	if(ymin.lt.1.) i = i - 1
	ny = ipwy - i
	if(ny.gt.ngry) goto 17
	ymin = 10.**i
	ygr = ny
	goto 18
17	ymin = 10.**(ipwy-ngry)
	ny = ngry
18	ny = -ny
20	intx  = nx
	inty  = ny
	xr(1) = xmin
	xr(2) = xmax
	yr(1) = ymin
	yr(2) = ymax
	return
	end

	subroutine pkg_plhist(x,y,ni)
	dimension x(*),y(*),xr(2),yr(2)
	il = 2
	if (ni.lt.0) il = -2
	n = abs(ni)
	if (n.lt.2) return
	do i = 1,n-1
	  xr(1) = x(i)
	  xr(2) = (x(i)+x(i+1)) / 2.
	  yr(1) = y(i)
	  yr(2) = y(i)
	  call pkg_pldatr(-1,xr,yr,il)
	  xr(1) = xr(2)
	  yr(1) = y(i)
	  yr(2) = y(i+1)
	  call pkg_pldatr(-1,xr,yr,il)
	  xr(2) = x(i+1)
	  yr(1) = y(i+1)
	  yr(2) = yr(1)
	  call pkg_pldatr(-1,xr,yr,il)
	enddo
	xr(2) = x(1)
	xr(1) = x(1) - abs(x(2) - x(1))/2.
	yr(1) = y(1)
	yr(2) = yr(1)
	call pkg_pldatr(-1,xr,yr,2)
	yr(1) = y(n)
	yr(2) = yr(1)
	xr(1) = x(n)
	xr(2) = x(n) + abs(x(n) - x(n-1)) / 2.
	call pkg_pldatr(-1,xr,yr,il)
	return
	end

      subroutine pkg_plcross(x,y,sigma,ni)
      dimension x(*),y(*),sigma(1),xr(2),yr(2)
      il = 2
      if (ni.lt.0) il = -2
      n = abs(ni)
      if (n.lt.2) return
      do i = 1,n
        xr(1) = x(i)
        xr(2) = xr(1)
        yr(1) = y(i) + sigma(i)
        yr(2) = y(i) - sigma(i)
        call pkg_pldatr(-1,xr,yr,il)
      enddo
      do i = 2,n-1
        yr(1) = y(i)
        yr(2) = yr(1)
        xr(1) = (x(i) + x(i-1)) / 2.
        xr(2) = (x(i) + x(i+1)) / 2.
        call pkg_pldatr(-1,xr,yr,il)
      enddo
c-- Artifically add the 'horizontal' bar for the first and the last bin.
      yr(1) = y(1)
      yr(2) = yr(1)
      xr(2) = (x(1) + x(2)) / 2.
      xr(1) = x(1) - abs(xr(2) - x(1))
      call pkg_pldatr(-1,xr,yr,il)
      yr(1) = y(n)
      yr(2) = yr(1)
      xr(1) = (x(n-1) + x(n)) / 2.
      xr(2) = x(n) + abs(xr(1) - x(n))
      call pkg_pldatr(-1,xr,yr,il)
      return
      end

      subroutine pkg_plpois(x,y,ni)
      dimension x(*),y(*),xr(2),yr(2)
      il = 2
      if (ni.lt.0) il = -2
      n = abs(ni)
      if (n.lt.2) return
      do i = 1,n
        xr(1) = x(i)
        xr(2) = xr(1)
        yr(1) = y(i) + sqrt(abs(y(i)))
        yr(2) = y(i) - sqrt(abs(y(i)))
        call pkg_pldatr(-1,xr,yr,il)
      enddo
      do i = 2,n-1
        yr(1) = y(i)
        yr(2) = yr(1)
        xr(1) = (x(i) + x(i-1)) / 2.
        xr(2) = (x(i) + x(i+1)) / 2.
        call pkg_pldatr(-1,xr,yr,il)
      enddo
c-- Artifically add the 'horizontal' bar for the first and the last bin.
      yr(1) = y(1)
      yr(2) = yr(1)
      xr(2) = (x(1) + x(2)) / 2.
      xr(1) = x(1) - abs(xr(2) - x(1))
      call pkg_pldatr(-1,xr,yr,il)
      yr(1) = y(n)
      yr(2) = yr(1)
      xr(1) = (x(n-1) + x(n)) / 2.
      xr(2) = x(n) + abs(xr(1) - x(n))
      call pkg_pldatr(-1,xr,yr,il)
      return
      end

	subroutine pkg_raster(ix,iy,xri,yri,idash)
	real*4 xri(2),yri(2)
	real*4 xd(2),yd(2),xr(2),yr(2)
	xr(1) = xri(1)
	xr(2) = xri(2)
	yr(1) = yri(1)
	yr(2) = yri(2)
	if (ix.lt.0) then
	  xr(2) = 10.**float(int(alog10(xr(2)))+1)
	  xr(1) = 10.**float(ix) * xr(2)
	else
	  xstep = float(ix)
	  call nice(xr(1),xr(2),xstep)
	endif
	xd(1) = xr(1)
	xd(2) = xr(2)
	if (iy.lt.0) then
	  yr(2) = 10.**float(int(alog10(yr(2)))+1)
	  yr(1) = 10.**float(iy) * yr(2)
	  ystep = yr(1)
	  yd(1) = yr(1)
	  do k = 1,abs(iy)
	    do i = 1,9
	      yd(1) = yd(1)  + ystep
	      yd(2) = yd(1)
	      call pkg_pldatr(-1,xd,yd,idash)
	    enddo
	    ystep = 10. * ystep
	  enddo
	else
	  ystep = float(iy)
	  call nice(yr(1),yr(2),ystep)
	  yd(1) = yr(1)
 3	  if (yd(1).le.yr(2)) then
	    yd(1) = yd(1) + ystep
	    yd(2) = yd(1)
	    call pkg_pldatr(-1,xd,yd,idash)
	    goto 3
	  endif
	endif
	xr(1) = xri(1)
	xr(2) = xri(2)
	yr(1) = yri(1)
	yr(2) = yri(2)
	if (iy.lt.0) then
	  yr(2) = 10.**float(int(alog10(yr(2)))+1)
	  yr(1) = 10.**float(iy) * yr(2)
	else
	  ystep = float(iy)
	  call nice(yr(1),yr(2),ystep)
	endif
	yd(1) = yr(1)
	yd(2) = yr(2)
	if (ix.lt.0) then
	  xr(2) = 10.**float(int(alog10(xr(2)))+1)
	  xr(1) = 10.**float(ix) * xr(2)
	  xstep = xr(1)
	  xd(1) = xr(1)
	  do k = 1,abs(ix)
	    do i = 1,9
	      xd(1) = xd(1)  + xstep
	      xd(2) = xd(1)
	      call pkg_pldatr(-1,xd,yd,idash)
	    enddo
	    xstep = 10. * xstep
	  enddo
	else
	  xstep = float(ix)
	  call nice(xr(1),xr(2),xstep)
	  xd(1) = xr(1)
 6	  if (xd(1).le.xr(2)) then
	    xd(1) = xd(1) + xstep
	    xd(2) = xd(1)
	    call pkg_pldatr(-1,xd,yd,idash)
	    goto 6
	  endif
	endif
	return
	end

c******************************************************** sron/rol *********
c
c  with this (improvised ad/hoc fred special) subroutine you can plot a
c  text block below your plots. in a normal rat = 1.0 , non-(p)reset plot
c  you can get up to 8 text lines.
c  
c  ctxt = text string to be plotted max 79 characters !!!!!!!! (%ref)
c  nrl  = line nr. for text line value : 1 <= nrl <= 8
c         n.b. line nr. 1 is closest to the frame and nr. 8 is the 
c              farthest.
c
c  written by : fred a. jansen
c               oct. 13, 1985    version : 1.0
c
c***************************************************************************
	subroutine pkg_pltextbl(ctxt,nrl)
	character*(*) ctxt
	logical swopen,swfrst,swxlog,swylog
	common /frparm/ pwidth,yfrac,yshft,x0,y0,yaxis,grdsz,symsz,dvx,
     *	                dtx,yhd1,yhd2,dvy,dty,chsz,dl,ndgt
	common /pkgpar/ swopen,swfrst,swxlog,swylog,xmin,xmax,ymin,ymax,
     *	                xfact,idset,redy
	xp = chsz
	yp = - dtx - (nrl+1) * chsz
	call symbol(xp,yp,chsz/2.,ctxt,ang,len(ctxt))
	return
	end

c**
c**	axtext -- plot routine to fully-annotate a single axis complete with a
c**	          centered title. the axis may be at any location and angle
c**	          on the plot.
c**
c**	calling sequence:
c**	   call axtext(xa,ya,axlen,ndiv,angle,ctit,tsiz,tdis,idt,
c**			sval,dval,vsiz,vdis)
c**
c**	parameters:
c**	   (all plot parameters are supposed to be in user units)
c**	   xa,ya - the coordinates of the start of the axis
c**	   axlen - the absolute value of axlen is the length of the axis.
c**	           if axlen is positive, all annotation appears on the clock-
c**	           wise side of the axis; if axlen is negative, it appears on
c**	           on the counterclockwise side.
c**	   ndiv  - the absolute value of ndiv is the number of axis divisions.
c**	           if ndiv is positive, the axis is a linear axis; if ndiv
c**	           is negative, the axis is a logarithmic axis.
c**	   angle - the angle in degrees that the axis line is to be rotated
c**	           about xa,ya measured counterclockwise from the horizontal.
c**	   ctit  - character variable with string to be placed along the axis
c**	           as a title. the title is parallel to and centered on the
c**		   axis.
c**	   tsiz  - the character height for the title
c**	   tdis  - the distance between the axis and the title
c**	   idt   - the absolute value of idt is the number of significant
c**	           digits to be plotted for each tick mark value. if idt is
c**	           positive, the values are written parallel to the axis; if
c**	           idt is negative, the values are written perpendicular on
c**	           the axis.
c**	   sval  - the value to be assigned to the first tick mark on the axis
c**	   dval  - the increment for linear values; meaningless for logarithmic
c**	           values
c**	   vsiz  - the character height for the annotation values
c**	   vdis  - the distance between the axis and the annotation
c**
c**
c**	if, for any value, the number of digits to be plotted exceeds the
c**	absolute value of idt, the title is extended with an offset and a
c**	factor, so that the annotation value can be represented as:
c**	            val = (title-offset) * 10**n
c**
c**
	subroutine axtext(xa,ya,axlen,ndiv,angle,ctit,tsiz,tdis,idt,
     *	                  sval,dval,vsiz,vdis)
c**
	character*(*) ctit
	character*1   csign,c123(3)
	character*4   cfmti
	character*8   cfmt
	character*80  cbuf
c**
	data cfmti/'(i*)'/,c123/'1','2','3'/,torad/.01745329/
c**
c**	initialization
c**
	xx    = xa
	yy    = ya
	axsiz = abs(axlen)
	side  = axsiz / axlen
	i     = abs(ndiv)
	xint  = axsiz / i
	ngr   = i + 1
	ang   = angle * torad
	sina  = sin(ang)
	cosa  = cos(ang)
	theta = angle
	if(idt.lt.0) theta = angle - 90.
	ndgt  = abs(idt)
	ntit  = len(ctit)
	cbuf  = '( '
	if(ndiv.lt.0) goto 10
c**
c******************************************************************************
c**		 l i n e a r   a x i s
c**
c**	build d-field of format statement:  format(fw.d)
c**	================================================
c**
	rmin = sval
	rmax = rmin + ndiv*dval
	amin = abs(rmin)
	amax = abs(rmax)
	ipow = log10(dval)
	if(dval.lt.1.) ipow = ipow - 1
	incr = dval / (10.**ipow) +.5
	if(incr.eq.10) ipow = ipow + 1
	cfmt = '(f**.0 )'
	id = 0
	if(ipow.lt.0) id = -ipow
	iw = id + 1
c**
c**	calculate maximum field length, factor and offset
c**	=================================================
c**
	isign = 0
	i1 = 0
	if(amin.ge.1.) i1 = log10(amin) + 1.
	if(rmin.lt.0.) i1 = i1 + 1
	ndgtmx = 0
	if(amax.ge.1.) ndgtmx = log10(amax) + 1.
	if(ndgtmx.lt.i1) then
	  isign = 1
	  ndgtmx = i1
	endif
	ndgtmx = ndgtmx + id
	if(id.ne.0) ndgtmx = ndgtmx + 1
	ipwf = 0
	ofs = 0.
	if(ndgtmx.gt.ndgt) then
	  if(ipow.ne.0) then
	    ipwf = -ipow
	    id = 0
	    iw = 1
	  endif
	  i1 = -ipow - ndgt + isign
	  tmp = 10.**i1
	  i1 = rmin * tmp
	  ofs = i1 / tmp
	endif
	write(cfmt(6:7),'(i2)') id
	fact = 10.**ipwf
	csign = '-'
	if(ofs.lt.0.) csign = '+'
c**
c**	build character string to plot title along the axis
c**	===================================================
c**
	nt = 0
c**
c**	build title
c**
	if(ntit.ne.0) then
	  if((ofs.ne.0.) .and. (ipwf.ne.0)) nt = 1
	  cbuf(nt+1:nt+ntit) = ctit
	  nt = nt + ntit
	endif
c**
c**	build offset
c**
	if(ofs.ne.0) then
	  if(ntit.ne.0) nt = nt + 1
	  nt = nt + 1
	  cbuf(nt:nt) = csign
	  if(ofs.gt.0.) nt = nt + 1
	  tmp = abs(ofs)
	  j = 0
	  if(tmp.ge.10.) j = log10(tmp)
	  ilng = iw + j + 1
	  if(ofs.lt.0.) ilng = ilng + 1
	  write(cfmt(3:4),'(i2)') ilng
	  nt1 = nt + 1
	  write(cbuf(nt1:nt+ilng),cfmt) ofs
	  if(ofs.lt.0.) cbuf(nt1:nt1) = ' '
	  nt = nt + ilng
	  if(id.eq.0) nt = nt - 1
	endif
	if(ipwf.eq.0) goto 15
c**
c**	build factor
c**
	if((ntit.ne.0) .and. (ofs.ne.0.)) then
	  nt = nt + 1
	  cbuf(nt:nt) = ')'
	endif
	if((ntit.ne.0) .or. (ofs.ne.0.)) nt = nt + 1
	cbuf(nt:nt+5) = '* 10**'
	nt = nt + 6
	if(ipwf.eq.1) then
	  nt = nt - 2
	  goto 15
	endif
	ilng = log10(float(iabs(ipwf))) + 1
	if(ipwf.lt.0) ilng = ilng + 1
	cfmti(3:3) = c123(ilng)
	nt1 = nt + 1
	write(cbuf(nt1:nt+ilng),cfmti) ipwf
	nt = nt + ilng
	goto 15
c**
c******************************************************************************
c**		 l o g a r i t h m i c   a x i s
c**
c**	calculate first value to be written
c**	===================================
c**
10	tmp = log10(sval)
	if(tmp.eq.0.) then
	  ipow = 0
	else
	  ipow = tmp + .5*abs(tmp)/tmp
	endif
	incr = 1
	indx = 1
	iw = 3
c**
c**	build character string to plot the title along the axis
c**	=======================================================
c**
	cbuf = 'log'
	nt = 3
	if(ntit.eq.0) goto 15
	nt = nt + 1
	cbuf(nt:nt) = '('
	cbuf(nt+1:nt+ntit) = ctit
	nt = nt + ntit + 1
	cbuf(nt:nt) = ')'
c**
c******************************************************************************
c**		 p l o t   t i t l e   a l o n g   t h e   a x i s
c**
15	if(nt.ne.0) then
	  ofs = nt/2. * tsiz
	  axszd2 = .5 * axsiz
	  xt = xx + axszd2 - ofs
	  yt = yy - side*tdis
	  if(side.gt.0.) yt = yt - tsiz
	  rx = xt*cosa - yt*sina
	  ry = yt*cosa + xt*sina
	  call symbol(rx,ry,tsiz,cbuf,angle,nt)
	endif
c**
c******************************************************************************
c**		 t i c k   m a r k   a n n o t a t i o n
c**
c**	calculate start position for first value
c**	========================================
c**
	yy = yy - side*vdis
	if(idt.lt.0) goto 18
c**
c**	parallel to axis
c**
	ijus = 0
	if(side.gt.0.) yy = yy - vsiz
	goto 20
c**
c**	perpendicular on axis
c**
18	ijus = -2. * side
	xx = xx - .5*vsiz
c**
c**	build character string to be plotted
c**	====================================
c**
20	i1 = -1
21	i1 = i1 + 1
	if(i1.eq.ngr) goto 99
	ilng = iw
	if(ndiv.gt.0) goto 22
c**
c**	build format for log value
c**	==========================
c**
	ival = ipow + i1*incr
	if(ival.ge.0) ilng = ilng - 1
	if(abs(ival).lt.10) ilng = ilng - 1
	cfmti(3:3) = c123(ilng)
c**
c**	encode value
c**
	write(cbuf,cfmti) ival
	goto 23
c**
c**	build w-field of format statement for lin value:  format(fw.d)
c**	===============================================================
c**
22	rval = (rmin + i1*dval) * fact
	tmp = abs(rval) + 10.**(-id-4)
	j = 0
	if(tmp.ge.10.) j = log10(tmp)
	ilng = ilng + j + 1
	if(rval.lt.0.) ilng = ilng + 1
	write(cfmt(3:4),'(i2)') ilng
c**
c**	encode value
c**
	write(cbuf,cfmt) rval
	if(id.eq.0) ilng = ilng - 1
23	indx = 1
	nch = ilng
	if(ilng.gt.ndgt) then
	  nch = ndgt
	  indx = ilng + 1 - nch
	  if(rval.lt.0.) cbuf(indx:indx) = cbuf(1:1)
	endif
c**
c**	plot annotation
c**	===============
c**
	xt = xx
	yt = yy
	ofs = nch * vsiz
	if(ijus) 28,26,27
26	xt = xt - .5*ofs
	goto 28
27	yt = yt + ofs
28	rx = xt*cosa - yt*sina
	ry = yt*cosa + xt*sina
	call symbol(rx,ry,vsiz,cbuf(indx:),theta,nch)
	xx = xx + xint
	goto 21
99	return
	end

      subroutine symbol(px,py,ph,istr,pa,n)
c======================================================================
c
c	plot symbol or string of symbols
c
c		px-	x-position (scaled by factor)
c		py-	y-position
c			position is coord. of bottom left, except for 
c			symbol nr 0-8, which are centered around coord.
c		ph-	height (scaled by factor)
c		istr-	character string
c			if n<0: nr of special char.
c			if n>0: ascii string
c		pa-	angle (degrees left from + x-axis)
c		n-	<0: istr is nr of special char.
c			>0: nr of ascii char. in istr
c
c
c method
c	each character is built up in a 7*7 matrix.
c	in a number of bytes in 'insmb1' resides the plotcode to form
c	these characters.
c	every byte gives the following information:
c	bits	0-3	the y-coordinate of a dot within a 7-7 matrix
c	bits	4-6	the x-coordinate of a dot within a 7-7 matrix
c	bit	7	the penposition,0=down,1=up
c
c	for instance	the character '+' is drawn with two lines
c			which are formed by four dots,c.q. four bytes.
c=======================================================================
	integer*2 ismbp(257),smbtab(1130)
	character*(*) istr
	integer*2 ikar,ismb
c  
	common /gsp /ox,oy,xo,xf,yo,yf,scale,iroute
	common /symb/xlast,ylast
c
	data pi/.0174533/
c
c...table pointers
c
	data ismbp/
     *	       1,   7,  17,  21,  25,  29,  35,  41,  45,  49,
     *	      55,  67,  75,  81,  83,  83,  83,  89,  95, 101,
     *	     107, 111, 117, 121, 127, 127, 127, 131, 137, 143,
     *	     149, 155, 161, 161, 165, 173, 181, 193, 205, 217,
     *	     221, 227, 233, 239, 243, 249, 251, 257, 259, 269,
     *	     275, 283, 297, 305, 315, 327, 331, 347, 359, 369,
     *	     381, 385, 389, 393, 403, 419, 427, 439, 447, 453,
     *	     461, 467, 477, 483, 489, 495, 501, 505, 509, 513,
     *	     525, 531, 543, 553, 565, 569, 575, 579, 585, 589,
     *	     595, 599, 603, 605, 609, 611, 613, 619, 629, 637,
     *	     645, 653, 663, 673, 685, 691, 695, 703, 709, 715,
     *	     725, 731, 741, 749, 759, 765, 775, 781, 789, 793,
     *	     803, 807, 811, 815, 827, 829, 841, 847, 847, 847,
     *	     847, 847, 847, 847, 847, 847, 847, 847, 847, 847,
     *	     847, 847, 847, 847, 847, 859, 859, 859, 859, 859,
     *	     859, 859, 859, 859, 859, 859, 859, 859, 859, 859,
     *	     859, 859, 859, 859, 859, 859, 859, 859, 859, 859,
     *	     859, 869, 883, 891, 891, 891, 891, 891, 891, 891,
     *	     891, 891, 891, 891, 891, 891, 891, 903, 913, 921,
     *	     921, 921, 921, 921, 921, 921, 921, 921, 921, 921,
     *	     921, 921, 921, 933, 937, 947, 953, 961, 971, 971,
     *	     971, 971, 971, 971, 971, 971, 971, 971, 971, 981,
     *	     987, 999,1007,1019,1019,1019,1027,1027,1027,1027,
     *	    1027,1027,1027,1027,1027,1031,1039,1051,1067,1079,
     *	    1089,1089,1089,1089,1089,1089,1089,1089,1089,1089,
     *	    1089,1091,1093,1103,1109,1125,1131/
c**
c**	arrays with data points
c**
	data smbtab/
     *	      164,  4,  0, 64, 68, 36,
     *	      164, 20,  3,  1, 16, 48, 65, 67, 52, 36,
     *	      164,  0, 64, 36,
     *	       36, 32,130, 66,
     *	       68,  0,132, 64,
     *	      164,  2, 32, 66, 36, 36,
     *	       36,  2, 66, 36, 32, 32,
     *	        0, 68,  4, 64,
     *	      132, 68,  0, 64,
     *	      132, 34, 68, 34, 32, 32,
     *	      196, 51, 19,  4, 19, 17,  0, 17, 49, 64, 49, 51,
     *	      132, 64,196,  0,130, 66,164, 32,
     *	       68,  4, 64,  0, 34, 34,
     *	       36, 32,
     *	      195, 37, 71, 37,117,117,
     *	      195,101, 71,101, 21, 21,
     *	      195, 72, 38, 72,102,102,
     *	      200, 67, 37, 67,101,101,
     *	      163, 71, 99, 99,
     *	      163, 99,229, 37,167,103,
     *	      167, 67,103,103,
     *	      165, 54, 66, 73,121,121,
     *	      163,103,167, 99,
     *	      163, 99,196, 72,166,102,
     *	      163, 99,228, 38,104,104,
     *	      178, 88,230, 38,164,100,
     *	      163, 99,164,102, 40, 40,
     *	      195, 68,198, 71,165,101,
     !	      194, 67,196, 73,
     "	      168, 56, 54, 40,216,104, 86, 88,
     #	       72,232, 66,228, 36,166,102,102,
     $	       82, 99,100, 85, 53, 38, 39, 56,104,201, 65, 65,
     %	      105,185, 40, 55, 72, 57,212, 67, 82, 99, 84, 84,
     &	      228, 66, 50, 35, 36, 87, 88, 73, 56, 55, 98, 98,
     '	      199, 73, 89, 71,
     (	      208, 65, 51, 54, 72, 89,
     )	      176, 65, 83, 86, 72, 57,
     *	      195, 71,230, 36,166,100,
     +	      195, 71,165,101,
     ,	      177, 66, 67, 51, 50, 66,
     -	      165,101,
     .	      178, 66, 67, 51, 50, 50,
     /	      105,105,
     o	      164, 39, 57, 89,103,100, 82, 50, 36, 36,
     1	      184, 73, 66, 50, 82, 82,
     2	      168, 57, 89,104,103, 35, 34, 98,
     3	      163, 50, 82, 99,101, 86, 54, 86,103,104, 89, 57, 40, 40,
     4	      169, 36,100, 84, 89, 82, 66, 98,
     5	      163, 50, 82, 99,101, 86, 38, 41,105,105,
     6	      166, 86,101, 99, 82, 50, 35, 40, 57, 89,104,104,
     7	      169,105, 67, 66,
     8	      214,103,104, 89, 57, 40, 39, 54, 86,101, 99, 82,
     8             50, 35, 37, 54,
     9	      163, 50, 82, 99,104, 89, 57, 40, 38, 53, 85,102,
     :	      179, 50, 66, 67, 51,182, 70, 71, 55, 54,
     ;	      177, 66, 67, 51, 50, 66,198, 71, 55, 54, 70, 70,
     <	      227, 37,103,103,
     =	      164,100,230, 38,
     >	      163,101, 39, 39,
     ?	      167, 40, 57, 89,104,103, 69, 68,195, 66,
     @	      228, 83, 51, 36, 39, 56, 88,103,101, 84, 68, 53,
     @             54, 71, 87,102,
     a	       40, 57, 89,104, 98,101, 37, 37,
     b	       41, 89,104,103, 86, 38, 86,101, 99, 82, 34, 34,
     c	      232, 89, 57, 40, 35, 50, 82, 99,
     d	       41, 89,104, 99, 82, 34,
     e	       98, 34, 38, 86, 38, 41,105,105,
     f	       38, 86, 38, 41,105,105,
     g	      232, 89, 57, 40, 35, 50, 82, 99,101, 85,
     h	       41, 38,102,105, 98, 98,
     i	      178, 82, 66, 73, 57, 89,
     j	      163, 50, 82, 99,105,105,
     k	       41, 37,105, 71, 98, 98,
     l	       41, 34, 98, 98,
     m	       41, 71,105, 98,
     n	       41, 98,105,105,
     o	      217, 57, 40, 35, 50, 82, 99,104, 89, 89, 89, 89,
     p	       41, 89,104,103, 86, 38,
     q	      210, 99,104, 89, 57, 40, 35, 50, 82,226, 68, 68,
     r	       41, 89,104,103, 86, 38, 86,101, 98, 98,
     s	      163, 50, 82, 99,101, 86, 54, 39, 40, 57, 89,104,
     t	      194, 73, 41,105,
     u	      169, 35, 50, 82, 99,105,
     v	      169, 66,105,105,
     w	       41, 34, 69, 98,105,105,
     x	      105,169, 98, 98,
     y	      194, 71, 41, 71,105,105,
     z	      169,105, 34, 98,
     [	      225, 49, 57,105,
     \	      169, 97,
     ]	      161, 81, 89, 41,
     *	      193, 73,
     *	      133,117,
     a	      165, 54, 69, 52, 37, 37,
     a	          226,101, 86, 54, 37, 35, 50, 82, 99,102,
     b	      169, 34, 82, 99,101, 86, 54, 37,
     c	      227, 82, 50, 35, 37, 54, 86,101,
     d	      233, 98, 50, 35, 37, 54,102,102,
     e	      164,100,101, 86, 54, 37, 35, 50, 82, 99,
     f	      178, 53, 37, 69, 53, 56, 73, 89,104,104,
     g	      161, 48, 80, 97,101, 86, 54, 37, 35, 50, 82, 99,
     h	       41, 37, 54, 86,101, 98,
     i	      194, 70,201, 73,
     j	      161, 48, 64, 81, 86,217, 89, 89,
     k	       41, 36,102, 69, 98, 98,
     l	      194, 73,128,128,128,128,
     m	       38, 37, 54, 69, 67, 69, 86,101, 98, 98,
     n	       38, 37, 54, 86,101, 98,
     o	      178, 82, 99,101, 86, 54, 37, 35, 50,50,
     p	       32, 38, 86,101, 99, 82, 50, 35,
     q	      227, 82, 50, 35, 37, 54, 86,101, 96,112,
     r	       38, 37, 54, 86,101,101,
     s	      163, 50, 82, 99, 84, 52, 37, 54, 86,101,
     t	      166,102,200, 67, 82, 99,
     u	      166, 35, 50, 82, 99, 98,102,102,
     v	      166, 66,102,102,
     w	      166, 35, 50, 67, 69, 67, 82, 99,102,102,
     x	      102,166, 98, 98,
     y	      166, 66, 48,102,
     z	      166,102, 34, 98,
     [	      233, 89, 72, 70, 53, 37, 53, 68, 66, 81, 97, 97,
     *	       32, 41,
     ]	      161, 49, 66, 68, 85,101, 85, 70, 72, 57, 41, 41,
     ^	      149, 38, 54, 84,100,117,
     *	      148, 37, 53, 83, 99,116,246,101, 85, 55, 39, 22,
     *	      246, 99, 82, 50, 35, 37, 54, 86,101,114,
     *	      160, 55, 72, 88,103,102, 85, 69, 85,100, 99, 82,
     *             66, 51,
     *	      167, 56, 72, 85,121, 85, 48, 48,
     *	      198, 54, 37, 35, 50, 66, 83, 85, 70, 56, 73,104,
     *	      210, 50, 35, 36, 84, 36, 37, 54, 86, 86,
     *	       53, 38, 53, 70, 86,101, 80, 80,
     *	      165,101,103, 88, 72, 55, 37, 35, 50, 66, 83,101,
     *	       70, 41, 70, 98,
     *	      160, 54, 51, 66, 82, 99,102, 99,114,114,
     *	      166, 54, 34, 67, 84,102,
     *	      165, 54, 70, 50, 70,118,102, 98,
     *	      160, 53, 70, 86,101, 99, 82, 66, 51, 51,
     *	      198, 85, 83, 66, 50, 35, 37, 54,103,103,
     *	      165, 54,118, 86, 66, 66,
     *	      161,105, 87, 55, 37, 36, 51, 83,101,102, 87, 87,
     *	      160, 89,183, 36, 51, 67, 84,119,
     *	      182, 37, 35, 50, 67, 69, 67, 82, 99,101, 86, 86,
     *	      161, 97, 81, 89,105, 41, 57, 49,
     *	       18, 72,114, 34,
     *	      226, 97, 33, 86, 41,105,104,104,
     *	      178, 82, 66, 70, 56, 41, 24,198, 88,105,120,120,
     *	      178, 82, 66, 73, 57, 89,215,102,101, 84, 52, 37,
     *             38, 55, 87, 87,
     *	      178, 82, 66, 73, 57, 89,167, 37, 52, 84,101,103,
     *	       50, 37, 39, 56, 88,103,101, 82, 98, 98,
     *	      145,113,
     *	      153,121,
     *	      177, 89,230, 87, 55, 38, 36, 51, 83,100,
     *	      164, 84,101,102, 87, 39,
     *	      197, 52, 36, 21, 22, 39, 55, 70, 69, 84,100,117,
     *            118,103, 87, 70,
     *         72, 98, 22,118, 34, 34/
c**
c**
	px1 = px
	py1 = py
	if (px .eq. 999.) px1=xlast
	if (py .eq. 999.) py1=ylast
c
	cx = px1
	cy = py1
	sz = ph/7.
	sa = sin(pi*pa)*sz
	ca = cos(pi*pa)*sz
	ipen=3
	if (n .le. -2) ipen=2
	ic=n
	if(n .le. 0) ic=1
c
c...to scope (iroute temp.set non-spooled)
c
	  do 5 k=1,ic
	     call plot(cx,cy,ipen)
	     ikar=ichar(istr(k:k))
	     ikar=ikar+1
	     ifs=ismbp(ikar)
	     ils=ismbp(ikar+1)-1
	     if (ifs .ge. ils) goto 32
c
		do 10 i=ifs,ils
		   ismb = smbtab(i)
		   ix=iand(ishft(ismb,-4) , 7 ) -2
		   iy=iand(ismb,15)-2
		   xt=cx+ca*ix-sa*iy
		   yt=cy+sa*ix+ca*iy
		   ipen=ishft(ismb,-7)+2
		   call plot(xt,yt,ipen)
10		continue
c
32	     cx=cx+7.*ca
	     cy=cy+7.*sa
	     ipen=3
5	  continue
c
c...determine pen position after string
c
999	xlast=px1+ic*7.*ca
	ylast=py1+ic*7.*sa
c
	call plot(xlast,ylast,3)
	return
	end

c
c-- Specially modified 'postscript only' version of plot
c
	subroutine plot(xin,yin,pen)
c
	real      xin, yin
	integer   pen
c
	real      xlow,xhigh,ylow,yhigh,x1,x2,y1,y2,fplot,fval,
     *            xl,xh,yl,yh,x_low,x_high,y_low,y_high,x_or,y_or,
     *            ax,bx,ay,by,x,y,xlen,ylen
	byte      init
	integer      plotrec(5)
	equivalence (plotrec(1),x1),(plotrec(2),y1),
     *              (plotrec(3),x2),(plotrec(4),y2)
c
	integer   lunpl,plunit
	character*(*) fname
c
	integer color, color_in
	integer ix1,ix2,iy1,iy2,i_x1,i_y1,ix_m,iy_m,i_xd,i_yd
	logical swinquire
c
	parameter (nvar=14)
c
	integer    num_context
	parameter (num_context=50)
	integer context(num_context),id_context
c-- used to be a structure, but g77 does not like it
	integer     pdsixlow(num_context)
	integer     pdsiylow(num_context)
	integer     pdsixhigh(num_context)
	integer     pdsiyhigh(num_context)
	real*4      pdsxlow(num_context)
	real*4      pdsxhigh(num_context)
	real*4      pdsylow(num_context)
	real*4      pdsyhigh(num_context)
	real*4      pdsxor(num_context)
	real*4      pdsyor(num_context)
	real*4      pdsfplot(num_context)
	real*4      pdsxlen(num_context)
	real*4      pdsylen(num_context)
	real*4      pdspwidth(num_context)
	real*4      pdsyfrac(num_context)
	real*4      pdsyshft(num_context)
	real*4      pdsparm(nvar,num_context)
	logical     pdsswopen(num_context)
	logical     pdsswfrst(num_context)
	logical     pdsswxlog(num_context)
	logical     pdsswylog(num_context)
	logical     pdsswrotateleft(num_context)
	real*4      pdsxmin(num_context)
	real*4      pdsxmax(num_context)
	real*4      pdsymin(num_context)
	real*4      pdsymax(num_context)
	real*4      pdsxfact(num_context)
	integer     pdsidset(num_context)
	real*4      pdsredy(num_context)
	integer     pdsnxd(num_context)
	integer     pdsnyd(num_context)
	integer     pdsiloc(num_context)
c-- end of former structure
c
	integer*4  nr_divisions
	parameter (nr_divisions=8)
	integer*4 plt_nx,plt_ny,plt_ixlow,plt_ixhigh,
     *            plt_iylow,plt_iyhigh,plt_nxd,plt_nyd,plt_csto,
     *            divisions(nr_divisions)
	real*4    plt_x_or,plt_y_or,x_click1,x_click2,y_click1,y_click2
	common /PS_plot/plt_nx,plt_ny,plt_ixlow,plt_ixhigh,
     *              plt_iylow,plt_iyhigh,plt_x_or,plt_y_or,
     *              plt_nxd,plt_nyd,divisions,plt_csto,
     *              x_click1,x_click2,y_click1,y_click2
c
	integer*4  n_colors, n_total
	parameter (n_colors=128,n_total=n_colors+6)
	byte        cmscolor(3,n_total)
	common /PS_color/cmscolor
c
	logical         swopen,swfrst,swxlog,swylog,swrotateleft
	real*4          pwidth,yfrac,yshft,parm
	real*4          xmin,xmax,ymin,ymax,xfact,redy
	integer*4       idset
	common /frparm/ pwidth,yfrac,yshft,parm(nvar)
	common /pkgpar/ swopen,swfrst,swxlog,swylog,xmin,xmax,ymin,ymax,
     *                  xfact,idset,redy
	data init/0/,fplot/1./
	data xlow/0./,xhigh/1./,ylow/0./,yhigh/1./
	data x_low/0./,y_low/0./,x_high/1./,y_high/1./
	data x_or/0./,y_or/0./
	data swrotateleft/.false./
	data plt_csto/0/
c
        save plunit
c
	if (init.eq.0) then
	  xlen     = parm(3) * xfact
	  ylen     = parm(3)
	  init     = 1
	  plt_x_or = 0.
	  plt_y_or = 0.
	  if (plt_csto.eq.0) plt_csto = n_colors + 1
	  if (.not.swrotateleft) then
	    call get_PS_WIDTH_LIMITS (plt_ixlow,plt_ixhigh)
	    call get_PS_HEIGHT_LIMITS(plt_iylow,plt_iyhigh)
	  else
	  endif
	endif
	ixlow  = plt_ixlow
	iylow  = plt_iylow
	ixhigh = plt_ixhigh
	iyhigh = plt_iyhigh
	ax     = float(ixlow-ixhigh) / (xlow-xhigh)
	bx     = float(ixlow) - ax * xlow
	ay     = float(iylow-iyhigh) / (ylow-yhigh)
	by     = float(iylow) - ay * ylow
	x_or   = plt_x_or
	y_or   = plt_y_or
	color  = plt_csto
	x      = xin * fplot
	y      = yin * fplot
	if (pen.eq.2) then
	  x2  = x + x_or
	  y2  = y + y_or
	  ix2 = int(ax*x2+bx)
	  iy2 = int(ay*y2+by)
	  ix1 = int(ax*x1+bx)
	  iy1 = int(ay*y1+by)
	  if (plunit.gt.0) then
	    if (.not.swrotateleft) then
	      call PS_Line(plunit,ix1       , iy1, ix2       , iy2, color)
	    else
	      call PS_Line(plunit,iyhigh-iy1, ix1, iyhigh-iy2, ix2, color)
	    endif
	  endif
	  x1 = x2
	  y1 = y2
	else if (pen.eq.3) then
	  x1 = x + x_or
	  y1 = y + y_or
	else if (pen.eq.-3) then
	  x_or     = x + x_or
	  y_or     = y + y_or
	  plt_x_or = x_or
	  plt_y_or = y_or
	  x1       = x
	  y1       = y
	else if (pen.eq.12) then
	  if (x.lt.x_low) x = x_low
	  if (y.lt.y_low) y = y_low
	  if (x.gt.x_high) x = x_high
	  if (y.gt.y_high) y = y_high
	  x2  = x + x_or
	  y2  = y + y_or
	  ix2 = int(ax*x2+bx)
	  iy2 = int(ay*y2+by)
	  ix1 = int(ax*x1+bx)
	  iy1 = int(ay*y1+by)
	  if (plunit.gt.0) then
	    if (.not.swrotateleft) then
	      call PS_Line(plunit,ix1       , iy1, ix2       , iy2, color)
	    else
	      call PS_Line(plunit,iyhigh-iy1, ix1, iyhigh-iy2, ix2, color)
	    endif
	  endif
	  x1 = x2
	  y1 = y2
	else if (pen.eq.13) then
	  if (x.lt.x_low) x = x_low
	  if (y.lt.y_low) y = y_low
	  if (x.gt.x_high) x = x_high
	  if (y.gt.y_high) y = y_high
	  x1 = x + x_or
	  y1 = y + y_or
	else if (pen.eq.999) then
	  init = 0
	  if (plunit.gt.0) then
	    call PS_epilogue(plunit)
	    close(unit=plunit)
	  endif
	  return
	else
	  write (*,'(a)') ' Unknown PEN parameter in PLOT call '
	endif
	return

	entry plots(fname,lunpl)
	plunit = lunpl
	open(unit=lunpl,file=fname)
	call PS_prologue(plunit)
	return

	entry bound(xl,xh,yl,yh)
	x_low = xl
	y_low = yl
	x_high = xh
	y_high = yh
	return

	entry factor(fval)
	fplot = fval
	return
c
c-- Newpen is slightly modified for PS use :
c-- call newpen(-2) results in a PS entry (2 setlinewidth)
c-- call newpen(1) still remains a color change of the 'pen'
c                  albeit through a pointer beyond the end of the 
c                  'regular' color table.
c
	entry newpen(color_in)
	if (color_in.lt.0) then
	  call PS_linewidth(abs(color_in))
	else
	  color    = max(min(color_in,n_total-n_colors),1) + n_colors
	  plt_csto = color
	endif
	return

	entry plot_convert(xin,yin,i_x1,i_y1)
	ixlow  = plt_ixlow
	iylow  = plt_iylow
	ixhigh = plt_ixhigh
	iyhigh = plt_iyhigh
	ax     = float(ixlow-ixhigh) / (xlow-xhigh)
	bx     = float(ixlow) - ax * xlow
	ay     = float(iylow-iyhigh) / (ylow-yhigh)
	by     = float(iylow) - ay * ylow
	x_or   = plt_x_or
	y_or   = plt_y_or
	csto   = plt_csto
	x      = xin * fplot
	y      = yin * fplot
	x2     = x + x_or
	y2     = y + y_or
	i_x1   = int(ax*x2+bx)
	i_y1   = int(ay*y2+by)
	return

	entry plot_add_offset(xin,yin)
	xin = xin + xlen
	yin = yin + ylen
	return

	entry plot_mouse_convert(ix_m,iy_m)
	if (swrotateleft) then
	i_xd = iy_m
c     i_yd = get_PS_WIDTH() - ix_m
	ix_m = i_xd
	iy_m = i_yd
	endif
	return

	entry plot_rotatedleft(swinquire)
	swinquire = .false.
	if (swrotateleft) swinquire = .true.
	return

	entry save_plot_context(id_context)
	i = 1
	do while (context(i).ne.0.and.i.le.num_context)
	  i = i + 1
	enddo
	if (i.le.num_context) then
	  id_context = i
	  context(i) = 1
	else
	  id_context = -1
	  return
	endif
c-- What about plt_nx, plt_ny, plt_csto
	pdsixlow(id_context)        = plt_ixlow
	pdsiylow(id_context)        = plt_iylow
	pdsixhigh(id_context)       = plt_ixhigh
	pdsiyhigh(id_context)       = plt_iyhigh
	pdsxlow(id_context)         = xlow
	pdsxhigh(id_context)        = xhigh
	pdsylow(id_context)         = ylow
	pdsyhigh(id_context)        = yhigh
	pdsxor(id_context)          = plt_x_or
	pdsyor(id_context)          = plt_y_or
	pdsfplot(id_context)        = fplot
	pdsxlen(id_context)         = xlen
	pdsylen(id_context)         = ylen
	pdspwidth(id_context)       = pwidth
	pdsyfrac(id_context)        = yfrac
	pdsyshft(id_context)        = yshft
	do k = 1,nvar
	  pdsparm(k,id_context)     = parm(k)
	enddo
	pdsswopen(id_context)       = swopen
	pdsswfrst(id_context)       = swfrst
	pdsswxlog(id_context)       = swxlog
	pdsswylog(id_context)       = swylog
	pdsswrotateleft(id_context) = swrotateleft
	pdsxmin(id_context)         = xmin
	pdsxmax(id_context)         = xmax
	pdsymin(id_context)         = ymin
	pdsymax(id_context)         = ymax
	pdsxfact(id_context)        = xfact
	pdsidset(id_context)        = idset
	pdsredy(id_context)         = redy
	pdsnxd(id_context)          = plt_nxd
	pdsnyd(id_context)          = plt_nyd
	pdsiloc(id_context)         = plt_iloc
	i_active_context                     = id_context
	return

	entry restore_plot_context(id_context)
	if (1.le.id_context.and.id_context.le.num_context) then
	  plt_ixlow        = pdsixlow(id_context)
	  plt_iylow        = pdsiylow(id_context)
	  plt_ixhigh       = pdsixhigh(id_context)
	  plt_iyhigh       = pdsiyhigh(id_context)
	  xlow             = pdsxlow(id_context)
	  xhigh            = pdsxhigh(id_context)
	  ylow             = pdsylow(id_context)
	  yhigh            = pdsyhigh(id_context)
	  plt_x_or         = pdsxor(id_context)
	  plt_y_or         = pdsyor(id_context)
	  fplot            = pdsfplot(id_context)
	  xlen             = pdsxlen(id_context)
	  ylen             = pdsylen(id_context)
	  pwidth           = pdspwidth(id_context)
	  yfrac            = pdsyfrac(id_context)
	  yshft            = pdsyshft(id_context)
	  do k = 1,nvar
	    parm(k)        = pdsparm(k,id_context)
	  enddo
	  swopen           = pdsswopen(id_context)
	  swfrst           = pdsswfrst(id_context)
	  swxlog           = pdsswxlog(id_context)
	  swylog           = pdsswylog(id_context)
	  swrotateleft     = pdsswrotateleft(id_context)
	  xmin             = pdsxmin(id_context)
	  xmax             = pdsxmax(id_context)
	  ymin             = pdsymin(id_context)
	  ymax             = pdsymax(id_context)
	  xfact            = pdsxfact(id_context)
	  idset            = pdsidset(id_context)
	  redy             = pdsredy(id_context)
	  plt_nxd          = pdsnxd(id_context)
	  plt_nyd          = pdsnyd(id_context)
	  plt_iloc         = pdsiloc(id_context)
	  i_active_context = id_context
	endif
	return

	entry delete_plot_context(id_context)
	if (1.le.id_context.and.id_context.le.num_context) then
	  context(id_context) = 0
	  id_context          = -1
	  i_active_context    = 0
	endif
	return

	entry get_active_plot_context_id(id_context)
	id_context = i_active_context
	return

	entry plot_rotate_left()
	swrotateleft = .true.
	return

	entry plot_rotate_none()
	swrotateleft = .false.
	return
	end

	logical function plot_translate(ix_mouse,iy_mouse,xmin,ymin,xmax,ymax,nx,ny,ix,iy,x_user,y_user)
	integer ix_mouse,iy_mouse,nx,ny,ix,iy
	real*4  x_user,y_user,xmin,ymin,xmax,ymax
c
c-- This function has two functions
c--  a) First determine whether a mouse position (ix,iy) is within 
c--     the plot frame limits
c--  b) If so, tranlate the ..,.. position pairs to user x,y such that 
c--     this represents a translated position in user coordinates.
c
c-- Inputs  : ix_mouse = mouse x in canvas
c--           iy_mouse = mouse y in canvas
c--           xmin     = user requested x minimum
c--           xmax     = user requested x maximum
c--           ymin     = user requested y minimum
c--           ymax     = user requested y maximum
c--           nx       = x division requested in plot
c--           ny       = y division requested in plot
c-- Outputs : ix       = x index in image for mouse position
c--           iy       = y index in image for mouse position
c--           x_user   = x position (mouse) in user coordinates
c--           y_user   = y position (mouse) in user coordinates
c
c-- Currently logarithmic axis are not properly accounted for (19980215)
c
	integer i_x1,i_y1,i_x2,i_y2,pltnx,pltny
	real*4  x,y,xr(2),yr(2),xfrac,yfrac
	logical swinquire
c
	plot_translate = .false.
	call plot_rotatedleft(swinquire)
	call plot_mouse_convert(ix_mouse,iy_mouse)
	x = 0.
	y = 0.
	call plot_convert(x,y,i_x1,i_y1)
	call plot_add_offset(x,y)
	call plot_convert(x,y,i_x2,i_y2)
c
	if (swinquire) then
	  ix = i_x1 - ix_mouse
	  iy = iy_mouse - i_y1
	  if ((1.le.ix.and.ix.le.(i_x1-i_x2)).and.(1.le.iy.and.iy.le.(i_y2-i_y1))) then
	    plot_translate = .true.
	    xr(1) = xmin
	    xr(2) = xmax
	    yr(1) = ymin
	    yr(2) = ymax
	    pltnx = nx
	    pltny = ny
	    call pkg_frame_actual_limits(pltnx,pltny,xr,yr)
	    xfrac = float(i_x1-ix_mouse) / float(i_x1-i_x2)
	    yfrac = float(iy_mouse-i_y1) / float(i_y2-i_y1)
	    if (pltnx.gt.0) then
	      x_user = xr(1) + xfrac*(xr(2)-xr(1))
	    else if (pltnx.lt.0) then
	      x_user = xr(2) * (10.** (float(pltnx) * (1.-xfrac)))
	    endif
	    if (pltny.gt.0) then
	      y_user = yr(1) + yfrac * (yr(2) - yr(1))
	    else if (pltny.lt.0) then
	      y_user = yr(2) * (10.** (float(pltny) * (1.-yfrac)))
	    endif
	  endif
	else
	  ix = ix_mouse - i_x1
	  iy = i_y1 - iy_mouse
	  if ((1.le.ix.and.ix.le.(i_x2-i_x1)).and.(1.le.iy.and.iy.le.(i_y1-i_y2))) then
	    plot_translate = .true.
	    xr(1) = xmin
	    xr(2) = xmax
	    yr(1) = ymin
	    yr(2) = ymax
	    pltnx = nx
	    pltny = ny
	    call pkg_frame_actual_limits(pltnx,pltny,xr,yr)
	    xfrac = float(ix_mouse-i_x1) / float(i_x2-i_x1)
	    yfrac = 1. - float(iy_mouse-i_y2) / float(i_y1-i_y2)
	    if (pltnx.gt.0) then
	      x_user = xr(1) + xfrac*(xr(2)-xr(1))
	    else if (pltnx.lt.0) then
	      x_user = xr(2) * (10.** (float(pltnx) * (1.-xfrac)))
	    endif
	    if (pltny.gt.0) then
	      y_user = yr(1) + yfrac * (yr(2) - yr(1))
	    else if (pltny.lt.0) then
	      y_user = yr(2) * (10.** (float(pltny) * (1.-yfrac)))
	    endif
	  endif
	endif
	return
	end

	subroutine PS_prologue(ilun)
	implicit none
	integer*4 ilun
c
	character*2  hexfmt(-128:127)
	character*24 cdate
	character*48 outhex
	character*56 outcol
	integer ix, iy, pointsize
	real*4  get_PS_SCALE
	integer i,j,i1,i2
c
	integer*4  n_colors, n_total
	parameter (n_colors=128,n_total=n_colors+6)
	byte       cmscolor(3,n_total)
	common /PS_color/cmscolor
c
	integer iflsh
	common /PSintern/iflsh
c
c-- Open PS flush file
c
	call get_lun(iflsh)
	open (unit=iflsh,file='psflsh.zzz',form='formatted',status='unknown')
c
	write (ilun ,'(''%!PS-Adobe-2.0'')')
	write (ilun ,'(''%%Creator: Fred Jansen ESTEC/SA S/W'')')
	write (ilun ,'(''%%Pages: 1'')')
	write (ilun ,'(''%%EndProlog'')')
	write (ilun ,'(''%%BeginSetup'')')
	write (ilun ,'(''%%/#copies 1 def'')')
	write (ilun ,'(''%%EndSetup'')')
	write (ilun ,'(''%%Page: '',i2,1x,i2)') 1,1
	call PS_prologue_page(0)
	return
	end

	subroutine PS_prologue_page(init)
	implicit none
c
	character*2  hexfmt(-128:127)
	character*24 cdate
	character*48 outhex
	character*56 outcol
	integer ix, iy, pointsize,init
	real*4  get_PS_SCALE
	integer i,j,i1,i2
c
	integer*4  n_colors, n_total
	parameter (n_colors=128,n_total=n_colors+6)
	byte       cmscolor(3,n_total)
	common /PS_color/cmscolor
c
	integer iflsh
	common /PSintern/iflsh
c
	pointsize = 12
	write (iflsh,'(''gsave'')')
	write (iflsh,'(''600 0 translate 90 rotate'')')
	write (iflsh,'(''/Helvetica findfont '',i3,'' scalefont setfont'')') pointsize
	call fdate(cdate)
	ix = 600 - 45 - pointsize
	iy = 60
c-- Temporarily edit out the dat stamp
c	write (iflsh,'(2i4,'' moveto ('',a,'') show'')') iy,ix,cdate
	write (iflsh,'(2f8.2,'' scale '')') get_PS_SCALE(), get_PS_SCALE()
	write (iflsh,'(''1 setlinewidth'')')
	write (iflsh,'(''true setstrokeadjust'')')
	write (iflsh,'(''/mt {moveto} bind def'')')
	write (iflsh,'(''/lt {lineto} bind def'')')
	write (iflsh,'(''/st {stroke} bind def'')')
	if (init.eq.0) call PS_init()
c
c-- Generate the color indexing table
c-- Write the 'hex' table
c
	do i = 0,255
	  j = i
	  if (i.ge.128) j = i - 256
	  write (hexfmt(j),'(z2)') i
	  if (hexfmt(j)(1:1).eq.' ') hexfmt(j)(1:1) = '0'
	enddo
c
c-- Write the 'indexed' color space !!
c
	do j = 1,len(outcol)
	  outcol(j:j) = ' '
	enddo
	i  = 1
	i2 = 0
	write (iflsh,'(''[/Indexed /DeviceRGB '',i3)') n_total - 1
	do i1 = 1,n_total
	  outcol(i  :i+1) = hexfmt(cmscolor(1,i1))
	  outcol(i+2:i+3) = hexfmt(cmscolor(2,i1))
	  outcol(i+4:i+5) = hexfmt(cmscolor(3,i1))
	  i = i + 7
	  if (i.gt.56) then
	    if (i2.eq.0) then
	      write (iflsh,'(a,a)') '<',outcol
	      i2 = 1
	    else
	      if (i1.eq.n_total) then
		write (iflsh,'(a,a)') outcol,'>'
	      else
		write (iflsh,'(a)') outcol
	      endif
	    endif
	    i = 1
	    do j = 1,len(outcol)
	      outcol(j:j) = ' '
	    enddo
	  endif
	enddo
	if (i.gt.1) then
	  write (iflsh,'(a,a)') outcol(1:i-1),'>'
	endif
	write (iflsh,'(''] setcolorspace'')')
	return
	end

	subroutine PS_epilogue(ilun)
	implicit none
	integer*4 ilun,lnblnk
	character*256 string
c
	integer iflsh
	common /PSintern/iflsh
c
	call PS_line_flush(ilun)
c
c-- Reinsert the flush file
c
	rewind (unit=iflsh)
	do while (.true.)
	  read (iflsh,'(a)',end=1) string
          write (ilun,'(a)') string(1:lnblnk(string))
	enddo
c
 1      write (ilun,'(''grestore'')')
	write (ilun,'(''showpage'')')
	write (ilun,'(''%%Trailer'')')
c
c-- Close PS flush file
c
	close (unit=iflsh,status='delete')
	call free_lun(iflsh)
c
	return
	end

	subroutine PS_line(ilun,ix1,iy1,ix2,iy2,color)
	implicit none
	integer*4 ilun,ix1,ix2,iy1,iy2,color,iline,iwdth
c
	integer*4  n_colors, n_total
	parameter (n_colors=128,n_total=n_colors+6)
	byte       cmscolor(3,n_total)
	common /PS_color/cmscolor
c
	logical nostroke, setwidth
	integer ix2_save,         iy2_save
	data    ix2_save/-99999/, iy2_save/-99999/, nostroke /.false./, setwidth/.false./
	save ix2_save, iy2_save, nostroke, setwidth, iwdth
c
	integer iflsh
	common /PSintern/iflsh
c
	if (ix1.eq.ix2_save.and.iy1.eq.iy2_save) then
	  write (iflsh,'(2i6, '' lt'')') ix2,iy2
	else
	  if (nostroke) then
	    write (iflsh,'(''st'')')
	    nostroke = .false.
	  endif
	  if (setwidth) then
	    write (iflsh,'(i5,'' setlinewidth'')') iwdth
	    setwidth = .false.
	  endif
c
	  write (iflsh,'(i4,'' setcolor '',2i5,'' mt '',2i5, '' lt'')') color-1,ix1,iy1,ix2,iy2
	  nostroke = .true.
	endif
	ix2_save = ix2
	iy2_save = iy2
	return
c
	entry PS_line_flush(ilun)
	if (nostroke) then
	  write (iflsh,'(''st'')')
	endif
	if (setwidth) then
	  write (iflsh,'(i5,'' setlinewidth'')') iwdth
	endif
	nostroke = .false.
	setwidth = .false.
	ix2_save = -99999
	iy2_save = -99999
	return
c
	entry PS_linewidth(iline)
	iwdth = iline
	setwidth = .true.
c
	return
	end

	subroutine PS_text(text,incolor)
	implicit none
c
c-- A mechanism needs to be added such that calling this subroutine when the PS 
c-- output is 'uninitialised' nothing happens, and the subroutine just returns 
c-- (possibly with a warning).
c
	character*(*) text
	integer incolor
c
	integer color,nout, ix, iy, iy_min, pointsize, lnblnk, j
	integer ix_in, iy_in, iy_min_in, pointsize_in, ix_ref, iy_ref, page_no
	logical swnewpage
	real*4  get_PS_SCALE
c
	integer iflsh
	common /PSintern/iflsh
c
	integer*4  n_colors, n_total
	parameter (n_colors=128,n_total=n_colors+6)
	byte       cmscolor(3,n_total)
	common /PS_color/cmscolor
c
	color = abs(incolor)
c
	nout = min(lnblnk(text),150)
	if (index(text(1:nout),')').ne.0) then
	  if (nout.lt.150) nout = nout + 1
	  j = index(text(1:nout),')')
	  if (text(j-1:j-1).ne.'\') then
	    text(j+1:nout) = text(j:nout-1)
	    text(j:j) = '\'
	  endif
	endif
	if (index(text(1:nout),'(').ne.0) then
	  if (nout.lt.150) nout = nout + 1
	  j = index(text(1:nout),'(')
	  if (text(j-1:j-1).ne.'\') then
	    text(j+1:nout) = text(j:nout-1)
	    text(j:j) = '\'
	  endif
	endif
c
	if (iy.ge.iy_min) then
	  if (incolor.ge.0) then
	    write (iflsh,'(i4,'' setcolor '',2i5,'' mt ('',a,'') show'')') color+n_colors,ix,iy,text(1:nout)
	  else
	    write (iflsh,'(''/Courier-Bold findfont '',i3,'' scalefont setfont'')') pointsize
	    write (iflsh,'(i4,'' setcolor '',2i5,'' mt ('',a,'') show'')') color+n_colors,ix,iy,text(1:nout)
	    write (iflsh,'(''/Courier findfont '',i3,'' scalefont setfont'')') pointsize
	  endif
	else
	  if (swnewpage) then
	    write (iflsh,'(''showpage'')')
c
c-- Repeat the entire sh'bang with gsave, colormap etc, for the next page.
c
	    page_no = page_no + 1
	    write (iflsh,'(''%%Page: '',i2,1x,i2)') page_no, page_no
	    call PS_prologue_page(1)
	    ix = ix_ref
	    iy = iy_ref
	    write (iflsh,'(''/Courier findfont '',i3,'' scalefont setfont'')') pointsize
	    write (iflsh,'(i4,'' setcolor '',2i5,'' mt ('',a,'') show'')') color+n_colors,ix,iy,text(1:nout)
	  endif
	endif
	iy = iy - pointsize
	return
c
	entry PS_text_init()
	pointsize = int(float(8) /get_PS_SCALE())
	iy        = int(float(540)/get_PS_SCALE())
	ix        = int(float(500)/get_PS_SCALE())
	iy_min    = int(float(30) /get_PS_SCALE())
	ix_ref    = ix
	iy_ref    = iy
	write (iflsh,'(''/Courier findfont '',i3,'' scalefont setfont'')') pointsize
	swnewpage = .false.
	page_no   = 1
	return
c
	entry PS_text_init_full(pointsize_in,ix_in,iy_in,iy_min_in)
	pointsize = int(float(pointsize_in) /get_PS_SCALE())
	iy        = int(float(iy_in)/get_PS_SCALE())
	ix        = int(float(ix_in)/get_PS_SCALE())
	iy_min    = int(float(iy_min_in) /get_PS_SCALE())
	ix_ref    = ix
	iy_ref    = iy
	write (iflsh,'(''/Courier findfont '',i3,'' scalefont setfont'')') pointsize
	swnewpage = .false.
	page_no   = 1
	return
c
	entry PS_text_allow_newpage()
	swnewpage = .true.
	return
	end

	subroutine PS_init()
	implicit none
	character*(*) prefix
c
	integer*4  n_colors, n_total, min, lnblnk
	parameter (n_colors=128,n_total=n_colors+6)
	byte       cmscolor(3,n_total)
	common /PS_color/cmscolor
c
	integer*4 ilun,i_table,i_request,i1,i2,i3,i,j,init
	character*200 filenm, path
c
	integer iflsh
	common /PSintern/iflsh
c
	data i_table/0/
	save i_table, path
	if (i_table.eq.0) i_table = 1
c
	if (0.le.i_table.and.i_table.le.15) then
	  write(filenm,'(a,z1,''.tab'')') path(1:lnblnk(path)),i_table
	endif
	call get_lun(ilun)
	open (unit=ilun,file=filenm,form='formatted',status='old')
	do i = 1,n_colors
	  read (ilun,*) i1,i2,i3
	  if (i1.gt.127) i1 = i1 - 256
	  if (i2.gt.127) i2 = i2 - 256
	  if (i3.gt.127) i3 = i3 - 256
	  cmscolor(1,i)      = i1
	  cmscolor(2,i)      = i2
	  cmscolor(3,i)      = i3
	  do j = 1,(256/n_colors) - 1
	    read (ilun,*) i1,i2,i3
	  enddo
	enddo
	close (unit=ilun)
	call free_lun(ilun)
c-- chartreuse
	cmscolor(1,n_colors+1) = '7f'x
	cmscolor(2,n_colors+1) = 'ff'x
	cmscolor(3,n_colors+1) = '00'x
c-- magenta
	cmscolor(1,n_colors+2) = 'ff'x
	cmscolor(2,n_colors+2) = '00'x
	cmscolor(3,n_colors+2) = 'ff'x
c-- spring green 1
	cmscolor(1,n_colors+3) = '00'x
	cmscolor(2,n_colors+3) = 'ff'x
	cmscolor(3,n_colors+3) = '7f'x
c-- salmon 1
	cmscolor(1,n_colors+4) = 'ff'x
	cmscolor(2,n_colors+4) = '8c'x
	cmscolor(3,n_colors+4) = '69'x
c-- orange red 1
	cmscolor(1,n_colors+5) = 'ff'x
	cmscolor(2,n_colors+5) = '45'x
	cmscolor(3,n_colors+5) = '00'x
c-- Open windows black (letters on B.G.)
	cmscolor(1,n_colors+6) = '00'x
	cmscolor(2,n_colors+6) = '00'x
	cmscolor(3,n_colors+6) = '00'x
c
	call PS_text_init()
	return
c
	entry PS_init_colourtable(i_request,prefix)
	i_table = i_request
	path    = prefix(1:min(lnblnk(prefix),len(path)))
	return
	end

	subroutine get_PS_WIDTH_LIMITS(ixlow, ixhigh)
	implicit none
	integer ixlow, ixhigh
	real*4  get_PS_SCALE
	logical swa4wide
	data swa4wide/.false./
c
	save swa4wide
c
	if (swa4wide) then
	  ixlow  = nint(float(60)  / get_PS_SCALE())
	  ixhigh = nint(float(800) / get_PS_SCALE())
	else
	  ixlow  = nint(float(60)  / get_PS_SCALE())
	  ixhigh = nint(float(550) / get_PS_SCALE())
	endif
c
	return
c
	entry set_PS_FULLPAGE()
	swa4wide = .true.
	return
c
	entry set_PS_NORMALPAGE()
	swa4wide = .false.
	return
	end

	subroutine get_PS_HEIGHT_LIMITS(iylow, iyhigh)
	implicit none
	integer iylow, iyhigh
	real*4 get_PS_SCALE
c
	iylow  = nint(float(60)  / get_PS_SCALE())
	iyhigh = nint(float(550) / get_PS_SCALE())
c
	return
	end

	real*4 function get_PS_SCALE()
	implicit none
	integer dpi
c
	dpi = 300
	get_PS_SCALE = float(72) / float(dpi)
c
	return
	end

c========================== S G P 4 Code ==============================
c========================== S G P 4 Code ==============================
c========================== S G P 4 Code ==============================
c========================== S G P 4 Code ==============================
c========================== S G P 4 Code ==============================
c========================== S G P 4 Code ==============================
c========================== S G P 4 Code ==============================
c========================== S G P 4 Code ==============================
c========================== S G P 4 Code ==============================
c========================== S G P 4 Code ==============================

*    this file contains a function to read two line element sets. while 
*    not formerly part of the sgp4 mathematical theory, it is 
*    required for practical implemenation.
*
*
*                           SUBROUTINE TWOLINE2RVSGP4
*
*  this function converts the two line element set character string data to
*    variables and initializes the sgp4 variables. several intermediate varaibles
*    and quantities are determined. note that the result is a "structure" so multiple
*    satellites can be processed simultaneously without having to reinitialize. the
*    verification mode is an important option that permits quick checks of any
*    changes to the underlying technical theory. this option works using a
*    modified tle file in which the start, stop, and delta time values are
*    included at the end of the second line of data. this only works with the
*    verification mode. the catalog mode simply propagates from -1440 to 1440 min
*    from epoch and is useful when performing entire catalog runs.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs        :
*    Numsats     - Number of satellites processed. It also becomes the record
*                  number for each satellite
*    typerun     - type of run                    verification 'V', catalog 'C', 
*                                                 manual 'M'
*    typeinput   - type of manual input           mfe 'M', epoch 'E', dayofyr 'D'
*    whichconst  - which set of constants to use  72, 84
*    opsmode   - type of manual input           afspc 'a', imporved 'i'
*
*  outputs       :
*    Code        - EOF indicator. Code = 999 when EOF reached
*
*  coupling      :
*    days2mdhms  - conversion of days to month, day, hour, minute, second
*    jday        - convert day month year hour minute second into julian date
*    sgp4init    - initialize the sgp4 variables
*
*  Files         :
*    Unit 10     - test.elm        input 2-line element set file
*    Unit 11     - test.bak        output file
*    Unit 15     - sgp4rec.bak     temporary file of record for 2 line element sets
*
*  references    :
*    norad spacetrack report #3
*    vallado, crawford, hujsak, kelso  2006
*------------------------------------------------------------------------------

      SUBROUTINE TwoLine2RVSGP4 ( NumSats, Typerun, typeinput, 
     &                            whichconst, Code )
        IMPLICIT NONE
        Character Typerun, typeinput
        Integer Code, NumSats, whichconst
        REAL*8 startmfe, stopmfe, deltamin

* ----------------------------  Locals  -------------------------------
        REAL*8 J2, mu, RadiusEarthKm,VKmPerSec, xke, tumin
        REAL*8 BC,EPDay, sec, xpdotp, j3, j4, j3oj2 
        REAL*8 startsec, stopsec, startdayofyr, stopdayofyr, jdstart, 
     &         jdstop
        INTEGER startyear, stopyear, startmon, stopmon, startday, 
     &          stopday, starthr, stophr, startmin, stopmin 
        INTEGER Yr,Mon,Day,Hr,Minute,  ICrdno,nexp,bexp, error
        CHARACTER Show
        Character*130 LongStr1,LongStr2

        COMMON /DebugHelp/ Help
        CHARACTER Help

        INCLUDE 'sgp4.cmn'
        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        Show = 'N'
        xpdotp        =  1440.0D0 / (2.0D0 * pi) ! 229.1831180523293

        CALL getgravconst( whichconst, tumin, mu, radiusearthkm, xke,  
     &       j2, j3, j4, j3oj2 );
        VKmPerSec     =  RadiusEarthKm * xke / 60.0D0

c        make sure the main program opens this file, otherwise do so here
c        ! store results in a temporary file of record
c        OPEN(115,FILE='Sgp4Rec.bak',ACCESS='DIRECT', FORM='UNFORMATTED',
c     &       RECL=1000,STATUS='UNKNOWN')

* ----------------- READ THE FIRST LINE OF ELEMENT SET ----------------
        Code = 0

        LongStr1 = ' '
   50   READ(110,'(a130)',END=999) LongStr1
        IF(LongStr1(1:1) .eq. '#') GOTO 50 ! Commented line of text, skip

        READ(LongStr1,500) ICRDNO,SatNum,SatName,EpochYr,EpDay,
     &                       NDot,NDDot,nexp,BStar,bexp,EPHTYP,ELNO
  500   FORMAT( I1,1X,I5,1X,A10,I2,D12.0,1X,D10.0,1X,
     &          F6.5,I2,1X,F6.5,I2,1X,I1,1X,I4 )

* ----------- READ THE SECOND LINE OF ELEMENT SET AND TIME ------------
        LongStr2 = ' '
   51   READ(110,'(a130)',END=999) LongStr2
        IF(LongStr2(1:1) .eq. '#') GOTO 51 ! Commented line of text, skip

        IF (Typerun.eq.'V') THEN
          READ(LongStr2,502) ICRDNO,Inclo,nodeo,Ecco,Argpo,Mo,No,REVI,
     &              startmfe, stopmfe, DeltaMin
         else
          READ(LongStr2,501) ICRDNO,Inclo,nodeo,Ecco,Argpo,Mo,No,REVI
         endif
  501   FORMAT( I1,7X,D8.0,1X,D8.0,1X,F7.7,1X,D8.0,1X,D8.0,1X,D11.0,I5)
  502   FORMAT( I1,7X,D8.0,1X,D8.0,1X,F7.7,1X,D8.0,1X,D8.0,1X,D11.0,I5,
     &          1X,F12.6,F12.6,F12.6 )

* ---------------------- CONVERT TO INTERNAL UNITS --------------------
* ---- RADIANS, DISTANCE IN EARTH RADII, AND VELOCITY IN ER/KEMIN) ----
        NDDot  = NDDot * 10.0D0**Nexp
        NDot   = NDot / (XPDOTP*1440)
        NDDot  = NDDot / (XPDOTP*1440*1440)
        BStar  = BStar * 10.0D0**Bexp

        No     = No / XPDOTP
        a      = (No*TUMin)**(-2.0D0/3.0D0)
        Inclo  = Inclo  * Deg2Rad
        nodeo  = nodeo * Deg2Rad
        Argpo  = Argpo * Deg2Rad
        Mo     = Mo   * Deg2Rad
                                                                        
        IF (DABS(Ecco-1.0D0) .gt. 0.000001D0) THEN
            Altp= (a*(1.0D0-Ecco))-1.0D0
            Alta= (a*(1.0D0+Ecco))-1.0D0
          ELSE
            Alta= 999999.9D0
            Altp= 2.0D0* (4.0D0/(No*No)**(1.0D0/3.0D0))
          ENDIF

        ! ---- Ballistic Coefficient ----
        IF (DABS(BStar) .gt. 0.00000001D0) THEN
            BC= 1.0D0/(12.741621D0*BStar)
          ELSE
            BC= 1.111111111111111D0
          ENDIF

        ! ----------------------------------------------------------------
        ! find sgp4epoch time of element set
        ! remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
        ! and minutes from the epoch (time)
        ! ----------------------------------------------------------------

        ! Temporary year fix
        IF (EpochYr.lt.57) THEN
            Yr = EpochYr + 2000
          ELSE
            Yr = EpochYr + 1900
          ENDIF

        CALL Days2MDHMS( Yr,EpDay, Mon,Day,Hr,Minute,Sec )
        CALL JDAY ( Yr,Mon,Day,Hr,Minute,Sec,  JDSatEpoch )

* ------------------- MAKE INITIAL PREDICTION AT EPOCH ----------------
        ! 2433281.5 - 2400000.5 = 33281.0, thus time from 1950
        CALL SGP4Init( whichconst,
     &                 SatNum,BStar, Ecco, JDSatEpoch-2433281.5D0,
     &                 Argpo,Inclo,Mo,No, nodeo, Error )

        ! ---- Write common block of data into file of record ----
        WRITE(115,Rec=NumSats) SatName,
     &          SatNum, ELNO  , EPHTYP, REVI  , EpochYr,
     &          BStar , Ecco  , Inclo , nodeo, Argpo , No    , Mo    ,
     &          NDot  , NDDot ,
     &          alta  , altp  , a     ,
     &          DeltaMin, JDSatEpoch, EpochDays,
     &          Isimp , Init  , Method, Opsmode,
     &          Aycof , CON41 , Cc1   , Cc4   , Cc5   , D2    , D3    ,
     &          D4    , Delmo , Eta   , ArgpDot,Omgcof, Sinmao,
     &          T2cof , T3cof , T4cof , T5cof , X1mth2, X7thm1, MDot  ,
     &          nodeDot,Xlcof, Xmcof , Xnodcf,
     &          D2201 , D2211 , D3210 , D3222 , D4410 , D4422 , D5220 ,
     &          D5232 , D5421 , D5433 , Dedt  , Del1  , Del2  , Del3  ,
     &          Didt  , Dmdt  , Dnodt , Domdt , E3    , Ee2   , Peo   ,
     &          Pgho  , Pho   , Pinco , Plo   , Se2   , Se3   , Sgh2  ,
     &          Sgh3  , Sgh4  , Sh2   , Sh3   , Si2   , Si3   , Sl2   ,
     &          Sl3   , Sl4   , GSTo  , Xfact , Xgh2  , Xgh3  , Xgh4  ,
     &          Xh2   , Xh3   , Xi2   , Xi3   , Xl2   , Xl3   , Xl4   ,
     &          Xlamo , Zmol  , Zmos  , Atime , Xli   , Xni   , IRez

        IF(Error .GT. 0) THEN
            WRITE( *,*) '# *** SGP4 Model Error ***',Error
          ENDIF

c      write tle output details
c      INCLUDE 'debug8.for'

        ! ---- Fix to indicate end-of-file
        GOTO 1000
  999   Code = 999
 1000   CONTINUE

       RETURN
       END  !       SUBROUTINE TwoLine2RVSGP4




c========================================================================
c========================================================================
c========================================================================
c========================================================================
c========================================================================

*     ----------------------------------------------------------------
*
*                               sgp4ext.for
*
*    this file contains extra routines needed for the main test program for sgp4.
*    these routines are derived from the astro libraries.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                    2007
*                              by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*               2 apr 07  david vallado
*                           misc updates for new baseline
*    changes :
*              14 aug 06  david vallado
*                           original baseline
*       ----------------------------------------------------------------



* ------------------------------------------------------------------------------
*
*                           SUBROUTINE MAG
*
*  This subroutine finds the magnitude of a vector.  The tolerance is set to
*    0.00000001D0, thus the 1.0D0E-16 for the squared test of underflows.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec       - Vector
*
*  OutPuts       :
*    Vec       - Answer stored in fourth component
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------  

      REAL*8 FUNCTION MAG    ( Vec )
        IMPLICIT NONE
        REAL*8 Vec(3)
* -----------------------------  Locals  ------------------------------
        Real*8 Temp

        ! --------------------  Implementation   ----------------------
        Temp= Vec(1)*Vec(1) + Vec(2)*Vec(2) + Vec(3)*Vec(3)

        IF ( DABS( Temp ) .ge. 1.0D-16 ) THEN
            MAG = DSQRT( Temp )
          ELSE
            MAG = 0.0D0
          ENDIF
      RETURN
      END  ! end mag


* ------------------------------------------------------------------------------
*
*                           SUBROUTINE CROSS
*
*  This subroutine crosses two vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec1        - Vector number 1
*    Vec2        - Vector number 2
*
*  OutPuts       :
*    OutVec      - Vector result of A x B
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------  

      SUBROUTINE CROSS       ( Vec1,Vec2, OutVec )
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3), OutVec(3)

        ! --------------------  Implementation   ----------------------
        OutVec(1)= Vec1(2)*Vec2(3) - Vec1(3)*Vec2(2)
        OutVec(2)= Vec1(3)*Vec2(1) - Vec1(1)*Vec2(3)
        OutVec(3)= Vec1(1)*Vec2(2) - Vec1(2)*Vec2(1)

      RETURN
      END  ! end cross


* ------------------------------------------------------------------------------
*
*                           FUNCTION DOT
*
*  This function finds the DOT product of two vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec1        - Vector number 1
*    Vec2        - Vector number 2
*
*  OutPuts       :
*    DOT         - Result
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------  

      REAL*8 FUNCTION DOT    ( Vec1,Vec2 )
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3)

        ! --------------------  Implementation   ----------------------
        DOT= Vec1(1)*Vec2(1) + Vec1(2)*Vec2(2) + Vec1(3)*Vec2(3)
      RETURN
      END  ! end dot


* ------------------------------------------------------------------------------
*
*                           SUBROUTINE ANGLE
*
*  This subroutine calculates the ANGLE between two vectors.  The output is
*    set to 999999.1D0 to indicate an undefined value.  Be SURE to check  
*    this at the output phase.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Vec1        - Vector number 1
*    Vec2        - Vector number 2
*
*  OutPuts       :
*    Theta       - ANGLE between the two vectors  -Pi to Pi
*
*  Locals        :
*    Temp        - Temporary REAL variable
*
*  Coupling      :
*    DOT           DOT Product of two vectors
*    DACOS         Arc Cosine FUNCTION
*
* ------------------------------------------------------------------------------  

      SUBROUTINE ANGLE       ( Vec1,Vec2, Theta )
        IMPLICIT NONE
        REAL*8 Vec1(3), Vec2(3), Theta, magvec1, magvec2
        EXTERNAL Dot, Mag
* -----------------------------  Locals  ------------------------------
        REAL*8 Temp, Dot, Mag
        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        magvec1 = MAG(vec1)
        magvec2 = MAG(vec2)
        IF ( magVec1*magVec2 .gt. Small**2 ) THEN
            Temp= DOT(Vec1,Vec2) / (magVec1*magVec2)
            IF ( DABS( Temp ) .gt. 1.0D0 ) THEN
                Temp= DSIGN(1.0D0, Temp)
              ENDIF
            Theta= DACOS( Temp ) 
          ELSE
            Theta= Undefined
          ENDIF
      RETURN
      END  ! end angle


* ------------------------------------------------------------------------------
*
*                           FUNCTION ASINH
*
*  This function evaluates the inverse hyperbolic sine.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    XVal        - ANGLE Value                                  any real
*
*  OutPuts       :
*    ASINH       - Result                                       any real
*
*  Locals        :
*    None.
*
*  Coupling      :
*    None.
*
* ------------------------------------------------------------------------------  

      REAL*8 FUNCTION ASINH( XVal )
        IMPLICIT NONE
        REAL*8 XVal

        ! --------------------  Implementation   ----------------------
        ASINH= DLOG( XVal + DSQRT( XVal*XVal + 1.0D0 ) )

      RETURN
      END  ! end asinh


* ------------------------------------------------------------------------------
*
*                           SUBROUTINE NEWTONNU
*
*  This subroutine solves Keplers equation when the true anomaly is known.
*    The Mean and Eccentric, parabolic, or hyperbolic anomaly is also found.
*    The parabolic limit at 168ø is arbitrary. The hyperbolic anomaly is also
*    limited. The hyperbolic sine is used because it's not double valued.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Ecc         - Eccentricity                   0.0D0 to
*    Nu          - True Anomaly                   -2Pi to 2Pi rad
*
*  Outputs       :
*    E0          - Eccentric Anomaly              0.0D0 to 2Pi rad       153.02 deg
*    M           - Mean Anomaly                   0.0D0 to 2Pi rad       151.7425 deg 
*
*  Locals        :
*    E1          - Eccentric Anomaly, next value  rad
*    SinE        - Sine of E
*    CosE        - Cosine of E
*    Ktr         - Index
*
*  Coupling      :
*    ASINH       - Arc hyperbolic sine
*    SINH        - Hyperbolic Sine
*
*  References    :
*    Vallado       2007, 85, Alg 5
*
* ------------------------------------------------------------------------------

      SUBROUTINE NEWTONNU    ( Ecc, Nu, E0, M )
        IMPLICIT NONE
        REAL*8 Ecc, Nu, E0, M
        EXTERNAL ASINH
* -----------------------------  Locals  ------------------------------
        REAL*8 SinE, CosE, ASINH

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        E0= 999999.9D0
        M = 999999.9D0
        ! --------------------------- Circular ------------------------
        IF ( DABS( Ecc ) .lt. 0.000001D0 ) THEN
            M = Nu
            E0= Nu 
          ELSE
            ! ---------------------- Elliptical -----------------------
            IF ( Ecc .lt. 0.999D0 ) THEN
                SinE= ( DSQRT( 1.0D0-Ecc*Ecc ) * DSIN(Nu) ) /
     &                ( 1.0D0+Ecc*DCOS(Nu) )
                CosE= ( Ecc + DCOS(Nu) ) / ( 1.0D0 + Ecc*DCOS(Nu) )
                E0  = DATAN2( SinE, CosE )
                M   = E0 - Ecc*DSIN(E0) 
              ELSE
                ! -------------------- Hyperbolic  --------------------
                IF ( Ecc .gt. 1.0001D0 ) THEN
                    IF ( ((Ecc .gt. 1.0D0) .and. (DABS(Nu)+0.00001D0
     &                     .lt. Pi-DACOS(1.0D0/Ecc)) ) ) THEN
                        SinE= ( DSQRT( Ecc*Ecc-1.0D0 ) * DSIN(Nu) ) /
     &                        ( 1.0D0 + Ecc*DCOS(Nu) )
                        E0  = ASINH( SinE )
                        M   = Ecc*DSINH(E0) - E0
                      ENDIF 
                  ELSE
                    ! ----------------- Parabolic ---------------------
                    IF ( DABS(Nu) .lt. 168.0D0/57.29578D0 ) THEN
                        E0= DTAN( Nu*0.5D0 )
                        M = E0 + (E0*E0*E0)/3.0D0 
                      ENDIF
                  ENDIF
              ENDIF
          ENDIF

        IF ( Ecc .lt. 1.0D0 ) THEN
            M = DMOD( M, 2.0D0*Pi )
            IF ( M .lt. 0.0D0 ) THEN
                M= M + 2.0D0*Pi 
              ENDIF
            E0 = DMOD( E0, 2.0D0*Pi )
          ENDIF 
      RETURN
      END  ! end newtonnu


* ------------------------------------------------------------------------------
*
*                           SUBROUTINE rv2coe
*
*  This subroutine finds the classical orbital elements given the Geocentric
*    Equatorial Position and Velocity vectors.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    R           - IJK Position vector            km
*    V           - IJK Velocity vector            km / s
*    mu          - gravitational parameter        km3 / s2
*
*  Outputs       :
*    P           - SemiLatus rectum               km
*    A           - semimajor axis                 km
*    Ecc         - Eccentricity
*    Incl        - inclination                    0.0D0 to Pi rad
*    Omega       - Longitude of Ascending Node    0.0D0 to 2Pi rad
*    Argp        - Argument of Perigee            0.0D0 to 2Pi rad
*    Nu          - True anomaly                   0.0D0 to 2Pi rad
*    M           - Mean anomaly                   0.0D0 to 2Pi rad
*    ArgLat      - Argument of Latitude      (CI) 0.0D0 to 2Pi rad
*    LamTrue     - True Longitude            (CE) 0.0D0 to 2Pi rad
*    LonPer      - Longitude of Periapsis    (EE) 0.0D0 to 2Pi rad
*
*  Locals        :
*    HBar        - Angular Momentum H Vector      km2 / s
*    EBar        - Eccentricity     E Vector
*    NBar        - Line of Nodes    N Vector
*    c1          - V**2 - u/R
*    RDotV       - R DOT V
*    Hk          - Hk norm vector
*    SME         - Specfic Mechanical Energy      km2 / s2
*    i           - index
*    E           - Eccentric, Parabolic,
*                  Hyperbolic Anomaly             rad
*    Temp        - Temporary variable
*    TypeOrbit   - Type of orbit                  EE, EI, CE, CI
*
*  Coupling      :
*    MAG         - Magnitude of a vector
*    CROSS       - CROSS product of two vectors
*    DOT         - DOT product of two vectors
*    ANGLE       - Find the ANGLE between two vectors
*    NEWTONNU    - Find the mean anomaly
*
*  References    :
*    Vallado       2007, 121, Alg 9, Ex 2-5
*
* ------------------------------------------------------------------------------

      SUBROUTINE rv2coe      ( R, V, mu, P, A, Ecc, Incl, Omega, Argp,
     &                         Nu, M, ArgLat, TrueLon, LonPer )
        IMPLICIT NONE
        REAL*8 R(3), V(3), mu, P, A, Ecc, Incl, Omega, Argp, Nu, M, 
     &         ArgLat, TrueLon, LonPer
        EXTERNAL DOT, MAG
* -----------------------------  Locals  ------------------------------
        REAL*8 c1, RDotV, hk, SME, Hbar(3), Ebar(3), Nbar(3),
     &         Dot, E, Temp, MAG, maghbar, magnbar, magr, magv
        INTEGER i
        CHARACTER*2 TypeOrbit

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------
        magr = MAG( R )
        magv = MAG( V )
        ! ------------------  Find H N and E vectors   ----------------
        CALL CROSS( R, V, HBar )
        maghbar = MAG(Hbar)
        IF ( maghbar .gt. Small ) THEN
            NBar(1)= -HBar(2)
            NBar(2)=  HBar(1)
            NBar(3)=   0.0D0
            magnbar = MAG( Nbar )
            c1 = magv**2 - mu/magr
            RDotV= DOT( R, V )
            DO i= 1 , 3
                EBar(i)= (c1*R(i) - RDotV*V(i))/mu
              ENDDO

            Ecc = MAG( EBar )

            ! ------------  Find a e and semi-Latus rectum   ----------
            SME= ( magv*magv*0.5D0 ) - ( mu/magr )
            IF ( DABS( SME ) .gt. Small ) THEN
                A= -mu / (2.0D0*SME)
              ELSE
                A= Infinite
              ENDIF
            P = maghbar*maghbar/mu

            ! -----------------  Find inclination   -------------------
            Hk= HBar(3)/maghbar
c            IF ( DABS( DABS(Hk) - 1.0D0 ) .lt. Small ) THEN
c                ! -------------  Equatorial Orbits   ------------------
c                IF ( DABS(HBar(3)) .gt. 0.0D0 ) THEN
c                    Hk= DSIGN(1.0D0, HBar(3))
c                  ENDIF
c              ENDIF
            Incl= DACOS( Hk ) 

            ! --------  Determine type of orbit for Later use  --------
            ! ------ Elliptical, Parabolic, Hyperbolic Inclined -------
            TypeOrbit= 'EI' 
            IF ( Ecc .lt. Small ) THEN
                ! ----------------  Circular Equatorial ---------------
                IF ( (Incl.lt.Small).or.(DABS(Incl-Pi).lt.Small) ) THEN
                    TypeOrbit= 'CE'
                  ELSE
                    ! --------------  Circular Inclined ---------------
                    TypeOrbit= 'CI'
                  ENDIF
              ELSE
                ! - Elliptical, Parabolic, Hyperbolic Equatorial --
                IF ( (Incl.lt.Small).or.(DABS(Incl-Pi).lt.Small) ) THEN
                    TypeOrbit= 'EE'
                  ENDIF
              ENDIF

            ! ----------  Find Longitude of Ascending Node ------------
            IF ( magnbar .gt. Small ) THEN
                Temp= NBar(1) / magnbar
                IF ( DABS(Temp) .gt. 1.0D0 ) THEN
                    Temp= DSIGN(1.0D0, Temp)
                  ENDIF
                Omega= DACOS( Temp ) 
                IF ( NBar(2) .lt. 0.0D0 ) THEN
                    Omega= TwoPi - Omega
                  ENDIF
              ELSE
                Omega= Undefined 
              ENDIF

            ! ---------------- Find Argument of perigee ---------------
            IF ( TypeOrbit .eq. 'EI' ) THEN
                CALL ANGLE( NBar, EBar, Argp )
                IF ( EBar(3) .lt. 0.0D0 ) THEN
                    Argp= TwoPi - Argp 
                  ENDIF
              ELSE
                Argp= Undefined 
              ENDIF

            ! ------------  Find True Anomaly at Epoch    -------------
            IF ( TypeOrbit(1:1) .eq. 'E' ) THEN
                CALL ANGLE( EBar, r, Nu )
                IF ( RDotV .lt. 0.0D0 ) THEN
                    Nu= TwoPi - Nu 
                  ENDIF
              ELSE
                Nu= Undefined 
              ENDIF

            ! ----  Find Argument of Latitude - Circular Inclined -----
            IF ( TypeOrbit .eq. 'CI' ) THEN
                CALL ANGLE( NBar, R, ArgLat )
                IF ( R(3) .lt. 0.0D0 ) THEN
                    ArgLat= TwoPi - ArgLat
                  ENDIF
              ELSE
                ArgLat= Undefined 
              ENDIF

            ! -- Find Longitude of Perigee - Elliptical Equatorial ----
            IF ( ( Ecc.gt.Small ) .and. (TypeOrbit.eq.'EE') ) THEN
                Temp= EBar(1)/Ecc
                IF ( DABS(Temp) .gt. 1.0D0 ) THEN
                    Temp= DSIGN(1.0D0, Temp)
                  ENDIF
                LonPer= DACOS( Temp ) 
                IF ( EBar(2) .lt. 0.0D0 ) THEN
                    LonPer= TwoPi - LonPer 
                  ENDIF
                IF ( Incl .gt. HalfPi ) THEN
                    LonPer= TwoPi - LonPer
                  ENDIF
              ELSE
                LonPer= Undefined
              ENDIF

            ! -------- Find True Longitude - Circular Equatorial ------
            IF ( ( magr.gt.Small ) .and. ( TypeOrbit.eq.'CE' ) ) THEN
                Temp= R(1)/magr
                IF ( DABS(Temp) .gt. 1.0D0 ) THEN
                    Temp= DSIGN(1.0D0, Temp)
                  ENDIF
                TrueLon= DACOS( Temp )
                IF ( R(2) .lt. 0.0D0 ) THEN
                    TrueLon= TwoPi - TrueLon
                  ENDIF
                IF ( Incl .gt. HalfPi ) THEN
                    TrueLon= TwoPi - TrueLon
                  ENDIF
              ELSE
                TrueLon= Undefined
              ENDIF

            ! ------------ Find Mean Anomaly for all orbits -----------
            CALL NEWTONNU(Ecc, Nu, E, M )

         ELSE
           P    = Undefined
           A    = Undefined
           Ecc  = Undefined
           Incl = Undefined
           Omega= Undefined 
           Argp = Undefined 
           Nu   = Undefined 
           M    = Undefined 
           ArgLat  = Undefined 
           TrueLon= Undefined 
           LonPer = Undefined 
         ENDIF 

      RETURN
      END  ! end rv2coe


* -----------------------------------------------------------------------------
*
*                           SUBROUTINE JDay
*
*  This subroutine finds the Julian date given the Year, Month, Day, and Time.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Year        - Year                           1900 .. 2100
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Universal Time Hour            0 .. 23
*    Min         - Universal Time Min             0 .. 59
*    Sec         - Universal Time Sec             0.0D0 .. 59.999D0
*    WhichType   - Julian .or. Gregorian calender   'J' .or. 'G'
*
*  Outputs       :
*    JD          - Julian Date                    days from 4713 BC
*
*  Locals        :
*    B           - Var to aid Gregorian dates
*
*  Coupling      :
*    None.
*
*  References    :
*    Vallado       2007, 189, Alg 14, Ex 3-14
* -----------------------------------------------------------------------------

      SUBROUTINE JDay        ( Year,Mon,Day,Hr,Min, Sec, JD )
        IMPLICIT NONE
        INTEGER Year, Mon, Day, Hr, Min
        REAL*8  Sec, JD

        ! --------------------  Implementation   ----------------------
        JD= 367.0D0 * Year
     &        - INT( (7* (Year+INT ( (Mon+9)/12.0) ) ) * 0.25D0 )
     &        + INT( 275*Mon / 9.0 )
     &        + Day + 1721013.5D0
     &        + ( (Sec/60.0D0 + Min ) / 60.0D0 + Hr ) / 24.0D0
*     &      - 0.5D0*DSIGN(1.0D0, 100.0D0*Year + Mon - 190002.5D0) + 0.5D0
      RETURN
      END  ! end jday


* -----------------------------------------------------------------------------
*
*                           SUBROUTINE DAYS2MDHMS
*
*  This subroutine converts the day of the year, days, to the equivalent month
*    day, hour, Minute and second.
*
*  Algorithm     : Set up array for the Number of days per month
*                  Find Leap Year - be sure to account for the 400 years
*                  Loop through a Temp value for WHILE the value is .lt. the days
*                  Perform INTEGER conversions to the correct day and month
*                  Convert remainder into H M S using type conversions
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    Year        - Year                          +1900 .. 2100+
*    Days        - Julian Day of the year         0.0D0  .. 366.0D0
*
*  OutPuts       :
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Hour                           0 .. 23
*    Min         - Minute                         0 .. 59
*    Sec         - Second                         0.0D0 .. 59.999D0
*
*  Locals        :
*    DayofYr     - Day of year
*    Temp        - Temporary REAL*8 values
*    IntTemp     - Temporary INTEGER value
*    i           - Index
*    LMonth[12]  - INTEGER Array containing the Number of days per month
*
*  Coupling      :
*    None.
* -----------------------------------------------------------------------------

      SUBROUTINE DAYS2MDHMS  ( Year,Days,  Mon,Day,Hr,Min,Sec )
        IMPLICIT NONE
        REAL*8 Days,Sec
        INTEGER Year, Mon, Day, Hr, Min
* ----------------------------  Locals  -------------------------------
        INTEGER IntTemp,i,DayofYr, LMonth(12)
        REAL*8 Temp

        ! --------------------  Implementation   ----------------------
        ! -------------- Set up array of days in month  ---------------
        DO i = 1,12
            LMonth(i) = 31
          ENDDO
        LMonth( 2) = 28
        LMonth( 4) = 30
        LMonth( 6) = 30
        LMonth( 9) = 30
        LMonth(11) = 30

        DayofYr= IDINT(Days )

        ! ---------------- Find month and Day of month ----------------
        IF (MOD(Year,4).eq.0) THEN
            LMonth(2)= 29
          ENDIF
        i= 1
        IntTemp= 0
        DO WHILE ( (DayofYr.gt.IntTemp + LMonth(i) ) .and. ( i.lt.12 ))
            IntTemp= IntTemp + LMonth(i)
            i= i+1
          ENDDO
        Mon= i
        Day= DayofYr - IntTemp

        ! ---------------- Find hours Minutes and seconds -------------
        Temp= (Days - DayofYr )*24.0D0
        Hr  = IDINT( Temp )
        Temp= (Temp-Hr) * 60.0D0
        Min = IDINT( Temp )
        Sec = (Temp-Min) * 60.0D0

        ! ---- Check for roundoff errors
c        IF (Sec .ge. 59.9999D0) THEN
c            Sec = 0.0D0
c            Min = Min + 1
c            IF (Min .gt. 59) THEN
c                Min = 0
c                Hr = Hr + 1
c                IF (Hr .gt. 23) THEN
c                    Hr = 0
c                    Day = Day + 1
c                  ENDIF
c              ENDIF
c          ENDIF
      RETURN
      END  ! end days2mdhms


* -----------------------------------------------------------------------------
*
*                           SUBROUTINE INVJDay
*
*  This subroutine finds the Year, month, day, hour, Minute and second
*  given the Julian date. TU can be UT1, TDT, TDB, etc.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date                    days from 4713 BC
*
*  OutPuts       :
*    Year        - Year                           1900 .. 2100
*    Mon         - Month                          1 .. 12
*    Day         - Day                            1 .. 28,29,30,31
*    Hr          - Hour                           0 .. 23
*    Min         - Minute                         0 .. 59
*    Sec         - Second                         0.0D0 .. 59.999D0
*
*  Locals        :
*    Days        - Day of year plus fractional
*                  portion of a day               days
*    Tu          - Julian Centuries from 0 h
*                  Jan 0, 1900
*    Temp        - Temporary real values
*    LeapYrs     - Number of Leap years from 1900
*
*  Coupling      :
*    DAYS2MDHMS  - Finds MD HMS given Days and Year
*
*  References    :
*    Vallado       2007, 208, Alg 22, Ex 3-13
* -----------------------------------------------------------------------------

      SUBROUTINE INVJDay     ( JD, Year,Mon,Day,Hr,Min, Sec )
        IMPLICIT NONE
        INTEGER Year, Mon, Day, Hr, Min
        REAL*8  Sec, JD
* ----------------------------  Locals  -------------------------------
        INTEGER LeapYrs
        REAL*8  Days, Tu, Temp

        ! --------------------  Implementation   ----------------------
        ! ---------------- Find Year and Days of the year -------------
        Temp   = JD-2415019.5D0
        Tu     = Temp / 365.25D0
        Year   = 1900 + IDINT( Tu )
        LeapYrs= IDINT( ( Year-1901 )*0.25D0 )
        Days   = Temp - ((Year-1900)*365.0D0 + LeapYrs )

        ! -------------- Check for case of beginning of a year --------
        IF ( Days .lt. 1.0D0 ) THEN
            Year   = Year - 1
            LeapYrs= IDINT( ( Year-1901 )*0.25D0 )
            Days   = Temp - ((Year-1900)*365.0D0 + LeapYrs )
          ENDIF

        ! ------------------ Find remaing data  -----------------------
        CALL DAYS2MDHMS( Year,Days, Mon,Day,Hr,Min,Sec )

      RETURN
      END

c=======================================================================
c=======================================================================
c=======================================================================
c=======================================================================
c=======================================================================
*   -------------------------------------------------------------------
*
*                               sgp4unit.for
*
*    this file contains the sgp4 procedures for analytical propagation
*    of a satellite. the code was originally released in the 1980 and 1986
*    spacetrack papers. a detailed discussion of the theory and history
*    may be found in the 2006 aiaa paper by vallado, crawford, hujsak,
*    and kelso.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                    2007
*                              by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              26 Aug 08  david vallado
*                           fix atime for faster operation in dspace
*                           add operationmode for afspc (a) or improved (i)
*                           performance mode
*    changes :
*              16 jun 08  david vallado
*                           update small eccentricity check
*              16 nov 07  david vallado
*                           misc fixes for better compliance
*               2 apr 07  david vallado
*                           misc fixes for constants
*              14 aug 06  david vallado
*                           chg lyddane choice back to strn3, constants,
*                           separate debug and writes, misc doc
*              26 jul 05  david vallado
*                           fixes for paper
*                           note that each fix is preceded by a
*                           comment with "sgp4fix" and an explanation of
*                           what was changed
*              10 aug 04  david vallado
*                           2nd printing baseline working
*              14 may 01  david vallado
*                           2nd edition baseline
*                     80  norad
*                           original baseline
*
*     *****************************************************************
*  Files         :
*    Unit 14     - sgp4test.dbg    debug output file


* -----------------------------------------------------------------------------
*
*                           SUBROUTINE DPPER
*
*  This Subroutine provides deep space long period periodic contributions
*    to the mean elements.  by design, these periodics are zero at epoch.
*    this used to be dscom which included initialization, but it's really a
*    recurring function.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    e3          -
*    ee2         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    se2 , se3 , Sgh2, Sgh3, Sgh4, Sh2, Sh3, Si2, Si3, Sl2, Sl3, Sl4 -
*    t           -
*    xh2, xh3, xi2, xi3, xl2, xl3, xl4 -
*    zmol        -
*    zmos        -
*    ep          - eccentricity                           0.0 - 1.0
*    inclo       - inclination - needed for lyddane modification
*    nodep       - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  outputs       :
*    ep          - eccentricity                           0.0 - 1.0
*    inclp       - inclination
*    nodep       - right ascension of ascending node
*    argpp       - argument of perigee
*    mp          - mean anomaly
*
*  locals        :
*    alfdp       -
*    betdp       -
*    cosip  , sinip  , cosop  , sinop  ,
*    dalf        -
*    dbet        -
*    dls         -
*    f2, f3      -
*    pe          -
*    pgh         -
*    ph          -
*    pinc        -
*    pl          -
*    sel   , ses   , sghl  , sghs  , shl   , shs   , sil   , sinzf , sis   ,
*    sll   , sls
*    xls         -
*    xnoh        -
*    zf          -
*    zm          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
*------------------------------------------------------------------------------

      SUBROUTINE DPPER( e3    , ee2   , peo   , pgho  , pho   , pinco ,
     &                  plo   , se2   , se3   , sgh2  , sgh3  , sgh4  ,
     &                  sh2   , sh3   , si2   , si3   , sl2   , sl3   ,
     &                  sl4   , T     , xgh2  , xgh3  , xgh4  , xh2   ,
     &                  xh3   , xi2   , xi3   , xl2   , xl3   , xl4   ,
     &                  zmol  , zmos  , inclo , init  ,
     &                  Eccp  , Inclp , nodep, Argpp , Mp,
     &                  operationmode )
        IMPLICIT NONE
        CHARACTER Init, operationmode
        REAL*8  e3    , ee2   , peo   , pgho  , pho   , pinco , plo   ,
     &          se2   , se3   , sgh2  , sgh3  , sgh4  , sh2   , sh3   ,
     &          si2   , si3   , sl2   , sl3   , sl4   , T     , xgh2  ,
     &          xgh3  , xgh4  , xh2   , xh3   , xi2   , xi3   , xl2   ,
     &          xl3   , xl4   , zmol  , zmos  , inclo ,
     &          Eccp  , Inclp , nodep, Argpp , Mp

* -------------------------- Local Variables --------------------------
        REAL*8  alfdp , betdp , cosip , cosop , dalf  , dbet  , dls   ,
     &          f2    , f3    , pe    , pgh   , ph    , pinc  , pl    ,
     &          sel   , ses   , sghl  , sghs  , shl   , shs   , sil   ,
     &          sinip , sinop , sinzf , sis   , sll   , sls   , xls   ,
     &          xnoh  , zf    , zm
        REAL*8  Zel   , Zes   , Znl   , Zns
        COMMON /DebugHelp/ Help
        CHARACTER Help
        INCLUDE 'ASTMATH.CMN'

* ----------------------------- Constants -----------------------------
        ZES  = 0.01675D0
        ZEL  = 0.05490D0
        ZNS  = 1.19459D-5
        ZNL  = 1.5835218D-4

* ------------------- CALCULATE TIME VARYING PERIODICS ----------------
        ZM   = ZMOS + ZNS*T

        IF (Init.eq.'y') ZM = ZMOS
        ZF   = ZM + 2.0D0*ZES*DSIN(ZM)
        SINZF= DSIN(ZF)
        F2   =  0.5D0*SINZF*SINZF - 0.25D0
        F3   = -0.5D0*SINZF*DCOS(ZF)
        SES  = SE2*F2 + SE3*F3
        SIS  = SI2*F2 + SI3*F3
        SLS  = SL2*F2 + SL3*F3 + SL4*SINZF
        SGHS = SGH2*F2 + SGH3*F3 + SGH4*SINZF
        SHS  = SH2*F2 + SH3*F3
        ZM   = ZMOL + ZNL*T

        IF (Init.eq.'y') ZM = ZMOL
        ZF   = ZM + 2.0D0*ZEL*DSIN(ZM)
        SINZF= DSIN(ZF)
        F2   =  0.5D0*SINZF*SINZF - 0.25D0
        F3   = -0.5D0*SINZF*DCOS(ZF)
        SEL  = EE2*F2 + E3*F3
        SIL  = XI2*F2 + XI3*F3
        SLL  = XL2*F2 + XL3*F3 + XL4*SINZF
        SGHL = XGH2*F2 + XGH3*F3 + XGH4*SINZF
        SHL  = XH2*F2 + XH3*F3
        PE   = SES + SEL
        PINC = SIS + SIL
        PL   = SLS + SLL
        PGH  = SGHS + SGHL
        PH   = SHS + SHL

        IF (Init.eq.'n') THEN
            PE    = PE   - PEO
            PINC  = PINC - PINCO
            PL    = PL   - PLO
            PGH   = PGH  - PGHO
            PH    = PH   - PHO
            Inclp = Inclp  + PINC
            Eccp  = Eccp   + PE
            SINIP = DSIN(Inclp)
            COSIP = DCOS(Inclp)

* ------------------------- APPLY PERIODICS DIRECTLY ------------------
c    sgp4fix for lyddane choice
c    strn3 used original inclination - this is technically feasible
c    gsfc used perturbed inclination - also technically feasible
c    probably best to readjust the 0.2 limit value and limit discontinuity
c    0.2 rad = 11.45916 deg
c    use next line for original strn3 approach and original inclination
c            IF (inclo.ge.0.2D0) THEN
c    use next line for gsfc version and perturbed inclination
            IF (Inclp.ge.0.2D0) THEN

                PH     = PH/SINIP
                PGH    = PGH - COSIP*PH
                Argpp  = Argpp + PGH
                nodep  = nodep + PH
                Mp     = Mp + PL
              ELSE

* ----------------- APPLY PERIODICS WITH LYDDANE MODIFICATION ---------
                SINOP  = DSIN(nodep)
                COSOP  = DCOS(nodep)
                ALFDP  = SINIP*SINOP
                BETDP  = SINIP*COSOP
                DALF   =  PH*COSOP + PINC*COSIP*SINOP
                DBET   = -PH*SINOP + PINC*COSIP*COSOP
                ALFDP  = ALFDP + DALF
                BETDP  = BETDP + DBET
                nodep = DMOD(nodep,TwoPi)
                ! sgp4fix for afspc written intrinsic functions
                ! nodep used without a trigonometric function ahead
                IF ((nodep .LT. 0.0D0) .and. (operationmode .eq. 'a')) 
     &                THEN
                    nodep = nodep + twopi
                  ENDIF
                XLS    = Mp + Argpp + COSIP*nodep
                DLS    = PL + PGH - PINC*nodep*SINIP
                XLS    = XLS + DLS
                XNOH   = nodep
                nodep  = DATAN2(ALFDP,BETDP)
                ! sgp4fix for afspc written intrinsic functions
                ! nodep used without a trigonometric function ahead
                IF ((nodep .LT. 0.0D0) .and. (operationmode .eq. 'a')) 
     &                THEN
                    nodep = nodep + twopi
                  ENDIF
                IF (DABS(XNOH-nodep) .GT. PI) THEN
                    IF(nodep .lt. XNOH) THEN
                        nodep = nodep+TWOPI
                      ELSE
                        nodep = nodep-TWOPI
                      ENDIF
                  ENDIF
                Mp   = Mp + PL
                Argpp=  XLS - Mp - COSIP*nodep
              ENDIF
          ENDIF

c        INCLUDE 'debug1.for'

      RETURN
      END  !  end dpper


* -----------------------------------------------------------------------------
*
*                           SUBROUTINE DSCOM
*
*  This Subroutine provides deep space common items used by both the secular
*    and periodics subroutines.  input is provided as shown. this routine
*    used to be called dpper, but the functions inside weren't well organized.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    epoch       -
*    ep          - eccentricity
*    argpp       - argument of perigee
*    tc          -
*    inclp       - inclination
*    nodep      - right ascension of ascending node
*    np          - mean motion
*
*  outputs       :
*    sinim  , cosim  , sinomm , cosomm , snodm  , cnodm
*    day         -
*    e3          -
*    ee2         -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    gam         -
*    peo         -
*    pgho        -
*    pho         -
*    pinco       -
*    plo         -
*    rtemsq      -
*    se2, se3         -
*    sgh2, sgh3, sgh4        -
*    sh2, sh3, si2, si3, sl2, sl3, sl4         -
*    s1, s2, s3, s4, s5, s6, s7          -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3         -
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4         -
*    nm          - mean motion
*    z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*    zmol        -
*    zmos        -
*
*  locals        :
*    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10         -
*    betasq      -
*    cc          -
*    ctem, stem        -
*    x1, x2, x3, x4, x5, x6, x7, x8          -
*    xnodce      -
*    xnoi        -
*    zcosg  , zsing  , zcosgl , zsingl , zcosh  , zsinh  , zcoshl , zsinhl ,
*    zcosi  , zsini  , zcosil , zsinil ,
*    zx          -
*    zy          -
*
*  coupling      :
*    none.
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
*------------------------------------------------------------------------------

      SUBROUTINE DSCOM( EPOCH , Eccp  , Argpp , Tc    , Inclp , nodep,
     &                  Np    ,
     &                  SNODM , CNODM , SINIM , COSIM , SINOMM, COSOMM,
     &                  DAY   , E3    , Ee2   , Eccm  , EMSQ  , GAM   ,
     &                  Peo   , Pgho  , Pho   , PInco , Plo   ,
     &                  RTemSq, Se2   , Se3   , Sgh2  , Sgh3  , Sgh4  ,
     &                  Sh2   , Sh3   , Si2   , Si3   , Sl2   , Sl3   ,
     &                  Sl4   , S1    , S2    , S3    , S4    , S5    ,
     &                  S6    , S7    , SS1   , SS2   , SS3   , SS4   ,
     &                  SS5   , SS6   , SS7   , SZ1   , SZ2   , SZ3   ,
     &                  SZ11  , SZ12  , SZ13  , SZ21  , SZ22  , SZ23  ,
     &                  SZ31  , SZ32  , SZ33  , Xgh2  , Xgh3  , Xgh4  ,
     &                  Xh2   , Xh3   , Xi2   , Xi3   , Xl2   , Xl3   ,
     &                  Xl4   , Xn    , Z1    , Z2    , Z3    , Z11   ,
     &                  Z12   , Z13   , Z21   , Z22   , Z23   , Z31   ,
     &                  Z32   , Z33   , Zmol  , Zmos )
        IMPLICIT NONE
        REAL*8  EPOCH , Eccp  , Argpp , Tc    , Inclp , nodep, Np    ,
     &          SNODM , CNODM , SINIM , COSIM , SINOMM, COSOMM, DAY   ,
     &          E3    , Ee2   , Eccm  , EMSQ  , GAM   , RTemSq, Se2   ,
     &          Peo   , Pgho  , Pho   , PInco , Plo   ,
     &          Se3   , Sgh2  , Sgh3  , Sgh4  , Sh2   , Sh3   , Si2   ,
     &          Si3   , Sl2   , Sl3   , Sl4   , S1    , S2    , S3    ,
     &          S4    , S5    , S6    , S7    , SS1   , SS2   , SS3   ,
     &          SS4   , SS5   , SS6   , SS7   , SZ1   , SZ2   , SZ3   ,
     &          SZ11  , SZ12  , SZ13  , SZ21  , SZ22  , SZ23  , SZ31  ,
     &          SZ32  , SZ33  , Xgh2  , Xgh3  , Xgh4  , Xh2   , Xh3   ,
     &          Xi2   , Xi3   , Xl2   , Xl3   , Xl4   , Xn    , Z1    ,
     &          Z2    , Z3    , Z11   , Z12   , Z13   , Z21   , Z22   ,
     &          Z23   , Z31   , Z32   , Z33   , Zmol  , Zmos

* -------------------------- Local Variables --------------------------
        REAL*8  c1ss  , c1L   , zcosis, zsinis, zsings, zcosgs,
     &          Zes   , zel
        INTEGER LsFlg
        REAL*8  a1    , a2    , a3    , a4    , a5    , a6    , a7    ,
     &          a8    , a9    , a10   , betasq, cc    , ctem  , stem  ,
     &          x1    , x2    , x3    , x4    , x5    , x6    , x7    ,
     &          x8    , xnodce, xnoi  , zcosg , zcosgl, zcosh , zcoshl,
     &          zcosi , zcosil, zsing , zsingl, zsinh , zsinhl, zsini ,
     &          zsinil, zx    , zy

        COMMON /DebugHelp/ Help
        CHARACTER Help
        INCLUDE 'ASTMATH.CMN'

* ------------------------------ Constants ----------------------------
        ZES    =  0.01675D0
        ZEL    =  0.05490D0
        C1SS   =  2.9864797D-6
        C1L    =  4.7968065D-7
        ZSINIS =  0.39785416D0
        ZCOSIS =  0.91744867D0
        ZCOSGS =  0.1945905D0
        ZSINGS = -0.98088458D0

* ----------------- DEEP SPACE PERIODICS INITIALIZATION ---------------
        XN     = Np
        Eccm   = Eccp
        SNODM  = DSIN(nodep)
        CNODM  = DCOS(nodep)
        SINOMM = DSIN(Argpp)
        COSOMM = DCOS(Argpp)
        SINIM  = DSIN(Inclp)
        COSIM  = DCOS(Inclp)
        EMSQ   = Eccm*Eccm
        BETASQ = 1.0D0-EMSQ
        RTEMSQ = DSQRT(BETASQ)

* --------------------- INITIALIZE LUNAR SOLAR TERMS ------------------
        PEO    = 0.0D0
        PINCO  = 0.0D0
        PLO    = 0.0D0
        PGHO   = 0.0D0
        PHO    = 0.0D0
        DAY    = EPOCH + 18261.5D0 + TC/1440.0D0
        XNODCE = DMOD(4.5236020D0 - 9.2422029D-4*DAY,TwoPi)
        STEM   = DSIN(XNODCE)
        CTEM   = DCOS(XNODCE)
        ZCOSIL = 0.91375164D0 - 0.03568096D0*CTEM
        ZSINIL = DSQRT(1.0D0 - ZCOSIL*ZCOSIL)
        ZSINHL = 0.089683511D0*STEM / ZSINIL
        ZCOSHL = DSQRT(1.0D0 - ZSINHL*ZSINHL)
        GAM    = 5.8351514D0 + 0.0019443680D0*DAY
        ZX     = 0.39785416D0*STEM/ZSINIL
        ZY     = ZCOSHL*CTEM + 0.91744867D0*ZSINHL*STEM
        ZX     = DATAN2(ZX,ZY)
        ZX     = GAM + ZX - XNODCE
        ZCOSGL = DCOS(ZX)
        ZSINGL = DSIN(ZX)

* ---------------------------- DO SOLAR TERMS -------------------------
        ZCOSG = ZCOSGS
        ZSING = ZSINGS
        ZCOSI = ZCOSIS
        ZSINI = ZSINIS
        ZCOSH = CNODM
        ZSINH = SNODM
        CC    = C1SS
        XNOI  = 1.0D0 / XN

        DO LSFlg = 1,2
            A1 =   ZCOSG*ZCOSH + ZSING*ZCOSI*ZSINH
            A3 =  -ZSING*ZCOSH + ZCOSG*ZCOSI*ZSINH
            A7 =  -ZCOSG*ZSINH + ZSING*ZCOSI*ZCOSH
            A8 =   ZSING*ZSINI
            A9 =   ZSING*ZSINH + ZCOSG*ZCOSI*ZCOSH
            A10=   ZCOSG*ZSINI
            A2 =   COSIM*A7 + SINIM*A8
            A4 =   COSIM*A9 + SINIM*A10
            A5 =  -SINIM*A7 + COSIM*A8
            A6 =  -SINIM*A9 + COSIM*A10

            X1 =  A1*COSOMM + A2*SINOMM
            X2 =  A3*COSOMM + A4*SINOMM
            X3 = -A1*SINOMM + A2*COSOMM
            X4 = -A3*SINOMM + A4*COSOMM
            X5 =  A5*SINOMM
            X6 =  A6*SINOMM
            X7 =  A5*COSOMM
            X8 =  A6*COSOMM

            Z31= 12.0D0*X1*X1 - 3.0D0*X3*X3
            Z32= 24.0D0*X1*X2 - 6.0D0*X3*X4
            Z33= 12.0D0*X2*X2 - 3.0D0*X4*X4
            Z1 =  3.0D0* (A1*A1 + A2*A2) + Z31*EMSQ
            Z2 =  6.0D0* (A1*A3 + A2*A4) + Z32*EMSQ
            Z3 =  3.0D0* (A3*A3 + A4*A4) + Z33*EMSQ
            Z11= -6.0D0*A1*A5 + EMSQ* (-24.0D0*X1*X7-6.0D0*X3*X5)
            Z12= -6.0D0* (A1*A6 + A3*A5) + EMSQ*
     &           ( -24.0D0*(X2*X7+X1*X8) - 6.0D0*(X3*X6+X4*X5) )
            Z13= -6.0D0*A3*A6 + EMSQ*(-24.0D0*X2*X8 - 6.0D0*X4*X6)
            Z21=  6.0D0*A2*A5 + EMSQ*(24.0D0*X1*X5-6.0D0*X3*X7)
            Z22=  6.0D0* (A4*A5 + A2*A6) + EMSQ*
     &           (  24.0D0*(X2*X5+X1*X6) - 6.0D0*(X4*X7+X3*X8) )
            Z23=  6.0D0*A4*A6 + EMSQ*(24.0D0*X2*X6 - 6.0D0*X4*X8)
            Z1 = Z1 + Z1 + BETASQ*Z31
            Z2 = Z2 + Z2 + BETASQ*Z32
            Z3 = Z3 + Z3 + BETASQ*Z33
            S3 = CC*XNOI
            S2 = -0.5D0*S3 / RTEMSQ
            S4 = S3*RTEMSQ
            S1 = -15.0D0*Eccm*S4
            S5 = X1*X3 + X2*X4
            S6 = X2*X3 + X1*X4
            S7 = X2*X4 - X1*X3

* ------------------------------ DO LUNAR TERMS -----------------------
            IF (LSFLG.eq.1) THEN
                SS1   = S1
                SS2   = S2
                SS3   = S3
                SS4   = S4
                SS5   = S5
                SS6   = S6
                SS7   = S7
                SZ1   = Z1
                SZ2   = Z2
                SZ3   = Z3
                SZ11  = Z11
                SZ12  = Z12
                SZ13  = Z13
                SZ21  = Z21
                SZ22  = Z22
                SZ23  = Z23
                SZ31  = Z31
                SZ32  = Z32
                SZ33  = Z33
                ZCOSG = ZCOSGL
                ZSING = ZSINGL
                ZCOSI = ZCOSIL
                ZSINI = ZSINIL
                ZCOSH = ZCOSHL*CNODM+ZSINHL*SNODM
                ZSINH = SNODM*ZCOSHL-CNODM*ZSINHL
                CC    = C1L
              ENDIF
          ENDDO

        ZMOL  = DMOD( 4.7199672D0 + 0.22997150D0*DAY-GAM,TwoPi )
        ZMOS  = DMOD( 6.2565837D0 + 0.017201977D0*DAY,TwoPi )

* ---------------------------- DO SOLAR TERMS -------------------------
        SE2 =   2.0D0*SS1*SS6
        SE3 =   2.0D0*SS1*SS7
        SI2 =   2.0D0*SS2*SZ12
        SI3 =   2.0D0*SS2*(SZ13-SZ11)
        SL2 =  -2.0D0*SS3*SZ2
        SL3 =  -2.0D0*SS3*(SZ3-SZ1)
        SL4 =  -2.0D0*SS3*(-21.0D0-9.0D0*EMSQ)*ZES
        SGH2=   2.0D0*SS4*SZ32
        SGH3=   2.0D0*SS4*(SZ33-SZ31)
        SGH4= -18.0D0*SS4*ZES
        SH2 =  -2.0D0*SS2*SZ22
        SH3 =  -2.0D0*SS2*(SZ23-SZ21)

* ---------------------------- DO LUNAR TERMS -------------------------
        EE2 =   2.0D0*S1*S6
        E3  =   2.0D0*S1*S7
        XI2 =   2.0D0*S2*Z12
        XI3 =   2.0D0*S2*(Z13-Z11)
        XL2 =  -2.0D0*S3*Z2
        XL3 =  -2.0D0*S3*(Z3-Z1)
        XL4 =  -2.0D0*S3*(-21.0D0-9.0D0*EMSQ)*ZEL
        XGH2=   2.0D0*S4*Z32
        XGH3=   2.0D0*S4*(Z33-Z31)
        XGH4= -18.0D0*S4*ZEL
        XH2 =  -2.0D0*S2*Z22
        XH3 =  -2.0D0*S2*(Z23-Z21)

c        INCLUDE 'debug2.for'

      RETURN
      END  !  dscom


* -----------------------------------------------------------------------------
*
*                           SUBROUTINE DSINIT
*
*  This Subroutine provides Deep Space contributions to Mean Motion Dot due
*    to geopotential resonance with half day and one day orbits.
*
*  Inputs        :
*    Cosim, Sinim-
*    Emsq        - Eccentricity squared
*    Argpo       - Argument of Perigee
*    S1, S2, S3, S4, S5      -
*    Ss1, Ss2, Ss3, Ss4, Ss5 -
*    Sz1, Sz3, Sz11, Sz13, Sz21, Sz23, Sz31, Sz33 -
*    T           - Time
*    Tc          -
*    GSTo        - Greenwich sidereal time                   rad
*    Mo          - Mean Anomaly
*    MDot        - Mean Anomaly dot (rate)
*    No          - Mean Motion
*    nodeo       - right ascension of ascending node
*    nodeDot     - right ascension of ascending node dot (rate)
*    XPIDOT      -
*    Z1, Z3, Z11, Z13, Z21, Z23, Z31, Z33 -
*    Eccm        - Eccentricity
*    Argpm       - Argument of perigee
*    Inclm       - Inclination
*    Mm          - Mean Anomaly
*    Xn          - Mean Motion
*    nodem       - right ascension of ascending node
*
*  Outputs       :
*    Eccm        - Eccentricity
*    Argpm       - Argument of perigee
*    Inclm       - Inclination
*    Mm          - Mean Anomaly
*    Xn          - Mean motion
*    nodem       - right ascension of ascending node
*    IRez        - Resonance flags              0-none, 1-One day,  2-Half day
*    Atime       -
*    D2201, D2211, D3210, D3222, D4410, D4422, D5220, D5232, D5421, D5433       -
*    Dedt        -
*    Didt        -
*    DMDT        -
*    DNDT        -
*    DNODT       -
*    DOMDT       -
*    Del1, Del2, Del3 -
*    Ses  , Sghl , Sghs , Sgs  , Shl  , Shs  , Sis  , Sls
*    THETA       -
*    Xfact       -
*    Xlamo       -
*    Xli         -
*    Xni
*
*  Locals        :
*    ainv2       -
*    aonv        -
*    cosisq      -
*    eoc         -
*    f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543        -
*    g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533        -
*    sini2       -
*    temp, temp1 -
*    Theta       -
*    xno2        -
*
*  Coupling      :
*    getgravconst-
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
*------------------------------------------------------------------------------

      SUBROUTINE DSINIT( whichconst,
     &                   Cosim , Emsq  , Argpo , S1    , S2    , S3    ,
     &                   S4    , S5    , Sinim , Ss1   , Ss2   , Ss3   ,
     &                   Ss4   , Ss5   , Sz1   , Sz3   , Sz11  , Sz13  ,
     &                   Sz21  , Sz23  , Sz31  , Sz33  , T     , Tc    ,
     &                   GSTo  , Mo    , MDot  , No    , nodeo ,nodeDot,
     &                   XPIDOT, Z1    , Z3    , Z11   , Z13   , Z21   ,
     &                   Z23   , Z31   , Z33   , Ecco  , EccSq ,
     &                   Eccm  , Argpm , Inclm , Mm    , Xn    , nodem,
     &                   IREZ  , Atime , D2201 , D2211 , D3210 , D3222 ,
     &                   D4410 , D4422 , D5220 , D5232 , D5421 , D5433 ,
     &                   Dedt  , Didt  , DMDT  , DNDT  , DNODT , DOMDT ,
     &                   Del1  , Del2  , Del3  , Xfact , Xlamo , Xli   ,
     &                   Xni )
        IMPLICIT NONE
        INTEGER  IRez, whichconst
        REAL*8   Cosim , Emsq  , Argpo , S1    , S2    , S3    , S4    ,
     &           S5    , Sinim , Ss1   , Ss2   , Ss3   , Ss4   , Ss5   ,
     &           Sz1   , Sz3   , Sz11  , Sz13  , Sz21  , Sz23  , Sz31  ,
     &           Sz33  , T     , Tc    , GSTo  , Mo    , MDot  , No    ,
     &           nodeo ,nodeDot,XPIDOT , Z1    , Z3    , Z11   , Z13   ,
     &           Z21   , Z23   , Z31   , Z33   , Eccm  , Argpm , Inclm ,
     &           Mm    , Xn    , nodem , Atime , D2201 , D2211 , D3210 ,
     &           D3222 , D4410 , D4422 , D5220 , D5232 , D5421 , D5433 ,
     &           Dedt  , Didt  , DMDT  , DNDT  , DNODT , DOMDT , Del1  ,
     &           Del2  , Del3  , Xfact , Xlamo , Xli   , Xni   , Ecco  ,
     &           Eccsq

* -------------------------- Local Variables --------------------------
        REAL*8  ainv2 , aonv  , cosisq, eoc   , f220  , f221  , f311  ,
     &          f321  , f322  , f330  , f441  , f442  , f522  , f523  ,
     &          f542  , f543  , g200  , g201  , g211  , g300  , g310  ,
     &          g322  , g410  , g422  , g520  , g521  , g532  , g533  ,
     &          ses   , sgs   , sghl  , sghs  , shs   , shl   , sis   ,
     &          sini2 , sls   , temp  , temp1 , Theta , xno2
        REAL*8  Q22   , Q31   , Q33   , ROOT22, ROOT44, ROOT54,
     &          RPTim , Root32, Root52, X2o3  , XKe   , Znl   ,
     &          Zns,  Emo, emsqo , tumin, mu, radiusearthkm, j2, j3, j4,
     &          j3oj2

        COMMON /DebugHelp/ Help
        CHARACTER Help
        INCLUDE 'ASTMATH.CMN'

        Q22    = 1.7891679D-6
        Q31    = 2.1460748D-6
        Q33    = 2.2123015D-7
        ROOT22 = 1.7891679D-6
        ROOT44 = 7.3636953D-9
        ROOT54 = 2.1765803D-9
        RPTim  = 4.37526908801129966D-3 ! this equates to 7.29211514668855e-5 rad/sec
        Root32 = 3.7393792D-7
        Root52 = 1.1428639D-7
        X2o3   = 2.0D0 / 3.0D0
        ZNL    = 1.5835218D-4
        ZNS    = 1.19459D-5

        ! sgp4fix identify constants and allow alternate values
        CALL getgravconst( whichconst, tumin, mu, radiusearthkm, xke,
     &       j2, j3, j4, j3oj2 )

* ------------------------ DEEP SPACE INITIALIZATION ------------------
        IREZ = 0
        IF ((XN.lt.0.0052359877D0).AND.(XN.GT.0.0034906585D0)) THEN
            IREZ = 1
          ENDIF
        IF ((XN.ge.8.26D-3).AND.(XN.LE.9.24D-3).AND.(Eccm.GE.0.5D0))THEN
            IREZ = 2
          ENDIF

* ---------------------------- DO SOLAR TERMS -------------------------
        SES  =  SS1*ZNS*SS5
        SIS  =  SS2*ZNS*(SZ11 + SZ13)
        SLS  = -ZNS*SS3*(SZ1 + SZ3 - 14.0D0 - 6.0D0*EMSQ)
        SGHS =  SS4*ZNS*(SZ31 + SZ33 - 6.0D0)
        SHS  = -ZNS*SS2*(SZ21 + SZ23)
c       sgp4fix for 180 deg incl
        IF ((Inclm.lt.5.2359877D-2).or.(Inclm.gt.pi-5.2359877D-2)) THEN
            SHS = 0.0D0
          ENDIF
        IF (SINIM.ne.0.0D0) THEN
            SHS = SHS/SINIM
          ENDIF
        SGS  = SGHS - COSIM*SHS

* ----------------------------- DO LUNAR TERMS ------------------------
        DEDT = SES + S1*ZNL*S5
        DIDT = SIS + S2*ZNL*(Z11 + Z13)
        DMDT = SLS - ZNL*S3*(Z1 + Z3 - 14.0D0 - 6.0D0*EMSQ)
        SGHL = S4*ZNL*(Z31 + Z33 - 6.0D0)
        SHL  = -ZNL*S2*(Z21 + Z23)
c       sgp4fix for 180 deg incl
        IF ((Inclm.lt.5.2359877D-2).or.(Inclm.gt.pi-5.2359877D-2)) THEN
            SHL = 0.0D0
          ENDIF
        DOMDT= SGS+SGHL
        DNODT= SHS
        IF (SINIM .ne. 0.0D0) THEN
            DOMDT = DOMDT-COSIM/SINIM*SHL
            DNODT = DNODT+SHL/SINIM
        ENDIF

* --------------- CALCULATE DEEP SPACE RESONANCE EFFECTS --------------
        DNDT  = 0.0D0
        THETA = DMOD(GSTo + TC*RPTIM,TwoPi)
        Eccm  = Eccm + DEDT*T
        emsq  = eccm**2
        Inclm = Inclm + DIDT*T
        Argpm = Argpm + DOMDT*T
        nodem = nodem + DNODT*T
        Mm    = Mm + DMDT*T
c   sgp4fix for negative inclinations
c   the following if statement should be commented out
c           IF(Inclm .lt. 0.0D0) THEN
c             Inclm  = -Inclm
c             Argpm  = Argpm-PI
c             nodem = nodem+PI
c           ENDIF

* ------------------ Initialize the resonance terms -------------------
        IF (IREZ .ne. 0) THEN
            AONV = (XN/XKE)**X2O3

* -------------- GEOPOTENTIAL RESONANCE FOR 12 HOUR ORBITS ------------
        IF (IREZ .eq. 2) THEN
            COSISQ = COSIM*COSIM
            emo    = Eccm
            emsqo  = emsq
            Eccm   = ecco
            emsq   = eccsq
            EOC    = Eccm*EMSQ
            G201   = -0.306D0-(Eccm-0.64D0)*0.440D0
            IF (Eccm.le.0.65D0) THEN
                G211 =   3.616D0 -  13.2470D0*Eccm +  16.2900D0*EMSQ
                G310 = -19.302D0 + 117.3900D0*Eccm - 228.4190D0*EMSQ +
     &                 156.591D0*EOC
                G322 = -18.9068D0+ 109.7927D0*Eccm - 214.6334D0*EMSQ +
     &                 146.5816D0*EOC
                G410 = -41.122D0 + 242.6940D0*Eccm - 471.0940D0*EMSQ +
     &                 313.953D0*EOC
                G422 =-146.407D0 + 841.8800D0*Eccm - 1629.014D0*EMSQ +
     &                1083.435D0*EOC
                G520 =-532.114D0 + 3017.977D0*Eccm - 5740.032D0*EMSQ +
     &                3708.276D0*EOC
              ELSE
                G211 =  -72.099D0 +  331.819D0*Eccm -  508.738D0*EMSQ +
     &                  266.724D0*EOC
                G310 = -346.844D0 + 1582.851D0*Eccm - 2415.925D0*EMSQ +
     &                 1246.113D0*EOC
                G322 = -342.585D0 + 1554.908D0*Eccm - 2366.899D0*EMSQ +
     &                 1215.972D0*EOC
                G410 =-1052.797D0 + 4758.686D0*Eccm - 7193.992D0*EMSQ +
     &                 3651.957D0*EOC
                G422 =-3581.690D0 + 16178.11D0*Eccm - 24462.77D0*EMSQ +
     &                12422.52D0*EOC
                IF (Eccm.gt.0.715D0) THEN
                    G520 =-5149.66D0 + 29936.92D0*Eccm -54087.36D0*EMSQ
     &                    + 31324.56D0*EOC
                  ELSE
                    G520 = 1464.74D0 -  4664.75D0*Eccm + 3763.64D0*EMSQ
                  ENDIF
              ENDIF
            IF (Eccm.lt.0.7D0) THEN
                G533 = -919.22770D0 + 4988.6100D0*Eccm-9064.7700D0*EMSQ
     &               + 5542.21D0*EOC
                G521 = -822.71072D0 + 4568.6173D0*Eccm-8491.4146D0*EMSQ
     &               + 5337.524D0*EOC
                G532 = -853.66600D0 + 4690.2500D0*Eccm-8624.7700D0*EMSQ
     &               + 5341.4D0*EOC
              ELSE
                G533 =-37995.780D0 + 161616.52D0*Eccm-229838.20D0*EMSQ+
     &              109377.94D0*EOC
                G521 =-51752.104D0 + 218913.95D0*Eccm-309468.16D0*EMSQ+
     &              146349.42D0*EOC
                G532 =-40023.880D0 + 170470.89D0*Eccm-242699.48D0*EMSQ+
     &              115605.82D0*EOC
              ENDIF
            SINI2 =  SINIM*SINIM
            F220  =  0.75D0* (1.0D0+2.0D0*COSIM+COSISQ)
            F221  =  1.5D0*SINI2
            F321  =  1.875D0*SINIM * (1.0D0-2.0D0*COSIM-3.0D0*COSISQ)
            F322  = -1.875D0*SINIM * (1.0D0+2.0D0*COSIM-3.0D0*COSISQ)
            F441  = 35.0D0*SINI2*F220
            F442  = 39.3750D0*SINI2*SINI2
            F522  =  9.84375D0*SINIM * (SINI2* (1.0D0-2.0D0*COSIM-
     &               5.0D0*COSISQ)+0.33333333D0 * (-2.0D0+4.0D0*COSIM+
     &               6.0D0*COSISQ) )
            F523  =  SINIM * (4.92187512D0*SINI2 * (-2.0D0-4.0D0*COSIM+
     &               10.0D0*COSISQ) + 6.56250012D0*
     &               (1.0D0+2.0D0*COSIM-3.0D0*COSISQ))
            F542  =  29.53125D0*SINIM * (2.0D0-8.0D0*COSIM+COSISQ*
     &               (-12.0D0+8.0D0*COSIM+10.0D0*COSISQ) )
            F543  = 29.53125D0*SINIM * (-2.0D0-8.0D0*COSIM+COSISQ*
     &               (12.0D0+8.0D0*COSIM-10.0D0*COSISQ) )

            XNO2   =  XN * XN
            AINV2  =  AONV * AONV
            TEMP1  =  3.0D0*XNO2*AINV2
            TEMP   =  TEMP1*ROOT22
            D2201  =  TEMP*F220*G201
            D2211  =  TEMP*F221*G211
            TEMP1  =  TEMP1*AONV
            TEMP   =  TEMP1*ROOT32
            D3210  =  TEMP*F321*G310
            D3222  =  TEMP*F322*G322
            TEMP1  =  TEMP1*AONV
            TEMP   =  2.0D0*TEMP1*ROOT44
            D4410  =  TEMP*F441*G410
            D4422  =  TEMP*F442*G422
            TEMP1  =  TEMP1*AONV
            TEMP   =  TEMP1*ROOT52
            D5220  =  TEMP*F522*G520
            D5232  =  TEMP*F523*G532
            TEMP   =  2.0D0*TEMP1*ROOT54
            D5421  =  TEMP*F542*G521
            D5433  =  TEMP*F543*G533
            XLAMO  =  DMOD(Mo+nodeo+nodeo-THETA-THETA,TwoPi)
            XFACT  = MDot + DMDT + 2.0D0 * (nodeDot+DNODT-RPTIM) - No

            Eccm = emo
            emsq = emsqo
          ENDIF

        IF (Irez .eq. 1) THEN
* -------------------- SYNCHRONOUS RESONANCE TERMS --------------------
            G200  = 1.0D0 + EMSQ * (-2.5D0+0.8125D0*EMSQ)
            G310  = 1.0D0 + 2.0D0*EMSQ
            G300  = 1.0D0 + EMSQ * (-6.0D0+6.60937D0*EMSQ)
            F220  = 0.75D0 * (1.0D0+COSIM) * (1.0D0+COSIM)
            F311  = 0.9375D0*SINIM*SINIM*
     &               (1.0D0+3.0D0*COSIM) - 0.75D0*(1.0D0+COSIM)
            F330  = 1.0D0+COSIM
            F330  = 1.875D0*F330*F330*F330
            DEL1  = 3.0D0*XN*XN*AONV*AONV
            DEL2  = 2.0D0*DEL1*F220*G200*Q22
            DEL3  = 3.0D0*DEL1*F330*G300*Q33*AONV
            DEL1  = DEL1*F311*G310*Q31*AONV
            XLAMO = DMOD(Mo+nodeo+Argpo-THETA,TwoPi)
            XFACT = MDot + XPIDOT - RPTIM + DMDT + DOMDT + DNODT - No
          ENDIF

* ---------------- FOR SGP4, INITIALIZE THE INTEGRATOR ----------------
         XLI   = XLAMO
         XNI   = No
         ATIME = 0.0D0
         XN    = No + DNDT
      ENDIF ! Ires non-zero

c        INCLUDE 'debug3.for'

      RETURN
      END  ! end dsinit


* -----------------------------------------------------------------------------
*
*                           SUBROUTINE DSPACE
*
*  This Subroutine provides deep space contributions to mean elements for
*    perturbing third body.  these effects have been averaged over one
*    revolution of the sun and moon.  for earth resonance effects, the
*    effects have been averaged over no revolutions of the satellite.
*    (mean motion)
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433       -
*    dedt        -
*    del1, del2, del3  -
*    didt        -
*    dmdt        -
*    dnodt       -
*    domdt       -
*    irez        - flag for resonance           0-none, 1-one day, 2-half day
*    argpo       - argument of perigee
*    argpdot     - argument of perigee dot (rate)
*    t           - time
*    tc          -
*    gsto        - gst
*    xfact       -
*    xlamo       -
*    no          - mean motion
*    atime       -
*    em          - eccentricity
*    ft          -
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         - mean motion
*    nodem       - right ascension of ascending node
*
*  outputs       :
*    atime       -
*    em          - eccentricity
*    argpm       - argument of perigee
*    inclm       - inclination
*    xli         -
*    mm          - mean anomaly
*    xni         -
*    nodem       - right ascension of ascending node
*    dndt        -
*    nm          - mean motion
*
*  locals        :
*    delt        -
*    ft          -
*    theta       -
*    x2li        -
*    x2omi       -
*    xl          -
*    xldot       -
*    xnddt       -
*    xndt        -
*    xomi        -
*
*  coupling      :
*    none        -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
*------------------------------------------------------------------------------

      SUBROUTINE DSPACE( IRez  , D2201 , D2211 , D3210 , D3222 , D4410 ,
     &                   D4422 , D5220 , D5232 , D5421 , D5433 , Dedt  ,
     &                   Del1  , Del2  , Del3  , Didt  , Dmdt  , Dnodt ,
     &                   Domdt , Argpo , ArgpDot, T    , TC    , GSTo  ,
     &                   Xfact , Xlamo , No    ,
     &                   Atime , Eccm  , Argpm , Inclm , Xli   , Mm  ,
     &                   XNi   , nodem, Dndt  , XN  )
        IMPLICIT NONE
        INTEGER  IRez
        Real*8   D2201 , D2211 , D3210 , D3222 , D4410 , D4422 , D5220 ,
     &           D5232 , D5421 , D5433 , Dedt  , Del1  , Del2  , Del3  ,
     &           Didt  , Dmdt  , Dnodt , Domdt , Argpo , ArgpDot,T     ,
     &           TC    , GSTo  , Xfact , Xlamo , No    , Atime , Eccm  ,
     &           Argpm , Inclm , Xli   , Mm    , Xni   , nodem, Dndt  ,
     &           XN

* -------------------------- Local Variables --------------------------
        INTEGER  iretn , iret
        REAL*8   Delt  , Ft    , theta , x2li  , x2omi , xl    , xldot ,
     &           xnddt , xndt  , xomi
        REAL*8   G22   , G32   , G44   , G52   , G54   , Fasx2 ,
     &           Fasx4 , Fasx6 , RPtim , Step2 , Stepn , Stepp

        COMMON /DebugHelp/ Help
        CHARACTER Help
        INCLUDE 'ASTMATH.CMN'

* ----------------------------- Constants -----------------------------
        FASX2 = 0.13130908D0
        FASX4 = 2.8843198D0
        FASX6 = 0.37448087D0
        G22   = 5.7686396D0
        G32   = 0.95240898D0
        G44   = 1.8014998D0
        G52   = 1.0508330D0
        G54   = 4.4108898D0
        RPTIM = 4.37526908801129966D-3
        STEPP =    720.0D0
        STEPN =   -720.0D0
        STEP2 = 259200.0D0

* --------------- CALCULATE DEEP SPACE RESONANCE EFFECTS --------------
        DNDT  = 0.0D0
        THETA = DMOD(GSTo + TC*RPTIM,TwoPi)
        Eccm  = Eccm + DEDT*T

        Inclm = Inclm + DIDT*T
        Argpm = Argpm + DOMDT*T
        nodem = nodem + DNODT*T
        Mm    = Mm + DMDT*T

c   sgp4fix for negative inclinations
c   the following if statement should be commented out
c        IF(Inclm .lt. 0.0D0) THEN
c            Inclm  = -Inclm
c            Argpm  = Argpm-PI
c            nodem = nodem+PI
c          ENDIF

c   sgp4fix for propagator problems
c   the following integration works for negative time steps and periods
c   the specific changes are unknown because the original code was so convoluted
c      sgp4fix take out atime = 0.0 and fix for faster operation
        Ft    = 0.0D0      ! Just in case - should be set in loops if used.

        IF (IREZ .ne. 0) THEN
* ----- UPDATE RESONANCES : NUMERICAL (EULER-MACLAURIN) INTEGRATION ---
* ---------------------------- EPOCH RESTART --------------------------
         ! sgp4fix streamline check
         IF ((atime .eq. 0.0D0) .or. (t * atime .le. 0.0D0) .or. 
     &       (dabs(t) .lt. dabs(atime)) ) THEN
               atime  = 0.0D0
               xni    = no
               xli    = xlamo
            ENDIF
           ! sgp4fix move check outside loop
           IF (t .gt. 0.0D0) THEN
               delt = stepp
             else
               delt = stepn
             ENDIF

            iretn = 381 ! added for do loop
            iret  =   0 ! added for loop
            DO WHILE (IRetn.eq.381)

* --------------------------- DOT TERMS CALCULATED --------------------
* ------------------- NEAR - SYNCHRONOUS RESONANCE TERMS --------------
            IF (IREZ .ne. 2) THEN
                XNDT  = DEL1*DSIN(XLI-FASX2) +
     &                  DEL2*DSIN(2.0D0*(XLI-FASX4)) +
     &                  DEL3*DSIN(3.0D0*(XLI-FASX6))
                XLDOT = XNI + XFACT
                XNDDT = DEL1*DCOS(XLI-FASX2) +
     &            2.0D0*DEL2*DCOS(2.0D0*(XLI-FASX4)) +
     &            3.0D0*DEL3*DCOS(3.0D0*(XLI-FASX6))
                XNDDT = XNDDT*XLDOT
              ELSE

* --------------------- NEAR - HALF-DAY RESONANCE TERMS ---------------
                XOMI = Argpo + ArgpDot*ATIME
                X2OMI= XOMI + XOMI
                X2LI = XLI + XLI
                XNDT = D2201*DSIN(X2OMI+XLI-G22) + D2211*DSIN(XLI-G22) +
     &                 D3210*DSIN( XOMI+XLI-G32) +
     &                 D3222*DSIN(-XOMI+XLI-G32) +
     &                 D4410*DSIN(X2OMI+X2LI-G44)+ D4422*DSIN(X2LI-G44)+
     &                 D5220*DSIN( XOMI+XLI-G52) +
     &                 D5232*DSIN(-XOMI+XLI-G52) +
     &                 D5421*DSIN( XOMI+X2LI-G54)+
     &                 D5433*DSIN(-XOMI+X2LI-G54)
                XLDOT = XNI+XFACT
                XNDDT = D2201*DCOS(X2OMI+XLI-G22) + D2211*DCOS(XLI-G22)+
     &                  D3210*DCOS( XOMI+XLI-G32) +
     &                  D3222*DCOS(-XOMI+XLI-G32) +
     &                  D5220*DCOS( XOMI+XLI-G52) +
     &                  D5232*DCOS(-XOMI+XLI-G52) +
     &                  2.0D0*(D4410*DCOS(X2OMI+X2LI-G44) +
     &                  D4422*DCOS(X2LI-G44) +
     &                  D5421*DCOS( XOMI+X2LI-G54) +
     &                  D5433*DCOS(-XOMI+X2LI-G54))
                XNDDT = XNDDT*XLDOT
              ENDIF

* ------------------------------- INTEGRATOR --------------------------
              !  sgp4fix move end checks to end of routine
              IF (DABS(T-ATIME).ge.STEPP) THEN
                  IRET  = 0
                  IRETN = 381
                ELSE
                  FT    = T-ATIME
                  IRETN = 0
                ENDIF

              IF (IRETN.EQ.381) THEN
                  XLI   = XLI + XLDOT*DELT + XNDT*STEP2
                  XNI   = XNI + XNDT*DELT + XNDDT*STEP2
                  ATIME = ATIME + DELT
                ENDIF

              ENDDO

            XN = XNI + XNDT*FT  + XNDDT*FT*FT*0.5D0
            XL = XLI + XLDOT*FT + XNDT*FT*FT*0.5D0
            IF(IREZ .ne. 1) THEN
                Mm   = XL-2.0D0*nodem+2.0D0*THETA
                DNDT = XN-No
              ELSE
                Mm   = XL-nodem-Argpm+THETA
                DNDT = XN-No
              ENDIF

            XN = No + DNDT
          ENDIF

c        INCLUDE 'debug4.for' 

      RETURN
      END  ! end dspace


* -----------------------------------------------------------------------------
*
*                           SUBROUTINE INITL
*
*  this subroutine initializes the spg4 propagator. all the initialization is
*    consolidated here instead of having multiple loops inside other routines.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    ecco        - eccentricity                           0.0 - 1.0
*    epoch       - epoch time in days from jan 0, 1950. 0 hr
*    inclo       - inclination of satellite
*    no          - mean motion of satellite
*    satn        - satellite number
*
*  outputs       :
*    ainv        - 1.0 / a
*    ao          - semi major axis
*    con41       -
*    con42       - 1.0 - 5.0 cos(i)
*    cosio       - cosine of inclination
*    cosio2      - cosio squared
*    eccsq       - eccentricity squared
*    method      - flag for deep space                    'd', 'n'
*    omeosq      - 1.0 - ecco * ecco
*    posq        - semi-parameter squared
*    rp          - radius of perigee
*    rteosq      - square root of (1.0 - ecco*ecco)
*    sinio       - sine of inclination
*    gsto        - gst at time of observation               rad
*    no          - mean motion of satellite
*
*  locals        :
*    ak          -
*    d1          -
*    del         -
*    adel        -
*    po          -
*
*  coupling      :
*    getgravconst-
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
*------------------------------------------------------------------------------

      SUBROUTINE INITL( Satn , whichconst, Ecco  , EPOCH , Inclo , No,
     &         Method, AINV  , AO    , CON41 , CON42 , COSIO , COSIO2,
     &         Eccsq , OMEOSQ, POSQ  , rp    , RTEOSQ, SINIO ,
     &         GSTo, operationmode )
        IMPLICIT NONE
        CHARACTER Method, operationmode
        INTEGER Satn, whichconst
        REAL*8 Ecco  , EPOCH , Inclo , No   ,
     &         AINV  , AO    , CON41 , CON42 , COSIO , COSIO2, 
     &         Eccsq , OMEOSQ, POSQ  , rp    , RTEOSQ, SINIO , GSTo

        COMMON /DebugHelp/ Help
        CHARACTER Help

* -------------------------- Local Variables --------------------------
c        sgp4fix use old way of finding gst
        Integer ids70
        REAL*8 ts70, ds70, tfrac, c1, thgr70, fk5r, c1p2p, thgr, thgro

        REAL*8  RadPerDay, Temp, TUT1
        REAL*8  ak, d1, del, adel, po
        REAL*8  X2o3, J2, XKE, tumin, mu, radiusearthkm, j3, j4, j3oj2
        INCLUDE 'ASTMATH.CMN'

* ------------------------ WGS-72 EARTH CONSTANTS ---------------------
        X2o3   = 2.0D0/3.0D0
        ! sgp4fix identify constants and allow alternate values
        CALL getgravconst( whichconst, tumin, mu, radiusearthkm, xke,
     &       j2, j3, j4, j3oj2 )

* ----------------- CALCULATE AUXILLARY EPOCH QUANTITIES --------------
        Eccsq  = Ecco*Ecco
        OMEOSQ = 1.0D0 - Eccsq
        RTEOSQ = DSQRT(OMEOSQ)
        COSIO  = DCOS(Inclo)
        COSIO2 = COSIO*COSIO

* ---------------------- UN-KOZAI THE MEAN MOTION ---------------------
        AK   =  (XKE/No)**X2O3
        D1   =  0.75D0*J2* (3.0D0*COSIO2-1.0D0) / (RTEOSQ*OMEOSQ)
        DEL  =  D1/(AK*AK)
        ADEL =  AK * ( 1.0D0 - DEL*DEL - DEL*
     &                 (1.0D0/3.0D0 + 134.0D0*DEL*DEL / 81.0D0) )
        DEL  =  D1/(ADEL*ADEL)
        No   =  No/(1.0D0 + DEL)

        AO   =  (XKE/No)**X2O3
        SINIO=  DSIN(Inclo)
        PO   =  AO*OMEOSQ
        CON42=  1.0D0-5.0D0*COSIO2
        CON41=  -CON42-COSIO2-COSIO2
        AINV =  1.0D0/AO
        POSQ =  PO*PO
        rp   =  AO*(1.0D0-Ecco)
        METHOD = 'n'

* ----------------- CALCULATE GREENWICH LOCATION AT EPOCH -------------
c       sgp4fix modern approach to finding sidereal time
        IF (operationmode .ne. 'a') THEN
            RadPerDay  = twopi * 1.002737909350795D0  !6.30038809866574D0
            Temp = Epoch + 2433281.5D0
            TUT1= ( DINT(Temp-0.5D0) + 0.5D0 - 2451545.0D0 ) / 36525.0D0
            Gsto= 1.75336855923327D0 + 628.331970688841D0*TUT1
     &             + 6.77071394490334D-06*TUT1*TUT1
     &             - 4.50876723431868D-10*TUT1*TUT1*TUT1
     &             + RadPerDay*( Temp-0.5D0-DINT(Temp-0.5D0) )
          ELSE
            ! sgp4fix use old way of finding gst
            ! count integer number of days from 0 jan 1970
           TS70  = EPOCH-7305.0D0
           IDS70 = TS70 + 1.0D-8
           TFRAC = TS70-IDS70
            ! find greenwich location at epoch
           C1    = 1.72027916940703639D-2
           THGR70= 1.7321343856509374D0
            FK5R  = 5.07551419432269442D-15
           C1P2P = C1+TWOPI
           gsto  = THGR70+C1*IDS70+C1P2P*TFRAC+TS70*TS70*FK5R
         ENDIF
         
        ! ------------------------ Check quadrants ---------------------
        Gsto = DMOD( Gsto,TwoPi )
        IF ( Gsto .lt. 0.0D0 ) THEN
            Gsto= Gsto + TwoPi
          ENDIF

c      write(*,*) Satn,'  gst delta ', gsto-gsto1

c        INCLUDE 'debug5.for' 

      RETURN
      END  ! end initl


* -----------------------------------------------------------------------------
*
*                             SUBROUTINE SGP4INIT
*
*  This subroutine initializes variables for SGP4.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    satn        - satellite number
*    bstar       - sgp4 type drag coefficient              kg/m2er
*    ecco        - eccentricity
*    epoch       - epoch time in days from jan 0, 1950. 0 hr
*    argpo       - argument of perigee (output if ds)
*    inclo       - inclination
*    mo          - mean anomaly (output if ds)
*    no          - mean motion
*    nodeo      - right ascension of ascending node
*
*  outputs       :
*    satrec      - common block values for subsequent calls
*    return code - non-zero on error.
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*  locals        :
*    CNODM  , SNODM  , COSIM  , SINIM  , COSOMM , SINOMM
*    Cc1sq  , Cc2    , Cc3
*    Coef   , Coef1
*    cosio4      -
*    day         -
*    dndt        -
*    em          - eccentricity
*    emsq        - eccentricity squared
*    eeta        -
*    etasq       -
*    gam         -
*    argpm       - argument of perigee
*    ndem        -
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    perige      - perigee
*    pinvsq      -
*    psisq       -
*    qzms24      -
*    rtemsq      -
*    s1, s2, s3, s4, s5, s6, s7          -
*    sfour       -
*    ss1, ss2, ss3, ss4, ss5, ss6, ss7         -
*    sz1, sz2, sz3
*    sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33        -
*    tc          -
*    temp        -
*    temp1, temp2, temp3       -
*    tsi         -
*    xpidot      -
*    xhdot1      -
*    z1, z2, z3          -
*    z11, z12, z13, z21, z22, z23, z31, z32, z33         -
*
*  coupling      :
*    getgravconst-
*    initl       -
*    dscom       -
*    dpper       -
*    dsinit      -
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
* ---------------------------------------------------------------------------- }

      SUBROUTINE SGP4Init ( whichconst,
     &                      Satn,   xBStar, xEcco,  Epoch, xArgpo,
     &                      xInclo, xMo,    xNo,    xnodeo, Error )
        IMPLICIT NONE
        INTEGER Satn, error, whichconst
        REAL*8  xBStar, xEcco, Epoch, xArgpo, xInclo, xMo, xNo, xnodeo
        REAL*8 T, r(3), v(3)

        INCLUDE 'sgp4.cmn'

        COMMON /DebugHelp/ Help
        CHARACTER Help

* -------------------------- Local Variables --------------------------

        REAL*8  Ao,ainv,con42,cosio,sinio,cosio2,Eccsq,omeosq,
     &          posq,rp,rteosq, CNODM , SNODM , COSIM , SINIM , COSOMM,
     &          SINOMM, Cc1sq ,
     &          Cc2   , Cc3   , Coef  , Coef1 , Cosio4, DAY   , Dndt  ,
     &          Eccm  , EMSQ  , Eeta  , Etasq , GAM   , Argpm , nodem,
     &          Inclm , Mm  , Xn    , Perige, Pinvsq, Psisq , Qzms24,
     &          RTEMSQ, S1    , S2    , S3    , S4    , S5    , S6    ,
     &          S7    , SFour , SS1   , SS2   , SS3   , SS4   , SS5   ,
     &          SS6   , SS7   , SZ1   , SZ2   , SZ3   , SZ11  , SZ12  ,
     &          SZ13  , SZ21  , SZ22  , SZ23  , SZ31  , SZ32  , SZ33  ,
     &          Tc    , Temp  , Temp1 , Temp2 , Temp3 , Tsi   , XPIDOT,
     &          Xhdot1, Z1    , Z2    , Z3    , Z11   , Z12   , Z13   ,
     &          Z21   , Z22   , Z23   , Z31   , Z32   , Z33 
        REAL*8  qzms2t, SS, mu, RadiusEarthKm, J2, j3oJ2,J4,X2o3,
     &          temp4, j3, xke, tumin
        INCLUDE 'ASTMATH.CMN'

* ---------------------------- INITIALIZATION -------------------------
        method = 'n'
c       clear sgp4 flag
        Error = 0

c      sgp4fix - note the following variables are also passed directly via sgp4 common. 
c      it is possible to streamline the sgp4init call by deleting the "x"
c      variables, but the user would need to set the common values first. we
c      include the additional assignment in case twoline2rv is not used. 
 
        bstar  = xbstar
        ecco   = xecco
        argpo  = xargpo
        inclo  = xinclo
        mo     = xmo
        no     = xno
        nodeo  = xnodeo

        ! sgp4fix identify constants and allow alternate values
        CALL getgravconst( whichconst, tumin, mu, radiusearthkm, xke,
     &       j2, j3, j4, j3oj2 )

        SS     = 78.0D0/RadiusEarthKm + 1.0D0
        QZMS2T = ((120.0D0-78.0D0)/RadiusEarthKm) ** 4
        X2o3   =  2.0D0 / 3.0D0
c     sgp4fix divisor for divide by zero check on inclination
c     the old check used 1.0D0 + cos(pi-1.0D-9), but then compared it to
c     1.5D-12, so the threshold was changed to 1.5D-12 for consistency
        temp4    =   1.5D-12

        Init = 'y'
        T = 0.0D0

        CALL INITL( Satn , whichconst, Ecco  , EPOCH , Inclo , No,
     &     Method, AINV  , AO    , CON41 , CON42 , COSIO , COSIO2,
     &     Eccsq , OMEOSQ, POSQ  , rp    , RTEOSQ, SINIO ,
     &     GSTo, Opsmode )

        IF(rp .lt. 1.0D0) THEN
c            Write(*,*) '# *** SATN',Satn,' EPOCH ELTS SUB-ORBITAL *** '
            Error = 5
          ENDIF

        IF(OMEOSQ .ge. 0.0D0 .OR. No .ge. 0.0D0) THEN
            ISIMP = 0
            IF (rp .lt. (220.0D0/RadiusEarthKm+1.0D0)) THEN
                ISIMP = 1
              ENDIF
            SFour  = SS
            QZMS24 = QZMS2T
            PERIGE = (rp-1.0D0)*RadiusEarthKm

* ----------- For perigees below 156 km, S and Qoms2t are altered -----
            IF(PERIGE .lt. 156.0D0) THEN
                SFour = PERIGE-78.0D0
                IF(PERIGE .le. 98.0D0) THEN
                    SFour = 20.0D0
                  ENDIF
                QZMS24 = ( (120.0D0-SFour)/RadiusEarthKm )**4
                SFour  = SFour/RadiusEarthKm + 1.0D0
              ENDIF
            PINVSQ = 1.0D0/POSQ

            TSI    = 1.0D0/(AO-SFour)
            ETA    = AO*Ecco*TSI
            ETASQ  = ETA*ETA
            EETA   = Ecco*ETA
            PSISQ  = DABS(1.0D0-ETASQ)
            COEF   = QZMS24*TSI**4
            COEF1  = COEF/PSISQ**3.5D0
            CC2    = COEF1*No* (AO* (1.0D0+1.5D0*ETASQ+EETA*
     &               (4.0D0+ETASQ) )+0.375D0*
     &         J2*TSI/PSISQ*CON41*(8.0D0+3.0D0*ETASQ*(8.0D0+ETASQ)))
            CC1    = BSTAR*CC2
            CC3    = 0.0D0
            IF(Ecco .GT. 1.0D-4) THEN
                CC3 = -2.0D0*COEF*TSI*J3OJ2*No*SINIO/Ecco
              ENDIF
            X1MTH2 = 1.0D0-COSIO2
            CC4    = 2.0D0*No*COEF1*AO*OMEOSQ*(ETA*(2.0D0+0.5D0*ETASQ)
     &              +Ecco*(0.5D0 + 2.0D0*ETASQ) - J2*TSI / (AO*PSISQ)*
     &              (-3.0D0*CON41*(1.0D0-2.0D0*
     &       EETA+ETASQ*(1.5D0-0.5D0*EETA))+0.75D0*X1MTH2*(2.0D0*ETASQ
     &       -EETA*(1.0D0+ETASQ))*DCOS(2.0D0*Argpo)))
            CC5    = 2.0D0*COEF1*AO*OMEOSQ* (1.0D0 + 2.75D0*
     &               (ETASQ + EETA) + EETA*ETASQ )
            COSIO4 = COSIO2*COSIO2
            TEMP1  = 1.5D0*J2*PINVSQ*No
            TEMP2  = 0.5D0*TEMP1*J2*PINVSQ
            TEMP3  = -0.46875D0*J4*PINVSQ*PINVSQ*No
            MDot   = No + 0.5D0*TEMP1*RTEOSQ*CON41 + 0.0625D0*TEMP2*
     &               RTEOSQ*(13.0D0 - 78.0D0*COSIO2 + 137.0D0*COSIO4)
            ArgpDot= -0.5D0*TEMP1*CON42 + 0.0625D0*TEMP2*
     &               (7.0D0 - 114.0D0*COSIO2 +
     &        395.0D0*COSIO4)+TEMP3*(3.0D0-36.0D0*COSIO2+49.0D0*COSIO4)
            XHDOT1 = -TEMP1*COSIO
            nodeDot = XHDOT1+(0.5D0*TEMP2*(4.0D0-19.0D0*COSIO2)+
     &                 2.0D0*TEMP3*(3.0D0 - 7.0D0*COSIO2))*COSIO
            XPIDOT = ArgpDot+nodeDot
            OMGCOF = BSTAR*CC3*DCOS(Argpo)
            XMCOF  = 0.0D0
            IF(Ecco .GT. 1.0D-4) THEN
                XMCOF = -X2O3*COEF*BSTAR/EETA
              ENDIF
            XNODCF = 3.5D0*OMEOSQ*XHDOT1*CC1
            T2COF  = 1.5D0*CC1
c           sgp4fix for divide by zero with xinco = 180 deg
            if (dabs(cosio+1.0).gt. 1.5d-12) THEN
                XLCOF  = -0.25D0*J3OJ2*SINIO*
     &                   (3.0D0+5.0D0*COSIO)/(1.0D0+COSIO)
              else
                XLCOF  = -0.25D0*J3OJ2*SINIO*
     &                   (3.0D0+5.0D0*COSIO)/temp4
              ENDIF
            AYCOF  = -0.5D0*J3OJ2*SINIO
            DELMO  = (1.0D0+ETA*DCOS(Mo))**3
            SINMAO = DSIN(Mo)
            X7THM1 = 7.0D0*COSIO2-1.0D0

* ------------------------ Deep Space Initialization ------------------
            IF ((TWOPI/No) .ge. 225.0D0) THEN
                METHOD = 'd'
                ISIMP  = 1
                TC     = 0.0D0
                Inclm  = Inclo
                CALL DSCOM( EPOCH     , Ecco  , Argpo , Tc    , Inclo ,
     &                  nodeo, No    ,
     &                  SNODM , CNODM , SINIM , COSIM , SINOMM, COSOMM,
     &                  DAY   , E3    , Ee2   , Eccm  , EMSQ  , GAM   ,
     &                  Peo   , Pgho  , Pho   , PInco , Plo   ,
     &                  RTemSq, Se2   , Se3   , Sgh2  , Sgh3  , Sgh4  ,
     &                  Sh2   , Sh3   , Si2   , Si3   , Sl2   , Sl3   ,
     &                  Sl4   , S1    , S2    , S3    , S4    , S5    ,
     &                  S6    , S7    , SS1   , SS2   , SS3   , SS4   ,
     &                  SS5   , SS6   , SS7   , SZ1   , SZ2   , SZ3   ,
     &                  SZ11  , SZ12  , SZ13  , SZ21  , SZ22  , SZ23  ,
     &                  SZ31  , SZ32  , SZ33  , Xgh2  , Xgh3  , Xgh4  ,
     &                  Xh2   , Xh3   , Xi2   , Xi3   , Xl2   , Xl3   ,
     &                  Xl4   , Xn    , Z1    , Z2    , Z3    , Z11   ,
     &                  Z12   , Z13   , Z21   , Z22   , Z23   , Z31   ,
     &                  Z32   , Z33   , Zmol  , Zmos )
                CALL DPPER( e3, ee2   , peo   , pgho  , pho   , pinco ,
     &                  plo   , se2   , se3   , sgh2  , sgh3  , sgh4  ,
     &                  sh2   , sh3   , si2   , si3   , sl2   , sl3   ,
     &                  sl4   , T     , xgh2  , xgh3  , xgh4  , xh2   ,
     &                  xh3   , xi2   , xi3   , xl2   , xl3   , xl4   ,
     &                  zmol  , zmos  , Inclm , init  ,
     &                  Ecco  , Inclo , nodeo, Argpo , Mo, Opsmode )

                Argpm  = 0.0D0 ! add for DS to work initial
                nodem  = 0.0D0
                Mm     = 0.0D0

                CALL DSINIT( whichconst,
     &                   Cosim ,Emsq, Argpo, S1    , S2    , S3    ,
     &                   S4    , S5    , Sinim , Ss1   , Ss2   , Ss3   ,
     &                   Ss4   , Ss5   , Sz1   , Sz3   , Sz11  , Sz13  ,
     &                   Sz21  , Sz23  , Sz31  , Sz33  , T     , Tc    ,
     &                   GSTo  , Mo    , MDot  , No    ,nodeo,nodeDot,
     &                   XPIDOT, Z1    , Z3    , Z11   , Z13   , Z21   ,
     &                   Z23   , Z31   , Z33   , ecco  , eccsq,
     &                   Eccm  , Argpm , Inclm , Mm    , Xn    , nodem,
     &                   IREZ  , Atime , D2201 , D2211 , D3210 , D3222 ,
     &                   D4410 , D4422 , D5220 , D5232 , D5421 , D5433 ,
     &                   Dedt  , Didt  , DMDT  , DNDT  , DNODT , DOMDT ,
     &                   Del1  , Del2  , Del3  , Xfact , Xlamo , Xli   ,
     &                   Xni )
            ENDIF

* ------------ Set variables if not deep space or rp < 220 -------------
            IF (ISIMP .ne. 1) THEN
                CC1SQ = CC1*CC1
                D2    = 4.0D0*AO*TSI*CC1SQ
                TEMP  = D2*TSI*CC1 / 3.0D0
                D3    = (17.0D0*AO + SFour) * TEMP
                D4    = 0.5D0*TEMP*AO*TSI*
     &                  (221.0D0*AO + 31.0D0*SFour)*CC1
                T3COF = D2 + 2.0D0*CC1SQ
                T4COF = 0.25D0* (3.0D0*D3+CC1*(12.0D0*D2+10.0D0*CC1SQ) )
                T5COF = 0.2D0* (3.0D0*D4 + 12.0D0*CC1*D3 + 6.0D0*D2*D2 +
     &                  15.0D0*CC1SQ* (2.0D0*D2 + CC1SQ) )
              ENDIF

          ENDIF ! ------ if nodeo and No are gtr 0

      init = 'n'

      CALL SGP4(whichconst, 0.0D0, r, v, error)

c        INCLUDE 'debug6.for'

      RETURN
      END  ! end sgp4init


* -----------------------------------------------------------------------------
*
*                             SUBROUTINE SGP4
*
*  this procedure is the sgp4 prediction model from space command. this is an
*    updated and combined version of sgp4 and sdp4, which were originally
*    published separately in spacetrack report #3. this version follows the
*    methodology from the aiaa paper (2006) describing the history and
*    development of the code.
*
*  author        : david vallado                  719-573-2600   28 jun 2005
*
*  inputs        :
*    satrec	 - initialised structure from sgp4init() call.
*    tsince	 - time eince epoch (minutes)
*
*  outputs       :
*    r           - position vector                     km
*    v           - velocity                            km/sec
*  return code - non-zero on error.
*                   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*                   2 - mean motion less than 0.0
*                   3 - pert elements, ecc < 0.0  or  ecc > 1.0
*                   4 - semi-latus rectum < 0.0
*                   5 - epoch elements are sub-orbital
*                   6 - satellite has decayed
*
*  locals        :
*    am          -
*    axnl, aynl        -
*    betal       -
*    COSIM   , SINIM   , COSOMM  , SINOMM  , Cnod    , Snod    , Cos2u   ,
*    Sin2u   , Coseo1  , Sineo1  , Cosi    , Sini    , Cosip   , Sinip   ,
*    Cosisq  , Cossu   , Sinsu   , Cosu    , Sinu
*    Delm        -
*    Delomg      -
*    Dndt        -
*    Eccm        -
*    EMSQ        -
*    Ecose       -
*    El2         -
*    Eo1         -
*    Eccp        -
*    Esine       -
*    Argpm       -
*    Argpp       -
*    Omgadf      -
*    Pl          -
*    R           -
*    RTEMSQ      -
*    Rdotl       -
*    Rl          -
*    Rvdot       -
*    Rvdotl      -
*    Su          -
*    T2  , T3   , T4    , Tc
*    Tem5, Temp , Temp1 , Temp2  , Tempa  , Tempe  , Templ
*    U   , Ux   , Uy    , Uz     , Vx     , Vy     , Vz
*    inclm       - inclination
*    mm          - mean anomaly
*    nm          - mean motion
*    nodem       - longi of ascending node
*    xinc        -
*    xincp       -
*    xl          -
*    xlm         -
*    mp          -
*    xmdf        -
*    xmx         -
*    xmy         -
*    nodedf     -
*    xnode       -
*    nodep      -
*    np          -
*
*  coupling      :
*    getgravconst-
*    dpper
*    dpspace
*
*  references    :
*    hoots, roehrich, norad spacetrack report #3 1980
*    hoots, norad spacetrack report #6 1986
*    hoots, schumacher and glover 2004
*    vallado, crawford, hujsak, kelso  2006
*------------------------------------------------------------------------------

      SUBROUTINE SGP4 ( whichconst, T, r, v, Error )
        IMPLICIT NONE
        INTEGER  Error, whichconst
        REAL*8   T, r(3), v(3)

        INCLUDE 'sgp4.cmn'

* -------------------------- Local Variables --------------------------
        REAL*8 AM    , Axnl  , Aynl  , Betal , COSIM , Cnod  ,
     &         Cos2u , Coseo1, Cosi  , Cosip , Cosisq, Cossu , Cosu  ,
     &         Delm  , Delomg, Eccm  , EMSQ  , Ecose , El2   , Eo1   ,
     &         Eccp  , Esine , Argpm , Argpp , Omgadf, Pl    ,
     &         Rdotl , Rl    , Rvdot , Rvdotl, SINIM ,
     &         Sin2u , Sineo1, Sini  , Sinip , Sinsu , Sinu  ,
     &         Snod  , Su    , T2    , T3    , T4    , Tem5  , Temp  ,
     &         Temp1 , Temp2 , Tempa , Tempe , Templ , U     , Ux    ,
     &         Uy    , Uz    , Vx    , Vy    , Vz    , Inclm , Mm  ,
     &         XN    , nodem , Xinc  , Xincp , Xl    , Xlm   , Mp  ,
     &         Xmdf  , Xmx   , Xmy   , Xnoddf, Xnode , nodep,
     &         Tc    , Dndt

        REAL*8 X2O3, J2,J3,XKE,J3OJ2, mr,mv,
     &         mu, RadiusEarthkm, VKmPerSec, temp4, tumin, j4
	INTEGER iter

        COMMON /DebugHelp/ Help
        CHARACTER Help
        INCLUDE 'ASTMATH.CMN'

* ------------------------ WGS-72 EARTH CONSTANTS ---------------------
* ---------------------- SET MATHEMATICAL CONSTANTS -------------------
      X2O3   = 2.0D0/3.0D0

c     Keep compiler ok for warnings on uninitialized variables
      mr = 0.0D0
      Coseo1 = 1.0D0
      Sineo1 = 0.0D0

      ! sgp4fix identify constants and allow alternate values
        CALL getgravconst( whichconst, tumin, mu, radiusearthkm, xke,
     &       j2, j3, j4, j3oj2 )
c     sgp4fix divisor for divide by zero check on inclination
c     the old check used 1.0D0 + cos(pi-1.0D-9), but then compared it to
c     1.5D-12, so the threshold was changed to 1.5D-12 for consistency
      temp4    =   1.5D-12
      VKmPerSec     =  RadiusEarthKm * xke/60.0D0

* ------------------------- CLEAR SGP4 ERROR FLAG ---------------------
      Error = 0

* ----------- UPDATE FOR SECULAR GRAVITY AND ATMOSPHERIC DRAG ---------
      XMDF   = Mo + MDot*T
      OMGADF = Argpo + ArgpDot*T
      XNODDF = nodeo + nodeDot*T
      Argpm  = OMGADF
      Mm     = XMDF
      T2     = T*T
      nodem  = XNODDF + XNODCF*T2
      TEMPA  = 1.0D0 - CC1*T
      TEMPE  = BSTAR*CC4*T
      TEMPL  = T2COF*T2
      IF (ISIMP .ne. 1) THEN
          DELOMG = OMGCOF*T
          DELM   = XMCOF*(( 1.0D0+ETA*DCOS(XMDF) )**3-DELMO)
          TEMP   = DELOMG + DELM
          Mm     = XMDF + TEMP
          Argpm  = OMGADF - TEMP
          T3     = T2*T
          T4     = T3*T
          TEMPA  = TEMPA - D2*T2 - D3*T3 - D4*T4
          TEMPE  = TEMPE + BSTAR*CC5*(DSIN(Mm) - SINMAO)
          TEMPL  = TEMPL + T3COF*T3 + T4*(T4COF + T*T5COF)
        ENDIF
      XN    = No
      Eccm  = Ecco
      Inclm = Inclo
      IF(METHOD .EQ. 'd') THEN
          TC     = T
          CALL DSPACE( IRez  , D2201 , D2211 , D3210 , D3222 , D4410 ,
     &                 D4422 , D5220 , D5232 , D5421 , D5433 , Dedt  ,
     &                 Del1  , Del2  , Del3  , Didt  , Dmdt  , Dnodt ,
     &                 Domdt , Argpo , ArgpDot, T    , TC    , GSTo ,
     &                 Xfact , Xlamo , No   ,
     &                 Atime , Eccm  , Argpm, Inclm , Xli   , Mm  ,
     &                 XNi   , nodem, Dndt  , XN  )
        ENDIF

c     mean motion less than 0.0
      IF(XN .LE. 0.0D0) THEN
          Error = 2
        ENDIF
      AM = (XKE/XN)**X2O3*TEMPA**2
      XN = XKE/AM**1.5D0
      Eccm = Eccm-TEMPE
c   fix tolerance for error recognition
      IF (Eccm .GE. 1.0D0 .or. Eccm.lt.-0.001D0 .or. AM .lt. 0.95) THEN
c	  write(6,*) '# Error 1, Eccm = ',  Eccm, ' AM = ', AM
          Error = 1
        ENDIF
c   sgp4fix change test condition for eccentricity   
      IF (Eccm .lt. 1.0D-6) Eccm = 1.0D-6
      Mm     = Mm+No*TEMPL
      XLM    = Mm+Argpm+nodem
      EMSQ   = Eccm*Eccm
      TEMP   = 1.0D0 - EMSQ
      nodem  = DMOD(nodem,TwoPi)
      Argpm  = DMOD(Argpm,TwoPi)
      XLM    = DMOD(XLM,TwoPi)
      Mm     = DMOD(XLM - Argpm - nodem,TwoPi)

* --------------------- COMPUTE EXTRA MEAN QUANTITIES -----------------
      SINIM  = DSIN(Inclm)
      COSIM  = DCOS(Inclm)

* ------------------------ ADD LUNAR-SOLAR PERIODICS ------------------
      Eccp   = Eccm
      XINCP  = Inclm
      Argpp  = Argpm
      nodep = nodem
      Mp     = Mm
      SINIP  = SINIM
      COSIP  = COSIM
      IF(METHOD .EQ. 'd') THEN
          CALL DPPER( e3    , ee2   , peo   , pgho  , pho   , pinco ,
     &                plo   , se2   , se3   , sgh2  , sgh3  , sgh4  ,
     &                sh2   , sh3   , si2   , si3   , sl2   , sl3   ,
     &                sl4   , T     , xgh2  , xgh3  , xgh4  , xh2   ,
     &                xh3   , xi2   , xi3   , xl2   , xl3   , xl4   ,
     &                zmol  , zmos  , Inclo , 'n'   ,
     &                Eccp  , XIncp , nodep, Argpp, Mp, Opsmode )
          IF(XINCP .lt. 0.0D0) THEN
              XINCP  = -XINCP
              nodep  = nodep + PI
              Argpp  = Argpp - PI
            ENDIF
          IF(Eccp .lt. 0.0D0 .OR. Eccp .GT. 1.0D0) THEN
              Error = 3
            ENDIF
        ENDIF

* ------------------------ LONG PERIOD PERIODICS ----------------------
      IF(METHOD .EQ. 'd') THEN
          SINIP =  DSIN(XINCP)
          COSIP =  DCOS(XINCP)
          AYCOF = -0.5D0*J3OJ2*SINIP
c         sgp4fix for divide by zero with xincp = 180 deg
          if (dabs(cosip+1.0).gt. 1.5d-12) THEN
              XLCOF  = -0.25D0*J3OJ2*SINIP*
     &                 (3.0D0+5.0D0*COSIP)/(1.0D0+COSIP)
            else
              XLCOF  = -0.25D0*J3OJ2*SINIP*
     &                 (3.0D0+5.0D0*COSIP)/temp4
            ENDIF
        ENDIF
      AXNL = Eccp*DCOS(Argpp)
      TEMP = 1.0D0 / (AM*(1.0D0-Eccp*Eccp))
      AYNL = Eccp*DSIN(Argpp) + TEMP*AYCOF
      XL   = Mp + Argpp + nodep + TEMP*XLCOF*AXNL

* ------------------------- SOLVE KEPLER'S EQUATION -------------------
      U    = DMOD(XL-nodep,TwoPi)
      EO1  = U
      ITER=0
c   sgp4fix for kepler iteration
c   the following iteration needs better limits on corrections
      Temp = 9999.9D0
      DO WHILE ((Temp.ge.1.0D-12).and.(ITER.lt.10))
          ITER=ITER+1
          SINEO1= DSIN(EO1)
          COSEO1= DCOS(EO1)
          TEM5  = 1.0D0 - COSEO1*AXNL - SINEO1*AYNL
          TEM5  = (U - AYNL*COSEO1 + AXNL*SINEO1 - EO1) / TEM5
          Temp  = DABS(Tem5)
          IF(Temp.gt.1.0D0) Tem5=Tem5/Temp ! Stop excessive correction
          EO1   = EO1+TEM5
        ENDDO

* ----------------- SHORT PERIOD PRELIMINARY QUANTITIES ---------------
      ECOSE = AXNL*COSEO1+AYNL*SINEO1
      ESINE = AXNL*SINEO1-AYNL*COSEO1
      EL2   = AXNL*AXNL+AYNL*AYNL
      PL    = AM*(1.0D0-EL2)
c     semi-latus rectum < 0.0
      IF ( PL .lt. 0.0D0 ) THEN
          Error = 4
        ELSE
          RL    = AM*(1.0D0-ECOSE)
          RDOTL = DSQRT(AM)*ESINE/RL
          RVDOTL= DSQRT(PL)/RL
          BETAL = DSQRT(1.0D0-EL2)
          TEMP  = ESINE/(1.0D0+BETAL)
          SINU  = AM/RL*(SINEO1-AYNL-AXNL*TEMP)
          COSU  = AM/RL*(COSEO1-AXNL+AYNL*TEMP)
          SU    = DATAN2(SINU,COSU)
          SIN2U = (COSU+COSU)*SINU
          COS2U = 1.0D0-2.0D0*SINU*SINU
          TEMP  = 1.0D0/PL
          TEMP1 = 0.5D0*J2*TEMP
          TEMP2 = TEMP1*TEMP

* ------------------ UPDATE FOR SHORT PERIOD PERIODICS ----------------
          IF(METHOD .EQ. 'd') THEN
              COSISQ = COSIP*COSIP
              CON41  = 3.0D0*COSISQ - 1.0D0
              X1MTH2 = 1.0D0 - COSISQ
              X7THM1 = 7.0D0*COSISQ - 1.0D0
            ENDIF
          mr   = RL*(1.0D0 - 1.5D0*TEMP2*BETAL*CON41) +
     &           0.5D0*TEMP1*X1MTH2*COS2U
          SU   = SU - 0.25D0*TEMP2*X7THM1*SIN2U
          XNODE= nodep + 1.5D0*TEMP2*COSIP*SIN2U
          XINC = XINCP + 1.5D0*TEMP2*COSIP*SINIP*COS2U
          mv   = RDOTL - XN*TEMP1*X1MTH2*SIN2U / XKE
          RVDOT= RVDOTL + XN*TEMP1* (X1MTH2*COS2U+1.5D0*CON41) / XKE

* ------------------------- ORIENTATION VECTORS -----------------------
          SINSU=  DSIN(SU)
          COSSU=  DCOS(SU)
          SNOD =  DSIN(XNODE)
          CNOD =  DCOS(XNODE)
          SINI =  DSIN(XINC)
          COSI =  DCOS(XINC)
          XMX  = -SNOD*COSI
          XMY  =  CNOD*COSI
          UX   =  XMX*SINSU + CNOD*COSSU
          UY   =  XMY*SINSU + SNOD*COSSU
          UZ   =  SINI*SINSU
          VX   =  XMX*COSSU - CNOD*SINSU
          VY   =  XMY*COSSU - SNOD*SINSU
          VZ   =  SINI*COSSU

* ----------------------- POSITION AND VELOCITY -----------------------
          r(1) = mr*UX * RadiusEarthkm
          r(2) = mr*UY * RadiusEarthkm
          r(3) = mr*UZ * RadiusEarthkm
          v(1) = (mv*UX + RVDOT*VX) * VKmPerSec
          v(2) = (mv*UY + RVDOT*VY) * VKmPerSec
          v(3) = (mv*UZ + RVDOT*VZ) * VKmPerSec
        ENDIF

* --------------------------- ERROR PROCESSING ------------------------
c     sgp4fix for decaying satellites
      if (mr .lt. 1.0D0) THEN
c          write(*,*) '# decay condition ',mr
          error = 6
        ENDIF

c        INCLUDE 'debug7.for'

      RETURN
      END  ! end sgp4

* -----------------------------------------------------------------------------
*
*                           FUNCTION GSTIME
*
*  This function finds the Greenwich SIDEREAL time.  Notice just the INTEGER
*    part of the Julian Date is used for the Julian centuries calculation.
*    We use radper Solar day because we're multiplying by 0-24 solar hours.
*
*  Author        : David Vallado                  719-573-2600    1 Mar 2001
*
*  Inputs          Description                    Range / Units
*    JD          - Julian Date                    days from 4713 BC
*
*  OutPuts       :
*    GSTIME      - Greenwich SIDEREAL Time        0 to 2Pi rad
*
*  Locals        :
*    Temp        - Temporary variable for reals   rad
*    TUT1        - Julian Centuries from the
*                  Jan 1, 2000 12 h epoch (UT1)
*
*  Coupling      :
*
*  References    :
*    Vallado       2007, 194, Eq 3-45
* -----------------------------------------------------------------------------

      REAL*8 FUNCTION GSTIME ( JD )
        IMPLICIT NONE
        REAL*8 JD
* ----------------------------  Locals  -------------------------------
        REAL*8 Temp, TUT1

        INCLUDE 'astmath.cmn'

        ! --------------------  Implementation   ----------------------

        TUT1= ( JD - 2451545.0D0 ) / 36525.0D0
        Temp= - 6.2D-6*TUT1*TUT1*TUT1
     &        + 0.093104D0*TUT1*TUT1
     &        + (876600.0D0*3600.0D0 + 8640184.812866D0)*TUT1
     &        + 67310.54841D0
        Temp= DMOD( Temp*Deg2Rad/240.0D0,TwoPi ) ! 360/86400 = 1/240, to deg, to rad

        ! ------------------------ Check quadrants --------------------
        IF ( Temp .lt. 0.0D0 ) THEN
            Temp= Temp + TwoPi
          ENDIF

        GSTIME= Temp

      RETURN
      END  ! end gstime


* -----------------------------------------------------------------------------
*
*                           function getgravconst
*
*  this function gets constants for the propagator. note that mu is identified to
*    facilitiate comparisons with newer models.
*
*  author        : david vallado                  719-573-2600   21 jul 2006
*
*  inputs        :
*    whichconst  - which set of constants to use  721, 72, 84
*
*  outputs       :
*    tumin       - minutes in one time unit
*    mu          - earth gravitational parameter
*    radiusearthkm - radius of the earth in km
*    xke         - reciprocal of tumin
*    j2, j3, j4  - un-normalized zonal harmonic values
*    j3oj2       - j3 divided by j2
*
*  locals        :
*
*  coupling      :
*
*  references    :
*    norad spacetrack report #3
*    vallado, crawford, hujsak, kelso  2006
*  ---------------------------------------------------------------------------- 

       SUBROUTINE getgravconst ( whichconst, tumin, mu, 
     &            radiusearthkm, xke, j2, j3, j4, j3oj2 )
       IMPLICIT NONE     
       REAL*8 radiusearthkm, xke, j2, j3, j4, j3oj2, mu, tumin
       INTEGER whichconst

       if (whichconst.eq.721) THEN
           ! -- wgs-72 low precision str#3 constants --
           radiusearthkm = 6378.135D0     ! km
           xke    = 0.0743669161D0
           mu     = 398600.79964D0            ! in km3 / s2
           tumin  = 1.0D0 / xke
           j2     =   0.001082616D0
           j3     =  -0.00000253881D0
           j4     =  -0.00000165597D0
           j3oj2  =  j3 / j2
         ENDIF
       if (whichconst.eq.72) THEN
           ! ------------ wgs-72 constants ------------
           mu     = 398600.8D0            ! in km3 / s2
           radiusearthkm = 6378.135D0     ! km
           xke    = 60.0D0 / dsqrt(radiusearthkm**3/mu)
           tumin  = 1.0D0 / xke
           j2     =   0.001082616D0
           j3     =  -0.00000253881D0
           j4     =  -0.00000165597D0
           j3oj2  =  j3 / j2
         ENDIF  
       if (whichconst.eq.84) THEN
           ! ------------ wgs-84 constants ------------
           mu     = 398600.5D0            ! in km3 / s2
           radiusearthkm = 6378.137D0     ! km
           xke    = 60.0D0 / dsqrt(radiusearthkm**3/mu)
           tumin  = 1.0D0 / xke
           j2     =   0.00108262998905D0
           j3     =  -0.00000253215306D0
           j4     =  -0.00000161098761D0
           j3oj2  =  j3 / j2
         ENDIF

       RETURN
       END  !  SUBROUTINE getgravconst



c=======================================================================
c=======================================================================
c=======================================================================
c=======================================================================
c=======================================================================
*  Files         :
*    Unit 110     - input elm file  input file for element sets
*    Unit 111     - sgp4test.out    output file
*    Unit 114     - sgp4test.dbg    debug output file
*    Unit 115     - sgp4rec.bak     temporary file of record for 2 line element sets
*
*  Uses object and include files:
*    Astmath.cmn,
*    Sgp4.cmn,
*    Sgp4ext,
*    Sgp4io,
*    Sgp4unit

      subroutine sgp4_init(ro, vo, long, lat)
        IMPLICIT NONE
        Character typerun, typeinput
        Character*12 InFileName

        Character*3 MonStr,Monthtitle(12)
        Integer Code, NumSats, TotalNumSats, k, error, whichconst
        Real*8 ro(3),vo(3), Tmfe

        REAL*8 p, ecc, incl, node, argp, nu, m,arglat,truelon,lonper

* ----------------------------  Locals  -------------------------------
        REAL*8 J2,TwoPi,Rad,mu, RadiusEarthKm,VKmPerSec, xke, thetaout, 
     &         de2ra, xpdotp, T, sec, JD, pi, j3, j4, j3oj2, tumin, theta, theta0g, long, lat, gstime, xyz(3), lla(3)
        INTEGER i,j, Year,yr,mon,day,hr,min

        save


        INCLUDE 'sgp4.cmn'

        COMMON /DebugHelp/ Help
        CHARACTER Help
        Help = 'N'

* ------------------------  Implementation   --------------------------
c
c-- improved sgp4 operation
        Opsmode    = 'i' ! improved sgp4 operation
c-- Select Manual as type of run
        typerun    = 'M'
c-- Select minutes from Epoch as the approach
        typeinput  = 'M'
c  'Input whichconst - 721, 72, 84'
        whichconst    = 84
        pi            =    4.0D0 * datan(1.0D0)  ! 3.14159265358979D0
        TwoPi         =    2.0D0 * pi    ! 6.28318530717959D0
        Rad           =   180.0D0 / pi   ! 57.29577951308230D0
        DE2RA         =    pi / 180.0D0  ! 0.01745329251994330D0
        xpdotp        =  1440.0 / (2.0 *pi)  ! 229.1831180523293D0
        ! sgp4fix identify constants and allow alternate values
        CALL getgravconst( whichconst, tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 )
        VKmPerSec     =  RadiusEarthKm * xke/60.0D0
        MonthTitle( 1)= 'Jan'
        MonthTitle( 2)= 'Feb'
        MonthTitle( 3)= 'Mar'
        MonthTitle( 4)= 'Apr'
        MonthTitle( 5)= 'May'
        MonthTitle( 6)= 'Jun'
        MonthTitle( 7)= 'Jul'
        MonthTitle( 8)= 'Aug'
        MonthTitle( 9)= 'Sep'
        MonthTitle(10)= 'Oct'
        MonthTitle(11)= 'Nov'
        MonthTitle(12)= 'Dec'
c
        InFileName = 'hrpt.tle'
        OPEN(110,FILE = InFileName ,STATUS='OLD',ACCESS = 'SEQUENTIAL' )
c-- 111 output file
        OPEN(111,FILE = 'tfor.out' ,STATUS='UNKNOWN',ACCESS = 'SEQUENTIAL' )
        OPEN(114,FILE = 'sgp4test.dbg' ,STATUS='UNKNOWN',ACCESS = 'SEQUENTIAL' )
c-- 115 temporary file of record for 2 line element sets ---
        OPEN(115,FILE = 'Sgp4Rec.bak', ACCESS = 'DIRECT',FORM = 'UNFORMATTED', RECL = 1100, STATUS = 'UNKNOWN' )

        ! ----------------- Test simple propagation -------------------
        NumSats = 0
        Numsats = NumSats + 1
        CALL TwoLine2RVSGP4 ( NumSats,typerun,typeinput,whichconst,Code )
c-- Initialisation done
        if (Code.eq.999) then
          stop ' ** TLE Read error - check hrpt.tle'
        endif
        T = 0.0D0
        CALL SGP4 ( whichconst, T, Ro, Vo, Error )
        return
c
      entry sgp4_run(Tmfe, ro, vo, long, lat, thetaout)
c-- now initialize time variables
      T  = Tmfe
      CALL SGP4 ( whichconst, T, Ro, Vo, Error )
c
      IF (Error .gt. 0 ) Write(*,*) '# Error in SGP4 .. ', Error
      IF ( error .eq. 0) THEN
        JD = JDSatEpoch + T/1440.0D0
        CALL INVJDAY( JD, Year,Mon,Day,Hr,Min, Sec )
        IF (Year.ge.2000) THEN
          Yr = Year - 2000
        ELSE
          Yr = Year - 1900
        ENDIF
        MonStr = MonthTitle( Mon )
c-- Convert x, y, z to Long, Lat
        theta  = datan2(ro(2),ro(1))
        if (theta.lt.0) theta = theta + twopi
        theta0g  = gstime(JD)   ! + frac(JD * 86400.0D0 * 7.29211510D-5)
        thetaout = theta0g
        long     = (theta-theta0g)*rad
        lat      = datan(ro(3)/dsqrt(ro(1)*ro(1)+ro(2)*ro(2)))*rad
c-- WGS84 check
        xyz(1) = ro(1) * 1000.0D0
        xyz(2) = ro(2) * 1000.0D0
        xyz(3) = ro(3) * 1000.0D0
        call xyz2lla(xyz, lla)
        if (lla(1).lt.0) lla(1) = lla(1) + twopi
        
c--
        WRITE( 111,'(F17.8,3F17.8,3F17.8,1x,I4,1x,A3,I3,I3,A1,I2,A1,F9.6,6F17.4)' ) t,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3),
     &                   Day,MonStr,Yr,Hr,':',Min,':',Sec, long, lat, JD, (lla(1)-theta0g)*rad, lla(2)*rad, lla(3)/1000.0D0
        long = (lla(1) - theta0g) * rad
        if (long.lt.-180.0D0) long = long + 360.0D0
        if (long.gt. 180.0D0) long = long - 360.0D0
        lat  =  lla(2) * rad
      ENDIF ! if error
      return
c
      entry sgp4_close()
      close (unit=110)
      close (unit=111)
      close (unit=114)
      close (unit=115)
      return
      END

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

c-- Original in Octave
cfunction xyz = lla2xyz(lla)
c  a=6378137;            % radius
c  e=8.1819190842622e-2; % eccentricity
c  asq=a*a;
c  esq=e*e;
c  lon=lla(1);
c  lat=lla(2);
c  alt=lla(3);
c  N=a/sqrt(1-esq*sin(lat)^2);
c  x=(N+alt)*cos(lat)*cos(lon);
c  y=(N+alt)*cos(lat)*sin(lon);
c  z=((1-esq)*N+alt)*sin(lat);
c  xyz=[x y z];
cendfunction;

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

c-- Original in Octave
cfunction lla=xyz2lla(xyz)
c  a=6378137;            % radius
c  e=8.1819190842622e-2; % eccentricity
c  asq=a*a;
c  esq=e*e;
c  x=xyz(1);
c  y=xyz(2);
c  z=xyz(3);
c  b=sqrt(asq*(1-esq));
c  bsq=b*b;
c  ep=sqrt((asq-bsq)/bsq);
c  p=sqrt(x*x+y*y);
c  th=atan2(a*z,b*p);
c  lon=atan2(y,x);
c  lat=atan2((z+ep*ep*b*sin(th)^3),(p-esq*a*cos(th)^3));
c  N=a/(sqrt(1-esq*sin(lat)^2));
c  %alt=p/cos(lat)-N; % Not needed here
c  g=lla2xyz([lon lat 0]);
c  gm=sqrt(g(1)*g(1)+g(2)*g(2)+g(3)*g(3));
c  am=sqrt(x*x+y*y+z*z);
c  alt=am-gm;
c  lla=[lon lat alt];
cendfunction;
