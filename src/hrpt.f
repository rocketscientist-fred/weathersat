c
c-- Dr Fred A.Jansen - ESA/ESTEC - written for personal use and in private time, using no ESA resources
c-- Use at your own risk !!!
c--
c-- HRPT reading s/w for NOAA, Metop, Meteor and FengYun .raw16, .hpt, .dat and .C10 images (also .bin output directly from XHRPT can be processed with binary switch on command line)
c-- The software has been extended to be able to work with AQUA Modis data - using J-L Milette deframed .bin files.
c
c-- It is not at all difficult to make the s/w crash - start simple and build on from that
c-- I am convinced there are a plethora of bugs left, but it can work :-)
c
c-- Compile with NO -O2 for plotsoft.f - use buildhrpt.sh and export MMODEL=LARGE - for including AQUA 250 m - or without:  export MMODEL=MEDIUM : 
c--
c-- gfortran -c -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double plotsoft.f
c-- gfortran -o hrpt.exe -cpp -DMEDIUM -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 hrpt.f plotsoft.o
c--
c-- or :
c--
c-- gfortran -c -ffixed-line-length-350 -fno-range-check -mcmodel=large -m64 -malign-double plotsoft.f
c-- gfortran -o hrpt.exe -cpp -DLARGE -ffixed-line-length-350 -fno-range-check -mcmodel=large -m64 -malign-double -O2 hrpt.f plotsoft.o
c--
c-- Verified to work on WSL, Ubuntu through VirtualBox on Windows and a recent 64-bit cygwin with gfortran and ImageMagick - cygwin and OS-X can only run the medium case - so no MODIS 250 m
c
c-- Or simpler, build with ./buildhrpt.sh
c-- Using ENV parameter MMODEL:
c--
c--   export MMODEL=LARGE
c
c-- Known problems:
c--   A pass crossing midnight will cause havoc and for Meteor anything between 00:00 and 03:00 will FAIL !
c--   In the opacity mode channel 10 (as an A on input) is not implemented yet
c--   The FY3 pixel offset for the channels >= 6 is NOT used for the projected/border images :-( Don know how to really do this
c--   All the rest of the bugs I have forgotten or have escaped me to date :-(
c--
c-- Todo:
c-- histcor=0.5,99.5 ? make the %-tages configurable
c-- sharpen=6x11+0.4-3 - add a string to pass on verbatim to the ImageMagick sharpen operation
c-- Use .bins from Peter's GRC script ? Can't find the GRC formats
c-- I might just as well remove the FY truecolor attempts and leave it to the user with opacity= etc etc ??
c
c-- Program usage:
c
c-- See help text at the end !
c
c-- And now for the code credits:
c--
c-- Top section - until get_lun and free_lun - all written fully by me using wikipedia, NOAA and Eumetsat public documentation
c-- Plot code - some 50% written in the early 80's by Hans Deutekom at the Cosmic Ray Working Group in Leiden - heavily modified by me for postscript
c--             the rest all written by me (pltextbl, plraster, plhist, plpois etc etc) - now in plotsoft.f (as -O2 as a compiler switch ruined the output)
c-- SGP4 code from the 70's written by David Vallado (see T.S Kelso's celectrak.com) in ancient style fortran and where needed adapted by me.
c-- xyz2lla and lla2xyz taken from https://ea4eoz.blogspot.com/2015/11/simple-wgs-84-ecef-conversion-functions.html
c--             with an implemented improvement from the footnote
c
      implicit none
c
c-- Remember to change the parameters 4096 and 4095 as hardcoded values when you start with more than 12 bits per colour !
c
      integer*4 nmaxlin, nmaxpix, nmaxch, nstereo, n_sza, ndivgrid, nmergemax, ncolmax
#ifdef LARGE
c      parameter (nmaxlin=24000, nmaxpix=5416, nmaxch=10, nstereo=270000000, n_sza=1800, ndivgrid=4, nmergemax=5, ncolmax=4096)
      parameter (nmaxlin=24000, nmaxpix=5416, nmaxch=10, nstereo=300000000, n_sza=1800, ndivgrid=10, nmergemax=5, ncolmax=4096)
#endif
#ifdef MEDIUM
      parameter (nmaxlin=11000, nmaxpix=2708, nmaxch=10, nstereo=66000000, n_sza=1800, ndivgrid=10, nmergemax=5, ncolmax=4096)
#endif
c
      integer*1 linevalid(nmaxlin), satch_valid(nmaxlin), stereo_count(nstereo,3)
      integer*2 satch_image(nmaxpix,nmaxch,nmaxlin), stereo_image(3,nstereo)
      integer*4 satch_linedoy(nmaxlin), linedoy(nmaxlin), szahist(n_sza)
      real*8    satch_linetime(nmaxlin), linetime(nmaxlin), longlat(3,nmaxpix+1,nmaxlin), satch_az(nmaxlin), satch_el(nmaxlin)
c
      integer*2 satch_valr(ndivgrid, ndivgrid), satch_valg(ndivgrid, ndivgrid), satch_valb(ndivgrid, ndivgrid)
      integer*2 satch_szar(ndivgrid, ndivgrid), satch_szag(ndivgrid, ndivgrid), satch_szab(ndivgrid, ndivgrid)
      integer*2 satch_valr_op(ndivgrid, ndivgrid), satch_valg_op(ndivgrid, ndivgrid), satch_valb_op(ndivgrid, ndivgrid)
      integer*4 ndivuse
      real*8    satch_ll(2, ndivgrid, ndivgrid)
c
      real*8    theta_merge(nmergemax), rlinedelay_merge(nmergemax)
      integer*4 nmerge, ilinedelay_in_merge(nmergemax), i_merge, nlines_gapcor_merge(nmergemax), nsat_pix_merge(nmergemax), lundump(3,nmergemax)
      logical   swtheta_merge, swlinedelay_merge, swnorth_merge(nmergemax)
c
      integer*1 lut(3,256,256), lutbuf(768), inbuf(27728)
      integer*2 channel_image(nmaxpix,10), imgbuf(nmaxpix), rgbbuf(3,nmaxpix), rgbbuf_correct(3,10000), index_correct(10000), pixcol(3), npos, channel_image_bin_metop(2071,5)
      integer*2 mr, mg, mb, fetch_pixel
      integer*4 lun, irec, i, j, k, l, lunimg, ios, lunchan(10), hist(ncolmax,10), lunlut, lunrgb(10), k1, k2, histrgb(ncolmax), j1, j2, histchimg(ncolmax)
      integer*4 visch(4), nvisch, rgb(3,4), histch(4,ncolmax), jch1(4), jch2(4), luntrue, xchoff(10), nsat_pix, ncorrect_pix, doyfile, doy_off
      integer*4 lun_correct, luntrue_correct, nchan, ir, ig, ib, nrecl, doy, i_arg, lunborder, nxs, nys, i_sza, ilinedelay, lunbatch, iosbatch, lunmerge
      integer*4 nblank, nlines_gapcor, k3, nwrite, i_off, irp, igp, ibp, nvalid, szach, npix1, npix2, perc1, perc2, hist1(ncolmax), hist2(ncolmax)
      integer*4 irps, igps, ibps, ndiv, opchr, opchg, opchb, opperc, longstep, latstep, mv, i_bcol, i_llcol, iclim, naltavg, gammalut(ncolmax)
      real*4    wl(4), rsum, gsum, bsum, sat_avg_alt, sat_fov, lat_start, lat_end, rclim
      real*4    dr, dg, db, dy, di, dq, sum1, sum2, szasum, wr1, wr2, wop1, wop2, altavg, rgamma, ggamma, bgamma, rcorrect, gcorrect, bcorrect
      real*8    timestamp, timstart, timprev, ro(3), vo(3), long, lat, long_stereo, lat_stereo, ro_copy(3), lla(3), fetch_ll
      real*8    theta0g, thetascan, theta_in, szalim, w1, w2, w3, szafact, rlinedelay, theta_start, theta_end, theta_step, etascan, eta_in, aqmtheta(2), bowgamma, rpos
      logical   swnorth, swfy, swnoaa, swmetop, swgapcor, swswath, swdebug, swproject, swborder, swgamma, swtheta, swszalim, swszamer, swlinedelay, swhybrid, swfast, swthetascan, swyiq
      logical   swbatch, swbatchnight, swhorizon, swsharp, swtruech, swbin, swopacity, swmonitor, swhistcor, swlonglat, swbcol, swllcol, swmn2, sweta, swaqmod, sw12bit
      logical   swreflong, swreadbindone, swborderhighres, swcorrect, swok, swbowt, swmerge, swlastmerge, swfirstmerge, swgammargb, swbowgamma, swgammargb_in, swprojecteq, swterra
      character*4   cfy, cmetop, cscan, cmn2, caqmod
      character*6   cnoaa
      character*9   cgamma
      character*10  creflong
      character*20  ctheta, cszalim, cszafact, clinedelay, csharp, ceta, cbowgamma, ccorrect
      character*30  cbowtie, cgammargb
      character*250 filedata, argstring, imgname, command, filedatain, outstring, mergestring
c
c-- I hate common blocks, but for the MODIS processing it avoids numerous arguments to subroutines.
c
      logical       swaqm1000, swaqm500, swaqm250, swaqmpan
      common /MODIS/swaqm1000, swaqm500, swaqm250, swaqmpan
c
      call get_command(command)
      write (*,'(a)') command(1:lnblnk(command))
      if (index(command,'--help').ne.0) then
        call get_help()
        stop
      endif
c
      swtheta_merge     = .false.
      swlinedelay_merge = .false.
      do i = 1, nmergemax
        theta_merge(i)         = 0.0D0
        rlinedelay_merge(i)    = 0.0D0
        ilinedelay_in_merge(i) = 0
      enddo
      swmerge         = .false.
      swfirstmerge    = .false.
      swlastmerge     = .false.
c-- Set the interpolation factor to ensure all the projected grid points are 'hit' - especially important towards the edge of the scan
      ndiv            = 12
      swnorth         = .false.
c-- Arguments from 3 onwards have no fixed order - set the initial values
      swgapcor        = .false.
      swdebug         = .false.
      swproject       = .false.
      swprojecteq     = .false.
      irp             = 2
      igp             = 2
      ibp             = 1
      swborder        = .false.
      swgamma         = .false.
      swtheta         = .false.
      swthetascan     = .false.
      swszalim        = .false.
      szalim          = -1.0D0
      swszamer        = .false.
      szach           = 4
      szafact         = 1.2
      swlinedelay     = .false.
      ilinedelay      = 0
      swhybrid        = .false.
      swfast          = .false.
      swyiq           = .false.
      swbatch         = .false.
      swbatchnight    = .false.
      swhorizon       = .false.
      swsharp         = .false.
      swtruech        = .false.
      swbin           = .false.
      swopacity       = .false.
      swmonitor       = .false.
      swhistcor       = .false.
      swlonglat       = .false.
      swbcol          = .false.
      swllcol         = .false.
      sweta           = .false.
      sw12bit         = .false.
      creflong        = '          '
      swreflong       = .false.
      swaqm1000       = .false.
      swaqm500        = .false.
      swaqm250        = .false.
      swaqmpan        = .false.
      swborderhighres = .false.
      swcorrect       = .false.
      swbowt          = .false.
      swgammargb      = .false.
      swgammargb_in   = .false.
      rgamma          = 1.0
      ggamma          = 1.0
      bgamma          = 1.0
      rcorrect        = 0.0
      gcorrect        = 0.0
      bcorrect        = 0.0
      bowgamma        = 0.0
      swbowgamma      = .false.
      swterra         = .false.
c-- Check if modis !! In that case, 12 bit colours. Also check the modis resolution required.
      do i_arg = 3,30
        call getarg(i_arg, argstring)
        if (index (argstring,'modis').ne.0) then
          sw12bit = .true.
c-- By default use the 1 km data
          swaqm1000 = .true.
          swaqm500  = .false.
          swaqm250  = .false.
          swaqmpan  = .false.
          if (index (argstring,'=correct').ne.0) swcorrect = .true. 
          if (index (argstring,'modis500').ne.0) then
            swaqm500  = .true.
            swaqm1000 = .false.
            swaqm250  = .false.
            swaqmpan  = .false.
          endif
          if (index (argstring,'modis250').ne.0) then
            swaqm250  = .true.
            swaqm1000 = .false.
            swaqm500  = .false.
            swaqmpan  = .false.
          endif
          if (index (argstring,'modispan').ne.0) then
            swaqmpan  = .true.
            swaqm1000 = .false.
            swaqm500  = .false.
            swaqm250  = .false.
          endif
        endif
      enddo
      iclim = 1023
      rclim = 1023.0
      if (sw12bit) then
        iclim = 4095
        rclim = 4095.0
      endif
c-- 
      do i_arg = 3,30
        call getarg(i_arg, argstring)
        if (index (argstring,'gapcor').ne.0) swgapcor = .true.
        if (index (argstring,'debug').ne.0)  swdebug  = .true.
        if (index (argstring,'project').ne.0) then
          swproject = .true.
          if (index(argstring,'=').ne.0) then
            i = index(argstring,'=')
            if (i.ne.0) swprojecteq = .true.
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
        if (index (argstring,'border').ne.0) then
          swborder  = .true.
          if (index (argstring,'=highres').ne.0) swborderhighres = .true.
        endif
        if (index (argstring,'gamma=').ne.0.and.index(argstring,'bowgamma=').eq.0) then
          swgamma  = .true.
          cgamma   = argstring(1:9)
          i = index(cgamma,'=')
          cgamma(i:i) = ' '
        endif
        if (index (argstring,'theta=').ne.0) then
          swtheta  = .true.
          ctheta   = argstring(1:lnblnk(argstring))
          if (index(ctheta,',').eq.0) then
            i = index(ctheta,'=')
            call string_to_r8(ctheta(i+1:lnblnk(ctheta))//' ',npos, theta_in)
            write (*,'(''Theta override     : '',F10.4)') theta_in
          else
            i = index(ctheta,'=')
            j = index(ctheta,',')
            call string_to_r8(ctheta(i+1:j-1)//' ',npos, theta_start)
            i = j
            if (index(ctheta(j+1:),',').eq.0) stop '** Improper theta step format **'
            j = index(ctheta(j+1:),',') + j
            call string_to_r8(ctheta(i+1:j-1)//' ',npos, theta_end)
            i = j
            call string_to_r8(ctheta(i+1:lnblnk(ctheta))//' ',npos, theta_step)
            write (*,'(''Theta step from    : '',F10.4,'' to : '',F10.4,'' with step size : '',F10.4)') theta_start, theta_end, theta_step
            swthetascan = .true.
          endif
        endif
        if (index (argstring,'szalim=').ne.0) then
c-- Requires swproject to be true
          swproject = .true.
          swszalim  = .true.
          cszalim  = argstring(1:lnblnk(argstring))
          i = index(cszalim,'=')
          call string_to_r8(cszalim(i+1:lnblnk(cszalim))//' ',npos, szalim)
          write (*,'(''SZA limit          : '',F10.4)') szalim
        endif
        if (index(argstring,'szamerge=').ne.0) then
          swszamer = .true.
          if (.not.swszalim) stop '** Define szalim= before szamerge= **'
          i = index(argstring,'=')
          if (ichar('0').le.ichar(argstring(i+1:i+1)).and.ichar(argstring(i+1:i+1)).le.ichar('9')) then
            read (argstring(i+1:i+1),'(i1)') szach
          else if (argstring(i+1:i+1).eq.'a'.or.argstring(i+1:i+1).eq.'A') then
            szach = 10
          endif
        endif
        if (index (argstring,'szafact=').ne.0) then
          cszafact = argstring(1:lnblnk(argstring))
          i = index(cszafact,'=')
          call string_to_r8(cszafact(i+1:lnblnk(cszafact))//' ',npos, szafact)
          write (*,'(''SZAfact override   : '',F10.4)') szafact
        endif
        if (index (argstring,'linedelay=').ne.0) then
          swlinedelay = .true.
          clinedelay = argstring(1:lnblnk(argstring))
          i = index(clinedelay,'=')
          call string_to_r8(clinedelay(i+1:lnblnk(clinedelay))//' ',npos, rlinedelay)
          ilinedelay = int(rlinedelay)
          write (*,'(''Line (=TLE) delay  : '',I10)') ilinedelay
        endif
        if (index (argstring,'hybrid').ne.0) swhybrid   = .true.
        if (index (argstring,'fast').ne.0) swfast       = .true.
        if (index (argstring,'YIQ').ne.0) swyiq         = .true.
        if (index (argstring,'horizon').ne.0) swhorizon = .true.
        if (index (argstring,'sharpen').ne.0) then
          swsharp = .true.
c                    12345678901234567890
          csharp  = '6x3+0.5+0           '
        endif
        if (index (argstring,'truech').ne.0) swtruech = .true.
        if (index (argstring,'binary').ne.0) swbin    = .true.
        if (index (argstring,'opacity').ne.0) then
          swopacity = .true.
          opchr  = 1
          opchg  = 9
          opchb  = 7
          opperc = 25
          i      = index(argstring,'=')
          if (i.ne.0) then
            read (argstring(i+1:i+1),'(I1)') opchr
            read (argstring(i+2:i+2),'(I1)') opchg
            read (argstring(i+3:i+3),'(I1)') opchb
            write (*,'(''Opacity RGB        : '',4x,3I2)') opchr, opchg, opchb
          endif
          i      = index(argstring,',')
          if (i.ne.0) then
            call string_to_i4(argstring(i+1:lnblnk(argstring))//' ',npos, opperc)
            write (*,'(''Opacity %          : '',I10)') opperc
          endif
        endif
        if (index(argstring,'monitor').ne.0) swmonitor = .true.
        if (index(argstring,'histcor').ne.0) swhistcor = .true.
        if (index (argstring,'longlat').ne.0) then
          swlonglat = .true.
          longstep  = 10
          latstep   = 10
          i         = index(argstring,'=')
          if (i.ne.0) then
            j      = index(argstring,',')
            call string_to_i4(argstring(i+1:j-1)//' ',npos, longstep)
            call string_to_i4(argstring(j+1:lnblnk(argstring))//' ',npos, latstep)
          endif
        endif
        if (index (argstring,'bcol=').ne.0) then
          swbcol = .true.
          i_bcol = i_arg
          i  = index(argstring,'=')
          j  = index(argstring,',')
          call string_to_i4(argstring(i+1:j-1)//' ',npos, mv)
          mr = min(max(mv,0),iclim)
          i  = j
          if (index(argstring(j+1:),',').eq.0) stop '** Improper bcol  RGB format **'
          j  = index(argstring(j+1:),',') + j
          call string_to_i4(argstring(i+1:j-1)//' ',npos, mv)
          mg = min(max(mv,0),iclim)
          i  = j
          call string_to_i4(argstring(i+1:lnblnk(argstring))//' ',npos, mv)
          mb = min(max(mv,0),iclim)
          write (*,'(''Border  RGB      R : '',I10,''  G : '',I10,''              B : '',I10)') mr, mg, mb
        endif
        if (index (argstring,'llcol=').ne.0) then
          swllcol = .true.
          i_llcol = i_arg
          i  = index(argstring,'=')
          j  = index(argstring,',')
          call string_to_i4(argstring(i+1:j-1)//' ',npos, mv)
          mr = min(max(mv,0),iclim)
          i  = j
          if (index(argstring(j+1:),',').eq.0) stop '** Improper llcol RGB format **'
          j  = index(argstring(j+1:),',') + j
          call string_to_i4(argstring(i+1:j-1)//' ',npos, mv)
          mg = min(max(mv,0),iclim)
          i  = j
          call string_to_i4(argstring(i+1:lnblnk(argstring))//' ',npos, mv)
          mb = min(max(mv,0),iclim)
          write (*,'(''Longlat RGB      R : '',I10,''  G : '',I10,''              B : '',I10)') mr, mg, mb
        endif
        if (index(argstring,'theta=').eq.0.and.index (argstring,'eta=').ne.0) then
          sweta  = .true.
          ceta   = argstring(1:lnblnk(argstring))
          i = index(ceta,'=')
          call string_to_r8(ceta(i+1:lnblnk(ceta))//' ',npos, eta_in)
          write (*,'(''Eta override       : '',F10.4)') eta_in
        endif
        if (index(argstring,'reflong=').ne.0) then
          swreflong = .true.
          i = index(argstring,'=')
          creflong = argstring(i+1:lnblnk(argstring))
        endif
        if (index(argstring,'bowtie=').ne.0) then
          swbowt = .true.
          do i = 1, len(cbowtie)
            cbowtie(i:i) = ' '
          enddo
          cbowtie = argstring(1:lnblnk(argstring))
        endif
c-- Check if merge and if so load files/parameters
        if (index(argstring,'merge').ne.0.and.index(argstring,'sza').eq.0) then
          swmerge      = .true.
          swfirstmerge = .true.
          swlastmerge  = .false.
          i_merge      = 1
c-- Read the merge file en store the parameters - also flush the filenames into the format/name used by the batch option, as merge uses this logic.
          call getarg(1, mergestring)
          call get_lun(lunmerge)
          open (unit=lunmerge, file=mergestring(1:lnblnk(mergestring)), form='formatted')
          call get_lun(lunbatch)
          open (unit=lunbatch,file='batch.txt',form='formatted')
          nmerge = 0
          ios    = 0
          do while (ios.eq.0) 
            read (lunmerge,'(a)',iostat=ios) mergestring
            if (ios.ne.0) goto 102
            i = 1
            do while (mergestring(i:i).eq.' '.and.i.le.len(mergestring))
              i = i + 1
            enddo
            k1 = i
            do while (mergestring(i:i).ne.' '.and.i.le.len(mergestring))
              i = i + 1
            enddo
            k2 = i - 1
            if (nmerge.ge.nmergemax) goto 102
            nmerge = nmerge + 1
            write (lunbatch,'(a)') mergestring(k1:k2)
c
            if (index(mergestring,'theta=').ne.0) then
              k1 = index(mergestring,'theta=') + 6
              i  = k1
              do while (mergestring(i:i).ne.' '.and.i.le.len(mergestring))
                i = i + 1
              enddo
              k2 = i - 1
              swtheta_merge  = .true.
              call string_to_r8(mergestring(k1:k2)//' ',npos, theta_merge(nmerge))
            endif
c
            if (index(mergestring,'linedelay=').ne.0) then
              k1 = index(mergestring,'linedelay=') + 10
              i  = k1
              do while (mergestring(i:i).ne.' '.and.i.le.len(mergestring))
                i = i + 1
              enddo
              k2 = i - 1
              swlinedelay_merge = .true.
              call string_to_r8(mergestring(k1:k2)//' ',npos, rlinedelay_merge(nmerge))
              ilinedelay_in_merge(nmerge) = int(rlinedelay_merge(nmerge))
            endif
          enddo
 102      close (unit=lunmerge)
          close (unit=lunbatch)
          call free_lun(lunmerge)
          call free_lun(lunbatch)
c-- End of reading merge file
        endif
        if (index (argstring,'gammargb=').ne.0) then
          swgammargb    = .true.
          swgammargb_in = .true.
          cgammargb     = argstring(1:lnblnk(argstring))
          i = index(cgammargb,'=')
          j = index(cgammargb,',')
          if (j.eq.0) stop '** Improper gammargb format **'
          call string_to_r4(cgammargb(i+1:j-1)//' ',npos, rgamma)
          i = j
          if (index(cgammargb(j+1:),',').eq.0) stop '** Improper gammargb format **'
          j = index(cgammargb(j+1:),',') + j
          call string_to_r4(cgammargb(i+1:j-1)//' ',npos, ggamma)
          i = j
          call string_to_r4(cgammargb(i+1:lnblnk(cgammargb))//' ',npos, bgamma)
        endif
        if (index (argstring,'bowgamma=').ne.0) then
          swbowgamma = .true.
          cbowgamma  = argstring(1:lnblnk(argstring))
          i = index(cbowgamma,'=')
          call string_to_r8(cbowgamma(i+1:lnblnk(cbowgamma))//' ',npos, bowgamma)
          write (*,'(''BowTie Gamma       : '',F10.4)') bowgamma
        endif
        if (index (argstring,'rcorrect=').ne.0) then
          ccorrect  = argstring(1:lnblnk(argstring))
          i = index(ccorrect,'=')
          call string_to_r4(ccorrect(i+1:lnblnk(ccorrect))//' ',npos, rcorrect)
        endif
        if (index (argstring,'gcorrect=').ne.0) then
          ccorrect  = argstring(1:lnblnk(argstring))
          i = index(ccorrect,'=')
          call string_to_r4(ccorrect(i+1:lnblnk(ccorrect))//' ',npos, gcorrect)
        endif
        if (index (argstring,'bcorrect=').ne.0) then
          ccorrect  = argstring(1:lnblnk(argstring))
          i = index(ccorrect,'=')
          call string_to_r4(ccorrect(i+1:lnblnk(ccorrect))//' ',npos, bcorrect)
        endif
      enddo
      if ((swborder.or.swlonglat).and.swsharp) write (*,'(''W A R N I N G : border/longlat and sharpen together may give different colour values for these !!'')')
c-- Borders require project
      if (swborder) swproject = .true.
c-- Polar stereograhic projection requires gapcor
      if (swproject) swgapcor = .true.
c-- If needed, rotate the polar stereograhic map, otherwise put back the polster_orig
      if (swreflong) then
        command = './polarstereo.exe '//creflong(1:lnblnk(creflong))
        write (*,'(a)') command(1:lnblnk(command))
        call system(command)
      else
        command = 'cp ./resource/polster_orig.rgb polster.rgb'
        write (*,'(a)') command(1:lnblnk(command))
        call system(command)
      endif
c-- Set the mission dependent variables
      swfy    = .false.
      swnoaa  = .false.
      swmetop = .false.
      swmn2   = .false.
      swaqmod = .false.
c
      call getarg(1, argstring)
      if (argstring(1:1).eq.'@') then
        swbatch = .true.
        if (index(argstring,'=').ne.0) then
          if (index(argstring,'=.C10').ne.0)   call system('ls -a1 '//argstring(2:index(argstring,'=')-1)//'*.C10   > batch.txt')
          if (index(argstring,'=.hpt').ne.0)   call system('ls -a1 '//argstring(2:index(argstring,'=')-1)//'*.hpt   > batch.txt')
          if (index(argstring,'=.raw16').ne.0) call system('ls -a1 '//argstring(2:index(argstring,'=')-1)//'*.raw16 > batch.txt')
          if (index(argstring,'=.bin').ne.0)   call system('ls -a1 '//argstring(2:index(argstring,'=')-1)//'*.bin   > batch.txt')
        else
          if (swbin) then
            if (swaqm250.or.swaqm500.or.swaqm1000.or.swaqmpan) then
              call system('ls -a1t '//argstring(2:lnblnk(argstring))//'AQ*.bin 2> /dev/null > batch.txt')
            else
              call system('ls -a1t '//argstring(2:lnblnk(argstring))//'*.bin '//argstring(2:lnblnk(argstring))//'*.raw16 2> /dev/null > batch.txt')
            endif
          else
            call system('ls -a1 '//argstring(2:lnblnk(argstring))//'*.C10 2> /dev/null > batch.txt')
            call system('ls -a1 '//argstring(2:lnblnk(argstring))//'*.hpt 2> /dev/null >> batch.txt')
            call system('ls -a1 '//argstring(2:lnblnk(argstring))//'*.raw16 2> /dev/null >> batch.txt')
          endif
        endif
        call get_lun(lunbatch)
        open (unit=lunbatch,file='batch.txt',form='formatted')
        iosbatch = 0
c-- Enable SZA calculation for automated night side switching !!
        swszalim = .true.
      endif
c-- Remember, merging files is implemented using the batch functionality, so swbatch is switched on here (not before as we do not need to call the ls command as above)
c-- The automated night time switching as is usual in batch mode, makes no sense for merging and is disabled !!!! So don't merge files from daytime and nighttime !!
      if (swmerge) then
        swbatch = .true.
        call get_lun(lunbatch)
        open (unit=lunbatch,file='batch.txt',form='formatted')
        iosbatch = 0
      endif
      if (swbatch.and.swhorizon.and.swgapcor) call map_azel_horizon_init()
 200  if (swbatch) then
        if (swmonitor.and.swbin) then
          call getarg(1, argstring)
          call monitor(argstring(2:lnblnk(argstring)), filedatain)
          write (*,'(''Processing monitor : '',a)') filedatain(1:lnblnk(filedatain))
        else
          read (lunbatch,'(a)',iostat=iosbatch) filedatain
          if (iosbatch.ne.0) goto 201
          write (*,'(''Processing queue   : '',a)') filedatain(1:lnblnk(filedatain))
        endif
        call save_string(filedatain(1:lnblnk(filedatain)))
        if (swbatchnight) then
          swbatchnight = .false.
          irp = irps
          igp = igps
          ibp = ibps
        endif
      else
        filedatain = argstring(1:lnblnk(argstring))
      endif
      j = 0
      i = lnblnk(filedatain)
      do while (i.ge.1.and.filedatain(i:i).ne.'/')
        i = i - 1
      enddo
      j = i + 1
      do i = 1, len(outstring)
        outstring(i:i) = ' '
      enddo
      call getenv('HRPTOUT',outstring)
      filedata = outstring(1:lnblnk(outstring))//filedatain(j:lnblnk(filedatain))
c-- Implement the modis PAN here. For the time being is is hardcoded to use ch 143 for RGB, with 1 being 250 m and 4 and 3 500 m - use i_merge to select
      if (swmerge.and.swaqmpan) then
        if (i_merge.eq.3) then
          swaqm1000 = .false.
          swaqm500  = .true.
          swaqm250  = .false.
        else if (i_merge.eq.2) then
          swaqm1000 = .false.
          swaqm500  = .true.
          swaqm250  = .false.
        else if (i_merge.eq.1) then
          swaqm1000 = .false.
          swaqm500  = .false.
          swaqm250  = .true.
        endif
      endif
c-- This is confusing a bit - some sats need gammargb for their default images (FY, AQUA), but they should not mess up the other sats in e.g. a batch processing
c-- So, here I determine if there was a command line defined gammargb, if not I reset swgammargb t.false.
      if (.not.swgammargb_in) swgammargb = .false.
c
      if (index(filedata,'FY').ne.0) then
        nchan   = 10
        swfy    = .true.
        swmetop = .false.
        swnoaa  = .false.
        swmn2   = .false.
        swaqmod = .false.
        cfy     = '.C10'
        if (swbin) cfy = '.frd'
c-- A good colour setting for FY appears to be RGB=197 with gammargb=0.93,0.95,0.90 - but, work out how to reset such values for batch processing
        ir      = 6
        ig      = 2
        ib      = 1
        if (.not.swprojecteq) irp = 1
        if (.not.swprojecteq) igp = 9
        if (.not.swprojecteq) ibp = 7
        if (.not.swgammargb_in) then
          swgammargb = .true.
          rgamma = 0.93
          bgamma = 0.93
          ggamma = 0.92
        endif
        nrecl = 27728
        if (swbin) nrecl = 26050
c-- Initialise the projection correction arrays
c-- Remember that the IFOV (= pixel size) is BIGGER than the sampling interval !
        sat_avg_alt  = 833.0
        sat_fov      = 55.37 * 2.0
        nsat_pix     = 2048
        ncorrect_pix = 0
        thetascan    = 0.00D0
        etascan      = 0.0D0
        if (index(filedata,'FY3C').ne.0.and.(.not.swlinedelay)) then
          swlinedelay = .true.
          ilinedelay  = -5
        endif
        if (index(filedata,'FY3B').ne.0.and.(.not.swlinedelay)) then
          swlinedelay = .true.
          ilinedelay  = -4
        endif
c-- Define the FY3B and FY3C channels to be used for trucolor
        nvisch   = 4
        visch(1) = 1
        visch(2) = 7
        visch(3) = 9
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
      if (index(filedata,'NOAA').ne.0) then
        nchan   = 5
        swnoaa  = .true.
        swfy    = .false.
        swmetop = .false.
        swmn2   = .false.
        swaqmod = .false.
        cnoaa   = '.raw16'
        ir      = 2
        ig      = 2
        ib      = 1
        nrecl   = 22180
c-- Initialise the projection correction arrays
c-- Remember that the IFOV (= pixel size) is BIGGER than the sampling interval !
        sat_avg_alt  = 850.0
        sat_fov      = 55.37 * 2.0
        nsat_pix     = 2048
        ncorrect_pix = 0
        nvisch       = 1
        if (index(filedata,'NOAA-18').ne.0.or.index(filedata,'NOAA_18').ne.0) then
          thetascan  = -1.40D0
        else
          thetascan  =  0.00D0
        endif
        etascan      = 0.0D0
        ilinedelay   = 0
c-- Define the NOAA x offsets ! 
        xchoff(1) =  0
        xchoff(2) =  0
        xchoff(3) =  0
        xchoff(4) =  0
        xchoff(5) =  0
      endif
      if (index(filedata,'Metop').ne.0.or.index(filedata,'M02').ne.0) then
        nchan   = 5
        swmetop = .true.
        swfy    = .false.
        swnoaa  = .false.
        swmn2   = .false.
        swaqmod = .false.
        cmetop  = '.hpt'
        if (swbin) cmetop = '.frd'
        ir      = 2
        ig      = 2
        ib      = 1
        nrecl   = 13864
        if (swbin) nrecl = 12960
c-- Initialise the projection correction arrays
c-- Remember that the IFOV (= pixel size) is BIGGER than the sampling interval !
        sat_avg_alt  = 805.0
        sat_fov      = 55.37 * 2.0
        nsat_pix     = 2048
        ncorrect_pix = 0
        nvisch       = 1
        thetascan    = -2.90D0
        etascan      = 0.0D0
        ilinedelay   = 0
c-- Define the Metop x offsets ! 
        xchoff(1) =  0
        xchoff(2) =  0
        xchoff(3) =  0
        xchoff(4) =  0
        xchoff(5) =  0
      endif
      if (index(filedata,'MN2').ne.0) then
        nchan   = 6
        swmn2   = .true.
        swfy    = .false.
        swnoaa  = .false.
        swmetop = .false.
        swaqmod = .false.
        cmn2    = '.dat'
        if (swbin) cmn2 = '.frd'
        if (.not.swbin) stop '** MN2 currently only supports .bin files **'
        ir      = 2
        ig      = 2
        ib      = 1
c-- Next line to be corrected if/when I support MN2 .dat files
        nrecl   = 0
        if (swbin) nrecl = 11850
c-- Initialise the projection correction arrays
c-- Remember that the IFOV (= pixel size) is BIGGER than the sampling interval !
        sat_avg_alt  = 833.0
        ilinedelay   = 0
        if (index(filedata,'MN2-2').ne.0) then
          sat_fov     = 55.25 * 2.0
          if (.not.swlinedelay) then
            swlinedelay = .true.
            ilinedelay  = 3
          endif
          etascan     = 0.0D0
        else
          sat_fov     = 55.25 * 2.0
          etascan     = 2.4D0
        endif
        nsat_pix     = 1572
        ncorrect_pix = 0
        nvisch       = 1
        thetascan    = 0.0D0
c-- Define the MN2 x offsets ! 
        xchoff(1) =  0
        xchoff(2) =  0
        xchoff(3) =  0
        xchoff(4) =  0
        xchoff(5) =  0
        xchoff(6) =  0
      endif
c-- Try AQUA MODIS ....... things like nr. of channels, RGB selection etc etc require loooooooooooooooots of work more
      if (index(filedata,'AQ_MODIS').ne.0) then
        nchan   = 5
        if (swaqm250) nchan = 2
        swmetop = .false.
        swfy    = .false.
        swnoaa  = .false.
        swmn2   = .false.
        swaqmod = .true.
        caqmod  = '.dat'
        ir      = 2
        ig      = 2
        ib      = 1
        nrecl   = 2708
c-- Initialise the projection correction arrays
c-- Remember that the IFOV (= pixel size) is BIGGER than the sampling interval !
        sat_avg_alt  = 705.0
        sat_fov      = 55.00 * 2.0
        if (swaqm1000) nsat_pix = 1354
        if (swaqm500)  nsat_pix = 2708
        if (swaqm250)  nsat_pix = 5416
        nrecl        = nsat_pix * 2
        ncorrect_pix =  0
        nvisch       =  1
c-- Experimentally determined in AQUA MODIS pass of Sep 29, 2019 - data from Jean-Luc
c-- Remember (!!!) thetascan is defined for a southern pass - the sign reversal for a Northern pass is done further down in the code !
        if (.not.swlinedelay) then
          ilinedelay   = -5
          swlinedelay  = .true.
        endif
        thetascan    =  0.100D0
        etascan      =  0.0D0
c-- Switch on the bowgamma => the best solution over a linear 0.18 - 0.42 bowtie slope
        if (.not.swbowgamma) then
          swbowgamma = .true.
          bowgamma   = 1.6
        endif
c-- For AQUA Modis modispan switch on the gammargb by default !!
        if ((.not.swgammargb_in).and.swaqmpan) then
          swgammargb = .true.
          rgamma = 0.95
          ggamma = 0.99
          bgamma = 0.93
        endif
        if ((.not.swprojecteq).and.swaqmpan) then
          irp = 1
          igp = 2
          ibp = 1
        endif
c-- Define the AQUA MODIS x offsets ! 
        xchoff(1) =  0
        xchoff(2) =  0
        xchoff(3) =  0
        xchoff(4) =  0
        xchoff(5) =  0
      endif
      if (index(filedata,'TE_MODIS').ne.0) then
        nchan   = 5
        if (swaqm250) nchan = 2
        swmetop = .false.
        swfy    = .false.
        swnoaa  = .false.
        swmn2   = .false.
        swaqmod = .true.
        swterra = .true.
        caqmod  = '.dat'
        ir      = 2
        ig      = 2
        ib      = 1
        nrecl   = 2708
c-- Initialise the projection correction arrays
c-- Remember that the IFOV (= pixel size) is BIGGER than the sampling interval !
        sat_avg_alt  = 705.0
        sat_fov      = 55.00 * 2.0
        if (swaqm1000) nsat_pix = 1354
        if (swaqm500)  nsat_pix = 2708
        if (swaqm250)  nsat_pix = 5416
        nrecl        = nsat_pix * 2
        ncorrect_pix =  0
        nvisch       =  1
c-- Copy from AQUA above
c-- Remember (!!!) thetascan is defined for a southern pass - the sign reversal for a Northern pass is done further down in the code !
        if (.not.swlinedelay) then
          ilinedelay   = -15
          swlinedelay  = .true.
        endif
        thetascan    =  0.100D0
        etascan      =  0.0D0
c-- Switch on the bowgamma => the best solution over a linear 0.18 - 0.42 bowtie slope
        if (.not.swbowgamma) then
          swbowgamma = .true.
          bowgamma   = 1.6
        endif
c-- For TERRA Modis modispan switch on the gammargb by default !!
        if ((.not.swgammargb_in).and.swaqmpan) then
          swgammargb = .true.
          rgamma = 0.95
          ggamma = 0.99
          bgamma = 0.93
        endif
        if ((.not.swprojecteq).and.swaqmpan) then
          irp = 1
          igp = 2
          ibp = 1
        endif
c-- Define the TERRA MODIS x offsets ! 
        xchoff(1) =  0
        xchoff(2) =  0
        xchoff(3) =  0
        xchoff(4) =  0
        xchoff(5) =  0
      endif
c-- In case of merge get the linedelay and/or theta (if defined in the merge file - remember, if it's there on 1 line it has to be there for all !!) and override the defaults here
      if (swmerge.and.swlinedelay_merge) then
        ilinedelay  = ilinedelay_in_merge(i_merge)
        swlinedelay = .true.
      endif
      if (swmerge.and.swtheta_merge) then
        theta_in  = theta_merge(i_merge)
        swtheta   = .true.
      endif
c
      call getarg(2,argstring)
      if (ichar('0').le.ichar(argstring(1:1)).and.ichar(argstring(1:1)).le.ichar('9').and.
     *    ichar('0').le.ichar(argstring(2:2)).and.ichar(argstring(2:2)).le.ichar('9').and.
     *    ichar('0').le.ichar(argstring(3:3)).and.ichar(argstring(3:3)).le.ichar('9')) then
        read (argstring(1:1),'(i1)') ir
        read (argstring(2:2),'(i1)') ig
        read (argstring(3:3),'(i1)') ib
      endif
      call correct_init(sat_avg_alt, sat_fov, nsat_pix, index_correct, 10000, ncorrect_pix,swaqmod)
c-- Calculate the 10 or 12 bit RGB for each of the visual channels - only used for FY3
      wl(1)    = (580.0 + 680.0) / 2.
      wl(2)    = (430.0 + 480.0) / 2.
      wl(3)    = (480.0 + 530.0) / 2.
      wl(4)    = (530.0 + 580.0) / 2.
      do i = 1, nvisch
        call wltoRGB(wl(i), rgb(1,i), rclim)
      enddo
c-- Calculate the sum of the R, G and B values to later use for normalisation
      rsum = 0.
      gsum = 0.
      bsum = 0.
      do i = 1, nvisch
        rsum = rsum + float(rgb(1,i))
        gsum = gsum + float(rgb(2,i))
        bsum = bsum + float(rgb(3,i))
      enddo
      rsum = rsum / rclim
      gsum = gsum / rclim
      bsum = bsum / rclim
      do i = 1, nvisch
        if (swfy) write (*,'(I5, 2x, F6.1, 2x, 3(Z4.4,3x))') i, wl(i), rgb(1,i), rgb(2,i), rgb(3,i)
      enddo
      if (swfy) write (*,'(15x, 3(F5.3,2x))') rsum, gsum, bsum
c      write (*,'(3(F6.2,2x))') rsum, gsum, bsum
c-- Read the RGB LUT - it is a 2D LUT in only 8 bits unfortunately
      if (swhybrid) then
        call get_lun(lunlut)
        open (unit=lunlut,file='./resource/LUT.rgb',form='unformatted',access='direct',recl=768)
        do i = 1,256
          read (lunlut,rec=i) lutbuf
          l = 1
          do j = 1,256
            do k = 1,3
              lut(k,i,j) = lutbuf(l)
              l = l + 1
            enddo
          enddo
        enddo
        close (unit=lunlut)
        call free_lun(lunlut)
      endif
c-- Start theta step loop here 
      if (swthetascan) thetascan = theta_start - theta_step
c-- This flag ensures the .bin to .frd conversion is only done once :-)
      swreadbindone = .false.
 100  continue
c-- Initialise all essential arrays to zero
      do i = 1, 4096
        histrgb(i) = 0
      enddo
      do i = 1, 4096
        do j = 1,4
          histch(j,i) = 0
        enddo
      enddo
      do i = 1,4096
        do j = 1,10
          hist(i,j) = 0
        enddo
      enddo
c
      do i = 1, nmaxlin
        satch_valid(i) = 0
        do j = 1,nmaxch
          do k = 1,nmaxpix
            satch_image(k, j, i) = 0
          enddo
        enddo
        do j = 1, nmaxpix + 1
          longlat(1,j,i) = 0.0D0
          longlat(2,j,i) = 0.0D0
          longlat(3,j,i) = 0.0D0
        enddo
      enddo
      do i = 1, n_sza
        szahist(i) = 0
      enddo
c
c-- For the case of processing binary .bin files replace the .bin with .frd in the filename as this is the filetype output by the call to readbin !
c-- further modifications made to support MODIS
      if (swbin.and.(.not.swnoaa).and.(.not.swreadbindone)) then
c-- This flag ensures the .bin to .frd conversion is only done once :-) - does this work with merge ??????? To be checked.
        swreadbindone = .true.
        i = index(filedata,'.bin')
        if (i.eq.0) stop 'No filename with .bin as an extension was provided'
        if (.not.swaqmod) then
          filedata(i:i+3) = '.frd'
          write (command,'(''./readbin.exe '',a,1x,a)') filedatain(1:lnblnk(filedatain)), filedata(1:lnblnk(filedata)) 
          write (*,'(a)') command(1:lnblnk(command))
          call system(command(1:lnblnk(command)))
          filedatain = filedata
        else
          filedata(i:i+10) = '-bandxx.dat'
          if (swcorrect) then
            if (swbowt) then
              write (command,'(''./readbin_modis.exe '',a,'' nopng correct '',a)') filedatain(1:lnblnk(filedatain)), cbowtie(1:lnblnk(cbowtie))
            else
              write (command,'(''./readbin_modis.exe '',a,'' nopng correct'')') filedatain(1:lnblnk(filedatain)) 
            endif
          else
            if (swbowt) then
              write (command,'(''./readbin_modis.exe '',a,'' nopng '',a)') filedatain(1:lnblnk(filedatain)), cbowtie(1:lnblnk(cbowtie))
            else
              write (command,'(''./readbin_modis.exe '',a,'' nopng'')') filedatain(1:lnblnk(filedatain)) 
            endif
          endif
          write (*,'(a)') command(1:lnblnk(command))
          call system(command(1:lnblnk(command)))
          filedatain = filedata
          call save_string(filedatain(1:lnblnk(filedatain)))
        endif
      endif
      call get_lun(lun)
      open (unit=lun,file=filedatain(1:lnblnk(filedatain)),form='unformatted',access='direct',recl=nrecl)
      if (swfy)    j = index(filedata,cfy)
      if (swnoaa)  j = index(filedata,cnoaa)
      if (swmetop) j = index(filedata,cmetop)
      if (swmn2)   j = index(filedata,cmn2)
      if (swaqmod) j = index(filedata,caqmod)
      if (.not.swfast) then
        do i = 1,nchan
          call get_lun(lunchan(i))
          write (imgname,'(a,a,z1,a)') filedata(1:j-1),'_ch',i,'.dat'
          open (unit=lunchan(i),file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=nsat_pix*2)
          close (unit=lunchan(i),status='delete', iostat=ios)
          open (unit=lunchan(i),file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=nsat_pix*2)
        enddo
      endif
      do i = 1,nchan
        if (i.eq.1.or.swhybrid) then
          call get_lun(lunrgb(i))
          write (imgname,'(a,a,z1,a)') filedata(1:j-1),'_ch1',i,'.rgb'
          open (unit=lunrgb(i),file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=nsat_pix*6)
          close (unit=lunrgb(i),status='delete', iostat=ios)
          open (unit=lunrgb(i),file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=nsat_pix*6)
        endif
      enddo
      if (swfy) then
        call get_lun(luntrue)
        write (imgname,'(a,a)') filedata(1:j-1),'_true.rgb'
        open (unit=luntrue,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=nsat_pix*6)
        close (unit=luntrue,status='delete', iostat=ios)
        open (unit=luntrue,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=nsat_pix*6)
c
        call get_lun(luntrue_correct)
        write (imgname,'(a,a)') filedata(1:j-1),'_true_correct.rgb'
        open (unit=luntrue_correct,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=2*ncorrect_pix*6)
        close (unit=luntrue_correct,status='delete', iostat=ios)
        open (unit=luntrue_correct,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=2*ncorrect_pix*6)
      endif
      call get_lun(lun_correct)
      write (imgname,'(a,a)') filedata(1:j-1),'_ch11_correct.rgb'
      open (unit=lun_correct,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=2*ncorrect_pix*6)
      close (unit=lun_correct,status='delete', iostat=ios)
      open (unit=lun_correct,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=2*ncorrect_pix*6)
c
c-- All the prerequisites have been fullfilled - start reading the data !!
c
c-- First loop through the image to retrieve the timestamp per scanline (needed to reject lines later)
c
      irec = 1
      do while (.true.)
        ios = 0
        if (swfy)    call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'FY', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
        if (swnoaa)  call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'NO', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
        if (swmetop) then
          if (swbin) then
            call get_scanline(lun, irec, ios, channel_image_bin_metop, 2071, nchan, nsat_pix, 'ME', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
            do i = 1, nchan
              do j = 1, nsat_pix
                channel_image(j,i) = channel_image_bin_metop(j+11,i)
              enddo
            enddo
          else
            call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'ME', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
          endif
        endif
        if (swmn2)   call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'M2', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
        if (swaqmod) call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'AM', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
        if (ios.ne.0) goto 2
        linetime(irec) = timestamp
        linedoy(irec)  = doy
        irec = irec + 1
      enddo
 2    irec = irec - 1
      call check_linetimes(linetime, linedoy, linevalid, irec, timstart, swdebug, swmn2, swaqmod)
c-- Allow for the TLE being off in time just a bit. Specified in # of lines.
      if (swlinedelay) then
        write (*,'(a,I8)') 'Line delay applied : ',ilinedelay
        do i = 1, irec
          if (linevalid(i).eq.1) then
c-- For MODIS the linedelay is defined in 250 m pixels - this is required for modispan to align the 500 m and 250 m data in truecolour images
            if (swaqmod) then
              linetime(i) = linetime(i) + dble(ilinedelay) * (1.477815D0 / 40.0D0)
            else
              linetime(i) = linetime(i) + dble(ilinedelay) * 0.1667D0
            endif
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
      if (swfy)                       call init_TLE('FY', filedata(1:lnblnk(filedata)), linedoy(i), doyfile)
      if (swnoaa)                     call init_TLE('NO', filedata(1:lnblnk(filedata)), linedoy(i), doyfile)
      if (swmetop)                    call init_TLE('ME', filedata(1:lnblnk(filedata)), linedoy(i), doyfile)
      if (swmn2)                      call init_TLE('M2', filedata(1:lnblnk(filedata)), linedoy(i), doyfile)
      if (swaqmod.and.(.not.swterra)) call init_TLE('AM', filedata(1:lnblnk(filedata)), linedoy(i), doyfile)
      if (swterra)                    call init_TLE('TE', filedata(1:lnblnk(filedata)), linedoy(i), doyfile)
c-- I have seen the doy in file on NOAA being wrong - correct with date in filename
      if (swnoaa) then
        if (doyfile.ne.linedoy(i)) then
          write (*,'(''ERROR doy in data  : '',I5,'' is NOT doy in filename : '',I5,'' using the latter ! '')') linedoy(i), doyfile
          doy_off = doyfile - linedoy(i)
          do j = 1, irec
            if (linevalid(j).eq.1) linedoy(j) = linedoy(j) + doy_off
          enddo
          call init_TLE('NO', filedata(1:lnblnk(filedata)), linedoy(i), doyfile)
        endif
      endif
      call sgp4_init(ro, vo, long, lat)
      if (.not.swmerge.or.swfirstmerge) call map_init()
      if (swgapcor.and.(.not.swmerge.or.swfirstmerge)) call map_stereo_init(swfast)
      if (swgapcor.and.(.not.swmerge.or.swfirstmerge)) call map_azel_init()
c-- Map out the subsatellite point and then the swath
c      swswath = .true.
      lat_start = 0.0
      lat_end   = 0.0
      altavg    = 0.0
      naltavg   = 0
      do i = 1, irec
        if (linevalid(i).eq.1) then
          call run_TLE(linetime(i),linedoy(i), ro, vo, long, lat, theta0g)
          call map_mark(long, lat, 0)
          call map_stereo_mark(long, lat, 0)
          if (lat_start.eq.0.0) lat_start = lat
          lat_end = lat
c-- Determine altitude 
          ro_copy(1) = ro(1) * 1000.0D0
          ro_copy(2) = ro(2) * 1000.0D0
          ro_copy(3) = ro(3) * 1000.0D0
          call xyz2lla(ro_copy, lla)
          altavg     = altavg  + lla(3)
          naltavg    = naltavg + 1
        endif
      enddo
      if (lat_start.lt.lat_end) swnorth = .true.
      if (lat_start.gt.lat_end) swnorth = .false.
      write (*,'(''Latitude start     : '',F8.2,'' Latitude end           : '',F8.2,'' swnorth : '',L1)') lat_start, lat_end, swnorth
      if (naltavg.gt.0) write (*,'(''Average altitude   : '',F8.2)') sngl(altavg)/float(naltavg)/1000.0
      if (swnorth.and.(.not.swthetascan)) thetascan = - thetascan
c-- Theta override
      if (swtheta.and.(.not.swthetascan)) thetascan = theta_in
      if (swthetascan) then
        thetascan = thetascan + theta_step
        if (thetascan.gt.theta_end) goto 101
      endif
      if (sweta) etascan = eta_in
c-- Second loop through the image(s) - write the channel images with or without gapcor and fill the satch_image (and related) taking gaps into account
      write (*,'(a,2F8.2)') 'Theta/Eta applied  : ', thetascan, etascan
      timprev       = 0.0D0
      nlines_gapcor = 0      
      irec          = 1
      do while (.true.)
        ios = 0
        if (swfy)    call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'FY', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
        if (swnoaa)  call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'NO', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
        if (swmetop) then
          if (swbin) then
            call get_scanline(lun, irec, ios, channel_image_bin_metop, 2071, nchan, nsat_pix, 'ME', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
            do i = 1, nchan
              do j = 1, nsat_pix
                channel_image(j,i) = channel_image_bin_metop(j+11,i)
              enddo
            enddo
          else
            call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'ME', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
          endif
        endif
        if (swmn2)   call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'M2', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
        if (swaqmod) call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'AM', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
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
                  nlines_gapcor = nlines_gapcor + 1
                  do j = 1,nchan
                    do i = 1,nsat_pix
                      imgbuf(i) = 0
                      satch_image(i,j,nlines_gapcor) = 0
                    enddo
c                    if (.not.swfast) write (lunchan(j),rec=nlines_gapcor) imgbuf
                    if (.not.swfast) call write_buffer(lunchan(j), nlines_gapcor, imgbuf, nsat_pix)
                  enddo
c-- Write the RGB image empty lines
                  if (swhybrid) then
                    do k = 2,nchan
                      do i = 1,nsat_pix
                        rgbbuf(1,i) = 0
                        rgbbuf(2,i) = 0
                        rgbbuf(3,i) = 0
                      enddo
c                      write (lunrgb(k),rec=nlines_gapcor) rgbbuf
                      call write_buffer(lunrgb(k), nlines_gapcor, rgbbuf, nsat_pix * 3)
                    enddo
                  endif
                  satch_linetime(nlines_gapcor) = timprev + dble(k3) * 0.1667D0
                  satch_linedoy(nlines_gapcor)  = linedoy(irec)
                  swswath                       = .false.
                  call run_TLE(satch_linetime(nlines_gapcor),satch_linedoy(nlines_gapcor), ro, vo, long, lat, theta0g)
                  aqmtheta(1) = 0.0D0
                  aqmtheta(2) = 0.0D0
                  if (.not.swaqmod) then
                    call scan_track(ro, vo, long, lat, satch_linetime(nlines_gapcor),satch_linedoy(nlines_gapcor), longlat(1,1,nlines_gapcor), swswath, theta0g, thetascan, etascan, swszalim, szalim, swnorth, swtheta, swdebug, swmn2, nsat_pix, sat_fov)
                  else
                    call scan_track_aqm(ro, vo, long, lat, satch_linetime(nlines_gapcor),satch_linedoy(nlines_gapcor), longlat(1,1,nlines_gapcor), swswath, theta0g, thetascan, etascan, swszalim, szalim, swnorth, swtheta, swdebug, swmn2, nsat_pix, sat_fov, aqmtheta, bowgamma)
                  endif
                enddo
              endif
              timprev = linetime(irec)
            endif
            nlines_gapcor = nlines_gapcor + 1
            do j = 1,nchan
              do i = 1,nsat_pix
                if (swnorth)      l = nsat_pix + 1 - i
                if (.not.swnorth) l = i
                imgbuf(l) = min(max(channel_image(i,j),0),iclim)
                hist(imgbuf(l)+1,j) = hist(imgbuf(l)+1,j) + 1
                satch_image(l,j,nlines_gapcor) = imgbuf(l)
              enddo
               if (.not.swfast) call write_buffer(lunchan(j), nlines_gapcor, imgbuf, nsat_pix)
            enddo
c-- Write the RGB image
            if (swhybrid) then
              do k = 2,nchan
                do i = 1,nsat_pix
                  k1 = 0
                  k2 = 0
                  do j = 1,nchan
                    if (j.eq.1) k1 = (channel_image(i,j) / 4) + 1
                    if (j.eq.k) k2 = (channel_image(i,j) / 4) + 1
                  enddo
                  if (swnorth)      l = nsat_pix + 1 - i
                  if (.not.swnorth) l = i
                  rgbbuf(1,l) = lut(1,k2,k1)
                  rgbbuf(2,l) = lut(2,k2,k1)
                  rgbbuf(3,l) = lut(3,k2,k1)
                enddo
c                write (lunrgb(k),rec=nlines_gapcor) rgbbuf
                call write_buffer(lunrgb(k), nlines_gapcor, rgbbuf, nsat_pix * 3)
              enddo
            endif
            satch_linetime(nlines_gapcor) = linetime(irec)
            satch_linedoy(nlines_gapcor)  = linedoy(irec)
            satch_valid(nlines_gapcor)    = 1
            swswath                       = .true.
            call run_TLE(satch_linetime(nlines_gapcor),satch_linedoy(nlines_gapcor), ro, vo, long, lat, theta0g)
            if (.not.swaqmod) then
              call scan_track(ro, vo, long, lat, satch_linetime(nlines_gapcor),satch_linedoy(nlines_gapcor), longlat(1,1,nlines_gapcor), swswath, theta0g, thetascan, etascan, swszalim, szalim, swnorth, swtheta, swdebug, swmn2, nsat_pix, sat_fov)
            else
              call scan_track_aqm(ro, vo, long, lat, satch_linetime(nlines_gapcor),satch_linedoy(nlines_gapcor), longlat(1,1,nlines_gapcor), swswath, theta0g, thetascan, etascan, swszalim, szalim, swnorth, swtheta, swdebug, swmn2, nsat_pix, sat_fov, aqmtheta, bowgamma)
            endif
            call azel(long, lat, ro, satch_az(nlines_gapcor), satch_el(nlines_gapcor))
            call map_azel_place(satch_az(nlines_gapcor), satch_el(nlines_gapcor))
            if (swbatch.and.swhorizon) call map_azel_horizon_place(satch_az(nlines_gapcor), satch_el(nlines_gapcor))
            call map_stereo_place_azel(satch_az(nlines_gapcor), satch_el(nlines_gapcor))
          endif
        else
          do j = 1,nchan
            do i = 1,nsat_pix
              if (swnorth)      l = nsat_pix + 1 - i
              if (.not.swnorth) l = i
              imgbuf(l) = min(max(channel_image(i,j),0),iclim)
c              if (0.le.imgbuf(l).and.imgbuf(l).le.(nsat_pix/2) - 1) then
c                imgbuf(l) = channel_image(i,j)
c              else
c                imgbuf(l) = 0
c              endif
              hist(imgbuf(l)+1,j) = hist(imgbuf(l)+1,j) + 1
            enddo
c            if (.not.swfast) write (lunchan(j),rec=irec) imgbuf
            if (.not.swfast) call write_buffer(lunchan(j), irec, imgbuf, nsat_pix)
          enddo
c-- Write the RGB image
          if (swhybrid) then
            do k = 2,nchan
              do i = 1,nsat_pix
                k1 = 0
                k2 = 0
                do j = 1,nchan
                  if (j.eq.1) k1 = (channel_image(i,j) / 4) + 1
                  if (j.eq.k) k2 = (channel_image(i,j) / 4) + 1
                enddo
                if (swnorth)      l = nsat_pix + 1 - i
                if (.not.swnorth) l = i
                rgbbuf(1,l) = lut(1,k2,k1)
                rgbbuf(2,l) = lut(2,k2,k1)
                rgbbuf(3,l) = lut(3,k2,k1)
              enddo
c              write (lunrgb(k),rec=irec) rgbbuf
              call write_buffer(lunrgb(k), irec, rgbbuf, nsat_pix * 3)
            enddo
          endif
        endif
c-- Create the histogram for an RGB image with R and G = Ch 2 and B = Ch 1 - 'The @petermeteor option'
        do i = 1,nsat_pix
          dr = float(channel_image(i,ir))
          dg = float(channel_image(i,ig))
          db = float(channel_image(i,ib))
c-- Transform to YIQ (The Y is used as Luminance in analogue colour TV) and build histogram on Y
          dy = 0.299 * dr + 0.587 * dg + 0.114 * db
          k = int(dy)
          k = min(max(k,0),iclim) + 1
          histrgb(k) = histrgb(k) + 1
        enddo
c-- Create the grayscale histograms per vis channel - For FY3 only
        if (swfy) then
          do i = 1,nsat_pix
            do j = 1, nchan
              do l = 1, nvisch
                if (j.eq.visch(l)) then
                  k = channel_image(i,j)
                  k = min(max(k,0),iclim) + 1
                  histch(l,k) = histch(l,k) + 1
                endif
              enddo
            enddo
          enddo
        endif
        irec = irec + 1
      enddo
 1    continue
      if (swdebug) then
        write (*,*) '** Azimuth & Elevation for observer **'
        do i = 1, nlines_gapcor
          if (satch_valid(i).eq.1) write (*,*) i, satch_linetime(i), longlat(1,(nsat_pix/2)+1,i), longlat(2,(nsat_pix/2)+1,i), satch_az(i), satch_el(i)
        enddo
      endif
c-- Determine the 0.5% and 99.5 % points in the histograms for the channels to be used in a truecolor image - For FY3 only
      if (swfy) then
        do j = 1, nvisch
          sum1 = 0.
          do k = 1,iclim+1
            sum1 = sum1 + float(histch(j,k))
          enddo
          sum2   = 0.
          jch1(j) = -1
          jch2(j) = -1
          do k = 1,iclim+1
            sum2 = sum2 + float(histch(j,k))
            if (sum2.ge.0.005*sum1.and.jch1(j).eq.-1) jch1(j) = k
            if (sum2.ge.0.995*sum1.and.jch2(j).eq.-1) jch2(j) = k
          enddo
          if (swdebug) write (*,*) jch1(j), jch2(j), histch(j,jch1(j)), histch(j,jch2(j)), visch(j)
        enddo
      endif
c-- Determine the 0.5% and 99.5 % points in the Y-histogram (Y from YIQ) for the RGB image created from ch 2, 2, 1
      sum1 = 0.
      do k = 1,iclim+1
        sum1 = sum1 + float(histrgb(k))
      enddo
      sum2 = 0.
      j1   = -1
      j2   = -1
      do k = 1,iclim+1
        sum2 = sum2 + float(histrgb(k))
        if (sum2.ge.0.005*sum1.and.j1.eq.-1) j1 = k
        if (sum2.ge.0.995*sum1.and.j2.eq.-1) j2 = k
      enddo
c
c-- Third loop through the image(s)
      timprev       = 0.0D0
      nlines_gapcor = 0      
      irec = 1
      do while (.true.)
        ios = 0
        if (swfy)    call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'FY', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
        if (swnoaa)  call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'NO', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
        if (swmetop) then
          if (swbin) then
            call get_scanline(lun, irec, ios, channel_image_bin_metop, 2071, nchan, nsat_pix, 'ME', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
            do i = 1, nchan
              do j = 1, nsat_pix
                channel_image(j,i) = channel_image_bin_metop(j+11,i)
              enddo
            enddo
          else
            call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'ME', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
          endif
        endif
        if (swmn2)   call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'M2', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
        if (swaqmod) call get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, nsat_pix, 'AM', inbuf, nrecl, timestamp, doy, swdebug, swbin, aqmtheta)
        if (ios.ne.0) goto 3
c
        if (swgapcor) then
          if (linevalid(irec).eq.1) then
            if (timprev.eq.0.0D0) then
              timprev = linetime(irec)
            else
              nblank = nint((linetime(irec)-timprev)/0.1667D0) - 1
              if (nblank.gt.0) then
                do k3 = 1, nblank
                  nlines_gapcor = nlines_gapcor + 1
                  do i = 1,nsat_pix
                    rgbbuf(1,i) = 0
                    rgbbuf(2,i) = 0
                    rgbbuf(3,i) = 0
                  enddo
c                  write (lunrgb(1),rec=nlines_gapcor) rgbbuf
                  call write_buffer(lunrgb(1), nlines_gapcor, rgbbuf, nsat_pix * 3)
c-- Handle the corrected CH11 image
                  call correct_apply(rgbbuf, 3, nsat_pix, rgbbuf_correct, 3, 2*ncorrect_pix, index_correct, ncorrect_pix, lun_correct, nlines_gapcor)
c-- Handle the truecolor image - For FY3 only
                  if (swfy) then
                    do i = 1,nsat_pix
                      rgbbuf(1,i) = 0
                      rgbbuf(2,i) = 0
                      rgbbuf(3,i) = 0
                    enddo
c                    write (luntrue,rec=nlines_gapcor) rgbbuf
                    call write_buffer(luntrue, nlines_gapcor, rgbbuf, nsat_pix * 3)
c-- Handle the corrected truecolor image
                   call correct_apply(rgbbuf, 3, nsat_pix, rgbbuf_correct, 3, 2*ncorrect_pix, index_correct, ncorrect_pix, luntrue_correct, nlines_gapcor)
                  endif
                enddo
              endif
              timprev = linetime(irec)
            endif
            nlines_gapcor = nlines_gapcor + 1
            do i = 1,nsat_pix
c-- Handle the offsets of the requested RGB channels - the formula is the same for N and S passes as it is an offset in scan (i) space
              k1 = i - xchoff(ir)
              k1 = max(min(nsat_pix,k1),1)
              dr = float(channel_image(k1,ir))
c
              k1 = i - xchoff(ig)
              k1 = max(min(nsat_pix,k1),1)
              dg = float(channel_image(k1,ig))
c
              k1 = i - xchoff(ib)
              k1 = max(min(nsat_pix,k1),1)
              db = float(channel_image(k1,ib))
c
              if (swnorth)      l = nsat_pix + 1 - i
              if (.not.swnorth) l = i
              rgbbuf(1,l) = min(max(int( (rclim + 1.0) * (dr-float(j1))/float(j2-j1+1)),0),iclim)
              rgbbuf(2,l) = min(max(int( (rclim + 1.0) * (dg-float(j1))/float(j2-j1+1)),0),iclim)
              rgbbuf(3,l) = min(max(int( (rclim + 1.0) * (db-float(j1))/float(j2-j1+1)),0),iclim)
            enddo
c            write (lunrgb(1),rec=nlines_gapcor) rgbbuf
            call write_buffer(lunrgb(1), nlines_gapcor, rgbbuf, nsat_pix * 3)
c-- Handle the corrected CH11 image
            call correct_apply(rgbbuf, 3, nsat_pix, rgbbuf_correct, 3, 2*ncorrect_pix, index_correct, ncorrect_pix, lun_correct, nlines_gapcor)
c-- Handle the truecolor image - For FY3 only
            if (swfy) then
              do i = 1,nsat_pix
                if (swnorth)      k = nsat_pix + 1 - i
                if (.not.swnorth) k = i
                rgbbuf(1,k) = 0
                rgbbuf(2,k) = 0
                rgbbuf(3,k) = 0
                do j = 1,nchan
                  if (swtruech) then
c-- Different approach possible R = ch1 G = (ch8 + ch9) / 2  B = ch7
                    l = 1
                    if (j.eq.visch(l)) then
                      k1 = i - xchoff(visch(l))
                      k1 = max(min(nsat_pix,k1),1)
                      dr = float(channel_image(k1,j) + 1)
                      if (jch1(l).le.int(dr).and.int(dr).le.jch2(l)) then
                        rgbbuf(1,k) = rgbbuf(1,k) + int(dr)
                      endif
                    endif
                    l = 2
                    if (j.eq.visch(l)) then
                      k1 = i - xchoff(visch(l))
                      k1 = max(min(nsat_pix,k1),1)
                      dr = float(channel_image(k1,j) + 1)
                      if (jch1(l).le.int(dr).and.int(dr).le.jch2(l)) then
                        rgbbuf(3,k) = rgbbuf(3,k) + int(dr)
                      endif
                    endif
                    l = 3
                    if (j.eq.visch(l)) then
                      k1 = i - xchoff(visch(l))
                      k1 = max(min(nsat_pix,k1),1)
                      dr = float(channel_image(k1,j) + 1)
                      if (jch1(l).le.int(dr).and.int(dr).le.jch2(l)) then
                        rgbbuf(2,k) = rgbbuf(2,k) + int(0.5 * dr)
                      endif
                    endif
                    l = 4
                    if (j.eq.visch(l)) then
                      k1 = i - xchoff(visch(l))
                      k1 = max(min(nsat_pix,k1),1)
                      dr = float(channel_image(k1,j) + 1)
                      if (jch1(l).le.int(dr).and.int(dr).le.jch2(l)) then
                        rgbbuf(2,k) = rgbbuf(2,k) + int(0.5 * dr)
                      endif
                    endif
c-- End of different approach
                  else
                    do l = 1, nvisch
                      if (j.eq.visch(l)) then
c-- The formula is the same for N and S passes as it is an offset in scan (i) space
                        k1 = i - xchoff(visch(l))
                        k1 = max(min(nsat_pix,k1),1)
                        dr = float(channel_image(k1,j) + 1)
                        if (jch1(l).le.int(dr).and.int(dr).le.jch2(l)) then
                          rgbbuf(1,k) = rgbbuf(1,k) + int(float(rgb(1,l)) * (dr-float(jch1(l)))/float(jch2(l)-jch1(l)+1))
                          rgbbuf(2,k) = rgbbuf(2,k) + int(float(rgb(2,l)) * (dr-float(jch1(l)))/float(jch2(l)-jch1(l)+1))
                          rgbbuf(3,k) = rgbbuf(3,k) + int(float(rgb(3,l)) * (dr-float(jch1(l)))/float(jch2(l)-jch1(l)+1))
                        endif
                      endif
                    enddo
                  endif
                enddo
                if (.not.swtruech) then
                  rgbbuf(1,k) = int(float(rgbbuf(1,k)) / rsum)
                  rgbbuf(2,k) = int(float(rgbbuf(2,k)) / gsum)
                  rgbbuf(3,k) = int(float(rgbbuf(3,k)) / bsum)
                endif
              enddo
c              write (luntrue,rec=nlines_gapcor) rgbbuf
              call write_buffer(luntrue, nlines_gapcor, rgbbuf, nsat_pix * 3)
c-- Handle the corrected truecolor image
              call correct_apply(rgbbuf, 3, nsat_pix, rgbbuf_correct, 3, 2*ncorrect_pix, index_correct, ncorrect_pix, luntrue_correct, nlines_gapcor)
            endif
          endif
        else
          do i = 1,nsat_pix
c-- Handle the offsets of the requested RGB channels - the formula is the same for N and S passes as it is an offset in scan (i) space
            k1 = i - xchoff(ir)
            k1 = max(min(nsat_pix,k1),1)
            dr = float(channel_image(k1,ir))
c
            k1 = i - xchoff(ig)
            k1 = max(min(nsat_pix,k1),1)
            dg = float(channel_image(k1,ig))
c
            k1 = i - xchoff(ib)
            k1 = max(min(nsat_pix,k1),1)
            db = float(channel_image(k1,ib))
c
            if (swnorth)      l = nsat_pix + 1 - i
            if (.not.swnorth) l = i
            rgbbuf(1,l) = min(max(int( (rclim + 1.0) * (dr-float(j1))/float(j2-j1+1)),0),iclim)
            rgbbuf(2,l) = min(max(int( (rclim + 1.0) * (dg-float(j1))/float(j2-j1+1)),0),iclim)
            rgbbuf(3,l) = min(max(int( (rclim + 1.0) * (db-float(j1))/float(j2-j1+1)),0),iclim)
          enddo
c          write (lunrgb(1),rec=irec) rgbbuf
          call write_buffer(lunrgb(1), irec, rgbbuf, nsat_pix * 3)
c-- Handle the corrected CH11 image
          call correct_apply(rgbbuf, 3, nsat_pix, rgbbuf_correct, 3, 2*ncorrect_pix, index_correct, ncorrect_pix, lun_correct, irec)
c-- Handle the truecolor image - For FY3 only
          if (swfy) then
            do i = 1,nsat_pix
              if (swnorth)      k = nsat_pix + 1 - i
              if (.not.swnorth) k = i
              rgbbuf(1,k) = 0
              rgbbuf(2,k) = 0
              rgbbuf(3,k) = 0
              do j = 1,nchan
                do l = 1, nvisch
                  if (j.eq.visch(l)) then
c-- The formula is the same for N and S passes as it is an offset in scan (i) space
                    k1 = i - xchoff(visch(l))
                    k1 = max(min(nsat_pix,k1),1)
                    dr = float(channel_image(k1,j) + 1)
                    if (jch1(l).le.int(dr).and.int(dr).le.jch2(l)) then
                      rgbbuf(1,k) = rgbbuf(1,k) + int(float(rgb(1,l)) * (dr-float(jch1(l)))/float(jch2(l)-jch1(l)+1))
                      rgbbuf(2,k) = rgbbuf(2,k) + int(float(rgb(2,l)) * (dr-float(jch1(l)))/float(jch2(l)-jch1(l)+1))
                      rgbbuf(3,k) = rgbbuf(3,k) + int(float(rgb(3,l)) * (dr-float(jch1(l)))/float(jch2(l)-jch1(l)+1))
                    endif
                  endif
                enddo
              enddo
              rgbbuf(1,k) = int(float(rgbbuf(1,k)) / rsum)
              rgbbuf(2,k) = int(float(rgbbuf(2,k)) / gsum)
              rgbbuf(3,k) = int(float(rgbbuf(2,k)) / bsum)
c              rgbbuf(1,k) = int(float(rgbbuf(1,k)) )
c              rgbbuf(2,k) = int(float(rgbbuf(2,k)) )
c              rgbbuf(3,k) = int(float(rgbbuf(2,k)) )
            enddo
            write (luntrue,rec=irec) rgbbuf
            call write_buffer(luntrue, irec, rgbbuf, nsat_pix * 3)
c-- Handle the corrected truecolor image
            call correct_apply(rgbbuf, 3, nsat_pix, rgbbuf_correct, 3, 2*ncorrect_pix, index_correct, ncorrect_pix, luntrue_correct, irec)
          endif
        endif
        irec = irec + 1
      enddo
 3    irec = irec - 1
c-- Release the logical Unit nrs = lun
      close (unit=lun)
      call free_lun(lun)
      do i = 1,nchan
        if (i.eq.1.or.swhybrid) then
          close (unit=lunrgb(i))
          call free_lun(lunrgb(i))
        endif
        if (.not.swfast) then
          close (unit=lunchan(i))
          call free_lun(lunchan(i))
        endif
      enddo
      close (unit=lun_correct)
      call free_lun(lun_correct)
      if (swfy) then
        close (unit=luntrue)
        call free_lun(luntrue)
        close (unit=luntrue_correct)
        call free_lun(luntrue_correct)
      endif
c
c-- Find the zoom parameters for the polar stereographic projection and execute the projection
c-- First zero out the image area
c
      do i = 1, nstereo
        stereo_image(1,i) = 0
        stereo_image(2,i) = 0
        stereo_image(3,i) = 0
        do j = 1, 3
          stereo_count(  i,j) = 0
        enddo
      enddo
c-- Do the hist normalisation before the merge jump !
c-- Try an on request histogram limits (user configurable) renormalisation of satch_image to 0 .. iclim
      if (swproject.and.swhistcor) then
        do j = 1, nchan
          sum1 = 0.
          do k = 1,iclim+1
            sum1 = sum1 + float(hist(k,j))
          enddo
          sum2 = 0.
          j1   = -1
          j2   = -1
          do k = 1,iclim+1
            sum2 = sum2 + float(hist(k,j))
            if (sum2.ge.0.005*sum1.and.j1.eq.-1) j1 = k
            if (sum2.ge.0.995*sum1.and.j2.eq.-1) j2 = k
          enddo
          j1 = j1 - 1
          j2 = j2 - 1
          do i = 1, nlines_gapcor
            if (satch_valid(i).eq.1) then
              do k = 1,nsat_pix
                satch_image(k,j,i) = min(max(int(float(satch_image(k,j,i) - j1) * float(iclim+1) / float(j2-j1+1)),0),iclim)
              enddo
            endif
          enddo
        enddo
      endif
c--
c-- For merging, start here with the disk stored files, but for all but the last loop skip this ! (after the last file the reloading starts)
c-- For merge, save nsat_pix, nlines_gapcor, satch_valid, swnorth, longlat, satch_image
      if (swmerge.and.(.not.swlastmerge)) then
        nsat_pix_merge(i_merge)      = nsat_pix
        nlines_gapcor_merge(i_merge) = nlines_gapcor
        swnorth_merge(i_merge)       = swnorth
        call dump_merge(i_merge, satch_image, satch_valid, longlat, nmaxpix, nmaxch, nlines_gapcor, lundump, nmergemax)
        goto 300
      endif
      if (swmerge.and.swlastmerge) then
        nsat_pix_merge(i_merge)      = nsat_pix
        nlines_gapcor_merge(i_merge) = nlines_gapcor
        swnorth_merge(i_merge)       = swnorth
      endif
      if (swproject) then
        call map_stereo_zoom_init(nsat_pix,irec,nxs,nys,swdebug,iclim)
        if (nxs*nys.gt.nstereo) then 
          write (*,*) ' ** Enlarge stereo Image in source code (nstereo) ', nxs*nys, ' at least required'
          stop
        endif
        write (*,'(''Projected nxs, nys : '',2I10)') nxs, nys
c-- If border colour or longlat colour are user defined, they can only be set after the map_stereo_zoom_init
c-- So repeat the getarg code from the top here and call the set_ routines
        if (swbcol) then
          call getarg(i_bcol,argstring)
          i  = index(argstring,'=')
          j  = index(argstring,',')
          call string_to_i4(argstring(i+1:j-1)//' ',npos, mv)
          mr = min(max(mv,0),iclim)
          i  = j
          j  = index(argstring(j+1:),',') + j
          call string_to_i4(argstring(i+1:j-1)//' ',npos, mv)
          mg = min(max(mv,0),iclim)
          i  = j
          call string_to_i4(argstring(i+1:lnblnk(argstring))//' ',npos, mv)
          mb = min(max(mv,0),iclim)
          call set_map_zoom_mark_colour(mr, mg, mb)
        endif
        if (swllcol) then
          call getarg(i_llcol,argstring)
          i  = index(argstring,'=')
          j  = index(argstring,',')
          call string_to_i4(argstring(i+1:j-1)//' ',npos, mv)
          mr = min(max(mv,0),iclim)
          i  = j
          j  = index(argstring(j+1:),',') + j
          call string_to_i4(argstring(i+1:j-1)//' ',npos, mv)
          mg = min(max(mv,0),iclim)
          i  = j
          call string_to_i4(argstring(i+1:lnblnk(argstring))//' ',npos, mv)
          mb = min(max(mv,0),iclim)
          call set_map_zoom_mark_longlat_colour(mr, mg, mb)
        endif
c-- End of border or lonlat colour user defined settings
      endif
c-- Create a Solar Zenith Angle histogram - this can only be done if projection is requested and swszalim = true
      if (swszalim) then
        do i = 1, nlines_gapcor
          if (satch_valid(i).eq.1) then
            do k = 1,nsat_pix
              i_sza = int(longlat(3,k,i)/0.1) + 1
              if (1.le.i_sza.and.i_sza.le.n_sza) then
                szahist(i_sza)   = szahist(i_sza)   + 1
              endif
            enddo
          endif
        enddo
        szasum = 0.
        i      = 0
        do k = 1, n_sza
          szasum = szasum + float(szahist(k)) * float(k-1) * 0.1
          i = i + szahist(k)
        enddo
        if (swbatch.and.i.gt.0) then
          if (szasum/float(i).gt.82.0) then
            write (*,'(''SZA average        : '',F10.2,'' switching to project=445'')') szasum/float(i)
            irps = irp
            igps = igp
            ibps = ibp
            irp  = 4
            igp  = 4
            ibp  = 5
            swbatchnight = .true.
          endif
        endif
      endif
      if (swproject) then
c-- Precursor logic for merge with IR channel in case of part of the swath 'in the dark'
c-- Caculate average luminosity in 2 deg before SZA limit
c-- Calculate average luminosity in 2 deg after SZA limit
c-- Use ratio as a normalisation
        if (swszalim.and.swszamer) then
          sum1  = 0.
          npix1 = 0
          perc1 = 0
          sum2  = 0.
          npix2 = 0
          perc2 = 0
          do i = 1, iclim+1
            hist1(i) = 0
            hist2(i) = 0
          enddo
          do i = 1, nlines_gapcor
            if (satch_valid(i).eq.1) then
              do k = 1,nsat_pix
c-- Daylight part. Translate RGB to YIQ and use Y
                if (szalim-5.0.le.longlat(3,k,i).and.longlat(3,k,i).le.szalim-0.5) then
                  dr = float(satch_image(k,irp,i))
                  dg = float(satch_image(k,igp,i))
                  db = float(satch_image(k,ibp,i))
                  dy = 0.299 * dr + 0.587 * dg + 0.114 * db
                  sum1  = sum1  + dy
                  npix1 = npix1 + 1
                  k1 = int(dy)
                  k1 = min(max(k1,0),iclim) + 1
                  hist1(k1) = hist1(k1) + 1
                endif
c-- Nightside part - use directly
                if (szalim+0.5.le.longlat(3,k,i).and.longlat(3,k,i).le.szalim+5.0) then
                  sum2  = sum2 + satch_image(k,szach,i)
                  npix2 = npix2 + 1
                  k2 = satch_image(k,szach,i)
                  k2 = min(max(k2,0),iclim) + 1
                  hist2(k2) = hist2(k2) + 1
                endif
              enddo
            endif
          enddo
          k = 0
          i = 1
          do while ((float(k)/float(npix1)).le.0.005.and.i.le.iclim+1)
           k = k + hist1(i)
           i = i + 1
          enddo
          perc1 = i
          k = 0
          i = 1
          do while ((float(k)/float(npix2)).le.0.005.and.i.le.iclim+1)
           k = k + hist2(i)
           i = i + 1
          enddo
          perc2 = i
c-- colour night = (satch(_image(k,szach,i) - (perc2-perc1)) * ( ((sum1 / float(npix1)) - perc1) / ((sum2 / float(npix2)) - perc2) )
          w3 =  szafact * ((sum1/float(npix1)) - float(perc1)) / ((sum2/float(npix2)) - float(perc2))
c          write (*,*) perc1, perc2, w3, ((sum1/float(npix1)) - float(perc1)), ((sum2/float(npix2)) - float(perc2)), npix1, npix2, szach, sum1/float(npix1), sum2/float(npix2)
        endif
        i_off = 1
        if (swyiq) call map_IQ_init()
        wop1 = 1.0 - (float(opperc) / 100.0)
        wop2 =       (float(opperc) / 100.0)
c-- This statement is used to merge back in the stored passes to be merged - in reverse order
 302    continue
        do i = 1, nlines_gapcor
          if (satch_valid(i).eq.1) then
            do k = 1,nsat_pix
c-- Calculate the fractional position (as seen from the center of the scanline
              rpos = abs((dble(2*k)/dble(nsat_pix)) - 1.0D0)
c-- Take care of some of the bowtie correction latitude sequence subtleties to prevent erroneous interpolation
              swok = .true.
c-- The rpos = 0.35 is an experimental value that closes the "gaps" at the 20 pixel boundaries in the 500 m image as a consequence of the bowtie latitude jump filter below
c-- Still don't understand it
              if (swaqmod) then
                if (swnorth       .and.(longlat(2,k,i  ).lt.longlat(2,k,i-1).or.longlat(2,k+1,i  ).lt.longlat(2,k+1,i-1)).and.rpos.gt.0.35) swok = .false.
                if (swnorth       .and.(longlat(2,k,i+1).lt.longlat(2,k,i  ).or.longlat(2,k+1,i+1).lt.longlat(2,k+1,i  )).and.rpos.gt.0.35) swok = .false.
                if ((.not.swnorth).and.(longlat(2,k,i  ).gt.longlat(2,k,i-1).or.longlat(2,k+1,i  ).gt.longlat(2,k+1,i-1)).and.rpos.gt.0.35) swok = .false.
                if ((.not.swnorth).and.(longlat(2,k,i+1).gt.longlat(2,k,i  ).or.longlat(2,k+1,i+1).gt.longlat(2,k+1,i  )).and.rpos.gt.0.35) swok = .false.
c-- The usual ndiv to be used for HRPT and inner part of AQUA Modis
                ndivuse = 4
c-- The ndiv to be used for the middle part of AQUA Modis
                if (rpos.gt.0.30) ndivuse = 6
c-- The ndiv to be used for the outermost part of AQUA Modis
                if (rpos.gt.0.80) ndivuse = 8
c-- The ndiv to be used for the outer-outermost part of AQUA Modis :-)
                if (rpos.gt.0.90) ndivuse = 10
              else
c-- The usual ndiv to be used for inner part of HRPT
                ndivuse = 4
c-- The ndiv to be used for the outer part of HRPT
                if (rpos.gt.0.70) ndivuse = 6
              endif
c-- Uncomment the following line if you simply want to see a projected image of the ORIGINAL pixels - without interpolation etc etc
c              ndivuse = 1
c-- Sanity check to prevent crashing if ndivuse > ndivgrid
              ndivuse = min(ndivuse, ndivgrid)
c
c-- First the code to take care of the merging of an IR channel based on Solar Zenith Angle
c
              if (swszalim.and.swszamer) then
                if (longlat(3,k,i).le.szalim-0.5) then
                  if (k.eq.1.or.k.eq.nsat_pix.or.i.eq.1.or.i.eq.nlines_gapcor) then
                  else
                    call grid_sample(satch_image(k-1,irp,i-1),satch_image(k  ,irp,i-1),satch_image(k+1,irp,i-1),
     *                               satch_image(k-1,irp,i)  ,satch_image(k  ,irp,i)  ,satch_image(k+1,irp,i)  ,
     *                               satch_image(k-1,irp,i+1),satch_image(k  ,irp,i+1),satch_image(k+1,irp,i+1),
     *                               longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valr, satch_ll, ndivuse)
                    call grid_sample(satch_image(k-1,igp,i-1),satch_image(k  ,igp,i-1),satch_image(k+1,igp,i-1),
     *                               satch_image(k-1,igp,i)  ,satch_image(k  ,igp,i)  ,satch_image(k+1,igp,i)  ,
     *                               satch_image(k-1,igp,i+1),satch_image(k  ,igp,i+1),satch_image(k+1,igp,i+1),
     *                               longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valg, satch_ll, ndivuse)
                    call grid_sample(satch_image(k-1,ibp,i-1),satch_image(k  ,ibp,i-1),satch_image(k+1,ibp,i-1),
     *                               satch_image(k-1,ibp,i)  ,satch_image(k  ,ibp,i)  ,satch_image(k+1,ibp,i)  ,
     *                               satch_image(k-1,ibp,i+1),satch_image(k  ,ibp,i+1),satch_image(k+1,ibp,i+1),
     *                               longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valb, satch_ll, ndivuse)
                    if (swok) then
                      do l = 1, ndivuse
                        do j = 1, ndivuse
                          if (swmerge.and.swaqmpan) then
                            if (i_merge.eq.3) then
                              pixcol(1) = 0
                              pixcol(2) = 0
                              pixcol(3) = fetch_pixel(satch_valb, ndivuse, j, l)
                            else if (i_merge.eq.2) then
                              pixcol(1) = 0
                              pixcol(2) = fetch_pixel(satch_valg, ndivuse, j, l)
                              pixcol(3) = 0
                            else if (i_merge.eq.1) then
                              pixcol(1) = fetch_pixel(satch_valr, ndivuse, j, l)
                              pixcol(2) = 0
                              pixcol(3) = 0
                            endif
                          else
                            pixcol(1) = fetch_pixel(satch_valr, ndivuse, j, l)
                            pixcol(2) = fetch_pixel(satch_valg, ndivuse, j, l)
                            pixcol(3) = fetch_pixel(satch_valb, ndivuse, j, l)
                          endif
                          long_stereo = fetch_ll(satch_ll, 2, ndivuse, 1, j, l)
                          lat_stereo  = fetch_ll(satch_ll, 2, ndivuse, 2, j, l)
                          if (swaqmod) then
                            if (swaqmpan) then
                              call map_stereo_zoom_place_aqm_pan(stereo_image,stereo_count,3,nxs,nys,long_stereo, lat_stereo, pixcol,i_merge)
                            else
                              call map_stereo_zoom_place_aqm(stereo_image,stereo_count,3,nxs,nys,long_stereo, lat_stereo, pixcol)
                            endif
                          else
                            call map_stereo_zoom_place(stereo_image,3,nxs,nys,long_stereo, lat_stereo, pixcol)
                          endif
                        enddo
                      enddo
                    endif
                  endif
                else
                  if (longlat(3,k,i).le.szalim+0.5) then
                    w1 = longlat(3,k,i) - (szalim-0.5)
                    w2 = 1.0D0 - w1
                    if (k.eq.1.or.k.eq.nsat_pix.or.i.eq.1.or.i.eq.nlines_gapcor) then
                    else
                      call grid_sample(satch_image(k-1,szach,i-1),satch_image(k  ,szach,i-1),satch_image(k+1,szach,i-1),
     *                                 satch_image(k-1,szach,i)  ,satch_image(k  ,szach,i)  ,satch_image(k+1,szach,i)  ,
     *                                 satch_image(k-1,szach,i+1),satch_image(k  ,szach,i+1),satch_image(k+1,szach,i+1),
     *                                 longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_szar, satch_ll, ndivuse)
                      call grid_sample(satch_image(k-1,szach,i-1),satch_image(k  ,szach,i-1),satch_image(k+1,szach,i-1),
     *                                 satch_image(k-1,szach,i)  ,satch_image(k  ,szach,i)  ,satch_image(k+1,szach,i)  ,
     *                                 satch_image(k-1,szach,i+1),satch_image(k  ,szach,i+1),satch_image(k+1,szach,i+1),
     *                                 longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_szag, satch_ll, ndivuse)
                      call grid_sample(satch_image(k-1,szach,i-1),satch_image(k  ,szach,i-1),satch_image(k+1,szach,i-1),
     *                                 satch_image(k-1,szach,i)  ,satch_image(k  ,szach,i)  ,satch_image(k+1,szach,i)  ,
     *                                 satch_image(k-1,szach,i+1),satch_image(k  ,szach,i+1),satch_image(k+1,szach,i+1),
     *                                 longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_szab, satch_ll, ndivuse)
                      call grid_sample(satch_image(k-1,irp,i-1),satch_image(k  ,irp,i-1),satch_image(k+1,irp,i-1),
     *                                 satch_image(k-1,irp,i)  ,satch_image(k  ,irp,i)  ,satch_image(k+1,irp,i)  ,
     *                                 satch_image(k-1,irp,i+1),satch_image(k  ,irp,i+1),satch_image(k+1,irp,i+1),
     *                                 longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valr, satch_ll, ndivuse)
                      call grid_sample(satch_image(k-1,igp,i-1),satch_image(k  ,igp,i-1),satch_image(k+1,igp,i-1),
     *                                 satch_image(k-1,igp,i)  ,satch_image(k  ,igp,i)  ,satch_image(k+1,igp,i)  ,
     *                                 satch_image(k-1,igp,i+1),satch_image(k  ,igp,i+1),satch_image(k+1,igp,i+1),
     *                                 longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valg, satch_ll, ndivuse)
                      call grid_sample(satch_image(k-1,ibp,i-1),satch_image(k  ,ibp,i-1),satch_image(k+1,ibp,i-1),
     *                                 satch_image(k-1,ibp,i)  ,satch_image(k  ,ibp,i)  ,satch_image(k+1,ibp,i)  ,
     *                                 satch_image(k-1,ibp,i+1),satch_image(k  ,ibp,i+1),satch_image(k+1,ibp,i+1),
     *                                 longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valb, satch_ll, ndivuse)
                      if (swok) then
                        do l = 1, ndivuse
                          do j = 1, ndivuse
                            if (swmerge.and.swaqmpan) then
                              if (i_merge.eq.3) then
                                pixcol(1) = 0
                                pixcol(2) = 0
                                pixcol(3) = w1 * (w3 * (fetch_pixel(satch_szab, ndivuse, j, l) - (perc2-perc1))) + w2 * fetch_pixel(satch_valb, ndivuse, j, l)
                              else if (i_merge.eq.2) then
                                pixcol(1) = 0
                                pixcol(2) = w1 * (w3 * (fetch_pixel(satch_szag, ndivuse, j, l) - (perc2-perc1))) + w2 * fetch_pixel(satch_valg, ndivuse, j, l)
                                pixcol(3) = 0
                              else if (i_merge.eq.1) then
                                pixcol(1) = w1 * (w3 * (fetch_pixel(satch_szar, ndivuse, j, l) - (perc2-perc1))) + w2 * fetch_pixel(satch_valr, ndivuse, j, l)
                                pixcol(2) = 0
                                pixcol(3) = 0
                              endif
                            else
                              pixcol(1) = w1 * (w3 * (fetch_pixel(satch_szar, ndivuse, j, l) - (perc2-perc1))) + w2 * fetch_pixel(satch_valr, ndivuse, j, l)
                              pixcol(2) = w1 * (w3 * (fetch_pixel(satch_szag, ndivuse, j, l) - (perc2-perc1))) + w2 * fetch_pixel(satch_valg, ndivuse, j, l)
                              pixcol(3) = w1 * (w3 * (fetch_pixel(satch_szab, ndivuse, j, l) - (perc2-perc1))) + w2 * fetch_pixel(satch_valb, ndivuse, j, l)
                            endif
                            long_stereo = fetch_ll(satch_ll, 2, ndivuse, 1, j, l)
                            lat_stereo  = fetch_ll(satch_ll, 2, ndivuse, 2, j, l)
                            if (swaqmod) then
                              if (swaqmpan) then
                                call map_stereo_zoom_place_aqm_pan(stereo_image,stereo_count,3,nxs,nys,long_stereo, lat_stereo, pixcol,i_merge)
                              else
                                call map_stereo_zoom_place_aqm(stereo_image,stereo_count,3,nxs,nys,long_stereo, lat_stereo, pixcol)
                              endif
                            else
                              call map_stereo_zoom_place(stereo_image,3,nxs,nys,long_stereo, lat_stereo, pixcol)
                            endif
                          enddo
                        enddo
                      endif
                    endif
                  else
                    if (k.eq.1.or.k.eq.nsat_pix.or.i.eq.1.or.i.eq.nlines_gapcor) then
                    else
                      call grid_sample(satch_image(k-1,szach,i-1),satch_image(k  ,szach,i-1),satch_image(k+1,szach,i-1),
     *                                 satch_image(k-1,szach,i)  ,satch_image(k  ,szach,i)  ,satch_image(k+1,szach,i)  ,
     *                                 satch_image(k-1,szach,i+1),satch_image(k  ,szach,i+1),satch_image(k+1,szach,i+1),
     *                                 longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valr, satch_ll, ndivuse)
                      call grid_sample(satch_image(k-1,szach,i-1),satch_image(k  ,szach,i-1),satch_image(k+1,szach,i-1),
     *                                 satch_image(k-1,szach,i)  ,satch_image(k  ,szach,i)  ,satch_image(k+1,szach,i)  ,
     *                                 satch_image(k-1,szach,i+1),satch_image(k  ,szach,i+1),satch_image(k+1,szach,i+1),
     *                                 longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valg, satch_ll, ndivuse)
                      call grid_sample(satch_image(k-1,szach,i-1),satch_image(k  ,szach,i-1),satch_image(k+1,szach,i-1),
     *                                 satch_image(k-1,szach,i)  ,satch_image(k  ,szach,i)  ,satch_image(k+1,szach,i)  ,
     *                                 satch_image(k-1,szach,i+1),satch_image(k  ,szach,i+1),satch_image(k+1,szach,i+1),
     *                                 longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valb, satch_ll, ndivuse)
                      if (swok) then
                        do l = 1, ndivuse
                          do j = 1, ndivuse
                            if (swmerge.and.swaqmpan) then
                              if (i_merge.eq.3) then
                                pixcol(1) = 0
                                pixcol(2) = 0
                                pixcol(3) = (w3 * (fetch_pixel(satch_valb, ndivuse, j, l) - (perc2-perc1)))
                              else if (i_merge.eq.2) then
                                pixcol(1) = 0
                                pixcol(2) = (w3 * (fetch_pixel(satch_valg, ndivuse, j, l) - (perc2-perc1)))
                                pixcol(3) = 0
                              else if (i_merge.eq.1) then
                                pixcol(1) = (w3 * (fetch_pixel(satch_valr, ndivuse, j, l) - (perc2-perc1)))
                                pixcol(2) = 0
                                pixcol(3) = 0
                              endif
                            else
                              pixcol(1) = (w3 * (fetch_pixel(satch_valr, ndivuse, j, l) - (perc2-perc1)))
                              pixcol(2) = (w3 * (fetch_pixel(satch_valg, ndivuse, j, l) - (perc2-perc1)))
                              pixcol(3) = (w3 * (fetch_pixel(satch_valb, ndivuse, j, l) - (perc2-perc1)))
                            endif
                            long_stereo = fetch_ll(satch_ll, 2, ndivuse, 1, j, l)
                            lat_stereo  = fetch_ll(satch_ll, 2, ndivuse, 2, j, l)
                            if (swaqmod) then
                              if (swaqmpan) then
                                call map_stereo_zoom_place_aqm_pan(stereo_image,stereo_count,3,nxs,nys,long_stereo, lat_stereo, pixcol,i_merge)
                              else
                                call map_stereo_zoom_place_aqm(stereo_image,stereo_count,3,nxs,nys,long_stereo, lat_stereo, pixcol)
                              endif
                            else
                              call map_stereo_zoom_place(stereo_image,3,nxs,nys,long_stereo, lat_stereo, pixcol)
                            endif
                          enddo
                        enddo
                      endif
                    endif
                  endif
                endif
c-- Then the code for the normal imaging (with the opacity part) - Lots of code repeat here !! (need to clean up)
              else
c-- Map to IQ colours to try and separate the land
                dr = float(satch_image(k,irp,i))
                dg = float(satch_image(k,igp,i))
                db = float(satch_image(k,ibp,i))
                dy = 0.299 * dr + 0.587 * dg + 0.114 * db
                di = 0.596 * dr - 0.274 * dg - 0.322 * db
                dq = 0.211 * dr - 0.523 * dg + 0.312 * db
                if (swyiq) call map_IQ_place(dy, di, dq)
                if (k.eq.1.or.k.eq.nsat_pix.or.i.eq.1.or.i.eq.nlines_gapcor) then
                else
                  call grid_sample(satch_image(k-1,irp,i-1),satch_image(k  ,irp,i-1),satch_image(k+1,irp,i-1),
     *                             satch_image(k-1,irp,i)  ,satch_image(k  ,irp,i)  ,satch_image(k+1,irp,i)  ,
     *                             satch_image(k-1,irp,i+1),satch_image(k  ,irp,i+1),satch_image(k+1,irp,i+1),
     *                             longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valr, satch_ll, ndivuse)
                  call grid_sample(satch_image(k-1,igp,i-1),satch_image(k  ,igp,i-1),satch_image(k+1,igp,i-1),
     *                             satch_image(k-1,igp,i)  ,satch_image(k  ,igp,i)  ,satch_image(k+1,igp,i)  ,
     *                             satch_image(k-1,igp,i+1),satch_image(k  ,igp,i+1),satch_image(k+1,igp,i+1),
     *                             longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valg, satch_ll, ndivuse)
                  call grid_sample(satch_image(k-1,ibp,i-1),satch_image(k  ,ibp,i-1),satch_image(k+1,ibp,i-1),
     *                             satch_image(k-1,ibp,i)  ,satch_image(k  ,ibp,i)  ,satch_image(k+1,ibp,i)  ,
     *                             satch_image(k-1,ibp,i+1),satch_image(k  ,ibp,i+1),satch_image(k+1,ibp,i+1),
     *                             longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valb, satch_ll, ndivuse)
                  if (swopacity.and.(.not.swbatchnight)) then
                    swok = .true.
                    call grid_sample(satch_image(k-1,opchr,i-1),satch_image(k  ,opchr,i-1),satch_image(k+1,opchr,i-1),
     *                               satch_image(k-1,opchr,i)  ,satch_image(k  ,opchr,i)  ,satch_image(k+1,opchr,i)  ,
     *                               satch_image(k-1,opchr,i+1),satch_image(k  ,opchr,i+1),satch_image(k+1,opchr,i+1),
     *                               longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valr_op, satch_ll, ndivuse)
                    call grid_sample(satch_image(k-1,opchg,i-1),satch_image(k  ,opchg,i-1),satch_image(k+1,opchg,i-1),
     *                               satch_image(k-1,opchg,i)  ,satch_image(k  ,opchg,i)  ,satch_image(k+1,opchg,i)  ,
     *                               satch_image(k-1,opchg,i+1),satch_image(k  ,opchg,i+1),satch_image(k+1,opchg,i+1),
     *                               longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valg_op, satch_ll, ndivuse)
                    call grid_sample(satch_image(k-1,opchb,i-1),satch_image(k  ,opchb,i-1),satch_image(k+1,opchb,i-1),
     *                               satch_image(k-1,opchb,i)  ,satch_image(k  ,opchb,i)  ,satch_image(k+1,opchb,i)  ,
     *                               satch_image(k-1,opchb,i+1),satch_image(k  ,opchb,i+1),satch_image(k+1,opchb,i+1),
     *                               longlat(1,k,i), longlat(1,k+1,i), longlat(1,k,i+1), longlat(1,k+1,i+1), satch_valb_op, satch_ll, ndivuse)
                  endif
                  if (swok) then
                    if (ndivuse.eq.1) then
                      if (swopacity.and.(.not.swbatchnight)) then
                        pixcol(1) = wop1 * satch_image(k  ,irp,i) + wop2 * satch_image(k  , opchr,i)
                        pixcol(2) = wop1 * satch_image(k  ,igp,i) + wop2 * satch_image(k  , opchg,i)
                        pixcol(3) = wop1 * satch_image(k  ,ibp,i) + wop2 * satch_image(k  , opchb,i)
                      else
                        if (swmerge.and.swaqmpan) then
                          if (i_merge.eq.3) then
                            pixcol(1) = 0
                            pixcol(2) = 0
                            pixcol(3) = satch_image(k  , ibp, i)
                          else if (i_merge.eq.2) then
                            pixcol(1) = 0
                            pixcol(2) = satch_image(k  , igp, i)
                            pixcol(3) = 0
                          else if (i_merge.eq.1) then
                            pixcol(1) = satch_image(k  , irp, i)
                            pixcol(2) = 0
                            pixcol(3) = 0
                          endif
                        else
                          pixcol(1) = satch_image(k  , irp, i)
                          pixcol(2) = satch_image(k  , igp, i)
                          pixcol(3) = satch_image(k  , ibp, i)
                        endif
                      endif
                      long_stereo = (longlat(1,k,i) + longlat(1,k+1,i)) / 2.0D0
                      lat_stereo  = (longlat(2,k,i) + longlat(2,k+1,i)) / 2.0D0
                      if (swaqmod) then
                        if (swaqmpan) then
                          call map_stereo_zoom_place_aqm_pan(stereo_image,stereo_count,3,nxs,nys,long_stereo, lat_stereo, pixcol,i_merge)
                        else
                          call map_stereo_zoom_place_aqm(stereo_image,stereo_count,3,nxs,nys,long_stereo, lat_stereo, pixcol)
                        endif
                      else
                        call map_stereo_zoom_place(stereo_image,3,nxs,nys,long_stereo, lat_stereo, pixcol)
                      endif
                    else
                      do l = 1, ndivuse
                        do j = 1, ndivuse
                          if (swopacity.and.(.not.swbatchnight)) then
                            pixcol(1) = wop1 * fetch_pixel(satch_valr, ndivuse, j, l) + wop2 * fetch_pixel(satch_valr_op, ndivuse, j, l)
                            pixcol(2) = wop1 * fetch_pixel(satch_valg, ndivuse, j, l) + wop2 * fetch_pixel(satch_valg_op, ndivuse, j, l)
                            pixcol(3) = wop1 * fetch_pixel(satch_valb, ndivuse, j, l) + wop2 * fetch_pixel(satch_valb_op, ndivuse, j, l)
c                          pixcol(1) = wop1 * satch_valr(j,l) + wop2 * satch_valr_op(j,l)
c                          pixcol(2) = wop1 * satch_valg(j,l) + wop2 * satch_valg_op(j,l)
c                          pixcol(3) = wop1 * satch_valb(j,l) + wop2 * satch_valb_op(j,l)
                          else
                            if (swmerge.and.swaqmpan) then
                              if (i_merge.eq.3) then
                                pixcol(1) = 0
                                pixcol(2) = 0
                                pixcol(3) = fetch_pixel(satch_valb, ndivuse, j, l)
                              else if (i_merge.eq.2) then
                                pixcol(1) = 0
                                pixcol(2) = fetch_pixel(satch_valg, ndivuse, j, l)
                                pixcol(3) = 0
                              else if (i_merge.eq.1) then
                                pixcol(1) = fetch_pixel(satch_valr, ndivuse, j, l)
                                pixcol(2) = 0
                                pixcol(3) = 0
                              endif
                            else
                              pixcol(1) = fetch_pixel(satch_valr, ndivuse, j, l)
                              pixcol(2) = fetch_pixel(satch_valg, ndivuse, j, l)
                              pixcol(3) = fetch_pixel(satch_valb, ndivuse, j, l)
                            endif
                          endif
                          long_stereo = fetch_ll(satch_ll, 2, ndivuse, 1, j, l)
                          lat_stereo  = fetch_ll(satch_ll, 2, ndivuse, 2, j, l)
                          if (swaqmod) then
                            if (swaqmpan) then
                              call map_stereo_zoom_place_aqm_pan(stereo_image,stereo_count,3,nxs,nys,long_stereo, lat_stereo, pixcol,i_merge)
                            else
                              call map_stereo_zoom_place_aqm(stereo_image,stereo_count,3,nxs,nys,long_stereo, lat_stereo, pixcol)
                            endif
                          else
                            call map_stereo_zoom_place(stereo_image,3,nxs,nys,long_stereo, lat_stereo, pixcol)
                          endif
                        enddo
                      enddo
                    endif
                  endif
                endif
              endif
            enddo
          endif
        enddo
c-- Reload the merge file for 1 to nmerge - 1 then jump back up to mark them on the stereo projection
        if (swlastmerge.and.i_merge.gt.1) then
          i_merge = i_merge - 1
          nsat_pix = nsat_pix_merge(i_merge)
          nlines_gapcor = nlines_gapcor_merge(i_merge)
          swnorth = swnorth_merge(i_merge)
          call load_merge(i_merge, satch_image, satch_valid, longlat, nmaxpix, nmaxch, nlines_gapcor, lundump, nmergemax)
c-- If the user wanted szamerge on the last pass to be merged now disable it to prevent ruining the other to be merged passes
          swszamer = .false.
          goto 302
        endif
        if (swlastmerge.and.i_merge.eq.1) then
          nsat_pix = nsat_pix_merge(nmerge)
          nlines_gapcor = nlines_gapcor_merge(nmerge)
          swnorth = swnorth_merge(nmerge)
        endif
c-- Average the projected zoom based on pixel hit count if necessary. Thereafter goto border loading and projection writing as necessary
        if (swaqmod) then
          if (swaqmpan) then
            call map_stereo_zoom_aver_aqm_pan(stereo_image,stereo_count,3,nxs,nys)
          else
            call map_stereo_zoom_aver_aqm(stereo_image,stereo_count,3,nxs,nys)
          endif
        endif
        if (swborder) then
          if (swborderhighres) then
            call get_borders(stereo_image,3,nxs,nys,swborderhighres, 1)
            call get_borders(stereo_image,3,nxs,nys,swborderhighres, 2)
            call get_borders(stereo_image,3,nxs,nys,swborderhighres, 3)
            call get_borders(stereo_image,3,nxs,nys,swborderhighres, 4)
          else
            call get_borders(stereo_image,3,nxs,nys,swborderhighres, 1)
          endif
c-- If required, draw the longlat grid after the borders
          if (swlonglat) then
            call map_stereo_zoom_mark_longlat(stereo_image,3,nxs,nys,longstep,latstep)
          endif
          call get_lun(lunborder)
          if (swfy)    k = index(filedata,cfy)
          if (swnoaa)  k = index(filedata,cnoaa)
          if (swmetop) k = index(filedata,cmetop)
          if (swmn2  ) k = index(filedata,cmn2)
          if (swaqmod) k = index(filedata,caqmod)
          if (swthetascan) then
            if (thetascan.lt.0.0) cscan(1:1) = 'M'
            if (thetascan.ge.0.0) cscan(1:1) = 'P'
            write (cscan(2:4),'(F3.1)') abs(thetascan)
            write (imgname,'(a,a)') filedata(1:k-1),'_'//cscan//'_border.rgb'
            open (unit=lunborder,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=3*nxs*2)
            close (unit=lunborder,status='delete', iostat=ios)
            open (unit=lunborder,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=3*nxs*2)
          else
            write (imgname,'(a,a)') filedata(1:k-1),'_border.rgb'
            open (unit=lunborder,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=3*nxs*2)
            close (unit=lunborder,status='delete', iostat=ios)
            open (unit=lunborder,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=3*nxs*2)
          endif
c-- If requested, carry out the gamma correction for separate RGB's
          if (swgammargb) then
            call gamma_correct(stereo_image, 3, nxs, nys, rgamma, ggamma, bgamma, gammalut, ncolmax, iclim, rcorrect, gcorrect, bcorrect)
          endif
c-- Write out rhe projected rgb file (with or without border !)
          do i = 1, nys
            call write_buffer(lunborder,i,stereo_image(1,i_off),3*nxs)
            i_off = i_off + nxs
          enddo
          close (unit=lunborder)
          call free_lun(lunborder)
c-- This projection always requires a flip !
          if (swthetascan) then
            if (thetascan.lt.0.0) cscan(1:1) = 'M'
            if (thetascan.ge.0.0) cscan(1:1) = 'P'
            write (cscan(2:4),'(F3.1)') abs(thetascan)
            if (swgamma) then
              if (swfy)    call create_png(filedata,cfy   ,'.rgb','_'//cscan//'_border.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swnoaa)  call create_png(filedata,cnoaa ,'.rgb','_'//cscan//'_border.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swmetop) call create_png(filedata,cmetop,'.rgb','_'//cscan//'_border.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swmn2)   call create_png(filedata,cmn2  ,'.rgb','_'//cscan//'_border.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swaqmod) call create_png(filedata,caqmod,'.rgb','_'//cscan//'_border.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
            else
              if (swfy)    call create_png(filedata,cfy   ,'.rgb','_'//cscan//'_border.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swnoaa)  call create_png(filedata,cnoaa ,'.rgb','_'//cscan//'_border.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swmetop) call create_png(filedata,cmetop,'.rgb','_'//cscan//'_border.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swmn2)   call create_png(filedata,cmn2  ,'.rgb','_'//cscan//'_border.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swaqmod) call create_png(filedata,caqmod,'.rgb','_'//cscan//'_border.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
            endif
          else
            if (swgamma) then
              if (swfy)    call create_png(filedata,cfy   ,'.rgb','_border.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swnoaa)  call create_png(filedata,cnoaa ,'.rgb','_border.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swmetop) call create_png(filedata,cmetop,'.rgb','_border.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swmn2)   call create_png(filedata,cmn2  ,'.rgb','_border.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swaqmod) call create_png(filedata,caqmod,'.rgb','_border.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
            else
              if (swfy)    call create_png(filedata,cfy   ,'.rgb','_border.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swnoaa)  call create_png(filedata,cnoaa ,'.rgb','_border.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swmetop) call create_png(filedata,cmetop,'.rgb','_border.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swmn2)   call create_png(filedata,cmn2  ,'.rgb','_border.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
              if (swaqmod) call create_png(filedata,caqmod,'.rgb','_border.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
            endif
          endif
        else
c-- Draw the longlat grid
          if (swlonglat) then
            call map_stereo_zoom_mark_longlat(stereo_image,3,nxs,nys,longstep,latstep)
          endif
          call get_lun(lunborder)
          if (swfy)    k = index(filedata,cfy)
          if (swnoaa)  k = index(filedata,cnoaa)
          if (swmetop) k = index(filedata,cmetop)
          if (swmn2)   k = index(filedata,cmn2)
          if (swaqmod) k = index(filedata,caqmod)
          write (imgname,'(a,a)') filedata(1:k-1),'_project.rgb'
          open (unit=lunborder,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=3*nxs*2)
          close (unit=lunborder,status='delete', iostat=ios)
          open (unit=lunborder,file=imgname(1:lnblnk(imgname)),form='unformatted', access='direct',recl=3*nxs*2)
c-- If requested, carry out the gamma correction for separate RGB's
          if (swgammargb) then
            call gamma_correct(stereo_image, 3, nxs, nys, rgamma, ggamma, bgamma, gammalut, ncolmax, iclim, rcorrect, gcorrect, bcorrect)
          endif
c-- Write out rhe projected rgb file (with or without border !)
          do i = 1, nys
            call write_buffer(lunborder,i,stereo_image(1,i_off),3*nxs)
            i_off = i_off + nxs
          enddo
          close (unit=lunborder)
          call free_lun(lunborder)
c-- This projection always requires a flip !
          if (swgamma) then
            if (swfy)    call create_png(filedata,cfy   ,'.rgb','_project.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
            if (swnoaa)  call create_png(filedata,cnoaa ,'.rgb','_project.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
            if (swmetop) call create_png(filedata,cmetop,'.rgb','_project.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
            if (swmn2)   call create_png(filedata,cmn2  ,'.rgb','_project.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
            if (swaqmod) call create_png(filedata,caqmod,'.rgb','_project.rgb',0, '-'//cgamma, .true., nxs, -nys, swbatch, swsharp, csharp)
          else
            if (swfy)    call create_png(filedata,cfy   ,'.rgb','_project.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
            if (swnoaa)  call create_png(filedata,cnoaa ,'.rgb','_project.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
            if (swmetop) call create_png(filedata,cmetop,'.rgb','_project.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
            if (swmn2)   call create_png(filedata,cmn2  ,'.rgb','_project.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
            if (swaqmod) call create_png(filedata,caqmod,'.rgb','_project.rgb',0,' ', .true., nxs, -nys, swbatch, swsharp, csharp)
          endif
        endif
      endif
c
c-- End of the projection/border part
c
 300  continue
      if (.not.swfast) then
        do i = 1,nchan
          nwrite = irec
          if (swgapcor) nwrite = nlines_gapcor
          if (swfy)    call create_png(filedata,cfy   ,'.dat','_ch',i,' ', swnorth, nsat_pix, nwrite, swbatch, swsharp, csharp)
          if (swnoaa)  call create_png(filedata,cnoaa ,'.dat','_ch',i,' ', swnorth, nsat_pix, nwrite, swbatch, swsharp, csharp)
          if (swmetop) call create_png(filedata,cmetop,'.dat','_ch',i,' ', swnorth, nsat_pix, nwrite, swbatch, swsharp, csharp)
          if (swmn2)   call create_png(filedata,cmn2  ,'.dat','_ch',i,' ', swnorth, nsat_pix, nwrite, swbatch, swsharp, csharp)
          if (swaqmod) call create_png(filedata,caqmod,'.dat','_ch',i,' ', swnorth, nsat_pix, nwrite, swbatch, swsharp, csharp)
        enddo
      endif
      if (swhybrid) then
        do i = 1,nchan
          nwrite = irec
          if (swgapcor) nwrite = nlines_gapcor
          if (swfy)    call create_png(filedata,cfy   ,'.rgb','_ch1',i,' ', swnorth, nsat_pix, nwrite, swbatch, swsharp, csharp)
          if (swnoaa)  call create_png(filedata,cnoaa ,'.rgb','_ch1',i,' ', swnorth, nsat_pix, nwrite, swbatch, swsharp, csharp)
          if (swmetop) call create_png(filedata,cmetop,'.rgb','_ch1',i,' ', swnorth, nsat_pix, nwrite, swbatch, swsharp, csharp)
          if (swmn2)   call create_png(filedata,cmn2  ,'.rgb','_ch1',i,' ', swnorth, nsat_pix, nwrite, swbatch, swsharp, csharp)
          if (swaqmod) call create_png(filedata,caqmod,'.rgb','_ch1',i,' ', swnorth, nsat_pix, nwrite, swbatch, swsharp, csharp)
        enddo
      endif
      if (swdebug) then
        do i = 1,iclim+1
          write (*,'(11I7)') i, (hist(i,k), k=1,nchan)
        enddo
      endif
      if (swfy)    k = index(filedata,cfy)
      if (swnoaa)  k = index(filedata,cnoaa)
      if (swmetop) k = index(filedata,cmetop)
      if (swmn2)   k = index(filedata,cmn2)
      if (swaqmod) k = index(filedata,caqmod)
      write (imgname,'(a,a)') filedata(1:k-1),'_ch'
      j = index(imgname,'_ch')
      if (.not.swfast) call hist_plot(hist, 4096, 10, iclim+1, nchan,imgname(1:j+2))
      write (imgname,'(a,a)') filedata(1:k-1),'_sza'
      j = index(imgname,'_sza')
      if (swszalim.and.(.not.swbatch)) then
        call sza_plot(szahist, n_sza, imgname(1:j+3))
      endif
c
      if (swfy) then
        nwrite = irec
        if (swgapcor) nwrite = nlines_gapcor
        if (swgamma) then
          call create_png(filedata,cfy,'.rgb','_true.rgb'        ,0, '-'//cgamma, swnorth, nsat_pix      , nwrite, swbatch, swsharp, csharp)
          call create_png(filedata,cfy,'.rgb','_true_correct.rgb',0, '-'//cgamma, swnorth, 2*ncorrect_pix, nwrite, swbatch, swsharp, csharp)
        else
          call create_png(filedata,cfy,'.rgb','_true.rgb'        ,0,'-gamma 1.4', swnorth, nsat_pix      , nwrite, swbatch, swsharp, csharp)
          call create_png(filedata,cfy,'.rgb','_true_correct.rgb',0,'-gamma 1.4', swnorth, 2*ncorrect_pix, nwrite, swbatch, swsharp, csharp)
        endif
      endif
c-- Write the CH11 corrected image
      nwrite = irec
      if (swgapcor) nwrite = nlines_gapcor
      if (swgamma) then
        if (swfy)    call create_png(filedata,cfy   ,'.rgb','_ch11_correct.rgb',0, '-'//cgamma, swnorth, 2*ncorrect_pix, nwrite, swbatch, swsharp, csharp)
        if (swnoaa)  call create_png(filedata,cnoaa ,'.rgb','_ch11_correct.rgb',0, '-'//cgamma, swnorth, 2*ncorrect_pix, nwrite, swbatch, swsharp, csharp)
        if (swmetop) call create_png(filedata,cmetop,'.rgb','_ch11_correct.rgb',0, '-'//cgamma, swnorth, 2*ncorrect_pix, nwrite, swbatch, swsharp, csharp)
        if (swmn2)   call create_png(filedata,cmn2  ,'.rgb','_ch11_correct.rgb',0, '-'//cgamma, swnorth, 2*ncorrect_pix, nwrite, swbatch, swsharp, csharp)
        if (swaqmod) call create_png(filedata,caqmod,'.rgb','_ch11_correct.rgb',0, '-'//cgamma, swnorth, 2*ncorrect_pix, nwrite, swbatch, swsharp, csharp)
      else
        if (swfy)    call create_png(filedata,cfy   ,'.rgb','_ch11_correct.rgb',0,'-gamma 1.4', swnorth, 2*ncorrect_pix, nwrite, swbatch, swsharp, csharp)
        if (swnoaa)  call create_png(filedata,cnoaa ,'.rgb','_ch11_correct.rgb',0,'-gamma 1.4', swnorth, 2*ncorrect_pix, nwrite, swbatch, swsharp, csharp)
        if (swmetop) call create_png(filedata,cmetop,'.rgb','_ch11_correct.rgb',0,'-gamma 1.4', swnorth, 2*ncorrect_pix, nwrite, swbatch, swsharp, csharp)
        if (swmn2)   call create_png(filedata,cmn2  ,'.rgb','_ch11_correct.rgb',0,'-gamma 1.4', swnorth, 2*ncorrect_pix, nwrite, swbatch, swsharp, csharp)
        if (swaqmod) call create_png(filedata,caqmod,'.rgb','_ch11_correct.rgb',0,'-gamma 1.4', swnorth, 2*ncorrect_pix, nwrite, swbatch, swsharp, csharp)
      endif
c
c
      if (swmerge.and.(.not.swlastmerge)) goto 301
      call map_write()
      if (swgapcor) call map_stereo_write(command)
      if (swgapcor) call map_azel_write(command)
      if (swgapcor.and.swbatch.and.swhorizon) call map_azel_horizon_write(command)
      if (swyiq) call map_IQ_write(command)
      if (swfy)    k = index(filedata,cfy)
      if (swnoaa)  k = index(filedata,cnoaa)
      if (swmetop) k = index(filedata,cmetop)
      if (swmn2)   k = index(filedata,cmn2)
      if (swaqmod) k = index(filedata,caqmod)
      write (imgname,'(a,a)') filedata(1:k-1),'_satmap.png'
      call system('mv satmap.png '//imgname(1:lnblnk(imgname)))
      write (*,'(a)') ' mv satmap.png '//imgname(1:lnblnk(imgname))
      if (swgapcor) then
        if (.not.swfast) then
          write (imgname,'(a,a)') filedata(1:k-1),'_satmap-stereo.png'
          write (*,'(a)') ' mv satmap-stereo.png '//imgname(1:lnblnk(imgname))
          call system('mv satmap-stereo.png '//imgname(1:lnblnk(imgname)))
        endif
        write (imgname,'(a,a)') filedata(1:k-1),'_satmap-stereo-zoom.png'
        write (*,'(a)') ' mv satmap-stereo-zoom.png '//imgname(1:lnblnk(imgname))
        call system('mv satmap-stereo-zoom.png '//imgname(1:lnblnk(imgname)))
        write (imgname,'(a,a)') filedata(1:k-1),'_azel.png'
        write (*,'(a)') ' mv satmap-azel.png '//imgname(1:lnblnk(imgname))
        call system('mv satmap-azel.png '//imgname(1:lnblnk(imgname)))
        if (swbatch.and.swgapcor.and.swhorizon) then
          write (imgname,'(a,a)') filedata(1:k-1),'_azelh.png'
          write (*,'(a)') ' mv satmap-azelh.png '//imgname(1:lnblnk(imgname))
          call system('mv satmap-azelh.png '//imgname(1:lnblnk(imgname)))
        endif
        if (swyiq) then
          write (imgname,'(a,a)') filedata(1:k-1),'_IQ.png'
          write (*,'(a)') ' mv satmap-IQ.png '//imgname(1:lnblnk(imgname))
          call system('mv satmap-IQ.png '//imgname(1:lnblnk(imgname)))
        endif
      endif
 301  continue
c
      nvalid = 0
      do i = 1, irec
        if (linevalid(i).eq.1) nvalid = nvalid + 1
      enddo
      write (*,'(''# of records : '',I8,'' of which : '',I8,'' in satellite data file and : '',I8,'' with a valid (in sequence) timestamp '')') nwrite, irec, nvalid
c
      call sgp4_end()
      if (swthetascan) goto 100
 101  continue
      if (swbatch) then
        if (swmerge) then
          if (swlastmerge) goto 201
          i_merge = i_merge + 1
          if (swfirstmerge) then
c-- Fix the zoom ration based on the first image in the sequence => is a test
            call map_stereo_rzoom(nsat_pix,irec)
            swfirstmerge = .false.
          endif
          if (i_merge.eq.nmerge) swlastmerge = .true.
        endif
        call get_scanline_reset()
        goto 200
      endif
 201  continue
      close (unit=lunbatch)
      call free_lun(lunbatch)
      stop
      end

      subroutine correct_init(sat_avg_alt, sat_fov, nsat_pix, index_correct, nmax, ncorrect_pix, swaqmod)
      implicit none
      integer*4 nsat_pix, ncorrect_pix, i, nmax
      integer*2 index_correct(nmax)
      real*4 sat_avg_alt, sat_fov
      logical swaqmod
c
      real*4 pix_ang, pix_siz, pix_pos, lat_off, sca_ang, rad_earth, cir_earth, pi
c
      save
c
      pi           = 2.0 * asin(1.0)
      ncorrect_pix = 0
      rad_earth    = 6371.0
      cir_earth    = pi * 2.0 * rad_earth
      pix_ang      = pi * (sat_fov / 180.0) / float(nsat_pix)
      pix_siz      = sat_avg_alt * tan(pix_ang)
c      write (*,*) pix_ang, pix_siz, (atan(pix_siz/rad_earth)/(2.0*pi))*cir_earth
c-- overwrite the calculated pixel size, as aparently (too much scan angle correction) it is wrong
c-- strange enough the pix_ang appears to be correct ..........
      if (.not.swaqmod) pix_siz = 1.08
      i = 1
      do while (.true.)
        pix_pos = 2.0 * pi * (((pix_siz/2.0) + (float(i-1) * pix_siz)) / cir_earth)
        lat_off = sin(pix_pos) * rad_earth
        sca_ang = atan(lat_off/(sat_avg_alt+rad_earth*(1.0-cos(pix_pos))))
c        write (*,*) i, nint(sca_ang/pix_ang)
        ncorrect_pix = ncorrect_pix + 1
        index_correct(i) = nint(sca_ang/pix_ang)
        if (nint(sca_ang/pix_ang).eq.nsat_pix/2) goto 1
        i = i + 1
      enddo
 1    continue
c      read (*,*)
c
      return
      end

      subroutine correct_apply(rgbbuf, nx, ny, rgbbuf_correct, nxc, nyc, index_correct, nxi, lun, irec)
      implicit none
      integer*4 nx, ny, nxc, nyc, nxi, lun, irec
      integer*2 rgbbuf(nx, ny), rgbbuf_correct(nxc, nyc), index_correct(nxi)
c
      integer*4 i, j, k, k1, k2
      if (nx.ne.3.or.nxc.ne.3) stop '** correct_apply - nx and nyc have to have a value of 3 (RGB only)'
      j = nyc / 2
      k = ny  / 2
      do i = 1,j
        k1 = index_correct(i) + k
        k2 = k - index_correct(i) + 1
        rgbbuf_correct(1, j + i)     = rgbbuf(1,k1)
        rgbbuf_correct(1, j - i + 1) = rgbbuf(1,k2)
        rgbbuf_correct(2, j + i)     = rgbbuf(2,k1)
        rgbbuf_correct(2, j - i + 1) = rgbbuf(2,k2)
        rgbbuf_correct(3, j + i)     = rgbbuf(3,k1)
        rgbbuf_correct(3, j - i + 1) = rgbbuf(3,k2)
      enddo
      call correct_write (lun, irec, rgbbuf_correct, nxc, nyc)
c
      return
      end

      subroutine correct_write (lun_correct, irec, histrgb_correct, nx, ny)
      implicit none
      integer*4 nx, ny, irec, lun_correct
      integer*2 histrgb_correct(nx, ny)
c
      write (lun_correct,rec=irec) histrgb_correct
c
      return
      end

      subroutine get_scanline(lun, irec, ios, channel_image, nmaxpix, nchan, npix, cmis, buf, nbuf, timestamp, doy, swdebug, swbin, aqmtheta)
      implicit none
      real*8    timestamp, aqmtheta(2)
      integer*4 lun, irec, ios, nmaxpix, npix, nchan, nbuf, doy
      integer*2 channel_image(nmaxpix, 10)
      integer*1 buf(nbuf)
      logical   swdebug, swbin
      character*2 cmis
c
      integer*1 buftim(28)
      integer*4 i, j, k, l, k1, k2, ihr, imn, init, i_time, iyr, imo, idy, dy2k, amdoy, lunam(10), luntim, initam, aqmnpix, nchan_local
      real*4    rsc
      real*8    offset, amtim, amtheta1, amtheta2
      character*2 cam1000(5), cam500(5), cam250(2)
      character*250 filename
c
      logical       swaqm1000, swaqm500, swaqm250, swaqmpan
      common /MODIS/swaqm1000, swaqm500, swaqm250, swaqmpan
c
      data init/0/, initam/0/, cam1000/'08','09','20','21','22'/, cam500/'03','04','05','06','07'/, cam250/'01','02'/
      equivalence (buftim(1), amtim), (buftim(9), amdoy), (buftim(13), amtheta1), (buftim(21), amtheta2)
c
      save
c
      aqmtheta(1) = 0.0D0
      aqmtheta(2) = 0.0D0
      if (initam.eq.0.and.cmis.eq.'AM') then
        call get_string(filename)
        j = index(filename,'xx')
        if (swaqm1000) then
          aqmnpix = 1354
          do i = 1,nchan
            call get_lun(lunam(i))
            filename(j:j+1) = cam1000(i)
            open (unit=lunam(i),file=filename(1:lnblnk(filename)),form='unformatted',access='direct',recl=2*aqmnpix)
          enddo
          call get_lun(luntim)
          open (unit=luntim,file=filename(1:j-1)//cam1000(1)//'.tim',form='unformatted',access='direct',recl=28)
        endif
        if (swaqm500) then
          aqmnpix = 2708
          do i = 1,nchan
            call get_lun(lunam(i))
            filename(j:j+1) = cam500(i)
            open (unit=lunam(i),file=filename(1:lnblnk(filename)),form='unformatted',access='direct',recl=2*aqmnpix)
          enddo
          call get_lun(luntim)
          open (unit=luntim,file=filename(1:j-1)//cam500(1)//'.tim',form='unformatted',access='direct',recl=28)
        endif
        if (swaqm250) then
          aqmnpix = 5416
          do i = 1,nchan
            call get_lun(lunam(i))
            filename(j:j+1) = cam250(i)
            open (unit=lunam(i),file=filename(1:lnblnk(filename)),form='unformatted',access='direct',recl=2*aqmnpix)
          enddo
          call get_lun(luntim)
          open (unit=luntim,file=filename(1:j-1)//cam250(1)//'.tim',form='unformatted',access='direct',recl=28)
        endif
        nchan_local = nchan
        initam = 1
      endif
      if (cmis.eq.'FY') then
        read (lun,rec=irec,iostat=ios) buf
        if (ios.ne.0) then
          init = 1
          goto 1
        endif
        k = 2001
        if (swbin) k = 437
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
        i_time = 11
        if (swbin) i_time = 26044
        call HRPT_time_FyMe(buf(i_time + 1),4,timestamp)
        if (timestamp.lt.43200.0D0) offset =  43200.0D0
        if (timestamp.ge.43200.0D0) offset = -43200.0D0
        timestamp = timestamp + offset
c-- Cowards approach for the time being. Cannot find the doy in the bit-shifted .bin files (is not before the time), so, use the date !
        if (swbin) then
          inquire (unit=lun,name=filename)
          j = index(filename,'FY3')
          read (filename(j+5 :j+8 ),'(I4)') iyr
          read (filename(j+10:j+11),'(I2)') imo
          read (filename(j+13:j+14),'(I2)') idy
          call datetodoy(iyr, imo, idy, doy)
        else
          doy = buf(i_time)
          if (doy.le.0) doy = doy + 256
          doy = doy * 2 + ishft(iand(buf(i_time + 1),'80'x),-7)
        endif
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
        if (swbin) k = 15
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
        i_time = 11
        if (swbin) i_time = 2
        call HRPT_time_FyMe(buf(i_time + 1),4,timestamp)
        if (swbin) then
          dy2k = buf(i_time)
          if (dy2k.le.0) dy2k = dy2k + 256
          dy2k = buf(i_time-1) * 256 + dy2k
          call dy2k_to_doy(dy2k,doy)
        else
          doy = buf(i_time)
          if (doy.le.0) doy = doy + 256
          doy = doy * 2 + ishft(iand(buf(i_time + 1),'80'x),-7)
        endif
        ihr = int (timestamp / 3600.0D0)
        imn = int((timestamp - dble(ihr)*3600.0D0)/60.0D0)
        rsc = timestamp - dble(ihr) * 3600.0D0 - dble(imn) * 60.0D0
        if (init.eq.0.and.swdebug) write (*,'(I4,2X,24(Z2.2,2x),F20.4,3I5,F8.3)') irec, (buf(j), buf(j+1), j = 1, 23, 2), timestamp, doy, ihr, imn, rsc
      endif
      if (cmis.eq.'M2') then
        read (lun,rec=irec,iostat=ios) buf
        if (ios.ne.0) then
          init = 1
          goto 4
        endif
        k = 0
        if (swbin) k = 51
        do i = 1, 1569, 4
          do j = 1,6
            k1 = buf(k)
            if (k1.lt.0) k1 = k1 + 256
            k1 = k1 * 4
            k2 = ishft(iand(buf(k+1),'C0'x),-6)
            channel_image(i  ,j) = k1 + k2
            k1 = iand(buf(k+1),'3F'x)
            k1 = k1 * 16
            k2 = ishft(iand(buf(k+2),'F0'x),-4)
            channel_image(i+1,j) = k1 + k2
            k1 = iand(buf(k+2),'0F'x)
            k1 = k1 * 64
            k2 = ishft(iand(buf(k+3),'FC'x),-2)
            channel_image(i+2,j) = k1 + k2
            k1 = iand(buf(k+3),'03'x)
            k1 = k1 * 256
            k2 = buf(k+4)
            if (k2.lt.0) k2 = k2 + 256
            channel_image(i+3,j) = k1 + k2
            k = k + 5
          enddo
        enddo
 4      continue
        if (swbin) then
          inquire (unit=lun,name=filename)
          if (index(filename,'MN2-2').ne.0.or.index(filename,'MN2_2').ne.0) then
            j = index(filename,'MN2') + 2
          else
            j = index(filename,'MN2')
          endif
          read (filename(j+9 :j+12 ),'(I4)') iyr
          read (filename(j+13:j+14),'(I2)')  imo
          read (filename(j+15:j+16),'(I2)')  idy
          call datetodoy(iyr, imo, idy, doy)
        else
        endif
c-- Meteor uses UTC + 3 - so correct here, but nightly runs before 03:00:00 will currently fail ! lower doy !
        ihr = buf( 9) - 3
        imn = buf(10)
        rsc = float(buf(12))
        if (rsc.lt.0.0) rsc = rsc + 256.0
c-- Sub-second interval is just a wild guess - replaced original guess by sum value derived from dumps
c        rsc = rsc * 0.00438657 + float(buf(11))
        rsc = rsc * 0.0043 + float(buf(11))
        timestamp = dble(ihr) * 3600.0D0 + dble(imn) * 60.0D0 + dble(rsc)
        if (init.eq.0.and.swdebug) write (*,'(I4,2X,24(Z2.2,2x),F20.4,3I5,F8.3)') irec, (buf(j), buf(j+1), j = 1, 23, 2), timestamp, doy, ihr, imn, rsc
      endif
      if (cmis.eq.'AM') then
        do j = 1, nchan
          read (lunam(j),rec=irec,iostat=ios) buf
          if (ios.ne.0) then
            init = 1
            goto 5
          endif
          k = 1
          do i = 1, aqmnpix
            k1 = buf(k)
            if (k1.lt.0) k1 = k1 + 256
            k2 = buf(k+1)
            if (k2.lt.0) k2 = k2 + 256
            channel_image(i  ,j) = (k1 + k2 * 256)
            k = k + 2
          enddo
        enddo
        read (luntim,rec=irec) buftim
c-- I do not understand the - 1.477815D0 , but apparently the MODIS time stamps label the last sample in a packt .... or something (20190929)
        timestamp   = amtim - 1.477815D0
        doy         = amdoy
        aqmtheta(1) = amtheta1
        aqmtheta(2) = amtheta2
        ihr         = int (timestamp / 3600.0D0)
        imn         = int((timestamp - dble(ihr)*3600.0D0)/60.0D0)
        rsc         = timestamp - dble(ihr) * 3600.0D0 - dble(imn) * 60.0D0
 5      continue
        if (init.eq.0.and.swdebug) write (*,'(I5,2X,24(Z2.2,2x),F20.4,3I5,3F8.3)') irec, (buf(j), buf(j+1), j = 1, 23, 2), timestamp, doy, ihr, imn, rsc, aqmtheta(1), aqmtheta(2)
      endif
c
      return
c
      entry get_scanline_reset()
      if (initam.eq.1) then
        do i = 1, nchan_local
          close (unit=lunam(i))
          call free_lun(lunam(i))
        enddo
        close (unit=luntim)
        call free_lun(luntim)
        initam = 0
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

      subroutine check_linetimes(linetime, linedoy, linevalid, irec, timstart, swdebug, swmn2, swaqmod)
c-- Determine the minimum and maximum valid time and check 'validity' ..... difficult algorithm and NOT bullet proof
c-- For Meteor satellites do a more relaxed interline filtering - appears to be not very precise - see the multiplier in get_scanline
      implicit none
      integer*4 irec
      integer*1 linevalid(irec)
      logical   swdebug, swmn2, swaqmod
      integer*4 linedoy(irec)
      real*8    linetime(irec), timstart, delta
c
      integer*4  i, j, timhist(864), i_tmax, v_tmax, i_min, i_max, ilinelast, doyhist(365), i_valid, irecval(3)
      real*8     delta_1, delta_2, timeval(3), lastvalid
c
      logical       swaqm1000, swaqm500, swaqm250, swaqmpan
      common /MODIS/swaqm1000, swaqm500, swaqm250, swaqmpan
c
      delta = 0.1666666D0
c-- For aqua we need a fix
      if (swaqm1000) delta = 0.1477815
      if (swaqm500)  delta = 0.1477815 / 2.0
      if (swaqm250)  delta = 0.1477815 / 4.0
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
          delta_1 = (linetime(i)   - linetime(i-1)) / delta
          delta_2 = (linetime(i+1) - linetime(i))   / delta
          if (.not.swmn2) then
            if (dabs(delta_1-1.0D0).lt.5.0D-3.and.dabs(delta_2-1.0D0).lt.5.0D-3) then
              linevalid(i-1) = 1
              linevalid(i)   = 1
              linevalid(i+1) = 1
            endif
          else
            if (dabs(delta_1-1.0D0).lt.5.0D-2.and.dabs(delta_2-1.0D0).lt.5.0D-2) then
              linevalid(i-1) = 1
              linevalid(i)   = 1
              linevalid(i+1) = 1
            endif
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
          if (.not.swmn2) then
            if (i.lt.i_valid) then
              delta_1 = (linetime(i_valid) - linetime(i)) / delta
              if (delta_1.gt.0.and.abs(delta_1-dble(nint(delta_1))).le.5.0D-2) then
                linevalid(i) = 1
              endif
            else
              delta_1 = (linetime(i) - linetime(i_valid)) / delta
              if (delta_1.gt.0.and.abs(delta_1-dble(nint(delta_1))).le.5.0D-2) then
                linevalid(i) = 1
              endif
            endif
          else
            if (i.lt.i_valid) then
              delta_1 = (linetime(i_valid) - linetime(i)) / delta
              if (delta_1.gt.0.and.abs(delta_1-dble(nint(delta_1))).le.8.0D-2) then
                linevalid(i) = 1
              endif
            else
              delta_1 = (linetime(i) - linetime(i_valid)) / delta
              if (delta_1.gt.0.and.abs(delta_1-dble(nint(delta_1))).le.8.0D-2) then
                linevalid(i) = 1
              endif
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
      integer*4 lun, iyr_tle, i, iyr_file, imn_file, ida_file, diy
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
        read (file(i:i+3),'(I4)') iyr_file
        read (file(i+4:i+5),'(I2)') imn_file
        read (file(i+6:i+7),'(I2)') ida_file
        call datetodoy(iyr_file, imn_file, ida_file, doyfile)
      endif
      if (csat.eq.'ME') then
        if (index(file,'_M02').ne.0) then
          i = index(file,'_M02') - 15
          read (file(i:i+3),'(I4)') iyr_file
        else
          i = index(file,'Metop') + 7
          read (file(i:i+3),'(I4)') iyr_file
        endif
      endif
      if (csat.eq.'M2') then
        if (index(file,'MN2-2').ne.0) then
          i = index(file,'MN2') + 11
        else
          i = index(file,'MN2') + 9
        endif
        read (file(i:i+3),'(I4)') iyr_file
      endif
      if (csat.eq.'AM') then
        i = index(file,'AQ_MODIS') + 9
        read (file(i:i+3),'(I4)') iyr_file
      endif
      if (csat.eq.'TE') then
        i = index(file,'TE_MODIS') + 9
        read (file(i:i+3),'(I4)') iyr_file
      endif
      write (*,'(''TLE Seek year      : '',I5)') iyr_file
      call doytodate(linedoy, iyr_file, cdate)
c      write (*,*) linedoy, cdate
      if (csat.ne.'AM'.and.csat.ne.'TE') then
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
      else
        foundtle = .false.
        call exist(tledir(1:lnblnk(tledir))//cdate//'-resource.tle', foundtle)
        if (foundtle) then
          call system('cp '//tledir(1:lnblnk(tledir))//cdate//'-resource.tle'//' ./weather.txt')
          write (*,'(''Using TLE file     : '',A)') tledir(1:lnblnk(tledir))//cdate//'-resource.tle'
        else
          call doytodate(linedoy-1, iyr_file, cdate)
          foundtle = .false.
          call exist(tledir(1:lnblnk(tledir))//cdate//'-resource.tle', foundtle)
          if (foundtle) then
            call system('cp '//tledir(1:lnblnk(tledir))//cdate//'-resource.tle'//' ./weather.txt')
            write (*,'(''Using TLE file     : '',A)') tledir(1:lnblnk(tledir))//cdate//'-resource.tle'
          else
            call doytodate(linedoy+1, iyr_file, cdate)
            foundtle = .false.
            call exist(tledir(1:lnblnk(tledir))//cdate//'-resource.tle', foundtle)
            if (foundtle) then
              call system('cp '//tledir(1:lnblnk(tledir))//cdate//'-resource.tle'//' ./weather.txt')
              write (*,'(''Using TLE file     : '',A)') tledir(1:lnblnk(tledir))//cdate//'-resource.tle'
            else
              write (*,*) 'TLE not found - using current weather.txt'
            endif
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
        if (file(index(file,'NOAA')+4:index(file,'NOAA')+4).eq.'-') then
          if (index(file,'NOAA-15').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 15" - > hrpt.tmp')
          if (index(file,'NOAA-18').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 18" - > hrpt.tmp')
          if (index(file,'NOAA-19').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 19" - > hrpt.tmp')
        else if (file(index(file,'NOAA')+4:index(file,'NOAA')+4).eq.'_') then
          if (index(file,'NOAA_15').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 15" - > hrpt.tmp')
          if (index(file,'NOAA_18').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 18" - > hrpt.tmp')
          if (index(file,'NOAA_19').ne.0) call system('cat weather.txt | grep -A 2 -e "NOAA 19" - > hrpt.tmp')
        endif
      endif
      if (csat.eq.'ME') then
        if (index(file,'opA').ne.0.or.index(file,'_M02').ne.0) call system('cat weather.txt | grep -A 2 -e "METOP-A" - > hrpt.tmp')
        if (index(file,'opB').ne.0) call system('cat weather.txt | grep -A 2 -e "METOP-B" - > hrpt.tmp')
        if (index(file,'opC').ne.0) call system('cat weather.txt | grep -A 2 -e "METOP-C" - > hrpt.tmp')
      endif
      if (csat.eq.'M2') then
        if (index(file,'MN2-2').ne.0) then
          call system('cat weather.txt | grep -A 2 -e "METEOR-M2 2" - > hrpt.tmp')
        else
          if (index(file,'MN2').ne.0) call system('cat weather.txt | grep -A 2 -e "METEOR-M 2" - > hrpt.tmp')
        endif
      endif
      if (csat.eq.'AM') then
        call system('cat weather.txt | grep -A 2 -e "AQUA" - > hrpt.tmp')
      endif
      if (csat.eq.'TE') then
        call system('cat weather.txt | grep -A 2 -e "TERRA" - > hrpt.tmp')
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
c-- Try to fix the "end-of-year" bug
      diy = 0
      if (iyr_file.ne.iyr_tle) then
        write (*,'(''End of year present: '',2I5)') iyr_tle, iyr_file
        diy = 365
        if (mod(iyr_tle,4).eq.0.and.iyr_tle.ne.2000) diy = diy + 1
      endif
c
      return
c
      entry run_tle(linetime,linedoy, ro, vo, long, lat, theta0g)
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

c-- This is a directory monitoring subroutine for batch processing binary files and waiting for new XHRPT .bin or .raw16 files to appear
      subroutine monitor(path, filename)
      implicit none
      character*(*) path, filename
c
      logical   swdone
      integer*1 buf(10000)
      integer*4 i, init, lunbatch, iosbatch, nstore, ios
      character*250 filedatain
      character*250 filedatastore(1000)
      data init/0/
c
      save
c-- Some good old fashioned spaghetti code :-( - As long as it works is the motto !
      if (init.eq.0) then
        call get_lun(lunbatch)
        nstore = 0
        init   = 1
      endif
 4    call system('ls -a1tr '//path(1:lnblnk(path))//'*.bin 2> /dev/null > monitor.txt')
      call system('ls -a1tr '//path(1:lnblnk(path))//'*.raw16 2> /dev/null >> monitor.txt')
      open (unit=lunbatch,file='monitor.txt',form='formatted')
 3    read (lunbatch,'(a)',iostat=iosbatch) filedatain
      if (iosbatch.ne.0) goto 1
      swdone = .false.
      do i = 1,nstore
        if (filedatain(1:lnblnk(filedatain)).eq.filedatastore(i)(1:lnblnk(filedatastore(i)))) then
          swdone = .true.
          goto 3
        endif
      enddo
c-- You ONLY get here when this file was not done - OR if nstore = 0 !
 2    nstore = nstore + 1
      filedatastore(nstore) = filedatain(1:lnblnk(filedatain))
      filename = filedatain(1:lnblnk(filedatain))
      close (unit=lunbatch)
c-- In a terminal, do a touch -f killhrpt in the data directory, and the program will stop when reaching that file
      if (index(filename,'killhrpt').ne.0) stop 'User triggered HRPT monitor stop'
c-- Check if the file can be opened ! If not, sleep !
 5    ios = 0
      open (unit=lunbatch,file=filename,form='unformatted',access='direct',recl=10000,iostat=ios)
      if (ios.ne.0) then
        close (unit=lunbatch)
        write (*,'(''Waiting for file to become openable'')')
        call system('sleep 10')
        goto 5
      endif
      ios = 0
      read (lunbatch,rec=1,iostat=ios)
      if (ios.ne.0) then
        close (unit=lunbatch)
        write (*,'(''Waiting for file to become readable'')')
        call system('sleep 10')
        goto 5
      endif
      close (unit=lunbatch)
      return
c-- Either all files were done or there was no file yet - wait and start again
 1    close (unit=lunbatch)
      call system('sleep 5')
      goto 4
      return
      end
      
      
c-- 17-Jan-2019 - modified for WGS84
      subroutine scan_track(ro, vo, long, lat, linetime, linedoy, longlat, swswath, theta0g, thetascan, etascan, swszalim, szalim, swnorth, swtheta, swdebug, swmn2, nsat_pix, sat_fov)
      implicit none
      real*4    sat_fov
      real*8    ro(3), vo(3), long, lat, linetime, longlat(3,5417), thetascan, etascan, szalim
      integer*4 linedoy, nsat_pix
      logical swswath, swszalim, swnorth, swtheta, swdebug, swmn2
c
      real*8 torad, todeg, x1, y1, z1, d1, d2, d3,  c1, c2, c3, a, b, c, xs1, ys1, zs1, angle_save
      real*8 xv1, yv1, zv1, xv2, yv2, zv2, angle,  longss, latss, angle_threshold, ystep, long_save, lat_save
      real*8 long_v, lat_v, c1000, a_r, b_r, xyz(3), lla(3), theta0g, twopi, phicor, angle1, angle2
      real*8 ro_copy(3), vo_copy(3), d3br, c21, c22, c31, c32, c33, h, onepi, d3i, scandist, depsilon
      real*8 aa(3), ba(3), ca(3), theta, long_int, lat_int, linetime_ref, erad
      real*8 lastlongplus, lastlatplus, lastlongmin, lastlatmin, satalt, eta
      integer*4 ind_pix, ind_pix_save, isgn, nstep, nstep_save, i, j, init, i_mid
c
      logical       swaqm1000, swaqm500, swaqm250, swaqmpan
      common /MODIS/swaqm1000, swaqm500, swaqm250, swaqmpan
c
      data init/0/
c
      save
c
      if (init.eq.0) call darkness_init()
c
c-- The FOV of 55.37 deg (defined through sat_fov) is the half width of the band scanned as seen from the satellite - 0.95 mrad sample step ! (ie .lt. IFOV)
c-- 20190805 - Just when I thought I never have to touch this routine again ..... What to do for Meteor ? 
c-- 20190806 - Made an attempt by using swmn2 and nsat_pix as argument and defining i_mid 
      angle_threshold = dble(sat_fov) / 2.0D0
      i_mid           = (nsat_pix / 2) + 1
      scandist        = 2.0D0 * angle_threshold / dble(nsat_pix)
      depsilon        = 1.0D0 * scandist
      a_r             = 6378137.0D0
      b_r             = 6356752.315D0
      c1000           = 1000.0D0
      torad           = 2.0D0 * dasin(1.0D0) / 180.0D0
      todeg           = 180.0D0 / (2.0D0 * dasin(1.0D0))
      twopi           = 4.0D0 * dasin(1.0D0)
      onepi           = 2.0D0 * dasin(1.0D0)
c
      do i = 1, 5417
        do j = 1,3
          longlat(j,i) = 0.0D0
        enddo
      enddo
      if (etascan.le.1.0D-6) then
        longlat(1,i_mid) = long
        longlat(2,i_mid) = lat
      endif
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
      isgn         = 1
      ystep        = 10.0D0
 2    continue
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
      satalt = lla(3)
      lla(3) = 0.0D0
      call lla2xyz(lla, xyz)
      x1     = xyz(1)
      y1     = xyz(2)
      z1     = xyz(3)
c-- Experimental - for MN2 - Rotate the vector from nadir point on earth to satellite around the vector d - angle is Eta
      eta    = etascan
      if (eta.gt.1.0D-6) then
        aa(1)  = d1
        aa(2)  = d2
        aa(3)  = d3
c-- ba is THE vector from Satellite to Nadir point ! this is the vector to compare with for binning the scanline
        ba(1)  = ro_copy(1) - x1
        ba(2)  = ro_copy(2) - y1
        ba(3)  = ro_copy(3) - z1
c--
        call rotate_vector_a_around_b(ba, aa, eta, ca)
        x1  = ro_copy(1) - ca(1)
        y1  = ro_copy(2) - ca(2)
        z1  = ro_copy(3) - ca(3)
c-- Put the x vector back on the Earth's surface by setting latitude (= lla(3) ) to zero
        xyz(1) = x1
        xyz(2) = y1
        xyz(3) = z1
        call xyz2lla(xyz, lla)
        lla(3) = 0.0D0
        call lla2xyz(lla, xyz)
        x1     = xyz(1)
        y1     = xyz(2)
        z1     = xyz(3)
      endif
c-- End of Eta rotation
c--
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
        if (swdebug) write (*,*) linetime_ref, linetime, theta, long, lat
      endif
c--
      call rotate_vector_a_around_b(aa, ba, theta, ca)
      d1     = ca(1)
      d2     = ca(2)
      d3     = ca(3)
c-- End of Theta rotation
c--
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
c-  Shift Y to find pixel boundaries - this should be OK as I rotate back the true satellite position to long = 0 and later rotate the solution back
c-- Silly construct to restart at original stepsize after sign reversal
      ystep        = dsign(10.0D0, ystep)
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
            if (init.eq.0.and.ind_pix_save.eq.-1.and.swdebug) then
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
              if (ind_pix.le.i_mid-1) then
                longlat(1,i_mid - isgn * ind_pix) = long_v
                longlat(2,i_mid - isgn * ind_pix) = lat_v
c                if (swswath) call darkness(longss, latss, linetime, linedoy, longlat(3,i_mid - isgn * ind_pix))
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
c--
c-- remember the 4 lines of code immediately below are for eta not equal zero !
      if (longlat(1,i_mid).eq.0.0D0.or.longlat(2,i_mid).eq.0.) then
        longlat(1,i_mid) = (longlat(1,i_mid-1) + longlat(1,i_mid+1)) / 2.0D0
        longlat(2,i_mid) = (longlat(2,i_mid-1) + longlat(2,i_mid+1)) / 2.0D0
      endif
c
      lastlongplus = longlat(1,i_mid)
      lastlatplus  = longlat(2,i_mid)
      lastlongmin  = longlat(1,i_mid)
      lastlatmin   = longlat(2,i_mid)
      do i = 1,i_mid - 1
        if (longlat(1,i_mid-i).eq.0.0D0) then
          longlat(1,i_mid-i) = lastlongmin
        else
          lastlongmin  = longlat(1,i_mid-i)
        endif
        if (longlat(1,i_mid+i).eq.0.0D0) then
          longlat(1,i_mid+i) = lastlongplus
        else
          lastlongplus = longlat(1,i_mid+i)
        endif
        if (longlat(2,i_mid-i).eq.0.0D0) then
          longlat(2,i_mid-i) = lastlatmin
        else
          lastlatmin  = longlat(2,i_mid-i)
        endif
        if (longlat(2,i_mid+i).eq.0.0D0) then
          longlat(2,i_mid+i) = lastlatplus
        else
          lastlatplus = longlat(2,i_mid+i)
        endif
      enddo
c-- Mark out the swath, calculate solar zenith angle if needed  etc
      if (swszalim) then
        do i = 1, nsat_pix
          call darkness((longlat(1,i)+longlat(1,i+1))/2.0D0, (longlat(2,i)+longlat(2,i+1))/2.0D0, linetime, linedoy, longlat(3,i))
        enddo
        if (szalim.ge.0.0D0) then
          do i = 1, nsat_pix
            if (longlat(3,i).le.szalim.and.szalim.lt.longlat(3,i+1)) then
              call map_stereo_mark((longlat(1,i)+longlat(1,i+1))/2.0D0, (longlat(2,i)+longlat(2,i+1))/2.0D0, 2)
            endif
            if (longlat(3,i).gt.szalim.and.szalim.ge.longlat(3,i+1)) then
              call map_stereo_mark((longlat(1,i)+longlat(1,i+1))/2.0D0, (longlat(2,i)+longlat(2,i+1))/2.0D0, 2)
            endif
          enddo
        endif
      endif
      if (swswath) then
        do i = 1, nsat_pix
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

c-- Separate AQM version of scan_track - with ascending theta over scan line
      subroutine scan_track_aqm(ro, vo, long, lat, linetime, linedoy, longlat, swswath, theta0g, thetascan, etascan, swszalim, szalim, swnorth, swtheta, swdebug, swmn2, nsat_pix, sat_fov, aqmtheta, bowgamma)
      implicit none
      real*4    sat_fov
      real*8    ro(3), vo(3), long, lat, linetime, longlat(3,5417), thetascan, etascan, szalim, aqmtheta(2), bowgamma
      integer*4 linedoy, nsat_pix
      logical swswath, swszalim, swnorth, swtheta, swdebug, swmn2
c
      real*8 torad, todeg, x1, y1, z1, d1, d2, d3,  c1, c2, c3, a, b, c, xs1, ys1, zs1, angle_save
      real*8 xv1, yv1, zv1, xv2, yv2, zv2, angle,  longss, latss, angle_threshold, ystep, long_save, lat_save
      real*8 long_v, lat_v, c1000, a_r, b_r, xyz(3), lla(3), theta0g, twopi, phicor, angle1, angle2
      real*8 ro_copy(3), vo_copy(3), d3br, c21, c22, c31, c32, c33, h, onepi, d3i, scandist, depsilon
      real*8 aa(3), ba(3), ca(3), theta, long_int, lat_int, linetime_ref, erad
      real*8 lastlongplus, lastlatplus, lastlongmin, lastlatmin, satalt, eta, thetastep
      integer*4 ind_pix, ind_pix_save, isgn, nstep, nstep_save, i, j, init, i_mid
c
      logical       swaqm1000, swaqm500, swaqm250, swaqmpan
      common /MODIS/swaqm1000, swaqm500, swaqm250, swaqmpan
c
      data init/0/
c
      save
c
      if (init.eq.0) call darkness_init()
c
c-- The FOV of 55.37 deg (defined through sat_fov) is the half width of the band scanned as seen from the satellite - 0.95 mrad sample step ! (ie .lt. IFOV)
c-- 20190805 - Just when I thought I never have to touch this routine again ..... What to do for Meteor ? 
c-- 20190806 - Made an attempt by using swmn2 and nsat_pix as argument and defining i_mid 
      angle_threshold = dble(sat_fov) / 2.0D0
      i_mid           = (nsat_pix / 2) + 1
      scandist        = 2.0D0 * angle_threshold / dble(nsat_pix)
      depsilon        = 1.0D0 * scandist
      a_r             = 6378137.0D0
      b_r             = 6356752.315D0
      c1000           = 1000.0D0
      torad           = 2.0D0 * dasin(1.0D0) / 180.0D0
      todeg           = 180.0D0 / (2.0D0 * dasin(1.0D0))
      twopi           = 4.0D0 * dasin(1.0D0)
      onepi           = 2.0D0 * dasin(1.0D0)
c
      do i = 1, 5417
        do j = 1,3
          longlat(j,i) = 0.0D0
        enddo
      enddo
      if (etascan.le.1.0D-6) then
        longlat(1,i_mid) = long
        longlat(2,i_mid) = lat
      endif
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
c-- If for bowtie correction I need to change the sign of the bowtie-theta between left side of scan and right side of scan the outer loop (label 2) should start here
c
      isgn         = 1
      ystep        = 10.0D0
 2    continue
c
c-- Set up the the bow tie correction theta for MODIS if needed .... let's see if it works
      if (swnorth) then
        if (swaqm1000.or.swaqm500.or.swaqm250.or.swaqmpan) then
          thetastep = abs(aqmtheta(2) - aqmtheta(1)) / (nsat_pix - i_mid)
          if (ystep.gt.0.0) then
            theta  = thetascan + aqmtheta(1)
          else
            theta  = thetascan - aqmtheta(1)
          endif
        else
          theta  = thetascan
        endif
      else
        if (swaqm1000.or.swaqm500.or.swaqm250.or.swaqmpan) then
          thetastep = - abs(aqmtheta(2) - aqmtheta(1)) / (nsat_pix - i_mid)
          if (ystep.gt.0.0) then
            theta  = thetascan + aqmtheta(1)
          else
            theta  = thetascan - aqmtheta(1)
          endif
        else
          theta  = thetascan
        endif
      endif
c-  Shift Y to find pixel boundaries - this should be OK as I rotate back the true satellite position to long = 0 and later rotate the solution back
c-- Silly construct to restart at original stepsize after sign reversal
      ystep        = dsign(10.0D0, ystep)
      ys1          = y1 - ystep
      angle        = 0.0D0
      ind_pix_save = -1
      nstep        = 0
      nstep_save   = 0
      do while (angle.le.angle_threshold + depsilon)
c-- TEMP !
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
        satalt = lla(3)
        lla(3) = 0.0D0
        call lla2xyz(lla, xyz)
        x1     = xyz(1)
        y1     = xyz(2)
        z1     = xyz(3)
c-- Experimental - for MN2 - Rotate the vector from nadir point on earth to satellite around the vector d - angle is Eta
        eta    = etascan
        if (eta.gt.1.0D-6) then
          aa(1)  = d1
          aa(2)  = d2
          aa(3)  = d3
c-- ba is THE vector from Satellite to Nadir point ! this is the vector to compare with for binning the scanline
          ba(1)  = ro_copy(1) - x1
          ba(2)  = ro_copy(2) - y1
          ba(3)  = ro_copy(3) - z1
c--
          call rotate_vector_a_around_b(ba, aa, eta, ca)
          x1  = ro_copy(1) - ca(1)
          y1  = ro_copy(2) - ca(2)
          z1  = ro_copy(3) - ca(3)
c-- Put the x vector back on the Earth's surface by setting latitude (= lla(3) ) to zero
          xyz(1) = x1
          xyz(2) = y1
          xyz(3) = z1
          call xyz2lla(xyz, lla)
          lla(3) = 0.0D0
          call lla2xyz(lla, xyz)
          x1     = xyz(1)
          y1     = xyz(2)
          z1     = xyz(3)
        endif
c-- End of Eta rotation
c--
c-- Experimental - Rotate the vector d around the vector from nadir point on earth to satellite
        aa(1)  = d1
        aa(2)  = d2
        aa(3)  = d3
c-- ba is THE vector from Satellite to Nadir point ! this is the vector to compare with for binning the scanline
        ba(1)  = ro_copy(1) - x1
        ba(2)  = ro_copy(2) - y1
        ba(3)  = ro_copy(3) - z1
c--
        call rotate_vector_a_around_b(aa, ba, theta, ca)
        d1     = ca(1)
        d2     = ca(2)
        d3     = ca(3)
c-- End of Theta rotation
c--
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
c-- TEMP !
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
            if (init.eq.0.and.ind_pix_save.eq.-1.and.swdebug) then
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
              if (ind_pix.le.i_mid-1) then
                longlat(1,i_mid - isgn * ind_pix) = long_v
                longlat(2,i_mid - isgn * ind_pix) = lat_v
c                if (swswath) call darkness(longss, latss, linetime, linedoy, longlat(3,i_mid - isgn * ind_pix))
              endif
              if (nstep-nstep_save.ge.10) then
                ystep = ystep * 1.414214
              endif
              nstep_save   = nstep
              ind_pix_save = ind_pix
              long_save    = longss
              lat_save     = latss
              angle_save   = angle
c-- Add the bow tie correction theta for MODIS if needed .... let's see if it works - do it here so it runs once per line !
              if (swaqm1000.or.swaqm500.or.swaqm250.or.swaqmpan) then
                if (bowgamma.gt.0.001) then
                  if (swnorth) then
                    if (ystep.gt.0.0) then
                      theta  = thetascan + aqmtheta(1) + (aqmtheta(2)-aqmtheta(1)) * (dabs((dble(2*(i_mid - isgn * ind_pix))/dble(nsat_pix)) - 1.0D0) ** (1.0D0/bowgamma))
                    else
                      theta  = thetascan - aqmtheta(1) - (aqmtheta(2)-aqmtheta(1)) * (dabs((dble(2*(i_mid - isgn * ind_pix))/dble(nsat_pix)) - 1.0D0) ** (1.0D0/bowgamma))
                    endif
                  else
                    if (ystep.gt.0.0) then
                      theta  = thetascan - aqmtheta(1) - (aqmtheta(2)-aqmtheta(1)) * (dabs((dble(2*(i_mid - isgn * ind_pix))/dble(nsat_pix)) - 1.0D0) ** (1.0D0/bowgamma))
                    else
                      theta  = thetascan + aqmtheta(1) + (aqmtheta(2)-aqmtheta(1)) * (dabs((dble(2*(i_mid - isgn * ind_pix))/dble(nsat_pix)) - 1.0D0) ** (1.0D0/bowgamma))
                    endif
                  endif
                else
                  if (ystep.gt.0.0) then
                    if (aqmtheta(1).gt.0.) then
                      theta  = theta + thetastep
                    else
                      theta  = theta - thetastep
                    endif
                  else
                    if (aqmtheta(1).gt.0.) then
                      theta  = theta - thetastep
                    else
                      theta  = theta + thetastep
                    endif
                  endif
                endif
              else
                theta  = thetascan
              endif
c-- Two check prints to check if theta behaves nicely
c              if (ind_pix.eq. 500) write (*,*) ind_pix, aqmtheta(1), aqmtheta(2), thetastep, thetascan, theta, ystep, longss, latss
c              if (ind_pix.eq.2500) write (*,*) ind_pix, aqmtheta(1), aqmtheta(2), thetastep, thetascan, theta, ystep, longss, latss
            else
              ind_pix_save = ind_pix
              long_save    = longss
              lat_save     = latss
              angle_save   = angle
            endif
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
c--
c-- remember the 4 lines of code immediately below are for eta not equal zero !
      if (longlat(1,i_mid).eq.0.0D0.or.longlat(2,i_mid).eq.0.) then
        longlat(1,i_mid) = (longlat(1,i_mid-1) + longlat(1,i_mid+1)) / 2.0D0
        longlat(2,i_mid) = (longlat(2,i_mid-1) + longlat(2,i_mid+1)) / 2.0D0
      endif
c
      lastlongplus = longlat(1,i_mid)
      lastlatplus  = longlat(2,i_mid)
      lastlongmin  = longlat(1,i_mid)
      lastlatmin   = longlat(2,i_mid)
      do i = 1,i_mid - 1
        if (longlat(1,i_mid-i).eq.0.0D0) then
          longlat(1,i_mid-i) = lastlongmin
        else
          lastlongmin  = longlat(1,i_mid-i)
        endif
        if (longlat(1,i_mid+i).eq.0.0D0) then
          longlat(1,i_mid+i) = lastlongplus
        else
          lastlongplus = longlat(1,i_mid+i)
        endif
        if (longlat(2,i_mid-i).eq.0.0D0) then
          longlat(2,i_mid-i) = lastlatmin
        else
          lastlatmin  = longlat(2,i_mid-i)
        endif
        if (longlat(2,i_mid+i).eq.0.0D0) then
          longlat(2,i_mid+i) = lastlatplus
        else
          lastlatplus = longlat(2,i_mid+i)
        endif
      enddo
c-- Mark out the swath, calculate solar zenith angle if needed  etc
      if (swszalim) then
        do i = 1, nsat_pix
          call darkness((longlat(1,i)+longlat(1,i+1))/2.0D0, (longlat(2,i)+longlat(2,i+1))/2.0D0, linetime, linedoy, longlat(3,i))
        enddo
        if (szalim.ge.0.0D0) then
          do i = 1, nsat_pix
            if (longlat(3,i).le.szalim.and.szalim.lt.longlat(3,i+1)) then
              call map_stereo_mark((longlat(1,i)+longlat(1,i+1))/2.0D0, (longlat(2,i)+longlat(2,i+1))/2.0D0, 2)
            endif
            if (longlat(3,i).gt.szalim.and.szalim.ge.longlat(3,i+1)) then
              call map_stereo_mark((longlat(1,i)+longlat(1,i+1))/2.0D0, (longlat(2,i)+longlat(2,i+1))/2.0D0, 2)
            endif
          enddo
        endif
      endif
      if (swswath) then
        do i = 1, nsat_pix
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
      
      subroutine create_png(file, ext1, ext2, ext3, n, option, swnorth, npix, nlines_in, swbatch, swsharp, csharp)
      implicit none
      logical   swnorth, swbatch, swsharp
      integer*4 n, npix, nlines_in
      character*(*) file, ext1, ext2, ext3, option, csharp
c
      logical swlevelrgb
      integer*4 j, nlines, i
      character*1 csq, cdq
      character*450 cstring2
      character*1000 imgname, command, cinput, cstring, cfilebatch
c
      if (nlines_in.lt.0) then
        call get_command(cinput)
        csq = char(ichar("'"))
        cdq = char(ichar('"'))
        cstring = ' -pointsize 20 -fill red -draw '//csq//'text 50,25 '//cdq//cinput(1:lnblnk(cinput))//cdq//csq
        if (swbatch) then
          call get_string(cstring2)
          do i = 1,len(cfilebatch)
            cfilebatch(i:i) = ' '
          enddo
          cfilebatch = ' -pointsize 20 -fill red -draw '//csq//'text 50,50 '//cdq//cstring2(1:lnblnk(cstring2))//cdq//csq
        else
          do i = 1,len(cstring2)
            cstring2(i:i) = ' '
          enddo
        endif
      else
        do i = 1, len(cstring)
          cstring(i:i) = ' '
        enddo
      endif
      swlevelrgb = .false.
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
          write (command,'(''convert -depth 16 -endian lsb -size '',I5.5,''x'',I5.5,1x,a,'' -equalize '',a,'' -flip -depth 16 '',a,'' '',a)') 
     *      npix, nlines, 'gray:'//imgname(1:lnblnk(imgname)), option(1:lnblnk(option)), cstring(1:lnblnk(cstring)), imgname(1:j-1)//'.png'
        else
          if (index(ext3,'true').eq.0) then
c-- The fix for a second line of text to identify the filename in 'batch' mode only needs to be applied here as for the polar stereographic projection swnorth is always true because of the flip required (up <-> down)
            if (nlines_in.lt.0.and.swbatch) then
              write (command,'(''convert -depth 16 -endian lsb -size '',I5.5,''x'',I5.5,1x,a,'' -equalize '',a,'' -flip -depth 16 '',a,'' '',a,'' '',a)') 
     *          npix, nlines, 'rgb:'//imgname(1:lnblnk(imgname)), option(1:lnblnk(option)), cstring(1:lnblnk(cstring)), cfilebatch(1:lnblnk(cfilebatch)), imgname(1:j-1)//'.png'
            else
              write (command,'(''convert -depth 16 -endian lsb -size '',I5.5,''x'',I5.5,1x,a,'' -equalize '',a,'' -flip -depth 16 '',a,'' '',a)') 
     *          npix, nlines, 'rgb:'//imgname(1:lnblnk(imgname)), option(1:lnblnk(option)), cstring(1:lnblnk(cstring)), imgname(1:j-1)//'.png'
            endif
          else
            write (command,'(''convert -depth 16 -endian lsb -size '',I5.5,''x'',I5.5,1x,a,'' -equalize '',a,'' -flip -depth 16 '',a,'' '',a)') 
     *        npix, nlines, 'rgb:'//imgname(1:lnblnk(imgname)), option(1:lnblnk(option)), cstring(1:lnblnk(cstring)), imgname(1:j-1)//'.png'
            if (index(ext3,'_correct').ne.0) swlevelrgb = .true.
          endif
        endif
      else
        if (index(ext2,'.rgb').eq.0) then
          write (command,'(''convert -depth 16 -endian lsb -size '',I5.5,''x'',I5.5,1x,a,'' -equalize '',a,'' -depth 16 '',a,'' '',a)') 
     *      npix, nlines, 'gray:'//imgname(1:lnblnk(imgname)), option(1:lnblnk(option)), cstring(1:lnblnk(cstring)), imgname(1:j-1)//'.png'
        else
          if (index(ext3,'true').eq.0) then
            write (command,'(''convert -depth 16 -endian lsb -size '',I5.5,''x'',I5.5,1x,a,'' -equalize '',a,'' -depth 16 '',a,'' '',a)') 
     *        npix, nlines, 'rgb:'//imgname(1:lnblnk(imgname)), option(1:lnblnk(option)), cstring(1:lnblnk(cstring)), imgname(1:j-1)//'.png'
          else
            write (command,'(''convert -depth 16 -endian lsb -size '',I5.5,''x'',I5.5,1x,a,'' -equalize '',a,'' -depth 16 '',a,'' '',a)') 
     *        npix, nlines, 'rgb:'//imgname(1:lnblnk(imgname)), option(1:lnblnk(option)), cstring(1:lnblnk(cstring)), imgname(1:j-1)//'.png'
            if (index(ext3,'_correct').ne.0) swlevelrgb = .true.
          endif
        endif
      endif
      write (*,*) command(1:lnblnk(command))
      call system(command)
c
      if (swsharp.and.(index(imgname(1:j-1),'_project').ne.0.or.index(imgname(1:j-1),'_border').ne.0)) then
        write (command,'(''convert '',a,'' -unsharp '',a,1x,a)') imgname(1:j-1)//'.png', csharp(1:lnblnk(csharp)), imgname(1:j-1)//'_sharpened'//'.png'
        write (*,*) command(1:lnblnk(command))
        call system(command)
      endif
      if (swlevelrgb) then
        write (command,'(''convert '',a,'' -channel RGB -separate '',a)') imgname(1:j-1)//'.png', 'channel-%d.png'
        write (*,*) command(1:lnblnk(command))
        call system(command)
        write (command,'(''convert channel-0.png -level 0.0%,80%,1.4 channel-0-mod.png'')')
        write (*,*) command(1:lnblnk(command))
        call system(command)
        write (command,'(''convert channel-1.png -level 0.0%,87%,1.2 channel-1-mod.png'')')
        write (*,*) command(1:lnblnk(command))
        call system(command)
        write (command,'(''convert channel-2.png -level 0.0%,95%,0.95 channel-2-mod.png'')')
        write (*,*) command(1:lnblnk(command))
        call system(command)
        write (command,'(''convert channel-0-mod.png channel-1-mod.png channel-2-mod.png -combine '',a)') imgname(1:j-1)//'-mod.png'
        write (*,*) command(1:lnblnk(command))
        call system(command)
        call system('rm channel*.png')
      endif
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
      character*11 csize
c
      save
c
      cmap = 'satmap' 
      write (csize,'(I5.5,''x'',I5.5)') nx, ny
c
      call get_lun(lun)
      open (unit=lun, file='./resource/worldmaphr.rgb', access='direct', form='unformatted', recl=3*nx)
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

      subroutine map_stereo_init(swfast)
      implicit none
      logical swfast, swdebug
c--
      integer*4 nxs, nys, nazel
      parameter (nxs=5000, nys=5000, nazel=1400)
c--
      logical   localfast
      integer*1 mapsline (3,nxs), maps(3,nxs,nys), masks(nxs, nys)
      integer*4 lun, i, j, ix, iy, itype, value, ixmin, iymin, ixmax, iymax, nmasked, nxsat, nysat, nzoom, nxz
      integer*4 ixminz, iyminz, ixmaxz, iymaxz, nx_map, ny_map, nc, nx_in, ny_in, ix1, iy1, ix2, iy2, ixdir, iydir
      integer*4 longstep, latstep, ilo1, ilo2, ila1, ila2, ncount, i4v1, i4v2, i4v3, iclim, i_merge
      integer*2 stereo_image(nc,nx_in,ny_in), pix(3), mark_r, mark_g, mark_b, mark_longlat_r, mark_longlat_g, mark_longlat_b, mr, mg, mb
      integer*1 stereo_count(nx_in, ny_in), stereo_count_rgb(nx_in, ny_in, 3)
      real*4    y_delta, x_delta, longs, lats, xs, ys, xr, yr, phi
      real*8    long, lat, torad, todeg, pi, pid4, long1, lat1, long2, lat2, rzoom, az, el, longmin, longmax, latmin, latmax, longval, rzoom_local
      character*1 csq, cdq
      character*6 cmap
      character*11 csize
      character*23 cextract
      character*50 argstring
      character*500 cstring
      character*(*) command
c
      save
c
      localfast = swfast
c
      cmap = 'satmap' 
      write (csize,'(I5.5,''x'',I5.5)') nxs, nys
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
      ncount  = 0
      rzoom_local = 0.0D0
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
c-- See if a rotation of the polar stereographic is required
      phi = 0.0
      do i = 3, 30
        call getarg(i,argstring)
        if (index(argstring,'reflong=').ne.0) then
          j = index(argstring,'=')
          read (argstring(j+1:lnblnk(argstring)),*) phi
          phi = phi * torad
        endif
      enddo
      return
c
      entry map_stereo_mark(long, lat, itype)
      longs = long * torad
      lats  =  lat * torad
      xs =  2.0 * float(nxs)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
      ys = -2.0 * float(nxs)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
      xr =   cos(phi) * xs - sin(phi) * ys
      yr =  (sin(phi) * xs + cos(phi) * ys)
      xs = xr
      ys = yr
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
      xr =   cos(phi) * xs - sin(phi) * ys
      yr =  (sin(phi) * xs + cos(phi) * ys)
      xs = xr
      ys = yr
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
      entry map_stereo_zoom_init(nxsat, nysat, nx_map, ny_map, swdebug, iclim)
      nmasked = 0
      do j = 1, nys
        do i = 1,nxs
          if (masks(i,j).eq.1) nmasked = nmasked + 1
        enddo
      enddo
c-- changed to the rzoom
      nzoom = nint(sqrt(float(nxsat*nysat)/float(nmasked)))
      rzoom =      sqrt(float(nxsat*nysat)/float(nmasked)) * 1.1
      if (rzoom_local.ne.0.0D0) then
        rzoom = rzoom_local
      endif
      if (swdebug) write (*,*) ixmin, ixmax, iymin, iymax, nmasked, nxsat*nysat, int(float(nxsat*nysat)/float(nmasked)), nzoom, rzoom
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
      longmin =  400.0D0
      longmax = -400.0D0
      latmin  =  100.0D0
      latmax  = -100.0D0
      mark_r  = iclim
      mark_g  = 0
      mark_b  = 0
      mark_longlat_r = 1012
      mark_longlat_g = 40
      mark_longlat_b = 692
c
      nx_map = ixmaxz - ixminz + 1
      ny_map = iymaxz - iyminz + 1
      if (swdebug) write (*,*) ixminz, ixmaxz, iyminz, iymaxz, nxz, nx_map*ny_map
      return
c
      entry map_stereo_rzoom(nxsat, nysat)
      nmasked = 0
      do j = 1, nys
        do i = 1,nxs
          if (masks(i,j).eq.1) nmasked = nmasked + 1
        enddo
      enddo
c-- changed to the rzoom
      rzoom_local = sqrt(float(nxsat*nysat)/float(nmasked)) * 1.1
      return

c
      entry map_stereo_zoom_mark(stereo_image, nc, nx_in, ny_in,long1, lat1, long2, lat2)
      longs = long1 * torad
      lats  =  lat1 * torad
      xs =  2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
      ys = -2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
      xr =   cos(phi) * xs - sin(phi) * ys
      yr =  (sin(phi) * xs + cos(phi) * ys)
      xs = xr
      ys = yr
      ix1 = min(max(int(xs)+(nxz/2),1),nxz)
      iy1 = min(max(int(ys)+(nxz/2),1),nxz)
      ix1 = ix1 - ixminz + 1
      iy1 = iy1 - iyminz + 1
      longs = long2 * torad
      lats  =  lat2 * torad
      xs =  2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
      ys = -2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
      xr =   cos(phi) * xs - sin(phi) * ys
      yr =  (sin(phi) * xs + cos(phi) * ys)
      xs = xr
      ys = yr
      ix2 = min(max(int(xs)+(nxz/2),1),nxz)
      iy2 = min(max(int(ys)+(nxz/2),1),nxz)
      ix2 = ix2 - ixminz + 1
      iy2 = iy2 - iyminz + 1
c-- Single Point
      if (ix1.eq.ix2.and.iy1.eq.iy2) then
        ix = ix1
        iy = iy1
        if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
          stereo_image(1,ix,iy) = mark_r
          stereo_image(2,ix,iy) = mark_g
          stereo_image(3,ix,iy) = mark_b
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
            stereo_image(1,ix,iy) = mark_r
            stereo_image(2,ix,iy) = mark_g
            stereo_image(3,ix,iy) = mark_b
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
            stereo_image(1,ix,iy) = mark_r
            stereo_image(2,ix,iy) = mark_g
            stereo_image(3,ix,iy) = mark_b
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
            stereo_image(1,ix,iy) = mark_r
            stereo_image(2,ix,iy) = mark_g
            stereo_image(3,ix,iy) = mark_b
          endif
        enddo
        goto 1
      else
        iydir = 1
        if (iy2.lt.iy1) iydir = -1
        do iy = iy1, iy2, iydir
          ix = nint(float(ix1) + ((float(iy-iy1)/float(iy2-iy1))*float(ix2-ix1)))
          if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
            stereo_image(1,ix,iy) = mark_r
            stereo_image(2,ix,iy) = mark_g
            stereo_image(3,ix,iy) = mark_b
          endif
        enddo
        goto 1
      endif
 1    continue
      return
c
      entry map_stereo_zoom_place(stereo_image,nc,nx_in,ny_in,long,lat,pix)
      longs = long * torad
      lats  =  lat * torad
      xs =  2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
      ys = -2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
      xr =   cos(phi) * xs - sin(phi) * ys
      yr =  (sin(phi) * xs + cos(phi) * ys)
      xs = xr
      ys = yr
      ix = min(max(int(xs)+(nxz/2),1),nxz)
      iy = min(max(int(ys)+(nxz/2),1),nxz)
      ix = ix - ixminz + 1
      iy = iy - iyminz + 1
      if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
        stereo_image(1,ix,iy) = pix(1)
        stereo_image(2,ix,iy) = pix(2)
        stereo_image(3,ix,iy) = pix(3)
      endif
c
      if (long.lt.longmin) longmin = long
      if (long.gt.longmax) longmax = long
      if (lat.lt.latmin)   latmin  = lat
      if (lat.gt.latmax)   latmax  = lat
      return
c
      entry map_stereo_zoom_place_aqm(stereo_image,stereo_count,nc,nx_in,ny_in,long,lat,pix)
      longs = long * torad
      lats  =  lat * torad
      xs =  2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
      ys = -2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
      xr =   cos(phi) * xs - sin(phi) * ys
      yr =  (sin(phi) * xs + cos(phi) * ys)
      xs = xr
      ys = yr
      ix = min(max(int(xs)+(nxz/2),1),nxz)
      iy = min(max(int(ys)+(nxz/2),1),nxz)
      ix = ix - ixminz + 1
      iy = iy - iyminz + 1
      if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
        i4v2 = stereo_count(ix,iy)
        if (i4v2.lt.0) i4v2 = i4v2 + 256
        if (i4v2.eq.0) then
          stereo_image(1,ix,iy) = pix(1)
          stereo_image(2,ix,iy) = pix(2)
          stereo_image(3,ix,iy) = pix(3)
        else if (i4v2.eq.1) then
          stereo_image(1,ix,iy) = (stereo_image(1,ix,iy) + pix(1)) / 2
          stereo_image(2,ix,iy) = (stereo_image(2,ix,iy) + pix(2)) / 2
          stereo_image(3,ix,iy) = (stereo_image(3,ix,iy) + pix(3)) / 2
        else if (i4v2.ge.2.and.i4v2.lt.255) then
          i4v1 = stereo_image(1,ix,iy)
          i4v2 = stereo_count(ix,iy)
          if (i4v2.lt.0) i4v2 = i4v2 + 256
          i4v3 = i4v2 * i4v1
          i4v1 = pix(1)
          i4v3 = i4v3 + i4v1
          i4v2 = i4v2 + 1
          i4v3 = i4v3 / i4v2
          stereo_image(1,ix,iy) = i4v3
          i4v1 = stereo_image(2,ix,iy)
          i4v2 = stereo_count(ix,iy)
          if (i4v2.lt.0) i4v2 = i4v2 + 256
          i4v3 = i4v2 * i4v1
          i4v1 = pix(2)
          i4v3 = i4v3 + i4v1
          i4v2 = i4v2 + 1
          i4v3 = i4v3 / i4v2
          stereo_image(2,ix,iy) = i4v3
          i4v1 = stereo_image(3,ix,iy)
          i4v2 = stereo_count(ix,iy)
          if (i4v2.lt.0) i4v2 = i4v2 + 256
          i4v3 = i4v2 * i4v1
          i4v1 = pix(3)
          i4v3 = i4v3 + i4v1
          i4v2 = i4v2 + 1
          i4v3 = i4v3 / i4v2
          stereo_image(3,ix,iy) = i4v3
c          stereo_image(1,ix,iy) = (stereo_image(1,ix,iy) * stereo_count(ix,iy) + pix(1)) / (stereo_count(ix,iy) + 1)
c          stereo_image(2,ix,iy) = (stereo_image(2,ix,iy) * stereo_count(ix,iy) + pix(2)) / (stereo_count(ix,iy) + 1)
c          stereo_image(3,ix,iy) = (stereo_image(3,ix,iy) * stereo_count(ix,iy) + pix(3)) / (stereo_count(ix,iy) + 1)
        endif
        i4v1 = stereo_count(ix,iy)
        if (i4v1.lt.0) i4v1 = i4v1 + 256
        if (i4v1.lt.255) then
          i4v1 = i4v1 + 1
          if (i4v1.gt.127) i4v1 = i4v1 - 256
          stereo_count(  ix,iy) = i4v1
        endif
      endif
c
      if (long.lt.longmin) longmin = long
      if (long.gt.longmax) longmax = long
      if (lat.lt.latmin)   latmin  = lat
      if (lat.gt.latmax)   latmax  = lat
      return
c
      entry map_stereo_zoom_aver_aqm(stereo_image,stereo_count,nc,nx_in,ny_in)
      do iy = 1, ny_in
        do ix = 1, nx_in
c          if (stereo_count(ix, iy).gt.0) then
c            stereo_image(1, ix, iy) = nint(float(stereo_image(1, ix, iy))/float(stereo_count(ix, iy)))
c            stereo_image(2, ix, iy) = nint(float(stereo_image(2, ix, iy))/float(stereo_count(ix, iy)))
c            stereo_image(3, ix, iy) = nint(float(stereo_image(3, ix, iy))/float(stereo_count(ix, iy)))
c          endif
          i4v1 = stereo_count(ix,iy)
          if (i4v1.lt.0) i4v1 = i4v1 + 256
          if (i4v1.gt.ncount) ncount = i4v1
        enddo
      enddo
      write (*,'(''Max pixel hit      : '',I10)') ncount
      return
c
      entry map_stereo_zoom_place_aqm_pan(stereo_image,stereo_count_rgb,nc,nx_in,ny_in,long,lat,pix,i_merge)
      longs = long * torad
      lats  =  lat * torad
      xs =  2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
      ys = -2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
      xr =   cos(phi) * xs - sin(phi) * ys
      yr =  (sin(phi) * xs + cos(phi) * ys)
      xs = xr
      ys = yr
      ix = min(max(int(xs)+(nxz/2),1),nxz)
      iy = min(max(int(ys)+(nxz/2),1),nxz)
      ix = ix - ixminz + 1
      iy = iy - iyminz + 1
      if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
        if (stereo_count_rgb(ix,iy,i_merge).eq.0) then
          if (i_merge.eq.1) then
            stereo_image(1,ix,iy) = pix(1)
          endif
          if (i_merge.eq.2) then
            stereo_image(2,ix,iy) = pix(2)
          endif
          if (i_merge.eq.3) then
            stereo_image(3,ix,iy) = pix(3)
          endif
        else if (stereo_count_rgb(ix,iy,i_merge).eq.1) then
          if (i_merge.eq.1) then
            stereo_image(1,ix,iy) = (stereo_image(1,ix,iy) + pix(1)) / 2
          endif
          if (i_merge.eq.2) then
            stereo_image(2,ix,iy) = (stereo_image(2,ix,iy) + pix(2)) / 2
          endif
          if (i_merge.eq.3) then
            stereo_image(3,ix,iy) = (stereo_image(3,ix,iy) + pix(3)) / 2
          endif
        else if (stereo_count_rgb(ix,iy,i_merge).ge.2) then
          if (i_merge.eq.1) then
            i4v1 = stereo_image(1,ix,iy)
            i4v2 = stereo_count_rgb(ix,iy,i_merge)
            if (i4v2.lt.0) i4v2 = i4v2 + 256
            i4v3 = i4v2 * i4v1
            i4v1 = pix(1)
            i4v3 = i4v3 + i4v1
            i4v2 = i4v2 + 1
            i4v3 = i4v3 / i4v2
            stereo_image(1,ix,iy) = i4v3
          endif
          if (i_merge.eq.2) then
            i4v1 = stereo_image(2,ix,iy)
            i4v2 = stereo_count_rgb(ix,iy,i_merge)
            if (i4v2.lt.0) i4v2 = i4v2 + 256
            i4v3 = i4v2 * i4v1
            i4v1 = pix(2)
            i4v3 = i4v3 + i4v1
            i4v2 = i4v2 + 1
            i4v3 = i4v3 / i4v2
            stereo_image(2,ix,iy) = i4v3
          endif
          if (i_merge.eq.3) then
            i4v1 = stereo_image(3,ix,iy)
            i4v2 = stereo_count_rgb(ix,iy,i_merge)
            if (i4v2.lt.0) i4v2 = i4v2 + 256
            i4v3 = i4v2 * i4v1
            i4v1 = pix(3)
            i4v3 = i4v3 + i4v1
            i4v2 = i4v2 + 1
            i4v3 = i4v3 / i4v2
            stereo_image(3,ix,iy) = i4v3
          endif
c          stereo_image(1,ix,iy) = (stereo_image(1,ix,iy) * stereo_count_rgb(ix,iy) + pix(1)) / (stereo_count_rgb(ix,iy) + 1)
c          stereo_image(2,ix,iy) = (stereo_image(2,ix,iy) * stereo_count_rgb(ix,iy) + pix(2)) / (stereo_count_rgb(ix,iy) + 1)
c          stereo_image(3,ix,iy) = (stereo_image(3,ix,iy) * stereo_count_rgb(ix,iy) + pix(3)) / (stereo_count_rgb(ix,iy) + 1)
        endif
        i4v1 = stereo_count_rgb(ix,iy,i_merge)
        if (i4v1.lt.0) i4v1 = i4v1 + 256
        i4v1 = i4v1 + 1
        if (i4v1.gt.127) i4v1 = i4v1 - 256
        stereo_count_rgb(ix,iy,i_merge) = i4v1
      endif
c
      if (long.lt.longmin) longmin = long
      if (long.gt.longmax) longmax = long
      if (lat.lt.latmin)   latmin  = lat
      if (lat.gt.latmax)   latmax  = lat
      return
c
      entry map_stereo_zoom_aver_aqm_pan(stereo_image,stereo_count_rgb,nc,nx_in,ny_in)
      ncount = 0
      do iy = 1, ny_in
        do ix = 1, nx_in
c          if (stereo_count(ix, iy).gt.0) then
c            stereo_image(1, ix, iy) = nint(float(stereo_image(1, ix, iy))/float(stereo_count(ix, iy)))
c            stereo_image(2, ix, iy) = nint(float(stereo_image(2, ix, iy))/float(stereo_count(ix, iy)))
c            stereo_image(3, ix, iy) = nint(float(stereo_image(3, ix, iy))/float(stereo_count(ix, iy)))
c          endif
          i4v1 = stereo_count_rgb(ix,iy,1)
          if (i4v1.lt.0) i4v1 = i4v1 + 256
          if (i4v1.gt.ncount) ncount = i4v1
        enddo
      enddo
      write (*,'(''Max pixel hit 1    : '',I10)') ncount
      ncount = 0
      do iy = 1, ny_in
        do ix = 1, nx_in
c          if (stereo_count(ix, iy).gt.0) then
c            stereo_image(1, ix, iy) = nint(float(stereo_image(1, ix, iy))/float(stereo_count(ix, iy)))
c            stereo_image(2, ix, iy) = nint(float(stereo_image(2, ix, iy))/float(stereo_count(ix, iy)))
c            stereo_image(3, ix, iy) = nint(float(stereo_image(3, ix, iy))/float(stereo_count(ix, iy)))
c          endif
          i4v1 = stereo_count_rgb(ix,iy,2)
          if (i4v1.lt.0) i4v1 = i4v1 + 256
          if (i4v1.gt.ncount) ncount = i4v1
        enddo
      enddo
      write (*,'(''Max pixel hit 2    : '',I10)') ncount
      ncount = 0
      do iy = 1, ny_in
        do ix = 1, nx_in
c          if (stereo_count(ix, iy).gt.0) then
c            stereo_image(1, ix, iy) = nint(float(stereo_image(1, ix, iy))/float(stereo_count(ix, iy)))
c            stereo_image(2, ix, iy) = nint(float(stereo_image(2, ix, iy))/float(stereo_count(ix, iy)))
c            stereo_image(3, ix, iy) = nint(float(stereo_image(3, ix, iy))/float(stereo_count(ix, iy)))
c          endif
          i4v1 = stereo_count_rgb(ix,iy,3)
          if (i4v1.lt.0) i4v1 = i4v1 + 256
          if (i4v1.gt.ncount) ncount = i4v1
        enddo
      enddo
      write (*,'(''Max pixel hit 3    : '',I10)') ncount
      return
c
      entry map_stereo_zoom_mark_longlat(stereo_image, nc, nx_in, ny_in, longstep, latstep)
      ilo1 = int(longmin/dble(longstep)) * longstep
      ilo2 = int(longmax/dble(longstep)) * longstep
      do i = ilo1, ilo2, longstep
        longs = dble(i) * torad
        lats  =  latmin * torad
        xs =  2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
        ys = -2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
        xr =   cos(phi) * xs - sin(phi) * ys
        yr =  (sin(phi) * xs + cos(phi) * ys)
        xs = xr
        ys = yr
        ix1 = min(max(int(xs)+(nxz/2),1),nxz)
        iy1 = min(max(int(ys)+(nxz/2),1),nxz)
        ix1 = ix1 - ixminz + 1
        iy1 = iy1 - iyminz + 1
        longs = dble(i) * torad
        lats  =  latmax * torad
        xs =  2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
        ys = -2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
        xr =   cos(phi) * xs - sin(phi) * ys
        yr =  (sin(phi) * xs + cos(phi) * ys)
        xs = xr
        ys = yr
        ix2 = min(max(int(xs)+(nxz/2),1),nxz)
        iy2 = min(max(int(ys)+(nxz/2),1),nxz)
        ix2 = ix2 - ixminz + 1
        iy2 = iy2 - iyminz + 1
c-- Single Point
        if (ix1.eq.ix2.and.iy1.eq.iy2) then
          ix = ix1
          iy = iy1
          if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
            stereo_image(1,ix,iy) = mark_longlat_r
            stereo_image(2,ix,iy) = mark_longlat_g
            stereo_image(3,ix,iy) = mark_longlat_b
          endif
          goto 2
        endif
c-- Vertical  Line
        if (ix1.eq.ix2) then
          ix = ix1
          iydir = 1
          if (iy2.lt.iy1) iydir = -1
          do iy = iy1, iy2, iydir
            if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
              stereo_image(1,ix,iy) = mark_longlat_r
              stereo_image(2,ix,iy) = mark_longlat_g
              stereo_image(3,ix,iy) = mark_longlat_b
            endif
          enddo
          goto 2
        endif
c-- Horizontal  Line
        if (iy1.eq.iy2) then
          iy = iy1
          ixdir = 1
          if (ix2.lt.ix1) ixdir = -1
          do ix = ix1, ix2, ixdir
            if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
              stereo_image(1,ix,iy) = mark_longlat_r
              stereo_image(2,ix,iy) = mark_longlat_g
              stereo_image(3,ix,iy) = mark_longlat_b
            endif
          enddo
          goto 2
        endif
c-- The Rest
        if (abs(ix2-ix1).gt.abs(iy2-iy1)) then
          ixdir = 1
          if (ix2.lt.ix1) ixdir = -1
          do ix = ix1, ix2, ixdir
            iy = nint(float(iy1) + ((float(ix-ix1)/float(ix2-ix1))*float(iy2-iy1)))
            if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
              stereo_image(1,ix,iy) = mark_longlat_r
              stereo_image(2,ix,iy) = mark_longlat_g
              stereo_image(3,ix,iy) = mark_longlat_b
            endif
          enddo
          goto 2
        else
          iydir = 1
          if (iy2.lt.iy1) iydir = -1
          do iy = iy1, iy2, iydir
            ix = nint(float(ix1) + ((float(iy-iy1)/float(iy2-iy1))*float(ix2-ix1)))
            if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
              stereo_image(1,ix,iy) = mark_longlat_r
              stereo_image(2,ix,iy) = mark_longlat_g
              stereo_image(3,ix,iy) = mark_longlat_b
            endif
          enddo
          goto 2
        endif
 2      continue
      enddo
      ila1 = int(latmin/dble(latstep)) * latstep
      ila2 = int(latmax/dble(latstep)) * latstep
      do i = ila1, ila2, latstep
        longval = dble(longmin) - 0.05D0
        do while (longval.lt.dble(longmax))
          longval = longval + 0.05D0
          longs = longval * torad
          lats  = dble(i) * torad
          xs =  2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
          ys = -2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
          xr =   cos(phi) * xs - sin(phi) * ys
          yr =  (sin(phi) * xs + cos(phi) * ys)
          xs = xr
          ys = yr
          ix1 = min(max(int(xs)+(nxz/2),1),nxz)
          iy1 = min(max(int(ys)+(nxz/2),1),nxz)
          ix1 = ix1 - ixminz + 1
          iy1 = iy1 - iyminz + 1
          longs = (longval + 0.05D0) * torad
          lats  = dble(i) * torad
          xs =  2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
          ys = -2.0 * float(nxz)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
          xr =   cos(phi) * xs - sin(phi) * ys
          yr =  (sin(phi) * xs + cos(phi) * ys)
          xs = xr
          ys = yr
          ix2 = min(max(int(xs)+(nxz/2),1),nxz)
          iy2 = min(max(int(ys)+(nxz/2),1),nxz)
          ix2 = ix2 - ixminz + 1
          iy2 = iy2 - iyminz + 1
c-- Single Point
          if (ix1.eq.ix2.and.iy1.eq.iy2) then
            ix = ix1
            iy = iy1
            if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
              stereo_image(1,ix,iy) = mark_longlat_r
              stereo_image(2,ix,iy) = mark_longlat_g
              stereo_image(3,ix,iy) = mark_longlat_b
            endif
            goto 3
          endif
c-- Vertical  Line
          if (ix1.eq.ix2) then
            ix = ix1
            iydir = 1
            if (iy2.lt.iy1) iydir = -1
            do iy = iy1, iy2, iydir
              if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
                stereo_image(1,ix,iy) = mark_longlat_r
                stereo_image(2,ix,iy) = mark_longlat_g
                stereo_image(3,ix,iy) = mark_longlat_b
              endif
            enddo
            goto 3
          endif
c-- Horizontal  Line
          if (iy1.eq.iy2) then
            iy = iy1
            ixdir = 1
            if (ix2.lt.ix1) ixdir = -1
            do ix = ix1, ix2, ixdir
              if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
                stereo_image(1,ix,iy) = mark_longlat_r
                stereo_image(2,ix,iy) = mark_longlat_g
                stereo_image(3,ix,iy) = mark_longlat_b
              endif
            enddo
            goto 3
          endif
c-- The Rest
          if (abs(ix2-ix1).gt.abs(iy2-iy1)) then
            ixdir = 1
            if (ix2.lt.ix1) ixdir = -1
            do ix = ix1, ix2, ixdir
              iy = nint(float(iy1) + ((float(ix-ix1)/float(ix2-ix1))*float(iy2-iy1)))
              if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
                stereo_image(1,ix,iy) = mark_longlat_r
                stereo_image(2,ix,iy) = mark_longlat_g
                stereo_image(3,ix,iy) = mark_longlat_b
              endif
            enddo
            goto 3
          else
            iydir = 1
            if (iy2.lt.iy1) iydir = -1
            do iy = iy1, iy2, iydir
              ix = nint(float(ix1) + ((float(iy-iy1)/float(iy2-iy1))*float(ix2-ix1)))
              if (1.le.ix.and.ix.le.nx_in.and.1.le.iy.and.iy.le.ny_in) then
                stereo_image(1,ix,iy) = mark_longlat_r
                stereo_image(2,ix,iy) = mark_longlat_g
                stereo_image(3,ix,iy) = mark_longlat_b
              endif
            enddo
            goto 3
          endif
 3        continue
        enddo
      enddo
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
      if (.not.localfast) then
        write (*,'(a)') ' convert -depth 8 -size '//csize//' rgb:'//cmap//'-stereo.rgb -flip '//' '//cmap//'-stereo.png'
        call system('convert -depth 8 -size '//csize//' rgb:'//cmap//'-stereo.rgb -flip '//' '//cmap//'-stereo.png')
c-- Add command and credit lines
        cstring = '-pointsize 20 -fill red -draw '//csq//'text 70,50 '//cdq//command(1:lnblnk(command))//cdq//csq
        write (*,'(a)') ' convert  '//cmap//'-stereo.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo.png'
        call system('convert  '//cmap//'-stereo.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo.png')
        cstring = '-pointsize 20 -fill red -draw '//csq//'text 70,80 '//cdq//'Image credit: Retro Stöckli, NASA Earth Observatory'//cdq//csq
        write (*,'(a)') ' convert  '//cmap//'-stereo.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo.png'
        call system('convert  '//cmap//'-stereo.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo.png')
      endif
c-- Zoom on shaded area (for twitter size image)
      write (cextract,'(I5.5,''x'',I5.5,''+'',I5.5,''+'',I5.5)') ixmax-ixmin+1, iymax-iymin+1, ixmin, iymin
      write (*,'(a)') ' convert -depth 8 -size '//csize//' -extract '//cextract//' rgb:'//cmap//'-stereo.rgb -flip '//' '//cmap//'-stereo-zoom.png'
      call system('convert -depth 8 -size '//csize//' -extract '//cextract//' rgb:'//cmap//'-stereo.rgb -flip '//' '//cmap//'-stereo-zoom.png')
c-- Add command and credit lines
      cstring = '-pointsize 10 -fill red -draw '//csq//'text 35,25 '//cdq//command(1:lnblnk(command))//cdq//csq
      write (*,'(a)') ' convert  '//cmap//'-stereo-zoom.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo-zoom.png'
      call system('convert  '//cmap//'-stereo-zoom.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo-zoom.png')
      cstring = '-pointsize 10 -fill red -draw '//csq//'text 35,40 '//cdq//'Image credit: Retro Stöckli, NASA Earth Observatory'//cdq//csq
      write (*,'(a)') ' convert  '//cmap//'-stereo-zoom.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo-zoom.png'
      call system('convert  '//cmap//'-stereo-zoom.png '//cstring(1:lnblnk(cstring))//' '//cmap//'-stereo-zoom.png')
      return
c
      entry set_map_zoom_mark_colour(mr, mg, mb)
      mark_r = mr
      mark_g = mg
      mark_b = mb
      return
c
      entry set_map_zoom_mark_longlat_colour(mr, mg, mb)
      mark_longlat_r = mr
      mark_longlat_g = mg
      mark_longlat_b = mb
      return
      end

      subroutine gamma_correct(stereo_image, nrgb, nxs, nys, rgamma_in, ggamma_in, bgamma_in, gammalut, ncol, iclim, rcorrect, gcorrect, bcorrect)
      implicit none
c
      real*4    rgamma_in, ggamma_in, bgamma_in, rcorrect, gcorrect, bcorrect
      integer*4 nrgb, nxs, nys, ncol, iclim
      integer*4 gammalut(ncol)
      integer*2 stereo_image(nrgb, nxs, nys)
c
      real*4 r, c, rpos, rgamma, ggamma, bgamma
      integer*4 i, j, k
      rgamma = rgamma_in
      ggamma = ggamma_in
      bgamma = bgamma_in
c-- R
      if (abs(rcorrect).gt.0.0001) then
        do i = 1, nxs
          rpos = abs((float(2*i)/float(nxs)) - 1.0)
          rgamma = rgamma_in + rpos * rcorrect
          do k = 1, iclim + 1
            r = float(k-1) / float(iclim)
            c = r ** (1.0/rgamma)
            gammalut(k) = nint(c*float(iclim))
          enddo
          do j = 1, nys
            if (0.lt.stereo_image(1,i,j).and.stereo_image(1,i,j).le.iclim) stereo_image(1,i,j) = gammalut(stereo_image(1,i,j) + 1)
          enddo
        enddo
      else
        do i = 1, iclim + 1
          r = float(i-1) / float(iclim)
          c = r ** (1.0/rgamma)
          gammalut(i) = nint(c*float(iclim))
        enddo
        do j = 1, nys
          do i = 1, nxs
            if (0.lt.stereo_image(1,i,j).and.stereo_image(1,i,j).le.iclim) stereo_image(1,i,j) = gammalut(stereo_image(1,i,j) + 1)
          enddo
        enddo
      endif
c-- G
      if (abs(gcorrect).gt.0.0001) then
        do i = 1, nxs
          rpos = abs((float(2*i)/float(nxs)) - 1.0)
          ggamma = ggamma_in + rpos * gcorrect
          do k = 1, iclim + 1
            r = float(k-1) / float(iclim)
            c = r ** (1.0/ggamma)
            gammalut(k) = nint(c*float(iclim))
          enddo
          do j = 1, nys
            if (0.lt.stereo_image(2,i,j).and.stereo_image(2,i,j).le.iclim) stereo_image(2,i,j) = gammalut(stereo_image(2,i,j) + 1)
          enddo
        enddo
      else
        do i = 1, iclim + 1
          r = float(i-1) / float(iclim)
          c = r ** (1.0/ggamma)
          gammalut(i) = nint(c*float(iclim))
        enddo
        do j = 1, nys
          do i = 1, nxs
            if (0.lt.stereo_image(2,i,j).and.stereo_image(2,i,j).le.iclim) stereo_image(2,i,j) = gammalut(stereo_image(2,i,j) + 1)
          enddo
        enddo
      endif
c-- B
      if (abs(bcorrect).gt.0.0001) then
        do i = 1, nxs
          rpos = abs((float(2*i)/float(nxs)) - 1.0)
          bgamma = bgamma_in + rpos * bcorrect
          do k = 1, iclim + 1
            r = float(k-1) / float(iclim)
            c = r ** (1.0/bgamma)
            gammalut(k) = nint(c*float(iclim))
          enddo
          do j = 1, nys
            if (0.lt.stereo_image(3,i,j).and.stereo_image(3,i,j).le.iclim) stereo_image(3,i,j) = gammalut(stereo_image(3,i,j) + 1)
          enddo
        enddo
      else
        do i = 1, iclim + 1
          r = float(i-1) / float(iclim)
          c = r ** (1.0/bgamma)
          gammalut(i) = nint(c*float(iclim))
        enddo
        do j = 1, nys
          do i = 1, nxs
            if (0.lt.stereo_image(3,i,j).and.stereo_image(3,i,j).le.iclim) stereo_image(3,i,j) = gammalut(stereo_image(3,i,j) + 1)
          enddo
        enddo
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
      write (*,'(a)') ' convert -depth 8 -size '//csize//' rgb:'//cmap//'-azel.rgb -flip '//' '//cmap//'-azel.png'
      call system('convert -depth 8 -size '//csize//' rgb:'//cmap//'-azel.rgb -flip '//' '//cmap//'-azel.png')
      return
      end

      subroutine map_azel_horizon_init()
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
      entry map_azel_horizon_place(az, el)
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
      entry map_azel_horizon_write(command)
      open (unit=lun, file=cmap//'-azelh.rgb', access='direct', form='unformatted', recl=3*nazel)
      do i = 1, nazel
        do j = 1, nazel
          mapsline(1,j) = maps(1,j,i)
          mapsline(2,j) = maps(2,j,i)
          mapsline(3,j) = maps(3,j,i)
        enddo
        write (lun,rec=i) mapsline
      enddo
      close (unit=lun)
c      call free_lun(lun)
      write (*,'(a)') ' convert -depth 8 -size '//csize//' rgb:'//cmap//'-azelh.rgb -flip '//' '//cmap//'-azelh.png'
      call system('convert -depth 8 -size '//csize//' rgb:'//cmap//'-azelh.rgb -flip '//' '//cmap//'-azelh.png')
      return
      end


      subroutine map_IQ_init()
      implicit none
c--
      integer*4 niq
      parameter (niq=1400)
c--
      integer*1 mapsline (3,niq), maps(3,niq,niq)
      integer*4 lun, i, j, ix, iy
      real*4    torad, todeg, ry, ri, rq, dy, di, dq, dr, dg, db
      character*6 cmap
      character*11 csize
      character*500 cstring
      character*(*) command
c
      save
c
      cmap    = 'satmap' 
      torad   = 2.0 * asin(1.0) / 180.0
      todeg   = 180.0 / 2.0 / asin(1.0)
      write (csize,'(I5.5,''x'',I5.5)') niq, niq
c
      do i = 1, niq
        do j = 1, niq
          maps(1,j,i) = 0
          maps(2,j,i) = 0
          maps(3,j,i) = 0
        enddo
      enddo
c-- Az-El test
      do i = 1, niq/2 - 1
        j  = (niq/2) - nint(float(i) / tan(asin(2.0*float(i)/float(niq))))
        ix = j + 1
        iy = (niq/2) - i
        ix = min(max(ix,1),niq)
        iy = min(max(iy,1),niq)
        dy = 0.5
        di = 2.0 * 0.5957 * float(ix) / float(niq) - 0.5957
        dq = 2.0 * 0.5226 * float(iy) / float(niq) - 0.5226
        dr = (1.000 * dy + 0.956 * di + 0.619 * dq) * 255.0
        dg = (1.000 * dy - 0.272 * di - 0.647 * dq) * 255.0
        db = (1.000 * dy - 1.106 * di + 1.703 * dq) * 255.0
        if (dr.gt.127.0) dr = dr - 256.0
        if (dg.gt.127.0) dg = dg - 256.0
        if (db.gt.127.0) db = db - 256.0
        maps(1,ix,iy) = nint(dr)
        maps(2,ix,iy) = nint(dg)
        maps(3,ix,iy) = nint(db)
        iy = (niq/2) + i
        iy = min(max(iy,1),niq)
        dy = 0.5
        di = 2.0 * 0.5957 * float(ix) / float(niq) - 0.5957
        dq = 2.0 * 0.5226 * float(iy) / float(niq) - 0.5226
        dr = (1.000 * dy + 0.956 * di + 0.619 * dq) * 255.0
        dg = (1.000 * dy - 0.272 * di - 0.647 * dq) * 255.0
        db = (1.000 * dy - 1.106 * di + 1.703 * dq) * 255.0
        if (dr.gt.127.0) dr = dr - 256.0
        if (dg.gt.127.0) dg = dg - 256.0
        if (db.gt.127.0) db = db - 256.0
        maps(1,ix,iy) = nint(dr)
        maps(2,ix,iy) = nint(dg)
        maps(3,ix,iy) = nint(db)
        j  = (niq/2) + nint(float(i) / tan(asin(2.0*float(i)/float(niq))))
        ix = j + 1
        iy = (niq/2) - i
        ix = min(max(ix,1),niq)
        iy = min(max(iy,1),niq)
        dy = 0.5
        di = 2.0 * 0.5957 * float(ix) / float(niq) - 0.5957
        dq = 2.0 * 0.5226 * float(iy) / float(niq) - 0.5226
        dr = (1.000 * dy + 0.956 * di + 0.619 * dq) * 255.0
        dg = (1.000 * dy - 0.272 * di - 0.647 * dq) * 255.0
        db = (1.000 * dy - 1.106 * di + 1.703 * dq) * 255.0
        if (dr.gt.127.0) dr = dr - 256.0
        if (dg.gt.127.0) dg = dg - 256.0
        if (db.gt.127.0) db = db - 256.0
        maps(1,ix,iy) = nint(dr)
        maps(2,ix,iy) = nint(dg)
        maps(3,ix,iy) = nint(db)
        iy = (niq/2) + i
        iy = min(max(iy,1),niq)
        dy = 0.5
        di = 2.0 * 0.5957 * float(ix) / float(niq) - 0.5957
        dq = 2.0 * 0.5226 * float(iy) / float(niq) - 0.5226
        dr = (1.000 * dy + 0.956 * di + 0.619 * dq) * 255.0
        dg = (1.000 * dy - 0.272 * di - 0.647 * dq) * 255.0
        db = (1.000 * dy - 1.106 * di + 1.703 * dq) * 255.0
        if (dr.gt.127.0) dr = dr - 256.0
        if (dg.gt.127.0) dg = dg - 256.0
        if (db.gt.127.0) db = db - 256.0
        maps(1,ix,iy) = nint(dr)
        maps(2,ix,iy) = nint(dg)
        maps(3,ix,iy) = nint(db)
      enddo
      do i = 1, niq
        iy = i
        ix = niq / 2
        ix = min(max(ix,1),niq)
        iy = min(max(iy,1),niq)
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
      enddo
      do i = 1, niq
        ix = i
        iy = niq / 2
        ix = min(max(ix,1),niq)
        iy = min(max(iy,1),niq)
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
      enddo
      do i = 1, 720
        ix = (niq/2) + int(dsin(dble(i/2)*torad)*dble(niq/6))
c iy is counterintuitive as the image is vertically flipped for display !
        iy = (niq/2) + int(dcos(dble(i/2)*torad)*dble(niq/6))
        ix = min(max(ix,1),niq)
        iy = min(max(iy,1),niq)
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
      enddo
      do i = 1, 720
        ix = (niq/2) + int(dsin(dble(i/2)*torad)*dble(niq/3))
c iy is counterintuitive as the image is vertically flipped for display !
        iy = (niq/2) + int(dcos(dble(i/2)*torad)*dble(niq/3))
        ix = min(max(ix,1),niq)
        iy = min(max(iy,1),niq)
        maps(1,ix,iy) = -1
        maps(2,ix,iy) = -1
        maps(3,ix,iy) = -1
      enddo
      return
c
      entry map_IQ_place(ry, ri, rq)
      ix = (niq/2) + int(float(niq)*ri/1024.0)
      iy = (niq/2) + int(float(niq)*rq/1024.0)
      ix = min(max(ix,1),niq)
      iy = min(max(iy,1),niq)
      dr = 1.000 * ry + 0.956 * ri + 0.621 * rq
      dg = 1.000 * ry - 0.272 * ri - 0.647 * rq
      db = 1.000 * ry - 1.106 * ri + 1.703 * rq
c      write (*,*) ri, rq, ix, iy, dr, dg, db
c      read (*,*)
      maps(1,ix,iy) = int(dr/8.0)
      maps(2,ix,iy) = int(dg/8.0)
      maps(3,ix,iy) = int(db/8.0)
      return
c
      entry map_IQ_write(command)
      open (unit=lun, file=cmap//'-IQ.rgb', access='direct', form='unformatted', recl=3*niq)
      do i = 1, niq
        do j = 1,niq
          mapsline(1,j) = maps(1,j,i)
          mapsline(2,j) = maps(2,j,i)
          mapsline(3,j) = maps(3,j,i)
        enddo
        write (lun,rec=i) mapsline
      enddo
      close (unit=lun)
      call free_lun(lun)
      write (*,'(a)') ' convert -depth 8 -size '//csize//' rgb:'//cmap//'-IQ.rgb -flip '//' '//cmap//'-IQ.png'
      call system('convert -depth 8 -size '//csize//' rgb:'//cmap//'-IQ.rgb -flip '//' '//cmap//'-IQ.png')
      return
      end

      subroutine dump_merge(i_merge, satch_image, satch_valid, longlat, nmaxpix, nmaxch, nmaxlin, lundump, nmergemax)
      implicit none
      integer*4 nmaxpix, nmaxch, nmaxlin, i_merge, nmergemax
      integer*1 satch_valid(nmaxlin)
      integer*2 satch_image(nmaxpix,nmaxch,nmaxlin)
      integer*4 lundump(3, nmergemax)
      real*8    longlat(3,nmaxpix+1,nmaxlin)
c
      integer*4 lunimg, lunval, lunll, i
      character*2 cmerge
      character*250 outstring
c
      save
c
      call get_lun(lundump(1,i_merge))
      call get_lun(lundump(2,i_merge))
      call get_lun(lundump(3,i_merge))
      call getenv('HRPTOUT',outstring)
      write (cmerge,'(I2.2)') i_merge
      open (unit=lundump(1,i_merge),file=outstring(1:lnblnk(outstring))//'img'//cmerge//'.dat',form='unformatted', access='direct', recl=nmaxpix*nmaxch*2)
      open (unit=lundump(2,i_merge),file=outstring(1:lnblnk(outstring))//'val'//cmerge//'.dat',form='unformatted', access='direct', recl=nmaxlin)
      open (unit=lundump(3,i_merge),file=outstring(1:lnblnk(outstring))//'ll_'//cmerge//'.dat',form='unformatted', access='direct', recl=(nmaxpix+1)*3*8)
      write (lundump(2,i_merge),rec=1) satch_valid
      do i = 1, nmaxlin
        call write_i2_buffer(lundump(1,i_merge),i,satch_image(1,1,i),nmaxpix*nmaxch)
        call write_r8_buffer(lundump(3,i_merge),i,longlat(1,1,i)    ,(nmaxpix+1)*3)
      enddo
      return
c
      entry load_merge(i_merge, satch_image, satch_valid, longlat, nmaxpix, nmaxch, nmaxlin, lundump, nmergemax)
      read (lundump(2,i_merge),rec=1) satch_valid
      do i = 1, nmaxlin
        call read_i2_buffer(lundump(1,i_merge),i,satch_image(1,1,i),nmaxpix*nmaxch)
        call read_r8_buffer(lundump(3,i_merge),i,longlat(1,1,i)    ,(nmaxpix+1)*3)
      enddo
      close (unit=lundump(1,i_merge),status='delete')
      close (unit=lundump(2,i_merge),status='delete')
      close (unit=lundump(3,i_merge),status='delete')
      call free_lun(lundump(1,i_merge))
      call free_lun(lundump(2,i_merge))
      call free_lun(lundump(3,i_merge))
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


      subroutine grid_sample( v1, v2, v3, v4, v5, v6, v7, v8, v9, ll1, ll2, ll3, ll4, val_img, val_ll, ndiv)
      implicit none
      integer*4 ndiv
      integer*2 v1, v2, v3, v4, v5, v6, v7, v8, v9, val_img(ndiv, ndiv), vy1, vy2, vy3
      real*8    ll1(2), ll2(2), ll3(2), ll4(2), val_ll(2, ndiv, ndiv), lly1, lly2
c
c The calling order of the corner points is 
c
c      ===> Scan direction
c
c     ll1--------ll2         |
c      |          |          |    Track direction - so for North ll3(2) > ll1(2) 
c      |          |          |
c      |          |          \/
c     ll3--------ll4
c
c-- There needs to be a check to see if the order of the latitudes is OK - this to prevent bow-tie correction interpolation errors for MODIS (see main code)
c
      integer*4 i, j
      real*4    wr1, wr2
c
      do i = 1, ndiv
        do j = 1,ndiv
          wr1 = abs(float(j-(ndiv/2)))/float(ndiv)
          wr2 = 1. - wr1
          if (j.le.(ndiv/2)) then
            vy1 = int(wr2 * float(v2) + wr1 * float(v1))
            vy2 = int(wr2 * float(v5) + wr1 * float(v4))
            vy3 = int(wr2 * float(v8) + wr1 * float(v7))
          else
            vy1 = int(wr2 * float(v2) + wr1 * float(v3))
            vy2 = int(wr2 * float(v5) + wr1 * float(v6))
            vy3 = int(wr2 * float(v8) + wr1 * float(v9))
          endif
          wr1 = abs(float(i-(ndiv/2)))/float(ndiv)
          wr2 = 1. - wr1
          if (i.le.(ndiv/2)) then
            val_img(j,i) = int(wr2 * float(vy2) + wr1 * float(vy1))
          else
            val_img(j,i) = int(wr2 * float(vy2) + wr1 * float(vy3))
          endif
        enddo
      enddo
      do i = 1, ndiv
        do j = 1, ndiv
          lly1 = ll1(1) + dble(j-1) * (ll2(1) - ll1(1)) / dble(ndiv)
          lly2 = ll3(1) + dble(j-1) * (ll4(1) - ll3(1)) / dble(ndiv)
          val_ll(1,j,i) = lly1 + dble(i-1) * (lly2 - lly1) / dble(ndiv)
          lly1 = ll1(2) + dble(j-1) * (ll2(2) - ll1(2)) / dble(ndiv)
          lly2 = ll3(2) + dble(j-1) * (ll4(2) - ll3(2)) / dble(ndiv)
          val_ll(2,j,i) = lly1 + dble(i-1) * (lly2 - lly1) / dble(ndiv)
        enddo
      enddo
      return
      end
      
      subroutine darkness_init()
      implicit none
      integer*4 linedoy
      real*8    linetimess, longss, latss, solzenang_out
c
      integer*4 dattim(8), monthdays(12), diy, doy, i, timezone
      real*4 days_in_year, pi, gamma, eqtime, decl, time_offset, tst, ha, long, lat, ha_rise, ha_set, sunrise, sunset, solzenang, rsc, tstval
      real*8 torad, todeg, linetime, linetimess_save
      character*5  c_zone
      character*8  c_date
      character*10 c_time
c
      data monthdays/31,0,31,30,31,30,31,31,30,31,30,31/
c
      save
c
      pi              = 2.0 * asin(1.0)
      torad           = 2.0D0 * dasin(1.0D0) / 180.0D0
      todeg           = 180.0D0 / (2.0D0 * dasin(1.0D0))
      linetimess_save = -1.0
c
      call date_and_time(c_date, c_time, c_zone, dattim)
      timezone = nint(float(dattim(4))/60.0)
c
      diy = 365
      if (mod(dattim(1),4).eq.0.and.dattim(1).ne.2000) diy = diy + 1
      days_in_year    = float(diy)
      return
c
      entry darkness(longss, latss, linetimess, linedoy, solzenang_out)
      long       = longss * torad
      lat        = latss  * torad
      doy        = linedoy
c- Satellite time is in UTC - convert to local (i.e. compensate for timezone or remove timezone from the code below)
      if (linetimess.ne.linetimess_save) then
        linetimess_save = linetimess
        linetime   = linetimess + dble(timezone) * 3600.0D0
        dattim(5)  = int (linetime / 3600.0D0)
        dattim(6)  = int((linetime - dble(dattim(5))*3600.0D0)/60.0D0)
        rsc        = linetime - dble(dattim(5)) * 3600.0D0 - dble(dattim(6)) * 60.0D0
c
        gamma      = ((2.0 * pi) / days_in_year) * (float(doy) - 1.0 + ((float(dattim(5))-12.0)/24.0) + float(dattim(6))/1440.)
c
        eqtime     = 229.18*(0.000075 + 0.001868*cos(gamma) - 0.032077*sin(gamma) - 0.014615*cos(2.0*gamma) - 0.040849*sin(2.0*gamma))
        decl       = 0.006918 - 0.399912*cos(gamma) + 0.070257*sin(gamma) - 0.006758*cos(2.0*gamma) + 0.000907*sin(2.0*gamma) - 0.002697*cos(3.0*gamma) + 0.00148*sin(3.0*gamma)
        tstval     = float(dattim(5)) * 60.0 + float(dattim(6)) + (rsc / 60.0)
      endif
c
      time_offset= eqtime + 4.0 * longss - 60.0 * float(timezone)
      tst        = tstval + time_offset
      ha         = ((tst / 4.0) - 180.0) * torad
      solzenang  = acos(sin(lat) * sin(decl) + cos(lat) * cos(decl) * cos(ha))
      solzenang  = solzenang * todeg
      solzenang_out = solzenang
c
      return
      end

      subroutine get_borders(stereo_image,nc,nx_in,ny_in,swborderhighres, i_highres)
      implicit none
      integer*4 nx_in, ny_in, nc, i_highres
      integer*2 stereo_image(nc,nx_in,ny_in)
      logical   swborderhighres
c
c      integer*1 buf(1048576)
      integer*1 buf(30000000)
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
      if (swborderhighres) then
        if (i_highres.eq.1) call get_shp_borders('./resource/border-L1-highres.shp',buf,30000000,n_actual,100)
        if (i_highres.eq.2) call get_shp_borders('./resource/border-L2-highres.shp',buf,30000000,n_actual,100)
        if (i_highres.eq.3) call get_shp_borders('./resource/shoreline-L1-highres.shp',buf,30000000,n_actual,100)
        if (i_highres.eq.4) call get_shp_borders('./resource/shoreline-L2-highres.shp',buf,30000000,n_actual,100)
      else
        call get_shp_borders('./resource/borders.shp',buf,30000000,n_actual,100)
      endif
      call move_bytes(buf(25),i1_fl(1),4)
      call move_bytes(buf(33),i1_shape(1),4)
      do while (.true.)
        if (swborderhighres) then
          if (i_highres.eq.1) call get_shp_borders('./resource/border-L1-highres.shp',buf,30000000,n_actual,52)
          if (i_highres.eq.2) call get_shp_borders('./resource/border-L2-highres.shp',buf,30000000,n_actual,52)
          if (i_highres.eq.3) call get_shp_borders('./resource/shoreline-L1-highres.shp',buf,30000000,n_actual,52)
          if (i_highres.eq.4) call get_shp_borders('./resource/shoreline-L2-highres.shp',buf,30000000,n_actual,52)
        else
          call get_shp_borders('./resource/borders.shp',buf,30000000,n_actual,52)
        endif
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
        if (swborderhighres) then
          if (i_highres.eq.1) call get_shp_borders('./resource/border-L1-highres.shp',buf,30000000,n_actual,i4_length*2+8-52)
          if (i_highres.eq.2) call get_shp_borders('./resource/border-L2-highres.shp',buf,30000000,n_actual,i4_length*2+8-52)
          if (i_highres.eq.3) call get_shp_borders('./resource/shoreline-L1-highres.shp',buf,30000000,n_actual,i4_length*2+8-52)
          if (i_highres.eq.4) call get_shp_borders('./resource/shoreline-L2-highres.shp',buf,30000000,n_actual,i4_length*2+8-52)
        else
          call get_shp_borders('./resource/borders.shp',buf,30000000,n_actual,i4_length*2+8-52)
        endif
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
c      save iptr, buffer, irec, lun, nbytes, eof
       save
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
c-- There is a difference here between my system at work and the system at home .......
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

      subroutine write_i2_buffer(lun,irec,array,n)
      implicit none
      integer*4 n, irec, lun
      integer*2 array(n)
c
      write (lun,rec=irec) array
      return
      end

      subroutine write_r8_buffer(lun,irec,array,n)
      implicit none
      integer*4 n, irec, lun
      real*8 array(n)
c
      write (lun,rec=irec) array
      return
      end

      subroutine read_i2_buffer(lun,irec,array,n)
      implicit none
      integer*4 n, irec, lun
      integer*2 array(n)
c
      read (lun,rec=irec) array
      return
      end

      subroutine read_r8_buffer(lun,irec,array,n)
      implicit none
      integer*4 n, irec, lun
      real*8 array(n)
c
      read (lun,rec=irec) array
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

      integer*2 function fetch_pixel(array, n, i, j)
      implicit none
      integer*4 n, i, j
      integer*2 array(n, n)
c
      fetch_pixel = array(i, j)
c
      return
      end
                  
      real*8 function fetch_ll(array, n1, n2, i, j, k)
      implicit none
      integer*4 n1, n2, i, j, k
      real*8    array(n1, n2, n2)
c
      fetch_ll = array(i, j, k)
c
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

      subroutine string_to_r4(string,npos,value)
      logical*1 swexp
      integer*2 npos,ilen,i_mantisse,i_main
      real*4    value
      character*1 cf
      character*10 cfmt
      character*(*) string
      ilen = len(string)
      do 1 i = 1,ilen
        if (string(i:i).ne.' ') then
          i1 = i
          do 2 j = i,ilen
            if (string(j:j).eq.' ') then
              i2 = j - 1
              goto 3
            endif
 2        continue
          i2 = ilen
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
        write (cfmt,'(''('',A1,I2,''.'',I1,'')'')') 
     *         cf,i_main,i_mantisse
      else
        write (cfmt,'(''('',A1,I1,''.'',I1,'')'')') 
     *         cf,i_main,i_mantisse
      endif
      read (string(i1:i2),cfmt) value
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

c-- Convert days since 1-Jan-2000 to day in current year
      subroutine dy2k_to_doy(dy2k, doy)
      implicit none
      integer*4 dy2k, doy
c
      integer*4 i_day, i_yr
c
      i_day = 0
      i_yr  = 2000
      do while (i_day.lt.dy2k)
        if (mod(i_yr,4).eq.0.and.i_yr.ne.2000) then
          i_day = i_day + 366
        else
          i_day = i_day + 365
        endif
        i_yr = i_yr + 1
      enddo
      i_yr = i_yr - 1
      if (mod(i_yr,4).eq.0.and.i_yr.ne.2000) then
        i_day = i_day - 366
      else
        i_day = i_day - 365
      endif
      doy = dy2k - i_day
      return
      end

      subroutine hist_plot(hist,nx,ny,nxin,nyin,filenm)
      implicit none
      integer*4 nx, ny, hist(nx, ny), nxin, nyin
      character*(*) filenm
c
      integer*4 lunplot, i, j
      real*4 xr(2), yr(2), rhisty(4096), rhistx(4096)
      character*79 ctxt
      call get_lun(lunplot)
      call PS_init_colourtable(1,'./resource/color')
      call set_PS_fullpage()
      call pkg_openpl(filenm(1:lnblnk(filenm))//'.ps', lunplot)
      call newpen(6)
      xr(1) = 0.
      xr(2) = float(nxin)
      yr(1) = hist(1,1)
      yr(2) = yr(1)
      do i = 1,nxin
        do j = 1,nyin
          if (hist(i,j).lt.yr(1)) yr(1) = hist(i,j)
          if (hist(i,j).gt.yr(2)) yr(2) = hist(i,j)
        enddo
      enddo
      call pkg_frame(11,-6,1.,xr,yr,'Bin #','Frequency','Channel histogram')
      do j = 1,nyin
        do i = 1,nxin
          rhisty(i) = float(hist(i,j))
          rhistx(i) = float(i)
        enddo
        call newpen(7-j)
        call pkg_plhist(rhistx,rhisty,nxin)
      enddo
      do j = 1,nyin
        write (ctxt,'(''Channel '',i1)') j
        call newpen(7-j)
        call pkg_pltextbl(ctxt(1:lnblnk(ctxt)),j)
      enddo
      call pkg_clospl()
      call free_lun(lunplot)
      return
      end

      subroutine sza_plot(szahist, n, filenm)
      implicit none
      integer*4 n, szahist(n)
      character*(*) filenm
c
      integer*4 lunplot, i, j
      real*4 xr(2), yr(2), rhisty(1800), rhistx(1800), sum
      character*79 ctxt
      call get_lun(lunplot)
      call PS_init_colourtable(1,'./resource/color')
      call set_PS_fullpage()
      call pkg_openpl(filenm(1:lnblnk(filenm))//'.ps', lunplot)
      call newpen(6)
      xr(1) = -1.0
      xr(2) = -1.0
      do i = 1, n
        if (szahist(i).ne.0.and.xr(1).eq.-1.0)     xr(1) = float(i-1) * 0.1
        if (szahist(n+1-i).ne.0.and.xr(2).eq.-1.0) xr(2) = float(n-i) * 0.1
      enddo
      sum = 0.
      do i = 1,n
        sum = sum + float(szahist(i))
        rhisty(i) = sum
        rhistx(i) = float(i-1) * 0.1
      enddo
      do i = 1,n
        rhisty(i) = rhisty(i) / sum
      enddo
      yr(1) = 0.
      yr(2) = 1.
      call pkg_frame(11,11,1.,xr,yr,'SZA','Cumul. Frequency','SZA histogram')
      call pkg_plhist(rhistx,rhisty,n)
      call newpen(2)
      call pkg_raster(11,11,xr,yr,-2)
      call pkg_clospl()
      call free_lun(lunplot)
      return
      end
      
      subroutine save_string(instring)
      implicit none
      character*(*) instring, outstring
c
      integer*4 i
      character*1000 savestring
c
      save savestring
c
      if (lnblnk(instring).gt.len(savestring)) then
        do i = 1,len(savestring)
          savestring(i:i) = ' '
        enddo
        savestring(1:len(savestring)) = instring(1:len(savestring))
      else
        do i = 1,len(savestring)
          savestring(i:i) = ' '
        enddo
        savestring(1:lnblnk(instring)) = instring(1:lnblnk(instring))
      endif
      return
c
      entry get_string(outstring)
      if (lnblnk(savestring).gt.len(outstring)) then
        do i = 1,len(outstring)
          outstring(i:i) = ' '
        enddo
        outstring(1:len(outstring)) = savestring(1:len(outstring))
      else
        do i = 1,len(outstring)
          outstring(i:i) = ' '
        enddo
        outstring(1:lnblnk(savestring)) = savestring(1:lnblnk(savestring))
      endif
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

      subroutine get_help()
      implicit none
c
      write (*,'(a)') 'Program usage:'
      write (*,'(a)') './hrpt.exe filename xyz'
      write (*,'(a)') '           filename is one of the FY, Metop or NOAA images in full - can have a directory in front of it, but it does need ./ at the beginning of the file name at the very least.'
      write (*,'(a)') '             The filename algorithm to separate directory and filename seeks for / backwards from the end of the string'
      write (*,'(a)') '           If filename starts with @ , the rest of the string is considered to be a directory (end has to be /) containing a batch set of data files to be processed'
      write (*,'(a)') '           If merge is used as an option, the filename is the file containing the files to be merged - and optional switches like theta'
      write (*,'(a)') '           This batch processing uses an average Solar Zenith Angle of 82 degrees to determine if a switch to the IR channels is needed'
      write (*,'(a)') ''
      write (*,'(a)') '           North or South is determined from the TLE calculations'
      write (*,'(a)') ''
      write (*,'(a)') '           xyz is a 3 digit code for generating the merged ch11 images with x for R, y for G and z for B - e.g. 221 is @petermeteor  '
      write (*,'(a)') '             standard (and is the default, except for FY where default is 621) - default is invoked e.g. when giving xxx'
      write (*,'(a)') ''
      write (*,'(a)') 'Parameters 3-25 are optional and have no specific order:'
      write (*,'(a)') ''
      write (*,'(a)') '  gapcor        ==> analyses the time codes of the lines and inserts gaps in the images for the missing lines - needed for border overplot'
      write (*,'(a)') '  debug         ==> Nr of diagnostic print outs - use >& hrpt.log to route to file'
      write (*,'(a)') '  project       ==> project the satellite image on the zoomed swath area :options: project=4 ==> ch 4 as channel in B/W image. project=221 use ch2 for R, ch 2 for G and ch 1 for B in projected image. '
      write (*,'(a)') '                    Only in case of FY satellites 0 => ch 10 ! Default w. project is project=221 - requires gapcor'        
      write (*,'(a)') '  border        ==> overplot country borders on projected satellite image - requires project (and project requires gapcor)'
      write (*,'(a)') '  gamma=x.y     ==> gamma correction for corrected image. Example : gamma=1.0 - default gamma=1.4'
      write (*,'(a)') '  theta=x.y     ==> Override the thetascan value below for testing purposes (x.y is free format)'
      write (*,'(a)') '                    theta=0.0 is a special option with a varying angle over latitude'
      write (*,'(a)') '                    theta=-2.91,1.89,0.2 means loop over theta starting from -2.91 up to 1.89 in 0.2 deg steps - avoid (!!!) theta becoming lower than 0.001 - this will trigger the special correction'
      write (*,'(a)') '                    Theta defines the rotation angle of the velocity vector around the Nadir vector (Nadir vector being the vector from the satellite to the sub-satellite point)'
      write (*,'(a)') '  eta=x.y       ==> Define the roration of the Nadir vector around the velocity vector (for now seen on MN2 only) default for MN2 : eta=2.4'
      write (*,'(a)') '  szalim=x.y    ==> Activate Solar Zenith Angle calculation (x.y is free format) - requires project (and project requires gapcor)'
      write (*,'(a)') '  szamerge=n    ==> Activate merging of IR channel with  the project=xyz image. Merge with channel n (1-digit - a or A is for the FY ch 10). Merge from szalim -0.5 to + 0.5 beyond that scaled IR'
      write (*,'(a)') '  szafact=x.y   ==> Normalisation of IR channel beyond szalim + 0.5 autoscaling. (default is 1.2) (x.y is free format)'
      write (*,'(a)') '  linedelay=x.y ==> # of lines of delay introduced on the linetime (can be negative) (x.y is free format but only int(x.y) is internally used) - for HRPT 166 msec per line, for MODIS see inside the code'
      write (*,'(a)') '  hybrid        ==> create the ch12, ch13 etc multi channel images using the GOES LUT'
      write (*,'(a)') '  fast          ==> minimise the output to what is really needed (no ch images, no worldmap no stereo full norther hemisphere, no PS plots etc)'
      write (*,'(a)') '  YIQ           ==> Experimental YIQ extensions'
      write (*,'(a)') '  horizon       ==> In batch mode, overplot AzEl plots to create a visibility mask'
      write (*,'(a)') '  sharpen       ==> sharpens the polar stereographic projected image (dont use with border or longlat enabled as it causes subtle border and longlat colour changes) - default : sharpen=6x3+0.5+0'
      write (*,'(a)') '  truech        ==> Different FY3 - truecolour approach, cf to default filter median wavelength translated to RGB and weighted - with now R = ch1 G = (ch8 + ch9) / 2  B = ch7'
      write (*,'(a)') '  binary        ==> The files (FY3 or Metop) to be analysed are in .bin format (direct XHRPT output)'
      write (*,'(a)') '  opacity       ==> Allow the introduction of a 2nd RGB to be merged at a certain opacity. Optional parameters: opacity=rgb,perc, like opacity=197,40 meaning merge 197 as RGB at 40% level'
      write (*,'(a)') '                    with the project= defined rgb image. Default is opacity=197,25. Just opacity=197 would add 197 at the default 25% level'
      write (*,'(a)') '                    Warning - this is largely only seful for the 10 ch FY3 data !!! Also, the >5 ch scanline offsets (see inside code) cause colour misery towards the edges'
      write (*,'(a)') '                    NB : Opacity is ignored for batch modes where SZA > 82.0 and the colour choices get changed to night mode'
      write (*,'(a)') '  monitor       ==> Intended for use with batch mode and binary to monitor the XHRPT output directory and auto-trigger the processing of .bin and .raw16 - terminate with a touch -f killhrpt.bin in the output directory'
      write (*,'(a)') '  histcor       ==> For the polar projection image, normalise the R, G and B histograms to their 0.5% to 99.5% range'
      write (*,'(a)') '  longlat       ==> Use longlat or longlat=n,m with n and m integer long and lat grid distances default: longlat=10,10'
      write (*,'(a)') '  bcol=         ==> Use: bcol=r,g,b  with r, g and b as 10 bit border color values (0-1023) default: bcol=1023,1023,0 - for Modis the 1023 changes to the 12 bit value 4095'
      write (*,'(a)') '  llcol=        ==> Use: llcol=r,g,b  with r, g and b as 10 bit longlat color values (0-1023) default: llcol=1012,40,692 - for Modis this changes to the 12 bit equivalent'
      write (*,'(a)') '  modis         ==> Switch the internal limits etc to MODIS values (like 12 bit colours) and use modis 1km data'
      write (*,'(a)') '                    modis250 or (!) modis500 as a parameter switch to the use of the 250 or 500m data respectively'
      write (*,'(a)') '                    If you add =correct to the end of one of the modis switches, it engages readbin_modis with a switch "correct" - this tries a linescan raster correction'
      write (*,'(a)') '                    REMEMBER ! for using 250 m data you need to have set export MMODEL=LARGE and have done ./buildhrpt.sh !!!'
      write (*,'(a)') '  reflong=x.y   ==> Rotate the polster.rgb to be at your reference longitude - takes about 1 minute extra. If parameter not present then polster_orig.rgb is copied to polster.rgb. reflong=90.0 for example is US centred'
      write (*,'(a)') '  merge         ==> A complex option to merge multiple images from different passes and different satellites. It also is required by the modispan option for Aqua MODIS truecolour processing'
      write (*,'(a)') '                ==> The first - filename - argument of the command line contains the files to be merged - or in case of panmodis, 3 (!!!! YES 3) times the same filename'
      write (*,'(a)') '                ==> In this text file after the filename you can specify non-standard theta and ilinedelay value - but if used once it needs (!!!!) to be present on all the lines'
      write (*,'(a)') '  modispan      ==> Used with the merge option to create truecolour modis images from (hardcoded RGB = ch 1, 4, 3 - so a mix of 250 m and 500 m data)'
      write (*,'(a)') '                ==> Example : ./hrpt.exe panmodis.txt 221 gapcor project=121 gammargb=0.92,0.98,0.92 fast histcor binary merge modispan=correct'
      write (*,'(a)') '  gammargb=x,y,z => Used to specify separate gammas for R, G and B example - gammargb=0.92,0.98,0.92 - only applied to the projected or the border image - note : it changes (slightly) the border and longlat colours'
      write (*,'(a)') ''
      write (*,'(a)') 'Inherited by other s/w'
      write (*,'(a)') '  bowtie=x.y,x.y ==> Passed on to readbin_modis.exe as part of the command string to set the bowtie correction parameters - default : bowtie=0.18,0.42 defined inside readbin_bowtie.f'
      write (*,'(a)') ''
      write (*,'(a)') 'Auxiliary files needed (in ./resource) :'
      write (*,'(a)') ''
      write (*,'(a)') '  LUT.rgb - a GOES lut for creating multicolor images'
      write (*,'(a)') '  worldmaphr.rgb - the MODIS 5400 by 2700 pixel map (credit: Retro Stöckli, NASA Earth Observatory) from https://visibleearth.nasa.gov/view.php?id=73580'
      write (*,'(a)') '  worldmapuhr.rgb - the MODIS 21600 by 10800 pixel map converted to polar stereographic (credit: Reto Stöckli, NASA Earth Observatory) from https://visibleearth.nasa.gov/view.php?id=73580'
      write (*,'(a)') '    used by polarstereo.f to create polster.rgb - if longref=x.y is not specified, polster_orig.rgb is copied to polster.rgb (File for Greenwich Meridian)'
      write (*,'(a)') '  polster_orig.rgb - see previous file'
      write (*,'(a)') '  color1.tab - color file for the plot software to generate the histograms'
      write (*,'(a)') '  color2.tab - color file for the plot software to generate the histograms'
      write (*,'(a)') '  color3.tab - color file for the plot software to generate the histograms'
      write (*,'(a)') '  color4.tab - color file for the plot software to generate the histograms'
      write (*,'(a)') '  borders.shp - low res borders file (For HRPT satellites)'
      write (*,'(a)') '  TLE files for the day of the dass, the day before or the day after - they arre searched for in directory the TLE environment variable points to'
      write (*,'(a)') ' '
      write (*,'(a)') 'For highres borders (for MODIS) :'
      write (*,'(a)') '  border-L1-highres.shp'
      write (*,'(a)') '  border-L2-highres.shp'
      write (*,'(a)') '  shoreline-L1-highres.shp'
      write (*,'(a)') '  shoreline-L2-highres.shp'
      write (*,'(a)') ''
      write (*,'(a)') 'Several exe files are needed as they possibly launch from hrpt.exe'
      write (*,'(a)') '  readbin.exe - for converting the .bin files from Joe s XHRPT directly to my format for the program (dont need Logs s/w - still has an occasional bug for FY)'
      write (*,'(a)') '  readbin_modis.exe - for converting the MODIS .bin files from Jean-Lucs format (deframed TM) to my format'
      write (*,'(a)') '  polarstereo.exe - for rotating the polarstereo.rgb to a different longitude'
      write (*,'(a)') ''
      write (*,'(a)') 'Several files are used in the compilation of the legacy SGP4 code (in ./src/) :'
      write (*,'(a)') '  astmath.cmn'
      write (*,'(a)') '  ASTMATH.CMN'
      write (*,'(a)') '  sgp4.cmn'
      write (*,'(a)') ''
      write (*,'(a)') 'Output files produced (in the directory defined by HRPTOUT) :'
      write (*,'(a)') ''
      write (*,'(a)') '  Per channel B/W images : for NOAA/Metop ch1.png through ch5.png 16 bit images using the full 10 bit dynamic range'
      write (*,'(a)') '                           for FengYun ch1.png through chA.png 16 bit images using the full 10 bit dynamic range'
      write (*,'(a)') '                           foe MODIS 250 m the two images for ch1 and ch 2. For 500 m Ch 3, 4, 5, 6 and 7. For 1 km the 5 channels (could be 10) selected in the s/w !!'
      write (*,'(a)') '  Largely deprecated : Two channel colour images using the GOES LUT (might suppress this in future)'
      write (*,'(a)') '  For FengYun only "true" colour images based on the four visible channel scans translated center wavelength to RGB'
      write (*,'(a)') '    and merged into one image - need to determine the weight'
      write (*,'(a)') '  RGB image ch11.png using for RGB the channels from the xyz command line parameter'
      write (*,'(a)') '  RGB image ch11_corrected.png using for RGB the channels from the xyz command line parameter and corrected to a fixed pixel size on the earths surface'
      write (*,'(a)') '  satmap.png image showing the image swath grayed out with a purple subsatellite line over the modis world map'
      write (*,'(a)') '  satmap_stereo-zoom.png - the greyed out area zoomed in on the worldmap'
      write (*,'(a)') '  Dependent on the parameters border and/or project - a ...border.png or ...project.png image - these can have nerged IR channels using szalim, szamerge and szafactor'
      write (*,'(a)') ''
      write (*,'(a)') 'New developments to be done:'
      write (*,'(a)') '  merging several passes into one big image (currently done in hrptmerge, but this is limited to a few file formats (.C10, .hpt and .raw16))'
      write (*,'(a)') ''
      write (*,'(a)') 'And now for the code:'
      write (*,'(a)') 'Top section - until get_lun and free_lun - all written fully by me using wikipedia, NOAA and Eumetsat public documentation'
      write (*,'(a)') 'Plot code - some 50% written in the early 80s by Hans Deutekom at the Cosmic Ray Working Group in Leiden - heavily modified by me for postscript'
      write (*,'(a)') '            the rest all written by me (pltextbl, plraster, plhist, plpois etc etc) - now in plotsoft.f (as -O2 as a compiler switch ruined the output)'
      write (*,'(a)') 'SGP4 code from the 70s written by David Vallado (see T.S Kelso celectrak.com) in ancient style fortran and where needed adapted by me.'
      write (*,'(a)') 'xyz2lla and lla2xyz taken from https://ea4eoz.blogspot.com/2015/11/simple-wgs-84-ecef-conversion-functions.html'
      write (*,'(a)') '            with an implemented improvement from the footnote'
      write (*,'(a)') ''
      write (*,'(a)') 'The s/w REQUIRES (!!!!) daily download TLE files - I do this on my NAS with embedded Linux and the CROTAB script : '
      write (*,'(a)') ''
      write (*,'(a)') '#!/bin/bash'
      write (*,'(a)') 'mytlefile=/share/Public/TLE/`(date +%Y%m%d)`-weather.tle'
      write (*,'(a)') 'wget --secure-protocol=auto https://www.celestrak.com/NORAD/elements/weather.txt'
      write (*,'(a)') 'mv weather.txt $mytlefile'
      write (*,'(a)') ''
      write (*,'(a)') 'NOTE !!!! : The TLE files are automatically selected from the directory pointed to by the TLE env variable. Expected format yyyymmdd-weather.tle'
      write (*,'(a)') ''
      write (*,'(a)') ''
      write (*,'(a)') 'Environment variables used :'
      write (*,'(a)') 'export TLE=/home/fjansen/psa/TLE/               ==> Points to the TLE repository. To analyse a file a TLE on the day itself, 1 before or 1 after is mandatory !!'
      write (*,'(a)') 'export HRPTOUT=/home/fjansen/Win7/SDR/Output/   ==> Points to the where the plots, images and data files are to be stored'
      write (*,'(a)') 'export MMODEL=LARGE  (for the HRPT, MODIS 1km and 500m and also doing MODIS 250 m data) or export MMODEL=MEDIUM for the HRPT, MODIS 1km and 500m'
      write (*,'(a)') ''
      write (*,'(a)') 'The s/w is built by the script ./buildhrpt.sh'
      write (*,'(a)') ''
      write (*,'(a)') 'Theta is a rotation of the scanline around the Nadir vector (satellite to subsatellite point)'
      write (*,'(a)') 'Eta   is a rotation of the scanline around the satellite velocity vector (effectively middle of scanline is not the subsatellite point)'
      write (*,'(a)') ''
      write (*,'(a)') 'The following defaults are built in :  '
      write (*,'(a)') ''
      write (*,'(a)') ' Satellite     theta (S)  theta(N)   linedelay    Eta (S)  Eta(N)   BowTie Theta on-axis  Bowtie-theta edge FOV   Bowtie Gamma'
      write (*,'(a)') ' NOAA-18          -1.4       1.4         0'
      write (*,'(a)') ' NOAA-19           0.0       0.0         0'
      write (*,'(a)') ' Metop-A          -2.9       2.9         0'
      write (*,'(a)') ' Metop-B          -2.9       2.9         0'
      write (*,'(a)') ' Metop-C          -2.9       2.9         0'
      write (*,'(a)') ' FY3B              0.0       0.0        -4'
      write (*,'(a)') ' FY3C              0.0       0.0        -5'
      write (*,'(a)') ' MN2               0.0       0.0         0         2.4                                                                           Use of Eta is experimental !!'
      write (*,'(a)') ' MN2-2             0.0       0.0         3                                                                                       (preliminary/untested)'
      write (*,'(a)') ' AQUA MODIS       -0.100     0.100      -5                                    0.00                0.42                 1.6'
      write (*,'(a)') ''
      write (*,'(a)') 'WSL example for Aqua Modis truecolour projected image generation - note the file naming convention below - MANDATORY'
      write (*,'(a)') ''
      write (*,'(a)') './hrpt.exe panmodis.txt 221 gapcor project=121 gammargb=0.95,0.99,0.93 histcor modispan=correct fast binary merge reflong=70.0 linedelay=0'
      write (*,'(a)') ''
      write (*,'(a)') 'uses an input file called panmodis.txt with 3 identical lines of text naming the required file eg : /mnt/f/AQ_MODIS-2019-09-29-1700.bin'
      write (*,'(a)') ''
      write (*,'(a)') ''
      write (*,'(a)') ''
      write (*,'(a)') ''
c
      return
      end
