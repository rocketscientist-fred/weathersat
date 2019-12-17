      implicit none
      integer*1 inbuf(27728), linevalid1(10000), linevalid2(10000)
      integer*4 lun1, lun2, lun3, nrecl, irec, ios,doy, nrec1, nrec2, nchan
      integer*2 channel_image(2048,10)
      integer*4 linedoy1(10000), linedoy2(10000), i1_low, i1_upp, i2_low, i2_upp, i, j
      real*8    timestamp, timstart, linetime1(10000), linetime2(10000)
      character*250 filedata1, filedata2, argstring
c
      call getarg(1, argstring)
      filedata1 = argstring(1:lnblnk(argstring))
      call getarg(2, argstring)
      filedata2 = argstring(1:lnblnk(argstring))
c
      nchan = 5
      nrecl = 13864
      call get_lun(lun1)
      open (unit=lun1,file=filedata1(1:lnblnk(filedata1)),form='unformatted',access='direct',recl=nrecl)
      irec = 1
      do while (.true.)
        ios = 0
        call get_scanline(lun1, irec, ios, channel_image, 2048, nchan, 'ME', inbuf, nrecl, timestamp, doy, .false.)
        if (ios.ne.0) goto 1
        linetime1(irec) = timestamp
        linedoy1(irec)  = doy
        irec = irec + 1
      enddo
 1    irec  = irec - 1
      nrec1 = irec
      call check_linetimes(linetime1, linedoy1, linevalid1, irec, timstart, .false.)
c
      call get_lun(lun2)
      open (unit=lun2,file=filedata2(1:lnblnk(filedata2)),form='unformatted',access='direct',recl=nrecl)
      irec = 1
      do while (.true.)
        ios = 0
        call get_scanline(lun2, irec, ios, channel_image, 2048, nchan, 'ME', inbuf, nrecl, timestamp, doy, .false.)
        if (ios.ne.0) goto 2
        linetime2(irec) = timestamp
        linedoy2(irec)  = doy
        irec = irec + 1
      enddo
 2    irec  = irec - 1
      nrec2 = irec
      call check_linetimes(linetime2, linedoy2, linevalid2, irec, timstart, .false.)
      if (linetime1(1).le.linetime2(1)) then
        if (linetime2(1).gt.linetime1(nrec1)) stop '** No overlap **'
        j = 1
        do while (j.le.nrec2)
          do i = 1, nrec1
            if (linetime1(i).eq.linetime2(j)) then
              i1_low = 1
              i1_upp = i - 1
              i2_low = j
              i2_upp = nrec2
              goto 3
            endif
          enddo
          j = j + 1
        enddo
      else
        if (linetime1(1).gt.linetime2(nrec2)) stop '** No overlap **'
        j = 1
        do while (j.le.nrec1)
          do i = 1, nrec2
            if (linetime2(i).eq.linetime1(j)) then
              i1_low = j
              i1_upp = nrec1
              i2_low = 1
              i2_upp = i - 1
              goto 3
            endif
          enddo
          j = j + 1
        enddo
      endif
 3    continue
      write (*,'(''Merging file 1 from line : '',I8,'' to : '',I8,'' with file 2 from line : '',I8,'' to : '',I8)') i1_low, i1_upp, i2_low, i2_upp
      irec = 1
      i    = index(filedata1,'.hpt')
      call get_lun(lun3)
      open (unit=lun3,file=filedata1(1:i-1)//'_Merged.hpt',form='unformatted',access='direct',recl=nrecl)
      do i = i1_low, i1_upp
        call get_file_buffer(lun1, i   , inbuf, nrecl)
        call put_file_buffer(lun3, irec, inbuf, nrecl)
        irec = irec + 1
      enddo
      do i = i2_low, i2_upp
        call get_file_buffer(lun2, i   , inbuf, nrecl)
        call put_file_buffer(lun3, irec, inbuf, nrecl)
        irec = irec + 1
      enddo
      close (unit=lun1)
      close (unit=lun2)
      close (unit=lun3)
      call free_lun(lun1)
      call free_lun(lun2)
      call free_lun(lun3)
      stop
      end

      subroutine get_file_buffer(lun,irec,array,n)
      implicit none
      integer*4 n, irec, lun
      integer*1 array(n)
c
      read (lun,rec=irec) array
      return
      end

      subroutine put_file_buffer(lun,irec,array,n)
      implicit none
      integer*4 n, irec, lun
      integer*1 array(n)
c
      write (lun,rec=irec) array
      return
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
