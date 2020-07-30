      implicit none
      logical swpng, swdebug, swcorrect, swnight, swfitlin, swterra
c
      integer*4 n_max, n_actual, n_req, n_valid
      parameter (n_max=1024, n_valid=2000000)
c
      integer*1 inbuf(642), timbuf(28), timbuf_zero(28), inbuf_night(276), buf(n_max), i_b1, i_b2, i_p1, i_p2, rotbuf(8), outbuf(276), rotout(284), pkt_type, irec_valid(n_valid)
c-- The image buffers - cater also for 13H and 14H
      integer*2 buffer_250m(216640,2), buffer_500m(54160,5), buffer_1km(13540,31), buffer_zero_250m(216640,2), buffer_zero_500m(54160,5), buffer_zero_1km(13540,31)
      integer*2 npos, index_correct(10000), rgbbuf(3, 5416), rgbbuf_correct(3,10000), bwbuf(5416,2)
      integer*4 lun, ios, irec, frame_count, mirror_side, packet_in_group, dy1958, doy, iyr, timhist(8640), yrhist(105), doyhist(366), packet_nr, pkt_length, msechist(8640000)
      integer*4 nlines(38), lunband(38), nrecl(38), irecband, luntim(38), doy_zero, imn, ida, i_hit, i_ptr, n_pkt, i_sav, j_ptr, maxlengap
      integer*4 i, i_max, f_min, i_doy, i_yr, i_sum, j, k, l, i_ext, i_dir, ib, it_ref, i_tref, i_arg, lunnight, i_tref_init, ncorrect_pix
      real*8    timestamp, t_min, t_max, timestamp_zero, theta1, theta_zero1, theta_pk_pk1, theta2, theta_zero2, theta_pk_pk2, timestamp_prev, timediff(n_valid), t_thresh
      character*2   cband
      character*250 argstring, outstring, command
c
      equivalence (timbuf(1), timestamp), (timbuf(9), doy), (timbuf(13), theta1), (timbuf(21), theta2)
      equivalence (timbuf_zero(1), timestamp_zero), (timbuf_zero(9), doy_zero), (timbuf_zero(13), theta_zero1), (timbuf_zero(21), theta_zero2)
c
      i_tref       = 8640001
      maxlengap    = 125
c-- The bow tie correction requires two different angles for left and right side correction  (at least for Ch 1) ........
c-- For a North direction pass, theta1 is the right hand side of the image :-) (which is -Y in hrsc.f .....)
c-- The bowgamma is defined in hrpt.f
      theta_pk_pk1 = 0.00D0
      theta_pk_pk2 = 0.42D0
c-- Value 2 is for the left hand side (as seen in the track direction - so left for a N pass); Value 1 for the right hand side
c
c-- Set all records to invalid
      do i = 1, n_valid
        irec_valid(i) = 1
        timediff(i)   = 0.0D0
      enddo
      t_thresh  = 5.0D0
      call get_lun(lun)
      swpng     = .true.
      swdebug   = .false.
      swcorrect = .false.
      swnight   = .false.
      swfitlin  = .false.
      swterra   = .false.
      do i_arg = 2,10
        call getarg(i_arg,argstring)
        if (index(argstring,'nopng').ne.0)   swpng     = .false.
        if (index(argstring,'debug').ne.0)   swdebug   = .true.
        if (index(argstring,'night').ne.0)   swnight   = .true.
        if (index(argstring,'correct').ne.0) swcorrect = .true.
        if (index(argstring,'bowtie=').ne.0) then
          i = index(argstring,'=') + 1
          j = index(argstring,',') - 1
          k = index(argstring,',') + 1
          l = lnblnk(argstring)
          if (i.ne.0.and.j.ne.0.and.k.ne.0.and.l.ne.0) then
            call string_to_r8(argstring(i:j)//' ',npos, theta_pk_pk1)
            write (*,'(''BowTie pk-pk1 override     : '',F10.4)') theta_pk_pk1
            call string_to_r8(argstring(k:l)//' ',npos, theta_pk_pk2)
            write (*,'(''BowTie pk-pk2 override     : '',F10.4)') theta_pk_pk2
          endif
        endif
        if (index(argstring,'fitlin').ne.0) swfitlin = .true.
      enddo
      call getarg(1,argstring)
      if (index(argstring,'TE_').ne.0) swterra = .true.
      n_actual = -1
      i_ext = index(argstring,'.bin')
      irec  = 1
      if (.not.swnight) then
        open (unit=lun,file=argstring(1:lnblnk(argstring)), access='direct',form='unformatted',recl=642)
      else
c-- Pre-filter to get rid of engineering packets
        call get_lun(lunnight)
        open (unit=lunnight,file=argstring(1:i_ext-1)//'.nig', access='direct',form='unformatted',recl=276)
        i_ptr = 0
        do i = lnblnk(argstring), 1, -1
          if (argstring(i:i).eq.'/'.and.i_ptr.eq.0) i_ptr = i
        enddo
        read (argstring(i_ptr+10:i_ptr+13),'(I4)')   iyr
        read (argstring(i_ptr+15:i_ptr+16),'(I2.2)') imn
        read (argstring(i_ptr+18:i_ptr+19),'(I2.2)') ida
        call datetodoy(iyr, imn, ida, doy)
        call doy_to_dy1958(dy1958, doy, iyr)
        i_b2  =       iand(dy1958,'000000FF'x)
        i_b1  = ishft(iand(dy1958,'0000FF00'x),-8)
        n_pkt = 269
        i_p2  =       iand(n_pkt ,'000000FF'x)
        i_p1  = ishft(iand(n_pkt ,'0000FF00'x),-8)
        i_hit = 0
        i_ptr = 0
        i_sav = 0
        j_ptr = 0
        do while (.true.)
          n_req = 1024
          call get_TM(argstring(1:lnblnk(argstring)), buf, n_max, n_actual, n_req)
          if (n_actual.eq.n_req) then
            do i = 1,n_actual
              do j = 1,7
                rotbuf(j) = rotbuf(j+1)
              enddo
              rotbuf(8) = buf(i)
              do j = 1,283
                rotout(j) = rotout(j+1)
              enddo
              rotout(284) = buf(i)
              i_ptr = i_ptr + 1
              if (buf(i).eq.i_p1.and.i_hit.eq.0) then
                i_hit = 1
                goto 10
              endif
              if (buf(i).eq.i_p2.and.i_hit.eq.1) then
                i_hit = 2
                goto 10
              endif
              if (buf(i).eq.i_b1.and.i_hit.eq.2) then
                i_hit = 3
                goto 10
              endif
              if (buf(i).eq.i_b2.and.i_hit.eq.3) then
c-- put the 8 bytes in a buffer and set j_ptr to 9
                write (*,'(I10,2x,I10,2x,8(Z2.2,1x))') i_ptr, i_ptr-i_sav, (rotbuf(k),k=1,8)
c                read (*,*)
                if (i_ptr-i_sav.eq.276) then
                  do k = 1, 276
                    outbuf(k) = rotout(k)
                  enddo
                  write (lunnight,rec=irec) outbuf
                  write (*,*) irec
                  irec = irec + 1
                endif
                i_sav = i_ptr
                i_hit = 0
                goto 10
              endif
              i_hit = 0
 10           continue
            enddo
c
c            if (pkt_type.eq.1.and.pkt_length.eq.269) then
c              n_req = 261
c              call get_TM(argstring(1:lnblnk(argstring)), buf(16), n_max, n_actual, n_req)
c            else if ((pkt_type.eq.2.or.pkt_type.eq.4).and.pkt_length.eq.635) then
c              n_req = 627
c              call get_TM(argstring(1:lnblnk(argstring)), buf(16), n_max, n_actual, n_req)
c            else
c              n_req = 261
c              call get_TM(argstring(1:lnblnk(argstring)), buf(16), n_max, n_actual, n_req)
c            endif
          else
            write (*,'(''Night file error : '',2I10)') n_actual, n_req
            goto 11
          endif
        enddo
 11      close (unit=lunnight)
        call free_lun(lunnight)
c
        open (unit=lun,file=argstring(1:i_ext-1)//'.nig', access='direct',form='unformatted',recl=276)
      endif
      i_dir = 0
      do i = lnblnk(argstring), 1, -1
        if (argstring(i:i).eq.'/'.and.i_dir.eq.0) i_dir = i
      enddo
      call getenv('HRPTOUT',outstring)
      do i = 1, 38
        nlines(i) = 0
      enddo
      do j = 1, 2
        do i = 1, 216640
          buffer_250m(i,j) = 0
          buffer_zero_250m(i,j) = 0
        enddo
      enddo
      do j = 1, 5
        do i = 1, 54160
          buffer_500m(i,j) = 0
          buffer_zero_500m(i,j) = 0
        enddo
      enddo
      do j = 1, 31
        do i = 1, 13540
          buffer_1km(i,j) = 0
          buffer_zero_1km(i,j) = 0
        enddo
      enddo
      do i = 1, 28
        timbuf_zero(i) = 0
      enddo
      theta1      = 0.0D0
      theta_zero1 = 0.0D0
      theta2      = 0.0D0
      theta_zero2 = 0.0D0
      if (.not.swnight) then
        do i = 1, 38
          call get_lun(lunband(i))
          call get_lun(luntim(i))
        enddo
        do i = 1, 2
          nrecl(i) = 5416
          write (cband,'(I2.2)') i
          open (unit=lunband(i),file=outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.dat', access='direct',form='unformatted',recl=nrecl(i)*2)
          open (unit=luntim(i) ,file=outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.tim', access='direct',form='unformatted',recl=28)
        enddo
        do i = 3, 7
          nrecl(i) = 2708
          write (cband,'(I2.2)') i
          open (unit=lunband(i),file=outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.dat', access='direct',form='unformatted',recl=nrecl(i)*2)
          open (unit=luntim(i) ,file=outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.tim', access='direct',form='unformatted',recl=28)
        enddo
        do i = 8, 38
          nrecl(i) = 1354
          write (cband,'(I2.2)') i
          open (unit=lunband(i),file=outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.dat', access='direct',form='unformatted',recl=nrecl(i)*2)
          open (unit=luntim(i) ,file=outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.tim', access='direct',form='unformatted',recl=28)
        enddo
      else
        do i = 20, 36
          call get_lun(lunband(i))
          call get_lun(luntim(i))
        enddo
        do i = 20, 36
          nrecl(i) = 1354
          write (cband,'(I2.2)') i
          open (unit=lunband(i),file=outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.dat', access='direct',form='unformatted',recl=nrecl(i)*2)
          open (unit=luntim(i) ,file=outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.tim', access='direct',form='unformatted',recl=28)
        enddo
      endif
c
c-- Determine valid time window, year and day of year by majority "voting" and use as criteria to filter bad records.
c
      do i = 1, 8640000
        msechist(i) = 0
      enddo
      do i = 1, 8640
        timhist(i) = 0
      enddo
      do i = 1, 105
        yrhist(i) = 0
      enddo
      do i = 1, 366
        doyhist(i) = 0
      enddo
      irec = 1
      do while (.true.)
        ios  = 0
        if (.not.swnight) then
          read (lun, rec=irec, iostat=ios) inbuf
        else
          read (lun, rec=irec, iostat=ios) inbuf_night
          do i = 1, 276
            inbuf(i) = inbuf_night(i)
          enddo
        endif
        if (ios.ne.0) goto 1
        dy1958 = inbuf(7)
        if (inbuf(7).lt.0) dy1958 = dy1958 + 256
        dy1958 = dy1958 * 256
        dy1958 = dy1958 + inbuf(8)
        if (inbuf(8).lt.0) dy1958 = dy1958 + 256
        call dy1958_to_doy(dy1958, doy, iyr)
        call MODIS_time(inbuf(9),6,timestamp)
        i = int(timestamp/10.0D0) + 1
        if (1.le.i.and.i.le.8640) timhist(i) = timhist(i) + 1
        i = int(timestamp*100.0D0) + 1
        if (1.le.i.and.i.le.8640000) msechist(i) = msechist(i) + 1
        i = iyr - 1995
        if (1.le.i.and.i.le.105) yrhist(i) = yrhist(i) + 1
        i = doy
        if (1.le.i.and.i.le.366) doyhist(i) = doyhist(i) + 1
        irec = irec + 1
      enddo
 1    i_max = -1
      if (swdebug) write (*,*) '# Records          : ', irec
      if (irec.gt.n_valid) stop "** Increase n_valid **"
      do i = 1, 8640
        if (timhist(i).ge.i_max) i_max = i
      enddo
c-- Take 10% of the maximum as a lower limit window
      f_min = timhist(i_max) / 10
      t_min = 0.0D0
      t_max = 0.0D0
      i_sum = 0
      do i = 1,8640
        if (timhist(i).ge.f_min) then 
          i_sum = i_sum + timhist(i)
        endif
      enddo
      if (swdebug) write (*,*) '# Valid timestamps : ', i_sum
      do i = 1,8640
        if (swdebug) write (*,*) i, timhist(i), f_min
        if (timhist(i).ge.f_min.and.t_min.eq.0.0D0) t_min = dble(i-1) * 10.0D0 + 5.0
      enddo
      do i = 8640, 1, -1
        if (timhist(i).ge.f_min.and.t_max.eq.0.0D0) t_max = dble(i+1) * 10.0D0 - 5.0  
      enddo
      write (*,'(''Time window validity     : '',2F15.2)') t_min, t_max
      i_max = -1
      do i = 1, 366
        if (doyhist(i).ge.i_max) then 
          i_max = doyhist(i)
          i_doy = i
        endif
      enddo
      write (*,'(''DOY validity             : '',I15)') i_doy
      i_max = -1
      do i = 1, 105
        if (yrhist(i).ge.i_max) then
          i_max = yrhist(i)
          i_yr  = i + 1995
        endif
      enddo
      write (*,'(''Year validity            : '',I15)') i_yr
c-- Determine which records inside the time window are valid (ie outlier elimination)
      irec = 1
      do while (.true.)
        ios  = 0
        if (.not.swnight) then
          read (lun, rec=irec, iostat=ios) inbuf
        else
          read (lun, rec=irec, iostat=ios) inbuf_night
          do i = 1, 276
            inbuf(i) = inbuf_night(i)
          enddo
        endif
        if (ios.ne.0) goto 12
        if (ishft(iand(inbuf(16),'80'x),-7).eq.0) then
          pkt_type    = ishft(iand(inbuf(15),'70'x),-4)
          dy1958 = inbuf(7)
          if (inbuf(7).lt.0) dy1958 = dy1958 + 256
          dy1958 = dy1958 * 256
          dy1958 = dy1958 + inbuf(8)
          if (inbuf(8).lt.0) dy1958 = dy1958 + 256
          call dy1958_to_doy(dy1958, doy, iyr)
          call MODIS_time(inbuf(9),6,timestamp)
          if (irec.gt.1) then
            if (t_min.le.timestamp.and.timestamp.le.t_max.and.iyr.eq.i_yr.and.doy.eq.i_doy.and.pkt_type.eq.0) then
c            if (t_min.le.timestamp.and.timestamp.le.t_max) then
              timediff(irec)   = timestamp - timestamp_prev
              irec_valid(irec) = 0
              timestamp_prev = timestamp
            endif
          endif
        endif
        irec = irec + 1
      enddo
 12   i = 1
      j = 1
      do while (j.le.irec-1)
        if (irec_valid(i).eq.0) then
          j = i + 1
          do while (irec_valid(j).ne.0.and.j.le.irec-1)
            j = j + 1
          enddo
          if (abs(timediff(i)).gt.t_thresh.and.(abs(timediff(i)+timediff(j)).lt.t_thresh)) irec_valid(i) = 1
          i = j
        else
          i = i + 1
        endif
      enddo
c-- End of outlier determination
c
c-- Test of 10 msec bin time hist as filter/trigger criteria for write
      i_max = -1
      do i = 1, 8640000
        if (msechist(i).ge.i_max) i_max = msechist(i)
      enddo
      if (swdebug) write (*,*) i_max
      j = 1
      it_ref = 0
      i_tref_init = -1
      do i = 1, 8640000
        if (msechist(i).ge.int(i_max/10)) then
          if (i_tref_init.lt.0.and.i.ge.nint(100.0*t_min)) i_tref_init = i
          if (swdebug) write (*,*) j, i, msechist(i), i - it_ref
          j      = j + 1
          it_ref = i
        endif
      enddo
c-- End of determining filter criteria
c-- Now prepare the corrections per idet/ifov - if requested
      if (swcorrect) then
        call modis_hist_init()
        irec = 1
c-- Separate the day and night passes - duplicates code, but simpler overview !
        if (.not.swnight) then
c-- Daytime data processing - is default
          i_tref = i_tref_init
          do while (.true.)
            ios  = 0
            read (lun, rec=irec, iostat=ios) inbuf
            if (ios.ne.0) goto 2
            if (ishft(iand(inbuf(16),'80'x),-7).eq.0) then
              packet_in_group = ishft(iand(inbuf(3),'C0'x),-6)
              packet_nr       =       iand(inbuf(3),'3F'x) * 256
              packet_nr       = packet_nr + inbuf(4)
              if (inbuf(4).lt.0) packet_nr = packet_nr + 256
              frame_count = iand(inbuf(16),'7F'x)
              frame_count = frame_count * 16 + ishft(iand(inbuf(17),'F0'x),-4)
              pkt_type    = ishft(iand(inbuf(15),'70'x),-4)
              mirror_side = iand(inbuf(15),'01'x)
              dy1958 = inbuf(7)
              if (inbuf(7).lt.0) dy1958 = dy1958 + 256
              dy1958 = dy1958 * 256
              dy1958 = dy1958 + inbuf(8)
              if (inbuf(8).lt.0) dy1958 = dy1958 + 256
              call dy1958_to_doy(dy1958, doy, iyr)
              call MODIS_time(inbuf(9),6,timestamp)
c-- Timestamp and doy go (modified) to the time buffer for the georeferencing SGP4 TLE calc
c              if (t_min.le.timestamp.and.timestamp.le.t_max.and.iyr.eq.i_yr.and.doy.eq.i_doy.and.pkt_type.eq.0) then
              if (irec_valid(irec).eq.0) then
                timestamp_prev = timestamp
c-- Check the timestamp to see whether a 'flush' is needed.
                i = int(timestamp*100.0D0) + 1
                if (msechist(i).gt.int(i_max/10).and.i.gt.i_tref) then
c-- Before we flush try to capture the wrong times that slip through - this could be done better by pre-building a tble to flag out of sequence times
                  if (nint((float(i)-float(i_tref))/147.7).lt.0.or.nint((float(i)-float(i_tref))/147.7).gt.maxlengap) goto 3
c-- Write band01 and band02
                  do ib = 1,2
                    call modis_hist(buffer_250m(1,ib), 216640, ib, swterra)
                  enddo
                  do j = 1,2
                    do i = 1, 216640
                      buffer_250m(i,j) = 0
                    enddo
                  enddo
c-- Write band03 through band07
                  do ib = 3,7
                    call modis_hist(buffer_500m(1,ib-2), 54160, ib, swterra)
                  enddo
                  do j = 1,5
                    do i = 1, 54160
                      buffer_500m(i,j) = 0
                    enddo
                  enddo
c-- Write band08 through band38
                  do ib = 8,38
                    call modis_hist(buffer_1km(1,ib-7), 13540, ib, swterra)
                  enddo
                  do j = 1,31
                    do i = 1, 13540
                      buffer_1km(i,j) = 0
                    enddo
                  enddo
                endif
                i = int(timestamp*100.0D0) + 1
                if (msechist(i).gt.int(i_max/10).and.(nint((float(i)-float(i_tref))/147.7).gt.0)) i_tref = i
                call band_image(inbuf, 642, packet_in_group, 1, buffer_250m(1,1), 4, 1354, 4, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group, 2, buffer_250m(1,2), 4, 1354, 4, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group, 3, buffer_500m(1,1), 2, 1354, 2, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group, 4, buffer_500m(1,2), 2, 1354, 2, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group, 5, buffer_500m(1,3), 2, 1354, 2, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group, 6, buffer_500m(1,4), 2, 1354, 2, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group, 7, buffer_500m(1,5), 2, 1354, 2, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group, 8, buffer_1km(1, 1) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group, 9, buffer_1km(1, 2) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,10, buffer_1km(1, 3) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,11, buffer_1km(1, 4) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,12, buffer_1km(1, 5) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,13, buffer_1km(1, 6) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,14, buffer_1km(1, 7) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,15, buffer_1km(1, 8) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,16, buffer_1km(1, 9) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,17, buffer_1km(1,10) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,18, buffer_1km(1,11) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,19, buffer_1km(1,12) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,20, buffer_1km(1,13) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,21, buffer_1km(1,14) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,22, buffer_1km(1,15) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,23, buffer_1km(1,16) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,24, buffer_1km(1,17) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,25, buffer_1km(1,18) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,26, buffer_1km(1,19) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,27, buffer_1km(1,20) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,28, buffer_1km(1,21) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,29, buffer_1km(1,22) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,30, buffer_1km(1,23) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,31, buffer_1km(1,24) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,32, buffer_1km(1,25) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,33, buffer_1km(1,26) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,34, buffer_1km(1,27) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,35, buffer_1km(1,28) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,36, buffer_1km(1,29) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,37, buffer_1km(1,30) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image(inbuf, 642, packet_in_group,38, buffer_1km(1,31) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              endif
            endif
 3          irec = irec + 1
          enddo
 2        irec = irec - 1
c-- Deal with the pending buffers
c-- Write band01 and band02
          do ib = 1,2
            call modis_hist(buffer_250m(1,ib), 216640, ib, swterra)
          enddo
          do j = 1,2
            do i = 1, 216640
              buffer_250m(i,j) = 0
            enddo
          enddo
c-- Write band03 through band07
          do ib = 3,7
            call modis_hist(buffer_500m(1,ib-2), 54160, ib, swterra)
          enddo
          do j = 1,5
            do i = 1, 54160
              buffer_500m(i,j) = 0
            enddo
          enddo
c-- Write band08 through band38
          do ib = 8,38
            call modis_hist(buffer_1km(1,ib-7), 13540, ib, swterra)
          enddo
          do j = 1,31
            do i = 1, 13540
              buffer_1km(i,j) = 0
            enddo
          enddo
        else
c-- Night data processing
          i_tref = i_tref_init
          do while (.true.)
            ios  = 0
            read (lun, rec=irec, iostat=ios) inbuf_night
            if (ios.ne.0) goto 4
            if (ishft(iand(inbuf_night(16),'80'x),-7).eq.0) then
              packet_in_group = ishft(iand(inbuf_night(3),'C0'x),-6)
              packet_nr       =       iand(inbuf_night(3),'3F'x) * 256
              packet_nr       = packet_nr + inbuf_night(4)
              if (inbuf_night(4).lt.0) packet_nr = packet_nr + 256
              frame_count = iand(inbuf_night(16),'7F'x)
              frame_count = frame_count * 16 + ishft(iand(inbuf_night(17),'F0'x),-4)
              pkt_type    = ishft(iand(inbuf_night(15),'70'x),-4)
              mirror_side = iand(inbuf_night(15),'01'x)
              dy1958 = inbuf_night(7)
              if (inbuf_night(7).lt.0) dy1958 = dy1958 + 256
              dy1958 = dy1958 * 256
              dy1958 = dy1958 + inbuf_night(8)
              if (inbuf_night(8).lt.0) dy1958 = dy1958 + 256
              call dy1958_to_doy(dy1958, doy, iyr)
              call MODIS_time(inbuf_night(9),6,timestamp)
c-- Timestamp and doy go (modified) to the time buffer for the georeferencing SGP4 TLE calc
c              if (t_min.le.timestamp.and.timestamp.le.t_max.and.iyr.eq.i_yr.and.doy.eq.i_doy.and.pkt_type.eq.0) then
              if (irec_valid(irec).eq.0) then
                timestamp_prev = timestamp
c-- Check the timestamp to see whether a 'flush' is needed.
                i = int(timestamp*100.0D0) + 1
                if (msechist(i).gt.int(i_max/10).and.i.gt.i_tref) then
c-- Before we flush try to capture the wrong times that slip through - this could be done better by pre-building a tble to flag out of sequence times
                  if (nint((float(i)-float(i_tref))/147.7).lt.0.or.nint((float(i)-float(i_tref))/147.7).gt.maxlengap) goto 5
c-- Write band20 through band36
                  do ib = 20,36
                    call modis_hist(buffer_1km(1,ib-7), 13540, ib, swterra)
                  enddo
                  do j = 1,31
                    do i = 1, 13540
                      buffer_1km(i,j) = 0
                    enddo
                  enddo
                endif
                i = int(timestamp*100.0D0) + 1
                if (msechist(i).gt.int(i_max/10).and.(nint((float(i)-float(i_tref))/147.7).gt.0)) i_tref = i
                call band_image_night(inbuf_night, 276, packet_in_group,20, buffer_1km(1,13) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,21, buffer_1km(1,14) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,22, buffer_1km(1,15) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,23, buffer_1km(1,16) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,24, buffer_1km(1,17) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,25, buffer_1km(1,18) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,26, buffer_1km(1,19) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,27, buffer_1km(1,20) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,28, buffer_1km(1,21) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,29, buffer_1km(1,22) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,30, buffer_1km(1,23) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,31, buffer_1km(1,24) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,32, buffer_1km(1,25) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,33, buffer_1km(1,26) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,34, buffer_1km(1,27) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,35, buffer_1km(1,28) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
                call band_image_night(inbuf_night, 276, packet_in_group,36, buffer_1km(1,29) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              endif
            endif
 5          irec = irec + 1
          enddo
 4        irec = irec - 1
c-- Deal with the pending buffers
c-- Write band20 through band36
          do ib = 20,36
            call modis_hist(buffer_1km(1,ib-7), 13540, ib, swterra)
          enddo
          do j = 1,31
            do i = 1, 13540
              buffer_1km(i,j) = 0
            enddo
          enddo
        endif
        call modis_hist_flush(swnight, swfitlin)
      endif
c--  End of correction preparation
      irec = 1
c-- Separate the day and night passes - duplicates code, but simpler overview !
      if (.not.swnight) then
c-- Daytime data processing - is default
        i_tref = i_tref_init
        do while (.true.)
          ios  = 0
          read (lun, rec=irec, iostat=ios) inbuf
          if (ios.ne.0) goto 6
          if (ishft(iand(inbuf(16),'80'x),-7).eq.0) then
            packet_in_group = ishft(iand(inbuf(3),'C0'x),-6)
            packet_nr       =       iand(inbuf(3),'3F'x) * 256
            packet_nr       = packet_nr + inbuf(4)
            if (inbuf(4).lt.0) packet_nr = packet_nr + 256
            frame_count = iand(inbuf(16),'7F'x)
            frame_count = frame_count * 16 + ishft(iand(inbuf(17),'F0'x),-4)
            pkt_type    = ishft(iand(inbuf(15),'70'x),-4)
            mirror_side = iand(inbuf(15),'01'x)
            dy1958 = inbuf(7)
            if (inbuf(7).lt.0) dy1958 = dy1958 + 256
            dy1958 = dy1958 * 256
            dy1958 = dy1958 + inbuf(8)
            if (inbuf(8).lt.0) dy1958 = dy1958 + 256
            call dy1958_to_doy(dy1958, doy, iyr)
            call MODIS_time(inbuf(9),6,timestamp)
c-- Timestamp and doy go (modified) to the time buffer for the georeferencing SGP4 TLE calc
c            if (t_min.le.timestamp.and.timestamp.le.t_max.and.iyr.eq.i_yr.and.doy.eq.i_doy.and.pkt_type.eq.0.and.irec_valid(irec).eq.0) then
            if (irec_valid(irec).eq.0) then
              if (swdebug) write (*,'(''Packet Count : '',I7,1x,''Packet # in group : '',I2,1x,''Mirror_side : '',I10,1x,''Frame Count : '',I10,'' Year, doy, time : '',2I7,F15.6,I15,F15.6,I10)') 
     *          packet_nr, packet_in_group, mirror_side + 1, frame_count, iyr, doy, timestamp, dy1958, timestamp - timestamp_prev, irec_valid(irec)
              timestamp_prev = timestamp
c-- Check the timestamp to see whether a 'flush' is needed.
              i = int(timestamp*100.0D0) + 1
              if (msechist(i).gt.int(i_max/10).and.i.gt.i_tref) then
c-- Before we flush try to capture the wrong times that slip through - this could be done better by pre-building a tble to flag out of sequence times
                if (swdebug) write (*,*) i, i_tref, nint((float(i)-float(i_tref))/147.7)
                if (nint((float(i)-float(i_tref))/147.7).lt.0.or.nint((float(i)-float(i_tref))/147.7).gt.maxlengap) goto 7
c-- Flush the required empty buffers
                if (nint((float(i)-float(i_tref))/147.7).gt.1) then
                  if (swdebug) write (*,'(''Flushing'',I6,'' empty buffers'')') nint((float(i)-float(i_tref))/147.7) - 1
                  do l = 1, nint((float(i)-float(i_tref))/147.7) - 1
c-- Write band01 and band02
                    doy_zero       = doy
                    do ib = 1,2
c-- write the timestamps
                      theta_zero1 = 0.0D0
                      theta_zero2 = 0.0D0
                      timestamp_zero = dble(i_tref)/100.0D0 + dble (l) * 1.477815
                      do k = nlines(ib)+1, nlines(ib) + 40
                        call write_buffer_i1(luntim(ib), k, timbuf_zero, 28)
                        timestamp_zero = timestamp_zero + (1.477815D0/40.0D0)
                      enddo
                      nlines(ib) = nlines(ib) + 40
                      irecband   = nlines(ib)
                      do i = 1, 10
                        do j = 1,4
                          k = ((i-1) * 4 + j - 1) * 5416 + 1
                          call write_buffer(lunband(ib), irecband, buffer_zero_250m(k,ib), nrecl(ib))
                          irecband = irecband - 1
                        enddo
                      enddo
                    enddo
c-- Write band03 through band07
                    do ib = 3,7
                      theta_zero1 = 0.0D0
                      theta_zero2 = 0.0D0
                      timestamp_zero = dble(i_tref)/100.0D0 + dble (l) * 1.477815
                      do k = nlines(ib)+1, nlines(ib) + 20
                        call write_buffer_i1(luntim(ib), k, timbuf_zero, 28)
                        timestamp_zero = timestamp_zero + (1.477815D0/20.0D0)
                      enddo
                      nlines(ib) = nlines(ib) + 20
                      irecband   = nlines(ib)
                      do i = 1, 10
                        do j = 1,2
                          k = ((i-1) * 2 + j - 1) * 2708 + 1
                          call write_buffer(lunband(ib), irecband, buffer_zero_500m(k,ib-2), nrecl(ib))
                          irecband = irecband - 1
                        enddo
                      enddo
                    enddo
c-- Write band08 through band38
                    do ib = 8,38
                      theta_zero1 = 0.0D0
                      theta_zero2 = 0.0D0
                      timestamp_zero = dble(i_tref)/100.0D0 + dble (l) * 1.477815
                      do k = nlines(ib)+1, nlines(ib) + 10
                        call write_buffer_i1(luntim(ib), k, timbuf_zero, 28)
                        timestamp_zero = timestamp_zero + (1.477815D0/10.0D0)
                      enddo
                      nlines(ib) = nlines(ib) + 10
                      irecband   = nlines(ib)
                      do i = 1, 10
                        do j = 1,1
                          k = ((i-1) * 1 + j - 1) * 1354 + 1
                          call write_buffer(lunband(ib), irecband, buffer_zero_1km(k,ib-7), nrecl(ib))
                          irecband = irecband - 1
                        enddo
                      enddo
                    enddo
                  enddo
                endif
c-- End flushing empty buffers, now flush the real one
c-- Write band01 and band02
                do ib = 1,2
                  if (swcorrect) call modis_hist_correct(buffer_250m(1,ib), 216640, ib, swfitlin)
                  timestamp_zero = timestamp
                  do k = nlines(ib) + 1, nlines(ib) + 40
                    theta1 = (theta_pk_pk1/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk1 / 39.0D0)
                    theta2 = (theta_pk_pk2/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk2 / 39.0D0)
                    call write_buffer_i1(luntim(ib), k, timbuf, 28)
                    timestamp = timestamp + (1.477815D0/40.0D0)
                  enddo
                  timestamp = timestamp_zero
                  nlines(ib) = nlines(ib) + 40
                  irecband   = nlines(ib)
                  do i = 1, 10
                    do j = 1,4
                      k = ((i-1) * 4 + j - 1) * 5416 + 1
                      call write_buffer(lunband(ib), irecband, buffer_250m(k,ib), nrecl(ib))
                      irecband = irecband - 1
                    enddo
                  enddo
                enddo
                do j = 1,2
                  do i = 1, 216640
                    buffer_250m(i,j) = 0
                  enddo
                enddo
c-- Write band03 through band07
                do ib = 3,7
                  if (swcorrect) call modis_hist_correct(buffer_500m(1,ib-2), 54160, ib, swfitlin)
                  timestamp_zero = timestamp
                  do k = nlines(ib) + 1, nlines(ib) + 20
                    theta1 = (theta_pk_pk1/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk1 / 19.0D0)
                    theta2 = (theta_pk_pk2/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk2 / 19.0D0)
                    call write_buffer_i1(luntim(ib), k, timbuf, 28)
                    timestamp = timestamp + (1.477815D0/20.0D0)
                  enddo
                  timestamp = timestamp_zero
                  nlines(ib) = nlines(ib) + 20
                  irecband   = nlines(ib)
                  do i = 1, 10
                    do j = 1,2
                      k = ((i-1) * 2 + j - 1) * 2708 + 1
                      call write_buffer(lunband(ib), irecband, buffer_500m(k,ib-2), nrecl(ib))
                      irecband = irecband - 1
                    enddo
                  enddo
                enddo
                do j = 1,5
                  do i = 1, 54160
                    buffer_500m(i,j) = 0
                  enddo
                enddo
c-- Write band08 through band38
                do ib = 8,38
                  if (swcorrect) call modis_hist_correct(buffer_1km(1,ib-7), 13540, ib, swfitlin)
                  timestamp_zero = timestamp
                  do k = nlines(ib) + 1, nlines(ib) + 10
                    theta1 = (theta_pk_pk1/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk1 /  9.0D0)
                    theta2 = (theta_pk_pk2/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk2 /  9.0D0)
                    call write_buffer_i1(luntim(ib), k, timbuf, 28)
                    timestamp = timestamp + (1.477815D0/10.0D0)
                  enddo
                  timestamp = timestamp_zero
                  nlines(ib) = nlines(ib) + 10
                  irecband   = nlines(ib)
                  do i = 1, 10
                    do j = 1,1
                      k = ((i-1) * 1 + j - 1) * 1354 + 1
                      call write_buffer(lunband(ib), irecband, buffer_1km(k,ib-7), nrecl(ib))
                      irecband = irecband - 1
                    enddo
                  enddo
                enddo
                do j = 1,31
                  do i = 1, 13540
                    buffer_1km(i,j) = 0
                  enddo
                enddo
              endif
              i = int(timestamp*100.0D0) + 1
              if (msechist(i).gt.int(i_max/10).and.(nint((float(i)-float(i_tref))/147.7).gt.0)) i_tref = i
              call band_image(inbuf, 642, packet_in_group, 1, buffer_250m(1,1), 4, 1354, 4, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group, 2, buffer_250m(1,2), 4, 1354, 4, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group, 3, buffer_500m(1,1), 2, 1354, 2, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group, 4, buffer_500m(1,2), 2, 1354, 2, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group, 5, buffer_500m(1,3), 2, 1354, 2, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group, 6, buffer_500m(1,4), 2, 1354, 2, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group, 7, buffer_500m(1,5), 2, 1354, 2, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group, 8, buffer_1km(1, 1) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group, 9, buffer_1km(1, 2) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,10, buffer_1km(1, 3) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,11, buffer_1km(1, 4) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,12, buffer_1km(1, 5) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,13, buffer_1km(1, 6) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,14, buffer_1km(1, 7) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,15, buffer_1km(1, 8) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,16, buffer_1km(1, 9) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,17, buffer_1km(1,10) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,18, buffer_1km(1,11) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,19, buffer_1km(1,12) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,20, buffer_1km(1,13) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,21, buffer_1km(1,14) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,22, buffer_1km(1,15) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,23, buffer_1km(1,16) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,24, buffer_1km(1,17) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,25, buffer_1km(1,18) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,26, buffer_1km(1,19) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,27, buffer_1km(1,20) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,28, buffer_1km(1,21) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,29, buffer_1km(1,22) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,30, buffer_1km(1,23) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,31, buffer_1km(1,24) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,32, buffer_1km(1,25) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,33, buffer_1km(1,26) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,34, buffer_1km(1,27) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,35, buffer_1km(1,28) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,36, buffer_1km(1,29) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,37, buffer_1km(1,30) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image(inbuf, 642, packet_in_group,38, buffer_1km(1,31) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
            endif
          endif
 7        irec = irec + 1
        enddo
 6      irec = irec - 1
        if (swdebug) write (*,*) irec
c-- Write the pending buffer
c-- Write band01 and band02
        do ib = 1,2
          if (swcorrect) call modis_hist_correct(buffer_250m(1,ib), 216640, ib, swfitlin)
          timestamp_zero = timestamp
          do k = nlines(ib) + 1, nlines(ib) + 40
            theta1 = (theta_pk_pk1/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk1 / 39.0D0)
            theta2 = (theta_pk_pk2/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk2 / 39.0D0)
            call write_buffer_i1(luntim(ib), k, timbuf, 28)
            timestamp = timestamp + (1.477815D0/40.0D0)
          enddo
          timestamp = timestamp_zero
          nlines(ib) = nlines(ib) + 40
          irecband   = nlines(ib)
          do i = 1, 10
            do j = 1,4
              k = ((i-1) * 4 + j - 1) * 5416 + 1
              call write_buffer(lunband(ib), irecband, buffer_250m(k,ib), nrecl(ib))
              irecband = irecband - 1
            enddo
          enddo
        enddo
c-- Write band03 through band07
        do ib = 3,7
          if (swcorrect) call modis_hist_correct(buffer_500m(1,ib-2), 54160, ib, swfitlin)
          timestamp_zero = timestamp
          do k = nlines(ib) + 1, nlines(ib) + 20
            theta1 = (theta_pk_pk1/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk1 / 19.0D0)
            theta2 = (theta_pk_pk2/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk2 / 19.0D0)
            call write_buffer_i1(luntim(ib), k, timbuf, 28)
            timestamp = timestamp + (1.477815D0/20.0D0)
          enddo
          timestamp = timestamp_zero
          nlines(ib) = nlines(ib) + 20
          irecband   = nlines(ib)
          do i = 1, 10
            do j = 1,2
              k = ((i-1) * 2 + j - 1) * 2708 + 1
              call write_buffer(lunband(ib), irecband, buffer_500m(k,ib-2), nrecl(ib))
              irecband = irecband - 1
            enddo
          enddo
        enddo
c-- Write band08 through band38
        do ib = 8,38
          if (swcorrect) call modis_hist_correct(buffer_1km(1,ib-7), 13540, ib, swfitlin)
          timestamp_zero = timestamp
          do k = nlines(ib) + 1, nlines(ib) + 10
            theta1 = (theta_pk_pk1/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk1 /  9.0D0)
            theta2 = (theta_pk_pk2/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk2 /  9.0D0)
            call write_buffer_i1(luntim(ib), k, timbuf, 28)
            timestamp = timestamp + (1.477815D0/10.0D0)
          enddo
          timestamp = timestamp_zero
          nlines(ib) = nlines(ib) + 10
          irecband   = nlines(ib)
          do i = 1, 10
            do j = 1,1
              k = ((i-1) * 1 + j - 1) * 1354 + 1
              call write_buffer(lunband(ib), irecband, buffer_1km(k,ib-7), nrecl(ib))
              irecband = irecband - 1
            enddo
          enddo
        enddo
        close (unit=lun)
        call free_lun(lun)
        do i = 1, 38
          close (unit=lunband(i))
          call free_lun(lunband(i))
          write (*,'('' # Lines in image band  '',I2.2,'' : '',I10)') i, nlines(i)
          close (unit=luntim(i))
          call free_lun(luntim(i))
        enddo
        if (swpng) then
c-- Create PNG's for band08 through band38 - remember, there are 36 bands, but 13 and 14 have a L and a H version, thus the total of 38
          do ib = 8, 38
            write (cband,'(I2.2)') ib
            write (command,'(''convert -depth 16 -endian lsb -size 1354x'',I5.5,'' gray:'',a,'' -equalize -depth 16 '',a)') nlines(ib), 
     *            outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.dat', 
     *            outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.png'
            write (*,*) command(1:lnblnk(command))
            call system(command(1:lnblnk(command)))
          enddo
c-- Create PNG's for band03 through band07
          do ib = 3, 7
            write (cband,'(I2.2)') ib
            write (command,'(''convert -depth 16 -endian lsb -size 2708x'',I5.5,'' gray:'',a,'' -equalize -depth 16 '',a)') nlines(ib), 
     *            outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.dat', 
     *            outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.png'
            write (*,*) command(1:lnblnk(command))
            call system(command(1:lnblnk(command)))
          enddo
c-- Create PNG's for band01 and band02
          do ib = 1, 2
            write (cband,'(I2.2)') ib
            write (command,'(''convert -depth 16 -endian lsb -size 5416x'',I5.5,'' gray:'',a,'' -equalize -depth 16 '',a)') nlines(ib), 
     *            outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.dat', 
     *            outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.png'
            write (*,*) command(1:lnblnk(command))
            call system(command(1:lnblnk(command)))
          enddo
c-- Create a quick-and-dirty earth curvature corrected RGB=221 image here ?
          if (swpng) then
            do i = 1, 3
              call get_lun(lunband(i))
            enddo
            do i = 1, 2
              nrecl(i) = 5416
              write (cband,'(I2.2)') i
              open (unit=lunband(i),file=outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.dat', access='direct',form='unformatted',recl=nrecl(i)*2)
            enddo
            call correct_init(710.0,110.0, 5416, index_correct, 10000, ncorrect_pix)
            open (unit=lunband(3),file=outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-221.rgb',form='unformatted', access='direct',recl=2*ncorrect_pix*6)
            close (unit=lunband(3),status='delete', iostat=ios)
            open (unit=lunband(3),file=outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-221.rgb',form='unformatted', access='direct',recl=2*ncorrect_pix*6)
            do i = 1, min(nlines(1), nlines(2))
              do j = 1, 2
                call read_buffer(lunband(j), i, bwbuf(1,j), nrecl(j))
              enddo
              do j = 1, 5416
                rgbbuf(1,j) = bwbuf(j,2)
                rgbbuf(2,j) = bwbuf(j,2)
                rgbbuf(3,j) = bwbuf(j,1)
              enddo
              call correct_apply(rgbbuf, 3, 5416, rgbbuf_correct, 3, 2 * ncorrect_pix, index_correct, ncorrect_pix, lunband(3), i)
            enddo
            do i = 1, 3
              close (unit=lunband(i))
              call free_lun(lunband(i))
            enddo
            write (command,'(''convert -depth 16 -endian lsb -size '',I5.5,''x'',I5.5,'' rgb:'',a,'' -equalize -depth 16 '',a)') 2*ncorrect_pix, min(nlines(1),nlines(2)), 
     *            outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-221.rgb', 
     *            outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-221'//'.png'
            write (*,*) command(1:lnblnk(command))
            call system(command(1:lnblnk(command)))
          endif
c-- End of quick-and-dirty earth curvature corrected RGB=221 image here
        endif
      else
c-- Night data processing
        i_tref = i_tref_init
        do while (.true.)
          ios  = 0
          read (lun, rec=irec, iostat=ios) inbuf_night
          if (ios.ne.0) goto 8
          if (ishft(iand(inbuf_night(16),'80'x),-7).eq.0) then
            packet_in_group = ishft(iand(inbuf_night(3),'C0'x),-6)
            packet_nr       =       iand(inbuf_night(3),'3F'x) * 256
            packet_nr       = packet_nr + inbuf_night(4)
            if (inbuf_night(4).lt.0) packet_nr = packet_nr + 256
            frame_count = iand(inbuf_night(16),'7F'x)
            frame_count = frame_count * 16 + ishft(iand(inbuf_night(17),'F0'x),-4)
            pkt_type    = ishft(iand(inbuf_night(15),'70'x),-4)
            mirror_side = iand(inbuf_night(15),'01'x)
            dy1958 = inbuf_night(7)
            if (inbuf_night(7).lt.0) dy1958 = dy1958 + 256
            dy1958 = dy1958 * 256
            dy1958 = dy1958 + inbuf_night(8)
            if (inbuf_night(8).lt.0) dy1958 = dy1958 + 256
            call dy1958_to_doy(dy1958, doy, iyr)
            call MODIS_time(inbuf_night(9),6,timestamp)
c-- Timestamp and doy go (modified) to the time buffer for the georeferencing SGP4 TLE calc
c            if (t_min.le.timestamp.and.timestamp.le.t_max.and.iyr.eq.i_yr.and.doy.eq.i_doy.and.pkt_type.eq.1) then
            if (irec_valid(irec).eq.0) then
              if (swdebug) write (*,'(''Packet Count : '',I7,1x,''Packet # in group : '',I2,1x,''Mirror_side : '',I10,1x,''Frame Count : '',I10,'' Year, doy, time : '',2I7,F15.6,I15,F15.6,I10)') 
     *          packet_nr, packet_in_group, mirror_side + 1, frame_count, iyr, doy, timestamp, dy1958, timestamp - timestamp_prev, irec_valid(irec)
              timestamp_prev = timestamp
c-- Check the timestamp to see whether a 'flush' is needed.
              i = int(timestamp*100.0D0) + 1
              if (msechist(i).gt.int(i_max/10).and.i.gt.i_tref) then
c-- Before we flush try to capture the wrong times that slip through - this could be done better by pre-building a tble to flag out of sequence times
                if (swdebug) write (*,*) i, i_tref, nint((float(i)-float(i_tref))/147.7)
                if (nint((float(i)-float(i_tref))/147.7).lt.0.or.nint((float(i)-float(i_tref))/147.7).gt.maxlengap) goto 9
c-- Flush the required empty buffers
                if (nint((float(i)-float(i_tref))/147.7).gt.1) then
                  if (swdebug) write (*,'(''Flushing'',I6,'' empty buffers'')') nint((float(i)-float(i_tref))/147.7) - 1
                  do l = 1, nint((float(i)-float(i_tref))/147.7) - 1
                    doy_zero       = doy
c-- Write band20 through band36
                    do ib = 20,36
                      theta_zero1 = 0.0D0
                      theta_zero2 = 0.0D0
                      timestamp_zero = dble(i_tref)/100.0D0 + dble (l) * 1.477815
                      do k = nlines(ib)+1, nlines(ib) + 10
                        call write_buffer_i1(luntim(ib), k, timbuf_zero, 28)
                        timestamp_zero = timestamp_zero + (1.477815D0/10.0D0)
                      enddo
                      nlines(ib) = nlines(ib) + 10
                      irecband   = nlines(ib)
                      do i = 1, 10
                        do j = 1,1
                          k = ((i-1) * 1 + j - 1) * 1354 + 1
                          call write_buffer(lunband(ib), irecband, buffer_zero_1km(k,ib-7), nrecl(ib))
                          irecband = irecband - 1
                        enddo
                      enddo
                    enddo
                  enddo
                endif
c-- End flushing empty buffers, now flush the real one
c-- Write band20 through band36
                do ib = 20,36
                  if (swcorrect) call modis_hist_correct(buffer_1km(1,ib-7), 13540, ib, swfitlin)
                  timestamp_zero = timestamp
                  do k = nlines(ib) + 1, nlines(ib) + 10
                    theta1 = (theta_pk_pk1/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk1 /  9.0D0)
                    theta2 = (theta_pk_pk2/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk2 /  9.0D0)
                    call write_buffer_i1(luntim(ib), k, timbuf, 28)
                    timestamp = timestamp + (1.477815D0/10.0D0)
                  enddo
                  timestamp = timestamp_zero
                  nlines(ib) = nlines(ib) + 10
                  irecband   = nlines(ib)
                  do i = 1, 10
                    do j = 1,1
                      k = ((i-1) * 1 + j - 1) * 1354 + 1
                      call write_buffer(lunband(ib), irecband, buffer_1km(k,ib-7), nrecl(ib))
                      irecband = irecband - 1
                    enddo
                  enddo
                enddo
                do j = 1,31
                  do i = 1, 13540
                    buffer_1km(i,j) = 0
                  enddo
                enddo
              endif
              i = int(timestamp*100.0D0) + 1
              if (msechist(i).gt.int(i_max/10).and.(nint((float(i)-float(i_tref))/147.7).gt.0)) i_tref = i
              call band_image_night(inbuf_night, 276, packet_in_group,20, buffer_1km(1,13) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,21, buffer_1km(1,14) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,22, buffer_1km(1,15) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,23, buffer_1km(1,16) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,24, buffer_1km(1,17) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,25, buffer_1km(1,18) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,26, buffer_1km(1,19) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,27, buffer_1km(1,20) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,28, buffer_1km(1,21) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,29, buffer_1km(1,22) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,30, buffer_1km(1,23) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,31, buffer_1km(1,24) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,32, buffer_1km(1,25) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,33, buffer_1km(1,26) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,34, buffer_1km(1,27) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,35, buffer_1km(1,28) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
              call band_image_night(inbuf_night, 276, packet_in_group,36, buffer_1km(1,29) , 1, 1354, 1, 10, mirror_side + 1, frame_count)
            endif
          endif
 9        irec = irec + 1
        enddo
 8      irec = irec - 1
        if (swdebug) write (*,*) irec
c-- Write the pending buffer
c-- Write band20 through band36
        do ib = 20, 36
          if (swcorrect) call modis_hist_correct(buffer_1km(1,ib-7), 13540, ib, swfitlin)
          timestamp_zero = timestamp
          do k = nlines(ib) + 1, nlines(ib) + 10
            theta1 = (theta_pk_pk1/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk1 /  9.0D0)
            theta2 = (theta_pk_pk2/2.0D0) - dble(k - nlines(ib) - 1) * (theta_pk_pk2 /  9.0D0)
            call write_buffer_i1(luntim(ib), k, timbuf, 28)
            timestamp = timestamp + (1.477815D0/10.0D0)
          enddo
          timestamp = timestamp_zero
          nlines(ib) = nlines(ib) + 10
          irecband   = nlines(ib)
          do i = 1, 10
            do j = 1,1
              k = ((i-1) * 1 + j - 1) * 1354 + 1
              call write_buffer(lunband(ib), irecband, buffer_1km(k,ib-7), nrecl(ib))
              irecband = irecband - 1
            enddo
          enddo
        enddo
        close (unit=lun)
        call free_lun(lun)
        do i = 20, 36
          close (unit=lunband(i))
          call free_lun(lunband(i))
          write (*,'('' # Lines in image band  '',I2.2,'' : '',I10)') i, nlines(i)
          close (unit=luntim(i))
          call free_lun(luntim(i))
        enddo
        if (swpng) then
c-- Create PNG's for band20 through band36
          do ib = 20, 36
            write (cband,'(I2.2)') ib
            write (command,'(''convert -depth 16 -endian lsb -size 1354x'',I5.5,'' gray:'',a,'' -equalize -depth 16 '',a)') nlines(ib), 
     *            outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.dat', 
     *            outstring(1:lnblnk(outstring))//argstring(i_dir+1:i_ext-1)//'-band'//cband//'.png'
            write (*,*) command(1:lnblnk(command))
            call system(command(1:lnblnk(command)))
          enddo
        endif
      endif
      stop
      end

      subroutine MODIS_time(buf,n,timestamp)
      implicit none
      integer*4 n
      integer*1 buf(n)
      real*8    timestamp
c
      integer*4 i, b
      integer*8 multiplier, time_msec, time_musec
c
      multiplier = 1
      time_msec  = 0
      time_musec = 0
      do i = 4, 1, -1
       b = buf(i)
       if (b.lt.0) b = b + 256
       time_msec  = time_msec + b * multiplier
       multiplier = multiplier * 256
      enddo
      multiplier = 1
      do i = 6, 5, -1
       b = buf(i)
       if (b.lt.0) b = b + 256
       time_musec  = time_musec + b * multiplier
       multiplier = multiplier * 256
      enddo
      timestamp = dble(time_msec) / 1000.0D0 + dble(time_musec)/1000000.0D0
c
      return
      end

c-- Convert days since 1-Jan-1958 to day in current year
      subroutine doy_to_dy1958(dy1958, doy, iyr)
      implicit none
      integer*4 dy1958, doy, iyr
c
      integer*4 i_day, i_yr
c
      dy1958 = 0
      do i_yr = 1958, iyr - 1
        if (mod(i_yr,4).eq.0.and.i_yr.ne.2000) then
          dy1958 = dy1958 + 366
        else
          dy1958 = dy1958 + 365
        endif
      enddo
      dy1958 = dy1958 + doy
      return
      end

c-- Convert days since 1-Jan-1958 to day in current year
      subroutine dy1958_to_doy(dy1958, doy, iyr)
      implicit none
      integer*4 dy1958, doy, iyr
c
      integer*4 i_day, i_yr
c
      i_day = 0
      i_yr  = 1958
      do while (i_day.lt.dy1958)
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
      iyr = i_yr
      doy = dy1958 - i_day
      return
      end

      subroutine band_image(inbuf, n, packet_in_group, nband, outbuf, nsamp, nfd, ndet, nfov, imir, ifd) 
      implicit none
      integer*4 n, packet_in_group, nband, nsamp, nfd, ndet, nfov, imir, ifd
      integer*2 outbuf(nsamp, nfd, ndet, nfov)
      integer*1 inbuf(n)
c
      integer*2 dn
      integer*4 isamp, idet, ifov, i_off, i_byte
      integer*4 index01(5), index03(5), index08(5)
      data index01/ 144,1140,2136,3132,4128/
      data index03/ 528,1524,2520,3516,4512/
      data index08/ 768,1764,2760,3756,4752/
c-- Sanity checks
      if (imir.eq.0.or.imir.eq.3) return
      if (ifd.lt.1.or.ifd.gt.1354) return
      if (packet_in_group.eq.0.or.packet_in_group.eq.3) return
c
c-- band01 thru band02 imagery
c
      if (nband.eq.1.or.nband.eq.2) then
        if (packet_in_group.eq.1) then
          do ifov = 1,5
            i_off = index01(ifov) + (nband - 1) * 16 * 12
            do isamp = 1,4
              do idet = 1,4
                i_byte = int(i_off / 8) + 1
                if (mod(i_off,8).eq.0) then
                  dn = inbuf(i_byte)
                  if (dn.lt.0) dn = dn + 256
                  dn = dn * 16 + ishft(iand(inbuf(i_byte+1),'F0'x),-4)
                else if (mod(i_off,8).eq.4) then
                  dn = iand(inbuf(i_byte),'0F'x) * 256
                  dn = dn + inbuf(i_byte + 1)
                  if (inbuf(i_byte + 1).lt.0) dn = dn + 256
                endif
                outbuf(isamp, ifd, idet, ifov) = dn
                i_off = i_off + 12
              enddo
            enddo
          enddo
        else if (packet_in_group.eq.2) then
          do ifov = 6,10
            i_off = index01(ifov-5) + (nband - 1) * 16 * 12
            do isamp = 1,4
              do idet = 1,4
                i_byte = int(i_off / 8) + 1
                if (mod(i_off,8).eq.0) then
                  dn = inbuf(i_byte)
                  if (dn.lt.0) dn = dn + 256
                  dn = dn * 16 + ishft(iand(inbuf(i_byte+1),'F0'x),-4)
                else if (mod(i_off,8).eq.4) then
                  dn = iand(inbuf(i_byte),'0F'x) * 256
                  dn = dn + inbuf(i_byte + 1)
                  if (inbuf(i_byte + 1).lt.0) dn = dn + 256
                endif
                outbuf(isamp, ifd, idet, ifov) = dn
                i_off = i_off + 12
              enddo
            enddo
          enddo
        endif
      endif
c
c-- band03 thru band07 imagery
c
      if (3.le.nband.and.nband.le.7) then
        if (packet_in_group.eq.1) then
          do ifov = 1,5
            i_off = index03(ifov) + (nband - 3) * 4 * 12
            do isamp = 1,2
              do idet = 1,2
                i_byte = int(i_off / 8) + 1
                if (mod(i_off,8).eq.0) then
                  dn = inbuf(i_byte)
                  if (dn.lt.0) dn = dn + 256
                  dn = dn * 16 + ishft(iand(inbuf(i_byte+1),'F0'x),-4)
                else if (mod(i_off,8).eq.4) then
                  dn = iand(inbuf(i_byte),'0F'x) * 256
                  dn = dn + inbuf(i_byte + 1)
                  if (inbuf(i_byte + 1).lt.0) dn = dn + 256
                endif
                outbuf(isamp, ifd, idet, ifov) = dn
                i_off = i_off + 12
              enddo
            enddo
          enddo
        else if (packet_in_group.eq.2) then
          do ifov = 6,10
            i_off = index03(ifov-5) + (nband - 3) * 4 * 12
            do isamp = 1,2
              do idet = 1,2
                i_byte = int(i_off / 8) + 1
                if (mod(i_off,8).eq.0) then
                  dn = inbuf(i_byte)
                  if (dn.lt.0) dn = dn + 256
                  dn = dn * 16 + ishft(iand(inbuf(i_byte+1),'F0'x),-4)
                else if (mod(i_off,8).eq.4) then
                  dn = iand(inbuf(i_byte),'0F'x) * 256
                  dn = dn + inbuf(i_byte + 1)
                  if (inbuf(i_byte + 1).lt.0) dn = dn + 256
                endif
                outbuf(isamp, ifd, idet, ifov) = dn
                i_off = i_off + 12
              enddo
            enddo
          enddo
        endif
      endif
c
c-- band08 thru band38 imagery
c
      if (8.le.nband.and.nband.le.38) then
        if (packet_in_group.eq.1) then
          do ifov = 1,5
            i_off = index08(ifov) + (nband - 8) * 1 * 12
            do isamp = 1,1
              do idet = 1,1
                i_byte = int(i_off / 8) + 1
                if (mod(i_off,8).eq.0) then
                  dn = inbuf(i_byte)
                  if (dn.lt.0) dn = dn + 256
                  dn = dn * 16 + ishft(iand(inbuf(i_byte+1),'F0'x),-4)
                else if (mod(i_off,8).eq.4) then
                  dn = iand(inbuf(i_byte),'0F'x) * 256
                  dn = dn + inbuf(i_byte + 1)
                  if (inbuf(i_byte + 1).lt.0) dn = dn + 256
                endif
                outbuf(isamp, ifd, idet, ifov) = dn
                i_off = i_off + 12
              enddo
            enddo
          enddo
        else if (packet_in_group.eq.2) then
          do ifov = 6,10
            i_off = index08(ifov-5) + (nband - 8) * 1 * 12
            do isamp = 1,1
              do idet = 1,1
                i_byte = int(i_off / 8) + 1
                if (mod(i_off,8).eq.0) then
                  dn = inbuf(i_byte)
                  if (dn.lt.0) dn = dn + 256
                  dn = dn * 16 + ishft(iand(inbuf(i_byte+1),'F0'x),-4)
                else if (mod(i_off,8).eq.4) then
                  dn = iand(inbuf(i_byte),'0F'x) * 256
                  dn = dn + inbuf(i_byte + 1)
                  if (inbuf(i_byte + 1).lt.0) dn = dn + 256
                endif
                outbuf(isamp, ifd, idet, ifov) = dn
                i_off = i_off + 12
              enddo
            enddo
          enddo
        endif
      endif
      return
      end

c-- Process the night TM - only 20 <= nband <= 36 (in this case these are the true channel nrs, not the channel nrs + 2 as for the daytime)
c-- I assume packet in group for nighttime is 3 and ignore it here ! (to prevent bit errors from causing unnecessary gaps)
      subroutine band_image_night(inbuf, n, packet_in_group, nband, outbuf, nsamp, nfd, ndet, nfov, imir, ifd) 
      implicit none
      integer*4 n, packet_in_group, nband, nsamp, nfd, ndet, nfov, imir, ifd
      integer*2 outbuf(nsamp, nfd, ndet, nfov)
      integer*1 inbuf(n)
c
      integer*2 dn
      integer*4 isamp, idet, ifov, i_off, i_byte
      integer*4 index08_1(5), index08_2(5)
      data index08_1/ 144, 348, 552, 756, 980/
      data index08_2/1164,1368,1572,1776,1980/
c-- Sanity checks
      if (imir.eq.0.or.imir.eq.3) return
      if (ifd.lt.1.or.ifd.gt.1354) return
c
c-- band08 thru band38 imagery
c
      if (20.le.nband.and.nband.le.36) then
        do ifov = 1,5
          i_off = index08_1(ifov) + (nband - 20) * 1 * 12
          do isamp = 1,1
            do idet = 1,1
              i_byte = int(i_off / 8) + 1
              if (mod(i_off,8).eq.0) then
                dn = inbuf(i_byte)
                if (dn.lt.0) dn = dn + 256
                dn = dn * 16 + ishft(iand(inbuf(i_byte+1),'F0'x),-4)
              else if (mod(i_off,8).eq.4) then
                dn = iand(inbuf(i_byte),'0F'x) * 256
                dn = dn + inbuf(i_byte + 1)
                if (inbuf(i_byte + 1).lt.0) dn = dn + 256
              endif
              outbuf(isamp, ifd, idet, ifov) = dn
              i_off = i_off + 12
            enddo
          enddo
        enddo
        do ifov = 6,10
          i_off = index08_2(ifov-5) + (nband - 20) * 1 * 12
          do isamp = 1,1
            do idet = 1,1
              i_byte = int(i_off / 8) + 1
              if (mod(i_off,8).eq.0) then
                dn = inbuf(i_byte)
                if (dn.lt.0) dn = dn + 256
                dn = dn * 16 + ishft(iand(inbuf(i_byte+1),'F0'x),-4)
              else if (mod(i_off,8).eq.4) then
                dn = iand(inbuf(i_byte),'0F'x) * 256
                dn = dn + inbuf(i_byte + 1)
                if (inbuf(i_byte + 1).lt.0) dn = dn + 256
              endif
              outbuf(isamp, ifd, idet, ifov) = dn
              i_off = i_off + 12
            enddo
          enddo
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

      subroutine write_buffer(lun,irec,array,n)
      implicit none
      integer*4 n, irec, lun
      integer*2 array(n)
c
      write (lun,rec=irec) array
      return
      end

      subroutine read_buffer(lun,irec,array,n)
      implicit none
      integer*4 n, irec, lun
      integer*2 array(n)
c
      read (lun,rec=irec) array
      return
      end

      subroutine write_buffer_i1(lun,irec,array,n)
      implicit none
      integer*4 n, irec, lun
      integer*1 array(n)
c
      write (lun,rec=irec) array
      return
      end

      subroutine modis_hist_init()
      implicit none
      integer*4 nbuf, band, ifov, idet
      integer*2 buffer(nbuf)
      logical swfitlin, swterra
c
      real*4    sum, histmin, histmax, dmaxmin, dhist
c-- The initial index 4 is to try and fix the apparent 4 pixel patter in Aqua and Terra band 1/2 images
      integer*4 nhist_250m(4, 40, 2), nhist_500m(2, 20, 5), nhist_1km(10, 31)
      real*4     hist_250m(4, 40, 2),  hist_500m(2, 20, 5),  hist_1km(10, 31)
      integer*4 i, j, k, l, lp
      logical   swnight
c-- See if this is permanently needed
      integer*4 hist_250m_full(4096, 40, 2),  hist_500m_full(4096, 20, 5),  hist_1km_full(4096, 10, 31), ih, ipix
      character*100 filenm, outstring, filinf
c
      save
c
      call getenv('HRPTOUT',outstring)
      do i = 1,2
        do j = 1,40
          do k = 1, 4
             hist_250m(k,j,i) = 0.
            nhist_250m(k,j,i) = 0
          enddo
          do k = 1, 4096
            hist_250m_full(k,j,i) = 0
          enddo
        enddo
      enddo
      do i = 1,5
        do j = 1,20
          do k = 1, 2
             hist_500m(k,j,i) = 0.
            nhist_500m(k,j,i) = 0
          enddo
          do k = 1, 4096
            hist_500m_full(k,j,i) = 0
          enddo
        enddo
      enddo
      do i = 1,31
        do j = 1,10
           hist_1km(j,i) = 0.
          nhist_1km(j,i) = 0
          do k = 1, 4096
            hist_1km_full(k,j,i) = 0
          enddo
        enddo
      enddo
c
      return
c
      entry modis_hist(buffer, nbuf, band, swterra)
c
      if (band.eq.1.or.band.eq.2) then
c-- First determine if the vertical strip has the right "amplitude" to use ....... and then determine its sum to get a mean
        dmaxmin = 175.0
        dhist   = 500.0
        if (swterra) then
          dmaxmin = 100.0
          dhist   = 300.0
        endif
        do l = 1,5416
          lp = mod(l,4) + 1
          histmin = 4096.
          histmax = -1.
          sum     = 0.
          do k = l, 39*5416 + l, 5416
            i = ((k-1)/5416) + 1
            if (min(max(buffer(k),0),4095).gt.histmax) histmax = min(max(buffer(k),0),4095)
            if (min(max(buffer(k),0),4095).lt.histmin) histmin = min(max(buffer(k),0),4095)
            sum = sum + min(max(buffer(k),0),4095)
          enddo
          if (histmax-histmin.lt.dmaxmin.and.histmin.le.dhist) then 
            do k = l, 39*5416 + l, 5416
              i = ((k-1)/5416) + 1
               hist_250m(lp,i,band) =  hist_250m(lp,i,band) + float(min(max(buffer(k),0),4095)) - (sum/40.0)
              nhist_250m(lp,i,band) = nhist_250m(lp,i,band) + 1
            enddo
          endif
          do k = l, 39*5416 + l, 5416
            i  = ((k-1)/5416) + 1
            ih = min(max(buffer(k),0),4095) + 1
            hist_250m_full(ih, i, band) = hist_250m_full(ih, i, band) + 1 
          enddo
        enddo
      endif
      if (3.le.band.and.band.le.7) then
c-- First determine if the vertical strip has the right "amplitude" to use ....... and then determine its sum to get a mean
        dmaxmin = 100.0
        dhist   = 500.0
        if (swterra) then
          dmaxmin = 100.0
          if (band.eq.3) dmaxmin = 150.0
          dhist   = 300.0
        endif
        do l = 1,2708
          lp = mod(l,2) + 1
          histmin = 4096.
          histmax = -1.
          sum     = 0.
          do k = l, 19*2708 + l, 2708
            i = ((k-1)/2708) + 1
            if (min(max(buffer(k),0),4095).gt.histmax) histmax = min(max(buffer(k),0),4095)
            if (min(max(buffer(k),0),4095).lt.histmin) histmin = min(max(buffer(k),0),4095)
            sum = sum + min(max(buffer(k),0),4095)
          enddo
          if (histmax-histmin.lt.dmaxmin.and.histmin.le.dhist) then 
            do k = l, 19*2708 + l, 2708
              i = ((k-1)/2708) + 1
               hist_500m(lp,i,band-2) =  hist_500m(lp,i,band-2) + float(min(max(buffer(k),0),4095)) - (sum/20.0)
              nhist_500m(lp,i,band-2) = nhist_500m(lp,i,band-2) + 1
            enddo
          endif
          do k = l, 19*2708 + l, 2708
            i  = ((k-1)/2708) + 1
            ih = min(max(buffer(k),0),4095) + 1
            hist_500m_full(ih, i, band-2) = hist_500m_full(ih, i, band-2) + 1 
          enddo
        enddo
      endif
      if (8.le.band.and.band.le.38) then
c-- First determine if the vertical strip has the right "amplitude" to use ....... and then determine its sum to get a mean
        do l = 1,1354
          histmin = 4096.
          histmax = -1.
          sum     = 0.
          do k = l, 9*1354 + l, 1354
            i = ((k-1)/1354) + 1
            if (min(max(buffer(k),0),4095).gt.histmax) histmax = min(max(buffer(k),0),4095)
            if (min(max(buffer(k),0),4095).lt.histmin) histmin = min(max(buffer(k),0),4095)
            sum = sum + min(max(buffer(k),0),4095)
          enddo
          if (histmax-histmin.lt.100.0) then 
            do k = l, 9*1354 + l, 1354
              i = ((k-1)/1354) + 1
               hist_1km(i,band-7) =  hist_1km(i,band-7) + float(min(max(buffer(k),0),4095)) - (sum/10.0)
              nhist_1km(i,band-7) = nhist_1km(i,band-7) + 1
            enddo
          endif
          do k = l, 9*1354 + l, 1354
            i  = ((k-1)/1354) + 1
            ih = min(max(buffer(k),0),4095) + 1
            hist_1km_full(ih, i, band-7) = hist_1km_full(ih, i, band-7) + 1 
          enddo
        enddo
      endif
      return
c
      entry modis_hist_correct(buffer, nbuf, band, swfitlin)
c
      if (swfitlin) then
        if (band.eq.1.or.band.eq.2) then
          do i = 1, 10
            do j = 1,4
              do k = ((i-1) * 4 + j - 1) * 5416 + 1, ((i-1) * 4 + j - 1) * 5416 + 5416
                lp = mod(k,4) + 1
                buffer(k) = buffer(k) - nint(hist_250m(lp,(i-1)*4+j,band))
              enddo
            enddo
          enddo
        endif
        if (3.le.band.and.band.le.7) then
          do i = 1, 10
            do j = 1,2
              do k = ((i-1) * 2 + j - 1) * 2708 + 1, ((i-1) * 2 + j - 1) * 2708 + 2708
                lp = mod(k,2) + 1
                buffer(k) = buffer(k) - nint(hist_500m(lp,(i-1)*2+j,band-2))
              enddo
            enddo
          enddo
        endif
        if (8.le.band.and.band.le.38) then
          do i = 1, 10
            do j = 1,1
              do k = ((i-1) * 1 + j - 1) * 1354 + 1, ((i-1) * 1 + j - 1) * 1354 + 1354
                buffer(k) = buffer(k) - nint(hist_1km((i-1)*1+j,band-7))
              enddo
            enddo
          enddo
        endif
      else
        if (band.eq.1.or.band.eq.2) then
          do i = 1, 10
            do j = 1,4
              do k = ((i-1) * 4 + j - 1) * 5416 + 1, ((i-1) * 4 + j - 1) * 5416 + 5416
                lp = mod(k,4) + 1
                buffer(k) = buffer(k) - nint(hist_250m(lp,(i-1)*4+j,band))
              enddo
            enddo
          enddo
        endif
        if (3.le.band.and.band.le.7) then
          do i = 1, 10
            do j = 1,2
              do k = ((i-1) * 2 + j - 1) * 2708 + 1, ((i-1) * 2 + j - 1) * 2708 + 2708
                lp = mod(k,2) + 1
                buffer(k) = buffer(k) - nint(hist_500m(lp,(i-1)*2+j,band-2))
              enddo
            enddo
          enddo
        endif
        if (8.le.band.and.band.le.38) then
          do i = 1, 10
            do j = 1,1
              do k = ((i-1) * 1 + j - 1) * 1354 + 1, ((i-1) * 1 + j - 1) * 1354 + 1354
                buffer(k) = buffer(k) - nint(hist_1km((i-1)*1+j,band-7))
              enddo
            enddo
          enddo
        endif
      endif
      return
c
      entry modis_hist_flush(swnight, swfitlin)
c
      if (.not.swnight) then
        do i = 1,2
          do j = 1,40
            do k = 1, 4
              if (nhist_250m(k,j,i).gt.0) then 
                hist_250m(k,j,i) = hist_250m(k,j,i) / float(nhist_250m(k,j,i))
c              write (*,*) i, j, hist_250m(j,i), nhist_250m(j,i)
              endif
            enddo
c            write (*,'(2I10,4(F10.3,1x))') i, j, (hist_250m(k,j,i),k=1,4)
c-- For the time being only plot the spectra here
            if (swfitlin) then
              if (j.ge.2) then
                write (filinf,'(a,''250m_b'',I2.2,''l''I2.2)') outstring(1:lnblnk(outstring)), i, j
                write (filenm,'(a,''250m_band_'',I2.2,''_line_''I2.2)') outstring(1:lnblnk(outstring)), i, j
                call hist_plot(hist_250m_full(1,j,i), hist_250m_full(1,1,i), 4096, 1, 4096, 1, filenm, ' ', filinf)
                call system('gmt psconvert '//filenm(1:lnblnk(filenm))//'.ps -A0.2c+s5c -Tg')
                call system('convert '//filenm(1:lnblnk(filenm))//'.png -rotate 90 '//filenm(1:lnblnk(filenm))//'.png')
              endif
            endif
c
          enddo
          if (swfitlin) call system('rm '//outstring(1:lnblnk(outstring))//'*.ps')
        enddo
        do i = 1,5
          do j = 1,20
            do k = 1, 2
              if (nhist_500m(k,j,i).gt.0) then 
                hist_500m(k,j,i) = hist_500m(k,j,i) / float(nhist_500m(k,j,i))
              endif
            enddo
c
            if (swfitlin) then
              if (j.ge.2) then
                write (filinf,'(a,''500m_b'',I2.2,''l''I2.2)') outstring(1:lnblnk(outstring)), i, j
                write (filenm,'(a,''500m_band_'',I2.2,''_line_''I2.2)') outstring(1:lnblnk(outstring)), i, j
                call hist_plot(hist_500m_full(1,j,i), hist_500m_full(1,1,i), 4096, 1, 4096, 1, filenm, ' ', filinf)
                call system('gmt psconvert '//filenm(1:lnblnk(filenm))//'.ps -A0.2c+s5c -Tg')
                call system('convert '//filenm(1:lnblnk(filenm))//'.png -rotate 90 '//filenm(1:lnblnk(filenm))//'.png')
              endif
            endif
c
          enddo
          if (swfitlin) call system('rm '//outstring(1:lnblnk(outstring))//'*.ps')
        enddo
        do i = 1,31
          do j = 1,10
            if (nhist_1km(j,i).gt.0.and.(.not.swfitlin)) then 
              hist_1km(j,i) = hist_1km(j,i) / float(nhist_1km(j,i))
c              write (*,*) i, j, hist_1km(j,i), nhist_1km(j,i)
            endif
          enddo
        enddo
      else
        do i = 13,29
          do j = 1,10
            if (nhist_1km(j,i).gt.0) then 
              hist_1km(j,i) = hist_1km(j,i) / float(nhist_1km(j,i))
c              write (*,*) i, j, hist_1km(j,i), nhist_1km(j,i)
            endif
          enddo
        enddo
      endif
c
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

      subroutine hist_plot(hist,hist_ref, nx,ny,nxin,nyin,filenm,ctxt, filinf)
      implicit none
      integer*4 nx, ny, hist(nx, ny), nxin, nyin, hist_ref(nx, ny)
      character*(*) filenm, ctxt, filinf
c
      integer*4 lunplot, i, j
      real*4 xr(2), yr(2), rhisty(4096), rhistx(4096)
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
c-- Fix the vertical scaling
      yr(2) = 99999.0
      j = len(filenm)
      do while (filenm(j:j).ne.'/')
        j = j - 1
      enddo
      call pkg_frame(11,-6,1.,xr,yr,'Bin #','Frequency',filinf(j+1:lnblnk(filinf)))
      do j = 1,nyin
        do i = 1,nxin
          rhisty(i) = float(hist(i,j))
          rhistx(i) = float(i)
        enddo
        call newpen(7-j)
        call pkg_plhist(rhistx,rhisty,nxin)
        do i = 1,nxin
          rhisty(i) = float(hist_ref(i,j))
          rhistx(i) = float(i)
        enddo
        call newpen(6-j)
        call pkg_plhist(rhistx,rhisty,nxin)
      enddo
      call newpen(6)
      call pkg_pltextbl(ctxt(1:lnblnk(ctxt)),1)
      call pkg_clospl()
      call free_lun(lunplot)
      return
      end

      subroutine get_TM(filenm,NCTRS_TM,n_max,n_actual,n_req)
      implicit none
c
      integer*4 n_max, n_actual, n_req
      integer*1 NCTRS_TM(n_max)
      character*(*) filenm
c
      integer*1 tm_buffer(1024)
c
      integer*4 iptr, irec, lun, n_used, i, nbytes, eof
c
      save iptr, tm_buffer, irec, lun, nbytes, eof
c
      n_used = 0
      if (n_actual.eq.-1) then
        call get_lun(lun)
        irec = 0
        open (unit=lun,file=filenm,status='old',form='unformatted',access='direct',recl=1024)
        inquire (file=filenm, size=nbytes)
        iptr = 0
        eof  = 0
        call get_file_buffer_i1(lun,1024,tm_buffer,irec,iptr)
      endif
      do i = 1, n_req
        iptr = iptr + 1
        if (iptr.gt.1024) then
          call get_file_buffer_i1(lun,1024,tm_buffer,irec,iptr)
          if (irec*1024.gt.nbytes) then
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
          if (iptr+(irec-1)*1024.gt.nbytes) then
            eof = 0
            close (unit=lun)
            call free_lun(lun)
            n_actual = -1
            return
          endif
        endif
        n_used = n_used + 1
        NCTRS_TM(n_used) = tm_buffer(iptr)
      enddo
      n_actual = n_used
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

      subroutine correct_init(sat_avg_alt, sat_fov, nsat_pix, index_correct, nmax, ncorrect_pix)
      implicit none
      integer*4 nsat_pix, ncorrect_pix, i, nmax
      integer*2 index_correct(nmax)
      real*4 sat_avg_alt, sat_fov
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
