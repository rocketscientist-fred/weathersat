      implicit none
      integer*4 nchan   , npix     , nrscans   , ngrey     , nfy
      parameter (nchan=5, npix=8192, nrscans=40, ngrey=4096, nfy=12390)
c
      logical   swcorrect, swdebug, swzero(3), swzerorec, swfix
      integer*1 filebuf(1024), nextbuf(1024), highres(nfy), shiftbuf(3), mask1(8), mask2(8), alignbuf(nfy), imir_buf(nrscans, nchan)
      integer*2 pixbuf(npix), channel_image(npix, nrscans, nchan), zero(npix), ival, offset(nchan, nrscans, 2)
      integer*2 rgbbuf(3, npix), rgbzero(3, npix), pixbuf2(npix)
      integer*4 n_actual, i, i1, i2, i3, j, k, l, i_sig, n_sig, iline, i_off, hist(50000), highres_offset(1000000)
      integer*4 lun(nchan), irec(nchan), iline_last(nchan), iline_last_all, i_channel, irec_max, n_wrong
      integer*4 n_actual_next, nbyte, nbit, jref, i_sync_bit, i_sync_byte, sync_hist(8), lunrgb, imir
      integer*4 chshift(nchan), lunout(nchan), iline_store(1000000), i_zero(3), dy2k, doy
      integer*4 imir_store(1000000), imir_rec(1000000, nchan), nrlines, imir_sum, imir_act
c
      character*2 cband
      character*5 cnpix
      character*500 filenm, command, argstring
c
      data mask1/'FF'x,'7F'x,'3F'x,'1F'x,'0F'x,'07'x,'03'x,'01'x/
      data mask2/'00'x,'80'x,'C0'x,'E0'x,'F0'x,'F8'x,'FC'x,'FE'x/
      data chshift/-8,0,-8,16,0/
c
      write (cnpix,'(I4.4,A1)') npix,'x'
      do i = 1, 50000
        hist(i) = 0
      enddo
      do i = 1, 1000000
        highres_offset(i) = 0
        iline_store(i)    = 0
        imir_store(i)     = 0
        do j = 1, nchan
          imir_rec(i,j) = 0
        enddo
      enddo
      do i = 1, nchan
        do j = 1, 40
          do k = 1, npix
            channel_image(k, j, i) = 0
          enddo
          imir_buf(j, i) = 0
        enddo
        iline_last(i) = -1
        irec(i)       =  0
        lun(i)        =  0
      enddo
      iline_last_all = -1
      do i = 1, npix
        zero(i) = 0
      enddo
      n_wrong = 0
      do i = 1, 8
        sync_hist(i) = 0
      enddo
      do i = 1, npix
        do j = 1, 3
          rgbzero(j, i) = 0
        enddo
      enddo
c-- Process all the command line arguments and take initialisation actions - if required
      call getarg(1, filenm)
      swcorrect = .false.
      swdebug   = .false.
      swfix     = .false.
      do i = 2,30
        call getarg(i, argstring)
        if (index(argstring(1:lnblnk(argstring)),'correct').ne.0) swcorrect = .true.
        if (index(argstring(1:lnblnk(argstring)),'debug').ne.0)   swdebug   = .true.
        if (index(argstring(1:lnblnk(argstring)),'fix').ne.0)     swfix     = .true.
      enddo
      if (swcorrect) call hist_correct_init(npix, nrscans, nchan, 2, ngrey)
c
c-- Find the sync marker alignment !
c
      n_actual = -1
      n_sig    = 0
      i_sig    = 1
      j        = 0
      jref     = 0
      k        = 0
      call get_file_bytes(filenm, filebuf, 1024, n_actual, 1024)
      do while (n_actual.gt.0)
        call get_file_bytes(filenm, nextbuf, 1024, n_actual_next, 1024)
        call  find_fy_mask(filebuf, nextbuf, 1024, nbyte, nbit)
        if (nbyte.gt.0) then
          if (j + nbyte - jref.le.50000) then
            hist(j + nbyte - jref) = hist(j + nbyte - jref) + 1
          endif
c-- High Res records
          if (j+nbyte-jref.eq.nfy) then
            k = k + 1
            highres_offset(k) = jref
            sync_hist(nbit + 1) = sync_hist(nbit + 1) + 1
          endif
c-- High Res records of double length - so damaged signature in the middle 
c-- For Alan's mersi.bin files this appears to point to garbage ..... so removed
c          if (j+nbyte-jref.eq.2*nfy) then
c            k = k + 1
c            highres_offset(k) = jref
c            k = k + 1
c            highres_offset(k) = jref + nfy
c          endif
          jref     = j + nbyte
        endif
        j        = j + n_actual
        n_actual = n_actual_next
        if (n_actual.gt.0) then
          do i = 1, n_actual
            filebuf(i) = nextbuf(i)
          enddo
        endif
      enddo
      do i = 1, 50000
        if (hist(i).gt.0.and.swdebug) write (*,*) i, hist(i)
      enddo
      i1 = sync_hist(1)
      i2 = 1
      do i = 2,8
        if (sync_hist(i).gt.i1) then
          i1 = sync_hist(i)
          i2 = i
        endif
      enddo
c
c-- It turns out that the sync marker I use ends up 1 bit below the start of the sync align in the data
c-- So where the sync marker akigns (nbit) is 1 bit below what value to use :-)
c
      i_sync_bit = i2 - 1
      write (*,'(''Sync marker aligned on bit (MSB=0, LSB=7, L -> R) : '', I3)') i_sync_bit
c
c-- End of sync marker search
c
c-- From inspecting the data, the MSB of the line counter comes 65 bits after the first bit of the sync marker (as I have defined it in find_fy_mask)
c
c-- For the lines below, I suspect that for i_sync_bit 4, 5, 6, 7 will need i_sync_byte = 10 and for 0, 1, 2, 3  i_sync_byte=9
      i_sync_bit = i_sync_bit + 1
      if (i_sync_bit.gt.7) then
        i_sync_bit  = i_sync_bit - 8
        i_sync_byte = 10
      else
        i_sync_byte = 9
      endif
      do i = 1, nchan
        call get_lun(lun(i))
        write (cband,'(I2.2)') i
        open (unit=lun(i),file=filenm(1:lnblnk(filenm))//'-band'//cband//'.dat',access='direct',form='unformatted',recl=16384)
        irec(i) = 0
      enddo
c-- Go through the file once to fix the line nrs out of sequence - a bit weird but could not think of another way ....
      n_actual = -2
      j        = 0
      l        = 1
      do while (n_actual.gt.0.or.n_actual.eq.-2)
        if (n_actual.eq.-2) n_actual = -1
        call get_file_bytes(filenm, filebuf, 1024, n_actual, 1024)
        if (n_actual.gt.0) then
          do i = 1, n_actual
            j = j + 1
            if (j.ge.highres_offset(l).and.j.le.highres_offset(l) + nfy - 1) then
              highres(j - highres_offset(l) + 1) = filebuf(i)
            endif
            if (j.eq.highres_offset(l)+nfy) then
              call bit_align_buffer(highres, nfy, alignbuf, nfy, i_sync_bit) 
              iline = alignbuf(9)
              if (iline.lt.0) iline = iline + 256
              imir  = ishft(iand(alignbuf(98),'F0'x),-4) + 1
c-- Alan's mersi.bin files contain mysteriously empty lines which confuse my buffer flushing, so remove artificially
              swzerorec = .false.
              if (iline.eq.0) then
                swzerorec = .true.
                i1 = 11
                do while (swzerorec.and.i1.le.100)
                  if (highres(i1).ne.0) swzerorec = .false.
                  i1 = i1 + 1
                enddo
              endif
              if (swzerorec) then
                iline_store(l) = -1
                imir_store(l)  = 1
              else
                iline = iline + 1
                iline_store(l) = iline
                imir_store(l)  = imir
                if (swdebug) write (*,'(100(z2.2,1x),4I8,1x)') (alignbuf(i1), i1 = 1,100), l, k,  iline, imir
              endif
              l = l + 1
              if (l.gt.k) goto 3
            endif
          enddo
        endif
      enddo
 3    continue
      do i = 2, k - 1
c-- fix the single wrong line
        if (iline_store(i).le.iline_store(i-1).and.iline_store(i+1).eq.iline_store(i-1)+2) then
          if (swdebug) write (*,'(''Fixed line counter : '',4I8)') iline_store(i-1), iline_store(i), iline_store(i+1), i
          iline_store(i) = iline_store(i-1) + 1
        endif
        if (iline_store(i).le.iline_store(i-1).and.iline_store(i+1).gt.iline_store(i-1)) then
          iline_store(i) = -1
        endif
        if (iline_store(i).lt.1.or.iline_store(i).gt.nchan*40) then
          iline_store(i) = -1
          n_wrong = n_wrong + 1
        endif
      enddo
c-- End of line fix
c
      n_actual = -2
      j        = 0
      l        = 1
      do while (n_actual.gt.0.or.n_actual.eq.-2)
        if (n_actual.eq.-2) n_actual = -1
        call get_file_bytes(filenm, filebuf, 1024, n_actual, 1024)
        if (n_actual.gt.0) then
          do i = 1, n_actual
            j = j + 1
            if (j.ge.highres_offset(l).and.j.le.highres_offset(l) + nfy - 1) then
              highres(j - highres_offset(l) + 1) = filebuf(i)
            endif
c-- Flush the highres buffer
            if (j.eq.highres_offset(l)+nfy) then
              iline = iline_store(l)
              imir  = imir_store(l)
              if (iline.lt.0) goto 2
              call bit_align_buffer(highres, nfy, alignbuf, nfy, i_sync_bit) 
              i_channel = ((iline - 1) / 40) + 1
              if (i_channel.gt.nchan) goto 2
              i_off = 91 + i_sync_byte
              do i2 = 1, 8191, 2
                i1 = i_off + ((i2 - 1) / 2) * 3
                shiftbuf(1) = alignbuf(i1  )
                shiftbuf(2) = alignbuf(i1+1)
                shiftbuf(3) = alignbuf(i1+2)
                ival = shiftbuf(1)
                if (ival.lt.0) ival = ival + 256
                ival = ival * 16
                ival = ival + ishft(iand(shiftbuf(2),'F0'x),-4)
                pixbuf(i2) = ival
                ival = iand(shiftbuf(2),'0F'x)
                ival = ival * 256
                ival = ival + shiftbuf(3)
                if (shiftbuf(3).lt.0) then
                  ival = ival + 256
                endif
                pixbuf(i2 + 1) = ival
              enddo
c-- The pixel buffer has been filled - store it in the channel_image buffer
c
c-- Check if there is a global wrap around. If so flush empty buffers to keep the 5 channel images sync'd
              if (iline.le.iline_last_all) then
                irec_max = irec(1)
                do i1 = 2, nchan
                  if (irec(i1).gt.irec_max) irec_max = irec(i1)
                enddo
                do i1 = 1, nchan
                  if (irec(i1).lt.irec_max) then
                    if (swdebug) write (*,'(''Flush empty ! '',3I8)') i1, irec(i1), irec_max
                    do i2 = 1, irec_max - irec(i1)
                      irec(i1) = irec(i1) + 1
                      call write_buffer(lun(i1), irec(i1), zero, npix)
                    enddo
                  endif
                enddo
c-- and then flush and clean the current buffers
                do i1 = 1, nchan
                  nrlines  = 0
                  imir_sum = 0
                  do i2 = 1, 40
                    if (imir_buf(i2, i1).ne.-1) then
                      nrlines  = nrlines + 1
                      imir_sum = imir_sum + imir_buf(i2, i1)
                    endif
                  enddo
                  imir_act = nint(float(imir_sum)/float(nrlines))
                  do i2 = 1, 40
                    irec(i1) = irec(i1) + 1
                    imir_rec(irec(i1), i1) = imir_act
                    call write_buffer(lun(i1), irec(i1), channel_image(1,i2,i1), npix)
c                    if (swcorrect) call hist_create(channel_image(1,i2,i1), npix, i2, i1)
                    do i3 = 1, npix
                      channel_image(i3,i2,i1) = 0
                    enddo
                    imir_buf(i2,i1) = -1
                  enddo
                enddo
              endif
              iline_last(i_channel) = iline
              iline_last_all        = iline
              i2 = mod(iline-1,40) + 1
              do i1 = 1, npix
                channel_image(i1,i2,i_channel) = pixbuf(i1)
              enddo
              imir_buf(i2, i_channel) = imir
c
 2            l = l + 1
              if (l.gt.k) goto 1
            endif
          enddo
        endif
      enddo
 1    continue
c-- Final sync align
      irec_max = irec(1)
      do i1 = 2, nchan
        if (irec(i1).gt.irec_max) irec_max = irec(i1)
      enddo
      do i1 = 1, nchan
        if (irec(i1).lt.irec_max) then
          do i2 = 1, irec_max - irec(i1)
            irec(i1) = irec(i1) + 1
            call write_buffer(lun(i1), irec(i1), zero, npix)
          enddo
        endif
      enddo
c-- If needed, temporarily fix Alan'1 "Line 0 problem" here
      if (swfix) then
        write (cband,'(I2.2)') 1
        open (unit=lun(1),file=filenm(1:lnblnk(filenm))//'-band'//cband//'.dat',access='direct',form='unformatted',recl=16384)
        i = 40
        do while (i.le.irec(1) - 2)
          read  (lun(1),rec=i)   pixbuf
          read  (lun(1),rec=i+2) pixbuf2
          do j = 1, npix
            pixbuf(j) = (pixbuf(j) + pixbuf2(j)) / 2
          enddo
          write (lun(1),rec=i+1) pixbuf
          i = i + 40
        enddo
      endif
c-- End of temporary fix
c-- Analyse the offset correection and add it to the channels
      if (swcorrect) then
        do i = 1, nchan
          do j = 1, irec(i)
            read (lun(i),rec=j) pixbuf
            call hist_create(pixbuf, npix, mod(j-1,40) + 1, i, imir_rec(j,i))
          enddo
        enddo
        call hist_analyse(offset, nchan, nrscans, 2, swdebug)
        do i = 1, nchan
          do j = 1, irec(i)
            imir_act = imir_rec(j, i)
            if (imir_act.eq.0) imir_act = 1
            read (lun(i),rec=j) pixbuf
            do k = 1, npix
              pixbuf(k) = min(max(pixbuf(k) + offset(i, mod(j-1,40) + 1, imir_act), 0), 4095)
            enddo
            write (lun(i),rec=j) pixbuf
          enddo
        enddo
      endif
c--
      do i = 1, nchan
        write (*,*) i, irec(i)
      enddo
      if (n_wrong.gt.0.and.swdebug) write (*,'(''Wrong line counters : '',I10)') n_wrong
c-- Convert data files to images
      do i = 1, nchan
        write (cband,'(I2.2)') i
        write (command,'(''convert -depth 16 -endian lsb -size 8192x'',I5.5,'' gray:'',a,'' -equalize -depth 16 -rotate 180 '',a)') irec(i), 
     *            filenm(1:lnblnk(filenm))//'-band'//cband//'.dat', 
     *            filenm(1:lnblnk(filenm))//'-band'//cband//'.png'
        write (*,*) command(1:lnblnk(command))
        call my_system(command(1:lnblnk(command)))
      enddo
      do i = 1, nchan
        write (cband,'(I2.2)') i
        open (unit=lun(i),file=filenm(1:lnblnk(filenm))//'-band'//cband//'.dat',access='direct',form='unformatted',recl=16384)
        call get_lun(lunout(i))
        write (cband,'(I2.2)') i
        open (unit=lunout(i),file=filenm(1:lnblnk(filenm))//'-band'//cband//'-shift.dat',access='direct',form='unformatted',recl=16384)
      enddo
      do i = 1, nchan
        do j = 1, irec(i)
          read (lun(i),rec=j) pixbuf
          if (chshift(i).ne.0) then
            if (chshift(i).gt.0) then
              do k = npix - chshift(i), 1, -1
                pixbuf(k+chshift(i)) = pixbuf(k)
              enddo
              do k = 1, chshift(i)
                pixbuf(k) = 0
              enddo
            else
              do k = abs(chshift(i)) + 1, npix
                pixbuf(k+chshift(i)) = pixbuf(k)
              enddo
              do k = npix + chshift(i) + 1, npix
                pixbuf(k) = 0
              enddo
            endif
          endif
          write (lunout(i),rec=j) pixbuf
        enddo
      enddo
      do i = 1, nchan
        close(unit=lun(i))
        call free_lun(lun(i))
        close(unit=lunout(i))
        call free_lun(lunout(i))
      enddo
      do i = 1, nchan
        write (cband,'(I2.2)') i
        write (command,'(''convert -depth 16 -endian lsb -size 8192x'',I5.5,'' gray:'',a,'' -equalize -depth 16 -rotate 180 '',a)') irec(i), 
     *            filenm(1:lnblnk(filenm))//'-band'//cband//'-shift.dat', 
     *            filenm(1:lnblnk(filenm))//'-band'//cband//'-shift.png'
        write (*,*) command(1:lnblnk(command))
        call my_system(command(1:lnblnk(command)))
      enddo
c-- Create a false colour composite from the shifted dat file for ch 3, 2 and 1 - in two different ways
      write (command,'(''convert -combine '',a,a,a,a)') 
     *            filenm(1:lnblnk(filenm))//'-band03-shift.png ', 
     *            filenm(1:lnblnk(filenm))//'-band02-shift.png ', 
     *            filenm(1:lnblnk(filenm))//'-band01-shift.png ', 
     *            filenm(1:lnblnk(filenm))//'-RGB-combine.png'
      write (*,*) command(1:lnblnk(command))
      call my_system(command(1:lnblnk(command)))
      write (command,'(''./deproject.exe '',a,a)') 
     *            filenm(1:lnblnk(filenm))//'-RGB-combine.png',' alt=848.0 sharpen bowtie=mersi1'
      write (*,*) command(1:lnblnk(command))
      call my_system(command(1:lnblnk(command)))
c
      stop
      end

      subroutine get_file_bytes(filenm,cadubuf,n_max,n_actual,n_req)
      implicit none
c
      integer*4 n_max, n_actual, n_req
      integer*1 cadubuf(n_max)
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
        cadubuf(n_used) = tm_buffer(iptr)
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

      subroutine my_system(command)
      implicit none
      character*(*) command
c
      integer*4 status
      status = 0
      call system(command(1:lnblnk(command)), status)
      if (status.ne.0) then
        write (*,'(''Call to execute : '',a,'' failed with status : '', I10)') command(1:lnblnk(command)), status
      endif
      return
      end

      subroutine FY_time(buf,n,timestamp)
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

c-- This subroutine searches for the FY3 line sync marker which is NOT necessarily aligned at byte level .... sigh
c-- So, search at bit level !
      subroutine find_fy_mask(buffer, nextbuf, n, nbyte, nbit)
      implicit none
      integer*4 n, nbyte, nbit
      integer*1 buffer(n), nextbuf(n)
c
      integer*1 mask(6), bytemask(7), masks(7), masks_reverse(7), bcheck
      integer*4 i, j, init
      logical   swmatch
c
      data init/0/
      save
c
      if (init.eq.0) then
        init    = 1
c-- The FY3 sync marker (my guess .....)
        mask(1) = '00'x
        mask(2) = '0F'x
        mask(3) = 'FF'x
        mask(4) = '00'x
        mask(5) = '0F'x
        mask(6) = 'FF'x
c-- The bit masks from LSB to MSB
        masks(1) = '01'x
        masks(2) = '03'x
        masks(3) = '07'x
        masks(4) = '0F'x
        masks(5) = '1F'x
        masks(6) = '3F'x
        masks(7) = '7F'x
c-- The bit masks from MSB to LSB
        masks_reverse(1) = '80'x
        masks_reverse(2) = 'C0'x
        masks_reverse(3) = 'E0'x
        masks_reverse(4) = 'F0'x
        masks_reverse(5) = 'F8'x
        masks_reverse(6) = 'FC'x
        masks_reverse(7) = 'FE'x
      endif
c
      nbyte = -1
      do j = 1, n
        do i = 0, 7
          if (i.eq.0) then
            bytemask(1) = mask(1)
            bytemask(2) = mask(2)
            bytemask(3) = mask(3)
            bytemask(4) = mask(4)
            bytemask(5) = mask(5)
            bytemask(6) = mask(6)
            swmatch = .true.
            if (j.le.n) then
              bcheck = buffer(j)
            else
              bcheck = nextbuf(j - n)
            endif
            if (bcheck.ne.bytemask(1)) swmatch = .false.
            if (j+1.le.n) then
              bcheck = buffer(j+1)
            else
              bcheck = nextbuf(j + 1 - n)
            endif
            if (bcheck.ne.bytemask(2)) swmatch = .false.
            if (j+2.le.n) then
              bcheck = buffer(j+2)
            else
              bcheck = nextbuf(j + 2 - n)
            endif
            if (bcheck.ne.bytemask(3)) swmatch = .false.
            if (j+3.le.n) then
              bcheck = buffer(j+3)
            else
              bcheck = nextbuf(j + 3 - n)
            endif
            if (bcheck.ne.bytemask(4)) swmatch = .false.
            if (j+4.le.n) then
              bcheck = buffer(j+4)
            else
              bcheck = nextbuf(j + 4 - n)
            endif
            if (bcheck.ne.bytemask(5)) swmatch = .false.
            if (j+5.le.n) then
              bcheck = buffer(j+5)
            else
              bcheck = nextbuf(j + 5 - n)
            endif
            if (bcheck.ne.bytemask(6)) swmatch = .false.
          else
            bytemask(1) =                                         ishft(mask(1),-i)
            bytemask(2) = ior  (ishft(iand(mask(1),masks(i)),8-i),ishft(mask(2),-i))
            bytemask(3) = ior  (ishft(iand(mask(2),masks(i)),8-i),ishft(mask(3),-i))
            bytemask(4) = ior  (ishft(iand(mask(3),masks(i)),8-i),ishft(mask(4),-i))
            bytemask(5) = ior  (ishft(iand(mask(4),masks(i)),8-i),ishft(mask(5),-i))
            bytemask(6) = ior  (ishft(iand(mask(5),masks(i)),8-i),ishft(mask(6),-i))
            bytemask(7) =       ishft(iand(mask(6),masks(i)),8-i)
            swmatch = .true.
            if (j.le.n) then
              bcheck = buffer(j)
            else
              bcheck = nextbuf(j - n)
            endif
            if (iand(bcheck,masks(8-i))      .ne.bytemask(1)) swmatch = .false.
            if (j+1.le.n) then
              bcheck = buffer(j+1)
            else
              bcheck = nextbuf(j + 1 - n)
            endif
            if (bcheck                       .ne.bytemask(2)) swmatch = .false.
            if (j+2.le.n) then
              bcheck = buffer(j+2)
            else
              bcheck = nextbuf(j + 2 - n)
            endif
            if (bcheck                       .ne.bytemask(3)) swmatch = .false.
c
            if (j+3.le.n) then
              bcheck = buffer(j+3)
            else
              bcheck = nextbuf(j + 3 - n)
            endif
            if (bcheck                       .ne.bytemask(4)) swmatch = .false.
            if (j+4.le.n) then
              bcheck = buffer(j+4)
            else
              bcheck = nextbuf(j + 4 - n)
            endif
            if (bcheck                       .ne.bytemask(5)) swmatch = .false.
            if (j+5.le.n) then
              bcheck = buffer(j+5)
            else
              bcheck = nextbuf(j + 5 - n)
            endif
            if (bcheck                       .ne.bytemask(6)) swmatch = .false.
            if (j+6.le.n) then
              bcheck = buffer(j+6)
            else
              bcheck = nextbuf(j + 6 - n)
            endif
            if (iand(bcheck,masks_reverse(i)).ne.bytemask(7)) swmatch = .false.
          endif
          if (swmatch) then
            nbyte = j
            nbit  = i
            return
          endif
        enddo
      enddo
      return
      end

      subroutine bit_align_buffer(inbuf, n, alignbuf, nalign, nbit) 
      implicit none
      integer*4 n, nalign, nbit
      integer*1 inbuf(n), alignbuf(nalign)
c
      integer*1 masks(7), masks_reverse(7)
      integer*4 i, init
      data init/0/
      save
c
      if (init.eq.0) then
        init    = 1
c-- The bit masks from LSB to MSB
        masks(1) = '01'x
        masks(2) = '03'x
        masks(3) = '07'x
        masks(4) = '0F'x
        masks(5) = '1F'x
        masks(6) = '3F'x
        masks(7) = '7F'x
c-- The bit masks from MSB to LSB
        masks_reverse(1) = '80'x
        masks_reverse(2) = 'C0'x
        masks_reverse(3) = 'E0'x
        masks_reverse(4) = 'F0'x
        masks_reverse(5) = 'F8'x
        masks_reverse(6) = 'FC'x
        masks_reverse(7) = 'FE'x
      endif
c
      if (nbit.eq.0) then
        do i = 1, n
          alignbuf(i) = inbuf(i)
        enddo
      else
        do i = 1, n-1
          alignbuf(i) = ior(ishft(iand(inbuf(i), masks(8-nbit)), nbit), ishft(iand(inbuf(i+1),masks_reverse(nbit)), nbit-8))
        enddo
      endif
      return
      end

      subroutine hist_correct_init(npix, nrscans, nchan, nmir, ngrey)
      implicit none
      integer*4 npix, nrscans, nchan, nmir, ngrey
c
      logical   swdebug
      integer*4 npix_int, nrscans_int, nchan_int, nmir_int, ngrey_int, n_pix, i, j, k, l, i_scan, i_channel, i_mir, i1, i2, i_histmax, histmax
      integer*2 channel_image(n_pix), offset(nchan, nrscans, nmir), zeropoint
      real*4    sumy, sumxy
      character*100 outstring, filenm, filinf
c
      integer*4, allocatable ::  hist(:, :, :, :)
      integer*4, allocatable :: nhist(:, :, :)
c
      save
c
      call getenv('HRPTOUT',outstring)
      npix_int    = npix
      nrscans_int = nrscans
      nchan_int   = nchan
      nmir_int    = nmir
      ngrey_int   = ngrey
      if (.not.allocated( hist)) allocate( hist(ngrey_int, nrscans_int, nchan_int, nmir_int))
      if (.not.allocated(nhist)) allocate(nhist(           nrscans_int, nchan_int, nmir_int))
      do i = 1, nmir_int
        do j = 1, nchan_int
          do k = 1, nrscans_int
            do l = 1, ngrey_int
              hist(l, k, j, i) = 0
            enddo
            nhist(k, j, i) = 0
          enddo
        enddo
      enddo
c
      return
c
      entry hist_create(channel_image, n_pix, i_scan, i_channel, i_mir)
      do i = 1, npix_int
        j = min(max(channel_image(i),0),4095) + 1
        hist(j, i_scan, i_channel, i_mir) = hist(j, i_scan, i_channel, i_mir) + 1
      enddo
      nhist(i_scan, i_channel, i_mir) = nhist(i_scan, i_channel, i_mir) + 1
c
      return
c
      entry hist_analyse(offset, nchan, nrscans, nmir, swdebug)
      do l = 1, nmir
        do i = 1, nchan_int
          do j = 1, nrscans_int
c            if(i.eq.1.and.j.eq.40) goto 2
            if (nhist(j, i, l).gt.0) then
              i_histmax = 2
              histmax   = hist(2, j, i, l)
              do k = 2, ngrey_int - 1
                if (hist(k, j, i, l).gt. histmax) then
                  i_histmax = k
                  histmax   = hist(k, j, i, l)
                endif
              enddo
              do i1 = i_histmax, 2, -1
                if (hist(i1, j, i, l).le.nint(0.3*histmax).and.hist(i1, j, i, l).gt.0) goto 1
              enddo
 1            continue
              i2 = i1 + 10
              if (i.eq.1) i2 = i2 + 10
              sumy  = 0.0
              sumxy = 0.0
              do k = max(i1, 1), min(i2,4096)
                sumy  = sumy  + float(hist(k, j, i, l))
                sumxy = sumxy + float(hist(k, j, i, l)) * float(k)
              enddo
              offset(i, j, l) = nint(sumxy/sumy)
              if (swdebug) write (*,*) l, i, j, nint(sumxy/sumy), i1, i2, histmax, int(0.8*histmax), i_histmax
c              do k = 1, 4096
c                rhist(k, j, i) = float(hist(k, j, i)) / float(nhist(k, j, i))
c              enddo
            endif
c-- For the time being only plot the spectra here
            if (j.ge.2.and.swdebug) then
              write (filinf,'(a,''250m_b'',I2.2,''l'',I2.2)') outstring(1:lnblnk(outstring)), i, j
              write (filenm,'(a,''250m_band_'',I2.2,''_line_'',I2.2,''_M_'',I1)') outstring(1:lnblnk(outstring)), i, j, l
              call hist_plot(hist(1,j,i,l), hist(1,1,i,1), 4096, 1, 4096, 1, filenm, ' ', filinf, i1, i2)
              call system('gmt psconvert '//filenm(1:lnblnk(filenm))//'.ps -A0.2c+s5c -Tg')
              call system('convert '//filenm(1:lnblnk(filenm))//'.png -rotate 90 '//filenm(1:lnblnk(filenm))//'.png')
            endif
c
          enddo
 2      continue
        enddo
        if (swdebug) call system('rm '//outstring(1:lnblnk(outstring))//'*.ps')
      enddo
      do i = 1, nchan_int
        zeropoint = offset(i,1,1)
        do j = 1, nrscans_int
          do l = 1, nmir
            offset(i,j,l) = zeropoint - offset(i,j,l)
            if (swdebug) write (*,*) l, i, j, offset(i,j,l), zeropoint
          enddo
        enddo
      enddo
      return
c
      entry hist_correct()
c
      return
c
      entry hist_correct_flush()
c
      return
      end

      subroutine hist_plot(hist,hist_ref, nx,ny,nxin,nyin,filenm,ctxt, filinf, i1, i2)
      implicit none
      integer*4 nx, ny, hist(nx, ny), nxin, nyin, hist_ref(nx, ny), i1, i2
      character*(*) filenm, ctxt, filinf
c
      integer*4 lunplot, i, j
      real*4 xr(2), yr(2), rhisty(4096), rhistx(4096), yrref(2)
      call get_lun(lunplot)
      call PS_init_colourtable(1,'./resource/color')
      call set_PS_fullpage()
      call pkg_openpl(filenm(1:lnblnk(filenm))//'.ps', lunplot)
      call newpen(6)
      xr(1) = 0.
      xr(2) = float(nxin)
      xr(1) = float(int(i1/100)       * 100)
      xr(2) = float((int(i2/100) + 1) * 100)
      yr(1) = hist(2,1)
      yr(2) = yr(1)
      do i = 2, nxin - 1
        do j = 1, nyin
          if (hist(i,j).lt.yr(1)) yr(1) = hist(i,j)
          if (hist(i,j).gt.yr(2)) yr(2) = hist(i,j)
        enddo
      enddo
c-- Check if the Y scaling for the reference hist is different and take the extremes
      yrref(1) = hist_ref(2,1)
      yrref(2) = yrref(1)
      do i = 2, nxin - 1
        do j = 1, nyin
          if (hist_ref(i,j).lt.yrref(1)) yrref(1) = hist_ref(i,j)
          if (hist_ref(i,j).gt.yrref(2)) yrref(2) = hist_ref(i,j)
        enddo
      enddo
      if (yrref(1).lt.yr(1)) yr(1) = yrref(1) 
      if (yrref(2).gt.yr(2)) yr(2) = yrref(2) 
c-- Fix the vertical scaling
c      yr(2) = 99999.0
      j = len(filenm)
      do while (filenm(j:j).ne.'/')
        j = j - 1
      enddo
      call pkg_frame(11,11,1.,xr,yr,'Bin #','Frequency',filinf(j+1:lnblnk(filinf)))
c      call pkg_frame(11,-6,1.,xr,yr,'Bin #','Frequency',filinf(j+1:lnblnk(filinf)))
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
      xr(1) = float(i1)
      xr(2) = float(i1)
      call pkg_pldatr(-1, xr, yr, -2)
      xr(1) = float(i2)
      xr(2) = float(i2)
      call pkg_pldatr(-1, xr, yr, -2)
      call pkg_pltextbl(ctxt(1:lnblnk(ctxt)),1)
      call pkg_clospl()
      call free_lun(lunplot)
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
