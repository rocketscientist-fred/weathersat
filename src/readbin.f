      implicit none
      integer*4 nbuf, ibuf, nbuffy
      parameter (nbuf=12960,nbuffy=26050)
      integer*1 buffer(nbuf), zerobuf(nbuf), buffy(nbuffy), zerobuffy(nbuffy)
c
      logical swvalidrec, swmetop, swfy, swdebug, swmn2
      integer*1 inbuf(1024), nextbuf(1024), b1, b2, b3, b4, alignbuf(1024)
      integer*4 lunfile, irec, ios, k, iptr, i4_scid, i4_vcid, i4_ccsds_ptr, i4_apid, i4_packet_count, i4_packet_length, i4_seq_flag, i4_packet_save
      integer*4 VCDU_data_length, VCDU_data_remaining, VCDU_full_records, VCDU_count_save, VCDU_full_left, VCDU_data_left, VCDU_next_header
      integer*4 nbyte, nbit, i, nbit_save, i_score, iline, nbit_last
      character*250 filedata, argstring
      character*500 carg
c
      integer*1 i1_vcdu_count(4)
      integer*4 i4_vcdu_count
c
      equivalence (i4_vcdu_count,i1_vcdu_count(1))
c
      do i = 1, nbuf
        zerobuf(i) = 0
        buffer(i)  = 0
      enddo
      do i = 1, nbuffy
        zerobuffy(i) = 0
        buffy(i)     = 0
      enddo
c
      swmetop = .false.
      swfy    = .false.
      swmn2   = .false.
      swdebug = .false.
      call getarg(1,argstring)
      filedata = argstring(1:lnblnk(argstring))
      if (index(filedata(1:lnblnk(filedata)),'Metop').ne.0) swmetop = .true.
      if (index(filedata(1:lnblnk(filedata)),'FY3').ne.0)   swfy    = .true.
      if (index(filedata(1:lnblnk(filedata)),'MN2').ne.0)   swmn2   = .true.
      call getarg(3,argstring)
      if (index(argstring,'debug').ne.0) swdebug = .true.
c
      if (swmn2) then 
        call scan_mn2_signature(filedata(1:lnblnk(filedata)), swdebug)
        goto 2
      endif
c
      VCDU_data_length = 886 - 2
      if (swmetop) VCDU_data_length = VCDU_data_length - 2
      if (swfy) VCDU_data_length = VCDU_data_length - 2
      VCDU_data_left   =  0
      VCDU_full_left   = -1
      VCDU_next_header = -1
      VCDU_count_save  = -1
      ibuf             =  1
c
      call get_lun(lunfile)
      open (unit=lunfile,file=filedata(1:lnblnk(filedata)),form='unformatted',access='direct',recl=1024)
      irec = 1
      b1   = '1A'x
      b2   = 'CF'x
      b3   = 'FC'x
      b4   = '1D'x
      do while (.true.)
        ios = 0
        read (lunfile,rec=irec,iostat=ios) inbuf
        if (ios.ne.0) goto 1
        ios = 0
        read (lunfile,rec=irec+1,iostat=ios) nextbuf
        swvalidrec = .false.
        i_score = 0
        if (inbuf(1).eq.b1) i_score = i_score + 1
        if (inbuf(2).eq.b2) i_score = i_score + 1
        if (inbuf(3).eq.b3) i_score = i_score + 1
        if (inbuf(4).eq.b4) i_score = i_score + 1
        if (i_score.ge.0) swvalidrec = .true.
        if (.not.swvalidrec) then
          write (*,'(I10,4(2x,z2.2))') irec, (inbuf(k), k=1,4)
        else
          i4_scid = iand(inbuf(5),'3F'x) * 4 + ishft(iand(inbuf(6),'B0'x),-6)
          i4_vcid = iand(inbuf(6),'3F'x)
          call move_bytes_reverse(inbuf(7),i1_vcdu_count(1),3)
          if (swdebug) write (*,*) i4_scid, i4_vcid, i4_vcdu_count
c-- As the entire analysis depens on the correct vcdu counter, do a simple protection against getting an erroneous one !
          if (abs(i4_vcdu_count-VCDU_count_save).ge.25000.and.VCDU_count_save.ge.0) goto 998
c
          iptr = 11
          if (swmetop.or.swfy) iptr = iptr + 2
          i4_ccsds_ptr = iand(inbuf(iptr),'07'x) * 256
          if (inbuf(iptr+1).lt.0) then
            i4_ccsds_ptr = i4_ccsds_ptr + inbuf(iptr+1) + 256
          else
            i4_ccsds_ptr = i4_ccsds_ptr + inbuf(iptr+1)
          endif
          if (swdebug.and.swmetop) write (*,'(''CCSDS Header Pointer : '',I10)') i4_ccsds_ptr
c-- One of the sanity checks to catch bad data
          if (swmetop.and.(i4_scid.lt.10.or.i4_scid.gt.14)) goto 998
c-- One of the sanity checks to catch bad data
          if (swmetop.and.i4_ccsds_ptr.gt.900.and.i4_ccsds_ptr.ne.2047) goto 998
          if (swmetop) then
            if (i4_ccsds_ptr.ne.2047) then
c-- Flush any outstanding data
              if (VCDU_full_left.eq.0) then
                call fill_buffer(buffer,nbuf,ibuf,inbuf(iptr+2),VCDU_data_left, iline)
              else
                call fill_buffer(buffer,nbuf,ibuf,zerobuf,VCDU_data_left, iline)
              endif
              VCDU_full_left  = 0
              VCDU_count_save = i4_vcdu_count
              VCDU_full_records = 0
              iptr = iptr + 2 + i4_ccsds_ptr
              i4_apid = iand(inbuf(iptr),'07'x) * 256
              if (inbuf(iptr+1).lt.0) then
                i4_apid = i4_apid + inbuf(iptr+1) + 256
              else
                i4_apid = i4_apid + inbuf(iptr+1)
              endif
              if (swdebug) write (*,*) i4_apid
              iptr = iptr + 2
              i4_seq_flag = ishft(iand(inbuf(iptr),'B0'x),-6)
              i4_packet_save  = i4_packet_count
              i4_packet_count = iand(inbuf(iptr),'3F'x) * 256
              if (inbuf(iptr+1).lt.0) then
                i4_packet_count = i4_packet_count + inbuf(iptr+1) + 256
              else
                i4_packet_count = i4_packet_count + inbuf(iptr+1)
              endif
c-- One of the sanity checks to catch bad data
              if ((abs(i4_packet_count-i4_packet_save).ge.2500).or.(i4_packet_count.le.i4_packet_save)) goto 998
              iptr = iptr + 2
              if (inbuf(iptr).lt.0) then
                i4_packet_length = (inbuf(iptr) + 256) * 256
              else
                i4_packet_length =  inbuf(iptr) * 256
              endif
              if (inbuf(iptr+1).lt.0) then
                i4_packet_length = i4_packet_length + inbuf(iptr+1) + 256
              else
                i4_packet_length = i4_packet_length + inbuf(iptr+1)
              endif
              i4_packet_length = i4_packet_length + 1
              if (swdebug) write (*,*) i4_seq_flag, i4_packet_count, i4_packet_length
c-- One of the sanity checks to catch bad data
              if (i4_packet_length.ne.12960) i4_packet_length = 12960
              VCDU_data_remaining = VCDU_data_length - i4_ccsds_ptr - 6
              VCDU_full_records   = int(float(i4_packet_length - VCDU_data_remaining)/float(VCDU_data_length))
              VCDU_full_left      = VCDU_full_records
              if (swdebug) write (*,*) 'Bytes after this header : ',VCDU_data_remaining,' Full records following : ', VCDU_full_records
              call fill_buffer(buffer,nbuf,ibuf,inbuf(iptr+2),VCDU_data_remaining, iline)
              VCDU_data_left = i4_packet_length - VCDU_data_remaining
            else
              if (VCDU_full_left - (i4_vcdu_count- VCDU_count_save).ge.0) then
                VCDU_full_left = VCDU_full_left - (i4_vcdu_count- VCDU_count_save)
                if (swdebug) write (*,*) 'VCDU_full_left : ',VCDU_full_left
                if (i4_vcdu_count- VCDU_count_save.eq.1) then
                  VCDU_count_save = i4_vcdu_count
                  VCDU_data_left  = VCDU_data_left - VCDU_data_length
                  if (swdebug) write (*,*) 'VCDU_data_left : ',VCDU_data_left
                  call fill_buffer(buffer,nbuf,ibuf,inbuf(iptr+2),VCDU_data_length, iline)
                else
                  if (VCDU_full_left.gt.0) then
c-- One or more records were missing, first write VCDU_data_length zeroes for every missing record !
                    do i = VCDU_count_save + 1, i4_vcdu_count - 1
                      VCDU_data_left  = VCDU_data_left - VCDU_data_length
                      if (swdebug) write (*,*) 'VCDU_data_left : ',VCDU_data_left, 'Written zero'
                      call fill_buffer(buffer,nbuf,ibuf,zerobuf,VCDU_data_length, iline)
                    enddo
                    VCDU_count_save = i4_vcdu_count
                    VCDU_data_left  = VCDU_data_left - VCDU_data_length
                    if (swdebug) write (*,*) 'VCDU_data_left : ',VCDU_data_left
                    call fill_buffer(buffer,nbuf,ibuf,inbuf(iptr+2),VCDU_data_length, iline)
                  else
c-- A few records were missing and the header of the next block was in there write the remaining zero bytes in the current 
c--   record and continue the search for a new header
c-- I don't think I should ever get here ?
                  endif
                endif
              else
                if (VCDU_data_left.gt.0) then
                if (swdebug) write (*,*) 'Gap detected : flushing',VCDU_data_left,' bytes'
c-- Enter the 'flush' routine here
                  call fill_buffer(buffer,nbuf,ibuf,zerobuf,VCDU_data_left, iline)
                VCDU_data_left = 0
                endif
              endif
            endif
          endif
          if (swfy) then
            iptr = iptr + 2
            call find_fy_mask(inbuf(iptr), nextbuf(iptr), VCDU_data_length, nbyte, nbit)
            if (nbyte.ge.0.and.swdebug) write (*,*) nbyte, nbit, VCDU_full_left, VCDU_full_records, i4_vcdu_count, VCDU_count_save, VCDU_next_header
            if (nbit.ge.0.and.nbit.ne.nbit_save.and.nbit_save.ge.0) nbit_last = nbit_save
            if (nbit.ge.0) nbit_save = nbit
            if (nbyte.ge.0.and.i4_vcdu_count.lt.VCDU_next_header) nbyte = -1
            if (nbyte.eq.-1) then
              i4_ccsds_ptr = 2047
            else
              i4_ccsds_ptr = nbyte
            endif
            if (swdebug) write (*,'(''Header Pointer : '',I10)') i4_ccsds_ptr
            if (i4_ccsds_ptr.ne.2047) then
c-- Flush any outstanding data
              if (VCDU_full_left.eq.0) then
c-- Catch cases I have seen of a corrupted header pointer - but be careful with the data_left = 882 case !
c-- Could also be a case of a wrongly aligned header - so use nbit_last !
                if (VCDU_data_left.lt.882.and.(i4_ccsds_ptr.ne.VCDU_data_left+1).and.VCDU_data_left.gt.0) then
                  i4_ccsds_ptr = VCDU_data_left + 1
                  call fill_buffy(buffy,nbuffy,ibuf,inbuf(iptr),VCDU_data_left,nbit_last, iline)
                else
                  call fill_buffy(buffy,nbuffy,ibuf,inbuf(iptr),VCDU_data_left,nbit_save, iline)
                endif
              else
                do i = 1, VCDU_full_left
                  call fill_buffy(buffy,nbuffy,ibuf,zerobuf,VCDU_data_length,nbit_save, iline)
                  VCDU_data_left  = VCDU_data_left - VCDU_data_length
                enddo
                if (i4_vcdu_count-VCDU_count_save.lt.28) then
                  call fill_buffy(buffy,nbuffy,ibuf,inbuf(iptr),VCDU_data_left,nbit_save, iline)
                else
                  call fill_buffy(buffy,nbuffy,ibuf,zerobuf,VCDU_data_left,nbit_save, iline)
                endif
              endif
              VCDU_full_left  = 0
              VCDU_count_save = i4_vcdu_count
c-- Packet length including the 3 sync bytes ?? (Kunstmaan article says 26047
              i4_packet_length = 26050
              VCDU_data_remaining = VCDU_data_length - i4_ccsds_ptr + 1
              VCDU_full_records   = int(float(i4_packet_length - VCDU_data_remaining)/float(VCDU_data_length))
              VCDU_full_left      = VCDU_full_records
              VCDU_next_header    = i4_vcdu_count + VCDU_full_records + 1
              if (swdebug) write (*,*) 'Bytes after this header : ',VCDU_data_remaining,' Full records following : ', VCDU_full_records, ' Next header in rec # : ', VCDU_next_header
              call fill_buffy(buffy,nbuffy,ibuf,inbuf(iptr+i4_ccsds_ptr-1),VCDU_data_remaining, nbit_save, iline)
              VCDU_data_left = i4_packet_length - VCDU_data_remaining
            else
              if (VCDU_full_left - (i4_vcdu_count- VCDU_count_save).ge.0) then
                VCDU_full_left = VCDU_full_left - (i4_vcdu_count- VCDU_count_save)
                if (swdebug) write (*,*) 'VCDU_full_left : ',VCDU_full_left
                if (i4_vcdu_count- VCDU_count_save.eq.1) then
                  VCDU_count_save = i4_vcdu_count
                  VCDU_data_left  = VCDU_data_left - VCDU_data_length
                  if (swdebug) write (*,*) 'VCDU_data_left : ',VCDU_data_left
                  call fill_buffy(buffy,nbuffy,ibuf,inbuf(iptr),VCDU_data_length, nbit_save, iline)
                else
                  if (VCDU_full_left.gt.0) then
c-- One or more records were missing, first write VCDU_data_length zeroes for every missing record !
                    do i = VCDU_count_save + 1, i4_vcdu_count - 1
                      VCDU_data_left  = VCDU_data_left - VCDU_data_length
                      if (swdebug) write (*,*) 'VCDU_data_left : ',VCDU_data_left, 'Written zero'
                      call fill_buffy(buffy,nbuffy,ibuf,zerobuf,VCDU_data_length, nbit_save, iline)
                    enddo
                    VCDU_count_save = i4_vcdu_count
                    VCDU_data_left  = VCDU_data_left - VCDU_data_length
                    if (swdebug) write (*,*) 'VCDU_data_left : ',VCDU_data_left
                    call fill_buffy(buffy,nbuffy,ibuf,inbuf(iptr),VCDU_data_length, nbit_save, iline)
                  else
c-- A few records were missing and the header of the next block was in there write the remaining zero bytes in the current 
c--   record and continue the search for a new header
c-- I don't think I should ever get here ?
                  endif
                endif
              else
                if (VCDU_data_left.gt.0) then
                  if (swdebug) write (*,*) 'Gap detected : flushing',VCDU_data_left,' bytes'
c-- Enter the 'flush' routine here
                  call fill_buffy(buffy,nbuffy,ibuf,zerobuf,VCDU_data_left, nbit_save, iline)
                  VCDU_data_left = 0
                endif
              endif
            endif
          endif
 999      continue
        endif
 998    continue
        irec = irec + 1
      enddo
 1    irec = irec - 1
      write (*,'('' Lines in image : '',I10)') iline
 2    continue
      stop
      end

      subroutine fill_buffer(buf,n,iptr,inbuf,nin, iline)
      implicit none
      integer*4 n, iptr, nin, iline
      integer*1 buf(n), inbuf(nin)
c
      integer*4 i, irec, init, lun
c
      real*8    timestamp
      integer*4 ios
      character*250 filename
c
      save
c
      if (init.eq.0) then
        init = 1
        irec = 0
        call getarg(2,filename)
        call get_lun(lun)
        open (unit=lun, file=filename(1:lnblnk(filename)), access='direct', form='unformatted', recl=n)
      endif
      do i = 1, nin
        if (iptr+i-1.le.n) buf(iptr+i-1) = inbuf(i)
      enddo
      iptr = iptr + nin
      if (iptr.gt.n) then
        irec = irec + 1
        write (lun,rec=irec) buf
c        write (*,*) 'Flushing buffer : ', irec
        iptr = 1
        do i = 1,n
          buf(i) = 0
        enddo
c
      endif
      iline = irec
      return
      end

      subroutine fill_buffy(buf,n,iptr,inbuf,nin, nbit, iline)
      implicit none
      integer*4 n, iptr, nin, nbit, iline
      integer*1 buf(n), inbuf(nin)
c
      integer*4 i, irec, init, lun, k
c
      real*8    timestamp
      integer*4 ios
      integer*1 alignbuf(26050)
      character*250 filename
c
      save
c
      if (init.eq.0) then
        init = 1
        irec = 0
        call getarg(2,filename)
        call get_lun(lun)
        open (unit=lun, file=filename(1:lnblnk(filename)), access='direct', form='unformatted', recl=n)
      endif
      do i = 1, nin
        if (iptr+i-1.le.n) buf(iptr+i-1) = inbuf(i)
      enddo
      iptr = iptr + nin
      if (iptr.gt.n) then
        call bit_align_buffer(buf, n, alignbuf, 26050, nbit) 
        nbit = -1
        irec = irec + 1
        write (lun,rec=irec) alignbuf
        iptr = 1
        do i = 1,n
          buf(i) = 0
        enddo
c
      endif
      iline = irec
      return
      end

c-- This subroutine searches for the FY3 line sync marker which is NOT necessarily aligned at byte level .... sigh
c-- So, search at bit level !
       subroutine find_fy_mask(buffer, nextbuf, n, nbyte, nbit)
      implicit none
      integer*4 n, nbyte, nbit
      integer*1 buffer(n), nextbuf(n)
c
      integer*1 mask(3), bytemask(4), masks(7), masks_reverse(7), bcheck
      integer*4 i, j, init
      logical   swmatch
c
      data init/0/
      save
c
      if (init.eq.0) then
        init    = 1
c-- The FY3 sync marker (from "De kunstmaan", 2017)
        mask(1) = '84'x
        mask(2) = '5B'x
        mask(3) = 'F5'x
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
          else
            bytemask(1) =                                         ishft(mask(1),-i)
            bytemask(2) = ior  (ishft(iand(mask(1),masks(i)),8-i),ishft(mask(2),-i))
            bytemask(3) = ior  (ishft(iand(mask(2),masks(i)),8-i),ishft(mask(3),-i))
            bytemask(4) =       ishft(iand(mask(3),masks(i)),8-i)
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
            if (j+3.le.n) then
              bcheck = buffer(j+3)
            else
              bcheck = nextbuf(j + 3 - n)
            endif
            if (iand(bcheck,masks_reverse(i)).ne.bytemask(4)) swmatch = .false.
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

      subroutine scan_mn2_signature(file, swdebug)
      implicit none
      character*(*) file
      logical swdebug
c
      logical   swsync
      integer*1 inbuf(1024), sigmn2(8), buffer(11850)
      integer*4 lun, irec, ios, i_mn2, mn2tmrange(2,4), i, j, k, i_ptr, lunout, iline
      character*250 fileout
c
      data sigmn2/'02'x,'18'x,'A7'x,'A3'x,'92'x,'DD'x,'9A'x,'BF'x/
c-- Where are the 948 bytes of MSU-MR TM in the 1024 byte block
      data mn2tmrange/  23, 260, 279, 516, 535, 772, 791, 1024/
c
      do k = 1, 11850
        buffer(k) = 0
      enddo
      call get_lun(lun)
      open (unit=lun,file=file(1:lnblnk(file)),form='unformatted',access='direct',recl=1024)
      call getarg(2,fileout)
      call get_lun(lunout)
      open (unit=lunout, file=fileout(1:lnblnk(fileout)), access='direct', form='unformatted', recl=11850)
      iline  = 1
      i_mn2  = 1
      irec   = 1
      swsync = .false.
      do while (.true.)
        ios = 0
        read (lun,rec=irec,iostat=ios) inbuf
        if (ios.ne.0) goto 1
        do i = 1,4
          do j = mn2tmrange(1,i), mn2tmrange(2,i)
            if (inbuf(j).eq.sigmn2(i_mn2)) then
              i_mn2 = i_mn2 + 1
              if (i_mn2.eq.9) then
                if (swdebug) write (*,*) 'Signature  ', irec, j
c-- New signature encountered before flushing buffer sets swsync to false. Meaning a few records were lost.
                if (swsync) then
                  write (lunout, rec=iline) buffer
                  do k = 1, 11850
                    buffer(k) = 0
                  enddo
                  iline = iline + 1
                  if (swdebug) write (*,*) 'Aborted flush', i_ptr, irec, j
                  swsync = .false.
                endif
                swsync = .true.
                i_mn2  = 1
                i_ptr  = 8
              endif
            else
              i_mn2  = 1
            endif
            if (swsync) then
              buffer(i_ptr) = inbuf(j)
              i_ptr = i_ptr + 1
c-- The buffer has been filled by all the data - flush and reset swsync
              if (i_ptr.gt.11850) then
                write (lunout, rec=iline) buffer
                do k = 1, 11850
                  buffer(k) = 0
                enddo
                iline = iline + 1
                if (swdebug) write (*,*) 'Regular flush', i_ptr, irec, j
                swsync = .false.
              endif
            else
            endif
          enddo
        enddo
        irec = irec + 1
      enddo
 1    iline = iline - 1
      close (unit=lun)
      call free_lun(lun)
      close (unit=lunout)
      call free_lun(lunout)
      write (*,'('' Lines in image : '',I10)') iline
c
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

