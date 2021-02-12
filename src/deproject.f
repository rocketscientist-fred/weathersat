      implicit none
c
      integer*4 nmax
      parameter (nmax=60000)
      integer*2 rgbbuf(3, nmax), index_correct(nmax), rgbbuf_correct(nmax)
c
      integer*2 npos
      integer*4 i_arg, luntmp, i, i1, i2, i3, nx, ny, lun, lunout, nrecl, ncorrect_pix, ios, nrscanlines
      real*4    alt, fov
      logical   swsharpen, swmersi1, swdebug
      character*1000 argstring, command
c
      alt         = 800.0
      fov         = 110.0
c-- Default AQUA MODIS Like with 40 scanlines per block
      nrscanlines = 40
      swsharpen = .false.
      swmersi1  = .false.
      swdebug   = .false.
      do i_arg = 2, 20
        call getarg(i_arg, argstring)
        if (index(argstring(1:lnblnk(argstring)),'alt=').ne.0) then
          call string_to_r4(argstring(index(argstring,'alt=')+4:lnblnk(argstring)), npos, alt)
          write (*,'(''Average altitude                         : '',F10.3)') alt
        endif
        if (index(argstring(1:lnblnk(argstring)),'fov=').ne.0) then
          call string_to_r4(argstring(index(argstring,'fov=')+4:lnblnk(argstring)), npos, fov)
          write (*,'(''Field of view (deg)                      : '',F10.3)') fov
        endif
        if (index(argstring(1:lnblnk(argstring)),'sharpen').ne.0) swsharpen = .true.
        if (index(argstring(1:lnblnk(argstring)),'bowtie=').ne.0) then
          if (index(argstring(1:lnblnk(argstring)),'mersi1').ne.0) swmersi1 = .true.
        endif
        if (index(argstring(1:lnblnk(argstring)),'nrscanlines=').ne.0) then
          call string_to_i4(argstring(index(argstring,'nrscanlines=')+12:lnblnk(argstring))//' ', npos, nrscanlines)
          write (*,'(''# of scanlines                           : '',I10)') nrscanlines
        endif
        if (index(argstring(1:lnblnk(argstring)),'debug').ne.0) swdebug = .true.
      enddo
      call getarg(1, argstring)
      if (index(argstring,'--help').ne.0) then
        call deproject_help()
        stop
      endif
      write (command,'(''identify '',a,'' >& ./tmp.txt'')') argstring(1:lnblnk(argstring))
      call my_system(command(1:lnblnk(command)))
      call get_lun(luntmp)
      open (unit=luntmp,file='./tmp.txt',access='sequential',form='formatted')
      read (luntmp,'(a)') argstring
      close (unit=luntmp)
      call free_lun(luntmp)
      i  = 1
      i1 = 0
      i2 = 0
      i3 = 0
      do while (i.le.lnblnk(argstring))
        if (argstring(i:i).eq.' '.and.i1.ne.0.and.i2.ne.0.and.i3.eq.0) i3 = i
        if (argstring(i:i).eq.' '.and.i1.ne.0.and.i2.eq.0)             i2 = i
        if (argstring(i:i).eq.' '.and.i1.eq.0)                         i1 = i
        i = i + 1
      enddo
      do i = i2 + 1, i3 - 1
        if (argstring(i:i).eq.'x') i1 = i
      enddo
      call string_to_i4(argstring(i2+1:i1-1)//' ', npos, nx)
      call string_to_i4(argstring(i1+1:i3-1)//' ', npos, ny)
      write (*,'(''nx : '',I5,'' ny : '',I5)') nx, ny
      call getarg(1, argstring)
      i  = lnblnk(argstring)
      i1 = 0
      do while (i1.eq.0.and.i.gt.0)
        if (argstring(i:i).eq.'.'.and.i1.eq.0) i1 = i
        i = i - 1
      enddo
      if (.not.swmersi1) then
c-- Convert the source image to RGB
        write (command,'(''convert '',a,'' -depth 16 rgb:'',a)') argstring(1:lnblnk(argstring)),argstring(1:i1)//'dat'
        write (*,*) command(1:lnblnk(command))
        call my_system(command(1:lnblnk(command)))
      endif
      if (swmersi1) then
        write (command,'(''convert '',a,'' -depth 16 -channel R -separate gray:'',a,a)') argstring(1:lnblnk(argstring)),argstring(1:i1)//'R.'//'dat'
        write (*,*) command(1:lnblnk(command))
        call my_system(command(1:lnblnk(command)))
        write (command,'(''convert '',a,'' -depth 16 -channel G -separate gray:'',a,a)') argstring(1:lnblnk(argstring)),argstring(1:i1)//'G.'//'dat'
        write (*,*) command(1:lnblnk(command))
        call my_system(command(1:lnblnk(command)))
        write (command,'(''convert '',a,'' -depth 16 -channel B -separate gray:'',a,a)') argstring(1:lnblnk(argstring)),argstring(1:i1)//'B.'//'dat'
        write (*,*) command(1:lnblnk(command))
        call my_system(command(1:lnblnk(command)))
        if (swdebug) then
c-- Back to single channel PNG's to check
          write (command,'(''convert -depth 16 -size '',I5.5,''x'',I5.5,1x,''gray:'',a,'' -depth 16 '',a)') nx, ny, argstring(1:i1)//'R.'//'dat',argstring(1:i1)//'R.'//'png'
          write (*,*) command(1:lnblnk(command))
          call my_system(command(1:lnblnk(command)))
          write (command,'(''convert -depth 16 -size '',I5.5,''x'',I5.5,1x,''gray:'',a,'' -depth 16 '',a)') nx, ny, argstring(1:i1)//'G.'//'dat',argstring(1:i1)//'G.'//'png'
          write (*,*) command(1:lnblnk(command))
          call my_system(command(1:lnblnk(command)))
          write (command,'(''convert -depth 16 -size '',I5.5,''x'',I5.5,1x,''gray:'',a,'' -depth 16 '',a)') nx, ny, argstring(1:i1)//'B.'//'dat',argstring(1:i1)//'B.'//'png'
          write (*,*) command(1:lnblnk(command))
          call my_system(command(1:lnblnk(command)))
        endif
c-- Bowtie correct (experimental)
        call bowtie_fix(argstring(1:i1)//'R.'//'dat', nx, nrscanlines)
        write (command,'(''convert -depth 16 -size '',I5.5,''x'',I5.5,1x,''gray:'',a,'' -depth 16 '',a)') nx, ny, argstring(1:i1)//'R.'//'dat.cor',argstring(1:i1)//'R.cor.'//'png'
        write (*,*) command(1:lnblnk(command))
        call my_system(command(1:lnblnk(command)))
        call bowtie_fix(argstring(1:i1)//'G.'//'dat', nx, nrscanlines)
        write (command,'(''convert -depth 16 -size '',I5.5,''x'',I5.5,1x,''gray:'',a,'' -depth 16 '',a)') nx, ny, argstring(1:i1)//'G.'//'dat.cor',argstring(1:i1)//'G.cor.'//'png'
        write (*,*) command(1:lnblnk(command))
        call my_system(command(1:lnblnk(command)))
        call bowtie_fix(argstring(1:i1)//'B.'//'dat', nx, nrscanlines)
        write (command,'(''convert -depth 16 -size '',I5.5,''x'',I5.5,1x,''gray:'',a,'' -depth 16 '',a)') nx, ny, argstring(1:i1)//'B.'//'dat.cor',argstring(1:i1)//'B.cor.'//'png'
        write (*,*) command(1:lnblnk(command))
        call my_system(command(1:lnblnk(command)))
        write (command,'(''convert '',3(a,1x),'' -combine '',a)') argstring(1:i1)//'R.cor.'//'png', argstring(1:i1)//'G.cor.'//'png', argstring(1:i1)//'B.cor.'//'png', argstring(1:i1)//'RGB.cor.'//'png'
        write (*,*) command(1:lnblnk(command))
        call my_system(command(1:lnblnk(command)))
c-- Convert the Bow-Tie corrected source image to RGB to be deprojected
        write (command,'(''convert '',a,'' -depth 16 rgb:'',a)') argstring(1:i1)//'RGB.cor.'//'png',argstring(1:i1)//'dat'
        write (*,*) command(1:lnblnk(command))
        call my_system(command(1:lnblnk(command)))
      endif
c--
      call get_lun(lun)
      call get_lun(lunout)
      nrecl = 6 * nx
      open (unit=lun,file=argstring(1:i1)//'dat', access='direct', form='unformatted', recl=nrecl)
      call correct_init(alt, fov, nx, index_correct, nmax, ncorrect_pix)
c
      open (unit=lunout,file=argstring(1:i1)//'rgb',form='unformatted', access='direct',recl=2*ncorrect_pix*6)
      close (unit=lunout,status='delete', iostat=ios)
      open (unit=lunout,file=argstring(1:i1)//'rgb',form='unformatted', access='direct',recl=2*ncorrect_pix*6)
c-- Force delete RGB file if already present, then re-open new file - bit stupid, but it works
      do i = 1, ny
        call read_buffer(lun, i, rgbbuf, nrecl/2)
        call correct_apply(rgbbuf, 3, nx, rgbbuf_correct, 3, 2 * ncorrect_pix, index_correct, ncorrect_pix, lunout, i)
      enddo
      if (swsharpen) then
        write (command,'(''convert -depth 16 -endian lsb -size '',I5.5,''x'',I5.5,'' rgb:'',a,'' -equalize -depth 16 -sharpen 6x3+0.5+0 '',a)') 2*ncorrect_pix, ny, 
     *            argstring(1:i1)//'rgb', argstring(1:i1-1)//'-deproject-sharpen.png'
      else
        write (command,'(''convert -depth 16 -endian lsb -size '',I5.5,''x'',I5.5,'' rgb:'',a,'' -equalize -depth 16 '',a)') 2*ncorrect_pix, ny, 
     *            argstring(1:i1)//'rgb', argstring(1:i1-1)//'-deproject.png'
      endif
      write (*,*) command(1:lnblnk(command))
      call system(command(1:lnblnk(command)))
      write (command,'(''rm '',a)') argstring(1:i1)//'dat'
      call system(command(1:lnblnk(command)))
      write (command,'(''rm '',a)') argstring(1:i1)//'rgb'
      call system(command(1:lnblnk(command)))
c--
      stop
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

      subroutine read_buffer(lun,irec,array,n)
      implicit none
      integer*4 n, irec, lun
      integer*2 array(n)
c
      read (lun,rec=irec) array
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

      subroutine string_to_r4(string,npos,value)
      logical*1 swexp
      integer*2 npos,ilen,i_mantisse,i_main
      real*4    value
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
 2        continue
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

      subroutine bowtie_fix(filenm, nmaxpix, nlines)
      implicit none
      integer*4 nmaxpix, nlines
      character*(*) filenm
c
      integer*4 nrscanlines
      real*4, allocatable :: displace(:, :)
c
      if (.not.allocated(displace)) allocate(displace(nmaxpix, nlines))
      call bowtie_displacement(displace, nmaxpix, nlines)
      call bowtie_lines(filenm, nmaxpix, nrscanlines)
      call bowtie_correct(filenm, displace, nmaxpix, nlines, nrscanlines)
      if (allocated(displace)) deallocate(displace)
      return
      end
      
      subroutine bowtie_correct(filenm, displace, nmaxpix, nlines, nrscanlines)
      implicit none
      integer*4 nmaxpix, nlines, nrscanlines
      real*4    displace(nmaxpix, nlines)
      character*(*) filenm
c
      integer*4 lunin, lunout, i, j, k
      real*4    w1, w2, d
c
      integer*2, allocatable :: inbuf(:)
      integer*4, allocatable :: rawimg(: , :)
      integer*4, allocatable :: corimg(: , :)
      real*4,    allocatable :: pixcount(: , :)
c
      call get_lun(lunin)
      if (.not.allocated(inbuf)) allocate(inbuf(nmaxpix))
      open (unit=lunin,file=filenm, access='direct',form='unformatted', recl=nmaxpix*2)
      if (.not.allocated(pixcount)) allocate(pixcount(nmaxpix, nrscanlines))
      if (.not.allocated(rawimg))   allocate(rawimg(nmaxpix, nrscanlines))
      if (.not.allocated(corimg))   allocate(corimg(nmaxpix, nrscanlines))
      do i = 1, nmaxpix
        do j = 1, nrscanlines
          pixcount(i,j) = 0.0
          rawimg(i,j)   = 0
          corimg(i,j)   = 0
        enddo
      enddo
      do i = 1, nrscanlines
        read (lunin,rec=i) inbuf
        do j = 1, nmaxpix
          if (inbuf(j).lt.0) then
            rawimg(j,i) = inbuf(j) + 65536
          else
            rawimg(j,i) = inbuf(j)
          endif
        enddo
      enddo
      close (unit=lunin)
      call free_lun(lunin)
      do i = 1, nmaxpix
        do j = 1, nrscanlines
          d = displace(i, nlines + 1 - (mod(j-1,nlines) + 1))
          w1 = abs(d - int(d))
          w2 = 1.0 - w1
          k = j + int(d)
          if (d.ge.0.0) then
            k = k + 1
            if (1.le.k.and.k.le.nrscanlines) then
              corimg(i,k) = corimg(i,k) + w1 * rawimg(i,j)
              pixcount(i,k) = pixcount(i,k) + w1
            endif
            k = k - 1
            if (1.le.k.and.k.le.nrscanlines) then
              corimg(i,k) = corimg(i,k) + w2 * rawimg(i,j)
              pixcount(i,k) = pixcount(i,k) + w2
            endif
          else
            k = k - 1
            if (1.le.k.and.k.le.nrscanlines) then
              corimg(i,k) = corimg(i,k) + w1 * rawimg(i,j)
              pixcount(i,k) = pixcount(i,k) + w1
            endif
            k = k + 1
            if (1.le.k.and.k.le.nrscanlines) then
              corimg(i,k) = corimg(i,k) + w2 * rawimg(i,j)
              pixcount(i,k) = pixcount(i,k) + w2
            endif
          endif
        enddo
      enddo
      do i = 1, nmaxpix
        do j = 1, nrscanlines
          if (pixcount(i,j).ne.0.0) corimg(i,j) = nint(float(corimg(i,j)) / pixcount(i,j))
c          if (pixcount(i,j).gt.0.2) corimg(i,j) = nint(float(corimg(i,j)) / pixcount(i,j))
        enddo
      enddo
      do i = 1, nmaxpix
        do j = 1, nrscanlines
          if (pixcount(i,j).eq.0.0.and.2.le.j.and.j.le.nrscanlines-1) corimg(i,j) = (corimg(i,j-1)+corimg(i,j+1)) / 2.0
c          if (pixcount(i,j).le.0.2.and.2.le.j.and.j.le.nrscanlines-1) corimg(i,j) = (corimg(i,j-1)+corimg(i,j+1)) / 2.0
        enddo
      enddo
      call get_lun(lunout)
      open (unit=lunout,file=filenm//'.cor', access='direct',form='unformatted', recl=nmaxpix*2)
      do i = 1, nrscanlines
        do j = 1, nmaxpix
          if (0.le.corimg(j,i).and.corimg(i,j).le.32767) then
            inbuf(j) = corimg(j,i)
          else
            inbuf(j) = corimg(j,i) - 65536
          endif
        enddo
        write (lunout,rec=i) inbuf
      enddo
c
      close (unit=lunout)
      call free_lun(lunout)
      if (allocated(inbuf)) deallocate(inbuf)
      if (allocated(pixcount)) deallocate(pixcount)
      if (allocated(rawimg))   deallocate(rawimg)
      if (allocated(corimg))   deallocate(corimg)
      return
      end

      subroutine bowtie_lines(filenm, nmaxpix, nrscanlines)
      implicit none
      integer*4 nmaxpix, nrscanlines
      character*(*) filenm
c
      integer*4 lun, ios
      integer*2, allocatable :: inbuf(:)
c
      call get_lun(lun)
      if (.not.allocated(inbuf)) allocate(inbuf(nmaxpix))
      open (unit=lun,file=filenm, access='direct',form='unformatted', recl=nmaxpix*2)
      nrscanlines = 0
      do while (.true.)
        read (lun,rec=nrscanlines+1, iostat=ios) inbuf
        if (ios.ne.0) goto 1
        nrscanlines = nrscanlines + 1
      enddo
 1    continue
c
      close (unit=lun)
      call free_lun(lun)
      if (allocated(inbuf)) deallocate(inbuf)
      return
      end

      subroutine bowtie_displacement(displace, nmaxpix, nlines)
      implicit none
      integer*4 nmaxpix, nlines
      real*4    displace(nmaxpix, nlines)
c
      integer*4 i, j, nsat_pix, i_arg, init
      integer*2 npos
      real*8 rpos, bowgamma, theta, theta1, theta2, torad, todeg, displace_y, theta_line_1, theta_line_2
      logical swdebug
      character*250 argstring
      data init/0/
c
      save init
c
      torad           = 2.0D0 * dasin(1.0D0) / 180.0D0
      todeg           = 180.0D0 / (2.0D0 * dasin(1.0D0))
      swdebug         = .false.
c
      bowgamma = 1.0D0
      nsat_pix = nmaxpix
      theta1   = 0.00D0
      theta2   = 0.80D0
c-- If need be override these parameters on the command line to allow finding optimum Bow-Tie parameters
      do i_arg = 2, 20
        call getarg(i_arg, argstring)
        if (index(argstring(1:lnblnk(argstring)),'bowgamma=').ne.0) then
          call string_to_r8(argstring(index(argstring,'bowgamma=')+9:lnblnk(argstring)), npos, bowgamma)
          write (*,'(''Bow-gamma                                : '',F10.3)') bowgamma
        endif
        if (index(argstring(1:lnblnk(argstring)),'theta1=').ne.0) then
          call string_to_r8(argstring(index(argstring,'theta1=')+7:lnblnk(argstring)), npos, theta1)
          write (*,'(''Theta1                                   : '',F10.3)') theta1
        endif
        if (index(argstring(1:lnblnk(argstring)),'theta2=').ne.0) then
          call string_to_r8(argstring(index(argstring,'theta2=')+7:lnblnk(argstring)), npos, theta2)
          write (*,'(''Theta2                                   : '',F10.3)') theta2
        endif
        if (index(argstring(1:lnblnk(argstring)),'debug').ne.0) swdebug = .true.
      enddo
c
      do j = 0, nlines - 1
        theta_line_1 = (theta1/2.0D0) - dble(j) * (theta1 / dble(nlines-1))
        theta_line_2 = (theta2/2.0D0) - dble(j) * (theta2 / dble(nlines-1))
c        write (*,*) j, theta_line_1, theta_line_2
        displace_y = 0.0
        do i = (nmaxpix / 2) - 1, 1, -1
          rpos = abs((dble(2*i)/dble(nsat_pix)) - 1.0D0)
          theta  = theta_line_1 + (theta_line_2 - theta_line_1) * (rpos ** (1.0D0/bowgamma))
          displace_y = displace_y + sin(theta * torad)
c          displace_y = sin(theta * torad) * float(abs(((nmaxpix/2)-i)))
          displace(i, j + 1) = displace_y
c          write (*,'(I5,2F10.3)') i, theta, displace
        enddo
        displace_y = 0.0
        do i = (nmaxpix / 2), nmaxpix
          rpos = abs((dble(2*i)/dble(nsat_pix)) - 1.0D0)
          theta  = theta_line_1 + (theta_line_2 - theta_line_1) * (rpos ** (1.0D0/bowgamma))
          displace_y = displace_y + sin(theta * torad)
c          displace_y = sin(theta * torad) * float(abs(((nmaxpix/2)-i)))
          displace(i, j + 1) = displace_y
c          write (*,'(I5,2F10.3)') i, theta, displace
        enddo
      enddo
      if (swdebug.and.init.eq.0) then
        init = 1
        call displace_plot(displace, nmaxpix, nlines)
      endif
c      do j = 1, nlines
c        write (*,*) j, displace(1,j), displace(nmaxpix/2,j), displace(nmaxpix,j)
c      enddo
      return
      end
      
      subroutine displace_plot(displace, nmaxpix, nlines)
      implicit none
      integer*4 nmaxpix, nlines
      real*4    displace(nmaxpix, nlines)
c      
      integer*4 lunplot, i, j, k
      real*4 xr(2), yr(2), yoff
      real*4,    allocatable :: x(:)
      real*4,    allocatable :: y(:)

c
      if (.not.allocated(x)) allocate(x(nmaxpix))
      if (.not.allocated(y)) allocate(y(nmaxpix))
c
      call get_lun(lunplot)
      call PS_init_colourtable(1,'./resource/color')
      call set_PS_fullpage()
      call pkg_openpl('BowTie_displace.ps', lunplot)
      call newpen(6)
      xr(1) = 0.
      xr(2) = float(nmaxpix)
      yr(1) = displace(1,1)
      yr(2) = displace(1,1)
      do i = 1, nlines
        do j = 1, nmaxpix
          if (displace(j,i).gt.yr(2)) yr(2) = displace(j,i)
          if (displace(j,i).lt.yr(1)) yr(1) = displace(j,i)
        enddo
      enddo
      if (mod(nlines,2).eq.0) then
        yr(1) = yr(1) - float(nlines/2) + 0.5
        yr(2) = yr(2) + float(nlines/2) - 0.5
      else
        yr(1) = yr(1) - float(nlines/2)
        yr(2) = yr(2) + float(nlines/2)
      endif
c-- Fix the vertical scaling
      call pkg_frame(11,11,1.,xr,yr,'Pixel','Displacement','BowTie Correction')
      call pkg_raster(11,11,xr,yr,-2)
      do k = 1, nlines
        i = nlines + 1 - k
        if (mod(nlines,2).eq.0) then
          yoff = float(i) - float(nlines/2) - 0.5
        else
          yoff = float(i) - float(nlines/2) - 1.0
        endif
        do j = 1, nmaxpix
          x(j) = float(j)
          y(j) = displace(j,k) + yoff
        enddo
        call pkg_pldatr(-1, x, y, nmaxpix)
      enddo
      call pkg_clospl()
      call free_lun(lunplot)
c
      if (allocated(x)) deallocate(x)
      if (allocated(y)) deallocate(y)
c
      return
      end
      
      subroutine deproject_help()
c
      write (*,'(''Run the program as follows: '')')
      write (*,'('' '')')
      write (*,'(''./deproject.exe /mnt/y/MERSI1-RGB-221.png alt=800.0 fov=110.0 '')')
      write (*,'('' '')')
      write (*,'('' First argument has to be the filename to be deprojected '')')
      write (*,'('' If the arguments alt= and/or fov= are not specified the defaults are alt=800.0 and fov=110.0 '')')
      write (*,'('' --help provides this help text'')')
      write (*,'('' The output ends in the same directory as the source file with -deproject added to the end of the filename (before the .png extension)'')')
      write (*,'('' '')')
      write (*,'('' Arguments:   '')')
      write (*,'('' '')')
      write (*,'(''     alt=810.0       define the average altitude of the satellite'')')
      write (*,'(''     fov=110.0       field-of-view of a scanline (deg)'')')
      write (*,'(''     bowtie=mersi1   Applies Bow-Tie correction. Currently only mersi1 supported, rest to be entered through expert parameters (below)'')')
      write (*,'(''     sharpen         create a sharpened (through unsharp masking) image'')')
      write (*,'('' '')')
      write (*,'('' Expert arguments (bowtie correction parameters:'')')
      write (*,'(''     theta1=0.0'')')
      write (*,'(''     theta2=0.8'')')
      write (*,'(''     bowgamma=1.0'')')
      write (*,'('' '')')
      write (*,'('' These are the bowtie parametrisation parameters (here for MERSI-1 (and MERSI-2 ??)'')')
      write (*,'('' For AUQA (in readbin_modis, these are bowgamma=1.6, theta1=0.0, theta2=1.05'')')
      write (*,'('' '')')
      write (*,'('' '')')
      return
      end
