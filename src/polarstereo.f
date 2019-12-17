      implicit none
c--
      integer*4 nx, ny, nxs, nys
      parameter (nx=21600, ny=10800, nxs=5000, nys=5000)
c--
      integer*1 mapline (3,nx), map(3,nx,ny)
      integer*1 mapsline (3,nxs), maps(3,nxs,nys)
      integer*4 lun, i, j, ix, iy, itype, value
      real*4    y_delta, x_delta, longs, lats, xs, ys, xr, yr, phi
      real*8    long, lat, torad, todeg, pi, pid4
      character*7 cmap
      character*9 csize
      character*100 argstring
c
      save
c
      call getarg(1,argstring)
      read (argstring(1:lnblnk(argstring)),*) phi
c
      cmap = 'polster'
      write (csize,'(I4.4,''x'',I4.4)') nxs, nys
c
      torad   = 2.0D0 * dasin(1.0D0) / 180.0
      todeg   = 180.0D0 / 2.0D0 / dasin(1.0D0)
      pi      = 2.0D0 * dasin(1.0D0)
      pid4    = pi / 4.0D0
      phi     = phi * torad
c
      call get_lun(lun)
      open (unit=lun, file='./resource/worldmapuhr.rgb', access='direct', form='unformatted', recl=3*nx)
      do i = 1, ny
        read (lun,rec=i) mapline
        do j = 1, nx
          map(1,j,i) = mapline (1,j)
          map(2,j,i) = mapline (2,j)
          map(3,j,i) = mapline (3,j)
        enddo
      enddo
      close (unit=lun)
      x_delta = 360.0 / float(nx)
      y_delta = 180.0 / float(ny)
      do i = 1, ny
        lats = 90.0 - (float(i-1) * y_delta) - y_delta
        if (lats.gt.20.0) then
          lats  = lats  * torad
          do j = 1, nx
            longs = -180.0 + (float(j-1) * x_delta) + x_delta
            longs = longs * torad
            xs =  2.0 * float(nxs)/2.81 * tan(pid4 - (lats/2.0) ) * sin(longs)
            ys = -2.0 * float(nxs)/2.81 * tan(pid4 - (lats/2.0) ) * cos(longs)
c-- Rotations here - remember in images, Y needs a flip - but that's done when creating the image - not here I'm currently guessing
            xr =   cos(phi) * xs - sin(phi) * ys
            yr =  (sin(phi) * xs + cos(phi) * ys)
            xs = xr
            ys = yr
c--
            ix = min(max(int(xs)+(nxs/2),1),nxs)
            iy = min(max(int(ys)+(nys/2),1),nys)
c-- Rotations here
            maps(1,ix,iy) = map(1,j,i)
            maps(2,ix,iy) = map(2,j,i)
            maps(3,ix,iy) = map(3,j,i)
c            write (*,'(2I5,2F10.3,2I7)') i, j, lats*(180.0/3.141592654), longs*(180.0/3.141592654), int(xs), int(ys)
          enddo
c          read (*,*)
        endif
      enddo
c
      open (unit=lun, file=cmap//'.rgb', access='direct', form='unformatted', recl=3*nxs)
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
      call system('convert -depth 8 -size '//csize//' rgb:'//cmap//'.rgb -flip '//cmap//'.png')
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

