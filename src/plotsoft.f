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
