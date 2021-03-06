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
*    Unit 110     - test.elm        input 2-line element set file
*    Unit 111     - test.bak        output file
*    Unit 114     - sgp4test.dbg    debug output file
*    Unit 115     - sgp4rec.bak     temporary file of record for 2 line element sets
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
c      CLOSE(111)
      return
c
      entry sgp4_end()
      close (unit=110)
      close (unit=111)
      close (unit=114)
      close (unit=115)
      END
