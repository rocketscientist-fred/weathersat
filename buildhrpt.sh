#!/bin/bash
if [ $MMODEL = "LARGE" ] 
then
  echo "gfortran -c -ffixed-line-length-350 -fno-range-check -mcmodel=large -m64 -malign-double ./src/plotsoft.f"
  gfortran -c -ffixed-line-length-350 -fno-range-check -mcmodel=large -m64 -malign-double ./src/plotsoft.f
  echo "gfortran -c -ffixed-line-length-350 -fno-range-check -mcmodel=large -m64 -malign-double -O2 ./src/sgp4tle.f"
  gfortran -c -ffixed-line-length-350 -fno-range-check -mcmodel=large -m64 -malign-double -O2 ./src/sgp4tle.f
  echo "gfortran -o hrpt.exe -cpp -DLARGE -ffixed-line-length-350 -fno-range-check -mcmodel=large -m64 -malign-double -O2 ./src/hrpt.f plotsoft.o sgp4tle.o"
  gfortran -o hrpt.exe -cpp -DLARGE -ffixed-line-length-350 -fno-range-check -mcmodel=large -m64 -malign-double -O2 ./src/hrpt.f plotsoft.o sgp4tle.o
  echo "gfortran -o tletrack.exe -ffixed-line-length-350 -fno-range-check -mcmodel=large -m64 -malign-double -O2 ./src/tletrack.f plotsoft.o sgp4tle.o"
  gfortran -o tletrack.exe -ffixed-line-length-350 -fno-range-check -mcmodel=large -m64 -malign-double -O2 ./src/tletrack.f plotsoft.o sgp4tle.o
fi
if [ $MMODEL = "MEDIUM" ] 
then
  echo "gfortran -c -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double ./src/plotsoft.f"
  gfortran -c -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double ./src/plotsoft.f
  echo "gfortran -c -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 ./src/sgp4tle.f"
  gfortran -c -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 ./src/sgp4tle.f
  echo "gfortran -o hrpt.exe -cpp -DMEDIUM -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 ./src/hrpt.f plotsoft.o sgp4tle.o"
  gfortran -o hrpt.exe -cpp -DMEDIUM -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 ./src/hrpt.f plotsoft.o sgp4tle.o
  echo "gfortran -o tletrack.exe -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 ./src/tletrack.f plotsoft.o sgp4tle.o"
  gfortran -o tletrack.exe -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 ./src/tletrack.f plotsoft.o sgp4tle.o
fi
echo "gfortran -o readbin.exe -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 ./src/readbin.f"
gfortran -o readbin.exe -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 ./src/readbin.f
echo "gfortran -o readbin_modis.exe -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 ./src/readbin_modis.f plotsoft.o"
gfortran -o readbin_modis.exe -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 ./src/readbin_modis.f plotsoft.o
echo "gfortran -o polarstereo.exe -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 ./src/polarstereo.f"
gfortran -o polarstereo.exe -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 ./src/polarstereo.f
echo "gfortran -o deproject.exe -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 ./src/deproject.f plotsoft.o"
gfortran -o deproject.exe -ffixed-line-length-350 -fno-range-check -mcmodel=medium -m64 -malign-double -O2 ./src/deproject.f plotsoft.o
rm *.o
