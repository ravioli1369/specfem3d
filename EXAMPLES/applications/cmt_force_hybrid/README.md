# README

This example provides a simple guide for the simultaneous use of point-force and moment-tensor sources in SEM simulations.


1. **Enable Hybrid `FORCESOLUTION` and `CMTSOLUTION`**
   In `DATA/Par_file`, set:
   ```fortran
   USE_CMT_AND_FORCE_SOURCE = .true.`
The program will then read both `FORCESOLUTION` and `CMTSOLUTION` from the `DATA` directory and run a simulation using both point force and CMT sources.

2. **Use Binary Source Files**
when
```fortran
USE_BINARY_SOURCE_FILE  = .true.
```
is set in `DATA/Par_file`, all source attributes and source time functions must be written to DATA/SOLUTION.bin.
  ```fortran
  access='stream'
  ```
  in the `write` statement.
- When both options are enabled, the program ignores `FORCESOLUTION` and `CMTSOLUTION` and instead reads from SOLUTION.bin.
- The format of SOLUTION.bin is illustrated in the Fortran snippet below:
```fortran
integer, parameter :: IO = 10

!!! source attributes are either integer/double precision

! open with C-style binary
open(IO, file='DATA/SOLUTION.bin', form='UNFORMATTED', action='write', access='stream')

! write number of CMT and force sources
write(IO) NS_CMT, NS_FORCE

! loop over each CMT source
do isource = 1, NS_CMT
  write(IO) shift_c(isource), hdur_c(isource), lat_c(isource), lon_c(isource), depth_c(isource)
  write(IO) mt(:,isource)    ! Mrr Mtt Mpp Mrt Mrp Mtp
  if (USE_EXTERNAL_SOURCE_FILE) write(IO) external_stf_cmt(:,isource)
enddo

! loop over each point force source
do isource = 1, NS_FORCE
  write(IO) shift_f(isource), hdur_f(isource), lat_f(isource), lon_f(isource), depth_f(isource)
  write(IO) force_stf(isource)
  write(IO) factor_f(isource), fe(isource), fn(isource), fzup(isource)
  if (USE_EXTERNAL_SOURCE_FILE) write(IO) external_stf_force(:,isource)
enddo
```
you can also refer to `write_binary_source.f90` for a small example.

3. **Use an External Source Time Function**
In `DATA/Par_file`, enable:
```fortran
USE_EXTERNAL_SOURCE_FILE = .true.
```