subroutine linspace(n,xmin,xmax,arr)
  implicit none
  integer,intent(in) :: n
  real,intent(in) :: xmin,xmax
  double precision,intent(inout) :: arr(n)

  !local
  integer :: i

  if (n == 1) then
    arr(1) = xmin
  else
    do i = 0,n-1
      arr(i + 1) = xmin + (xmax - xmin) / (n-1) * i
    enddo
  endif

end subroutine linspace

program main
  implicit none

  !!!!!!!! ATTENTION HERE !!!!!!!!!!!!!!!!!!!!!
  ! all parameters are double precision except:
  ! integer :: NS_CMT NS_FORCE  no. of CMT/FORCE sources
  ! integer :: force_stf(NS_Force) source time function flag in FORCESOLUTION
  ! real ::  external_stf_force/cmt, with shape(NSTEP,NS_FORCE/CMT)
  ! in binary file, you should first dump CMT, and then force

  ! no. of cmt and force sources
  integer,parameter :: NS_CMT = 10, NS_FORCE = 5
  integer,parameter :: IO = 123, IO_TXT = 12114 ! FILE number
  integer,parameter ::  NSTEP= 3000
  integer :: isource,it,i
  logical :: USE_EXTERNAL_SOURCE = .true.

  ! allocate space for force
  integer,dimension(NS_FORCE) :: force_stf ! source time function flag
  double precision,dimension(NS_FORCE) :: lat_f,lon_f,depth_f,hdur_f,shift_f, &
                                          factor_f,fe,fn,fzup
  real :: external_stf_force(NSTEP,NS_FORCE)

  ! allocate space for cmt
  double precision,dimension(NS_CMT) :: lat_c,lon_c,depth_c,hdur_c,shift_c
  double precision, dimension(6,NS_CMT) :: mt ! Mrr Mtt Mpp Mrt Mrp Mtp
  real :: external_stf_cmt(NSTEP,NS_CMT)

  !filename
  character(len=128) :: filename

  ! set value for FORCE, same as FORCESOLUTION
  force_stf = 0
  factor_f = 2.0d12
  fe = -1.0; fn = 2.0; fzup = 1.0
  lat_f = 0; depth_f = 1.
  hdur_f = 0.5; shift_f = 0.
  call linspace(NS_FORCE,-512.54,541.5,lon_f)

  do isource = 1,NS_FORCE
    write(filename,'(a,i0,a)') './DATA/stf.force.',isource-1,'.txt'
    open(IO_TXT,file=trim(filename))
    do it = 1,NSTEP
      read(IO_TXT,*) external_stf_force(it,isource)
    enddo;
    close(IO_TXT)
  enddo

  ! set value for CMT, same as CMTSOLUTION
  mt(:,1) = (/1.0,1.0,0.0,1.,0.,0./)
  mt = mt * 1.0d23
  lon_c = 0.; depth_c = 0.5
  hdur_c = 5.0; shift_c = 0.
  call linspace(NS_CMT,-500.,500.,lat_c)

  do i = 1,NS_CMT
    lat_c(i) = -523.131 + 1000. / (NS_CMT - 1) * (i-1)
  enddo

  do isource = 1,NS_CMT
    write(filename,'(a,i0,a)') './DATA/stf.cmt.',isource-1,'.txt'
    open(IO_TXT,file=trim(filename))
    do it = 1,NSTEP
      read(IO_TXT,*) external_stf_cmt(it,isource)
    enddo;
    close(IO_TXT)
  enddo

  ! write to binary file and txt file
  open(IO,file='DATA/SOLUTION.bin',form='UNFORMATTED',action='write',access='stream')
  write(IO) NS_CMT,NS_FORCE

  ! write cmt
  open(IO_TXT,file='./DATA/CMTSOLUTION',action='write')
  do isource = 1,NS_CMT
    write(IO) shift_c(isource),hdur_c(isource),lat_c(isource),lon_c(isource),depth_c(isource)
    write(IO) mt(:,isource) !  Mrr Mtt Mpp Mrt Mrp Mtp
    if (USE_EXTERNAL_SOURCE) write(IO) external_stf_cmt(:,isource)


    ! txt file
    write(IO_TXT,*)'PDE  1999 01 01 00 00 00.00  67000 67000 -25000 4.2 4.2 hom_explosion'
    write(IO_TXT,*) 'event name:       hom_explosion'
    write(IO_TXT,*) 'time shift:', shift_c(isource)
    write(IO_TXT,*) 'half duration:',   hdur_c(isource)
    write(IO_TXT,*) 'latorUTM:',       lat_c(isource)
    write(IO_TXT,*) 'longorUTM:',      lon_c(isource)
    write(IO_TXT,*) 'depth:',           depth_c(isource)
    write(IO_TXT,*) 'Mrr:',  mt(1,isource)
    write(IO_TXT,*) 'Mtt:',  mt(2,isource)
    write(IO_TXT,*) 'Mpp:',  mt(3,isource)
    write(IO_TXT,*) 'Mrt:',  mt(4,isource)
    write(IO_TXT,*) 'Mrp:',  mt(5,isource)
    write(IO_TXT,*) 'Mtp:',  mt(6,isource)
    write(filename,'(a,i0,a)') './DATA/stf.cmt.',isource-1,'.txt'
    write(IO_TXT,*) trim(filename)
  enddo
  close(IO_TXT)

  ! write force
  open(IO_TXT,file='./DATA/FORCESOLUTION',action='write')
  do isource = 1,NS_FORCE
    write(IO) shift_f(isource),hdur_f(isource),lat_f(isource),lon_f(isource),depth_f(isource)
    write(IO) force_stf(isource)
    write(IO) factor_f(isource),fe(isource),fn(isource),fzup(isource)
    if (USE_EXTERNAL_SOURCE) write(IO) external_stf_force(:,isource)

    !txt file
    write(IO_TXT,*) 'FORCE  001'
    write(IO_TXT,*) 'time shift:',    shift_f(isource)
    write(IO_TXT,*) 'hdurorf0:', hdur_f(isource)
    write(IO_TXT,*) 'latorUTM:',      lat_f(isource)
    write(IO_TXT,*) 'longorUTM:',  lon_f(isource)
    write(IO_TXT,*) 'depth:',         depth_f(isource)
    write(IO_TXT,*) 'source time function:',            force_stf(isource)
    write(IO_TXT,*) 'factor force source:',             factor_f(isource)
    write(IO_TXT,*) 'component dir vect source E:',     fe(isource)
    write(IO_TXT,*) 'component dir vect source N:',     fn(isource)
    write(IO_TXT,*) 'component dir vect source Z_UP:',  fzup(isource)
    write(filename,'(a,i0,a)') './DATA/stf.force.',isource-1,'.txt'
    write(IO_TXT,*) trim(filename)
  enddo
  close(IO)
  close(IO_TXT)

end program main