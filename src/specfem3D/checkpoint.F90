!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

module checkpoint_par

  use constants, only: CUSTOM_REAL

  implicit none

  integer, parameter :: IOUT_CKPT = 245
  integer, parameter :: CKPT_MAGIC = 314159265
  integer, parameter :: CKPT_MAGIC_END = 271828182
  integer, parameter :: CKPT_VERSION = 1

  ! alternating slots so a kill during a checkpoint write leaves the previous one intact
  integer :: ckpt_slot = 0

  ! resume state, set by checkpoint_check_restart()
  logical :: RESUME_FROM_CHECKPOINT = .false.
  integer :: it_resume = 0

  ! the device pressure seismogram buffer is packed differently from the host array
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: ckpt_seismo_p_gpu
  integer :: ckpt_size_dva = 0, ckpt_size_p = 0

end module checkpoint_par

!
!-------------------------------------------------------------------------------------------------
!

  subroutine checkpoint_filename(islot,fname)

  use constants, only: MAX_STRING_LEN
  use specfem_par, only: prname

  implicit none

  integer,intent(in) :: islot
  character(len=MAX_STRING_LEN),intent(out) :: fname

  write(fname,'(a,a,i1.1,a)') prname(1:len_trim(prname)),'checkpoint_',islot,'.bin'

  end subroutine checkpoint_filename

!
!-------------------------------------------------------------------------------------------------
!

  subroutine checkpoint_okname(islot,fname)

  use constants, only: MAX_STRING_LEN
  use specfem_par, only: prname

  implicit none

  integer,intent(in) :: islot
  character(len=MAX_STRING_LEN),intent(out) :: fname

  write(fname,'(a,a,i1.1,a)') prname(1:len_trim(prname)),'checkpoint_',islot,'.ok'

  end subroutine checkpoint_okname

!
!-------------------------------------------------------------------------------------------------
!

  subroutine checkpoint_check_supported()

! refuses configurations whose per-step state this checkpoint does not capture

  use constants, only: INJECTION_TECHNIQUE_IS_DSM,INJECTION_TECHNIQUE_IS_AXISEM
  use specfem_par
  use specfem_par_coupling, only: do_save_coupling_wavefield
  use pml_par, only: NSPEC_CPML
  use checkpoint_par

  implicit none

  integer :: ier

  if (NTSTEP_BETWEEN_CHECKPOINTS <= 0) return

  ! GPU seismogram buffer sizes, matching the device allocations
  if (GPU_MODE .and. do_save_seismograms .and. nrec_local > 0) then
    ckpt_size_dva = NDIM * nrec_local * nlength_seismogram
    ckpt_size_p = nrec_local * nlength_seismogram * NB_RUNS_ACOUSTIC_GPU
    allocate(ckpt_seismo_p_gpu(ckpt_size_p),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating checkpoint GPU pressure buffer')
    ckpt_seismo_p_gpu(:) = 0._CUSTOM_REAL
  endif

  if (UNDO_ATTENUATION_AND_OR_PML) &
    call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported with UNDO_ATTENUATION_AND_OR_PML')
  if (LTS_MODE) &
    call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported with LTS_MODE')
  if (EXACT_UNDOING_TO_DISK) &
    call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported with EXACT_UNDOING_TO_DISK')
  if (SIMULATION_TYPE == 2) &
    call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported for SIMULATION_TYPE == 2')
  if (SIMULATION_TYPE == 3 .and. READ_ADJSRC_ASDF) &
    call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported for SIMULATION_TYPE == 3 with READ_ADJSRC_ASDF')
  if (SIMULATION_TYPE == 3 .and. SAVE_MOHO_MESH) &
    call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported for SIMULATION_TYPE == 3 with SAVE_MOHO_MESH')
  if (NOISE_TOMOGRAPHY /= 0) &
    call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported with NOISE_TOMOGRAPHY')
  if (ASDF_FORMAT) &
    call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported with ASDF_FORMAT seismograms')
  if (do_save_seismograms) then
    if (USE_BINARY_FOR_SEISMOGRAMS) &
      call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS needs ASCII seismograms, set USE_BINARY_FOR_SEISMOGRAMS = .false.')
    if (SAVE_ALL_SEISMOS_IN_ONE_FILE) &
      call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported with SAVE_ALL_SEISMOS_IN_ONE_FILE')
    if (HDF5_FORMAT) &
      call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported with HDF5_FORMAT seismograms')
  endif
  if (ADIOS_ENABLED) &
    call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported with ADIOS_ENABLED')
  if (POROELASTIC_SIMULATION .and. GPU_MODE) &
    call exit_MPI(myrank,'poroelastic GPU runs are not supported at all')
  if (USE_LDDRK .and. ELASTIC_SIMULATION .and. GPU_MODE) &
    call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported with USE_LDDRK on GPU')
  if (PML_CONDITIONS .and. GPU_MODE) &
    call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported with PML_CONDITIONS on GPU')
  if (RECIPROCITY_AND_KH_INTEGRAL) &
    call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported with RECIPROCITY_AND_KH_INTEGRAL')
  if (do_save_coupling_wavefield) &
    call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported while saving a coupling wavefield')

  ! DSM carries it_dsm plus two sequential file positions; AxiSEM reads its boundary file sequentially
  if (COUPLE_WITH_INJECTION_TECHNIQUE) then
    if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_DSM) &
      call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported for DSM injection')
    if (INJECTION_TECHNIQUE_TYPE == INJECTION_TECHNIQUE_IS_AXISEM) &
      call exit_MPI(myrank,'NTSTEP_BETWEEN_CHECKPOINTS not supported for AxiSEM injection')
  endif

  end subroutine checkpoint_check_supported

!
!-------------------------------------------------------------------------------------------------
!

  subroutine checkpoint_check_restart()

! scans both slots, picks the newest one all ranks agree on

  use constants, only: MAX_STRING_LEN,IMAIN
  use specfem_par
  use checkpoint_par

  implicit none

  integer :: islot,ier,islot_best,irank,icand,it_cand,it_newest
  integer :: mysteps(2)
  integer, dimension(:,:), allocatable :: allsteps
  integer :: ihead(4),itail(2)
  character(len=MAX_STRING_LEN) :: fname,okname
  logical :: exists,common_to_all

  RESUME_FROM_CHECKPOINT = .false.
  it_resume = 0
  ckpt_slot = 0

  if (NTSTEP_BETWEEN_CHECKPOINTS <= 0) return

  mysteps(:) = 0
  islot_best = -1

  do islot = 0,1
    call checkpoint_filename(islot,fname)
    inquire(file=trim(fname),exist=exists)
    if (.not. exists) cycle

    open(unit=IOUT_CKPT,file=trim(fname),status='old',form='unformatted',action='read',iostat=ier)
    if (ier /= 0) cycle

    ! header: magic, version, nproc, it
    read(IOUT_CKPT,iostat=ier) ihead
    if (ier /= 0) then
      close(IOUT_CKPT)
      cycle
    endif
    close(IOUT_CKPT)

    if (ihead(1) /= CKPT_MAGIC) cycle
    if (ihead(2) /= CKPT_VERSION) cycle
    if (ihead(3) /= NPROC) cycle

    ! the .ok marker is written only after the payload file is closed
    call checkpoint_okname(islot,okname)
    inquire(file=trim(okname),exist=exists)
    if (.not. exists) cycle
    open(unit=IOUT_CKPT,file=trim(okname),status='old',form='unformatted',action='read',iostat=ier)
    if (ier /= 0) cycle
    read(IOUT_CKPT,iostat=ier) itail
    close(IOUT_CKPT)
    if (ier /= 0) cycle
    if (itail(1) /= CKPT_MAGIC_END) cycle
    if (itail(2) /= ihead(4)) cycle

    mysteps(islot+1) = ihead(4)
  enddo

  ! a kill can land between one rank writing its marker and another writing its own,
  ! so resume from the newest step every rank still holds rather than giving up
  allocate(allsteps(2,0:NPROC-1),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating checkpoint step table')
  call gather_all_all_i(mysteps,2,allsteps,2,NPROC)

  it_resume = 0
  it_newest = maxval(allsteps)

  do icand = 1,2
    it_cand = mysteps(icand)
    if (it_cand <= 0) cycle
    if (it_cand <= it_resume) cycle
    common_to_all = .true.
    do irank = 0,NPROC-1
      if (allsteps(1,irank) /= it_cand .and. allsteps(2,irank) /= it_cand) then
        common_to_all = .false.
        exit
      endif
    enddo
    if (common_to_all) then
      it_resume = it_cand
      islot_best = icand - 1
    endif
  enddo
  deallocate(allsteps)

  if (it_resume <= 0) then
    if (myrank == 0 .and. it_newest > 0) then
      write(IMAIN,*) 'Checkpoint: no step is present on every rank (newest seen ',it_newest,'), starting from scratch'
      call flush_IMAIN()
    endif
    it_resume = 0
    return
  endif

  if (it_resume >= NSTEP) then
    if (myrank == 0) then
      write(IMAIN,*) 'Checkpoint: found a completed run at step ',it_resume,', nothing to do'
      call flush_IMAIN()
    endif
    it_resume = 0
    return
  endif

  if (myrank == 0 .and. it_resume < it_newest) then
    write(IMAIN,*) 'Checkpoint: newest step seen is ',it_newest,' but only ',it_resume,' is on every rank'
    call flush_IMAIN()
  endif

  RESUME_FROM_CHECKPOINT = .true.
  ckpt_slot = islot_best

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Checkpoint: resuming from step ',it_resume,' of ',NSTEP
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine checkpoint_check_restart


!
!-------------------------------------------------------------------------------------------------
!

  subroutine checkpoint_open_abs_w(fid,fname,flen,filesize)

! absorbing-record files must not be truncated when resuming

  use checkpoint_par, only: RESUME_FROM_CHECKPOINT

  implicit none

  integer,intent(in) :: fid,flen
  character(len=*),intent(in) :: fname
  integer(kind=8),intent(in) :: filesize

  if (RESUME_FROM_CHECKPOINT) then
    call open_file_abs_rw(fid,fname,flen,filesize)
  else
    call open_file_abs_w(fid,fname,flen,filesize)
  endif

  end subroutine checkpoint_open_abs_w

!
!-------------------------------------------------------------------------------------------------
!

  subroutine checkpoint_gpu_to_host()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use checkpoint_par

  implicit none

  if (.not. GPU_MODE) return

  if (ACOUSTIC_SIMULATION) &
    call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic,potential_dot_acoustic, &
                                        potential_dot_dot_acoustic,Mesh_pointer)

  if (ELASTIC_SIMULATION) then
    call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ,veloc,accel,Mesh_pointer)
    if (ATTENUATION) then
      call transfer_rmemory_from_device(Mesh_pointer,R_xx,R_yy,R_xy,R_xz,R_yz,R_trace,size(R_xx,kind=4))
      call transfer_strain_from_device(Mesh_pointer,epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                       epsilondev_xz,epsilondev_yz,epsilondev_trace,size(epsilondev_xx,kind=4))
    endif
  endif

  ! backward wavefields and kernel accumulators live only on the device during a kernel run
  if (SIMULATION_TYPE == 3) then
    if (ACOUSTIC_SIMULATION) then
      call transfer_b_fields_ac_from_device(NGLOB_AB,b_potential_acoustic,b_potential_dot_acoustic, &
                                            b_potential_dot_dot_acoustic,Mesh_pointer)
      call transfer_kernels_ac_to_host(Mesh_pointer,rho_ac_kl,kappa_ac_kl,NSPEC_AB)
      if (APPROXIMATE_HESS_KL) &
        call transfer_kernels_hess_ac_tohost(Mesh_pointer,hess_ac_kl,hess_rho_ac_kl,hess_kappa_ac_kl,NSPEC_AB)
    endif
    if (ELASTIC_SIMULATION) then
      call transfer_b_fields_from_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel,Mesh_pointer)
      call transfer_b_strain_from_device(Mesh_pointer,b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy, &
                                         b_epsilondev_xz,b_epsilondev_yz,b_epsilondev_trace, &
                                         size(b_epsilondev_xx,kind=4))
      call transfer_b_eps_trace_from_device(Mesh_pointer,b_epsilon_trace_over_3, &
                                            size(b_epsilon_trace_over_3,kind=4))
      if (ATTENUATION) &
        call transfer_b_rmemory_from_device(Mesh_pointer,b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz, &
                                            b_R_trace,size(b_R_xx,kind=4))
      call transfer_kernels_el_to_host(Mesh_pointer,rho_kl,mu_kl,kappa_kl,cijkl_kl,NSPEC_AB)
      if (APPROXIMATE_HESS_KL) &
        call transfer_kernels_hess_el_tohost(Mesh_pointer,hess_kl,hess_rho_kl,hess_kappa_kl,hess_mu_kl,NSPEC_AB)
    endif
  endif

  ! seismogram samples accumulate on the device between flushes
  if (ckpt_size_dva > 0) &
    call transfer_seismograms_from_device(Mesh_pointer,seismograms_d,seismograms_v,seismograms_a, &
                                          ckpt_seismo_p_gpu,ckpt_size_dva,ckpt_size_p)

  end subroutine checkpoint_gpu_to_host

!
!-------------------------------------------------------------------------------------------------
!

  subroutine checkpoint_host_to_gpu()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use checkpoint_par

  implicit none

  if (.not. GPU_MODE) return

  if (ACOUSTIC_SIMULATION) &
    call transfer_fields_ac_to_device(NGLOB_AB,potential_acoustic,potential_dot_acoustic, &
                                      potential_dot_dot_acoustic,Mesh_pointer)

  if (ELASTIC_SIMULATION) then
    call transfer_fields_el_to_device(NDIM*NGLOB_AB,displ,veloc,accel,Mesh_pointer)
    if (ATTENUATION) then
      call transfer_rmemory_to_device(Mesh_pointer,R_xx,R_yy,R_xy,R_xz,R_yz,R_trace,size(R_xx,kind=4))
      call transfer_strain_to_device(Mesh_pointer,epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                     epsilondev_xz,epsilondev_yz,epsilondev_trace,size(epsilondev_xx,kind=4))
    endif
  endif

  if (SIMULATION_TYPE == 3) then
    if (ACOUSTIC_SIMULATION) then
      call transfer_b_fields_ac_to_device(NGLOB_AB,b_potential_acoustic,b_potential_dot_acoustic, &
                                          b_potential_dot_dot_acoustic,Mesh_pointer)
      call transfer_kernels_ac_to_device(Mesh_pointer,rho_ac_kl,kappa_ac_kl,NSPEC_AB)
      if (APPROXIMATE_HESS_KL) &
        call transfer_kernels_hess_ac_todevice(Mesh_pointer,hess_ac_kl,hess_rho_ac_kl,hess_kappa_ac_kl,NSPEC_AB)
    endif
    if (ELASTIC_SIMULATION) then
      call transfer_b_fields_to_device(NDIM*NGLOB_AB,b_displ,b_veloc,b_accel,Mesh_pointer)
      call transfer_b_strain_to_device(Mesh_pointer,b_epsilondev_xx,b_epsilondev_yy,b_epsilondev_xy, &
                                       b_epsilondev_xz,b_epsilondev_yz,b_epsilondev_trace, &
                                       size(b_epsilondev_xx,kind=4))
      call transfer_b_eps_trace_to_device(Mesh_pointer,b_epsilon_trace_over_3, &
                                          size(b_epsilon_trace_over_3,kind=4))
      if (ATTENUATION) &
        call transfer_b_rmemory_to_device(Mesh_pointer,b_R_xx,b_R_yy,b_R_xy,b_R_xz,b_R_yz, &
                                          b_R_trace,size(b_R_xx,kind=4))
      call transfer_kernels_el_to_device(Mesh_pointer,rho_kl,mu_kl,kappa_kl,cijkl_kl,NSPEC_AB)
      if (APPROXIMATE_HESS_KL) &
        call transfer_kernels_hess_el_todevice(Mesh_pointer,hess_kl,hess_rho_kl,hess_kappa_kl,hess_mu_kl,NSPEC_AB)
    endif
  endif

  if (ckpt_size_dva > 0) &
    call transfer_seismograms_to_device(Mesh_pointer,seismograms_d,seismograms_v,seismograms_a, &
                                        ckpt_seismo_p_gpu,ckpt_size_dva,ckpt_size_p)

  end subroutine checkpoint_host_to_gpu

!
!-------------------------------------------------------------------------------------------------
!

  subroutine checkpoint_save()

! writes the full restartable state of this rank at the end of time step it

  use constants, only: MAX_STRING_LEN,IMAIN
  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use pml_par
  use gravity_perturbation, only: GRAVITY_SIMULATION,gravity_checkpoint_save
  use checkpoint_par

  implicit none

  integer :: ier,islot
  character(len=MAX_STRING_LEN) :: fname,okname
  logical :: exists

  ! alternate slots so the previous checkpoint survives a kill during this write
  islot = 1 - ckpt_slot

  call checkpoint_gpu_to_host()

  ! invalidate the slot we are about to overwrite
  call checkpoint_okname(islot,okname)
  inquire(file=trim(okname),exist=exists)
  if (exists) then
    open(unit=IOUT_CKPT,file=trim(okname),status='old',iostat=ier)
    if (ier == 0) close(IOUT_CKPT,status='delete')
  endif

  call checkpoint_filename(islot,fname)
  open(unit=IOUT_CKPT,file=trim(fname),status='replace',form='unformatted',action='write',iostat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error opening checkpoint file for writing')

  write(IOUT_CKPT) CKPT_MAGIC,CKPT_VERSION,NPROC,it
  write(IOUT_CKPT) NSTEP,NGLOB_AB,NSPEC_AB,SIMULATION_TYPE,seismo_offset,seismo_current,nlength_seismogram,nrec_local, &
                   ckpt_size_p,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC,num_abs_boundary_faces

  call checkpoint_body(.true.)

  close(IOUT_CKPT)

  ! marker written only once the payload is on disk and closed
  open(unit=IOUT_CKPT,file=trim(okname),status='replace',form='unformatted',action='write',iostat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error opening checkpoint marker file for writing')
  write(IOUT_CKPT) CKPT_MAGIC_END,it
  close(IOUT_CKPT)

  ! only flip once the new slot is complete
  call synchronize_all()
  ckpt_slot = islot

  if (myrank == 0) then
    write(IMAIN,*) 'Checkpoint written at time step ',it
    call flush_IMAIN()
  endif

  end subroutine checkpoint_save

!
!-------------------------------------------------------------------------------------------------
!

  subroutine checkpoint_restore()

! restores the state written by checkpoint_save() and sets it_begin accordingly

  use constants, only: MAX_STRING_LEN
  use specfem_par
  use checkpoint_par

  implicit none

  integer :: ier
  integer :: ihead(4),ibody(12)
  character(len=MAX_STRING_LEN) :: fname

  if (.not. RESUME_FROM_CHECKPOINT) return

  call checkpoint_filename(ckpt_slot,fname)
  open(unit=IOUT_CKPT,file=trim(fname),status='old',form='unformatted',action='read',iostat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error opening checkpoint file for reading')

  read(IOUT_CKPT) ihead
  read(IOUT_CKPT) ibody

  if (ihead(4) /= it_resume) call exit_MPI(myrank,'Checkpoint step mismatch on restore')
  if (ibody(1) /= NSTEP) call exit_MPI(myrank,'Checkpoint NSTEP mismatch: rerun with the same Par_file')
  if (ibody(2) /= NGLOB_AB) call exit_MPI(myrank,'Checkpoint NGLOB_AB mismatch: mesh changed')
  if (ibody(3) /= NSPEC_AB) call exit_MPI(myrank,'Checkpoint NSPEC_AB mismatch: mesh changed')
  if (ibody(4) /= SIMULATION_TYPE) call exit_MPI(myrank,'Checkpoint SIMULATION_TYPE mismatch')
  if (ibody(7) /= nlength_seismogram) call exit_MPI(myrank,'Checkpoint nlength_seismogram mismatch')
  if (ibody(8) /= nrec_local) call exit_MPI(myrank,'Checkpoint nrec_local mismatch')
  if (ibody(9) /= ckpt_size_p) call exit_MPI(myrank,'Checkpoint GPU_MODE/seismogram layout mismatch')
  if (ibody(10) /= nadj_rec_local) call exit_MPI(myrank,'Checkpoint nadj_rec_local mismatch')
  if (ibody(11) /= NTSTEP_BETWEEN_READ_ADJSRC) call exit_MPI(myrank,'Checkpoint NTSTEP_BETWEEN_READ_ADJSRC mismatch')
  ! the Stacey record grid must match the forward run that wrote absorb_*.bin
  if (ibody(12) /= num_abs_boundary_faces) call exit_MPI(myrank,'Checkpoint num_abs_boundary_faces mismatch')

  seismo_offset = ibody(5)
  seismo_current = ibody(6)

  call checkpoint_body(.false.)

  close(IOUT_CKPT)

  call checkpoint_host_to_gpu()

  call checkpoint_trim_seismograms()

  it_begin = it_resume + 1

  end subroutine checkpoint_restore

!
!-------------------------------------------------------------------------------------------------
!

  subroutine checkpoint_trim_seismograms()

! a flush may have landed after the last checkpoint, or been cut short by the kill;
! drop everything on disk past the checkpointed sample cursor so the next flush appends cleanly

  use constants, only: MAX_STRING_LEN,IMAIN
  use specfem_par, only: myrank,seismo_offset,do_save_seismograms,OUTPUT_FILES,SU_FORMAT

  implicit none

  character(len=MAX_STRING_LEN) :: cmd
  character(len=16) :: nstr
  integer :: ier

  if (.not. do_save_seismograms) return

  ! SU writes are positioned by absolute byte offset, so a resumed run overwrites in place
  if (SU_FORMAT) return

  if (myrank == 0) then
    write(nstr,'(i16)') seismo_offset
    write(cmd,'(a)') 'for f in '//trim(OUTPUT_FILES)//'/*.sem*; do [ -f "$f" ] || continue; '// &
                     'head -n '//trim(adjustl(nstr))//' "$f" > "$f.ckpttmp" && mv "$f.ckpttmp" "$f"; done'
    call execute_command_line(trim(cmd),exitstat=ier)
    if (ier /= 0) then
      write(IMAIN,*) 'Checkpoint: warning, could not trim seismogram files (status ',ier,')'
      call flush_IMAIN()
    endif
  endif

  call synchronize_all()

  end subroutine checkpoint_trim_seismograms

!
!-------------------------------------------------------------------------------------------------
!

  subroutine checkpoint_body(do_write)

! single description of the payload, so the write and read paths cannot drift apart

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic
  use pml_par
  use gravity_perturbation, only: GRAVITY_SIMULATION,gravity_checkpoint_save,gravity_checkpoint_restore
  use checkpoint_par, only: IOUT_CKPT,ckpt_seismo_p_gpu,ckpt_size_p

  implicit none

  logical,intent(in) :: do_write

  if (ACOUSTIC_SIMULATION) then
    call ckpt_arr(do_write,potential_acoustic,size(potential_acoustic))
    call ckpt_arr(do_write,potential_dot_acoustic,size(potential_dot_acoustic))
    call ckpt_arr(do_write,potential_dot_dot_acoustic,size(potential_dot_dot_acoustic))
    if (USE_LDDRK) then
      call ckpt_arr(do_write,potential_acoustic_lddrk,size(potential_acoustic_lddrk))
      call ckpt_arr(do_write,potential_dot_acoustic_lddrk,size(potential_dot_acoustic_lddrk))
    endif
  endif

  if (ELASTIC_SIMULATION) then
    call ckpt_arr(do_write,displ,size(displ))
    call ckpt_arr(do_write,veloc,size(veloc))
    call ckpt_arr(do_write,accel,size(accel))
    if (ATTENUATION) then
      call ckpt_arr(do_write,R_trace,size(R_trace))
      call ckpt_arr(do_write,R_xx,size(R_xx))
      call ckpt_arr(do_write,R_yy,size(R_yy))
      call ckpt_arr(do_write,R_xy,size(R_xy))
      call ckpt_arr(do_write,R_xz,size(R_xz))
      call ckpt_arr(do_write,R_yz,size(R_yz))
      call ckpt_arr(do_write,epsilondev_trace,size(epsilondev_trace))
      call ckpt_arr(do_write,epsilondev_xx,size(epsilondev_xx))
      call ckpt_arr(do_write,epsilondev_yy,size(epsilondev_yy))
      call ckpt_arr(do_write,epsilondev_xy,size(epsilondev_xy))
      call ckpt_arr(do_write,epsilondev_xz,size(epsilondev_xz))
      call ckpt_arr(do_write,epsilondev_yz,size(epsilondev_yz))
    endif
    if (USE_LDDRK) then
      call ckpt_arr(do_write,displ_lddrk,size(displ_lddrk))
      call ckpt_arr(do_write,veloc_lddrk,size(veloc_lddrk))
      if (ATTENUATION) then
        call ckpt_arr(do_write,R_trace_lddrk,size(R_trace_lddrk))
        call ckpt_arr(do_write,R_xx_lddrk,size(R_xx_lddrk))
        call ckpt_arr(do_write,R_yy_lddrk,size(R_yy_lddrk))
        call ckpt_arr(do_write,R_xy_lddrk,size(R_xy_lddrk))
        call ckpt_arr(do_write,R_xz_lddrk,size(R_xz_lddrk))
        call ckpt_arr(do_write,R_yz_lddrk,size(R_yz_lddrk))
      endif
    endif
  endif

  if (POROELASTIC_SIMULATION) then
    call ckpt_arr(do_write,displs_poroelastic,size(displs_poroelastic))
    call ckpt_arr(do_write,velocs_poroelastic,size(velocs_poroelastic))
    call ckpt_arr(do_write,accels_poroelastic,size(accels_poroelastic))
    call ckpt_arr(do_write,displw_poroelastic,size(displw_poroelastic))
    call ckpt_arr(do_write,velocw_poroelastic,size(velocw_poroelastic))
    call ckpt_arr(do_write,accelw_poroelastic,size(accelw_poroelastic))
  endif

  ! PML convolution memory variables
  if (PML_CONDITIONS .and. NSPEC_CPML > 0) then
    if (allocated(rmemory_displ_elastic)) &
      call ckpt_arr(do_write,rmemory_displ_elastic,size(rmemory_displ_elastic))
    if (allocated(rmemory_dux_dxl_x)) &
      call ckpt_arr(do_write,rmemory_dux_dxl_x,size(rmemory_dux_dxl_x))
    if (allocated(rmemory_dux_dyl_x)) &
      call ckpt_arr(do_write,rmemory_dux_dyl_x,size(rmemory_dux_dyl_x))
    if (allocated(rmemory_dux_dzl_x)) &
      call ckpt_arr(do_write,rmemory_dux_dzl_x,size(rmemory_dux_dzl_x))
    if (allocated(rmemory_duy_dxl_x)) &
      call ckpt_arr(do_write,rmemory_duy_dxl_x,size(rmemory_duy_dxl_x))
    if (allocated(rmemory_duy_dyl_x)) &
      call ckpt_arr(do_write,rmemory_duy_dyl_x,size(rmemory_duy_dyl_x))
    if (allocated(rmemory_duz_dxl_x)) &
      call ckpt_arr(do_write,rmemory_duz_dxl_x,size(rmemory_duz_dxl_x))
    if (allocated(rmemory_duz_dzl_x)) &
      call ckpt_arr(do_write,rmemory_duz_dzl_x,size(rmemory_duz_dzl_x))
    if (allocated(rmemory_dux_dxl_y)) &
      call ckpt_arr(do_write,rmemory_dux_dxl_y,size(rmemory_dux_dxl_y))
    if (allocated(rmemory_dux_dyl_y)) &
      call ckpt_arr(do_write,rmemory_dux_dyl_y,size(rmemory_dux_dyl_y))
    if (allocated(rmemory_duy_dxl_y)) &
      call ckpt_arr(do_write,rmemory_duy_dxl_y,size(rmemory_duy_dxl_y))
    if (allocated(rmemory_duy_dyl_y)) &
      call ckpt_arr(do_write,rmemory_duy_dyl_y,size(rmemory_duy_dyl_y))
    if (allocated(rmemory_duy_dzl_y)) &
      call ckpt_arr(do_write,rmemory_duy_dzl_y,size(rmemory_duy_dzl_y))
    if (allocated(rmemory_duz_dzl_y)) &
      call ckpt_arr(do_write,rmemory_duz_dzl_y,size(rmemory_duz_dzl_y))
    if (allocated(rmemory_duz_dyl_y)) &
      call ckpt_arr(do_write,rmemory_duz_dyl_y,size(rmemory_duz_dyl_y))
    if (allocated(rmemory_dux_dxl_z)) &
      call ckpt_arr(do_write,rmemory_dux_dxl_z,size(rmemory_dux_dxl_z))
    if (allocated(rmemory_dux_dzl_z)) &
      call ckpt_arr(do_write,rmemory_dux_dzl_z,size(rmemory_dux_dzl_z))
    if (allocated(rmemory_duy_dyl_z)) &
      call ckpt_arr(do_write,rmemory_duy_dyl_z,size(rmemory_duy_dyl_z))
    if (allocated(rmemory_duy_dzl_z)) &
      call ckpt_arr(do_write,rmemory_duy_dzl_z,size(rmemory_duy_dzl_z))
    if (allocated(rmemory_duz_dxl_z)) &
      call ckpt_arr(do_write,rmemory_duz_dxl_z,size(rmemory_duz_dxl_z))
    if (allocated(rmemory_duz_dyl_z)) &
      call ckpt_arr(do_write,rmemory_duz_dyl_z,size(rmemory_duz_dyl_z))
    if (allocated(rmemory_duz_dzl_z)) &
      call ckpt_arr(do_write,rmemory_duz_dzl_z,size(rmemory_duz_dzl_z))
    if (allocated(PML_displ_old)) &
      call ckpt_arr(do_write,PML_displ_old,size(PML_displ_old))
    if (allocated(PML_displ_new)) &
      call ckpt_arr(do_write,PML_displ_new,size(PML_displ_new))
    if (allocated(rmemory_potential_acoustic)) &
      call ckpt_arr(do_write,rmemory_potential_acoustic,size(rmemory_potential_acoustic))
    if (allocated(rmemory_dpotential_dxl)) &
      call ckpt_arr(do_write,rmemory_dpotential_dxl,size(rmemory_dpotential_dxl))
    if (allocated(rmemory_dpotential_dyl)) &
      call ckpt_arr(do_write,rmemory_dpotential_dyl,size(rmemory_dpotential_dyl))
    if (allocated(rmemory_dpotential_dzl)) &
      call ckpt_arr(do_write,rmemory_dpotential_dzl,size(rmemory_dpotential_dzl))
    if (allocated(PML_potential_acoustic_old)) &
      call ckpt_arr(do_write,PML_potential_acoustic_old,size(PML_potential_acoustic_old))
    if (allocated(PML_potential_acoustic_new)) &
      call ckpt_arr(do_write,PML_potential_acoustic_new,size(PML_potential_acoustic_new))
    if (allocated(rmemory_coupling_ac_el_displ)) &
      call ckpt_arr(do_write,rmemory_coupling_ac_el_displ,size(rmemory_coupling_ac_el_displ))
    if (allocated(rmemory_coupling_el_ac_potential)) &
      call ckpt_arr(do_write,rmemory_coupling_el_ac_potential,size(rmemory_coupling_el_ac_potential))
    if (allocated(rmemory_coupling_el_ac_potential_dot_dot)) &
      call ckpt_arr(do_write,rmemory_coupling_el_ac_potential_dot_dot,size(rmemory_coupling_el_ac_potential_dot_dot))
  endif

  ! partially filled seismogram buffers; dummy 1x1x1 when the component is not saved
  call ckpt_arr(do_write,seismograms_d,size(seismograms_d))
  call ckpt_arr(do_write,seismograms_v,size(seismograms_v))
  call ckpt_arr(do_write,seismograms_a,size(seismograms_a))
  call ckpt_arr(do_write,seismograms_p,size(seismograms_p))
  call ckpt_arr(do_write,seismograms_eps,size(seismograms_eps))
  if (ckpt_size_p > 0) call ckpt_arr(do_write,ckpt_seismo_p_gpu,ckpt_size_p)


  ! SIMULATION_TYPE 3: backward wavefields, kernel accumulators and the adjoint-source chunk
  if (SIMULATION_TYPE == 3) then
    if (allocated(source_adjoint)) &
      call ckpt_arr(do_write,source_adjoint,size(source_adjoint))
    if (allocated(b_displ)) &
      call ckpt_arr(do_write,b_displ,size(b_displ))
    if (allocated(b_veloc)) &
      call ckpt_arr(do_write,b_veloc,size(b_veloc))
    if (allocated(b_accel)) &
      call ckpt_arr(do_write,b_accel,size(b_accel))
    if (allocated(b_epsilondev_xx)) &
      call ckpt_arr(do_write,b_epsilondev_xx,size(b_epsilondev_xx))
    if (allocated(b_epsilondev_yy)) &
      call ckpt_arr(do_write,b_epsilondev_yy,size(b_epsilondev_yy))
    if (allocated(b_epsilondev_xy)) &
      call ckpt_arr(do_write,b_epsilondev_xy,size(b_epsilondev_xy))
    if (allocated(b_epsilondev_xz)) &
      call ckpt_arr(do_write,b_epsilondev_xz,size(b_epsilondev_xz))
    if (allocated(b_epsilondev_yz)) &
      call ckpt_arr(do_write,b_epsilondev_yz,size(b_epsilondev_yz))
    if (allocated(b_epsilondev_trace)) &
      call ckpt_arr(do_write,b_epsilondev_trace,size(b_epsilondev_trace))
    if (allocated(b_epsilon_trace_over_3)) &
      call ckpt_arr(do_write,b_epsilon_trace_over_3,size(b_epsilon_trace_over_3))
    if (allocated(b_R_trace)) &
      call ckpt_arr(do_write,b_R_trace,size(b_R_trace))
    if (allocated(b_R_xx)) &
      call ckpt_arr(do_write,b_R_xx,size(b_R_xx))
    if (allocated(b_R_yy)) &
      call ckpt_arr(do_write,b_R_yy,size(b_R_yy))
    if (allocated(b_R_xy)) &
      call ckpt_arr(do_write,b_R_xy,size(b_R_xy))
    if (allocated(b_R_xz)) &
      call ckpt_arr(do_write,b_R_xz,size(b_R_xz))
    if (allocated(b_R_yz)) &
      call ckpt_arr(do_write,b_R_yz,size(b_R_yz))
    if (allocated(b_displ_lddrk)) &
      call ckpt_arr(do_write,b_displ_lddrk,size(b_displ_lddrk))
    if (allocated(b_veloc_lddrk)) &
      call ckpt_arr(do_write,b_veloc_lddrk,size(b_veloc_lddrk))
    if (allocated(b_R_trace_lddrk)) &
      call ckpt_arr(do_write,b_R_trace_lddrk,size(b_R_trace_lddrk))
    if (allocated(b_R_xx_lddrk)) &
      call ckpt_arr(do_write,b_R_xx_lddrk,size(b_R_xx_lddrk))
    if (allocated(b_R_yy_lddrk)) &
      call ckpt_arr(do_write,b_R_yy_lddrk,size(b_R_yy_lddrk))
    if (allocated(b_R_xy_lddrk)) &
      call ckpt_arr(do_write,b_R_xy_lddrk,size(b_R_xy_lddrk))
    if (allocated(b_R_xz_lddrk)) &
      call ckpt_arr(do_write,b_R_xz_lddrk,size(b_R_xz_lddrk))
    if (allocated(b_R_yz_lddrk)) &
      call ckpt_arr(do_write,b_R_yz_lddrk,size(b_R_yz_lddrk))
    if (allocated(rho_kl)) &
      call ckpt_arr(do_write,rho_kl,size(rho_kl))
    if (allocated(mu_kl)) &
      call ckpt_arr(do_write,mu_kl,size(mu_kl))
    if (allocated(kappa_kl)) &
      call ckpt_arr(do_write,kappa_kl,size(kappa_kl))
    if (allocated(cijkl_kl)) &
      call ckpt_arr(do_write,cijkl_kl,size(cijkl_kl))
    if (allocated(hess_kl)) &
      call ckpt_arr(do_write,hess_kl,size(hess_kl))
    if (allocated(hess_rho_kl)) &
      call ckpt_arr(do_write,hess_rho_kl,size(hess_rho_kl))
    if (allocated(hess_mu_kl)) &
      call ckpt_arr(do_write,hess_mu_kl,size(hess_mu_kl))
    if (allocated(hess_kappa_kl)) &
      call ckpt_arr(do_write,hess_kappa_kl,size(hess_kappa_kl))
    if (allocated(b_potential_acoustic)) &
      call ckpt_arr(do_write,b_potential_acoustic,size(b_potential_acoustic))
    if (allocated(b_potential_dot_acoustic)) &
      call ckpt_arr(do_write,b_potential_dot_acoustic,size(b_potential_dot_acoustic))
    if (allocated(b_potential_dot_dot_acoustic)) &
      call ckpt_arr(do_write,b_potential_dot_dot_acoustic,size(b_potential_dot_dot_acoustic))
    if (allocated(rho_ac_kl)) &
      call ckpt_arr(do_write,rho_ac_kl,size(rho_ac_kl))
    if (allocated(kappa_ac_kl)) &
      call ckpt_arr(do_write,kappa_ac_kl,size(kappa_ac_kl))
    if (allocated(hess_ac_kl)) &
      call ckpt_arr(do_write,hess_ac_kl,size(hess_ac_kl))
    if (allocated(hess_rho_ac_kl)) &
      call ckpt_arr(do_write,hess_rho_ac_kl,size(hess_rho_ac_kl))
    if (allocated(hess_kappa_ac_kl)) &
      call ckpt_arr(do_write,hess_kappa_ac_kl,size(hess_kappa_ac_kl))
    if (allocated(b_displs_poroelastic)) &
      call ckpt_arr(do_write,b_displs_poroelastic,size(b_displs_poroelastic))
    if (allocated(b_velocs_poroelastic)) &
      call ckpt_arr(do_write,b_velocs_poroelastic,size(b_velocs_poroelastic))
    if (allocated(b_accels_poroelastic)) &
      call ckpt_arr(do_write,b_accels_poroelastic,size(b_accels_poroelastic))
    if (allocated(b_displw_poroelastic)) &
      call ckpt_arr(do_write,b_displw_poroelastic,size(b_displw_poroelastic))
    if (allocated(b_velocw_poroelastic)) &
      call ckpt_arr(do_write,b_velocw_poroelastic,size(b_velocw_poroelastic))
    if (allocated(b_accelw_poroelastic)) &
      call ckpt_arr(do_write,b_accelw_poroelastic,size(b_accelw_poroelastic))
    if (allocated(b_epsilonsdev_xx)) &
      call ckpt_arr(do_write,b_epsilonsdev_xx,size(b_epsilonsdev_xx))
    if (allocated(b_epsilonsdev_yy)) &
      call ckpt_arr(do_write,b_epsilonsdev_yy,size(b_epsilonsdev_yy))
    if (allocated(b_epsilonsdev_xy)) &
      call ckpt_arr(do_write,b_epsilonsdev_xy,size(b_epsilonsdev_xy))
    if (allocated(b_epsilonsdev_xz)) &
      call ckpt_arr(do_write,b_epsilonsdev_xz,size(b_epsilonsdev_xz))
    if (allocated(b_epsilonsdev_yz)) &
      call ckpt_arr(do_write,b_epsilonsdev_yz,size(b_epsilonsdev_yz))
    if (allocated(b_epsilonwdev_xx)) &
      call ckpt_arr(do_write,b_epsilonwdev_xx,size(b_epsilonwdev_xx))
    if (allocated(b_epsilonwdev_yy)) &
      call ckpt_arr(do_write,b_epsilonwdev_yy,size(b_epsilonwdev_yy))
    if (allocated(b_epsilonwdev_xy)) &
      call ckpt_arr(do_write,b_epsilonwdev_xy,size(b_epsilonwdev_xy))
    if (allocated(b_epsilonwdev_xz)) &
      call ckpt_arr(do_write,b_epsilonwdev_xz,size(b_epsilonwdev_xz))
    if (allocated(b_epsilonwdev_yz)) &
      call ckpt_arr(do_write,b_epsilonwdev_yz,size(b_epsilonwdev_yz))
    if (allocated(b_epsilons_trace_over_3)) &
      call ckpt_arr(do_write,b_epsilons_trace_over_3,size(b_epsilons_trace_over_3))
    if (allocated(b_epsilonw_trace_over_3)) &
      call ckpt_arr(do_write,b_epsilonw_trace_over_3,size(b_epsilonw_trace_over_3))
    if (allocated(rhot_kl)) &
      call ckpt_arr(do_write,rhot_kl,size(rhot_kl))
    if (allocated(rhof_kl)) &
      call ckpt_arr(do_write,rhof_kl,size(rhof_kl))
    if (allocated(sm_kl)) &
      call ckpt_arr(do_write,sm_kl,size(sm_kl))
    if (allocated(eta_kl)) &
      call ckpt_arr(do_write,eta_kl,size(eta_kl))
    if (allocated(mufr_kl)) &
      call ckpt_arr(do_write,mufr_kl,size(mufr_kl))
    if (allocated(B_kl)) &
      call ckpt_arr(do_write,B_kl,size(B_kl))
    if (allocated(C_kl)) &
      call ckpt_arr(do_write,C_kl,size(C_kl))
    if (allocated(M_kl)) &
      call ckpt_arr(do_write,M_kl,size(M_kl))
  endif

  ! gravity time series are held in memory for the whole run
  if (GRAVITY_SIMULATION) then
    if (do_write) then
      call gravity_checkpoint_save(IOUT_CKPT)
    else
      call gravity_checkpoint_restore(IOUT_CKPT)
    endif
  endif

  end subroutine checkpoint_body

!
!-------------------------------------------------------------------------------------------------
!

  subroutine ckpt_arr(do_write,a,n)

! external on purpose: sequence association lets one routine handle every rank

  use constants, only: CUSTOM_REAL
  use checkpoint_par, only: IOUT_CKPT

  implicit none

  logical,intent(in) :: do_write
  integer,intent(in) :: n
  real(kind=CUSTOM_REAL),intent(inout) :: a(n)

  if (do_write) then
    write(IOUT_CKPT) a
  else
    read(IOUT_CKPT) a
  endif

  end subroutine ckpt_arr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine checkpoint_cleanup()

! removes both slots after a completed run, so a rerun in the same directory starts from scratch

  use constants, only: MAX_STRING_LEN
  use specfem_par, only: NTSTEP_BETWEEN_CHECKPOINTS
  use checkpoint_par

  implicit none

  integer :: islot,ier
  character(len=MAX_STRING_LEN) :: fname
  logical :: exists

  if (NTSTEP_BETWEEN_CHECKPOINTS <= 0) return

  do islot = 0,1
    call checkpoint_okname(islot,fname)
    inquire(file=trim(fname),exist=exists)
    if (exists) then
      open(unit=IOUT_CKPT,file=trim(fname),status='old',iostat=ier)
      if (ier == 0) close(IOUT_CKPT,status='delete')
    endif
    call checkpoint_filename(islot,fname)
    inquire(file=trim(fname),exist=exists)
    if (exists) then
      open(unit=IOUT_CKPT,file=trim(fname),status='old',iostat=ier)
      if (ier == 0) close(IOUT_CKPT,status='delete')
    endif
  enddo

  end subroutine checkpoint_cleanup
