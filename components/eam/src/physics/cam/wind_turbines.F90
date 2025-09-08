
module wt_drag

!--------------------------------------------------------------------------
! Wind Turbine Drag parameterizations
!--------------------------------------------------------------------------
! Purpose:
!
! Module to compute the atmospheric forcing due to wind turbines
!
! Author: Chris Chartrand
!
!--------------------------------------------------------------------------

  use shr_kind_mod,  only: r8 => shr_kind_r8
  use ppgrid,        only: pcols, pver, pverp
  use constituents,  only: pcnst
  use physics_types, only: physics_state, physics_ptend, physics_ptend_init
  use spmd_utils,    only: iam, npes, masterproc
  use cam_history,   only: outfld
  use cam_logfile,   only: iulog
  use cam_abortutils,only: endrun

  use physconst,     only: pi, rearth

! Typical module header
  implicit none
  private
  save

!
! PUBLIC: interfaces
!
  public :: wt_read                  ! Read turbine parameters
  public :: wt_init                  ! Initialization
  public :: wt_tend                  ! Forcing tendency

!
! PRIVATE: Rest of the data and interfaces are private to this module
!

  logical,public :: l_wt_drag
  logical :: l_wt_debug

  !! turbine characteristics
  integer ntypes, nturbs
  integer,allocatable :: ncount(:)
  real*8,allocatable :: height(:), radius(:)
  real*8,allocatable :: wtlat(:), wtlon(:)
  integer,allocatable :: turbtype(:)

  !! power/thrust/Ctke curve table values
  integer :: tottablelines
  integer,allocatable :: ntablelines(:)
  integer,allocatable :: istart(:), istop(:)
  real*8,allocatable :: wtvel(:),wtct(:),wtcp(:),wtctke(:)

  integer, allocatable :: wtowner(:), wtcol(:), wtchnk(:)
  integer, allocatable :: wtnper(:,:)

!==========================================================================
contains
!==========================================================================

subroutine wt_init()

use cam_history,      only: addfld, add_default, horiz_only
use phys_grid,        only: get_ncols_p, phys_grid_initialized 

#if ( defined SPMD )
  use mpishorthand
#endif

    integer :: itype 
    character(len=16) :: suffix

    !! TODO remove this when done debugging
    l_wt_debug = .False.
    !l_wt_debug = .True.

    if (masterproc) then
      write(iulog,*) 'Performing wind turbine initialization'
      if (l_wt_debug) write(*,*) 'DEBUG preread l_wt_drag = ',l_wt_drag
      call wt_read()
      if (l_wt_debug) write(*,*) 'DEBUG postread l_wt_drag = ',l_wt_drag
    endif

    ! if read found turbine data, broadcast it
#ifdef SPMD
    call mpibcast (l_wt_drag, 1, mpilog, 0, mpicom)
#endif
    if (.not. l_wt_drag) then
      return
    endif

    call wt_bcast()

    ! Total velocity tendency output.
    call addfld ('UTWT',(/ 'lev' /), 'A','m/s2', &
         'u tendency - wind turbine drag')
    call addfld ('VTWT',(/ 'lev' /), 'A','m/s2', &
         'v tendency - wind turbine drag')
    !! turbulence tendency
    call addfld ('DTKEWT',(/ 'lev' /), 'A','m2/s3', &
         'TKE tendency - wind turbine')

    call add_default('UTWT', 1, ' ')
    call add_default('VTWT', 1, ' ')
    call add_default('DTKEWT ', 1, ' ')


    do itype = 1,ntypes
        write(suffix,'(i0.0)') itype

        !! total power produced
        call addfld ('PTWT'//trim(suffix),horiz_only, 'A','W', &   !! 1/2 Cp A U^3 * Nt
             'Wind Turbine Power absorption total')

        !! instantanious fields for debugging
        call addfld ('PowWT'//trim(suffix),horiz_only, 'I','W', &   !! 1/2 Cp A U^3 
             'Wind Turbine Power absorption per')
        call addfld ('FXWT'//trim(suffix),horiz_only, 'I','m4/s2', &   !! turbine force
             'x Force/density - wind turbine drag')
        call addfld ('FYWT'//trim(suffix),horiz_only, 'I','m4/s2', &
             'y Force/density - wind turbine drag')
        call addfld ('KTWT'//trim(suffix),horiz_only, 'I','m2/s3', &    !! tke source
             'tke - wind turbine production')

        call addfld ('CTWT'//trim(suffix),horiz_only, 'I','-', &
             'Thrust Coefficient wind turbine drag')
        call addfld ('CPWT'//trim(suffix),horiz_only, 'I','-', &
             'Power Coefficient wind turbine drag')
        call addfld ('CTKEWT'//trim(suffix),horiz_only, 'I','-', &
             'TKE Coefficient wind turbine drag')

        call addfld ('AreaWT'//trim(suffix),horiz_only, 'I','m2', &
             'Disk Area wind turbine drag')
        call addfld ('WTUD'//trim(suffix),horiz_only, 'I','m/s', &
             'WT Disk Average Velocity')
        call addfld ('WTVD'//trim(suffix),horiz_only, 'I','m/s', &
             'WT Disk Average Velocity')

        call add_default('PTWT'//trim(suffix), 1, ' ')
        call add_default('PowWT'//trim(suffix), 1, ' ')
        call add_default('FXWT'//trim(suffix), 1, ' ')
        call add_default('FYWT'//trim(suffix), 1, ' ')
        call add_default('KTWT'//trim(suffix), 1, ' ')

        call add_default('CTWT'//trim(suffix), 1, ' ')
        call add_default('CPWT'//trim(suffix), 1, ' ')
        call add_default('CTKEWT'//trim(suffix), 1, ' ')

        call add_default('AreaWT'//trim(suffix), 1, ' ')
        call add_default('WTUD'//trim(suffix), 1, ' ')
        call add_default('WTVD'//trim(suffix), 1, ' ')
    enddo


    !! get grid columns for turbines
    if ( phys_grid_initialized() ) then
      call wt_getcols()
    else
      call endrun( "Wind Turbine error: phys_grid not called yet" )
    endif

end subroutine wt_init

!!******************************************************************************
  !-----------------------------------------------------------------------
  ! Time independent initialization for wind turbine parameters
  !-----------------------------------------------------------------------

subroutine wt_read()

  use ioFileMod,        only: getfil
  use units,           only: getunit, freeunit

#if ( defined SPMD )
  use mpishorthand
#endif

  ! VARIABLES NEEDED TO READ IN TABLE OF SPECTRA FROM FILE

  ! Local variables
  integer :: unitn, ierr
  integer, allocatable :: tableunitn(:)
  real*8 :: latdum, londum, diam, drag0, nompow, ctke0
  integer :: itype, iturb, ictkemodel
  integer :: i,j,k,ii,nn,il,iu,idx

  character(len=*), parameter :: subname = 'wt_read'

  character(len=256) :: msgstring

  !! io files
  character(len=256) :: wt_loc_file = 'windturbines.txt'
  character(len=256) :: wt_drag_file ! drag coefficient file
  character(len=256) :: wt_drag_file_loc ! local filepath of wt_drag_file
  character(len=256) :: wt_loc_file_loc ! local filepath of wt_loc_file
  character(len=256) :: stringdum

  logical :: lexists
  !----------------------------------------------------------------------
  ! read in turbine parameters
  !-----------------------------------------------------------------------

    !! TODO put this in the namelist file
    l_wt_drag=.True.

    ! search for turbine location file, returns as 'wt_loc_file_loc'
    call getfil(wt_loc_file, wt_loc_file_loc, 1, lexists)

    if ( .not. lexists) then
        write(iulog,*) 'No wind turbine locations to read'
        l_wt_drag=.False.
        return
    endif

    ! if it exists, read it
    unitn = getunit()
    open( unitn, file=trim(wt_loc_file_loc), status='old', iostat=ierr )
    if (ierr /= 0) then
        call endrun(subname // ':: ERROR opening turbine location list')
        l_wt_drag=.False.
        return
    endif

    ! count the number of lines (location entries)
    nturbs=0
    ntypes=0
    do while (.true.)
        !! TODO enforce formatting check real, real, int
        read(unitn, *, iostat=ierr) latdum, londum, itype
        if (ierr /= 0) then
           exit
        endif
        nturbs = nturbs+1
        ntypes = max(ntypes,itype)
    enddo

    !! allocate arrays
    allocate(ncount(ntypes))
    allocate(ntablelines(ntypes))
    allocate(height(ntypes))
    allocate(radius(ntypes))
    allocate(istart(ntypes))
    allocate(istop(ntypes))

    allocate(wtlat(nturbs))
    allocate(wtlon(nturbs))
    allocate(turbtype(nturbs))

    allocate(tableunitn(ntypes))

    write(iulog,*) nturbs,'total turbines defined'
    write(iulog,*) ntypes,'types of turbines defined'

    !! now populate the location arrays
    ncount(:) = 0
    rewind(unitn)
    do iturb=1,nturbs
        read(unitn, *, iostat=ierr) wtlat(iturb), wtlon(iturb), turbtype(iturb)
        if (ierr /= 0) then
           write(msgstring,*) ':: ERROR ',ierr,' while reading turbine locations at line number: ',iturb,' of ',nturbs
           call endrun(subname // msgstring)
        end if
        if (wtlon(iturb) < 0.0) then
            wtlon(iturb) = wtlon(iturb) + 360.0
        endif      
        itype = turbtype(iturb)
        ncount(itype) = ncount(itype) +1
    enddo
    close(unitn)
    call freeunit(unitn)

    write(iulog,*) 'Finished reading wind turbine locations'
    do itype=1,ntypes
        if (ncount(itype)==0) cycle
        write(iulog,*) ncount(itype),'turbines of type',itype
    enddo

    ! read in wtct, wtcp, wtctke tables from wind-turbine-XX.tbl
    tottablelines=0
    ntablelines(:) = 0
    do itype=1,ntypes
        if (ncount(itype) == 0) then  !! skip any types which weren't defined
            ntablelines(itype)=0
            istart(itype) = tottablelines
            istop(itype) = tottablelines
        else
            unitn = getunit()
            tableunitn(itype) = unitn

            ! search for turbine drag file, returns local path as 'wt_drag_file_loc'
            write(wt_drag_file,'(a,i0.0,a)') 'wind-turbine-',itype,'.tbl'
            call getfil(wt_drag_file, wt_drag_file_loc, 1, lexists)

            open(unitn, file=trim(wt_drag_file_loc), status='old', iostat=ierr )

            read(unitn, *, iostat=ierr)  ntablelines(itype)

            if (ierr /= 0) then
                write(*,*) ':: ERROR reading turbine drag table size'
                call exit(1)
            endif
            write(*,*) ntablelines(itype)
            ntablelines(itype) = ntablelines(itype) + 1  !! add one for the 0th entry

            istart(itype) = tottablelines +1   !! first index in coefficient tables for each type
            tottablelines = tottablelines + ntablelines(itype)
            istop(itype) = tottablelines  !! last index of each table
        endif
    enddo

    ! allocate curve arrays
    allocate(wtvel(tottablelines))
    allocate(wtct(tottablelines))
    allocate(wtcp(tottablelines))
    allocate(wtctke(tottablelines))

    ! populate curve arrays
    do itype=1,ntypes
        if (ncount(itype) == 0) cycle !! skip any types which weren't defined

        unitn=tableunitn(itype)

        il = istart(itype)
        iu = istop(itype)

        !! check number of entries to see if Ctke is specified
        read(unitn, '(a)')  stringDum  !! read first line
        read(stringDum,*,iostat=ierr) height(itype), diam, drag0, nompow, ctke0  !! check for 5 floats in the line
        if (ierr == 0) then  !! found Ctke
            !write(*,*)  ':: Model type 1'
            ictkeModel=1
        else  !! read error means Ctke was NOT specified
            read(stringDum,*,iostat=ierr) height(itype), diam, drag0, nompow  !! try converting just 4 floats
            if (ierr /= 0) then  !! still errors, then bad format
                write(msgstring,*)  ':: ERROR reading turbine drag table data', ierr
                call endrun(subname // msgstring)
            endif
            !write(*,*)  ':: Model type 2'
            ictkeModel=2   !! fitch model, Ctke = Ct-Cp
            ctke0 = drag0 !  Ctke0 = Ct0-Cp0
        endif      
        radius(itype) = diam/2.0

        !! set 1st value of table arrays
        wtvel(il) = 0.0
        wtcp(il)  = 0.0
        wtct(il)  = drag0
        wtctke(il) = ctke0

        wtct(il) = max(wtct(il),0.0)  !! just to be safe
        wtctke(il) = max(wtctke(il),0.0)  !! just to be safe
        
        !! read in the rest of the table
        do ii = 1, ntablelines(itype)-1
           idx = il+ii
           if ( ictkeModel == 1) then   !! defined ctke
               read(unitn, *, iostat=ierr)  wtvel(idx), wtcp(idx), wtct(idx), wtctke(idx)
           else if ( ictkeModel == 2) then   !! Fitch model
               read(unitn, *, iostat=ierr)  wtvel(idx), wtcp(idx), wtct(idx)
               wtctke(idx) = wtct(idx) - wtcp(idx)
           endif
           wtct(idx) = max(wtct(idx),0.0)  !! just to be safe
           wtcp(idx) = max(wtcp(idx),0.0)  !! just to be safe
           wtctke(idx) = max(wtctke(idx),0.0)  !! just to be safe
           if (ierr /= 0) then
             write(msgstring,*)  ':: ERROR ',ierr,' while reading coef table type',itype,'at line number: ',ii+1
                call endrun(subname // msgstring)
           end if
        enddo

        write(iulog,*) 'Finished reading wind turbine data', itype, ntypes 
        close(unitn)
        call freeunit(unitn)
    enddo  !! ntypes

#ifdef SPMD
#endif


end subroutine wt_read

!!******************************************************************************

subroutine wt_bcast()

#if ( defined SPMD )
  use mpishorthand
#endif

#ifdef SPMD
    call mpibcast (nturbs,        1,         mpiint, 0, mpicom)
    call mpibcast (ntypes,        1,         mpiint, 0, mpicom)

    call mpibcast (tottablelines, 1,         mpiint, 0, mpicom)
#endif
    ! already allocated on master
    if (.not. masterproc) then
        allocate(ncount(ntypes))
        allocate(ntablelines(ntypes))
        allocate(height(ntypes))
        allocate(radius(ntypes))
        allocate(istart(ntypes))
        allocate(istop(ntypes))

        allocate(wtlat(nturbs))
        allocate(wtlon(nturbs))
        allocate(turbtype(nturbs))

        allocate(wtvel(tottablelines))
        allocate(wtct(tottablelines))
        allocate(wtcp(tottablelines))
        allocate(wtctke(tottablelines))
    endif
#ifdef SPMD
    call mpibcast (ncount,       ntypes,  mpiint, 0, mpicom)
    call mpibcast (ntablelines,  ntypes,  mpiint, 0, mpicom)

    call mpibcast (height,   ntypes,    mpir8,  0, mpicom)
    call mpibcast (radius,   ntypes,    mpir8,  0, mpicom)
    call mpibcast (istart,   ntypes,    mpiint, 0, mpicom)
    call mpibcast (istop,    ntypes,    mpiint, 0, mpicom)

    call mpibcast (wtlat,    nturbs,    mpir8,  0, mpicom)
    call mpibcast (wtlon,    nturbs,    mpir8,  0, mpicom)
    call mpibcast (turbtype, nturbs,    mpiint, 0, mpicom)

    call mpibcast (wtvel,    tottablelines, mpir8,  0, mpicom)
    call mpibcast (wtct,     tottablelines, mpir8,  0, mpicom)
    call mpibcast (wtcp,     tottablelines, mpir8,  0, mpicom)
    call mpibcast (wtctke,   tottablelines, mpir8,  0, mpicom)
#endif

end subroutine wt_bcast


!!******************************************************************************
subroutine wt_crossec( icol, itype, nlines,                    & 
                       u, v, zi,                               &
                       wtvel, wtct, wtcp, wtctke,              &
                       Cd, Cp, Ctke, udisk, vdisk, areatot,    &
                       UUA, UVA, U3A                           &
                     )

    use ppgrid,        only:  pver, pverp

    integer :: icol,itype,nlines
    real*8 :: u(pver), v(pver), zi(pverp)
    real*8 :: wtvel(nlines), wtct(nlines), wtcp(nlines), wtctke(nlines)
    real*8 :: cd, cp, ctke, udisk, vdisk, areatot
    real*8 :: UUA(pver), UVA(pver), U3A(pver)

    integer :: i,j,k
    real*8 :: ttop, tbot, R2, d_R, pi
    real*8 :: area(pverp)
    real*8 :: wtblockage, umagloc, umagdisk, interp

    !! top and bottom of the turbine disc
    tbot = height(itype)-radius(itype)
    ttop = height(itype)+radius(itype)

    ! R^2
    R2 = radius(itype)**2
    pi = acos(-1.0)

    !! layer interface loop, get a_k:
    !!     area of turbine cross-section above layer interface k
    do k=1,pverp
        !! zi goes from upper to lower
        if (zi(k) <= tbot ) then
          area(k) = pi*R2
        elseif (zi(k) >= ttop ) then
          area(k) = 0.0
        else
          d_R = (zi(k)-height(itype))/radius(itype)
          !! area of cross section above layer k
          area(k) = R2*( acos(d_R) - d_R*sqrt(1.0-d_R**2) )
        endif
    enddo

    !! get disk velocity
    udisk=0.0
    vdisk=0.0
    areatot=0.0
    do k=1,pver
        !! cross sectional area within layer is difference of
        !!   areas above upper and lower boundary
        !! This is Aijk in Fitch 2012
        wtblockage = area(k+1)-area(k)
        areatot = areatot + wtblockage

        !! get area averaged disk velocity
        udisk = udisk + u(k)*wtblockage
        vdisk = vdisk + v(k)*wtblockage

        !! save U*A for each face
        umagLoc = sqrt(u(k)**2+v(k)**2)

        !! this is |U|U*A in Fitch
        UUA(k) = umagLoc*u(k)*wtblockage
        UVA(k) = umagLoc*v(k)*wtblockage

        U3A(k) = umagLoc**3*wtblockage
    enddo
    !! are averaged disk velocity
    udisk = udisk / areatot
    vdisk = vdisk / areatot

    !! magnitude of disk velocity
    umagDisk = sqrt(udisk**2 + vdisk**2)

    !! use average disk velocity to get Cd, Cp
    if ( umagDisk > wtvel(nlines) ) then
        Cd = wtct(nlines)
        Cp = wtcp(nlines)
    else if ( umagDisk < wtvel(1) ) then
        Cd = wtct(1)
        Cp = wtcp(1)
    else
        do i=2,nlines
            !! at the first value higher than umagDisk, interpolate, break the loop
            if (wtvel(i) > umagDisk ) then
                !! linear interpolant between table velocity values
                interp = (umagDisk - wtvel(i-1) ) / ( wtvel(i) - wtvel(i-1) )
                Cd   = interp * (wtct(i)-wtct(i-1)) + wtct(i-1)
                Cp   = interp * (wtcp(i)-wtcp(i-1)) + wtcp(i-1)
                Ctke = interp * (wtctke(i)-wtctke(i-1)) + wtctke(i-1)
                exit
            endif
        enddo
    endif

end subroutine wt_crossec

!!******************************************************************************
!!  subroutine to locate column of each turbine
!!******************************************************************************
subroutine wt_getcols()
    use phys_grid, only: phys_grid_find_col

#if ( defined SPMD )
  use mpishorthand
#endif

    integer :: i,j,k
    real*8, allocatable :: rlat(:), rlon(:), idyn_dist(:)

    allocate(wtowner(nturbs))     ! rank of chunk owner
    allocate(wtchnk(nturbs))      ! local chunk index
    allocate(wtcol(nturbs))       ! column index within the chunk

    if (l_wt_debug) then
    if (masterproc) then
        write(iulog,*) 'Turbine DEBUG', nturbs
    endif
    endif

    do i=1,nturbs
        call phys_grid_find_col( wtlat(i), wtlon(i), wtowner(i), wtchnk(i), wtcol(i) )
!        write(iulog,*) i, wtowner(i), wtcol(i), wtchnk(i)
    enddo

    if (l_wt_debug) then
    if (masterproc) then
        write(iulog,*) 'End Turbine DEBUG'
        write(iulog,*) 'Finished binning turbines'
    endif
    endif

end subroutine wt_getcols

!!******************************************************************************

subroutine wt_tend(state, sgh, pbuf, dt, ptend, cam_in)
  !-----------------------------------------------------------------------
  ! Interface for wind turbine drag parameterization.
  !-----------------------------------------------------------------------
  use physics_types,  only: physics_state_copy !, set_dry_to_wet
  use physics_buffer, only: physics_buffer_desc, pbuf_get_index, pbuf_get_field, pbuf_set_field
  use camsrfexch, only: cam_in_t
  use phys_grid,  only: get_area_p
  use phys_grid, only: get_gcol_p
  !------------------------------Arguments--------------------------------
  type(physics_state), intent(in) :: state      ! physics state structure
  ! Standard deviation of orography.
  real(r8), intent(in) :: sgh(pcols)
  type(physics_buffer_desc), pointer :: pbuf(:) ! Physics buffer
  real(r8), intent(in) :: dt                    ! time step
  ! Parameterization net tendencies.
  type(physics_ptend), intent(out):: ptend
  type(cam_in_t), intent(in) :: cam_in

  type(physics_state) :: state1     ! Local copy of state variable

  integer :: lchnk                  ! chunk identifier
  integer :: ncol                   ! number of atmospheric columns

  integer :: i, k                   ! loop indices
  integer :: icol, iturb, itype, gcol

  !! variables for calculating thrust, etc.
  real(r8) :: Cd(state%ncol,ntypes)
  real(r8) :: Cp(state%ncol,ntypes)
  real(r8) :: Ctke(state%ncol,ntypes)
  real(r8) :: areaTot(state%ncol,ntypes)
  real(r8) :: UUA(state%ncol,ntypes,pver)
  real(r8) :: UVA(state%ncol,ntypes,pver)
  real(r8) :: U3A(state%ncol,ntypes,pver)
  real(r8) :: wtud(state%ncol,ntypes), wtvd(state%ncol,ntypes)
  real(r8) :: interp, vol, totvol, volinv

!  ! Which constituents are being affected by diffusion.
  logical  :: lq(pcnst)

  ! Contiguous copies of state arrays.
  real(r8) :: u(state%ncol,pver)
  real(r8) :: v(state%ncol,pver)
  real(r8) :: zm(state%ncol,pver)
  real(r8) :: zi(state%ncol,pver+1)

  real(r8) :: dTKEdt(state%ncol,pver)

  real(r8) :: wtpowTot(state%ncol,ntypes)

  real(r8) :: wtpower(state%ncol,ntypes)
  real(r8) :: forcex(state%ncol,ntypes)
  real(r8) :: forcey(state%ncol,ntypes)
  real(r8) :: tkesrc(state%ncol,ntypes)

  real(r8), pointer :: tke(:,:)
  real(r8), pointer :: up2(:,:)
  real(r8), pointer :: vp2(:,:)
  real(r8), pointer :: wp2(:,:)
  integer :: tke_idx, up2_idx, vp2_idx, wp2_idx

  real(r8) :: fxtmp, fytmp, tketmp, powtmp
  character(len=16) :: suffix

!  !------------------------------------------------------------------------
!
  !! TODO remove this, unnecessary
  ! Make local copy of input state.
  call physics_state_copy(state, state1)

  if (l_wt_debug) then
  if (masterproc) then
      write(iulog,*) 'DEBUG Entering Turbine Tendendcy'
  endif
  endif

  lchnk = state1%lchnk
  ncol  = state1%ncol

  !! probably don't need to copy state in the first place
  u = state1%u(:ncol,:)
  v = state1%v(:ncol,:)
  zm = state1%zm(:ncol,:)
  zi = state1%zi(:ncol,:)

  !! grab tke
  tke_idx  = pbuf_get_index('tke')
  up2_idx  = pbuf_get_index('UP2_nadv')
  vp2_idx  = pbuf_get_index('VP2_nadv')
  wp2_idx  = pbuf_get_index('WP2_nadv')
  call pbuf_get_field(pbuf, tke_idx,     tke)
  call pbuf_get_field(pbuf, up2_idx,     up2)
  call pbuf_get_field(pbuf, vp2_idx,     vp2)
  call pbuf_get_field(pbuf, wp2_idx,     wp2)

  lq = .true.
  call physics_ptend_init(ptend, state1%psetcols, "wind_turb_drag", &
       ls=.true., lu=.true., lv=.true., lq=lq)

  !! find intersected column layers
  if (l_wt_debug) then
  if (masterproc) then
      write(iulog,*) 'DEBUG Calculating Crossection Geometry'
  endif
  endif

  if ( .not. allocated(wtnper) ) then
      allocate(wtnper(ncol,ntypes))
  endif

  !! loop turbines
  if (l_wt_debug) then
  if (masterproc) then
      write(iulog,*) 'DEBUG Calculating turbine crossection'
  endif
  endif


  !! ncolumn, ntypes
  Cd(:,:) = 0.0
  Cp(:,:) = 0.0
  Ctke(:,:) = 0.0
  wtnper(:,:) = 0
  areatot(:,:) = 0.0


  !! loop turbines, bin to columns, and calculate coefficients
  do iturb = 1, nturbs
      if (iam==wtowner(iturb) .and.  lchnk==wtchnk(iturb))  then
          !! column the turbine is in
          icol=wtcol(iturb)
          itype=turbtype(iturb)
          !! number of turbines per column
          wtnper(icol,itype) = wtnper(icol,itype) + 1

          if ( wtnper(icol,itype) == 1) then  !! first turbine in column
              !! need to calculate values
              call wt_crossec( icol, itype, ntablelines(itype), &   !! column #, type #, length of coefficicent curve table
                               u( icol,: ), v( icol,: ),        &   !! pass in velocity at each layer in column
                               zi( icol,: ),                    &   !! layer interface heights
                               wtvel(  istart(itype): ),        &   !! velocity axis
                               wtct(   istart(itype): ),        &   !! Ct curve for this turbine type
                               wtcp(   istart(itype): ),        &   !! Ct curve for this turbine type
                               wtctke( istart(itype): ),        &   !! Ctke curve for this turbine type
                               !! start of return values
                               Cd(   icol,itype ),              &   !! Cd interpolated from curves using disk averaged velocity
                               Cp(   icol,itype ),              &
                               Ctke( icol,itype ),              &
                               wtud( icol,itype ),              &   !! disk velocity averaged over all intersected layers
                               wtvd( icol,itype ),              &
                               areatot( icol,itype ),           &   !! occluded area in this column by turbine type
                               UUA( icol, itype, :),            &   !! U*Umag*Area
                               UVA( icol, itype, :),            &   !! V*Umag*Area
                               U3A( icol, itype, :)             &   !! Umag^3 *Area
                        )
          endif
      endif
  enddo

  !! loop columns and turbine types, sum up forces and sources
  forcex(:,:)   = 0.0
  forcey(:,:)   = 0.0
  tkesrc(:,:)   = 0.0
  wtpower(:,:)  = 0.0
  wtpowTot(:,:) = 0.0

  ptend%u(:,:)  = 0.0
  ptend%v(:,:)  = 0.0
  dTKEdt(:,:)   = 0.0
  do icol = 1, ncol
      do itype=1,ntypes
          if ( wtnper(icol,itype) == 0 ) cycle  !! nothing to do here

          totvol=0.0
          do k = 1, pver
               vol = get_area_p(lchnk,icol)  !! rads^2
               vol = vol*rearth**2        !! meters**2
               vol = vol*(zi(icol,k)-zi(icol,k+1))  !! meters**3
               volinv = 1.0/vol              !! inverse
               totvol = totvol + vol
               if (l_wt_debug) then
                   write(*,*) 'DEBUG Volume Calc:', get_area_p(lchnk,icol),rearth,(zi(icol,k)-zi(icol,k+1)),vol
               endif

               !! force, source, and power production of one turbine at level k
               fxtmp  =  -1.0/2.0*Cd(icol,itype)  *UUA(icol,itype,k)   !! X-force of a single turbine in leyer k
               fytmp  =  -1.0/2.0*Cd(icol,itype)  *UVA(icol,itype,k)
               tketmp =   1.0/2.0*Ctke(icol,itype)*U3A(icol,itype,k)   !! turbulent source from a single turbine in layer k
               powtmp =   1.0/2.0*Cp(icol,itype)  *U3A(icol,itype,k)   !! power production from a single turbine in layer k

               !! forces and power are summed over level to get toal per turbine
               forcex(icol,itype)  = forcex(icol,itype)  + fxtmp       !! x-force of a single turbine, summed over all layers
               forcey(icol,itype)  = forcey(icol,itype)  + fytmp
               tkesrc(icol,itype)  = tkesrc(icol,itype)  + tketmp      !! tke source of a single turbine, summed over all layers
               wtpower(icol,itype) = wtpower(icol,itype) + powtmp      !! power of a single turbine, summed over all layers

               !! update tendency
               !! tendencies are summed over types to get total contribution to cell
               ptend%u(icol,k) = ptend%u(icol,k) +fxtmp*wtnper(icol,itype)*volinv
               ptend%v(icol,k) = ptend%v(icol,k) +fytmp*wtnper(icol,itype)*volinv
               dTKEdt(icol,k)  = dTKEdt(icol,k)  +tketmp*wtnper(icol,itype)*volinv
          end do
          !! turbine power is scaled by number/column to get total production
          wtpowtot(icol,itype)  = wtpower(icol,itype)*wtnper(icol,itype)  !! total power from all turbines in this column of this type
      end do
      !! after summing tendencies over all turbine types, update turbulence
      do k = 1, pver
           tke(icol,k) = tke(icol,k) + dTKEdt(icol,k)*dt
           up2(icol,k) = up2(icol,k) + 2.0/3.0*dTKEdt(icol,k)*dt
           vp2(icol,k) = vp2(icol,k) + 2.0/3.0*dTKEdt(icol,k)*dt
           wp2(icol,k) = wp2(icol,k) + 2.0/3.0*dTKEdt(icol,k)*dt
      enddo
  end do
  ! set new turbulence values
  call pbuf_set_field(pbuf, tke_idx,     tke)
  call pbuf_set_field(pbuf, up2_idx,     up2)
  call pbuf_set_field(pbuf, vp2_idx,     vp2)
  call pbuf_set_field(pbuf, wp2_idx,     wp2)

  ! Write output fields to history file
  call outfld('UTWT', ptend%u,  ncol, lchnk)   !! U-tendency @ level, column
  call outfld('VTWT', ptend%v,  ncol, lchnk)
  call outfld('DTKEWT', dTKEdt,  ncol, lchnk)  !! TKE tendency @ level, column

  do itype = 1,ntypes
      write(suffix,'(i0.0)') itype
      call outfld('PTWT'//trim(suffix),   wtpowTot(:,itype),  ncol, lchnk)  !! total power @ column by type

      call outfld('PowWT'//trim(suffix),  wtpower(:,itype),  ncol, lchnk)  !! wind turbine power @ column by type
      call outfld('FXWT'//trim(suffix),   forcex(:,itype),  ncol, lchnk)   !! force of EACH turbine of itype in column 
      call outfld('FYWT'//trim(suffix),   forcey(:,itype),  ncol, lchnk)
      call outfld('KTWT'//trim(suffix),   tkesrc(:,itype),  ncol, lchnk)

      call outfld('CTWT'//trim(suffix),   cd(:,itype),  ncol, lchnk)      !! coefficents used in caluclations, interpolated form disk velocity
      call outfld('CPWT'//trim(suffix),   cp(:,itype),  ncol, lchnk)
      call outfld('CTKEWT'//trim(suffix), ctke(:,itype),  ncol, lchnk)

      call outfld('AreaWT'//trim(suffix), areatot(:,itype),  ncol, lchnk) !! wind turbine disc area @ column by type
      call outfld('WTUD'//trim(suffix),   wtud(:,itype),  ncol, lchnk)    !! disk average velocity @ column by type
      call outfld('WTVD'//trim(suffix),   wtvd(:,itype),  ncol, lchnk)
  enddo
    
end subroutine wt_tend


!==========================================================================

end module wt_drag
