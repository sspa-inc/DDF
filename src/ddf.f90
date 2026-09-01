program ddf
  use m_river,   only :   nriv, xriv, yriv
  use m_grid,    only :   lfac, ncell, xcell, ycell, make_contour, hasNodeName, nodenames, create_grid
  use m_time,    only :   ntime,times
  use m_well,    only :   nwel, well, calc_grid_parameters, onereach
  use m_io,      only :   readdata,configfile,outfile,write_vals,write_vals_rowname
  use m_spatial, only :   pdist
  implicit none
  ! This program estimates the drawdown caused by a pumping well adjacent to a stream network.
  ! the weights of imaginary wells are determined by a linear system composing zero head change at all stream nodes: sum(weights * dd_riv_image) = dd_riv_pump


  ! local
  integer                    :: iwel, itime
  real, allocatable          :: leak(:), leaknet(:,:)               ! positive from stream to aquifer; negative from aquifer to stream
  real, allocatable          :: ddoff(:), dd(:), ddnet(:,:)         ! postive is head drawdown, negative is head raise

  character*999              :: str
  real                       :: dtoff, dt
  logical                    :: pumpoff


  call get_command_argument(1, configfile)
  call get_command_argument(2, str)
  onereach = trim(str)=="1"

  call readdata()
  if (lfac>0.0) then
    call calc_grid_parameters()
    call create_grid()
  end if

  allocate(ddnet(ncell, ntime), ddoff(ncell), dd(ncell))
  if (nriv>0) then
    allocate(leaknet(nriv, ntime),leak(nriv))
    leaknet = 0.0
    leak = 0.0
  else
    allocate(leaknet(1, ntime),leak(1))
  end if
  ddnet = 0.0

  do iwel=1, nwel
    if (abs(well(iwel)%q)<1e-10) cycle
    if(nriv>0) call well(iwel)%filter_riv()
    call well(iwel)%set_dist()
    if (nriv>0) call well(iwel)%filter_node()
  end do

  do itime = 1, ntime
    dt = times(itime)
    print*, "Calculating drawdown at Time ", dt
    do iwel=1, nwel
      if (abs(well(iwel)%q)<1e-10) cycle
      dtoff = 0.0
      if (well(iwel)%timeoff>0 .and. dt>well(iwel)%timeoff)  dtoff = dt - well(iwel)%timeoff
      ! calculate the constant pumping drawdown
      call well(iwel)%calc_dd(times(itime), dd, leak)
      if (well(iwel)%nriv>0) leaknet(:, itime) = leaknet(:, itime) + leak
      ! calculate if pumping is off
      if(dtoff>0) then
        call well(iwel)%calc_dd(dtoff, ddoff, leak)
        dd = dd - ddoff
        if (well(iwel)%nriv>0) leaknet(:, itime) = leaknet(:, itime) - leak*dtoff/dt
      end if
      ! else
      !    call well(iwel)%calc_dd(times(itime), dd)
      ! end if
      ddnet(:, itime) = ddnet(:, itime) + dd
    end do
    if (lfac>0.0) call make_contour(ddnet(:,itime), trim(outfile)//"_cc_"//real2str(times(itime)))
  end do
  ! print*, xcell(3281),ycell(3281)
  ! print*, minloc(r_wel_cel), active_cel(minloc(r_wel_cel))
  ! print*, "ddwel", ddwel(1741:1745)
  ! print*, "ddimg", ddimg(1741:1745)
  ! print*, "dd   ", dd(1742)
  ! print*, "idx  ", active_cel(1742)
  ! print*, "ddnet", ddnet(3200,:)
  ! print*, "theis", theis(rimg, times)
  if (hasNodeName) then
    call write_vals_rowname(trim(outfile)//"_dd.csv", ncell, ntime, xcell, ycell, "x,y"//strtime(ntime), ddnet, nodenames)
  else
    call write_vals(trim(outfile)//"_dd.csv", ncell, ntime, xcell, ycell, "x,y"//strtime(ntime), ddnet)
  end if
  if (nriv>0) call write_vals(trim(outfile)//"_dl.csv", nriv , ntime, xriv , yriv , "x,y"//strtime(ntime), leaknet)

  print*, "DDF Finished!"
  contains

  !> @brief Format the first nt values of the global times array as a comma-separated string.
  !>
  !> @param[in] nt Number of time values to include (1..ntime).
  !> @return Character string containing the formatted times.
  function strtime(nt)
  character*2028  :: strtime
  integer         :: nt
  write(strtime, "(*(:',',G0))") times(:nt)
  end function

  !> @brief Convert a real to a trimmed string with three decimal places.
  !>
  !> @param[in] r Real value to format.
  !> @return Character(20) string containing r formatted with F0.3 and left-adjusted.
  function real2str(r)
  character*20    :: real2str
  real            :: r
  write(real2str, "(F0.3)") r
  real2str = adjustl(real2str)
  end function

end program
