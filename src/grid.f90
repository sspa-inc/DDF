module m_grid
  implicit none
  ! cell location
  integer              :: ncell, nx, ny
  real, allocatable    :: xcell(:)          ! grid x coordinates
  real, allocatable    :: ycell(:)          ! grid y coordinates
  real, allocatable    :: xx(:), yy(:)      ! coordinates saved for contouring
  real                 :: cellsize1, cellmutl, contourlevel, lfac, gridrmax
  real                 :: xcenter, ycenter
  logical              :: hasNodeName = .false.
  character(len=64), allocatable    :: nodenames(:)      ! cell names

  contains



  !> @brief Create a structured, symmetric grid centered on the wells.
  !>
  !> Builds a square grid of cells around the wells using geometric spacing
  !> controlled by `cellsize1` and `cellmutl`. The grid extent is chosen based
  !> on maximum Theis influence and `gridrmax`.
  !>
  !> @note Populates module-level arrays `xx`, `yy`, `xcell`, `ycell`, `ncell`, `nx`, `ny`.
  subroutine create_grid()
  !
  integer               :: nx2, ii, ij
  real                  :: xtmp, ytmp, csize

  if (cellmutl - 1.0 < 1e-3) then
    cellmutl = 1.0
    nx2 = int(gridrmax / cellsize1) + 1
  else
    nx2 = int(log((gridrmax * (cellmutl - 1.0) / cellsize1 + 1)) / log(cellmutl))
  end if
  nx = nx2 * 2 + 1
  ny = nx
  ! print*, rmax, nx, nx2, nx2+2
  allocate(xx(-nx2:nx2), yy(-nx2:nx2))
  allocate(xcell(nx*ny), ycell(nx*ny))

  xx(0) = 0.0
  csize = 1.0
  do ii = 1, nx2
    xx(ii) = xx(ii-1) + csize
    csize = csize * cellmutl
  end do
  xx(-nx2:-1) = - xx(nx2:1:-1)
  yy = xx

  xx = xcenter + xx * cellsize1 / abs(lfac)
  yy = ycenter + yy * cellsize1 / abs(lfac)
  ncell = 0
  do ii=-nx2, nx2
    ytmp = yy(ii)
    cloop: do ij=-nx2, nx2
      ncell = ncell + 1
      xcell(ncell) = xx(ij)
      ycell(ncell) = ytmp
    end do cloop
  end do
  print*, "Grid created!"
  end subroutine create_grid


  !> @brief Create contour lines (GeoJSON) from an nx-by-ny grid of values.
  !>
  !> Writes contours to <outfile>.geojson using the global `contourlevel` and
  !> the grid arrays `xx` and `yy` for coordinates.
  !>
  !> @param[in] dd 2D array of values (nx,ny) to contour.
  !> @param[in] outfile Base filename (no extension) for the output GeoJSON.
  subroutine make_contour(dd, outfile)
  real              :: dd(nx, ny)
  character(*)      :: outfile

  integer           :: ifile, ii, imin, imax, nl
  real              :: dmin, dmax
  real, allocatable :: levels(:)
  if (contourlevel<=0.0) return
  dmin = minval(dd)
  imin = int(dmin/contourlevel)
  dmax = maxval(dd)
  imax = int(dmax/contourlevel)
  if (dmin>=0.0) imin = imin + 1
  if (dmax< 0.0) imax = imax - 1
  open(newunit=ifile, file=trim(outfile)//".geojson", status="replace")

  levels = [(ii,ii=imin,imax)] * contourlevel
  nl = imax-imin+1
  print*, "Writing", nl, "contour levels to ", trim(outfile)//".geojson"
  call conrec(dd,1,nx,1,ny,xx,yy,nl,levels,ifile)

  end subroutine

end module

