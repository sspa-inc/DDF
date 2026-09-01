module m_io
  integer                    :: ioerr
  character*999              :: configfile
  character*999              :: outfile
  contains


  !> @brief Replace characters from a set in a string with a target character.
  !>
  !> Replaces any character found in `charset` within `string` with
  !> `target_char`. If `charset` is omitted it defaults to comma and tab.
  !> If `target_char` is omitted it defaults to a single space.
  !>
  !> @param[in] string Input string.
  !> @param[in] charset Optional array of characters to replace.
  !> @param[in] target_char Optional replacement character.
  !> @return New string with replacements applied.
  pure function replace_line(string, charset, target_char) result(res)
  character(*), intent(in) :: string
  character, intent(in), optional :: charset(:), target_char
  character(len(string)) :: res
  character, allocatable :: real_charset(:)
  character              :: real_target_char
  integer :: i
  res = string
  if (present(charset)) then
    real_charset = charset
  else
    real_charset = [",", char(9)]
  end if
  if (present(target_char)) then
    real_target_char = target_char
  else
    real_target_char = " "
  end if
  do i = 1, len(string)
    if (any(string(i:i) == real_charset)) then
      res(i:i) = real_target_char
    end if
  end do
  end function replace_line

  !> @brief Read an unlimited length line from unit number lun into a deferred-
  !! length character string (line).
  !!
  !! Tack on a single space to the end so that routines like URWORD continue to
  !! function as before.
  !<
  subroutine get_line(lun, line, iostat)
  ! -- dummy
  integer, intent(in) :: lun
  character(len=:), intent(out), allocatable :: line
  integer, intent(out) :: iostat
  ! -- local
  integer, parameter :: buffer_len = 1024
  character(len=buffer_len) :: buffer
  character(len=:), allocatable :: linetemp
  integer :: size_read, linesize
  !
  ! -- initialize
  line = ''
  linetemp = ''
  !
  ! -- process
  do
    read (lun, '(A)', iostat=iostat, advance='no', size=size_read) buffer
    if (is_iostat_eor(iostat)) then
      linesize = len(line)
      deallocate (linetemp)
      allocate (character(len=linesize) :: linetemp)
      linetemp(:) = line(:)
      deallocate (line)
      allocate (character(len=linesize + size_read + 1) :: line)
      line(:) = linetemp(:)
      line(linesize + 1:) = buffer(:size_read)
      linesize = len(line)
      line(linesize:linesize) = ' '
      iostat = 0
      exit
    else if (iostat == 0) then
      linesize = len(line)
      deallocate (linetemp)
      allocate (character(len=linesize) :: linetemp)
      linetemp(:) = line(:)
      deallocate (line)
      allocate (character(len=linesize + size_read) :: line)
      line(:) = linetemp(:)
      line(linesize + 1:) = buffer(:size_read)
    else
      exit
    end if
  end do
  end subroutine get_line

  !> @brief Count the number of lines in a text file.
  !>
  !> Opens the file and counts non-empty read attempts to determine the
  !> number of lines. Prints a detection message to standard output.
  !>
  !> @param[in] file Path to the file.
  !> @return Number of lines detected in the file.
  function linecount(file)
  integer         :: linecount
  character(*)    :: file
  open(newunit=ifile, file=trim(file), status='old')
  linecount = 0
  do
    read(ifile,*, iostat=ioerr)
    if (ioerr/=0) EXIT
    linecount = linecount + 1
  end do
  close(ifile)
  print*, "Detect ", linecount, " lines in "//trim(file)
  end function

  !> @brief Detect number of columns in a delimited text file using the first line.
  !>
  !> Reads the first line using `get_line`, normalizes delimiters to spaces, then
  !> counts the number of whitespace-delimited tokens.
  !>
  !> @param[in] file Path to the file.
  !> @return Detected number of columns on the first data line.
  function columncount(file)
  integer         :: columncount
  character(*)    :: file
  character(len=:), allocatable :: line
  integer         :: i0

  open(newunit=ifile, file=trim(file), status='old')
  call get_line(ifile, line, ioerr)
  close(ifile)
  line = replace_line(adjustl(line))

  columncount = 1
  do i = 2, len_trim(line)
    if (line(i:i) /= " " .and. line(i-1:i-1) == " ") columncount = columncount + 1
  end do
  print*, "Detect ", columncount, " columns in "//trim(file)
  end function

  !> @brief Read configuration and data files (grid, wells, river, times).
  !>
  !> Reads settings from `configfile` (global), then reads the referenced
  !> grid, well, river and time files. Populates global module variables such
  !> as `well`, `xcell`, `xriv`, `times`, etc.
  !>
  !> @note This routine uses module globals and opens files based on paths read
  !> from the configuration file.
  subroutine readdata()
  use m_grid
  use m_river
  use m_time
  use m_well, only            : nwel, well, onereach, ddtor, dptor, rimg
  character*999              :: gridfile=""
  character*999              :: rivfile=""
  character*999              :: welfile=""
  character*999              :: segfile=""
  character*999              :: timefile=""
  real                       :: nodefac, tmp(5)
  integer                    :: n, pos

  character(len=:), allocatable :: line

  ! read main input
  open(newunit=ifile, file=trim(configfile), status='old')
  do
    call get_line(ifile, line, ioerr)
    if (ioerr /= 0) exit
    line = adjustl(line)
    if (line(1:1) /= '#' .and. line(1:1) /= '!') exit
  end do

  do pos = 1, 5
    read(line,*,iostat=ioerr) tmp(1:pos)  ! it automatically fails if non-numeric characters are encountered
    if (ioerr==0) then
      n=pos
    else
      exit
    end if
  end do
  if (n>=1) then; lfac    = tmp(1); else; lfac    = 1.0  ; end if
  if (n>=2) then; ddtor   = tmp(2); else; ddtor   = 1e-3 ; end if
  if (n>=3) then; dptor   = tmp(3); else; dptor   = 1e-10; end if
  if (n>=4) then; rimg    = tmp(4); else; rimg    = 0.5  ; end if
  if (n>=5) then; nodefac = tmp(5); else; nodefac = 0.0  ; end if
  ! sanity check
  if (nodefac<0) nodefac = 0
  if (nodefac>1) nodefac = 1
  if (ddtor  <  0.0) error stop "ddtor   must be positive or zero."
  if (dptor  <  0.0) error stop "dptor   must be positive or zero."
  if (rimg   <= 0.0) error stop "rimg    must be positive."
  
  if (lfac < 0) then
    read(ifile,*) gridfile
  else
    read(ifile,*) cellsize1, cellmutl, gridrmax
  end if
  read(ifile,*) contourlevel
  read(ifile,*) welfile
  read(ifile,*) rivfile
  ! read(ifile,*) segfile
  read(ifile,*) timefile
  read(ifile,*) outfile
  close(ifile)

  ! read grid coordinates
  if (lfac < 0) then
    print*, "Reading grid coordiantes in "//trim(gridfile)
    ncell = linecount(gridfile) - 1
    hasNodeName = columncount(gridfile)>2
    allocate(xcell(ncell), ycell(ncell))
    open(newunit=ifile, file=trim(gridfile), status='old')
    if (hasNodeName) then
      allocate(nodenames(0:ncell))
      read(ifile,*) nodenames(0),nodenames(0),nodenames(0)
    else
      read(ifile,*)
    end if
    if (hasNodeName) then
      do ii=1, ncell
        read(ifile, *) xcell(ii), ycell(ii), nodenames(ii)
      end do
    else
      do ii=1, ncell
        read(ifile, *) xcell(ii), ycell(ii)
      end do
    end if
    close(ifile)
  end if

  ! read well file
  print*, "Reading wells in "//trim(welfile)
  nwel = linecount(welfile) - 1
  allocate(well(nwel))
  open(newunit=ifile, file=trim(welfile), status='old')
  read(ifile,*)
  do iwel=1, nwel
    read(ifile,*) well(iwel)%x, well(iwel)%y, well(iwel)%trans, well(iwel)%sto, well(iwel)%q, well(iwel)%timeoff
    well(iwel)%iwel  = iwel
    well(iwel)%qsign = sign(1.0, well(iwel)%q)
    well(iwel)%q = abs(well(iwel)%q)
    well(iwel)%nriv = 0
  end do
  close(ifile)

  ! read riv coordinates
  print*, "Reading stream nodes in "//trim(rivfile)
  nriv = linecount(rivfile) - 1
  if (nriv>0) then
    if (nriv==1) onereach = .true.
    allocate(xriv(nriv), yriv(nriv), x1riv(nriv), y1riv(nriv), x2riv(nriv), y2riv(nriv), iseg(nriv), lriv(nriv))
    open(newunit=ifile, file=trim(rivfile), status='old')
    read(ifile,*)
    do iriv=1, nriv
      read(ifile, *) x1riv(iriv), y1riv(iriv), x2riv(iriv), y2riv(iriv), iseg(iriv)
    end do
    xriv = x1riv * nodefac + (1.0 - nodefac) * x2riv
    yriv = y1riv * nodefac + (1.0 - nodefac) * y2riv
    lriv = sqrt((x1riv - x2riv) ** 2 + (y1riv - y2riv) ** 2)
    nseg = maxval(iseg)
    close(ifile)
  end if

  ! read time intervals
  print*, "Reading times in "//trim(timefile)
  ntime = linecount(timefile) - 1
  allocate(times(ntime))
  open(newunit=ifile, file=trim(timefile), status='old')
  read(ifile,*)
  do itime=1, ntime
    read(ifile, *) times(itime)
  end do
  close(ifile)

  end subroutine readdata


  !> @brief Write a CSV-style file of x,y coordinates and optional values.
  !>
  !> Writes a header with column names and then each row with index, x, y and
  !> optional numeric columns `vals` when provided.
  !>
  !> @param[in] filename Output file path.
  !> @param[in] nrow Number of rows to write.
  !> @param[in] ncol Number of additional numeric columns per row.
  !> @param[in] x,y Coordinate arrays of length nrow.
  !> @param[in] columns Column header text appended to the header line.
  !> @param[in,optional] vals Optional real array (nrow x ncol) with values.
  subroutine write_vals(filename, nrow, ncol, x, y, columns, vals)
  character(len=*)        :: filename
  integer                 :: nrow, ncol
  character(len=*)        :: columns
  real                    :: x(:), y(:)
  real, optional          :: vals(nrow,ncol)
  ! print*, "weights", weights
  open(newunit=ifile, file=trim(filename), status="replace")
  write(ifile, "(A)") "idx,"//trim(columns)
  do ii=1, nrow
    if (ncol>0 .and. present(vals)) then
      write(ifile, "(I0,',',G0.8,',',G0.8,*(:',',G0.8))") ii, x(ii), y(ii), vals(ii, :)
    else
      write(ifile, "(I0,',',G0.8,',',G0.8,*(:',',G0.5))") ii, x(ii), y(ii)
    end if
  end do
  close(ifile)
  end subroutine


  !> @brief Write a CSV with a row name column plus coordinates and values.
  !>
  !> Similar to `write_vals` but includes a row-name string per row.
  !>
  !> @param[in] filename Output file path.
  !> @param[in] nrow Number of rows to write.
  !> @param[in] ncol Number of numeric columns per row.
  !> @param[in] x,y Coordinate arrays of length nrow.
  !> @param[in] columns Column header text.
  !> @param[in] vals Numeric array (nrow x ncol) with values.
  !> @param[in] rownames String array (0:nrow) containing the column header at index 0 and each row name.
  subroutine write_vals_rowname(filename, nrow, ncol, x, y, columns, vals, rownames)
  character(len=*)        :: filename
  integer                 :: nrow, ncol
  character(len=*)        :: columns
  character(len=*)        :: rownames(0:nrow)
  real                    :: x(:), y(:)
  real                    :: vals(nrow,ncol)
  ! print*, "weights", weights
  open(newunit=ifile, file=trim(filename), status="replace")
  write(ifile, "(A)") "idx,"//trim(rownames(0))//","//trim(columns)
  do ii=1, nrow
    write(ifile, "(I0,',',A,',',G0.8,',',G0.8,*(:',',G0.5))") ii, trim(rownames(ii)), x(ii), y(ii), vals(ii, :)
  end do
  close(ifile)
  end subroutine


end module
