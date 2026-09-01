module m_well
  implicit none
  real, parameter         :: pi = 3.1415926535897932384626433
  real, parameter         :: pi2 = sqrt(pi)
  real, parameter         :: gamma = 0.5772156649015328606
  real                    :: rimg = 0.5         ! minimum distance (well radius)
  real                    :: ddtor = 1e-3       ! minimum drawdown induced by pumping well at the stream node to be included in calculation
  Real                    :: dptor = 1e-10      ! minimum leakage rate from a stream node to be included in calculation
  real, parameter         :: extention=5280*100 ! extend one hundred miles
  logical                 :: onereach
  integer                 :: maxriv = -999      ! maximum number of river segments that can be used by a well
  type t_well
    integer              :: iwel
    real                 :: x
    real                 :: y
    real                 :: q                 ! positive for pumping; negative for injection
    real                 :: timeoff
    real                 :: qsign
    real                 :: sto
    real                 :: trans

    ! below are used during computation and only need to be calculated once
    real, allocatable    :: r_riv(:),r_cel(:),r_img_cel(:,:),r_img_riv(:,:),l_riv(:)

    integer              :: ncell
    logical, allocatable :: mask_cel(:)
    integer, allocatable :: cidx(:)

    integer              :: nriv
    logical, allocatable :: mask_riv(:)
    integer,allocatable  :: ridx(:)
    logical, allocatable :: mask_node(:)

  contains
  procedure            :: filter_riv
  procedure            :: filter_node
  procedure            :: filter_side
  procedure            :: riv_radius
  procedure            :: set_dist
  procedure            :: theis
  procedure            :: solve_weight
  procedure            :: calc_dd
  end type t_well
  integer                 :: nwel
  type(t_well),allocatable:: well(:)
  contains

  !> @brief Convert a logical mask array to a 1-based index array of true elements.
  !>
  !> Packs the indices i (1..size(mask)) where mask(i) is true into an
  !> integer, allocatable array and returns it.
  !>
  !> @param[in] mask Logical array to convert.
  !> @return Integer array of indices where mask is true.
  function mask2index(mask)
  integer,allocatable  :: mask2index(:)
  logical,intent(in)   :: mask(:)
  integer              :: ii
  mask2index = pack([(ii,ii=1,size(mask))], mask)
  end function

  !> @brief Theis solution for drawdown at radius r and time t for this well.
  !>
  !> This elemental function evaluates the Theis well solution using the
  !> well's pumping rate and aquifer properties stored in the t_well object.
  !>
  !> @param[in] self Class instance providing q, sto, trans, etc.
  !> @param[in] r Distance from well (length units consistent with trans/sto).
  !> @param[in] t Time since pumping started.
  !> @return Drawdown value (same units as head).
  elemental function theis(self, r, t)
  class(t_well), intent(in)     :: self
  real                          :: theis
  real, intent(in)              :: r, t
  theis = self%q * wellfuc1(calc_u(r, t, self%sto, self%trans)) / (4.0*pi*self%trans)
  end function

  !> @brief Compute dimensionless parameter u used in well functions.
  !>
  !> @param[in] r Radius.
  !> @param[in] t Time.
  !> @param[in] sto Storage coefficient.
  !> @param[in] trans Transmissivity.
  !> @return Dimensionless u = r^2 * sto / (4 * trans * t).
  elemental function calc_u(r, t, sto, trans)
  ! well function 1
  real              :: calc_u
  real, intent(in)  :: r, t, sto, trans
  calc_u = r*r*sto/(4.0*trans*t)
  end function

  !> @brief Approximation of the exponential integral / well function for Theis.
  !>
  !> Provides an approximation for the function used in the Theis solution.
  !>
  !> @param[in] U Dimensionless argument (calc_u result).
  !> @return Value of the well function used by the Theis solution.
  elemental function wellfuc1(U)
  ! well function 1
  real              :: wellfuc1
  real, intent(in)  :: U
  real              :: U2,U3,U4
  real              :: a
  integer           :: n
  U2 = U*U
  if (U .GT. 1.0) then
    wellfuc1=(exp(-U)/U)*(U2+2.334733*U+.250621)/(U2+3.330657*U + 1.681534)
  else
    U3 = U2*U
    U4 = U3*U
    wellfuc1=-log(U)-gamma+.99999193*U-.24991055*U2+.05519968*U3-.00976004*U4+.00107857*U4*U
  endif
  ! wellfuc1 = -gamma-log(U)+U
  ! a = U
  ! do n=1, 20
  !   a = -a*U*n/(n+1)**2
  !   wellfuc1 = wellfuc1+a
  ! end do
  end function

  !> @brief Estimate an effective radius for a river reach that the drawdown at this radius matches the average drawdown of the reach.
  !>
  !> Computes an approximate radius at which the drawdown equals the average
  !> drawdown of the reach. Currently the routine returns a minimum radius
  !> `rimg` due to accuracy considerations.
  !>
  !> @param[in] self Class instance with aquifer parameters.
  !> @param[in] l Length of the stream reach.
  !> @param[in] t Time since pumping started.
  !> @return Effective radius to use for river image-well computations.
  elemental function riv_radius(self, l, t)
  class(t_well), intent(in)  :: self
  real, intent(in)           :: l, t ! stream reach length and time
  real                       :: riv_radius

  real                       :: dd_avg, fpiT, ll, log_term
  riv_radius = rimg
  return ! experiment shows that using average drawdown decrease accuracy
  ll = l / 2.0
  fpiT = 4.0 * pi * self%trans
  dd_avg = (self%q / fpiT) * ( -gamma + log( (4.0 * self%trans * t) / (self%sto * ll**2) ) + 1.0 )
  log_term = (fpiT / self%q) * dd_avg + gamma
  riv_radius = max(sqrt((4.0 * self%trans * t / self%sto) * exp(-log_term)), rimg)
  end function


  !> @brief Extend a river segment endpoints upstream and downstream by an incremental length.
  !>
  !> Moves the two endpoints (x1,y1),(x2,y2) outward along the segment direction
  !> by a total of `increase` (split equally upstream/downstream).
  !>
  !> @param[inout] x1,y1 Coordinates of the first endpoint (modified in place).
  !> @param[inout] x2,y2 Coordinates of the second endpoint (modified in place).
  !> @param[in] increase Length to extend the segment by.
  subroutine extend_riv(x1, y1, x2, y2, increase)
  ! extend the river section upstream and downstream by a length of `increase`.
  real, intent(inout)        :: x1, y1, x2, y2
  real, intent(in   )        :: increase
  real        :: dx, dy, dl
  dx = x2 - x1
  dy = y2 - y1
  dl = sqrt(dx**2 + dy**2)
  dx = dx * increase / dl
  dy = dy * increase / dl
  x1 = x1 - dx
  y1 = y1 - dy
  x2 = x2 + dx
  y2 = y2 + dy
  end subroutine

  !> @brief check if any of the stream node of a segment is on the same side of well
  !>
  !> if all stream nodes of a segment are on different side, this segment will be excluded from calculation.
  subroutine filter_riv(self)
  ! filter river
  use m_spatial, only   : cdist, doIntersect
  use m_river
  ! use m_io,      only   : write_vals
  class(t_well)        :: self
  ! local
  integer              :: i,j,k
  real                 :: r_wel_riv(nriv)
  logical              :: sameside(nriv)

  r_wel_riv = cdist(self%x, self%y, xriv, yriv)
  ! step 1 find useful riv points
  allocate(self%mask_riv(nriv))
  if (onereach) then
    self%mask_riv = .false.
    self%mask_riv(minloc(r_wel_riv, 1)) = .true.
  else
    sameside = .true.
    sloop: do k = 1, nseg
      rloop: do i=1, nriv
        ! check if any of the RIV node of the current segment is on the same side of well,
        ! if so, this segment needs to be included
        if (iseg(i) /= k) cycle rloop
        jloop: do j=1, nriv
          if (iseg(j) == k) cycle jloop ! do not check using its own reach
          if (doIntersect(x1riv(j),y1riv(j),x2riv(j),y2riv(j),self%x,self%y,xriv(i),yriv(i))) then
            sameside(i) = .false.
            ! print*, "Well", self%iwel, " excludes Segment", k, " through Node", i, " by Reach", j
            exit jloop ! found this river node is on the other side; move to next river node
          end if
        end do jloop
      end do rloop
      ! if (all(pack()))
    end do sloop
    self%mask_riv = sameside
  end if

  self%nriv = count(self%mask_riv)
  self%ridx = mask2index(self%mask_riv)
  self%r_riv = r_wel_riv(self%ridx)
  self%l_riv = lriv(self%ridx)
  if (onereach) then
    ! extend the segment like infinitely long
    do i=1, self%nriv
      j = self%ridx(i)
      call extend_riv(x1riv(j),y1riv(j),x2riv(j),y2riv(j),extention)
    end do
  end if
  ! call write_vals("riv_active.csv", self%nriv, 0, xriv(self%ridx), yriv(self%ridx), "x,y")
  end subroutine filter_riv


  !> @brief Return true if point (xs,ys) is on the same side of the river as the well.
  !>
  !> Uses river segment intersection tests to decide whether the provided point
  !> lies on the same side as the well instance.
  !>
  !> @param[in] self Well instance to compare against.
  !> @param[in] xs,ys Point coordinates to test.
  !> @return .true. when the point is on the same side; .false. otherwise.
  elemental function filter_side(self, xs, ys)
  ! check if xs,ys is on the same side of the river as the well
  use m_river
  use m_spatial, only         : doIntersect
  class(t_well),intent(in)   :: self
  logical                    :: filter_side
  real, intent(in)           :: xs, ys
  integer                    :: k, iriv
  filter_side = .false.
  do iriv=1, self%nriv
    k = self%ridx(iriv)
    if (doIntersect(x1riv(k),y1riv(k),x2riv(k),y2riv(k),self%x,self%y,xs,ys)) return
  end do
  filter_side = .true.
  end function filter_side

  !> @brief Compute distances between the well, cells, and stream nodes.
  !>
  !> Populates r_cel, r_riv, r_img_cel and r_img_riv arrays and related index/mask
  !> data needed for later drawdown calculations.
  !>
  !> @param[inout] self Well object; its distance and index fields are updated.
  subroutine set_dist(self)
  ! calculate all distance used later
  use m_spatial, only   : cdist
  use m_river
  use m_grid
  ! use m_io,      only   : write_vals
  class(t_well)        :: self

  ! local
  integer              :: iimg
  real, allocatable    :: ximg(:),yimg(:),dimg(:)

  ! filter cells and calculate distance between well and cells

  ! call write_vals("active.csv", self%ncell, 0, xcell(self%cidx), ycell(self%cidx), "x,y")
  ! place the (imaginary) wells
  if (self%nriv>0) then
    self%mask_cel = self%filter_side(xcell, ycell)
    self%ncell = count(self%mask_cel)
    self%cidx = mask2index(self%mask_cel)
    self%r_cel = cdist(self%x, self%y, xcell(self%cidx), ycell(self%cidx))
    self%r_cel = abs(lfac)*self%r_cel
    allocate(ximg(self%nriv),yimg(self%nriv),dimg(self%nriv))
    if (onereach) then
      ximg = 2.0 * xriv(self%ridx) - self%x
      yimg = 2.0 * yriv(self%ridx) - self%y
    else
      ximg = xriv(self%ridx) - self%x
      yimg = yriv(self%ridx) - self%y
      dimg = sqrt(ximg**2 + yimg**2)
      ximg = xriv(self%ridx) + rimg * ximg/dimg
      yimg = yriv(self%ridx) + rimg * yimg/dimg
      ! ximg = xriv(self%ridx)
      ! yimg = yriv(self%ridx)
    end if
    ! call write_vals("image.csv", self%nriv, 0, ximg, yimg, "x,y")

    ! calculate the distance from riv/cell to the (imaginary) wells
    allocate(self%r_img_riv(self%nriv, self%nriv))
    allocate(self%r_img_cel(self%ncell, self%nriv))
    do iimg = 1, self%nriv
      self%r_img_riv(:, iimg) = cdist(ximg(iimg), yimg(iimg), xriv (self%ridx), yriv (self%ridx))
      self%r_img_cel(:, iimg) = cdist(ximg(iimg), yimg(iimg), xcell(self%cidx), ycell(self%cidx))
    end do

    ! apply length unit conversion factor
    self%r_riv    =abs(lfac)*self%r_riv
    self%r_img_riv=abs(lfac)*self%r_img_riv
    self%r_img_cel=abs(lfac)*self%r_img_cel
    ! make sure there is no distance smaller than radius `rmig`
    where(self%r_cel<rimg)    self%r_cel=rimg
    where(self%r_riv<rimg)     self%r_riv=rimg
    where(self%r_img_riv<rimg) self%r_img_riv=rimg
    where(self%r_img_cel<rimg) self%r_img_cel=rimg
  else
    allocate(self%mask_cel(ncell)); self%mask_cel=.true.
    self%ncell = ncell
    self%cidx = mask2index(self%mask_cel)
    self%r_cel = cdist(self%x, self%y, xcell, ycell)
    self%r_cel = abs(lfac)*self%r_cel
    where(self%r_cel<rimg)    self%r_cel=rimg
  end if

  end subroutine

  !> @brief Solve for weights/leakages of stream nodes that balance river drawdown.
  !>
  !> Solves a linear system to find non-negative weights on stream nodes so that
  !> the resulting drawdown at river nodes matches the drawdown from the pumping
  !> well. Negative weights are set to zero and renormalized.
  !>
  !> @param[in] self Well instance with precomputed distances.
  !> @param[in] t Time at which to evaluate drawdown and weights.
  !> @return Array of weights for river nodes (length self%nriv).
  function solve_weight(self, t)
  class(t_well)           :: self
  real                    :: t
  real                    :: solve_weight(self%nriv)
  ! local
  logical                 :: pumpoff
  real                    :: wtot,ddwel(self%nriv)
  integer                 :: ii, neq, info, irivmin
  real   , allocatable    :: b(:),ddriv(:,:),matA(:,:),ddmaxriv(:)
  integer, allocatable    :: idx(:), ipiv(:)
  logical                 :: ddmask(self%nriv)

  solve_weight = 0.0
  if (onereach) then
    solve_weight(1) = 1.0
    return
  end if

  pumpoff = self%timeoff>0 .and. t>self%timeoff
  ddwel = self%theis(self%r_riv, t)
  ! if (pumpoff) ddwel = ddwel - self%theis(self%r_riv, t-self%timeoff)

  ! mask of stream nodes that should be included
  ddmask = (ddwel>=ddtor) .and. (self%mask_node)

  ! filter the `maxriv` most influential river segments
  neq = count(ddmask)
  if (maxriv>0) then
      if (neq>maxriv) then
        allocate(ddmaxriv(maxriv))
        irivmin = 0
        do ii = 1, self%nriv
            if (ddmask(ii)) then
                irivmin = irivmin + 1
                ddmaxriv(irivmin) = ddwel(ii)
            end if
            if (irivmin==maxriv) exit
        end do
        irivmin = minloc(ddmaxriv, 1)
        do ii = maxriv+1, self%nriv
            if (ddwel(ii)>ddmaxriv(irivmin)) then
              ddmaxriv(irivmin) = ddwel(ii)
              irivmin = minloc(ddmaxriv, 1)
            end if
        end do
        ddmask = ddwel >= ddmaxriv(irivmin)
        neq = count(ddmask)
      end if
  end if
  if (neq==0) return
  allocate(ipiv(neq), matA(neq,neq), ddriv(neq,neq), b(neq))
  idx = mask2index(ddmask)
  b(1:neq) = ddwel(idx)

  ! calculate weights for each riv segment that could cancel the effect of pumping well
  do ii = 1, neq
    ! self%r_img_riv(idx(ii),idx(ii)) = self%riv_radius(self%l_riv(idx(ii)), t)
    ddriv(:, ii) = self%theis(self%r_img_riv(idx,idx(ii)), t)
    ! if (pumpoff) ddriv(:, ii) = ddriv(:, ii) - self%theis(self%r_img_riv(idx,idx(ii)), t-self%timeoff)
  end do
  matA = max(0.0, ddriv)
  call sgesv(neq, 1, matA, neq, ipiv, b(1:neq), neq, info)
  ! call write_matrix(neq, ddriv, ddwel(idx))
  if (info==0) then
    ! weight correction by removing negative weights
    ! print*, "        ", "positive", sum(pack(b, b>0))
    ! print*, "        ", "negative", sum(pack(b, b<0))
    ! print*, "        ", "total   ", sum(b)
    wtot = sum(b)
    where(b<dptor) b=0.0
    b = b * wtot / sum(b)
    if (sum(b) /= sum(b)) b=0.0
    do ii = 1, neq
      solve_weight(idx(ii)) = b(ii)
    end do
    return
  else
    print*, "Well", self%iwel, "sgesv error at time", t, neq, self%nriv, info
    call write_matrix(neq, ddriv, ddwel(idx))
    stop
  end if

  end function

  !> @brief Return true if late time (t) depletion at a node is greater than or equal to dptor.
  !>
  !> Uses solve_weight function to calculate stream node weights at a late time period.
  !> Compares the node weights to dptor, and creates a true/false mask
  !>
  !> @param[in] self Well instance to compare against.
  !> @param[in] t Late-time period (in days) used to calculate weights.
  !> @return .true. when the node weight is greater than or equal to dptor; .false. otherwise.
  subroutine filter_node(self)

  use m_river
  class(t_well)        :: self
  real                 :: late_weights(self%nriv)
  real                 :: tmax            ! time for drawdown to reach the furthest stream nodex
  real                 :: rmax            ! maximum distance from active stream node to well center

  allocate(self%mask_node(self%nriv))
  self%mask_node = .true.

  if (dptor<=0.0) return
  ! filter the stream nodes positioned on the further side of active reaches from the well
  rmax = maxval(self%r_riv)

  ! back calculate the time to reach rmax using Cooper and Jacob approximation; 1.3 is an overshoot
  tmax = 1.3 * rmax*rmax * self%sto / (4.0 * self%trans) / exp(- ddtor * 4.0 * pi * self%trans / abs(self%q) + gamma)

  late_weights = self%solve_weight(tmax)
  self%mask_node = late_weights>=dptor
  end subroutine filter_node

  !> @brief Dump linear system matrices to disk for debugging.
  !>
  !> Writes matrix A and right-hand side B to simple text files for inspection.
  !>
  !> @param[in] neq Dimension of the square matrix.
  !> @param[in] matA Matrix of size (neq,neq) to write.
  !> @param[in] rhsB Right-hand side vector of length neq.
  subroutine write_matrix(neq, matA, rhsB)
  integer, intent(in) :: neq
  real, intent(in)    :: matA(neq, neq), rhsB(neq)

  integer             :: ii, ifile
  open(newunit=ifile, file="matA.dat", status="replace")
  do ii = 1, neq
    write(ifile, "(*(G10.3))") matA(:, ii)
  end do
  close(ifile)
  !open(newunit=ifile, file="imgD.dat", status="replace")
  !do ii = 1, neq
  !   write(ifile, "(*(G10.3))") self%r_img_riv(idx,idx(ii))
  !end do
  !close(ifile)
  open(newunit=ifile, file="rhsB.dat", status="replace")
  write(ifile, "(1G10.3)") rhsB
  close(ifile)
  end subroutine

  !> @brief Analytical cumulative Stream Depletion Factor (SDF) using erfc.
  !>
  !> Returns the cumulative stream depletion factor for a given transmissivity
  !> T, storativity S, reach radius r and time.
  !>
  !> @param[in] r Radius.
  !> @param[in] t Time.
  !> @param[in] sto Storage coefficient.
  !> @param[in] trans Transmissivity.
  !> @return Cumulative SDF (dimensionless).
  elemental function cum_sdf(r, t, sto, trans)
  ! calculate the Stream Depletion Factor (cumulative stream depletion) by integrting error function
  real, intent(in)           :: r, t, sto, trans
  real                       :: cum_sdf

  ! local
  real                       :: u2, u

  u = calc_u(r, t, sto, trans)
  u2 = sqrt(u)

  cum_sdf = (1.0 + 2.0*u) * erfc(u2) - 2.0 / pi2 * u2 * exp(-u)
  end function

  !> @brief Numerical integration of cumulative Stream Depletion Factor (SDF).
  !>
  !> Performs an adaptive trapezoidal integration of the error-function based
  !> SDF when the analytical expression is unsuitable or unstable.
  !>
  !> @param[in] T Transmissivity.
  !> @param[in] S Storativity.
  !> @param[in] r Reach effective radius.
  !> @param[in] time Time over which to integrate.
  !> @return Numerical estimate of cumulative SDF.
  function cum_sdf_numerical(T, S, r, time)
  ! calculate the Stream Depletion Factor (cumulative stream depletion) by integrting error function
  real, intent(in)           :: T, S, r, time
  real                       :: cum_sdf_numerical
  ! local
  real, parameter            :: dt_fac=1.05 ! time step increase factor
  real                       :: tmin        ! Small lower bound to avoid division by zero
  real                       :: alpha, tt, dt, sdr1, sdr2
  integer                    :: n

  alpha = (S * r**2) / (4.0 * T)
  tmin = alpha / 200.0
  n = 0
  cum_sdf_numerical = 0.0
  sdr2 = erfc(sqrt(alpha / tmin))
  dt = tmin
  tt = tmin
  do while (tt < time)
    dt = min(dt, time-tt)
    tt = tt + dt
    sdr1 = sdr2
    sdr2 = erfc(sqrt(alpha / tt))
    cum_sdf_numerical = cum_sdf_numerical + 0.5 * (sdr1 + sdr2) * dt
    if (sdr2/sdr1<dt_fac .or. dt/time<0.0001) dt = dt * dt_fac ! increase step size when function plateau or dt is small
    ! Uncomment below for debugging:
    ! n = n + 1
    ! print *, 'Step:', n, 'Time:', tt, 'dt:', dt, 'Value:', sdr2
  end do
  cum_sdf_numerical = cum_sdf_numerical / time
  end function

  !> @brief Calculate drawdown contribution from this well at time t.
  !>
  !> Computes drawdown at all grid cells due to this well, accounting for
  !> stream nodes leakage. Optionally returns per-river-node leakage.
  !>
  !> @param[in] self Well instance with precomputed distances and weights.
  !> @param[in] t Time at which to compute drawdown.
  !> @param[out] ddnet Array (ncell) filled with drawdown values.
  !> @param[out, optional] leakage If present, filled with leakage per river node.
  subroutine calc_dd(self, t, ddnet, leakage)
  use m_river, only           : nriv
  use m_grid , only           : ncell
  class(t_well),intent(in)    :: self
  real, intent(in)            :: t
  real, intent(out), optional :: leakage(:)
  real, intent(out)           :: ddnet(ncell)

  ! local
  real                        :: weights(self%nriv)
  real                        :: ddwel(self%ncell)
  real                        :: ddimg(self%ncell)
  real                        :: dd(self%ncell)
  integer                     :: iriv

  if (nriv>0) leakage = 0.0
  ddnet = 0.0
  if (self%nriv>0) then
    weights = self%solve_weight(t)
  else
    weights = 0.0
  end if
  ! if (itime==30) call write_vals("riv_weight.csv", nriv, 1, xriv, yriv, ["x                   ","y                   ","weight"], weights)
  ddwel = self%theis(self%r_cel, t)
  ddimg = 0.0
  do iriv = 1, self%nriv
    ! if (weights(iriv)<dptor) then
    !   weights(iriv) = 0
    !   cycle
    ! end if
    dd = self%theis(self%r_img_cel(:,iriv), t)
    ddimg = ddimg + weights(iriv) * dd
  end do
  dd = ddwel - ddimg
  ! print*, t, self%r_cel(1414), ddwel(1414), ddimg(1414), dd(1414)
  ! print*, self%r_img_cel(1414,:)
  ! stop
  where(dd<ddtor) dd=0.0
  if (self%nriv>0) then
    if (onereach) then
      weights(1) = cum_sdf(self%r_riv(1), t, self%sto, self%trans)
      ! weights(1) = cum_sdf_numerical(self%trans, self%sto, self%r_riv(1), t)
    end if
    leakage = unpack(weights * self%qsign * self%q, self%mask_riv, 0.0)
  end if
  ddnet = unpack(dd * self%qsign, self%mask_cel, 0.0)
  end subroutine calc_dd

  !> @brief Calculate the well center and maximum radius for well influence.
  !>
  subroutine calc_grid_parameters()
    use m_spatial, only    : pdist
    use m_grid   , only    : xcenter, ycenter, gridrmax, lfac
    use m_time   , only    : times
    real                   :: tmax, rmax
    tmax = max(maxval(times), maxval(well%timeoff))
    rmax = maxval(sqrt(exp(-(ddtor / abs(well%q) * (4.0*pi*well%trans) + gamma)) * 4.0 * well%trans * tmax / well%sto))*1.3 ! Cooper and Jacob approximation
    if (nwel>1) rmax = rmax + maxval(pdist(well%x, well%y)*abs(lfac))
    gridrmax = max(rmax, gridrmax)
    xcenter = sum(well%x) / nwel
    ycenter = sum(well%y) / nwel
  end subroutine calc_grid_parameters
end module
