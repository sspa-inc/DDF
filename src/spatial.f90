module m_spatial
  implicit none
  contains

  !> @brief Compute Euclidean distance between two points.
  !>
  !> @param[in] x1,y1 Coordinates of the first point.
  !> @param[in] x2,y2 Coordinates of the second point.
  !> @return Euclidean distance between the two points.
  elemental function cdist(x1,y1,x2,y2) result(distances)
  real, intent(in)  :: x1,y1,x2,y2
  real              :: distances
  distances = sqrt((x1 - x2)**2 + (y1 - y2)**2)
  end function cdist

  !> @brief Pairwise distances between points in arrays xs, ys.
  !>
  !> Returns an allocated array containing distances for each unique pair
  !> (i,j) with i<j in lexicographic order.
  !>
  !> @param[in] xs,xs Arrays of x and y coordinates of length n.
  !> @return Allocated real array of length n*(n-1)/2 with pairwise distances.
  function pdist(xs, ys) result(distances)
  real, intent(in)  :: xs(:), ys(:)
  real, allocatable :: distances(:)
  integer           :: n, k

  integer           :: ii, ij
  n = size(xs)
  allocate(distances((n-1)*n/2))
  k = 0
  do ii = 1, n-1
    do ij = ii + 1, n
      k = k + 1
      distances(k) = sqrt((xs(ii) - xs(ij))**2 + (ys(ii) - ys(ij))**2)
    end do
  end do
  end function pdist

  ! https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
  ! Given three collinear points p, q, r, the function checks if
  ! point q lies on line segment 'pr'
  !> @brief Check whether point q lies on segment pr (collinear case).
  !>
  !> @param[in] px,py Coordinates of p.
  !> @param[in] qx,qy Coordinates of q.
  !> @param[in] rx,ry Coordinates of r.
  !> @return .true. if q is on segment pr (inclusive endpoints).
  elemental function onSegment(px, py, qx, qy, rx, ry)
  logical                 :: onSegment
  real, intent(in)        :: px, py, qx, qy, rx, ry
  if ( (qx <= max(px, rx)) .and. (qx >= min(px, rx)) .and. &
    (qy <= max(py, ry)) .and. (qy >= min(py, ry))) onSegment = .true.
  onSegment = .false.
  end function onSegment

  !> @brief Determine orientation of three ordered points (p,q,r).
  !>
  !> Returns 0 for collinear, 1 for clockwise, and 2 for counter-clockwise.
  !>
  !> @param[in] px,py Coordinates of p.
  !> @param[in] qx,qy Coordinates of q.
  !> @param[in] rx,ry Coordinates of r.
  !> @return 0 (collinear), 1 (clockwise), or 2 (counterclockwise).
  elemental function  orientation(px, py, qx, qy, rx, ry)
  integer                 :: orientation
  real, intent(in)        :: px, py, qx, qy, rx, ry
  ! to find the orientation of an ordered triplet (p,q,r)
  ! function returns the following values:
  ! 0 : Collinear points
  ! 1 : Clockwise points
  ! 2 : Counterclockwise

  ! See https://www.geeksforgeeks.org/orientation-3-ordered-points/amp/
  ! for details of below formula.
  real                    :: val

  val = ((qy - py) * (rx - qx)) - ((qx - px) * (ry - qy))
  if (val > 0) then

    ! Clockwise orientation
    orientation = 1
  else if (val < 0) then

    ! Counterclockwise orientation
    orientation = 2
  else
    ! Collinear orientation
    orientation = 0
  end if
  end function  orientation

  !> @brief Return .true. when two points are closer than a tiny tolerance.
  !>
  !> Useful to avoid false intersections caused by numerical noise.
  !>
  !> @param[in] p1x,p1y First point.
  !> @param[in] p2x,p2y Second point.
  !> @return .true. if points are within tolerance.
  elemental function close_point(p1x,p1y,p2x,p2y)
  logical                 :: close_point
  real, intent(in)        :: p1x,p1y,p2x,p2y
  real, parameter         :: tor_touch=1e-8

  close_point = abs(p1x-p2x)<tor_touch .and. abs(p1y-p2y)<tor_touch

  end function close_point


  !> @brief Determine whether two line segments p1q1 and p2q2 intersect.
  !>
  !> Uses orientation tests and special-case checks for collinearity to
  !> determine whether two segments intersect (including endpoints).
  !>
  !> @param[in] p1x,p1y,q1x,q1y Endpoints of the first segment.
  !> @param[in] p2x,p2y,q2x,q2y Endpoints of the second segment.
  !> @return .true. if the segments intersect, .false. otherwise.
  elemental function doIntersect(p1x,p1y,q1x,q1y,p2x,p2y,q2x,q2y)

  ! Find the 4 orientations required for
  ! the general and special cases
  logical                 :: doIntersect
  real, intent(in)        :: p1x,p1y,q1x,q1y,p2x,p2y,q2x,q2y
  integer                 :: o1,o2,o3,o4

  doIntersect = .false.

  ! check bounding box
  if (p1x<p2x .and. p1x<q2x .and. q1x<p2x .and. q1x<q2x) return
  if (p1x>p2x .and. p1x>q2x .and. q1x>p2x .and. q1x>q2x) return
  if (p1y<p2y .and. p1y<q2y .and. q1y<p2y .and. q1y<q2y) return
  if (p1y>p2y .and. p1y>q2y .and. q1y>p2y .and. q1y>q2y) return

  o1 = orientation(p1x, p1y, q1x, q1y, p2x, p2y)
  o2 = orientation(p1x, p1y, q1x, q1y, q2x, q2y)
  o3 = orientation(p2x, p2y, q2x, q2y, p1x, p1y)
  o4 = orientation(p2x, p2y, q2x, q2y, q1x, q1y)

  ! General case
  if ((o1 /= o2) .and. (o3 /= o4)) then
    doIntersect = .true.
    return
  end if

  ! Special Cases

  ! p1 , q1 and p2 are collinear and p2 lies on segment p1q1
  if ((o1 == 0) .and. onSegment(p1x, p1y, p2x, p2y, q1x, q1y)) then
    doIntersect = .true.
    return
  end if

  ! p1 , q1 and q2 are collinear and q2 lies on segment p1q1
  if ((o2 == 0) .and. onSegment(p1x, p1y, q2x, q2y, q1x, q1y)) then
    doIntersect = .true.
    return
  end if

  ! p2 , q2 and p1 are collinear and p1 lies on segment p2q2
  if ((o3 == 0) .and. onSegment(p2x, p2y, p1x, p1y, q2x, q2y)) then
    doIntersect = .true.
    return
  end if

  ! p2 , q2 and q1 are collinear and q1 lies on segment p2q2
  if ((o4 == 0) .and. onSegment(p2x, p2y, q1x, q1y, q2x, q2y)) then
    doIntersect = .true.
    return
  end if
  end function doIntersect

end module
