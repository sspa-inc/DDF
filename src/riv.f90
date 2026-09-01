module m_river
  implicit none
  integer              :: nriv, nseg                ! number of stream reaches/nodes, number of stream segments
  integer, allocatable :: iseg(:)                   ! segment index for each stream reach
  real, allocatable    :: xriv(:),x1riv(:),x2riv(:) ! x coordinates of a stream reach mid point, start point, end point
  real, allocatable    :: yriv(:),y1riv(:),y2riv(:) ! y coordinates of a stream reach mid point, start point, end point
  real, allocatable    :: lriv(:)                   ! length of a stream reach
end module

