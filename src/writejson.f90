!> @brief Write the opening JSON FeatureCollection header to unit iout.
!>
!> @param[in] iout Fortran unit number to write the JSON text to.
subroutine jsonhead(iout)
integer     :: iout
write(iout, "(A)") '{'
write(iout, "(A)") '  "type": "FeatureCollection",'
write(iout, "(A)") '  "features": ['
end subroutine

!> @brief Append a LineString feature (with Level property) to the JSON features array.
!>
!> This writes a single LineString feature representing a river segment or line between
!> (x1,y1) and (x2,y2) with a numeric property "Level" set to z. If z equals 0.0,
!> nothing is written.
!>
!> @param[in] iout Fortran unit number to write the JSON text to.
!> @param[in] x1,x2   Coordinates of the first and second point (X).
!> @param[in] y1,y2   Coordinates of the first and second point (Y).
!> @param[in] z       Numeric Level property; if zero the feature is skipped.
subroutine jsonadd(iout,x1,y1,x2,y2,z)
integer     :: iout
real        :: x1,y1,x2,y2,z
if(z==0.0) return
! if (first) then
!     write(iout, "(A)")                 '      {'
!     first = .false.
! else
!     write(iout, "(A)")                 '      , {'
! end if
! write(iout, "(A)")                 '        "type": "Feature",'
! write(iout, "(A)")                 '        "geometry": {'
! write(iout, "(A)")                 '          "type": "LineString",'
! write(iout, "(A)")                 '          "coordinates": ['
! write(iout, "(A,G0.8,',',G0.8,A)") '            [',x1,y1,'],'
! write(iout, "(A,G0.8,',',G0.8,A)") '            [',x2,y2,']'
! write(iout, "(A)")                 '          ]'
! write(iout, "(A)")                 '        },'
! write(iout, "(A,G0.5,A)")          '        "properties": {"Level": ',z,'}'
! write(iout, "(A)")                 '      }'

write(iout, "(A,4(G0.8,A),G0.5,A)") &
  '      {"type": "Feature","geometry": {"type": "LineString","coordinates": [[',&
  x1,',',y1,'],[',x2,',',y2,']]},"properties": {"Level": ',z,'}},'
end subroutine

!> @brief Write the terminating JSON elements for the FeatureCollection.
!>
!> This closes the features array and the outer JSON object. The additional
!> coordinate arguments are currently unused but kept for API compatibility
!> with other json write helpers.
!>
!> @param[in] iout Fortran unit number to write the JSON text to.
!> @param[in] x1,x2,y1,y2,z Unused coordinates/level parameters (kept for compatibility).
subroutine jsonend(iout)
integer     :: iout
write(iout, "(A)") '      {}'
write(iout, "(A)") '  ]'
write(iout, "(A)") '}'
end subroutine