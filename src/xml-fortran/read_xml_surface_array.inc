! Part of XML-Fortran library:
!
! $Id: read_xml_array.inc,v 1.3 2007/02/26 20:33:38 arjenmarkus Exp $
!
   type(XML_PARSE), intent(inout)                :: info
   character(len=*), intent(in)                  :: tag
   logical, intent(inout)                        :: endtag
   character(len=*), dimension(:,:), intent(in)  :: attribs
   integer, intent(in)                           :: noattribs
   character(len=*), dimension(:), intent(in)    :: data
   integer, intent(in)                           :: nodata
   logical, intent(inout)                        :: has_var

   character(len=len(attribs(1,1)))              :: buffer
   integer                                       :: idx
   integer                                       :: ierr
   integer                                       :: i,j

   !
   ! The big trick:
   ! A string long enough to hold all data strings
   !
   character(len=nodata*(len(data(1))+1))        :: bufferd
   character(len=11)                             :: op
   integer                                       :: start


   !
   ! The value can be stored in an attribute values="..." or in
   ! the data
   !
   has_var = .false.
   idx = xml_find_attrib( attribs, noattribs, 'values', buffer )
   if ( idx .gt. 0 ) then
      call read_from_buffer( buffer, var, ierr )
      if ( buffer .ne. ' ' ) then
         has_var = .true.
      endif
   else
      bufferd = ' '
      start   = 1
      j = 1
      do idx = 1,nodata
         if ( data(idx) .ne. ' ' ) then
            ! add to buffer char by char, allowing for preprocessing of
            ! openmc operators into their corresponding integer ids
            !write(*,*),"data ",data(idx)
            do i=1,len(data(idx))
                if (data(idx)(i:i) == '(') then
                    write(op,'(I11)') huge(0)
                    bufferd(j:j+11) = op
                    j = j + 11
                else if (data(idx)(i:i) == ')') then
                    write(op,'(I11)') huge(0)-1
                    bufferd(j:j+11) = op
                    j = j + 11
                else if (data(idx)(i:i) == ':') then
                    write(op,'(I11)') huge(0)-2 ! union
                    bufferd(j:j+11) = op
                    j = j + 11
                else if (data(idx)(i:i) == '#') then
                    write(op,'(I11)') huge(0)-3 ! difference
                    bufferd(j:j+11) = op
                    j = j + 11
                else
                    bufferd(j:j) = data(idx)(i:i)
                    j = j + 1
                end if
            end do
         endif
      enddo
      call read_from_buffer( bufferd, var, ierr )
      if ( bufferd .ne. ' ' ) then
         has_var = .true.
      endif
   endif

   if ( ierr .ne. 0 ) then
      write(*,*) 'Error reading variable - tag = ', trim(tag)
      has_var = .false.
   endif
