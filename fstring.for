*====================================================================
* Module origname:        FString
* Module specs:           Text file reading functions
* Project:                Units
* Functions:              5
* Style warnings:         0
*====================================================================
* TODO:                   0
*====================================================================

* LINK WITH: -


*--------------------------------------------------------------------
*	int iReadNonEmptyStringEx(int ifs, char *buf,
*					char *comment, int ierr);
*	int iReadNonEmptyString(int ifs, char *buf);
*
* iReadNonEmptyStringEx() reads line from input file stream 'ifs' and
* stores it into an external buffer pointed to by 'buf'. String parts
* starting from token 'comment' to end of line are ignored.
* Error messages are redirected to the output stream 'ierr'.
*
* Return Value: length of resulting string (written into 'buf')
*
* Errors: unexpected EOF or other file error. On error, a message
* is written to the error stream, empty string is returned in 'buf'
* and 0 for the function return value.
*
* Bug: comment substring will be considered as comment start even if
* placed inside quotes.
*
* iReadNonEmptyString() is a shortcut with comment='#' and ierr=6.
*--------------------------------------------------------------------
	integer function iReadNonEmptyString(ifs,buf)
	integer ifs	! input stream
c	integer ibuf_maxl
	character buf*(*)	!(ibuf_maxl)
	iReadNonEmptyString=iReadNonEmptyStringEx(ifs,buf,'#',6)
	return
	end
*--------------------------------------------------------------------
	integer function iReadNonEmptyStringEx(ifs,buf,comment,ierr)

	integer ifs	! input stream
c	integer ibuf_maxl
* I don't know why, but it doesn't work with strings longer
* than 256 bytes, hanging up the machine (!)
	character buf*(256) !(ibuf_maxl)
	character comment*(*)

	buf=' '
	dowhile (lentrim(buf).eq.0)
	   read(ifs,'(a256)',err=201,end=202) buf
	   if (lentrim(buf).eq.0) continue
d	   write(ierr,*) 'have read "',buf(1:lentrim(buf)),'"'
	   m=index(buf,comment)
	   if (m.gt.1) buf=buf(1:m-1)
	   if (m.eq.1) buf=' '
	enddo
  204	format(a256)

	! Remove nonprinting characters (TAB) here.
	dowhile (lentrim(buf).gt.0.and.buf(1:1).le.char(32))
	  buf=buf(2:)
	enddo
	do i=1,lentrim(buf)
	  if (buf(i:i).le.char(32)) buf=buf(1:i-1)//char(32)//buf(i+1:)
	enddo
* "buf(m:n)" selects n characters starting from m-th position

	goto 203

  201	if (ierr.gt.0) write(ierr,*) 'Error reading input stream ',ifs
	buf=' '
c	exit
	goto 203
  202	if (ierr.gt.0) write(ierr,*) 'Unexpected end of file ',ifs
	buf=' '
c	exit
	goto 203

  203	iReadNonEmptyStringEx=lentrim(buf)
d	write(ierr,*) 'returning "',buf(1:lentrim(buf)),'"'
	return
	end

*--------------------------------------------------------------------
*       void WritePercent(int istr, int idone, int itotal);
*
* Writes the value 100*idone/itotal with '%' sign in the beginning
* of the line.
*--------------------------------------------------------------------
	subroutine WritePercent(istr,idone,itotal)

	write(istr,401) char(13),
     .	idnint(100.0d0*dble(idone)/dble(itotal))
  401	format(a1,i3,'%',\)

	return
	end

*====================================================================

	integer function iFirstWordLen(buf)

	integer ret
	character buf*(256)

	call BTrimEx(buf,256)
	ret=0
	dowhile (ret.lt.256.and.buf(ret+1:ret+2).gt.char(32))
	  ret=ret+1
	enddo

d	write(6,*) 'iFW: "',buf(1:ret),'"',ret
	iFirstWordLen=ret
	return
	end

*====================================================================
*       int iIniFileOpenSection(char *file_name, char *section_name,
*				int file_number, int stop_if_failed);
*
* Opens the file specified by  file_name  and positions the file
* pointer in the beginning of  section_name  section.
*
* The number of Fortran stream used to work with the file is set by
* file_number and returned upon success.
*
* If the requested section has not been opened (no file, etc.) and
* stop_if_failed was set to a nonzero, the program terminates.
*
* IniFile is a text file with data entries "x=value" arranged into
* sections:
*       [section1]
*       twopi=6.28
*       infile=my_data.dat
*       [section2]
*       m=12
*       .....
*--------------------------------------------------------------------
	integer function iIniFileOpenSection(file_name,section_name,
     +	ifile_number,ifail)

	integer ifs,ret,ifail
	character file_name*(*),section_name*(*)
	character buf*(256),secn*(256)

	ret=0
	ifs=ifile_number
	secn=section_name
d	write(6,*) 'file "',file_name,'", section "',section_name,'"'
d	write(6,*) 'conv. section "',secn(1:lentrim(secn)),'"'
	open(ifs,file=file_name,form='formatted',err=203)
	dowhile (iReadNonEmptyStringEx(ifs,buf,';',0).gt.0)
d	   write(6,*) 'want "','['//secn(1:lentrim(secn))//']',
d     . '", got "',buf(1:lentrim(buf)),'"'
	   if (buf(1:lentrim(buf)).eq.
     . '['//secn(1:lentrim(secn))//']') then
d		write(6,*) 'positioned after "',buf(1:lentrim(buf)),'"'
		ret=ifs
		exit
	   endif
	enddo
	if (ret.eq.0) close(ifs)

  203	iIniFileOpenSection=ret

c Terminate the application if the successful opening of the section
c was marked as critical by a nonzero ifail:
	if (ret.eq.0.and.ifail.ne.0) then
	write(*,*) 'IniRead: Can''t open section [',section_name,
c     +	secn(1:lentrim(secn)),
     +	'] in project file ',file_name	!(1:lentrim(file_name))
	stop ' '
	endif

	return
	end


*--------------------------------------------------------------------
*       iIniFileReadEntryString(int ifs, char *buf);
*
* Reads next entry string from IniFile pointed to by  ifs  opened on
* correct section.
* Returns the length of the resulting string.
* On EOF or section end  buf='' and 0 is returned.
*--------------------------------------------------------------------
	integer function iIniFileReadEntryString(ifs,buf)

	integer ifs,ret
	character buf*(256)

	ret=0
	if (ifs.ne.0) ret=iReadNonEmptyStringEx(ifs,buf,';',0) ! output stream =0
d	write(6,*) 'IFReadEntryString got "',buf(1:lentrim(buf)),'"',ret
	if (ret.eq.0.or.index(buf,'[').eq.1) then
	   buf=' '
	   ret=0
d	   write(6,*) 'closing stream ',ifs
	   if (ifs.ne.0) close(ifs)     !!!
	endif

	iIniFileReadEntryString=ret
	return
	end

*====================================================================

	integer function iStringToDoubleArray(buf,dest,dest_maxl,ierr)

	character buf*(*)
	double precision dest(*)
	integer dest_maxl ! maximum number of items read
	integer ierr	! error output stream

	integer m,n

d	write(ierr,*) 'iS2DA("',buf(1:lentrim(buf)),'"):'

	n=0
	m=iFirstWordLen(buf)
	dowhile (m.gt.0.and.n.lt.dest_maxl)
	   dest(n+1)=0.0d0
	   read(buf(1:m),*,err=204) dest(n+1)
	   goto 206
  204	if (ierr.gt.0) write(ierr,*)
     +	'Error reading double precision value'
  206	buf=buf(m+1:)
d	   write(ierr,*) n,dest(n+1)
	   n=n+1
	   m=iFirstWordLen(buf)
	enddo

c	goto 207
c  205	if (ierr.gt.0) write(ierr,*) 'Unexpected end of input stream'
  207	iStringToDoubleArray=n
d	write(ierr,*) 'iS2DA returned ',n,' values:'
d	write(ierr,*) (dest(i),i=1,n)
	return
	end

*--------------------------------------------------------------------
	integer function iStringToIntArray(buf,dest,dest_maxl,ierr)

	character buf*(*)
	integer dest(*)
	integer dest_maxl ! maximum number of items read
	integer ierr	! error output stream

	integer m,n

d	write(ierr,*) 'iS2DA("',buf(1:lentrim(buf)),'"):'

	n=0
	m=iFirstWordLen(buf)
	dowhile (m.gt.0.and.n.lt.dest_maxl)
	   read(buf(1:m),*,err=204) dest(n+1)
	   goto 206
  204	if (ierr.gt.0) write(ierr,*)
     +	'Error reading integer value'
  206	buf=buf(m+1:)
d	   write(ierr,*) n,dest(n+1)
	   n=n+1
	   m=iFirstWordLen(buf)
	enddo

c	goto 207
c  205	if (ierr.gt.0) write(ierr,*) 'Unexpected end of input stream'
  207	iStringToIntArray=n
	return
	end

*--------------------------------------------------------------------
	logical function GetLogical(s,defaultval)
	character*(*) s
	character m
	logical ret,defaultval

	m=s(1:1)
d	print *,m,index('YyÑ§Ttà®1',m),index('Nnç≠Ffã´0',m)

	ret=defaultval
	if (index('YyÑ§Ttà®1',m).gt.0) then
	  ret=.TRUE.
	elseif (index('Nnç≠Ffã´0',m).gt.0) then
	  ret=.FALSE.
	else
	  print *,'GetLogical: §Æ´¶≠Æ °Î´Æ °Î‚Ï Æ§≠Æ ®ß YÑTà1/NçFã0; ',
     +	'Ø‡®≠Ô‚Æ ',defaultval
	endif

d	print *,'GetLogical returns ',ret
	GetLogical=ret
	return
	end
