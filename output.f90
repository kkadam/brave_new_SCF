!*************************************************************************
!*
!*  OUTPUT
!*
!*************************************************************************
subroutine output(frnum, out_template, data_array)
implicit none
include 'runscf.h'
!*************************************************************************
!*
!*  Global Variables


!*
!************************************************************************* 
!*
!*   Local variables

character(len=54) :: out_file
character :: numfr*(4)                     
!*
!*************************************************************************
!*
!*   Subroutine Arguments

integer :: frnum, i, j

character(len=50) :: out_template

real, dimension(numr,numz,numphi) :: data_array

!*
!*************************************************************************
! Initialize local variables
!out_template   = 'frame'

! create the filenames for the files every pe is going
!write(out_file,'(a,i4)') trim(out_template),frnum

!    write (numfr,900) frnum
!    900 format (I4,F8.3)	
	
	
out_file="invisible_cows"
print*,out_file
open(unit=50,file=trim(out_file),form='unformatted',convert='BIG_ENDIAN',status='unknown') 

write(50) data_array

close(50)

       open(unit=10,file="star1")
         do j=1,numz
           do i=1,numr  
             write(10,*) i,j,data_array(i,j,1) 
           enddo
           write(10,*)
         enddo
       close(10)    
       
       
       open(unit=10,file="star2")
         do j=1,numz
           do i=1,numr  
             write(10,*) i,j,data_array(i,j,numphi/2) 
           enddo
           write(10,*)
         enddo
       close(10)        
       
       
       open(unit=15,file='density',form='unformatted',status='unknown')
       write(15) data_array
       close(15) 	
	
return
end subroutine output
