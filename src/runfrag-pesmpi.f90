program main
!use mpi
implicit none
include "mpif.h"

integer(kind=4)::ierr, nid, num_procs, ist, i
character(len=3)::jobname(21)
character(len=7):: jobname2(9)
character(len=80)::fline
logical::frag,frag2

! initialize mpi 
! ierr for fortran only
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, nid, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)

jobname = (/ 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', &
            'GLY', 'HIE', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', &
            'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'cap' /)
jobname2 = (/ 'ALA-ALA', 'ALA-GLY', 'ALA-PRO', 'GLY-ALA', 'GLY-GLY', &
            'GLY-PRO', 'PRO-ALA', 'PRO-GLY', 'PRO-PRO' /)

do i = nid+1, 30, num_procs
frag=.false.
frag2=.false.
open(11,file='pesfrag')
if(i<22)then
    open(22,file='input.'//jobname(i))
    do while(.true.)
        read(11,'(A80)',iostat=ist)fline
        if(ist/=0)exit
        if(fline(1:8)=='fragment'.and.fline(9:10)/='2b'.and.fline(12:14)==jobname(i))then
            write(22,'(a)')fline(1:11)//'.xyz'
            frag = .true.
        endif
        if(fline(1:3)=='cap'.and.jobname(i)=='cap')then
            write(22,'(a)')fline(1:6)//'.xyz'
            frag = .true.
        endif
    enddo
    close(22)
    if(frag)then
        call system('./package/'//jobname(i)//'/nn.sh')
    endif
else
    open(22,file='input.'//jobname2(i-21))
    write(*,*)'input.'//jobname2(i-21)
    do while(.true.)
        read(11,'(A80)',iostat=ist)fline
        if(ist/=0)exit
        if(fline(1:10)=='fragment2b'.and.fline(14:20)==jobname2(i-21))then
            write(22,'(a)')fline(1:13)//'.xyz'
            frag2 = .true.
        endif
    enddo
    close(22)
    if(frag2)then
        call system('./package/'//jobname2(i-21)//'/nn.sh')
    endif
endif
close(11)
enddo
call MPI_Barrier(MPI_COMM_WORLD,ierr)
end
