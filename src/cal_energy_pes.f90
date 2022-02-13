! cal pes energy , prepare input.pqr *.pqr *.log, change MAXRES MAXAT if necessary
! method if DFT fot default, if use different method , change reading of SCF Done
! by HYQ,2019
module abc
integer,parameter::MAXAT=90000,MAXRES=6000
integer firstprotres,lastprotres,resno(MAXRES)
integer resatom(MAXRES),presatom(MAXRES,100),presatom_AB(MAXRES,100),resatom_cap(MAXRES),presatom_cap(MAXRES,100)
real(kind=8)::efragment(MAXRES),efragment_chg(MAXRES),econcap(MAXRES),econcap_chg(MAXRES)
real(kind=8)::epair(MAXRES),epair_chg(MAXRES),eresidue(MAXRES),eresidue_chg(MAXRES)
real(kind=8)::Ecap_A_Cap,Ecap_Cap,E_pair,Eres(MAXRES),E_qkql,E_qiqj,E_total,Ecap_AB_Cap
real(kind=8)::force(MAXAT,3),temp(3)
character(len=80):: systemname
integer natom,nres,modum, ist
character(len=4)::atomname(MAXAT)
character(len=3)::residue(MAXAT),residuename(MAXAT)
character(len=2)::element(MAXAT)
real(kind=8)::coord(3,MAXAT)
real(kind=8)::qmcharge(MAXAT),rad(MAXAT),tcharge
logical::normal,countatom(MAXRES,100),countatom_cap(MAXRES,100)

endmodule

Program main
use abc
implicit none
systemname = 'input.pqr' ! pqr filename of the full system, cant be protein.pqr
call fileRead
call cap_AB_Cap
call cap_A_Cap
!call concap
call mmpart
call cal_E_total
endprogram main

!==================sum energy of cap_A_Cap================================
!
! read pqr file
!
subroutine fileRead
use abc
implicit none
integer i,j,k, nptemp, ttnumber
character*80 line
character*6 sn(MAXAT)
double precision chargef(MAXRES)

! 读取时原子序号是自动累加的，不使用文件中的,
! 残基序号必须是连续的
open(701,file=trim(adjustl(systemname)))
chargef = 0.d0
nptemp = 0
i = 0
do k=1,MAXAT
read(701,'(a80)',end = 101) line
    if(line(1:5) .eq. 'ATOM ') then
        i = i + 1
        read(line,100)sn(i),ttnumber,atomname(i),residue(i),resno(i),&
            (coord(j,i),j=1,3),qmcharge(i),rad(i),element(i)
100 format(a6,i5,1x,a4,1x,a3,2x,i4,4x,3f8.3,f8.4,f8.3,7x,a2)
        nptemp=resno(i)
        chargef(nptemp) = chargef(nptemp) + qmcharge(i)
        if(nptemp.gt.MAXRES) then
            write(0,*) 'too many residues: ', nptemp, MAXRES
            stop
        endif
        residuename(nptemp)=residue(i)
    endif
enddo
101 natom = i
firstprotres=resno(1)
lastprotres=resno(natom)  !the number of last protein residue, not total number of residues
nres = resno(natom) - resno(1) +1  !total number of res
close(701)

do i=1,natom
    if(resno(i) == firstprotres)then
        resno(i)=firstprotres+1
    endif
    if(resno(i) == lastprotres)then
        resno(i)=lastprotres-1
    endif
    if(element(i).eq.'  ') then
        element(i)(1:1) = atomname(i)(2:2)
    endif
enddo

firstprotres = firstprotres +1
lastprotres = lastprotres -1
end subroutine fileRead

subroutine cap_AB_Cap
use abc
implicit none 
integer i,j,k
character*80 filei,fileiline
Ecap_AB_Cap = 0.0
force = 0.0d0
countatom = .false.
do i=firstprotres,lastprotres-1
    filei = 'fragment2b'//char(48+i/100) &
        //char(48+(i-i/100*100)/10)//char(48+(i-i/10*10))
    open(41,file=filei(1:13)//'.pqr')
    j=0
    do while (.true.)
        read(41,'(a80)',iostat=ist) fileiline
        if(ist/=0)exit
        if((fileiline(1:4)=='ATOM'))then
            j=j+1
            read(fileiline,'(6x,i5)')presatom_AB(i,j)
            if((fileiline(18:20)/='MOD'))then
                countatom(i,j) = .true.
            endif
        endif
    enddo
    close(41)
    resatom(i)=j

    open(40,file=filei(1:13)//'.log')
    normal = .true.

        do while(.true.)
            read(40,'(a80)',iostat=ist) fileiline
            if(ist/=0)exit
            if(fileiline(1:7).eq.'energy:')then
                read(40,*) efragment(i)
            endif
            if(fileiline(1:6).eq.'force:')then
                do j = 1, resatom(i)
                    read(40,*) (temp(k),k=1,3)
                    if(countatom(i,j))then
                        !write(409,'(i6,i8,3f15.9)') i,presatom_AB(i,j),(temp(k),k=1,3)
                        force(presatom_AB(i,j),1) = force(presatom_AB(i,j),1) + temp(1)
                        force(presatom_AB(i,j),2) = force(presatom_AB(i,j),2) + temp(2)
                        force(presatom_AB(i,j),3) = force(presatom_AB(i,j),3) + temp(3)
                    endif
                enddo
                normal = .false.
                exit
            endif
        enddo
    ! endif

    close(40)
    if(normal)then
        write(233,*)'error in ',filei
    endif

    write(6,2001)'fragment2b',i,'QM energy=',efragment(i)
2001 format(a10,i3,1x,a16,f16.8)
    Ecap_AB_Cap = Ecap_AB_Cap + efragment(i)
enddo
open(99,file='force.dat0')
do i = 1, natom
    write(99,'(I5,3x,3f15.9)')i,(force(i,k),k=1,3)
enddo
close(99)
end subroutine cap_AB_Cap

subroutine cap_A_Cap
use abc
implicit none 
integer i,j,k
character*80 filei,fileiline
Ecap_A_Cap = 0.0
countatom = .false.
do i=firstprotres+1,lastprotres-1
    filei = 'fragment'//char(48+i/100) &
        //char(48+(i-i/100*100)/10)//char(48+(i-i/10*10))
    open(41,file=filei(1:11)//'.pqr')
    j=0
    do while (.true.)
        read(41,'(a80)',iostat=ist) fileiline
        if(ist/=0)exit
        if((fileiline(1:4)=='ATOM'))then
            j=j+1
            read(fileiline,'(6x,i5)')presatom(i,j)
            if((fileiline(18:20)/='MOD'))then
                countatom(i,j) = .true.
            endif
        endif
    enddo
    close(41)
    resatom(i)=j

    open(40,file=filei(1:11)//'.log')
    normal = .true.
    do while(.true.)
        read(40,'(a80)',iostat=ist) fileiline
        if(ist/=0)exit
        if(fileiline(1:7).eq.'energy:')then
            read(40,*) efragment(i)
        endif
        if(fileiline(1:6).eq.'force:')then
            do j = 1, resatom(i)
                read(40,*) (temp(k),k=1,3)
                if(countatom(i,j))then
                    !write(409,'(i6,i8,3f15.9)') i,presatom(i,j),(temp(k),k=1,3)
                    force(presatom(i,j),1) = force(presatom(i,j),1) - temp(1)
                    force(presatom(i,j),2) = force(presatom(i,j),2) - temp(2)
                    force(presatom(i,j),3) = force(presatom(i,j),3) - temp(3)
                endif
            enddo
            normal = .false.
            exit
        endif
    enddo


    close(40)
    if(normal)then
        write(233,*)'error in ',filei
    endif

    write(6,2000)'fragment',i,'QM energy=',efragment(i)
2000 format(a8,i3,1x,a16,f16.8)
    Ecap_A_Cap = Ecap_A_Cap + efragment(i)
enddo
open(99,file='force.dat00')
do i = 1, natom
    write(99,'(I5,3x,3f15.9)')i,(force(i,k),k=1,3)
enddo
close(99)
end subroutine cap_A_Cap
!
!========================MM region====================================
!
subroutine mmpart
use abc
implicit none
integer i,j,k
real(kind=8)::dis, x1,y1,z1, x2,y2,z2, fele(MAXAT,3)
character(len=9) aa,bb
open(80,file='fort.501')
read(80,'(a9,f18.10)') aa,E_qkql
write(*,*)'MM energy = ',0.529158477*E_qkql
open(81,file='fort.203')
fele = 0.0d0
do while(.true.)
    read(81,*,iostat=ist) i,j
    if(ist/=0)exit
    x1 = coord(i,1)
    y1 = coord(i,2)
    z1 = coord(i,3)
    x2 = coord(j,1)
    y2 = coord(j,2)
    z2 = coord(j,3)
    dis = dsqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    fele(i,1) = fele(i,1) + qmcharge(i)*qmcharge(i)/(dis**3)*(x1-x2)
    fele(i,2) = fele(i,2) + qmcharge(i)*qmcharge(i)/(dis**3)*(y1-y2)                        
    fele(i,3) = fele(i,3) + qmcharge(i)*qmcharge(i)/(dis**3)*(z1-z2)
    fele(j,1) = fele(j,1) + qmcharge(j)*qmcharge(j)/(dis**3)*(x2-x1)
    fele(j,2) = fele(j,2) + qmcharge(j)*qmcharge(j)/(dis**3)*(y2-y1)                        
    fele(j,3) = fele(j,3) + qmcharge(j)*qmcharge(j)/(dis**3)*(z2-z1)
enddo
close(81)
    !Unit conversion
fele = fele*0.5291771d0/627.509474d0

!open(88,file='force.dat1')
!open(99,file='force.dat11')
do i = 1, natom
!    write(88,'(I5,3x,3f15.9)')i,(force(i,k),k=1,3)
!    write(99,'(I5,3x,3f15.9)')i,(fele(i,k),k=1,3)
    force(i,1)=force(i,1)+fele(i,1)
    force(i,2)=force(i,2)+fele(i,2)
    force(i,3)=force(i,3)+fele(i,3)
enddo
!close(88)
!close(99)
endsubroutine mmpart
!
!=====================EE-GMFCC formula===============================
!
!E_total = EcapAcap-Econcap+Epair-mmpart
subroutine cal_E_total
use abc
open(99,file='force.dat')
do i = 1, natom
    write(99,'(I5,3x,3f15.9)')i,(force(i,k),k=1,3)
enddo
close(99)
E_total = 0.d0
E_total = Ecap_AB_Cap - Ecap_A_Cap + 0.529158477*E_qkql
open(88,file='energy.dat')
write(88,'(a16,f20.10)') ' Total energy = ', E_total
write(*,*) 'Total energy = ', E_total
write(88,'(a16,f20.10)') 'onebody energy= ',Ecap_AB_Cap-Ecap_A_Cap
write(88,'(a16,f20.10)') '    MM energy = ',0.529158477*E_qkql
close(88)
endsubroutine cal_E_total

    
    
    