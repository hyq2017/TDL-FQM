module comanmr
implicit none
integer,parameter::MAXAT=10000,MAXRES=6000
integer lastprotatom,firstprotres,lastprotres,resno(MAXAT),charge(MAXRES)
integer cout1(MAXAT,MAXAT),cout2(MAXAT,MAXAT)
integer natom,nres,modum
character(len=4)::atomname(MAXAT)
character(len=5)::dlabel(MAXAT)
character(len=3)::residue(MAXAT),residuename(MAXAT)
character(len=2)::element(MAXAT)
real(kind=8)::coord(3,MAXAT)
real(kind=8)::qmcharge(MAXAT),rad(MAXAT),tcharge
logical::connect(MAXRES,MAXRES)
logical::atomsign(MAXAT),mmsign(MAXAT)
character*80 filek,filej,filel,filem, gmemory, gnproc, gmethod,systemname
character(len=4) basename
logical::ter(0:MAXRES), wg09, wxyz, wprotein
integer qselectC(MAXRES),qselectCA(MAXRES),qselectCB(MAXRES),qselectN(MAXRES)
end module comanmr

Program main
use comanmr

!============== user change part ==========================
    systemname = 'input.pqr' ! pqr filename of the full system, cant be protein.pqr
    gmemory = '%mem=3GB'
    gnproc = '%nprocshared=4'
    gmethod = '#p wB97XD/6-31g* nosymm SCF=tight force'
    wg09 = .true.  !whether to write com and pqr file of frags
    wxyz = .true. !whether to write xyz file of frags
    wprotein = .true. !whether to write com file of the protein 

!==================== end ============================

    if(wg09)then
        open(901,file='g09frag')
    endif
    if(wxyz)then
        open(902,file='pesfrag')
    endif
    call fileRead
    if(wprotein)then
        call protwrite !写出全原子体系输入文件
    endif
    call fragselect
    call fragCreate
    if(wg09)then
        write(901,*)
        close(901)
    endif
    if(wxyz)then
        write(902,*)
        close(902)
    endif
endprogram main
!==========================subroutines====================
!
! read pqr file
!
subroutine fileRead
use comanmr
implicit none
integer i,j,k
integer nptemp
!integer ist
integer ttnumber
character*80 line
character*6 sn(MAXAT)
double precision chargef(MAXRES)

! 读取时原子序号是自动累加的，不使用文件中的,
! 残基序号必须是连续的
100 format(a6,i5,1x,a4,1x,a3,2x,i4,4x,3f8.3,f8.4,f8.3,7x,a2)

ter = .false.
chargef = 0.d0
charge = 0
!lastprotatom = 0 ! 暂时搁置
nptemp = 0
i = 0
open(701,file=trim(adjustl(systemname)))
do k=1,MAXAT
    read(701,'(a80)',end = 101) line
    if(line(1:5) .eq. 'ATOM ') then
        i = i + 1
        read(line,100)sn(i),ttnumber,atomname(i),residue(i),resno(i),&
            (coord(j,i),j=1,3),qmcharge(i),rad(i),element(i)
        nptemp=resno(i)
        chargef(nptemp) = chargef(nptemp) + qmcharge(i)
        if(nptemp.gt.MAXRES) then
            write(0,*) 'too many residues: ', nptemp, MAXRES
            stop
        endif
        residuename(nptemp)=residue(i)
        !if(resno(i) .eq. lastprotres) then
        !    lastprotatom = i
        !endif
    else if (line(1:3).eq.'TER') then
        ter(nptemp) = .true.
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
    if(element(i).eq.'  ')then
        element(i)(1:1) = atomname(i)(2:2)
    endif
enddo

do i=firstprotres, lastprotres
    charge(i) = nint(chargef(i))
    !write(201,*) i,charge(i),natom

    do j=firstprotres, lastprotres
        connect(i,j)=.false.
    enddo
enddo

firstprotres = firstprotres +1
lastprotres = lastprotres -1
charge(firstprotres)=charge(firstprotres)+charge(firstprotres-1)
charge(lastprotres)=charge(lastprotres)+charge(lastprotres+1)
end subroutine fileRead

subroutine protwrite
use comanmr
implicit none
integer cfrag,i,iqm
cfrag=0
do i=firstprotres, lastprotres
    cfrag = cfrag + charge(i)
enddo

open(88,file='protein.com')
!write(901,'(a)') 'g09 protein.com'
open(89,file='protein.pqr')
write(88,'(a)') '%mem=30GB'
write(88,'(a)') trim(adjustl(gnproc))
write(88,'(a)') trim(adjustl(gmethod))
write(88,*)
write(88,'(a)') 'full protein '
write(88,*)
write(88,'(1x,i3,2x,i2)') cfrag,1
iqm = 0
do i=1,natom
    iqm = iqm + 1
    call addatom(i,iqm,88,89)
enddo
close(89)
write(88,*)
close(88)
end subroutine protwrite 

!
! prepare C and CA for making frag
!
subroutine fragselect
use comanmr
implicit none
integer i

do i=1,natom
    if(resno(i).le.lastprotres) then
        if(atomname(i).eq." C  ") qselectC(resno(i))=i
        if(atomname(i).eq." CA ") qselectCA(resno(i))=i
        if(atomname(i).eq." CB ") qselectCB(resno(i))=i
        if(atomname(i).eq." N  ") qselectN(resno(i))=i
!   else
!       if(resno(i-1).ne.resno(i)) qselectC(resno(i))=i
    endif
enddo
qselectN(lastprotres) = qselectCA(lastprotres)-2
end subroutine fragselect

subroutine fragCreate
use comanmr
implicit none
integer i,j,k,l,l1,l2,ll,iqm,cfrag,kk,iqmca,jj,kkatom,n1,n2,iiqm,llatom
integer kstart,kfinal,ktemp,iresm1,kbstart,kbfinal,kkstart,kkfinal
double precision x,y,z,dis,E_qkql,E_qiqj, E_qkql2

dis = 0.d0
E_qkql = 0.d0
E_qkql2 = 0.d0
E_qiqj = 0.d0
!==========================residues with caps=============================
do k=firstprotres+1,lastprotres-1
    iqm=0
    filek = 'fragment'//char(48+k/100) &
        //char(48+(k-k/100*100)/10)//char(48+(k-k/10*10))
    if(wg09)then
        write(901,'(a)') 'g09 '//filek(1:11)//'.com'
        open(30,file=filek(1:11)//'.com')
        open(31,file=filek(1:11)//'.pqr')
        write(30,'(a)') trim(adjustl(gmemory))
        write(30,'(a)') trim(adjustl(gnproc))
        write(30,'(a)') trim(adjustl(gmethod))
        write(30,*)
        write(30,'(a,i4)') 'GMFCC fragment ',k
    endif
    if(wxyz)then
        write(902,'(a)') filek(1:11)//residuename(k)//'.xyz'
        open(32,file=filek(1:11)//'.xyz')
    endif
    
    do i=1,natom
        atomsign(i)=.false.
    enddo

    cfrag=0
    kstart = qselectCA(k-1)
    kkstart = qselectC(k-1)
    kkfinal = qselectCA(k+1)+1
    if(k==firstprotres)then
        kkstart = 1
    elseif(k==lastprotres)then
        kkfinal = qselectC(k)+2
        if(residuename(k+1)=='NME')then
            kkfinal = qselectC(k)+7
        endif
    elseif(ter(k))then
        kkfinal = qselectC(k)+2
    elseif(ter(k-1))then
        kkstart = qselectC(k-1)+3
    endif

    if(residuename(k+1) == 'PRO')then
        kkfinal = qselectC(k)+2
        kfinal = qselectCA(k+1)
    endif

    call getcharge(kkstart,kkfinal,qmcharge,tcharge)
    if (k/=firstprotres)then
        tcharge = tcharge + qmcharge(kstart) + qmcharge(kstart+1)
    endif
    if(residuename(k+1) == 'PRO')then
        tcharge = tcharge + qmcharge(kfinal) +qmcharge(kfinal+1)
    endif
    cfrag = nint(tcharge)

    if(wg09)then
        write(30,*)
        write(30,'(1x,i3,2x,i2)') cfrag,1
    endif
!
!   QM region (one residue with caps)
!
    if(k/=firstprotres) then
        n1=qselectCA(k-1)
        n2=qselectN(k-1)
        call xyzchange(coord(1,n2),coord(2,n2),coord(3,n2), &
            coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
        iqm = iqm + 1
        if(wg09)then
            call addH(iqm,x,y,z,30,31)
        endif
        if(wxyz)then
            write(32,'(3f12.5)')x,y,z
        endif
    endif
    if(k/=firstprotres)then
        do kk=kstart,kstart+1
            atomsign(kk) = .true.
            iqm = iqm + 1
            if(wg09)then
                call addatom(kk,iqm,30,31)
            endif
            if(wxyz)then
                write(32,'(3f12.5)')(coord(j,kk),j=1,3)
            endif
        enddo
    endif
    if(k/=firstprotres) then
        n1=qselectCA(k-1)
        if(residuename(k-1) == 'PRO')then
            n2=qselectCB(k-1)
        else
            n2=qselectCA(k-1)+2
        endif
        call xyzchange(coord(1,n2),coord(2,n2),coord(3,n2), &
            coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
        iqm = iqm + 1
        if(wg09)then
            call addH(iqm,x,y,z,30,31)
        endif
        if(wxyz)then
            write(32,'(3f12.5)')x,y,z
        endif
    endif
    do kk=kkstart,kkfinal
        atomsign(kk) = .true.
        iqm = iqm + 1
        if(wg09)then
            call addatom(kk,iqm,30,31)
        endif
        if(wxyz)then
            write(32,'(3f12.5)')(coord(j,kk),j=1,3)
        endif
    enddo
    if(residuename(k+1) == 'PRO')then
        n1=qselectC(k)+2
        n2=qselectC(k)+3
        call xyzchangeN(coord(1,n2),coord(2,n2),coord(3,n2), &
            coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
        iqm = iqm + 1
        if(wg09)then
            call addH(iqm,x,y,z,30,31)
        endif
        if(wxyz)then
            write(32,'(3f12.5)')x,y,z
        endif
    endif
    if(residuename(k+1) == 'PRO')then
        do kk=kfinal,kfinal+1
            atomsign(kk) = .true.
            iqm = iqm + 1
            if(wg09)then
                call addatom(kk,iqm,30,31)
            endif
            if(wxyz)then
                write(32,'(3f12.5)')(coord(j,kk),j=1,3)
            endif
        enddo
    endif
!
!   add H
!
    if(k/=lastprotres) then
        n1=qselectCA(k+1)
        n2=qselectC(k+1)
        call xyzchange(coord(1,n2),coord(2,n2),coord(3,n2), &
            coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
        iqm = iqm + 1
        if(wg09)then
            call addH(iqm,x,y,z,30,31)
        endif
        if(wxyz)then
            write(32,'(3f12.5)')x,y,z
        endif
        n1=qselectCA(k+1)
        if(residuename(k+1) == 'PRO')then
            n2=qselectCB(k+1)
        else
            n2=qselectCA(k+1)+2
        endif
        call xyzchange(coord(1,n2),coord(2,n2),coord(3,n2), &
            coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
        iqm = iqm + 1
        if(wg09)then
            call addH(iqm,x,y,z,30,31)
        endif
        if(wxyz)then
            write(32,'(3f12.5)')x,y,z
        endif
    endif

    if(wg09)then
        close(31)
        write(30,*)
        close(30)    
    endif
    if(wxyz)then
        close(32)
    endif
!
!  mm region
!
    !
    ! to subtract duplicate mm interaction (will be used in energy calculation)
    !
    ! 1. count method
    do i=1,natom
        if(atomsign(i).eqv..true.) then
            do j=i+1,natom
                if(atomsign(j).eqv..true.) then
                    cout2(i,j) = cout2(i,j) + 1
                endif
            enddo
        endif
    enddo

enddo  !big loop 

!=============================neighbored 2body with caps===================
do k = firstprotres, lastprotres-1
    iqm=0
    filek = 'fragment2b'//char(48+k/100) &
        //char(48+(k-k/100*100)/10)//char(48+(k-k/10*10))
    if(wg09)then
        write(901,'(a)') 'g09 '//filek(1:13)//'.com'
        open(30,file=filek(1:13)//'.com')
        open(31,file=filek(1:13)//'.pqr')
        write(30,'(a)') trim(adjustl(gmemory))
        write(30,'(a)') trim(adjustl(gnproc))
        write(30,'(a)') trim(adjustl(gmethod))
        write(30,*)
        write(30,'(a,i4)') 'GMFCC fragment2b ',k
    endif
    if(wxyz)then
        write(902,'(a)') filek(1:13)//residuename(k)//"-"//residuename(k+1)//'.xyz'
        open(32,file=filek(1:13)//'.xyz')
    endif

    do i=1,natom
        atomsign(i)=.false.
    enddo

    cfrag=0
    kstart = qselectCA(k-1)
    kkstart = qselectC(k-1)
    kkfinal = qselectCA(k+2)+1

    if(k==firstprotres)then
        kkstart = 1
    elseif(k==lastprotres-1)then
        kkfinal = qselectC(k+1)+2
        if(residuename(k+2)=='NME')then
            kkfinal = qselectC(k+1)+7
        endif
    elseif(ter(k+1))then
        kkfinal = qselectC(k+1)+2
    elseif(ter(k-1))then
        kkstart = qselectC(k-1)+3
    endif

    if(residuename(k+2) == 'PRO')then
        kkfinal = qselectC(k+1)+2
        kfinal = qselectCA(k+2)
    endif

    call getcharge(kkstart,kkfinal,qmcharge,tcharge)
    if (k/=firstprotres)then
        tcharge = tcharge + qmcharge(kstart) + qmcharge(kstart+1)
    endif
    if(residuename(k+1) == 'PRO')then
        tcharge = tcharge + qmcharge(kfinal) +qmcharge(kfinal+1)
    endif
    cfrag = nint(tcharge)

    if(wg09)then
        write(30,*)
        write(30,'(1x,i3,2x,i2)') cfrag,1
    endif

!
!   QM region (one residue with caps)
!
    if(k/=firstprotres) then
        n1=qselectCA(k-1)
        n2=qselectN(k-1)
        call xyzchange(coord(1,n2),coord(2,n2),coord(3,n2), &
            coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
        iqm = iqm + 1
        if(wg09)then
            call addH(iqm,x,y,z,30,31)
        endif
        if(wxyz)then
            write(32,'(3f12.5)')x,y,z
        endif
    endif
    if(k/=firstprotres)then
        do kk=kstart,kstart+1
            atomsign(kk) = .true.
            iqm = iqm + 1
            if(wg09)then
                call addatom(kk,iqm,30,31)
            endif
            if(wxyz)then
                write(32,'(3f12.5)')(coord(j,kk),j=1,3)
            endif
        enddo
    endif
    if(k/=firstprotres) then
        n1=qselectCA(k-1)
        if(residuename(k-1) == 'PRO')then
            n2=qselectCB(k-1)
        else
            n2=qselectCA(k-1)+2
        endif
        call xyzchange(coord(1,n2),coord(2,n2),coord(3,n2), &
            coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
        iqm = iqm + 1
        if(wg09)then
            call addH(iqm,x,y,z,30,31)
        endif
        if(wxyz)then
            write(32,'(3f12.5)')x,y,z
        endif
    endif
    do kk=kkstart,kkfinal
        atomsign(kk) = .true.
        iqm = iqm + 1
        if(wg09)then
            call addatom(kk,iqm,30,31)
        endif
        if(wxyz)then
            write(32,'(3f12.5)')(coord(j,kk),j=1,3)
        endif
    enddo
    if(residuename(k+2) == 'PRO')then
        n1=qselectC(k+1)+2
        n2=qselectC(k+1)+3
        call xyzchangeN(coord(1,n2),coord(2,n2),coord(3,n2), &
            coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
        iqm = iqm + 1
        if(wg09)then
            call addH(iqm,x,y,z,30,31)
        endif
        if(wxyz)then
            write(32,'(3f12.5)')x,y,z
        endif
    endif
    if(residuename(k+2) == 'PRO')then
        do kk=kfinal,kfinal+1
            atomsign(kk) = .true.
            iqm = iqm + 1
            if(wg09)then
                call addatom(kk,iqm,30,31)
            endif
            if(wxyz)then
                write(32,'(3f12.5)')(coord(j,kk),j=1,3)
            endif
        enddo
    endif
!
!   add H
!
    if(k+1/=lastprotres) then
        n1=qselectCA(k+2)
        n2=qselectC(k+2)
        call xyzchange(coord(1,n2),coord(2,n2),coord(3,n2), &
            coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
        iqm = iqm + 1
        if(wg09)then
            call addH(iqm,x,y,z,30,31)
        endif
        if(wxyz)then
            write(32,'(3f12.5)')x,y,z
        endif

        n1=qselectCA(k+2)
        if(residuename(k+2) == 'PRO')then
            n2=qselectCB(k+2)
        else
            n2=qselectCA(k+2)+2
        endif
        call xyzchange(coord(1,n2),coord(2,n2),coord(3,n2), &
            coord(1,n1),coord(2,n1),coord(3,n1),x,y,z)
        iqm = iqm + 1
        if(wg09)then
            call addH(iqm,x,y,z,30,31)
        endif
        if(wxyz)then
            write(32,'(3f12.5)')x,y,z
        endif
    endif

    if(wg09)then
        close(31)
        write(30,*)
        close(30)    
    endif
    if(wxyz)then
        close(32)
    endif

!
!  mm region
!
    !
    ! to subtract duplicate mm interaction (will be used in energy calculation)
    !
    ! 1. count method
    do i=1,natom
        if(atomsign(i).eqv..true.) then
            do j=i+1,natom
                if(atomsign(j).eqv..true.) then
                    cout1(i,j) = cout1(i,j) + 1
                endif
            enddo
        endif
    enddo

    !2. k-l method
    if(k<lastprotres-1)then
        if(residuename(k+2)=='PRO')then
            l1 = kfinal - 1
            l2 = kfinal + 2
        else
            l1 = kkfinal + 1
            l2 = kkfinal + 2
        endif
        if(k>firstprotres)then
            do i=kstart, kstart+1
                do j =  kkfinal+1, l1
                    write(203,*) i,j
                    dis = dsqrt((coord(1,i)-coord(1,j))**2 &
                        +(coord(2,i)-coord(2,j))**2+(coord(3,i)-coord(3,j))**2)
                    E_qkql2 =  E_qkql2 + qmcharge(i)*qmcharge(j)/dis
                enddo
                do j =  l2, natom
                    write(203,*) i,j
                    dis = dsqrt((coord(1,i)-coord(1,j))**2 &
                        +(coord(2,i)-coord(2,j))**2+(coord(3,i)-coord(3,j))**2)
                    E_qkql2 =  E_qkql2 + qmcharge(i)*qmcharge(j)/dis
                enddo
            enddo
        endif
        n1 = qselectCA(k)-1
        n2 = qselectC(k)-1
        do i=kkstart, n1
            do j =  kkfinal+1, l1
                write(203,*) i,j
                dis = dsqrt((coord(1,i)-coord(1,j))**2 &
                    +(coord(2,i)-coord(2,j))**2+(coord(3,i)-coord(3,j))**2)
                E_qkql2 =  E_qkql2 + qmcharge(i)*qmcharge(j)/dis
            enddo
            do j =  l2, natom
                write(203,*) i,j
                dis = dsqrt((coord(1,i)-coord(1,j))**2 &
                    +(coord(2,i)-coord(2,j))**2+(coord(3,i)-coord(3,j))**2)
                E_qkql2 =  E_qkql2 + qmcharge(i)*qmcharge(j)/dis
            enddo
        enddo
        do i=n1+3, n2
            do j =  kkfinal+1, l1
                write(203,*) i,j
                dis = dsqrt((coord(1,i)-coord(1,j))**2 &
                    +(coord(2,i)-coord(2,j))**2+(coord(3,i)-coord(3,j))**2)
                E_qkql2 =  E_qkql2 + qmcharge(i)*qmcharge(j)/dis
            enddo
            do j =  l2, natom
                write(203,*) i,j
                dis = dsqrt((coord(1,i)-coord(1,j))**2 &
                    +(coord(2,i)-coord(2,j))**2+(coord(3,i)-coord(3,j))**2)
                E_qkql2 =  E_qkql2 + qmcharge(i)*qmcharge(j)/dis
            enddo
        enddo
    endif
enddo  !big loop k=firstprotres+1,lastprotres-1
write(501,'(a9,f18.10)') 'E_qkql = ',E_qkql2

end subroutine fragCreate

!=====================functions===========================
subroutine xyzchange(xold,yold,zold,xzero,yzero,zzero,xnew,ynew,znew)
implicit none
double precision xold,yold,zold,xzero,yzero,zzero,xnew,ynew,znew,grad

grad=dsqrt(1.09d0**2/((xold-xzero)**2+(yold-yzero)**2 &
    +(zold-zzero)**2))
xnew=xzero+grad*(xold-xzero)
ynew=yzero+grad*(yold-yzero)
znew=zzero+grad*(zold-zzero)
return   
end subroutine xyzchange

subroutine xyzchangeN(xold,yold,zold,xzero,yzero,zzero,xnew,ynew,znew)
implicit none
double precision xold,yold,zold,xzero,yzero,zzero,xnew,ynew,znew,grad

grad=dsqrt(1.00d0**2/((xold-xzero)**2+(yold-yzero)**2 &
    +(zold-zzero)**2))
xnew=xzero+grad*(xold-xzero)
ynew=yzero+grad*(yold-yzero)
znew=zzero+grad*(zold-zzero)
return   
end subroutine xyzchangeN

subroutine addatom(jj,jiqm,jfile1,jfile2)
use comanmr
implicit none
integer jj,jiqm,j,jfile1,jfile2
    if(jiqm.lt.10)then
        dlabel(jiqm)(1:1) = element(jj)(1:1)
        write(dlabel(jiqm)(2:2),'(i1)') jiqm
        dlabel(jiqm)(3:5) = '   '
    elseif(jiqm.lt.100) then
        dlabel(jiqm)(1:1) = element(jj)(1:1)
        write(dlabel(jiqm)(2:3),'(i2)') jiqm
        dlabel(jiqm)(4:5) = '  '
    elseif(jiqm.lt.1000) then
        dlabel(jiqm)(1:1) = element(jj)(1:1)
        write(dlabel(jiqm)(2:4),'(i3)') jiqm
        dlabel(jiqm)(5:5) = ' '
    elseif(jiqm.lt.10000) then
        dlabel(jiqm)(1:1) = element(jj)(1:1)
        write(dlabel(jiqm)(2:5),'(i4)') jiqm
    else
        write(0,*) 'Error: more than 100000 atoms'
        stop
    endif
    write(jfile1,'(a,2x,3f12.5)') dlabel(jiqm)(1:1),(coord(j,jj),j=1,3)
    write(jfile2,'(a,i5,1x,a4,1x,a3,i6,4x,3f8.3,f8.4,f8.3)') 'ATOM  ', &
        jj,atomname(jj),residue(jj),resno(jj),(coord(j,jj),j=1,3), &
        qmcharge(jj),rad(jj)
end subroutine addatom

subroutine addH(iqm,x,y,z,kfile1,kfile2)
    use comanmr
    implicit none
    integer iqm,kfile1,kfile2
    real(kind=8) x,y,z
    dlabel(iqm)(1:1) = 'H'
    if(iqm.lt.10) then
        write(dlabel(iqm)(2:2),'(i1)') iqm
        dlabel(iqm)(3:5) = '   '
    elseif(iqm.lt.100) then
        write(dlabel(iqm)(2:3),'(i2)') iqm
        dlabel(iqm)(4:5) = '  '
    elseif(iqm.lt.1000) then
        write(dlabel(iqm)(2:4),'(i3)') iqm
        dlabel(iqm)(5:5) = ' '
    elseif(iqm.lt.1000) then
        write(dlabel(iqm)(2:5),'(i4)') iqm
    else
        write(0,*) 'Error: more than 10000 atoms'
        stop
    endif
    write(kfile1,'(a,2x,3f12.5)') dlabel(iqm)(1:1),x,y,z
    modum = modum+1
    if(modum.lt.10) then
        write(kfile2,'(a,i5,2x,a,i1,a,3f8.3,f8.4,f8.3)') 'ATOM  ', &
            iqm,'H',modum,'  MOD  9999    ',x,y,z,0.0,1.2
    else
        write(kfile2,'(a,i5,2x,a,i2,a,3f8.3,f8.4,f8.3)') 'ATOM  ', &
            iqm,'H',modum,' MOD  9999    ',x,y,z,0.0,1.2
    endif
    return 
end subroutine addH

subroutine getcharge(jstart,jfinal,charge1,tcharge1)
    implicit none
    integer jstart,jfinal
    real(kind=8)::charge1(10000),tcharge1
    integer::i
    tcharge1 = 0.0d0
    do i = jstart, jfinal
        tcharge1 = tcharge1 + charge1(i)
    enddo
    return
end subroutine getcharge
