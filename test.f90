program test
implicit none
integer :: n, z, ZDnrange, NDzrange
integer, allocatable :: zdripndep(:)
integer, allocatable :: ndripzdep(:)
real,    allocatable :: BindingEarr(:,:)
open(unit=20,  file="Zdrip.dat",                status="unknown")
open(unit=21,  file="Zdripsearch.dat",          status="unknown")
open(unit=22,  file="Ndrip.dat",                status="unknown")
open(unit=23,  file="Ndripsearch.dat",          status="unknown")
open(unit=24,  file="BindingEnergy.dat",        status="unknown")

!Test files
open(unit=25,  file="BETEST.dat",               status="unknown")
open(unit=26,  file="SEPTEST.dat",              status="unknown")



write(*, "(A)", advance="NO") "desired n range for protron drip: nfinz=" 
read *, ZDnrange



allocate (zdripndep(ZDnrange))
!call ProtonDripLine(ZDnrange, zdripndep)
call testProtonDripLine(ZDnrange, zdripndep)
!call ProntonDripSeparationBindingEnergy(ZDnrange, zdripndep)


!NDzrange=zdripndep(ZDnrange)
!allocate (ndripzdep(NDzrange))
!call NeutronDripLine(NDzrange, ndripzdep)


!call BindingEnergy(ZDnrange)


end program test

subroutine ProntonDripSeparationBindingEnergy(ZDnrange, zdripndep)
implicit none
integer :: n, z
integer :: ZDnrange, ZDmaxzrange
integer :: zdripndep(ZDnrange)
real    :: bez(ZDnrange)
real    :: sez
ZDmaxzrange=ZDnrange

do n=1, ZDnrange
    do z=1, ZDmaxzrange

        call BindingE(z, n, bez(z))
        call BindingE(z+1, n, bez(z+1))
        write(25, *) n, bez(z) 
        if(bez(z+1)>0 .and. bez(z)>0) then    
        sez=bez(z+1)-bez(z)
        if(sez>0)then
        write(26, *) n, sez
        end if
        end if

    end do 
end do
end subroutine ProntonDripSeparationBindingEnergy





subroutine BindingEnergy(ZDnrange)
implicit none 
integer :: n, z, ZDnrange, NDzrange
integer, allocatable :: zdripndep(:)
integer, allocatable :: ndripzdep(:)
real,    allocatable :: BindingEarr(:,:)

allocate (zdripndep(ZDnrange))
call ProtonDripLine(ZDnrange, zdripndep)
NDzrange=zdripndep(ZDnrange)
allocate (ndripzdep(NDzrange))
call NeutronDripLine(NDzrange, ndripzdep)
allocate (BindingEarr(NDzrange, ZDnrange))

do n=1, ZDnrange
   do z=1, NDzrange
       if(z<=zdripndep(n) .and. n<=ndripzdep(z)) then
       write(24, *) n, z !Binding Energy test
        call BindingE(z, n, BindingEarr(z,n))
        end if
    end do
end do

do n=1, ZDnrange
    do z=1, NDzrange

    if(z<=zdripndep(n) .and. n<=ndripzdep(z)) then
    print *, "n=", n, "z=", z, "BE=", BindingEarr(z, n)
    end if

   end do
end do


print *, "ZDnrange=", ZDnrange, "NDzrange=", NDzrange


end subroutine BindingEnergy






subroutine ProtonDripLine(ZDnrange, zdripndep)
implicit none
integer :: n, z
integer :: ZDnrange, ZDzrange
integer :: zdripndep(ZDnrange)
real    :: bez(4*ZDnrange)
real    :: sez


print *, zdripndep


do n=1, ZDnrange

    do z=1, zdripndep(n-1)+10
        call BindingE(z,n,bez(z))
        write(21, *) n, z
    end do


    do z=1, zdripndep(n-1)+10
        sez=bez(z+1)-bez(z)
        if (sez<0.) then 
        zdripndep(n)=z
        exit       
        end if 
    end do  

end do
do n=1, ZDnrange !SEC 1 
    write(20,*) n, zdripndep(n)
    print *,    "Proton Drip:   ", "n=", n, " z=", zdripndep(n)
end do 

end subroutine ProtonDripLine




subroutine testProtonDripLine(ZDnrange, zdripndep)
implicit none
integer :: n, z
integer :: ZDnrange
integer :: initialz, finalz, ZDmaxzrange
integer :: zdripndep(ZDnrange)
real    :: befirst, besecond
real    :: sepE 

initialz=0
finalz=20
ZDmaxzrange=4*ZDnrange

do n=1, ZDnrange
do z=initialz, ZDmaxzrange

call BindingE(z, n, befirst)
call BindingE(z+1, n, besecond)
sepE=besecond-befirst

write(21, *) n, z


!print *, "n=", n, " z=", z, " sepE=", sepE
if (sepE<0 .and. besecond>0 .and. sepE<0)then
zdripndep(n)=z
initialz=z-5
finalz=z+5
exit
end if

end do
end do

do n=1, ZDnrange
print *, n, zdripndep(n)

write(20, *) n, zdripndep(n)
end do




end subroutine
















subroutine NeutronDripLine(NDzrange, ndripzdep)
implicit none
integer :: n, z
integer :: NDnrange, NDzrange
integer :: ndripzdep(NDzrange)
real    :: ben(4*NDzrange)
real    :: sen

do z=1, NDzrange

    do n=1, ndripzdep(z-1)+10
        call BindingE(z,n,ben(n))
        write(23, *) n, z
    end do


    do n=1,  ndripzdep(z-1)+10
        sen=ben(n+1)-ben(n)
        if (sen<0.) then 
        ndripzdep(z)=n
        exit       
        end if 
    end do  

end do

do z=1, NDzrange !SEC 1 
    write(22,*) ndripzdep(z), z
    print *,    "Neutron Drip: ", " z=", z, " n=", ndripzdep(z)
end do 

end subroutine NeutronDripLine










!subroutine to find the binding energy
subroutine BindingE(z,n,be)
implicit none
real:: be, p1, p2, p3, p4, p5
integer:: z, n, a
a=z+n
!-----------------------------
!----------LQDM---------------
!-----------------------------
p1 =15.5*(A)
p2=  -16.8*(A**(2./3.))
p3=-(0.72*(Z*(Z-1.))*(A**(-1./3.)))
p4= -(23.*((Z-N)**2.))/A
!||||||||||||||||||||||||||||||||||
if(    (mod(z,2)==0)  ) then
    if(    (mod(n,2)==0)  ) then
        p5=34.*(a**(-3./4.))
        !print *, "p5 z,n= even"
    end if
    if(    (mod(n,2)==1)  ) then
        p5=0
        !a is odd
    end if
end if !main if

if(    (mod(z,2)==1)  ) then
    if(    (mod(n,2)==1)  ) then
        p5=-34.*(a**(-3./4.))
        !print *, "p5 z,n= odd"
    end if
    if(    (mod(n,2)==0)  ) then
        p5=0
        !a is odd
    end if
end if !main if
!p5=0 !ignoring P5
!-----------------------------
BE= p1+p2+p3+p4+p5 !creating array of binding energy dependent on N
!-----------------------------
!-----------------------------
!-----------------------------
end subroutine bindingE
