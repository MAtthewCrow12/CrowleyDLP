
program LiquidDM
implicit none
integer              ::  ZDnrange, NDzrange
integer, allocatable :: zdripndep(:)
integer, allocatable :: ndripzdep(:)



!write(*, "(A)", advance="NO") "desired n range: "
!read *, ZDnrange
ZDnrange=177


allocate (zdripndep(ZDnrange))
call ProtonDripLine(ZDnrange, zdripndep)
!NDzrange=zdripndep(ZDnrange)
!allocate (ndripzdep(NDzrange))
!call NeutronDripLine(NDzrange, ndripzdep)
!call ProntonDripSeparationBindingEnergy(ZDnrange, zdripndep)
!call BindingEnergy(ZDnrange)


end program LiquidDM


!    |||//   ||||||\
!      //    ||    ||
!     //     ||    ||
!    //|||   ||||||/  
subroutine ProtonDripLine(ZDnrange, zdripndep)
implicit none
integer :: n, z
integer :: ZDnrange
integer :: initialz, finalz, ZDmaxzrange
integer :: zdripndep(ZDnrange)
real    :: befirst, besecond
real    :: sepE 


open(unit=20,  file="Zdrip.dat",                status="unknown")
open(unit=21,  file="Zdripsearch.dat",          status="unknown")


initialz=-5
finalz=20
ZDmaxzrange=4*ZDnrange
do n=2, ZDnrange

    do z=initialz, finalz
        call BindingE(z, n, befirst)
        call BindingE(z+1, n, besecond)
        sepE=besecond-befirst
        write(21, *) n, z
    
        print *, 'n: ', n, 'z:', z
        !print *, 'initialz:', initialz
        !print *, 'finalz:', finalz
        print *, 'first: ', befirst
        print *, 'second: ', besecond
        print *, 'separation: ', sepE
        print *, ' '
        

        if (sepE<0 .and. besecond>0 .and. befirst>0)then
            print *, 'IF STATEMENT ACTIVATED !!!!!!!!!!!'
            print *, ' '
            print *, ' '
            zdripndep(n)=z
            initialz=z-10
            finalz=z+10
            exit
        !else
            !zdripndep(n)=0
        end if
    end do

    print *, ' '
end do


do n=2, ZDnrange
    !print *, "n=", n, "zdripndep(n)=", zdripndep(n)
    write(20, *) n, zdripndep(n)
end do


end subroutine ProtonDripLine


!  |\ || ||||||\
!  |\\|| ||    ||
!  ||\\| ||    ||
!  || \| ||||||/    
subroutine NeutronDripLine(NDzrange, ndripzdep)
integer :: n, z
integer :: NDzrange
integer :: ndripzdep(NDzrange)
integer :: initaln, finaln
real    :: befirst, besecond, NDmaznrange
real    :: sepE


open(unit=22,  file="Ndrip.dat",                status="unknown")
open(unit=23,  file="Ndripsearch.dat",          status="unknown")


initialn=-5
finaln=20
NDmaznrange=4*NDzrange
do z=1, NDzrange
    do n=initialn, finaln
        call BindingE(z, n, befirst)
        call BindingE(z, n+1, besecond)
        sepE=besecond-befirst
        write(23, *) n, z
        if (sepE<0 .and. besecond>0 .and. befirst>0)then
            ndripzdep(z)=n
            initialn=n-10
            finaln=n+10
            exit
        else
            ndripzdep(z)=0
        end if
    end do
end do


do z=1, NDzrange
    !print *, "z=", z, "ndripzdep(z)=", ndripzdep(z)
    write(22, *)  ndripzdep(z), z
end do


end subroutine NeutronDripLine


! |||||\\    |||||||||
! ||||||\\   ||
! ||  |||||  ||
! ||||||//   |||||||||
! ||||||\\   |||||||||
! ||  |||||  ||
! ||||||//   ||
! |||||//    |||||||||
subroutine BindingEnergy(ZDnrange)
implicit none 
integer :: i, n, z, a, ZDnrange, NDzrange, lastn, lastz
integer, allocatable :: zdripndep(:)
integer, allocatable :: ndripzdep(:)
real,    allocatable :: theBindingEnergy(:,:)
real,    allocatable :: expBindingEnergy(:,:)
integer              :: lowerbound(177)
integer              :: upperbound(177)
real                 :: expBindEnergy
integer              :: dummy1
character            :: dummy2

real                 :: diff

!open statements
open(unit=24,  file="theBindingEnergy.dat",     status="unknown")
open(unit=25,  file="EXPERIMENT_AME2016.dat",   status="unknown")
open(unit=26,  file="expBindingEnergy.dat",     status="unknown")
open(unit=27,  file="expupperdripline.dat",     status="unknown")
open(unit=28,  file="explowerdripline.dat",     status="unknown")
open(unit=29,  file="expthediff.dat",           status="unknown")


!Finding theoretical Proton/Neutron drip line
allocate (zdripndep(ZDnrange))
call ProtonDripLine(ZDnrange, zdripndep)
NDzrange=zdripndep(ZDnrange)
allocate (ndripzdep(NDzrange))
call NeutronDripLine(NDzrange, ndripzdep)
!allocate (theBindingEnergy(NDzrange, ZDnrange))
allocate (theBindingEnergy(125, 180))
allocate (expBindingEnergy(125, 180))


!creating matrix of Theoretical Bnding Energies WITHIN THEORETICAL BINDING ENERGY DRIP LINES
!do n=1, ZDnrange  
   !do z=1, NDzrange
       !if(z<=zdripndep(n) .and. n<=ndripzdep(z)) then
           !call BindingE(z, n, theBindingEnergy(z,n))
       !end if
   !end do
!end do


!creating matriz of all Experimental Binding Energies
Do i=1, 3433
    read(25, *) dummy1, dummy2, a, n, z, expBindEnergy
    expBindingEnergy(z,n)=expBindEnergy
    !call BindingE(z, n, theBindingEnergy(z,n))
    !write(24,*) n, theBindingEnergy(z,n)
    !write(26,*) n, -expBindingEnergy(z,n)
    !write(29,*) n, theBindingEnergy(z,n)+expBindingEnergy(z,n)
end do


!finding the borders of the of the experimental Binding Energies Data
do n=1, 177
    do z=1, 120
        if(expBindingEnergy(z-1,n)==0 .and. expBindingEnergy(z,n) /=0)then
            lowerbound(n)=z
            write(28,*) n, lowerbound(n)
        end if 
        if(expBindingEnergy(z-1,n) /=0 .and. expBindingEnergy(z,n)==0)then
            upperbound(n)=z-1
            write(27,*) n, upperbound(n)
        end if
    end do
end do

do n=1, 177
    do z=lowerbound(n), upperbound(n)
        if (z<93 .and. z>19) then 
        call BindingE(z, n, theBindingEnergy(z,n))
        write(24,*) n, theBindingEnergy(z,n)
        write(26,*) n, -expBindingEnergy(z,n)
        write(29,*) n, z, theBindingEnergy(z,n)+expBindingEnergy(z,n)
        end if
    end do
end do




end subroutine BindingEnergy




























!}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}} currently this subroutine calculates the Proton Drip line as I used to.
!}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}} I have since modified and placed in the regular proton drip line 
!}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
subroutine testProtonDripLine(ZDnrange, zdripndep)
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

end subroutine testProtonDripLine
!}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
!}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
!}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
subroutine testNeutronDripLine(NDzrange, ndripzdep)
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

end subroutine testNeutronDripLine
!}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
!}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
!}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}










!subroutine to find the binding energy
subroutine BindingE(z,n,be)
implicit none
real:: be, p1, p2, p3, p4, p5
integer:: z, n, a
a=z+n
!-----------------------------
!----------LQDM---------------
!-----------------------------
p1 =15.76*(A)
p2=  -17.81*(A**(2./3.))
!p3= -(0.711*(Z*(Z-1.))*(A**(-1./3.)))
!!!!!
p3=- 0.711 * (z**2.) * a**(-1./3.)     !BETTER          !-(0.711*(Z*(Z-1.))*(A**(-1./3.)))
!!!!!

p4= -(23.702*((Z-N)**2.))/A


!p1 =15.5*(A)
!p2=  -16.8*(A**(2./3.))
!p3=-(0.72*(Z*(Z-1.))*(A**(-1./3.)))
!p4= -(23.*((Z-N)**2.))/A
!||||||||||||||||||||||||||||||||||
!if(    (mod(z,2)==0)  ) then
    !if(    (mod(n,2)==0)  ) then
    !    p5=34.*(a**(-3./4.))
        !print *, "p5 z,n= even"
    !end if
    !if(    (mod(n,2)==1)  ) then
     !   p5=0
        !a is odd
    !end if
!end if !main if

!if(    (mod(z,2)==1)  ) then
!    if(    (mod(n,2)==1)  ) then
!        p5=-34.*(a**(-3./4.))
!        !print *, "p5 z,n= odd"
!    end if
!    if(    (mod(n,2)==0)  ) then
!        p5=0
        !a is odd
!    end if
!end if 
if(    (mod(z,2)==0) .and.  (mod(n,2)==0)  ) then
    
    p5=34.*(a**(-3./4.))
    
else if (    (mod(z,2)==1) .and.  (mod(n,2)==1)  ) then
        
    p5=-34.*(a**(-3./4.))

else

    p5=0 

end if




!-----------------------------
BE= p1+p2+p3+p4+p5 !creating array of binding energy dependent on N
!-----------------------------
!-----------------------------
!-----------------------------
end subroutine bindingE













