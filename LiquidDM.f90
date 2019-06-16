program LiquidDM
implicit none

real::    sen, sez           
integer:: z, n, a, i, temp, j=0,f=0, g=0
integer:: nfinn, zfinn, nfinz, zfinz

real,    allocatable :: ben(:), bez(:)
!integer, allocatable :: ndripzdep(:)
!integer, allocatable :: ndripndep(:)   !array of z values for neutron drip the first neutron dependent the second proton dependent
integer, allocatable :: zdripndep(:)
!real,    allocatable :: BindEM(:,:)

open(unit=32,  file="ProtonDripline.dat",   status="unknown")
!open(unit=33,  file="NeutronDripline.dat",  status="unknown")

open(unit=43,  file="ZDripTest.dat",        status="unknown")
!open(unit=44,  file="NDripTest.dat",        status="unknown")
!open(unit=36,  file="BEMattest.dat",        status="unknown")


!/////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////
write(*, "(A)", advance="NO") "desired n range for protron drip: nfinz="      !outer loop limit value should be smaller than the inner loop limit value 
read *, nfinz                        
!zfinz=4*nfinz
allocate(zdripndep(nfinz))

call ProtonDripline(nfinz, zdripndep)                  
!/////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////


end program LiquidDM

!=================================================================
!===========SEC 1 Proton drip line================================
!=================================================================
subroutine ProtonDripline(nfinz, zdripndep)
implicit none 
integer::zfinz, nfinz, n, z, f=0
integer, allocatable :: zdripndep(:)
real,    allocatable :: bez(:)
real :: sez

zfinz=4*nfinz


allocate(bez(zfinz))



nlpA: do n=1, nfinz


zlpB: do z=zdripndep(n-1), zdripndep(n-1)+5
call BindingE(z,n,bez(z))
f=f+1
write(43, *) n, z
end do zlpB

sepEzB: do z=zdripndep(n-1),   zfinz

sez=bez(z+1)-bez(z)
if (sez<0.) then 

zdripndep(n)=z
goto 111       
end if 

end do sepEzB  
111 continue 


end do nlpA
!=================================================================
!=========SEC 3 Neutron drip line=================================
!=================================================================    

!zfinn=zdripndep(nfinz)
!allocate(ndripzdep(zfinn))
!nfinn=4*zfinn                   
!allocate(ben(nfinn))

!zloopA: do z=1, zfinn 


!nloopB: do n=ndripzdep(z-1), ndripzdep(z-1)+5    
!call BindingE(z, n, ben(n))
!write(44, *) n, z
!g=g+1
!end do nloopB
!sepEnB: do n=ndripzdep(z-1), nfinn

!sen=ben(n+1)-ben(n)
!if (sen<0.) then 
!ndripzdep(z)=n 
!goto 333
!end if

!end do sepEnB   
!333 continue


!end do zloopA
!=================================================================
!========SEC 5 printing/writing===================================
!=================================================================
prn1A: do n=1, nfinz !SEC 1 
write(32,*) n, zdripndep(n)
end do prn1A

!prn4A: do z=1, zfinn  !SEC 3
!write(33,*) ndripzdep(z), z
!end do prn4A 

!=================================================================
!=================================================================
!===================finding binding energies======================
!=================================================================
!=================================================================

!allocate(  BindEM( zfinn, nfinz )  )

!do z=1, zfinn
!do n=1, nfinz !ndripzdep(zfinn)

!if(z<=zdripndep(n) .and. n<=ndripzdep(z)) then
!write(36, *) n, z !Binding Energy test
!call BindingE(z, n, BindEM(z,n))
!j=j+1
!end if
!end do
!end do

!write (*, 32) zfinn, nfinz

!32 format ( /, "binding energy values have been found in a ", I3, " by ", I3, " matrix." /, /)

write (*, 43) f
!write (*, 44) g
!write (*, 45) j
!write (*, 46) j+g+f
43 format ("Binding Energy function called: ", I9, " times by Proton Drip line." /)
!44 format ("Binding Energy function called: ", I9, " times by Neutron Drip Line."  /)
!45 format ("Binding Energy function called: ", I9, " times by Binding Energy Matrix." /)
!46 format ("Total Binding Energy calls: ",     I9 /)
!=================================================================
!=================================================================
end subroutine ProtonDripline





!-------------------------------------
!-------------------------------------
!-------------------------------------
!-------------------------------------
!-------------------------------------
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
