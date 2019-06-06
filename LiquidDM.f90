program LiquidDM
implicit none



real::    sen, sez           
integer:: z, n, a, i, temp
integer:: nfinn, zfinn, nfinz, zfinz

real,    allocatable :: ben(:), bez(:)
integer, allocatable :: ndripzdep(:)
integer, allocatable :: ndripndep(:)   !array of z values for neutron drip the first neutron dependent the second proton dependent
integer, allocatable :: zdripndep(:)
real,    allocatable :: BindEM(:,:)

open(unit=32, file="ProtonDripline.dat", status="unknown")
open(unit=33, file="NeutronDripline.dat", status="unknown")
open(unit=36, file="test.dat", status="unknown")


!/////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////
write(*, "(A)", advance="NO") "desired n range for protron drip: nfinz="      !outer loop limit value should be smaller than the inner loop limit value 
read *, nfinz                        
zfinz=4*nfinz

write(*, "(A)", advance="NO") "desired z range for neutron drip: zfinn="      !outer loop limit value should be smaller than the inner loop limit value      
read *, zfinn 
nfinn=4*zfinn                   
!/////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////
allocate(bez(zfinz))
allocate(zdripndep(nfinz))

allocate(ben(nfinn))
allocate(ndripzdep(zfinn))

allocate(ndripndep(zfinn))


!=================================================================
!===========SEC 1 Proton drip line================================
!=================================================================
nlpA: do n=1, nfinz


zlpB: do z=zdripndep(n-1), zfinz
call BindingE(z,n,bez(z))
end do zlpB

!finding the separation energy: the difference between binding energies incremented z 
sepEzB: do z=zdripndep(n-1), zfinz

sez=bez(z+1)-bez(z)
if (sez<0.) then 
!}}}}}}}}}}}}}}}
zdripndep(n)=z!}  !this array is the proton drip line z values dependent on n
!}}}}}}}}}}}}}}}
goto 111       
end if 

end do sepEzB  
111 continue 


end do nlpA
!=================================================================
!=========SEC 3 Neutron drip line=================================
!=================================================================    
zloopA: do z=1, zfinn


nloopB: do n=ndripzdep(z-1), nfinn
call BindingE(z, n, ben(n))
end do nloopB
        
!finding the separation energy: the difference between binding energies incremented n
sepEnB: do n=ndripzdep(z-1), nfinn

sen=ben(n+1)-ben(n)
if (sen<0.) then 
!}}}}}}}}}}}}}}}}
ndripzdep(z)=n!}  this array is the NEUTRON DRIP LINE N VALUES DEPENDENT ON Z 
!}}}}}}}}}}}}}}}} 
goto 333
end if

end do sepEnB   
333 continue


end do zloopA
!=================================================================
!========SEC 5 printing/writing===================================
!=================================================================
prn1A: do n=1, nfinz !SEC 1 
!||||||||||||||||||||||||||||||||
!print *,  "n=", n, "z=", zdripndep(n),  "proton drip line value"
!printing the proton drip line z values dependent on n
write(32,*) n, zdripndep(n)
!||||||||||||||||||||||||||||||||
end do prn1A

prn4A: do z=1, zfinn  !SEC 3
!|||||||||||||||||||||||||||||||||||
!print *, "z=", z, "n=", ndripzdep(z), "Neutron drip line value"
!printing the Neutron drip line n values dependent on z
write(33,*) ndripzdep(z), z
!|||||||||||||||||||||||||||||
end do prn4A 

!=================================================================
!=================================================================
!===================finding binding energy========================
!=================================================================
!=================================================================

allocate(BindEM(zfinn,ndripzdep(zfinn)))

Do A=2, 100
do z=1, A-1
N=A-z
print *, ' '
if(z<=zdripndep(n)) then
if(n<=ndripzdep(z)) then
write(36, *) n, z
call BindingE(z, n, BindEM(z,n))
print *, a, z, n, BindEM(z,n)
end if
end if

end do
end do



!do n=2, ndripzdep(zfinn)


!do z=1, zdripndep(n)
!if(n==ndripzdep(z))then
!ndripndep(n)=z
!goto 345
!end if

!end do 
!345 continue
!end do 


!do n=2, ndripzdep(zfinn) 
!i=1
!do while (ndripndep(n)==0)
!ndripndep(n)=ndripndep(n+i)
!i=i+1
!end do
!write(36,*) n, ndripndep(n)
!print *, 'n=', n, 'ndripndep(n)=', ndripndep(n), 'zdripndep(n)=', zdripndep(n)
!end do

!do n=2, ndripzdep(zfinn)
!do z=ndripndep(n), zdripndep(n)
!write(36,*) n, z
!end do
!end do
!=================================================================
!=================================================================
!=================================================================
!=================================================================


end program LiquidDM
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
