program testsemf
implicit none
integer :: z, n, a
real    ::be
open(unit=20,  file="check_py_ldm.dat",                status="unknown")

do n=1, 60
    !print *, "n=", n
    print  *, "!!!!!!"
    do z=1, 60

        call BindingE(z,n,be)
        !print  *, "a=", n+z, "n=", n, "z=", z, "be=", be
        write(20, *) z, n, be

end do
    print  *, "!!!!!!"
end do




end program testsemf


!subroutine to find the binding energy
subroutine BindingE(z,n,be)
implicit none
real:: be, p1, p2, p3, p4, p5, p5test
integer:: z, n, a
a=z+n
!-----------------------------
!----------LQDM---------------
!-----------------------------
p1 =15.76*(A)
p2=  -17.81*(A**(2./3.))

p3=- 0.711 * (z**2.) * a**(-1./3.)                !-(0.711*(Z*(Z-1.))*(A**(-1./3.)))

p4= -(23.702*((Z-N)**2.))/A
!print *, p1
!print *, p2
!print *, "!!!"
!print *, p3
!print *, "!!!"
!print *, p4

!p1 =15.5*(A)
!p2=  -16.8*(A**(2./3.))
!p3=-(0.72*(Z*(Z-1.))*(A**(-1./3.)))
!p4= -(23.*((Z-N)**2.))/A
!||||||||||||||||||||||||||||||||||
!if(    (mod(z,2)==0)  ) then
!    if(    (mod(n,2)==0)  ) then
!        p5=34.*(a**(-3./4.))
!        print *, "p5 z,n= even"
!    end if
!    if(    (mod(n,2)==1)  ) then
!        p5=0
!        print*, "a is odd"
!    end if
!end if !main if

!if(    (mod(z,2)==1)  ) then
!    if(    (mod(n,2)==1)  ) then
!        p5=-34.*(a**(-3./4.))
!        print *, "p5 z,n= odd"
!    end if
!    if(    (mod(n,2)==0)  ) then
!        p5=0
!       !print *, "a is odd"
!    end if
!end if 

!print *, '!!!!!'
if(    (mod(z,2)==0) .and.  (mod(n,2)==0)  ) then
    
    p5=34.*(a**(-3./4.))
    !print *, "n/z EVEN  p5 !=0 & +"
    
else if (    (mod(z,2)==1) .and.  (mod(n,2)==1)  ) then
        
    p5=-34.*(a**(-3./4.))
    !print *, "n/z ODD  p5 !=0 & -"

else

    p5=0 
    !print *, "A odd p5=0"

end if


!-----------------------------
BE= p1+p2+p3+p4+p5 !creating array of binding energy dependent on N
!-----------------------------
!-----------------------------
!-----------------------------
!print *, 'n=', n
!print *, 'z=', z
!print *, 'a=',  a
!print *, 'p1=',  p1
!print *, 'p2=',  p2
!print *, 'p3=',  p3
!print *, 'p4=',  p4
!print *, 'p5=',  p5test
print *, 'N,Z', n,z, 'BE=',  BE

end subroutine bindingE
