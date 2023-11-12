function MPD_AL(eps) result(Y)
implicit none
real(8)::eps, x0=-3.668221137346 
real(8) a, b, x, f, c
real(8):: Y(2)
integer i

f(x) = 2d0*x**4 - 9d0*x**3 - 60d0*x**2+1d0
a = -4d0
b = -3d0
i = 0d0

do while (b - a > eps)
c = (a + b)/2d0
i = i + 1

if (f(c)>0) then
a = c
else
b = c
endif

end do

y(1) = dabs(c-x0)
y(2) = real(i)

end function MPD_AL

function MPD_TR(eps) result(Y)
implicit none
real(8)::eps, x0=2.5061841466899
real(8) a, b, x, f, c
real(8):: Y(2)
integer i

f(x) = x*dlog10(x) - 1
a = 1d0
b = 5d0
i = 0d0

do while (b - a > eps)
c = (a + b)/2
i = i + 1

if (f(c)>0) then
b = c
else
a = c
endif

end do

Y(1) = dabs(c-x0)
Y(2) = real(i)


end function MPD_TR

function MMN_AL(eps) result(Y)
implicit none
real(8)::eps, x0=-3.668221137346 
real(8) x, f, k1, f1, x_prev
real(8):: Y(2)
integer i

f(x) = 2d0*x**4 - 9d0*x**3 - 60d0*x**2+1d0
f1(x) = 8d0*x**3 - 27d0*x**2 - 120d0*x
x = -4d0

k1 = f1(x)

do while (dabs(x - x_prev)>1.512d0*eps)
i = i + 1
x_prev = x
x = x_prev - f(x_prev)/k1
enddo

Y(1)= dabs(x -x0)
Y(2)= real(i)

end function MMN_AL

function MMN_TR(eps) result(Y)
implicit none
real(8)::eps, x0=2.5061841466899 
real(8) x, f, k1, f1, x_prev
real(8):: Y(2)
integer i

f(x) = x*dlog10(x) - 1d0
f1(x) = (dlog10(x) + 1d0)/dlog(10d0)
x=5d0

k1 = f1(x)

do while (dabs(x - x_prev)>2d0*eps)
i = i + 1
x_prev = x
x = x_prev - f(x_prev)/k1
enddo

Y(1)= dabs(x - x0)
Y(2)= real(i)

end function MMN_TR

program test_methods
implicit none
real(8) eps
real(8):: f1(2), f2(2)
integer i

interface

function MPD_AL(eps) result(Y)
real(8):: eps, Y(2)
end function MPD_AL

function MPD_TR(eps) result(Y)
real(8):: eps, Y(2)
end function MPD_TR

function MMN_AL(eps) result(Y)
real(8):: eps, Y(2)
end function MMN_AL

function MMN_TR(eps) result(Y)
real(8):: eps, Y(2)
end function MMN_TR

end interface

open(11, file='MPD.txt', status='replace')

do i=1, 7, 1
eps = 0.1d0**(i)

f1 = MPD_AL(eps)
f2 = MPD_TR(eps)

write(11, '(f9.7 a f11.9 a i3 a f11.9 a i3)') eps, ' ', f1(1), ' ', int(f1(2)), ' ', f2(1), ' ', int(f2(2)) 

enddo

open(12, file='MMN.txt', status='replace')

do i=1, 7, 1
eps = 0.1d0**(i)
f1 = MMN_AL(eps)
f2 = MMN_TR(eps)

write(12, '(f9.7 a f11.9 a i3 a f11.9 a i3)') eps, ' ', f1(1), ' ', int(f1(2)), ' ', f2(1), ' ', int(f2(2)) 

enddo
end program test_methods
