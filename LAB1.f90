program bisection_method
implicit none
integer i
real(8) a, b, c, eps, x, f, g
real(8):: a0=-4, b0=-3, a1=1, b1=5

f(x) = 2*x**4 - 9*x**3 - 60*x**2+1
g(x) = x*log10(x) - 1
read(*,*) eps

a = a0
b = b0

write(*,'(a)') '                ALG FUNC'
write(*,'(a)') '-----------------------------------------'
write(*,'(a)') '| i |   Left    |   Right   |   Root    |'
write(*,'(a)') '-----------------------------------------'
do while (b - a > eps)
c = (a + b)/2
i = i + 1
write(*,'(a i2 a f11.7 a f11.7 a f11.7 a)') '| ', i, '|', a, '|', b, '|', c, '|'
write(*,'(a)') '-----------------------------------------'
if (f(c)>0) then
a = c
else
b = c
endif

enddo
i = 0
write(*,*)
a = a1
b = b1
write(*,'(a)') '              TRANSC FUNC'
write(*,'(a)') '-----------------------------------------'
write(*,'(a)') '| i |   Left    |   Right   |   Root    |'
write(*,'(a)') '-----------------------------------------'
do while (b - a > eps)
c = (a + b)/2
i = i + 1
write(*,'(a i2 a f11.7 a f11.7 a f11.7 a)') '| ', i, '|', a, '|', b, '|', c, '|'
write(*,'(a)') '-----------------------------------------'
if (g(c)>0) then
b = c
else
a = c
endif

enddo

end program bisection_method
