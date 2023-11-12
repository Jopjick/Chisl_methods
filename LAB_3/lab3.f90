module funcs
implicit none

contains

function norm_v(a, k) result(n)
implicit none
real(8):: a(20)
real(8) n
integer i
integer:: k
n=0d0
do i=1,k
n = n + a(i)**2
enddo
n = n**0.5
end function norm_v

function norm1_v(x, n) result(norm)
implicit none
real(8):: x(20)
integer:: n
integer i
real(8) norm

norm=0d0
do i=1,n
norm = norm + dabs(x(i))
enddo

end function norm1_v

function norm1_ar(A, n) result(norm)
implicit none
real(8):: A(20,20), temp_norm(20)
real(8) temp
integer:: n
integer i, j
real(8) norm

temp=0d0
do i=1,n
do j=1,n
temp = temp + dabs(a(i,j))
enddo
temp_norm(i) = temp
temp = 0d0
enddo
norm = maxval(temp_norm)

end function norm1_ar

function Read_ar(k1,k) result(a)
implicit none
integer, intent(in):: k, k1
integer i, j
real(8):: a(20,20)
do i=1,k
read(k1, *) (a(i,j), j=1,k)
enddo

end function Read_ar

function multiply_ar_ar(ar1,ar2, n) result(res)
implicit none
real(8):: ar1(20,20), ar2(20,20), res(20,20)
integer:: n
integer i, j, k
real(8) temp

temp = 0d0
do i=1, n
do j=1, n

do k =1, n
temp = temp + ar1(i,k)*ar2(k,j)
enddo
res(i,j) = temp
temp = 0d0
enddo
enddo

end function multiply_ar_ar

function multiply_v_ar(v,ar, n) result(res)
implicit none
real(8):: v(20), ar(20,20), res(20)
integer:: n
real(8) temp
integer i, j

temp = 0d0
do i=1,n
do j=1,n
temp = temp + v(j)*ar(j,i)
enddo
res(i) = temp
temp = 0d0
enddo

end function multiply_v_ar

function multiply_v_scalar(v1,v2, n) result(res)
implicit none
real(8):: v1(20), v2(20), res
integer:: n
integer i
real(8) temp

temp = 0d0
do i = 1, n
temp = temp + v1(i)*v2(i)
enddo
res = temp

end function multiply_v_scalar

function multiply_v_v(v1,v2, n) result(res)
implicit none
real(8):: v1(20), v2(20), res(20,20)
integer:: n
integer i, j

do i = 1, n
do j = 1, n
res(i,j) = v1(i)*v2(j)
enddo
enddo


end function multiply_v_v

function make_b(A, x, k) result(b)
implicit none
real(8):: A(20,20), x(20), b(20)
integer i, j
real(8) temp
integer , intent(in):: k

temp = 0d0
do i=1,k
do j=1,k
temp = temp + A(i,j)*x(j)
enddo
b(i) = temp
temp = 0d0
enddo

end function make_b

function ar_reflection(x,y,n) result(P)
implicit none
real(8):: w(20), P(20,20), x(20), y(20), e(20,20)
integer:: n
integer i
real(8) b 

w = x - y
P = multiply_v_v(w,w,n)

do i=1,n
e(i,i) = 1d0
enddo

b = 2d0/norm_v(w,n)**2

P = E - b*P


end function ar_reflection

function make_r(A, b, n) result(res)
implicit none
real(8):: A(20,20), res(20,20), b(20)
integer:: n
real(8):: x(20), y(20), P(20,20), P_temp(20,20)
integer k, i, j

res = A

do k = 0, n-1

do i=1,n
x(i)=0d0
y(i)=0d0
do j=1,n
P(i,j)=0d0
enddo
enddo

do i=1, k

P(i,i) = 1d0
enddo


do i=k+1,n
x(i - k) = res(i,k + 1)
enddo

if (x(1).EQ.0d0) then
y(1) = norm_v(x,n)
else
y(1) = -x(1)/dabs(x(1))*norm_v(x,n)
endif

P_temp = ar_reflection(x, y, n-k)
do i=k + 1,n
do j=k + 1,n
P(i,j) = P_temp(i - k,j - k)

enddo
enddo

res = multiply_ar_ar(P, res, n)
b = make_b(P,b, n)
enddo


end function make_r

function solve_r(R, b, n) result(x)
implicit none
real(8):: R(20,20), b(20), x(20)
integer:: n
integer i, j
real(8) temp
temp = 0d0

do i=1,n
x(i)=0d0
enddo

do i =n, 1, -1

do j= n, i+1, -1
temp = temp + x(j)*R(i,j)
enddo

x(i) = (b(i) - temp)/R(i,i)
temp = 0d0
enddo

end function solve_r

subroutine print_ar(A, n)
implicit none
real(8):: A(20,20)
integer:: n
integer i

do i =1, n
write(*,*) A(i,1:n)
enddo

end subroutine print_ar

end module

program lab_3
use funcs
implicit none
real(8):: U(20,20), A(20,20), Q(20,20), Q_T(20,20), R(20,20)
real(8):: x_i(20), b_i(20), b(20), x(20), b_j(20), d_A(20,20), A1(20,20), R1(20,20), x1(20), b1(20), d_a_new(20,20)
real(8):: H(20,20), C(20,20)
integer i, n, j, k
real(8) cond, dx,dA, nA, nd
n=12

open(100, file='A.txt', status='old')
open(101, file='Ort_T.txt', status='old')
open(102, file='Ort.txt', status='old')
U = Read_ar(100,n)
Q_T = Read_ar(101,n)
Q = Read_ar(102,n)
close(100)
close(101)
close(102)
open(100, file='results.txt', status='replace')

open(101, file='d_A.txt', status='old')
d_A = Read_ar(101, n)
close(101)

do j=1, 20
cond =5d0**j
U(1,1) = cond
A = multiply_ar_ar(Q_T,U,n)
A = multiply_ar_ar(A,Q,n)
A1 = A+d_A

do i=1,n
x_i(i)=i*(-1d0)**i
enddo

b_i = make_b(A, x_i, n)
b = b_i
b1 = b_i
R = make_r(A, b, n)
x = solve_r(R,b,n)

R1 = make_r(A1, b1, n)
x1 = solve_r(R1,b1,n)

b_j = make_b(R, x_i, n)

dx = norm1_v(x1-x, n)/norm1_v(x, n)
nA = norm1_ar(A, n)
nd = norm1_ar(d_A, n)
dA = norm1_ar(d_A, n)/nA

write(100,*) cond, norm_v(x-x_i, n), norm_v(b-b_j, n), dx/dA, cond/(1-cond/norm1_ar(A, n)*norm1_ar(d_A, n))

enddo
close(100)

do k=1, 1
open(110, file='pogr1.txt', status='replace')
cond = 200d0
U(1,1) = cond
A = multiply_ar_ar(Q_T,U,n)
A = multiply_ar_ar(A,Q,n)

do j=1,10
d_A_new = (1.5d0**j)*d_A
A1 = A + d_A_new

b_i = make_b(A, x_i, n)
b = b_i
b1 = b_i

R = make_r(A, b, n)
x = solve_r(R,b,n)

R1 = make_r(A1, b1, n)
x1 = solve_r(R1,b1,n)

write(110,*) norm1_v(x-x1, n)/norm1_v(x, n), norm1_ar(d_A_new, n)/norm1_ar(A, n), cond/(1d0-cond/norm1_ar(A, n)*norm1_ar(d_A_new, n))
enddo
close(110)


open(110, file='pogr2.txt', status='replace')
cond = 100000000d0
U(1,1) = cond
A = multiply_ar_ar(Q_T,U,n)
A = multiply_ar_ar(A,Q,n)

do j=1,10

d_A_new = (1.5d0**j)*d_A
A1 = A + d_A_new

b_i = make_b(A, x_i, n)
b = b_i
b1 = b_i

R = make_r(A, b, n)
x = solve_r(R,b,n)

R1 = make_r(A1, b1, n)
x1 = solve_r(R1,b1,n)

write(110,*) norm1_v(x-x1, n)/norm1_v(x, n), norm1_ar(d_A_new, n)/norm1_ar(A, n), cond/(1d0-cond/norm1_ar(A, n)*norm1_ar(d_A_new, n))
enddo
close(110)
enddo


open(105, file='hilbert.txt', status='replace')
do k=2,12


do i=1,20
x(i) = 0d0
enddo

n = k
do i=1,k
x_i(i)=i*(-1d0)**i
do j=1,k
H(i,j) = 1d0/(i+j-1)
R(i,j) = 0d0
R1(i,j) = 0d0
enddo
enddo

b_i = make_b(H, x_i, n)
b = b_i
b1 = b_i
R = make_r(H, b, n)
x = solve_r(R,b,n)

d_A = (0.01d0)*d_A
A1=H+d_A

R1 = make_r(A1, b1, n)
x1 = solve_r(R1,b1,n)

b_j = make_b(R, x_i, n)

dx = norm1_v(x1-x, n)/norm1_v(x, n)
nA = norm1_ar(A, n)
nd = norm1_ar(d_A, n)
dA = norm1_ar(d_A, n)/nA

write(105,*) k, norm_v(x-x_i, n), norm_v(b-b_j, n), norm1_v(x1-x, n)/norm1_v(x, n)/dA
enddo
close(105)


U(1,1) = 0
A = multiply_ar_ar(Q_T,U,n)
A = multiply_ar_ar(A,Q,n)
C = A


b_i = make_b(C, x_i, n)
b = b_i
b1 = b_i
R = make_r(C, b, n)
x = solve_r(R,b,n)
A1 = C + d_A
R1 = make_r(A1, b1, n)
x1 = solve_r(R1,b1,n)

b_j = make_b(R, x_i, n)

dx = norm1_v(x1-x, n)/norm1_v(x, n)
nA = norm1_ar(C, n)
nd = norm1_ar(d_A, n)
dA = norm1_ar(d_A, n)/nA
write(*,'(f8.3 f8.3 f8.3 f8.3 f8.3 f8.3 f8.3 f8.3 f8.3 f8.3 f8.3 f8.3)') x(1:12)
write(*,*) '      ||Xi-X||                      ||AX-B||'
write(*,*) norm_v(x-x_i, n), norm_v(b-b_j, n)
write(*,*)

C = 10d0*A

b_i = make_b(C, x_i, n)
b = b_i
b1 = b_i
R = make_r(C, b, n)
x = solve_r(R,b,n)
A1 = C + d_A
R1 = make_r(A1, b1, n)
x1 = solve_r(R1,b1,n)

b_j = make_b(R, x_i, n)

dx = norm1_v(x1-x, n)/norm1_v(x, n)
nA = norm1_ar(C, n)
nd = norm1_ar(d_A, n)
dA = norm1_ar(d_A, n)/nA
write(*,'(f8.3 f8.3 f8.3 f8.3 f8.3 f8.3 f8.3 f8.3 f8.3 f8.3 f8.3 f8.3)') x(1:12)
write(*,*) '      ||Xi-X||                      ||AX-B||'
write(*,*) norm_v(x-x_i, n), norm_v(b-b_j, n)

end program lab_3
