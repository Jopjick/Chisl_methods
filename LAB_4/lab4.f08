module funcs

contains

function Read_ar(k1,k) result(a)
implicit none
integer, intent(in):: k, k1
integer i, j
real(8):: a(k,k)
do i=1,k
read(k1, *) (a(i,j), j=1,k)
enddo
end function Read_ar

function multiply_ar_v(ar,v, n) result(res)
implicit none
real(8), allocatable:: v(:), ar(:,:), res(:)
integer:: n
real(8) temp
integer i, j

allocate(res(n))

temp = 0d0
do i=1,n
do j=1,n
temp = temp + v(j)*ar(i,j)
enddo
res(i) = temp
temp = 0d0
enddo

end function multiply_ar_v

function vect_norm(v) result(res)
implicit none
real(8), allocatable:: v(:)
real(8) res

res = maxval(v)

end function vect_norm

function arr_norm(A, n) result(res)
implicit none
real(8), allocatable:: A1(:), A(:,:)
real(8) res, temp
integer:: n
integer i, j

allocate(A1(n))
temp = 0d0
do i=1, n
do j=1, n
temp = temp + dabs(A(i,j))
enddo
A1(i) = temp
temp = 0d0
enddo
res = maxval(A1)

end function arr_norm

function next_x(x, A, b, alpha, n) result(res)
implicit none
real(8), allocatable:: x(:), A(:,:), b(:), res(:)
real(8):: alpha
integer:: n

allocate(res(n))

res = multiply_ar_v(A, x, n)
res = res-b
res = alpha*res
res = x - res

end function next_x

function multiply_ar_ar(ar1,ar2, n) result(res)
implicit none
real(8), allocatable:: ar1(:,:), ar2(:,:), res(:,:)
integer:: n
integer i, j, k
real(8) temp

allocate(res(n,n))

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

function create_A(D, n) result(A)
implicit none
real(8), allocatable:: D(:,:), A(:,:), Q(:,:), Q_t(:,:)
integer:: n

allocate(A(n,n), Q(n,n), Q_t(n,n))

open(100, file='Ort.txt', status='old')
open(101, file='Ort_T.txt', status='old')
Q = Read_ar(100, n)
Q_t = Read_ar(101, n)
close(100)
close(101)

A = multiply_ar_ar(Q_t, D, n)
A = multiply_ar_ar(A, Q, n)

end function create_A 

function create_alpha(D, n) result(alpha)
implicit none
real(8), allocatable:: D(:,:), diag(:)
real(8):: alpha(2)
integer:: n
real(8) d1, d2
integer i

allocate(diag(n))

do i=1,n
diag(i) = D(i,i)
enddo

d1 = minval(diag)
d2 = maxval(diag)

alpha(1) = 4/(2d0*(d1+d2)-(d2-d1)*2d0**0.5)
alpha(2) = 4/(2d0*(d1+d2)+(d2-d1)*2d0**0.5)

end function create_alpha

subroutine print_ar(A, n)
implicit none
real(8), allocatable:: A(:,:)
integer:: n
integer i


do i =1, n
write(*,*) A(i,1:n)
enddo

end subroutine print_ar

end module

program lab4
use funcs
implicit none
real(8), allocatable:: D(:,:), A(:,:), x_i(:), b(:), x(:), d_x(:), n_A(:), x_prev(:), e(:,:), c(:,:)
real(8):: alpha(2), alph1(2)
real(8) eps, det
integer i, n, k, j

n= 12

allocate(D(n,n), A(n,n), x_i(n), b(n), x(n), d_x(n), n_A(n), x_prev(n), e(n,n), c(n,n))

do i=1,12
D(i,i) = i
e(i,i) = 1
enddo

alpha = create_alpha(D, n)

A = create_A(D, 12)
do i=1,n
x_i(i)=i*(-1d0)**i
enddo

b = multiply_ar_v(A, x_i, n)

C = E-alpha(1)*A
alph1(1) = arr_norm(C, n)
C = E-alpha(2)*A
alph1(2) = arr_norm(C, n)

!Ошибка и невязка от точности

open(100, file='l4_o_n.txt', status='replace')
do i=1, 14

eps =0.1d0**i

do j=1,n
x(j)= b(j)/A(j,j)
x_prev(j) = 0d0
enddo

d_x = x - x_prev
j=0
k=1
do while (vect_norm(d_x).gt.(dabs(eps*(1-alph1(k))/alph1(k))))
j = j+1
k = mod(j,2)
if (k.eq.0) then
k=2
endif

x_prev = x
x=next_x(x,A,b,alpha(k),n)
d_x = x - x_prev

enddo
n_A = multiply_ar_v(A,x,n)
n_A = n_A-b
write(100,*) eps, vect_norm(d_x), vect_norm(n_A)
enddo
close(100)

!Число итераций от определителя

open(101, file='l4_i_d.txt', status='replace')
eps = 0.1d0**14

do i=n,2,-1
D(i,i) = D(i,i)/6d0
enddo

do i =n, 1,-1

if (i.lt.n) then
D(i+1,i+1) = D(i+1,i+1)/10d0
endif

A = create_A(D,n)
alpha = create_alpha(D,n)
det = 1d0
do j=1,n
x(j)= b(j)/A(j,j)
det = det*D(j,j)
enddo

b = multiply_ar_v(A, x_i, n)

d_x = x-x_i
j=0
k=1
do while (vect_norm(d_x).gt.(dabs(eps*(1-alph1(k))/alph1(k))))
j = j+1
k = mod(j,2)
if (k.eq.0) then
k=2
endif

x=next_x(x,A,b,alpha(k),n)
d_x = x-x_i
enddo

write(*,*) j, det
enddo

close(101)
end program lab4
