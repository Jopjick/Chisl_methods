function Read_ar(k) result(a)
implicit none
integer, intent(in):: k
integer i, j
real(8):: a(12,12)
do i=1,12
read(k, *) (a(i,j), j=1,12)
enddo

end function Read_ar

function Read_v(k) result(a)
implicit none
integer, intent(in):: k
integer i
real(8):: a(12)

read(k, *) (a(i), i=1,12)

end function Read_v

function solve_l(l, b, k) result(z)
implicit none
real(8):: l(12,12), b(12)
real(8):: z(12)
integer i, j
integer:: k
real(8) temp

do i=1,k
temp = 0d0
do j=1,i-1
temp = temp + l(i,j)*z(j)
enddo
z(i) = b(i)-temp
temp = 0d0
enddo
end function solve_l

function solve_d(d, z, k) result(y)
implicit none
real(8):: d(12,12), z(12)
real(8):: y(12)
integer i
integer:: k
do i=1,k
y(i) = z(i)/d(i,i)
enddo

end function solve_d

function solve_r(r, y, k) result(x)
implicit none
real(8):: r(12,12), y(12)
real(8):: x(12)
real(8) temp
integer i, j
integer:: k
do i=k,1,-1
temp = 0d0
do j=12,i+1,-1
temp = temp + r(i,j)*x(j)
enddo
x(i) = y(i) - temp
temp = 0d0
enddo

end function solve_r 

function norm(a, k) result(n)
implicit none
real(8):: a(12)
real(8) n
integer i
integer:: k
n=0d0
do i=1,k
n = n + a(i)**2
enddo
n = n**0.5
end function norm

function make_matrix(mat1,mat2,mat3) result(mat4)
implicit none
real(8):: mat1(12,12), mat2(12,12), mat3(12,12), mat_temp(12,12), mat4(12,12)
integer i, j ,k
real(8) temp
temp = 0d0
do i=1,12
do j=1,12

do k =1,12
temp = temp + mat1(i,k)*mat2(k,j)
enddo
mat_temp(i,j) = temp
temp = 0d0
enddo
enddo
temp = 0d0
do i=1,12
do j=1,12

do k =1,12
temp = temp + mat_temp(i,k)*mat3(k,j)
enddo
mat4(i,j) = temp
temp = 0d0
enddo
enddo
end function make_matrix

function make_b(A, x, k) result(b)
implicit none
real(8):: A(12,12), x(12), b(12)
integer i, j
real(8) temp
integer:: k

temp = 0d0
do i=1,k
do j=1,k
temp = temp + A(i,j)*x(j)
enddo
b(i) = temp
temp = 0d0
enddo

end function make_b

function solve_matrix(cond, b1, A1, x_new) result(x)
implicit none
real(8):: U(12,12), l(12,12), d(12,12), r(12,12), b(12), x_new(12)
real(8):: q(12,12), q_t(12,12), A(12,12), A1(12,12)
real(8):: x(12), y(12), z(12), x_i(12), b1(12), b_d(12), b_new(12)
integer i1, j1, i
real(8) cond

interface

function Read_ar(k) result(a)
integer, intent(in):: k
real(8):: a(12,12)
end function Read_ar

function make_matrix(q,a,q_t) result(b)
real(8):: q(12,12), q_t(12,12), a(12,12), b(12,12)
end function make_matrix

function Read_v(k) result(a)
integer, intent(in):: k
real(8):: a(12)
end function Read_v

function solve_l(l,b, k) result(z)
real(8):: l(12,12), b(12)
real(8):: z(12)
integer:: k
end function solve_l

function solve_d(d,z, k) result(y)
real(8):: d(12,12), z(12)
real(8):: y(12)
integer:: k
end function solve_d

function solve_r(r, y, k) result(x)
real(8):: r(12,12), y(12)
real(8):: x(12)
integer:: k
end function

function make_b(A, x, k) result(b)
real(8):: A(12,12), x(12), b(12)
integer:: k
end function make_b

end interface

open(101, file='Ort.txt', status='old')
open(102, file='Ort_T.txt', status='old')
open(103, file='A.txt', status='old')

q = Read_ar(101)
q_t = Read_ar(102)
U = Read_ar(103)
close(101)
close(102)
close(103)
U(1,1) = cond
A = make_matrix(q_t,U,q)
A1=A
do i1=1,12
do j1=1,12
l(i1,j1)=0
d(i1,j1)=0
r(i1,j1)=0
enddo
enddo

call factor(l,d,r,A)

do i=1,12
x_i(i) = (i*(-1)**i)*1d0
enddo

open(107, file='db.txt', status='old')
read(107,*) (b_d(i), i=1,12)
close(107)

b = make_b(A,x_i, 12)

do i=1,12
b_new(i) = b(i) + b_d(i)
enddo

z = solve_l(l,b_new, 12)
y = solve_d(d,z, 12)
x_new = solve_r(r,y, 12)

do i1=1,12
y(i1)=0
z(i1)=0

enddo

z = solve_l(l,b, 12)
y = solve_d(d,z, 12)
x = solve_r(r,y, 12)

b1 = b
end function

function hilbert(b1, A, k) result(x)
implicit none
real(8) b1(12), A(12,12), x(12)
integer:: k
real(8) l(12,12), d(12,12), r(12,12)
real(8) x_i(12), b(12), y(12), z(12)

interface

function solve_l(l,b, k) result(z)
real(8):: l(12,12), b(12)
real(8):: z(12)
integer:: k
end function solve_l

function solve_d(d,z, k) result(y)
real(8):: d(12,12), z(12)
real(8):: y(12)
integer:: k
end function solve_d

function solve_r(r, y, k) result(x)
real(8):: r(12,12), y(12)
real(8):: x(12)
integer:: k
end function

function make_b(A, x, k) result(b)
real(8):: A(12,12), x(12), b(12)
integer:: k
end function make_b

end interface


b = make_b(A, x_i, k)

call factor_h(l,d,r,A,k)

z = solve_l(l,b, k)
y = solve_d(d,z, k)
x = solve_r(r,y, k)

b1=b

end function hilbert

program lab2
implicit none
real(8):: x(12), x_i(12), x_d(12), b1(12), A1(12,12), b2(12), b_d(12), x_new(12), x_d_new(12), H(12,12)
real(8) norm
integer i, j, k

interface

function hilbert(b1, A, k) result(x)
real(8) b1(12), A(12,12), x(12)
integer:: k
end function hilbert

function solve_matrix(cond, b1, A1, x_new) result(x)
real(8) cond
real(8):: x(12), b1(12), A1(12,12), x_new(12)
end function

function Read_v(k) result(a)
integer, intent(in):: k
real(8):: a(12)
end function Read_v

function make_b(A, x, k) result(b)
real(8):: A(12,12), x(12), b(12)
integer:: k
end function make_b

end interface

do i=1,12
x_i(i) = (i*(-1)**i)*1d0
enddo

open(205, file='osh.txt', status='replace')
open(207, file='cond.txt', status='replace')
do i=1, 1000
x= solve_matrix(i*10d0, b1, A1, x_new)
b2 = make_b(A1,x, 12)
do j=1,12
b_d(j) = b2(j) - b1(j) 
x_d(j) = x(j)-x_i(j)
x_d_new(j) = x(j)-x_new(j)
enddo
write(205,*) i*10d0, norm(x_d, 12), norm(b_d, 12)
write(207,*) i*10d0, (norm(x_d, 12)/norm(x, 12))/(norm(b_d, 12)/norm(b1, 12))
enddo
close(205)
close(207)

open(111, file='hilbert.txt', status='replace')

do k=2,12

do i=2,k
x_i(i)=i*(-1d0)*i
do j=1,k
H(i,j)=1d0/(i+j-1d0)
enddo
enddo

x = hilbert(b1, H, i)
b2 = make_b(H,x, k)
do j=1,k
b_d(j) = b2(j) - b1(j) 
x_d(j) = x(j)-x_i(j)
x_d_new(j) = x(j)-x_new(j)
enddo
write(111,*) k, norm(x_d, k), norm(b_d, k)
enddo
close(111)
end program

subroutine factor(l,d,r,a)
implicit none
integer i, j, k, m
real(8):: a(12,12)
real(8):: l(12,12), d(12,12), r(12,12)
real(8) temp

do m=1,12
temp = 0d0

do k=1,m-1
temp = temp + l(m,k)*d(k,k)*r(k,m)
enddo
d(m,m) = a(m,m) - temp
temp = 0d0

do j=m,12

do k=1,m-1
temp = temp + l(m,k)*d(k,k)*r(k,j)
enddo
r(m,j) = (a(m,j)-temp)/d(m,m)
temp = 0d0
enddo

do i=m,12

do k=1,m-1
temp = temp + l(i,k)*d(k,k)*r(k,m)
enddo
l(i,m) = (a(i,m)- temp)/d(m,m)
temp = 0d0
enddo

enddo

end subroutine

subroutine factor_h(l,d,r,a,k1)
implicit none
integer i, j, k, m
integer:: k1
real(8):: a(12,12)
real(8):: l(12,12), d(12,12), r(12,12)
real(8) temp

do m=1,k1
temp = 0d0

do k=1,m-1
temp = temp + l(m,k)*d(k,k)*r(k,m)
enddo
d(m,m) = a(m,m) - temp
temp = 0d0

do j=m,12

do k=1,m-1
temp = temp + l(m,k)*d(k,k)*r(k,j)
enddo
r(m,j) = (a(m,j)-temp)/d(m,m)
temp = 0d0
enddo

do i=m,12

do k=1,m-1
temp = temp + l(i,k)*d(k,k)*r(k,m)
enddo
l(i,m) = (a(i,m)- temp)/d(m,m)
temp = 0d0
enddo

enddo

end subroutine

