PROGRAM Escape_tent
!---------------------------------------------------------------------
!
!  This program computes the escape_tent Map 
!
!  Physica A 367 (2006) 158-172
!  
!  Programmer Danilo Szezech - 22/08/2016
!---------------------------------------------------------------------
IMPLICIT NONE
! Declare local variables
INTEGER :: i,j,tfinal,P,k,rede
Real,dimension(0:100) :: x,xn
real, external :: f
Real :: sigma,pi,soma,r
pi=4.*atan(1.0)
OPEN (UNIT=1,FILE='mapa.dat',STATUS='UNKNOWN')
!

! Parametros e CIs
rede=100

! Lendo condições iniciais
OPEN (UNIT=3, FILE='condicoes_iniciais.dat', STATUS='OLD')
read(3,*) (x(k), k=1, rede) 
close(3)


! Lendo parâmetros de rede
OPEN (UNIT=4,FILE='parametros_rede.dat',STATUS='OLD')
read(4,*) tfinal
read(4,*) sigma
read(4,*) r 
close(4)

P=int(rede*r)

!
! Mapa logistico
!
do i=1,tfinal
!atualiza os vetores com o mapa
    do j=1,rede
        x(j)=f(x(j))
    enddo
!
!
! Finite Range Coupling
!
    do j=1,rede
        soma=0.0
        do k=-P,P
            soma=soma+x(modulo(j+K,rede+1))-x(j)
        enddo
    xn(j)=x(j)+(sigma/(2.*P))*soma
    enddo
    do j = 1,rede
        x(j)=xn(j)
enddo

! Write out values of the map
!if(i > tfinal - 400 .and. mod(i,2) == 0) write(1,100) (x(k),k=0,rede)
write(1,100) (x(k),k=1,rede)
enddo
!
!-------------------------------------------------------------------
!  Call system("gnuplot < commands2.txt")
! Call system('rm -f saida.dat')
OPEN(10,ACCESS='SEQUENTIAL',FILE='gp.txt')
write(10,*)'set term png size 800,640'
write(10,*)'set output "Sigma ',sigma,'.png"'
write(10,*)'set grid'
write(10,*)'set style data points'
write(10,*)'unset key'
write(10,*)'set title "Sigma ',sigma,'"'
write(10,*)'set ylabel "<---Y--->"'
write(10,*)'set xlabel "<---X--->"'
write(10,*)'set view 0,0'
write(10,*)'splot "mapa.dat" matrix w pm3d'
write(10,*)'quit'
close(10)
Call system("gnuplot < gp.txt")

!Call system('rm -f gp.txt ')
!Call system('rm -f saida.dat ')
66 format(a21,I2,a5)
67 format(a,a,a)
100 format(101f11.6)
close(1)
!
END PROGRAM Escape_tent

real function f(x)
    implicit none
    real, intent(in) :: x
    real, parameter :: r=3.80

if(x > 0. .and. x < 1.) f=r*x*(1.-x)
end function f