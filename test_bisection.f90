include "utils.f90"
program test
use utils
implicit none
integer :: i
double precision,allocatable :: outr(:),outPl(:), outPr(:), outP(:), outex(:), outex2(:)
double precision ::  E_c,C, Ab, Bb, moment1, moment2,var

m_elect=0.51099895000d0 

!Values to tweak
Z=8
A=16
nn=3
l=2
! 1=Coulomb potential, 2=Finite-size potential, 3=Wood-Saxon + Finite-size
V_selector = 3
!Potential Well depth (in hartree units)
V0 = -10d0*10**6/27.211
V0 = -20d0*10**6/27.211
!-----

m_nucl = A*931.5d0/m_elect
! m_nucl = (Z*938.27208816d0 + (A-Z)*939.56542052d0)/m_elect 
m_xi = 1321.71d0/m_elect
mu = m_xi*m_nucl/(m_xi+m_nucl)

!Calculate variance of r to get an idea of the endpoint in which wavefunction converges to 0
moment1 = 1.d0/(2*Z*mu)*(3*nn**2-l*(l+1))
moment2 = 1.d0/(2*Z**2*mu**2)*nn**2*(5*nn**2+1-3*l*(l+1))
var = sqrt(moment2-moment1**2)
write(*,*) var


open(1, file="test_bisection.txt")
        E_c = -1.d0*Z**2*mu/(2*nn**2)
        write(*,*) "MU", mu
        write(*,*) "COLOUMB ENERGY", E_c
        Ab = E_c-5000.d0
        Bb = E_c+5000.d0
call bisection(20*var,Ab,Bb,0.0000001d0,C,outr, outPl,outPr,outP)
write(*,*) "Energy in keV", C*27.211386245988d0/1000
write(*,*) "Shift with respect to Coulombian in  keV", (C-E_c)*27.211386245988d0/1000

! Generate data predicted by analytical methods for Coloumb potential
call func_data(0.0d0, 20*var, fdx, Pnl, outex, outex2)
write(*,*) quadnorm(outr,outP), quadnorm(outex,outex2), RN()

! Export data
write(1,'(F20.12, 2X,F20.12,2X,F20.12,2X,F20.12, 2X, F20.12,2X, F20.12)')&
                (outr(i), outPl(i),outPr(i),outP(i),outex(i), outex2(i), i=1, min(size(outPl),size(outPr),size(outex)))
! Execute plotting
! CALL execute_command_line('python "plotter.py"')
CALL execute_command_line('gnuplot "plot.gpi"')

contains


! Generates data from a defined function
subroutine func_data(xi,xf,dx,func, out_xdata,out_ydata)
        double precision,intent(in) :: xi, xf
        double precision,intent(out),allocatable :: out_xdata(:), out_ydata(:)
        integer :: i, N
        double precision :: x, func, dx

        x = xi
        N = 1

        do while (x < xf)
        x = x + dx(x)
        N = N + 1 
        end do
        
        x = xi
        allocate(out_xdata(N))
        allocate(out_ydata(N))
        do i = 1, N
                if (i == N-1) then
                        x = xf 
                else
                        x = x + dx(x)
                endif
                out_xdata(i) = x
                out_ydata(i) = func(x)
        end do
end

!Analytical function for radial wave function of hydrogen-like atoms
double precision function Pnl(r)
        double precision :: N,r, ro
        ro = 2.d0*Z*mu/nn
        N = ro**(3.d0/2) * ((gamma(dble(nn-l)))/(2*nn*((gamma(dble(nn+l+1)))**3)))**(1.d0/2)
        Pnl = N * r* (ro*r)**dble(l) * exp(-ro*r/2) *Lf(ro*r)
        if (Lf(ro*r) > HUGE(r)) then
                write(*,*) "Too big"
        endif
 end

 ! generalized laguerre polynomials
 double precision function Lf(x)
         double precision :: x, Lf0, Lf1, num, den, temp
         integer :: i, k, alfa
         k = nn-l-1
         alfa = 2*l+1
         temp=0.d0
         

         Lf0 = 1.d0
         Lf1 = 1.d0 + alfa - x

         do i=0, k
                 num = (-1.d0)**i * (gamma(dble(nn+l+1)))**2 * x**i
                 den = gamma(dble(i+1))*gamma(dble(nn-l-i))*gamma(dble(2*l+2+i))
                 temp = temp +num/den
         end do
         Lf = temp
         

end
end program

