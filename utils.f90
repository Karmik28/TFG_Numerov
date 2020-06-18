module utils
implicit none
integer :: Z, A, nn, l, V_selector
double precision :: E, mu, m_nucl,m_elect, m_xi, V0
double precision, parameter :: a0 = 5.292d-11, r0 = 1.2d-15

contains

double precision function fdx(x)
        double precision :: x
                if (x < RN()*10) then
                        fdx = 0.00001d0/mu
                        ! fdx = 0.00002d0
                else
                        fdx = 0.001d0/mu
                endif

end
double precision function numerov(dx, y0,y1, g0,g1,g2, s0,s1,s2)
        double precision :: dx, y0,y1, g0,g1,g2, s0,s1,s2, den, num
        
        num = 2*y1*(1 - 5*g1*dx**2/12) - y0*(1 + g0*dx**2/12) + (s2 + 10*s1 + s0)*dx**2/12
        den = 1 + g2*dx**2/12
        numerov = num/den
end


double precision function numerov_varstep(dx1,dx2, y0,y1, g0,g1,g2)
        double precision :: dx1, dx2, y0,y1,g0,g1,g2
        double precision :: num, den, c
        c = dx2/dx1
        num = (1+c)*y1 - c*y0 - dx1**2/12 * ((1 + 4*c +4*c**2 +c**3)*g1*y1 + (c + c**2 -c**3)*g0*y0)
        den = 1 + dx1**2/12 * (c**2 +c -1)*g2

        ! num = 2*y1 -y0 -dx1**2/12 *(10*g1 + g0)
        ! den = 1+dx1**2/12 * g2
        numerov_varstep = num/den

end 

! g function
double precision function g(r)
        double precision :: r 
        g = 2*mu*(E-V(r)) - (dble(l)*(dble(l)+1))/r**2
end

!Nucleus Radius
double precision function RN()
        RN = r0 *dble(A)**(1.d0/3) /a0
end



! Wood-Saxon potential
double precision function V(r)
        double precision :: r, a, Vsph
        ! a_mu takes into account value of rest mass cascade particle
        a = 0.6d-15
        ! a = 0.6d-15/2586

        if (V_selector == 1) then
                  V = -1.d0*Z/r

        else if (V_selector == 2) then

                if (r < RN()) then 
                        V = -1.d0*Z/2 * (3*RN()**2-r**2)/RN()**3
                
                else
                        V = -1.d0*Z/r
                endif
        else if (V_selector == 3) then

                if (r < RN()) then 
                        Vsph = -1.d0*Z/2 * (3*RN()**2-r**2)/RN()**3
                
                else
                        Vsph = -1.d0*Z/r
                endif
                V = Vsph + V0/(1+exp((r-RN())/(a/a0)))
        endif
end



!Integrates square norm of function with trapezoidal rule
double precision function quadnorm(x, y)
        integer :: i
        double precision :: x(:), y(:), counter
        counter = 0.d0
        do i=1, size(y) - 1 
                counter = counter + (y(i)**2+y(i+1)**2) * abs(x(i+1)-x(i))
        enddo
        quadnorm =  sqrt(counter/2)
end

subroutine solver_numerov(rf, dr, r_l, P_l, direction)
        integer :: i, N
        double precision, intent(out),allocatable :: r_l(:), P_l(:)
        double precision, intent(in) :: rf
        double precision :: dr, r
        character(1), intent(in):: direction

        r = dr(epsilon(rf))
        N = 1


! Determine size of the arrays to be allocated
        do while (r < rf)
                r = r +dr(r)
                N = N + 1
        end do
       
        allocate(r_l(N))
        allocate(P_l(N))

        r_l(1) = dr(epsilon(rf))
        do i=1, size(r_l)-1
        r_l(i+1) = r_l(i) + dr(r_l(i))
        end do

! Set initial conditions depending on side we integrate from         
        if (direction == "L") then
                ! P_l(1) = r_l(1)**(l+1)
                ! P_l(2) = r_l(2)**(l+1)
                P_l(1) =r_l(1)**(l+1)
                P_l(2) =r_l(2)**(l+1)

        else if (direction == "R") then
                r_l = r_l(N:1:-1)
                P_l(1) = exp(-sqrt(-2*mu*E)*r_l(1))
                P_l(2) = exp(-sqrt(-2*mu*E)*r_l(2))
        else
                write(*,*) direction,"  Error, integration direction unspecified"
                stop
        endif

        ! Apply numerov's algorithm to the whole array
        do i=1, size(P_l)-2
                P_l(i+2) = numerov_varstep(dr(r_l(i)),dr(r_l(i+1)), P_l(i),P_l(i+1), g(r_l(i)), g(r_l(i+1)), &
                             g(r_l(i+2)) )
        end do


        P_l = P_l / quadnorm(r_l, P_l)
        if (direction == "R") then 
                P_l = P_l(N:1:-1) * (-1)**(nn-l+1)
                r_l = r_l(N:1:-1)
        endif

end


integer function derivs(xl, yl, xr,yr,ic)
        double precision, dimension(:), intent(in) :: xl,yl,xr,yr
        double precision :: difl,difr,diff, Vlist(size(xl))
        integer :: ic, i
        do i=1, size(xl) 
         Vlist(i) = V(xl(i))
        end do
        ic = minloc(abs(Vlist-E),DIM=1)
        difl = abs((yl(ic)-yl(ic-1)) / ( (xl(ic)-xl(ic-1))*yl(ic)))
        difr = abs((yr(ic)-yr(ic-1)) / ( (xr(ic)-xr(ic-1))*yr(ic)))
        diff = difr-difl
        if (diff > 0.d0) then
                 derivs = 1
         else if (diff < 0.d0) then
                 derivs = 0
         endif
         write(*,*) "DIFF", diff, xl(ic), xr(ic)
end
        

subroutine bisection(xf,A,B, eps, C,r,pl,pr,P)
        double precision, intent(in) :: A,B,eps, xf
        double precision ::AA, BB
        double precision, allocatable,intent(out) :: r(:), P(:)
        double precision, allocatable :: rl(:), pl(:), rr(:), pr(:)
        double precision, intent(out) :: C
        integer :: N, i,condV,ic, sz
        N = nint(log(abs(A-B)/eps)/log(2.d0))+1
        
        AA = A ; BB = B

        do i=1, N 
                E = (AA+BB)/2
                call solver_numerov(xf,fdx,rl,pl,"L")
                call solver_numerov(xf,fdx,rr,pr,"R")
                condV = derivs(rl,pl,rr,-pr,ic)
                if (condV == 0) then
                        BB = E

                else if ( condV == 1) then
                        AA = E
                endif
                write(*,*) "Bisection iteration",i,"E:",E
        enddo 
        sz = size(rl)
        allocate(r(sz))
        allocate(P(sz))
        r = rl
        P(1:ic) = Pl(1:ic)
        P(ic+1:sz) = Pl(ic)/Pr(ic)*Pr(ic+1:sz)
        P = P/quadnorm(r,P)
        C = E
end


end module
