!Last change:  YKO  May 14 2018    3:26 pm_v1
module param
DOUBLE PRECISION :: sigma, beta, r, b_inf, b_sup, b_inf_LC, b_sup_LC, y_inf, y_sup, eps_inf, eps_sup, std_eps, pi_inf, pi_sup, &
                    rho, mean_y, std_y, width, prob_excl_end, coupon, delta, gamma, d0, d1, cdf_inf, alpha0, alpha1, &
                    cdf_sup, b_recov, kappa, absortion, slope_ss, prob_ss, premium_grid_value, recovery, theta_cost, zero
                    
parameter (beta = 0.92, sigma = 2d+0, std_eps = 0.0280d+0, rho = 0.69d+0, mean_y = -0.5d+0*std_eps**2, width = 1.5d+0, cdf_inf = 3.167d-5, cdf_sup = 1d+0 - 3.167d-5,  &
          gamma = 0.0d+0, r = 0.04, kappa = 0.69d+0, absortion = 0.0d+0, prob_ss = 0.15, recovery = 0.5d+0, alpha0=0d+0, alpha1=0d+0, theta_cost = 2000.20d+0, zero = 0d+0,&
          d0 = -1.85d+0, d1 = 2.242d+0, b_recov = 0d+0, dur_excl = 1.5d+0, delta = 0.155d+0, premium_grid_value = 0d+0, slope_ss = 0, prob_excl_end = 1d+0/ dur_excl)

!sigma = coefficient of relative risk aversion in private consumption
!r = international interest rate
!rho = coefficient of autocorrelation in income
!std_eps = standard deviation of innovation to income shocks
!mean_y = unconditional mean of the growth trend shockq
!std_y = standard deviation of growth trend shock
!b_recov = number of bond extended
!slope_ss = SENSITIVITY OF SS PROBABILITY WRT CURRENT INCOME
!prob_ss = PARAMETER THAT REGULATES THE LONG RUN PROBABILITY OF SWITCHES FROM LOW TO HIGH RISK PREMIUM.
INTEGER :: b_num_LC, b_num_long, y_num, quad_num, premium_num, nout, i_y_global, i_excl_global, theta_num, &
           i_b_LC_global, i_b_long_global, i_default_global, indicator_putazo, i_b_LC_next, i_b_long_next,&
           b_num_LC_finer, b_num_long_finer, y_num_finer, i_premium_global, index_premium_initial,&
           num_LC_next, y_num_coarse, b_num_long_coarse, b_num_LC_coarse, net_num

parameter (pi_num = 3, b_num_LC = 3, b_num_long = 25, y_num =25, quad_num = 50, premium_num = 1, &
           b_num_LC_finer = 3, b_num_long_finer = 50, y_num_finer =50, num_LC_next = 3, theta_num =1,&
           y_num_coarse =15, b_num_long_coarse = 30, b_num_LC_coarse = 30, net_num = 100, num_debt = 40, num_pi = pi_num)

DOUBLE PRECISION :: b_grid_long(1:b_num_long), b_grid_LC(1:b_num_LC), premium_grid(2), y_grid(1:y_num), default_grid(1:2), indicator_tirar,&
                    vector_v_global(y_num), y_initial, g_initial, b_LC_initial, b_long_initial, b_global,counter, b_long_global, &
                    b_LC_global, q_LC_global, quad_w(1:quad_num), quad_x(1:quad_num), pi_num_grid(num_pi),&
                    b_LC_global_optimize, b_grid_long_finer(1:b_num_long_finer), b_grid_LC_finer(1:b_num_LC_finer), &
                    y_grid_finer(1:y_num_finer), b_next_LC_grid(num_LC_next), &
                    b_grid_long_coarse(1:b_num_long), b_grid_LC_coarse(1:b_num_LC), y_grid_coarse(1:y_num), indicator_global_search, cons_global, &
                    b_next_LC_global, b_next_long_global, c_global, g_global, output_global,&
                    net_grid(net_num), net_global, q_global, b_LC_next_global, util_global, ev_global,pi_global_optimize, pi_global

!y_initial = growth rate state. USED TO SOLVE OPTIMAL SAVINGS RULE AT EVERY ITERATION
!b_initial = debt state. USED TO SOLVE OPTIMAL SAVINGS RULE AT EVERY ITERATION
DOUBLE PRECISION, DIMENSION(b_num_LC, b_num_long, y_num, premium_num) :: v_matrix
INTEGER, DIMENSION(b_num_LC, b_num_long, y_num, 2) :: default_decision
DOUBLE PRECISION, DIMENSION(b_num_LC, b_num_long, y_num, premium_num,2) :: b_next_matrix
DOUBLE PRECISION, DIMENSION(b_num_LC, b_num_long, y_num, premium_num) :: v0_matrix, b0_next_LC_matrix, b1_next_LC_matrix, b1_next_long_matrix,&
 q_paid_matrix, b0_next_long_matrix, v1_matrix, q_nodef_matrix, q_nodef_LC_matrix, q_def_matrix, g0_matrix, g1_matrix, q_LC_def_matrix, pi_matrix
DOUBLE PRECISION, DIMENSION(premium_num, premium_num) :: trans_matrix
DOUBLE PRECISION, DIMENSION (b_num_LC_finer, b_num_long_finer, y_num_finer, 2) :: EV_matrix, q_menu_matrix, q_LC_menu_matrix, EV_excl_matrix, q_menu_def_matrix, q_LC_menu_def_matrix

!needed to read files created with linear interpolation
INTEGER :: b_num_LC_linear, b_num_long_linear, y_num_linear
PARAMETER(b_num_LC_linear = 20, b_num_long_linear = 20, y_num_linear = 25)
double precision :: b_grid_LC_linear(b_num_LC_linear), b_grid_long_linear(b_num_long_linear), y_grid_linear(y_num_linear)

 
!needed to use BS2IN - spline interpolation in 2 dimensions
INTEGER :: KORDER
PARAMETER( KORDER=3 )

INTEGER :: dim_auxiliar_long, dim_auxiliar_LC, dim_auxiliar_long_finer, dim_auxiliar_LC_finer
PARAMETER(dim_auxiliar_long = KORDER + b_num_long, dim_auxiliar_LC = KORDER + b_num_LC,&
          dim_auxiliar_long_finer = KORDER + b_num_long_finer, dim_auxiliar_LC_finer = KORDER + b_num_LC_finer )

DOUBLE PRECISION, DIMENSION(b_num_LC_finer, y_num_finer, premium_num) :: break_matrix_b_limit
DOUBLE PRECISION, DIMENSION(4, b_num_LC_finer, y_num_finer, premium_num) :: coeff_matrix_b_limit


DOUBLE PRECISION, DIMENSION (dim_auxiliar_long) :: b_long_knots
DOUBLE PRECISION, DIMENSION (dim_auxiliar_LC) :: b_LC_knots

DOUBLE PRECISION, DIMENSION (dim_auxiliar_long_finer) :: b_long_knots_finer
DOUBLE PRECISION, DIMENSION (dim_auxiliar_LC_finer) :: b_LC_knots_finer

DOUBLE PRECISION, DIMENSION(b_num_LC, b_num_long, y_num, premium_num) :: coeff_matrix_v0, coeff_matrix_q, coeff_matrix_q_LC, coeff_matrix_q_def, coeff_matrix_v1, coeff_matrix_q_LC_def, coeff_matrix_pi
DOUBLE PRECISION, DIMENSION(b_num_LC_finer, b_num_long_finer, y_num_finer, premium_num) :: coeff_matrix_ev, coeff_matrix_q_menu, coeff_matrix_q_LC_menu, coeff_matrix_ev_excl, coeff_matrix_q_LC_menu_def, coeff_matrix_q_menu_def

end module




!SPECIFY GRID VALUES
subroutine compute_grid
USE param
DOUBLE PRECISION :: dos, DNORDF, delta_y, prob_mass, b_inf0, b_sup0, b_inf1, b_sup1, y_left, &
prob_vector(y_num), y_right, min_output, max_debt_min_revenue
INTEGER :: b_num_half, i,j
EXTERNAL :: DNORDF

coupon = 1d+0 - exp(-r)*(1-delta) 
std_y = std_eps/ SQRT(1 - rho**2)

y_inf = mean_y - 5*std_y
y_sup = mean_y + 5*std_y
delta_y = 5*std_y / (y_num - 1d+0)

b_inf = 0.00001 !MIN(0.55d+0*EXP(y_sup), .99*EXP(y_inf))
b_sup =  min(1.8d+0, exp(y_inf)/coupon-0.05d+0)

b_inf_LC = 0.00001d+0 !MIN(0.55d+0*EXP(y_sup), .99*EXP(y_inf))
b_sup_LC =  0.00002d+0

pi_inf = 0.000001d+0
pi_sup = 0.5d+0

open (10, FILE='graphs\b_grid_LC.txt',STATUS='replace')
open (11, FILE='graphs\y_grid.txt',STATUS='replace')
open (12, FILE='graphs\b_grid_long.txt',STATUS='replace')
open (13, FILE='graphs\bounds.txt',STATUS='replace')


WRITE(13, '(F15.11, X, F15.11)') b_inf, b_sup
WRITE(13, '(F15.11, X, F15.11)') b_inf_LC, b_sup_LC
WRITE(13, '(F15.11, X, F15.11)') y_inf, y_sup


do i=1,b_num_LC
   b_grid_LC(i) = b_inf_LC + (b_sup_LC - b_inf_LC)*(i-1d+0)/(b_num_LC - 1d+0)
end do

do i=1,b_num_long
   b_grid_long(i) = b_inf + (b_sup - b_inf)*(i-1d+0)/(b_num_long - 1d+0)
end do


do i=1,b_num_LC
 WRITE(10, '(F12.8, X, F12.8)') b_grid_LC(i)
end do

do i=1,b_num_long
 WRITE(12, '(F12.8, X, F12.8)') b_grid_long(i)
end do

do i=1,y_num
   y_grid(i) = y_inf + (y_sup - y_inf) * (i-1) / (y_num - 1)
   WRITE(11, '(F12.8)') y_grid(i)
end do



close(10)
CLOSE(11)
CLOSE(12)
CLOSE(13)


open (16, FILE='graphs\trans_matrix.txt',STATUS='replace')
do i=1,premium_num
   do j=1,premium_num
       WRITE(16, '(F12.8)') trans_matrix(i,j)
   end do
end do
CLOSE(16)


dos = 2d+00


default_grid(1) = 0.0d+0
default_grid(2) = 1.0d+0

eps_inf = -4*std_eps
eps_sup = 4*std_eps

open (20, FILE='graphs\b_LC_finer.txt',STATUS='replace')
open (22, FILE='graphs\b_long_finer.txt',STATUS='replace')

do i=1,b_num_LC_finer
   b_grid_LC_finer(i) = b_inf_LC + (b_sup_LC - b_inf_LC)*(i-1d+0)/(b_num_LC_finer - 1d+0)
   WRITE(20, '(F12.8)') b_grid_LC_finer(i)
end do

do i=1,b_num_long_finer
   b_grid_long_finer(i) = b_inf + (b_sup - b_inf)*(i-1d+0)/(b_num_long_finer - 1d+0)
   WRITE(22, '(F12.8)') b_grid_long_finer(i)
end do

close(20)
close(22)
do i=1,y_num_finer
   y_grid_finer(i) = y_inf + (y_sup - y_inf) * (i-1) / (y_num_finer - 1)
end do

pi_num_grid(1) = pi_inf
pi_num_grid(2) = pi_inf + 0.001d+0
if (num_pi > 5) then
    do i=2,num_pi-3
       pi_num_grid(i+1) = pi_inf + (pi_sup - pi_inf) * (i-1)/ (num_pi - 3)
    end do
    pi_num_grid(num_pi-1) = pi_sup-0.001d+0
    pi_num_grid(num_pi) = pi_sup
end if

b_next_LC_grid(1) = b_inf_LC
b_next_LC_grid(2) = b_sup_LC !+ 0.001d+0
do i=2,num_LC_next-3
   b_next_LC_grid(i+1) = b_inf_LC + (b_sup_LC - b_inf_LC) * (i-1)/ (num_LC_next - 3)
end do
b_next_LC_grid(num_LC_next-1) = b_sup_LC-0.001d+0
b_next_LC_grid(num_LC_next) = b_sup_LC

open (16, FILE='graphs\pdf_y.txt',STATUS='replace')
!Aproximate unconditional p.d.f. of y
y_left  = (y_grid(1) + delta_y)/ std_y
y_right = (y_grid(y_num) - delta_y)/ std_y

prob_vector(1)     = DNORDF(y_left)
prob_vector(y_num) = 1d+0 - DNORDF(y_right)
write(16, '(F12.8)') prob_vector(1)
   do j=2,y_num-1
     y_right = (y_grid(j) + delta_y )/ std_y
     y_left  = (y_grid(j) - delta_y )/ std_y
     prob_vector(j) = (DNORDF(y_right) - DNORDF(y_left) ) !/prob_may_num
!     WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)')  y_grid(j), y_left, y_right, prob_vector(j)
     write(16, '(F12.8)') prob_vector(j)
   end do
write(16, '(F12.8)') prob_vector(y_num)
close (16)

end subroutine



!COMPUTE QUADRATURE POINTS AND WEIGHTS USING LEGENDRE AND GAUSSIAN QUADRATURE RULES
!THEY ARE USED TO COMPUTE NUMERICAL INTEGRALS
subroutine quadrature
USE param
INTEGER :: N_v, N_q, N, IWEIGH, IWEIGH4, NFIX
parameter(N = quad_num)
DOUBLE PRECISION:: ALFA, BETA1, QW(1:N), QX(1:N), QXFIX(2)
PARAMETER(ALFA=0, BETA1=0, IWEIGH=1, IWEIGH4=4)
external DGQRUL

!quad_w = weights using Gauss legendre quadrature rule
!quad_x = points using Gauss legendre quadrature rule

NFIX = 0
CALL DGQRUL(N,IWEIGH, ALFA, BETA1, NFIX, QXFIX, QX, QW)

quad_w = QW
quad_x = QX

end subroutine



!COMPUTES THE DIFFERENCE BETWEEN THE VALUE FUNCTION UNDER NO DEFAULT AND DEFAULT
!NEED b_ylobal TO BE DEFINED BEFORE.
DOUBLE PRECISION function dif_fun(y)
USE param
DOUBLE PRECISION, INTENT(IN) :: y
DOUBLE PRECISION :: v0_fun, v1_fun
EXTERNAL v0_fun, v1_fun


dif_fun = v0_fun(b_LC_global, b_long_global,  y, i_premium_global) - v1_fun(b_LC_global, b_long_global, y, i_premium_global)
!WRITE(nout, '(F12.8, X, F12.8, X, F12.8)') y, v0_fun(b_global, y), v1_fun(b_global, y)
end

!FUNCTION USED TO COMPUTE THE MINIMUM INCOME FOR WHICH THE CURRENT PARTY IN POWER DOES NOT DEFAULT
!b = outstanding debt
DOUBLE PRECISION function y_fun(b_LC, b_long, index_premium)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_LC, b_long
INTEGER, INTENT(IN) :: index_premium
INTEGER :: MAXFN, num_tirar
PARAMETER(num_tirar = 100)
DOUBLE PRECISION :: dif_fun, ERRABS, ERRREL, left, right, y_max, y_min, dif_right, dif_left,&
                    y_tirar(num_tirar), tirar
EXTERNAL dif_fun, DZBREN

MAXFN=1000
ERRREL = 1D-10
ERRABS = 1D-10


y_max = rho*y_initial + (1-rho)*mean_y + width*eps_sup
y_min = rho*y_initial + (1-rho)*mean_y + width*eps_inf


b_LC_global = b_LC     !b_global IS USED IN dif_fun TO FIX THE VALUE OF b (OUTSTANDING DEBT)
b_long_global = b_long     !b_global IS USED IN dif_fun TO FIX THE VALUE OF b (OUTSTANDING DEBT)
i_premium_global = index_premium
!open (120, FILE='tirar1.txt',STATUS='replace')
!do i=1,num_tirar
!    tirar = y_min + (y_max - y_min) * (i-1d+0)/(num_tirar-1)
!    left = dif_fun(tirar)
!end do
!CLOSE(120)


dif_right = dif_fun(y_max)
dif_left  = dif_fun(y_min)

!1) DETERMINE WHETHER THERE IS AN INTERIOR ROOT OR NOT
!   a) IF THERE IS AN INTERIOR ROOT, SPECIFY WHETHER THE difference function IS INCREASING IN g OR NOT.
!   b) IF THERE IS NO INTERIOR ROOT, DETERMINE WHETHER THE GOV'T ALWAYS OR NEVER DEFAULTS

!WRITE(nout, *) y_min, y_max, dif_right, dif_left


if (dif_right*dif_left<0) then  !THERE IS AN INTERIOR ROOT
   left =  y_min
   right = y_max
   CALL DZBREN (dif_fun, ERRABS, ERRREL, left, right, MAXFN)
   y_fun = right
else if (dif_left>0) then
      y_fun = y_min
else !dif_right<0
      y_fun = y_max
end if

end


DOUBLE PRECISION function objective_function(X)
USE param
DOUBLE PRECISION, INTENT (IN) :: X(3)
INTEGER :: other_type, i,j,h, t, d
DOUBLE PRECISION :: acum, exp_v_next, value_next, y_next, y_threshold, y_next1, acum2,acum1, q, DNORDF, v0_fun, v1_fun, &
                    borrowing, cdf, cdf1, cdf_threshold, DNORIN, output, b_LC, b_long, EV_fun, q_vector(2), output_cost_ss,&
                    q_menu_fun, penalty, consum, q_LC, util, ind_payment_suspension, q_LC_menu_fun, borrowing_LC

EXTERNAL v0_fun, v1_fun, DNORIN, EV_fun, q_menu_fun, q_LC_menu_fun

b_LC = MAX(MIN(X(1), b_sup_LC), b_inf_LC)
b_long  = MAX(MIN(X(2), b_sup), b_inf)
pi  = MAX(MIN(X(3), pi_sup), pi_inf)
penalty = 10d+0*(max(pi_inf-X(3),0d+0)**2d+0 + max(X(3)-pi_sup,0d+0)**2d+0 + max(b_inf-X(2),0d+0)**2d+0 + max(X(2)-b_sup,0d+0)**2d+0 + max(b_inf_LC-X(1),0d+0)**2d+0 + max(X(1)-b_sup_LC,0d+0)**2d+0)

!CALCULATE EXPECTED CONTINUATION VALUE
exp_v_next = EV_fun(b_LC, b_long, y_initial, index_premium_initial)

q = q_menu_fun(b_LC, b_long, y_initial, index_premium_initial)
q_LC = q_LC_menu_fun(b_LC, b_long, y_initial, index_premium_initial)
borrowing =  b_long - b_long_initial*(1d+0 - delta)
borrowing_LC =  b_LC - b_LC_initial*(1-delta)/(1d+0 + pi)
output = EXP(y_initial) 

consum = ((1-pi)**theta_cost)*output - coupon*b_long_initial + borrowing*q + borrowing_LC*q_LC - b_LC_initial*coupon/(1d+0 + pi) - absortion
util = (max(consum, 1d-6))**(1-sigma)/(1-sigma)  !COMPUTE CURRENT UTILITY FLOW
if (borrowing_LC < 0d+0 ) then !this imposes no-buyback constraint
    objective_function = 10d+16
elseif(borrowing_LC > 0d+0 .AND. q < q_global) then !this imposes the borrowing constraint and q constraint
    objective_function = 10d+16
else
    objective_function =  - util - beta*exp_v_next + penalty
end if

c_global = consum
util_global = util
ev_global = exp_v_next

end function

DOUBLE PRECISION function objective_function_excl(X)
USE param
DOUBLE PRECISION, INTENT(IN) :: X(2)
INTEGER :: other_type, i,j,h, t, d
DOUBLE PRECISION :: q_fun, acum, exp_v_excl, value_next, y_next, y_fun, y_threshold, y_next1, acum2,acum1, q, v0_fun, v1_fun, &
                    value_next_tirar, acum_int, borrowing, cdf, cdf1, cdf_threshold, output, b_LC, b_long, q_LC, exp_v_next,&
                    EV_excl_fun, util, interpolate, consum

EXTERNAL q_fun, y_fun, v0_fun, v1_fun, EV_excl_fun, interpolate

b_LC = X(1)
b_long = X(2)

output = min(EXP(y_initial)*(1d+0 - d0) - d1*EXP(y_initial)**2d+0, exp(y_initial)) !EXP(y_initial)
consum = output  - absortion
util = (max(consum, 1d-6))**(1-sigma)/(1-sigma)  !COMPUTE CURRENT UTILITY FLOW

exp_v_next = EV_excl_fun(b_LC*exp(r), b_long*exp(r), y_initial, index_premium_initial) !COMPUTE CONTINUATION VALUE UNDER DEFAULT
objective_function_excl = util + beta * exp_v_next 
!write(6, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)') b_LC, b_long, y_initial, output, consum, exp_v_next, objective_function_excl

c_global = consum
output_global = output
end


DOUBLE PRECISION function objective_function_debt(debt)
USE param
DOUBLE PRECISION, INTENT(IN) :: debt
DOUBLE PRECISION :: objective_function, b_LC, b_long, pi
EXTERNAL objective_function

b_LC = b_LC_global_optimize
b_long  = MAX(MIN(debt, b_sup), b_inf)
pi = pi_global_optimize

objective_function_debt = objective_function((/b_LC, b_long, pi/))

end

    
subroutine optimize(b_next, v_value)
USE param
integer :: t, d, index_opt, index_vector(1), num, num_tirar, N, IPARAM(7), IBTYPE,IERSVR, IPACT, ISACT,&
           index_LC, index_long, i_premium, i, ITER, NP, num_long, index_opt_pi, index_opt_LC
parameter (num_long=100, num_tirar = 200, N=3, NP = 3)
DOUBLE PRECISION :: b_next(3), v_value, new_value, q_fun, old_value, b, q,vector(b_num_long),&
                    b_next_LC, b_next_long, b_next_guess(3), b_tirar(2), b_next_inf, b_next_sup, &
                    difference, XSCALE(N), FSCALE, RPARAM(7), FVALUE, XUB(N), XLB(N), fun_right, fun_left,&
                    derivative, eps,  y_threshold, y_fun, y_current_type_1, y_current_type_2,DNORDF, q_menu_fun,&
                    objective_function, b_LC_next_inf, b_LC_next_sup, FTOL, matrix_direction(N,N), &
                    b_next_graph(2), pi

EXTERNAL q_fun, DUVMIF, objective_function, ERSET, y_fun, DNORDF, q_menu_fun

!if (indicator_tirar>0) then
!   open (117, FILE='tirar1.txt',STATUS='replace')
!end if



!IERSVR = 0!4
!IPACT = 0!1
!ISACT = 0!1

!call ERSET(IERSVR, IPACT, ISACT)

b_next_inf = b_inf
b_next_sup = b_sup
index_opt_pi=1
index_opt_LC=1
index_opt = 1
old_value = 10d+15
if (indicator_global_search > 0) then  !SEARCH FOR BEST INITIAL GUESS USING ALL GRID RANGE
    do j = 1, num_pi
        pi = pi_num_grid(j)
        pi_global_optimize = pi
        do i=1,num_LC_next
            b_next_LC = b_next_LC_grid(i)
            b_LC_global_optimize = b_next_LC
            !FIND THE OPTIMAL DEBT ACCUMULATION AND IMPLIED CONTINUATION VALUE FOR A CHOICE OF b_next_LC and pi
            call optimize_debt(b_next_long, new_value)
            if (new_value <= old_value) then
                index_opt_LC = i 
                index_opt_pi = j
                b_next_guess(2) = b_next_long
                b_next_guess(3) = pi
                b_next_guess(1) = b_next_LC
                old_value = new_value
                index_opt_LC = i
            end if
            !if (indicator_tirar > 1) then
            !write(nout, '(F10.6, X, I3, X, G15.6, X, F10.6, X, F10.6, X, F10.6)') b_next_LC, index_opt_LC, b_next_long, pi, new_value!, &
            !end if
            !b_next_guess(1), b_next_guess(2)
        end do
    end do
else !FOCUS GRID SEARCH IN A NEIGHBORHOOD OF PREVIOUS OPTIMAL VALUE
    if (i_default_global < 2) then
        b_next_guess(1) = b0_next_LC_matrix(i_b_LC_global, i_b_long_global, i_y_global, index_premium_initial)
        b_next_guess(2) = b0_next_long_matrix(i_b_LC_global, i_b_long_global, i_y_global, index_premium_initial)
        b_next_guess(3) = pi_matrix(i_b_LC_global, i_b_long_global, i_y_global, index_premium_initial)
    else
        b_next_guess(1) = b1_next_LC_matrix(i_b_LC_global, i_b_long_global, i_y_global, index_premium_initial)
        b_next_guess(2) = b1_next_long_matrix(i_b_LC_global, i_b_long_global, i_y_global, index_premium_initial)
        b_next_guess(3) = pi_matrix(i_b_LC_global, i_b_long_global, i_y_global, index_premium_initial)
    end if
    index_opt_LC = 2 !SET IT TO A NUMBER > 1 & < num_LC_next SO THAT THE OPTIMIZATION ROUTINE IS INVOKED
    index_opt_pi = 2
end if

XSCALE = 1d+0
FSCALE = 1d+0
IPARAM(1) = 0d+0
IBTYPE = 0
XLB(1) = b_inf_LC
XLB(2) = b_inf
XLB(3) = pi_inf
XUB(1) = b_sup_LC
XUB(2) = b_sup
XUB(3) = pi_sup

if (index_opt_LC ==1) then
    b_next = b_next_guess
    v_value = -objective_function(b_next)
elseif(index_opt_LC ==num_LC_next ) then
    b_next = b_next_guess
    v_value = -objective_function(b_next)
else
    matrix_direction = 0d+0
    matrix_direction(1,1) = 1d+0
    matrix_direction(2,2) = 1d+0
    matrix_direction(3,3) = 1d+0
    FTOL = 1d-8

    !call POWELL(objective_function, b_next_guess, matrix_direction, N, NP, FTOL, ITER, FVALUE)
    call POWELL(b_next_guess, matrix_direction, N, NP, FTOL, ITER, FVALUE)
    b_next(1) = MAX(MIN(b_next_guess(1), b_sup_LC), b_inf_LC)
    b_next(2)  = MAX(MIN(b_next_guess(2), b_sup), b_inf)
    b_next(3)  = MAX(MIN(b_next_guess(3), pi_sup), pi_inf)
    pi_global = b_next(3)
    v_value = -FVALUE

end if

10 end subroutine


subroutine optimize_debt(b_next, v_value)
USE param
integer :: MAXFN, t, d, index_opt, index_vector(1), num, i
parameter (num=15)
DOUBLE PRECISION :: b_next, v_value, objective_function_debt, new_value, q_fun, old_value, b, q, q_menu_fun,&
                    b_next_grid(num_debt), STEP, BOUND, XACC, b_next_guess, b_tirar_sup, b_tirar_inf, b_next_inf, b_next_sup,&
                    ancho, debt_inf, debt_sup
EXTERNAL objective_function_debt, DUVMIF, q_menu_fun

!open(100, FILE ='tirar.txt',STATUS='replace')

b_next_inf = b_inf
b_next_sup = b_sup

if (indicator_global_search > 0) then  !SEARCH FOR BEST INITIAL GUESS USING ALL GRID RANGE
  old_value = 10d+15
  index_opt = 1
   do i=1,num_debt
    !DEFINE THE GRID IN A NEIGHBORHOOD OF THE PREVIOUS OPTIMAL VALUE
!    b_next_grid(i) = debt_inf + (debt_sup - debt_inf) * (i-1)/ (num_debt - 1) !DEFINE GRID BASED ON PREVIOUS OPTIMAL VALUE
     b_next_grid(i) = b_next_inf + (b_next_sup - b_next_inf) * (i-1)/ (num_debt - 1)
     new_value = objective_function_debt(b_next_grid(i))
     if (new_value <= old_value) then
        index_opt = i
        b_next_guess = b_next_grid(i)
        old_value = new_value
     end if
   end do
else !FOCUS GRID SEARCH IN A NEIGHBORHOOD OF PREVIOUS OPTIMAL VALUE


   ancho = 0.25d+0
!IF POLICY FUNCTIONS HAVE CONVERGED, USE PREVIOUS OPTIMAL VALUE FOR DEBT TO NARROW GRID SEARCH
 debt_inf = MAX(b_next_inf, b0_next_long_matrix(i_b_LC_global, i_b_long_global, i_y_global, index_premium_initial) - ancho)
 debt_sup = MIN(b_next_sup, b0_next_long_matrix(i_b_LC_global, i_b_long_global, i_y_global, index_premium_initial) + ancho)

  old_value = 10d+15
  index_opt = 1
  do i=1,num_debt
    !DEFINE THE GRID IN A NEIGHBORHOOD OF THE PREVIOUS OPTIMAL VALUE
    b_next_grid(i) = debt_inf + (debt_sup - debt_inf) * (i-1)/ (num_debt - 1) !DEFINE GRID BASED ON PREVIOUS OPTIMAL VALUE
    new_value = objective_function_debt(b_next_grid(i))
    if (new_value <= old_value) then
       index_opt = i
       b_next_guess = b_next_grid(i)
       old_value = new_value
    end if
!write(nout, '(I3, X, F12.8, X, I3, X, F12.8, X, F12.8)') i, b_next_grid(i), index_opt, new_value, &
!q_menu_fun(b_LC_global_optimize, b_next_grid(i), y_initial, index_premium_initial)

!write(100, '(F12.8, X, F12.8, X, F12.8)') b_next_grid(i), MIN(new_value, 100)
  end do
end if

!CLOSE(100)
!pause
STEP = (b_next_grid(2) - b_next_grid(1))*0.5
BOUND = 10*EXP(y_sup)
XACC = 1d-5
MAXFN = 1000



b_tirar_sup = b_next_sup - 1.0d-4
b_tirar_inf = b_next_inf + 1.0d-4

if (index_opt == 1 .AND. objective_function_debt(b_tirar_inf) >= objective_function_debt(b_next_grid(index_opt))) then
   b_next = b_next_guess
   v_value = old_value

elseif(index_opt == 1 .AND. objective_function_debt(b_tirar_inf) < objective_function_debt(b_next_grid(index_opt))) then
       !THE MAXIMUM MAY BE CLOSE TO MINIMUM GRID POINT BUT NOT AT MINIMUM GRID POINT
       b_next_guess = b_tirar_inf
       call DUVMIF(objective_function_debt, b_next_guess, STEP, BOUND, XACC, MAXFN, b_next)
       v_value = objective_function_debt(b_next)
elseif (index_opt == num_debt .AND. objective_function_debt(b_tirar_sup) >= objective_function_debt(b_next_grid(index_opt))) then
       b_next = b_next_guess
       v_value = old_value
elseif(index_opt == num_debt .AND. objective_function_debt(b_tirar_sup) < objective_function_debt(b_next_grid(index_opt))) then
       !THE MAXIMUM MAY BE CLOSE TO MAXIMUM GRID POINT BUT NOT AT MAXIMUM GRID POINT
       b_next_guess = b_tirar_sup
       call DUVMIF(objective_function_debt, b_next_guess, STEP, BOUND, XACC, MAXFN, b_next)
       v_value = objective_function_debt(b_next)
else
   if (ABS(old_value - objective_function_debt(b_next_grid(index_opt-1)))/(b_next_grid(2)-b_next_grid(1)) < 1d-10 .OR.  &
       ABS(old_value - objective_function_debt(b_next_grid(index_opt+1)))/(b_next_grid(2)-b_next_grid(1)) < 1d-10) then
       b_next = b_next_guess
       v_value = old_value
   else
      call DUVMIF(objective_function_debt, b_next_guess, STEP, BOUND, XACC, MAXFN, b_next)
       v_value = objective_function_debt(b_next)
   end if
end if

10 end subroutine


    
DOUBLE PRECISION function v0_fun(b_LC, b_long, y, index_premium)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_LC, b_long, y
INTEGER, INTENT(IN) :: index_premium
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl, interpolate3
INTEGER ::  NINTV, index_y
EXTERNAL DBS2VL, interpolate3

bs=  min(max(b_inf_LC+1d-6, b_LC),b_sup_LC-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)

v0_fun = interpolate3(bs, bl, y, v0_matrix(:,:,:, index_premium))
end function

   

DOUBLE PRECISION function v1_fun(b_LC, b_long, y, index_premium)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_LC, b_long, y
INTEGER, INTENT(IN) :: index_premium
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl, interpolate3
INTEGER :: index_y, NINTV
EXTERNAL DBS2VL, interpolate3

bs=  min(max(b_inf_LC+1d-6, b_LC),b_sup_LC-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)
v1_fun = interpolate3(bs, bl, y, v1_matrix(:,:,:, index_premium))

end function



DOUBLE PRECISION function interpolate3(bs, bl, y, matrix)
USE param
DOUBLE PRECISION, INTENT(IN) :: bs, bl, y, matrix(b_num_LC, b_num_long, y_num)
DOUBLE PRECISION ::  slope, weight, acum, ratio_bs, ratio_bl, ratio_y, eps
INTEGER :: index_bs, index_bl, index_y, i_bs, i_bl, i_y

eps = 1d-8 !NEEDED TO ADD EPSILON BECAUSE IF y = y_grid(i) IT MAY CHOOSE index_y = i-1 DUE TO APPROX ERROR IN THE ROUNDING
index_y = (MAX(MIN(INT((y_num-1)*(y+eps-y_grid(1))/(y_grid(y_num)-y_grid(1)))+1,y_num-1), 1))
index_bs = (MAX(MIN(INT((b_num_LC-1)*(bs+eps-b_grid_LC(1))/(b_grid_LC(b_num_LC)-b_grid_LC(1)))+1,b_num_LC-1), 1))
index_bl = (MAX(MIN(INT((b_num_long -1)*(bl+eps-b_grid_long(1)) /(b_grid_long(b_num_long)-b_grid_long(1)))+1,b_num_long-1), 1))



ratio_bs = (bs - b_grid_LC(index_bs)) / (b_grid_LC(index_bs+1) - b_grid_LC(index_bs))
ratio_bl = (bl - b_grid_long(index_bl)) / (b_grid_long(index_bl+1) - b_grid_long(index_bl))
ratio_y = (y - y_grid(index_y)) / (y_grid(index_y+1) - y_grid(index_y))
acum=0

do i_bs=0,1
   do i_bl = 0,1
      do i_y=0,1
        weight = ((1-i_bs)*(1 - ratio_bs) + i_bs*ratio_bs) * ((1-i_bl)*(1 - ratio_bl) + i_bl*ratio_bl) * &
                 ((1-i_y)*(1 - ratio_y) + i_y*ratio_y )
        acum = acum + matrix(index_bs + i_bs, index_bl + i_bl, index_y + i_y)*weight
!        if (indicator_tirar>0) then
!            WRITE(nout, '(F12.8, X, I4, X, F12.8)')&
!            i_y*ratio_y, i_y,ratio_y
!        end if
      end do
   end do
end do

!f_value = acum
interpolate3 = acum
end function
    
DOUBLE PRECISION function interpolate3_finer(bs, bl, y, matrix)
USE param
DOUBLE PRECISION, INTENT(IN) :: bs, bl, y, matrix(b_num_LC_finer, b_num_long_finer, y_num_finer)
DOUBLE PRECISION ::  slope, weight, acum, ratio_bs, ratio_bl, ratio_y, eps
INTEGER :: index_bs, index_bl, index_y, i_bs, i_bl, i_y

eps = 1d-8 !NEEDED TO ADD EPSILON BECAUSE IF y = y_grid(i) IT MAY CHOOSE index_y = i-1 DUE TO APPROX ERROR IN THE ROUNDING
index_y = (MAX(MIN(INT((y_num_finer-1)*(y+eps-y_grid_finer(1))/(y_grid_finer(y_num_finer)-y_grid_finer(1)))+1,y_num_finer-1), 1))
index_bs = (MAX(MIN(INT((b_num_LC_finer-1)*(bs+eps-b_grid_LC_finer(1))/(b_grid_LC_finer(b_num_LC_finer)-b_grid_LC_finer(1)))+1,b_num_LC_finer-1), 1))
index_bl = (MAX(MIN(INT((b_num_long_finer -1)*(bl+eps-b_grid_long_finer(1)) /(b_grid_long_finer(b_num_long_finer)-b_grid_long_finer(1)))+1,b_num_long_finer-1), 1))



ratio_bs = (bs - b_grid_LC_finer(index_bs)) / (b_grid_LC_finer(index_bs+1) - b_grid_LC_finer(index_bs))
ratio_bl = (bl - b_grid_long_finer(index_bl)) / (b_grid_long_finer(index_bl+1) - b_grid_long_finer(index_bl))
ratio_y = (y - y_grid_finer(index_y)) / (y_grid_finer(index_y+1) - y_grid_finer(index_y))
acum=0

do i_bs=0,1
   do i_bl = 0,1
      do i_y=0,1
        weight = ((1-i_bs)*(1 - ratio_bs) + i_bs*ratio_bs) * ((1-i_bl)*(1 - ratio_bl) + i_bl*ratio_bl) * &
                 ((1-i_y)*(1 - ratio_y) + i_y*ratio_y )
        acum = acum + matrix(index_bs + i_bs, index_bl + i_bl, index_y + i_y)*weight
!        if (indicator_tirar>0) then
!            WRITE(nout, '(F12.8, X, I4, X, F12.8)')&
!            i_y*ratio_y, i_y,ratio_y
!        end if
      end do
   end do
end do

!f_value = acum
interpolate3_finer = acum
end function    


DOUBLE PRECISION function q_paid_fun(b_LC, b_long, y, index_premium)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_LC, b_long, y
INTEGER, INTENT(IN) :: index_premium
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl, interpolate3
INTEGER :: index_y, NINTV
EXTERNAL DBS2VL, interpolate3


!!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
bs=  min(max(b_inf_LC+1d-6, b_LC),b_sup_LC-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)

q_paid_fun = interpolate3(bs, bl, y, q_nodef_matrix(:,:,:, index_premium))

end

DOUBLE PRECISION function q_LC_paid_fun(b_LC, b_long, y, index_premium)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_LC, b_long, y
INTEGER, INTENT(IN) :: index_premium
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl, interpolate3
INTEGER :: index_y, NINTV
EXTERNAL DBS2VL, interpolate3


!!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
bs=  min(max(b_inf_LC+1d-6, b_LC),b_sup_LC-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)

q_LC_paid_fun = interpolate3(bs, bl, y, q_nodef_LC_matrix(:,:,:, index_premium))

end
    
    
DOUBLE PRECISION function q_paid_def_fun(b_LC, b_long, y, index_premium)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_LC, b_long, y
INTEGER, INTENT(IN) :: index_premium
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl, interpolate3
INTEGER :: index_y, NINTV
EXTERNAL DBS2VL, interpolate3


!!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
bs=  min(max(b_inf_LC+1d-6, b_LC),b_sup_LC-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)

q_paid_def_fun = interpolate3(bs, bl, y, q_def_matrix(:,:,:, index_premium))
end

DOUBLE PRECISION function q_LC_def_fun(b_LC, b_long, y, index_premium)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_LC, b_long, y
INTEGER, INTENT(IN) :: index_premium
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl, interpolate3
INTEGER :: index_y, NINTV
EXTERNAL DBS2VL, interpolate3


!!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
bs=  min(max(b_inf_LC+1d-6, b_LC),b_sup_LC-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)

q_LC_def_fun = interpolate3(bs, bl, y, q_LC_def_matrix(:,:,:, index_premium))
end
    
    
DOUBLE PRECISION function EV_fun(b_LC, b_long, y, index_premium)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_LC, b_long, y
INTEGER, INTENT(IN) :: index_premium
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl, interpolate3_finer
INTEGER :: index_b, NINTV, index_y
EXTERNAL DBS2VL, interpolate3_finer

!!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
bs=  min(max(b_inf_LC+1d-6, b_LC),b_sup_LC-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)

EV_fun = interpolate3_finer(bs, bl, y, EV_matrix(:,:,:, index_premium))
end function


    
    
DOUBLE PRECISION function EV_excl_fun(b_LC, b_long, y, index_premium)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_LC, b_long, y
INTEGER, INTENT(IN) :: index_premium
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl, interpolate3_finer
INTEGER :: index_b, NINTV, index_y
EXTERNAL DBS2VL, interpolate3_finer

!!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
bs=  min(max(b_inf_LC+1d-6, b_LC),b_sup_LC-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)

EV_excl_fun = interpolate3_finer(bs, bl, y, EV_excl_matrix(:,:,:, index_premium))
end function

    


DOUBLE PRECISION function q_menu_fun(b_LC, b_long, y, index_premium)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_LC, b_long, y
INTEGER, INTENT(IN) :: index_premium
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl, interpolate3_finer
INTEGER :: index_b, NINTV, index_y
EXTERNAL DBS2VL, interpolate3_finer

!!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
bs=  min(max(b_inf_LC+1d-6, b_LC),b_sup_LC-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)
q_menu_fun = interpolate3_finer(bs, bl, y, q_menu_matrix(:,:,:, index_premium))
end


DOUBLE PRECISION function q_LC_menu_fun(b_LC, b_long, y, index_premium)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_LC, b_long, y
INTEGER, INTENT(IN) :: index_premium
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl, interpolate3_finer
INTEGER :: index_b, NINTV, index_y
EXTERNAL DBS2VL, interpolate3_finer

!!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
bs=  min(max(b_inf_LC+1d-6, b_LC),b_sup_LC-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)
q_LC_menu_fun = interpolate3_finer(bs, bl, y, q_LC_menu_matrix(:,:,:, index_premium))
end 
    
DOUBLE PRECISION function q_menu_def_fun(b_LC, b_long, y, index_premium)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_LC, b_long, y
INTEGER, INTENT(IN) :: index_premium
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl, interpolate3_finer
INTEGER :: index_b, NINTV, index_y
EXTERNAL DBS2VL, interpolate3_finer

!!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
bs=  min(max(b_inf_LC+1d-6, b_LC),b_sup_LC-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)
q_menu_def_fun = interpolate3_finer(bs, bl, y, q_menu_def_matrix(:,:,:, index_premium))
end

DOUBLE PRECISION function q_LC_menu_def_fun(b_LC, b_long, y, index_premium)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_LC, b_long, y
INTEGER, INTENT(IN) :: index_premium
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl, interpolate3_finer
INTEGER :: index_b, NINTV, index_y
EXTERNAL DBS2VL, interpolate3_finer

!!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
bs=  min(max(b_inf_LC+1d-6, b_LC),b_sup_LC-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)
q_LC_menu_def_fun = interpolate3_finer(bs, bl, y, q_LC_menu_def_matrix(:,:,:, index_premium))
end

DOUBLE PRECISION function pi_fun(b_LC_int, b_int, y, index_premium)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_LC_int, b_int, y
INTEGER, INTENT(IN) :: index_premium
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, b_LC, b, interpolate3
INTEGER :: index_y, NINTV
EXTERNAL DBS2VL, interpolate3
!NINTV = b_num - 1
index_y = (MAX(MIN(INT((y_num-1)*(y-y_grid(1))/(y_grid(y_num)-y_grid(1)))+1,y_num-1), 1))

!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
b_LC=  min(max(b_inf_LC+1d-6, b_LC_int),b_sup_LC-1d-6)
b=  min(max(b_inf+1d-6, b_int),b_sup-1d-6)

pi_fun = interpolate3(b_LC, b, y, pi_matrix(:,:,:, index_premium))
end
    
subroutine compute_knots
USE param
!INTEGER :: KORDER
!PARAMETER(KORDER=3)
!DOUBLE PRECISION, DIMENSION (b_num_long + KORDER) :: b_long_knots
!DOUBLE PRECISION, DIMENSION (b_num_LC + KORDER) :: b_LC_knots
EXTERNAL DBSNAK

CALL DBSNAK (b_num_long, b_grid_long, 3, b_long_knots)
CALL DBSNAK (b_num_LC, b_grid_LC,3, b_LC_knots)

CALL DBSNAK (b_num_long_finer, b_grid_long_finer, 3, b_long_knots_finer)
CALL DBSNAK (b_num_LC_finer, b_grid_LC_finer,3, b_LC_knots_finer)
! CALL DBSNAK (b_num_long, b_grid_long, KORDER, b_long_knots)
! CALL DBSNAK (b_num_LC, b_grid_LC, KORDER, b_LC_knots)
end subroutine


double precision function b_limit_fun(bnext)
use param
double precision, intent(in) :: bnext
double precision :: q_menu_fun
external q_menu_fun

!write(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8)') b_LC_next_global, bnext, y_initial
b_limit_fun =  q_menu_fun(b_LC_next_global, bnext, y_initial, index_premium_initial) - q_global
!write(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8)') bnext, q_menu_fun(b_LC_next_global, bnext, y_initial, i_m_initial), q_global, b_limit_fun
end function  
    
    
!COMPUTE MATRICES AND SPLINE COEFFICIENTS FOR q(bs, bl, y) and E[V(bs, bl, y') | y)]
subroutine compute_q_ev
USE param
INTEGER :: i_b_LC, i_b_long, i_y, i, ILEFT, IRIGHT, i_premium, index_premium_next, KXORD, KYORD, LDF, MAXFN
parameter (KXORD =3, KYORD = 3)
DOUBLE PRECISION :: b_LC, b_long, y_threshold, cdf_threshold, y_fun, DNORDF, exp_v_next, DNORIN, &
                    acum1, acum2, acum(2), cdf, cdf1, y_next, y_next1, value_next, v0_fun, v1_fun, exp_v_excl, acum_excl(2),&
                    exp_v_excl_ends, FDATA(b_num_LC_finer), DLEFT, DRIGHT, BREAK_GRID(b_num_LC_finer), m_next_LC,&
                    CSCOEF(4,b_num_LC_finer), DCSVAL, tirar, FDATA_2D_finer(b_num_LC_finer, b_num_long_finer), &
                    BSCOEF_finer(b_num_LC_finer, b_num_long_finer), q_menu_fun, q_vector(2), gamma_t,&
                    left, right, ERRREL, ERRABS, b_limit_fun, b_limit_fun_spline, exp_m(premium_num), exp_q(premium_num), exp_q_def(premium_num),&
                    exp_m_LC(premium_num), exp_q_LC(premium_num), ind_suspension_next, m_next, q_next_LC, q_paid_fun, q_LC_paid_fun, q_LC_def_fun, q_paid_def_fun,&
                    exp_q_def_LC(premium_num), b_LC_next_post_def, b_next_post_def, b_next_def, b_LC_next_def, b_LC_remain_def, b_remain_def, pi_fun

DOUBLE PRECISION, DIMENSION (b_num_LC_finer, y_num_finer, 2) :: blong_limit_matrix

EXTERNAL y_fun, DNORDF, v0_fun, v1_fun, DCSDEC, DNORIN, DBS2IN, q_menu_fun, DZBREN, q_paid_fun, q_LC_paid_fun, q_LC_def_fun, q_paid_def_fun, pi_fun



do i_premium = 1,premium_num
    index_premium_initial = i_premium
    do i_y = 1,y_num_finer
        !WRITE(nout, *) i_y
        y_initial = y_grid_finer(i_y)
        do i_b_LC = 1, b_num_LC_finer
            b_LC = b_grid_LC_finer(i_b_LC)
            do i_b_long = 1, b_num_long_finer
                b_long = b_grid_long_finer(i_b_long)
                exp_m = 0d+0
                exp_q = 0d+0
                exp_m_LC = 0d+0
                exp_q_LC = 0d+0
                do index_premium_next = 1,premium_num
                    gamma_t = alpha0 + alpha1*EXP(y_initial) !USED TO COMPUTE THE PRICING KERNEL
                    !COMPUTE EXPECTED VALUE FUNCTION FOR THE NEXT PERIOD IF IT CHOOSES A PORTFOLIO (b_LC, b_long)
                    !WHEN CURRENT INCOME = y_initial
                    y_threshold   = y_fun(b_LC, b_long, index_premium_next)
                    cdf_threshold = DNORDF((y_threshold - rho*y_initial - (1-rho)* mean_y)/std_eps)
                    exp_v_next = 0d+0
                    acum1=0d+0
                    acum2=0d+0
                    if (cdf_threshold <= cdf_inf) then
                        do i=1,quad_num
                            cdf = 0.5*(quad_x(i) + 1)*(cdf_sup - cdf_inf) + cdf_inf
                            y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
                            m_next = EXP(-r - premium_grid(i_premium) * (y_next - ( (1-rho)*mean_y + rho* y_initial)) - 0.5d+0*(premium_grid(i_premium)**2d+0)*std_eps**2d+0)
                            m_next_LC =  EXP(-r)*1d+0/(1d+0 + pi_fun(b_LC, b_long, y_next, index_premium_next))

                            value_next = v0_fun(b_LC, b_long, y_next, index_premium_next)
                            acum1 = acum1 + 0.5*(cdf_sup - cdf_inf)*quad_w(i) * value_next * trans_matrix(i_premium, index_premium_next)
                    
                            exp_m(index_premium_next) = exp_m(index_premium_next) + m_next*coupon*0.5* (cdf_sup - cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            exp_q(index_premium_next) = exp_q(index_premium_next) + (1d+0-delta)*m_next*q_paid_fun(b_LC, b_long, y_next, index_premium_next)*0.5*(cdf_sup-cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            exp_m_LC(index_premium_next) = exp_m_LC(index_premium_next) + m_next_LC*coupon* 0.5* (cdf_sup - cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)

                            q_next_LC = (1d+0 - delta)*q_LC_paid_fun(b_LC, b_long, y_next, index_premium_next)
                            exp_q_LC(index_premium_next) = exp_q_LC(index_premium_next) + m_next_LC*q_next_LC* 0.5* (cdf_sup - cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            
                        end do
                    elseif(cdf_threshold >= cdf_sup) then
                        do i=1,quad_num
                            cdf = 0.5*(quad_x(i) + 1)*(cdf_sup - cdf_inf) + cdf_inf
          	                y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
                            m_next = EXP(-r - premium_grid(i_premium) * (y_next - ( (1-rho)*mean_y + rho* y_initial)) - 0.5d+0*(premium_grid(i_premium)**2d+0)*std_eps**2d+0)
                            m_next_LC =  EXP(-r)*1d+0/(1d+0 + pi_fun(b_LC, b_long, y_next, index_premium_next))

          	                value_next = v1_fun(b_LC, b_long, y_next, index_premium_next)
             	            acum1 = acum1 + 0.5*(cdf_sup - cdf_inf)*quad_w(i) * value_next*trans_matrix(i_premium, index_premium_next)
                            
                            exp_q(index_premium_next) = exp_q(index_premium_next) + m_next* q_paid_def_fun(b_LC, b_long, y_next, index_premium_next) *0.5*(cdf_sup-cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            q_next_LC = q_LC_def_fun(b_LC, b_long, y_next, index_premium_next)
                            exp_q_LC(index_premium_next) = exp_q_LC(index_premium_next) + m_next_LC*q_next_LC* 0.5* (cdf_sup - cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                        end do
                    else
                        do i=1,quad_num
                            cdf = 0.5*(quad_x(i) + 1)*(cdf_threshold - cdf_inf) + cdf_inf
                            y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
                            value_next = v1_fun(b_LC, b_long, y_next, index_premium_next)
                            acum1 = acum1 + 0.5*(cdf_threshold - cdf_inf)*quad_w(i) * value_next*trans_matrix(i_premium, index_premium_next)
                            m_next = EXP(-r - premium_grid(i_premium) * (y_next - ( (1-rho)*mean_y + rho* y_initial)) - 0.5d+0*(premium_grid(i_premium)**2d+0)*std_eps**2d+0)
                            m_next_LC =  EXP(-r)*1d+0/(1d+0 + pi_fun(b_LC, b_long, y_next, index_premium_next))
                            
                            exp_q(index_premium_next) = exp_q(index_premium_next) + m_next *q_paid_def_fun(b_LC, b_long, y_next, index_premium_next) * 0.5* (cdf_threshold-cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            q_next_LC = q_LC_def_fun(b_LC, b_long, y_next, index_premium_next)
                            exp_q_LC(index_premium_next) = exp_q_LC(index_premium_next) + m_next_LC*q_next_LC* 0.5* (cdf_threshold-cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            
                            cdf1 = 0.5*(quad_x(i) + 1)*(cdf_sup - cdf_threshold) + cdf_threshold
                            y_next1 = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf1) * std_eps 
                            value_next = v0_fun(b_LC, b_long, y_next1, index_premium_next)
                            acum2 = acum2 + 0.5*(cdf_sup - cdf_threshold)*quad_w(i) * value_next*trans_matrix(i_premium, index_premium_next)
                   
                            m_next = EXP(-r - premium_grid(i_premium) * (y_next1 - ( (1-rho)*mean_y + rho* y_initial)) - 0.5d+0*(premium_grid(i_premium)**2d+0)*std_eps**2d+0)
                            m_next_LC =  EXP(-r)*1d+0/(1d+0 + pi_fun(b_LC, b_long, y_next1, index_premium_next))
                            exp_m(index_premium_next) = exp_m(index_premium_next) + m_next*coupon*0.5* (cdf_sup - cdf_threshold)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            exp_q(index_premium_next) = exp_q(index_premium_next) + (1d+0-delta)*m_next*q_paid_fun(b_LC, b_long, y_next1, index_premium_next)*0.5*(cdf_sup - cdf_threshold)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            exp_m_LC(index_premium_next)   = exp_m_LC(index_premium_next) + m_next_LC*coupon* 0.5* (cdf_sup - cdf_threshold)*quad_w(i)*trans_matrix(i_premium, index_premium_next)

                            q_next_LC = (1d+0 - delta)*q_LC_paid_fun(b_LC, b_long, y_next1, index_premium_next)
                            exp_q_LC(index_premium_next) = exp_q_LC(index_premium_next) + m_next_LC*q_next_LC* 0.5* (cdf_sup - cdf_threshold)*quad_w(i)*trans_matrix(i_premium, index_premium_next)

                        end do
                    end if
                    acum(index_premium_next) =  (acum1+acum2)/(cdf_sup - cdf_inf)
                    exp_m(index_premium_next) = exp_m(index_premium_next)/(cdf_sup - cdf_inf)
                    exp_q(index_premium_next) = exp_q(index_premium_next)/(cdf_sup - cdf_inf)
                    exp_m_LC(index_premium_next) = exp_m_LC(index_premium_next)/(cdf_sup - cdf_inf)
                    exp_q_LC(index_premium_next) = exp_q_LC(index_premium_next)/(cdf_sup - cdf_inf)
                end do
                EV_matrix(i_b_LC, i_b_long, i_y, i_premium) = SUM(acum)
                q_menu_matrix(i_b_LC, i_b_long, i_y, i_premium) = MAX(SUM(exp_m + exp_q), 0d+0) 
                q_LC_menu_matrix(i_b_LC, i_b_long, i_y, i_premium) = MAX(SUM(exp_m_LC + exp_q_LC), 0d+0)

                exp_q_def = 0d+0
                exp_q_def_LC = 0d+0
                exp_m = 0d+0
                exp_q = 0d+0
                exp_m_LC = 0d+0
                exp_q_LC = 0d+0
                acum = 0d+0
                acum_excl = 0d+0
                do index_premium_next = 1,premium_num
                    y_threshold = y_fun((1+r)*recovery*b_LC, (1+r)*recovery*b_long, index_premium_next)
                    cdf_threshold = DNORDF((y_threshold - rho*y_initial - (1-rho)* mean_y)/std_eps)
                    b_next_post_def = recovery*b_long*(1d+0 +r)
                    b_LC_next_post_def = recovery*b_LC*(1d+0 +r)
                    b_next_def = b_next_post_def !*(1d+0 +r)
                    b_LC_next_def = b_LC_next_post_def !*(1d+0 +r)
                    b_remain_def = (1+r)*b_long
                    b_LC_remain_def = (1+r)*b_LC
                    !write(6, '(I4, X, I4, X, I4, X, F12.8, X, F12.8, X, F12.8)') i_b_LC, i_b_long, i_y, b_next_post_def, b_next_def, b_remain_def
                    if (cdf_threshold <= cdf_inf) then !THRESHOLD TOMORROW IS LOW
                        do i=1,quad_num
                            cdf = 0.5*(quad_x(i) + 1)*(cdf_sup - cdf_inf) + cdf_inf
                            y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
                            m_next = EXP(-r - premium_grid(i_premium) * ( y_next - ( (1-rho)*mean_y + rho* y_initial)) - 0.5d+0*(premium_grid(i_premium)**2d+0)*std_eps**2d+0)
                            m_next_LC =  EXP(-r)*1d+0/(1d+0 + pi_fun(b_LC_next_post_def, b_next_post_def, y_next, index_premium_next))
                            !GET b_recov_local BONDS FOR EACH BOND IN DEFAULT
                            exp_m(index_premium_next) = exp_m(index_premium_next) + m_next *recovery*(1+r)*coupon * 0.5* (cdf_sup - cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            exp_q(index_premium_next) = exp_q(index_premium_next) + (1d+0-delta) * m_next*recovery*(1+r)*&
                                q_paid_fun(b_LC_next_post_def, b_next_post_def, y_next, index_premium_next) *0.5*(cdf_sup-cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                
                            exp_m_LC(index_premium_next) = exp_m_LC(index_premium_next) + m_next_LC *recovery*(1+r)*coupon*0.5* (cdf_sup - cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            q_next_LC = (1d+0 - delta)*q_LC_paid_fun(b_LC_next_post_def, b_next_post_def, y_next, index_premium_next)
                            exp_q_LC(index_premium_next) = exp_q_LC(index_premium_next) + m_next_LC*recovery*(1+r)*q_next_LC* 0.5* (cdf_sup - cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                
                            exp_q_def(index_premium_next) = exp_q_def(index_premium_next) + m_next* (1d+0 +r)* &
                            q_paid_def_fun(b_LC_remain_def, b_remain_def, y_next, index_premium_next) *0.5*(cdf_sup-cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            m_next_LC =  EXP(-r)*1d+0/(1d+0 + pi_fun(b_LC_remain_def, b_remain_def, y_next, index_premium_next))
                            q_next_LC = q_LC_def_fun(b_LC_remain_def, b_remain_def, y_next, index_premium_next)
                            exp_q_def_LC(index_premium_next) = exp_q_def_LC(index_premium_next) + m_next_LC* (1d+0 +r)*q_next_LC* 0.5* (cdf_sup - cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)                        
                        end do
                    elseif (cdf_threshold > cdf_sup) then !THRESHOLD TOMORROW IS HIGH
                        do i=1,quad_num
                            cdf = 0.5*(quad_x(i) + 1)*(cdf_sup - cdf_inf) + cdf_inf
                            y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
                            m_next = EXP(-r - premium_grid(i_premium) * ( y_next - ( (1-rho)*mean_y + rho* y_initial)) - 0.5d+0*(premium_grid(i_premium)**2d+0)*std_eps**2d+0)
                            m_next_LC =  EXP(-r)*1d+0/(1d+0 + pi_fun(b_LC_next_def, b_next_def, y_next, index_premium_next))
                            
                            exp_q(index_premium_next) = exp_q(index_premium_next) +m_next*recovery*(1d+0 +r)*&
                                q_paid_def_fun(b_LC_next_def, b_next_def, y_next, index_premium_next) *0.5*(cdf_sup-cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                
                            q_next_LC = q_LC_def_fun(b_LC_next_def, b_next_def, y_next, index_premium_next)
                            exp_q_LC(index_premium_next) = exp_q_LC(index_premium_next) + m_next_LC*recovery*(1d+0 +r)*q_next_LC* 0.5* (cdf_sup - cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                
                            exp_q_def(index_premium_next) = exp_q_def(index_premium_next) + m_next* (1d+0 +r)* &
                            q_paid_def_fun(b_LC_remain_def, b_remain_def, y_next, index_premium_next) *0.5*(cdf_sup-cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            
                            m_next_LC =  EXP(-r)*1d+0/(1d+0 + pi_fun(b_LC_remain_def, b_remain_def, y_next, index_premium_next))
                            q_next_LC = q_LC_def_fun(b_LC_remain_def, b_remain_def, y_next, index_premium_next)
                            exp_q_def_LC(index_premium_next) = exp_q_def_LC(index_premium_next) + m_next_LC* (1d+0 +r)*q_next_LC* 0.5* (cdf_sup - cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                        end do
                    else
                        do i=1,quad_num
                            !compute integral for y < y_threshold
                            cdf = 0.5*(quad_x(i) + 1)*(cdf_threshold - cdf_inf) + cdf_inf
                            y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
                            m_next = EXP(-r - premium_grid(i_premium) * ( y_next - ( (1-rho)*mean_y + rho* y_initial)) - 0.5d+0*(premium_grid(i_premium)**2d+0)*std_eps**2d+0)
                            m_next_LC =  EXP(-r)*1d+0/(1d+0 + pi_fun(b_LC_next_def, b_next_def, y_next, index_premium_next))
                            exp_q(index_premium_next) = exp_q(index_premium_next) +m_next*recovery*(1d+0 +r)*&
                                q_paid_def_fun(b_LC_next_def, b_next_def, y_next, index_premium_next) *0.5*(cdf_threshold - cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                
                            q_next_LC = q_LC_def_fun(b_LC_next_def, b_next_def, y_next, index_premium_next)
                            exp_q_LC(index_premium_next) = exp_q_LC(index_premium_next) +m_next_LC*recovery*(1d+0 +r)*q_next_LC* 0.5* (cdf_threshold - cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            
                            exp_q_def(index_premium_next) = exp_q_def(index_premium_next) + m_next* (1d+0 +r)* &
                            q_paid_def_fun(b_LC_remain_def, b_remain_def, y_next, index_premium_next) *0.5*(cdf_threshold - cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            
                            m_next_LC =  EXP(-r)*1d+0/(1d+0 + pi_fun(b_LC_remain_def, b_remain_def, y_next, index_premium_next))
                            q_next_LC = q_LC_def_fun(b_LC_remain_def, b_remain_def, y_next, index_premium_next)
                            exp_q_def_LC(index_premium_next) = exp_q_def_LC(index_premium_next) + m_next_LC* (1d+0 +r)*q_next_LC* 0.5* (cdf_threshold - cdf_inf)*quad_w(i)*trans_matrix(i_premium, index_premium_next)

                            !compute integral for y > y_threshold
                            cdf = 0.5*(quad_x(i) + 1)*(cdf_sup - cdf_threshold) + cdf_threshold
                            y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
                            m_next = EXP(-r - premium_grid(i_premium) * ( y_next - ( (1-rho)*mean_y + rho* y_initial)) - 0.5d+0*(premium_grid(i_premium)**2d+0)*std_eps**2d+0)
                            m_next_LC =  EXP(-r)*1d+0/(1d+0 + pi_fun(b_LC_next_def, b_next_def, y_next, index_premium_next))
                            exp_m(index_premium_next) = exp_m(index_premium_next) + m_next *recovery*(1+r)*coupon * 0.5* (cdf_sup - cdf_threshold)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            exp_q(index_premium_next) = exp_q(index_premium_next) + (1d+0-delta) * m_next*recovery*(1+r)*&
                                q_paid_fun(b_LC_next_post_def, b_next_post_def, y_next, index_premium_next) *0.5*(cdf_sup - cdf_threshold)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                
                            exp_m_LC(index_premium_next) = exp_m_LC(index_premium_next) + m_next_LC *recovery*(1+r)*coupon* 0.5* (cdf_sup - cdf_threshold)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            q_next_LC = (1d+0 - delta)*q_LC_paid_fun(b_LC_next_post_def, b_next_post_def, y_next, index_premium_next) 
                            exp_q_LC(index_premium_next) = exp_q_LC(index_premium_next) + m_next_LC*recovery*(1d+0 +r)*q_next_LC* 0.5* (cdf_sup - cdf_threshold)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                
                            !COMPUTE NEXT PERIOD BOND PRICE WHEN ECONOMY KEEP EXCLUDED
                            exp_q_def(index_premium_next) = exp_q_def(index_premium_next) + m_next* (1d+0 +r)* &
                            q_paid_def_fun(b_LC_remain_def, b_remain_def, y_next, index_premium_next) *0.5*(cdf_sup-cdf_threshold)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            
                            m_next_LC =  EXP(-r)*1d+0/(1d+0 + pi_fun(b_LC_remain_def, b_remain_def, y_next, index_premium_next))
                            q_next_LC = q_LC_def_fun(b_LC_remain_def, b_remain_def, y_next, index_premium_next) 
                            exp_q_def_LC(index_premium_next) = exp_q_def_LC(index_premium_next) + m_next_LC* (1d+0 +r)*q_next_LC* 0.5* (cdf_sup - cdf_threshold)*quad_w(i)*trans_matrix(i_premium, index_premium_next)
                            !write(6, '(F12.8, X, F12.8, X, F12.8, X, F12.8)') y_next, m_next, q_paid_fun(b_LC_next, b_long_next, y_next, index_premium_next), exp_q(index_premium_next)
                        end do
                    end if
                    !NEED TO ADJUST EXPECTATION TO TAKE INTO ACCOUNT THAT
                    !THE ACTUAL PROBABILITY MASS BETWEEN y_threshold and infty IS UNDERESTIMATED USING LEBESGE-..
                    exp_q(index_premium_next) =  exp_q(index_premium_next) / (cdf_sup - cdf_inf)
                    exp_m(index_premium_next) = exp_m(index_premium_next) / (cdf_sup - cdf_inf)
                    exp_q_LC(index_premium_next) = exp_q_LC(index_premium_next) / (cdf_sup - cdf_inf)
                    exp_m_LC(index_premium_next) = exp_m_LC(index_premium_next) / (cdf_sup - cdf_inf)
                    exp_q_def(index_premium_next) = exp_q_def(index_premium_next)/ (cdf_sup - cdf_inf)
                    exp_q_def_LC(index_premium_next) = exp_q_def_LC(index_premium_next)/ (cdf_sup - cdf_inf)
                    !write(6, '(I3, X, F12.8, X, F12.8, X, F12.8, X, F12.8)') index_premium_next, exp_q(index_premium_next), exp_q_LC(index_premium_next),  exp_q_def(index_premium_next), acum(index_premium_next)
                    !end do
                    exp_v_excl      = 0d+0
                    exp_v_excl_ends = 0d+0
                    do i=1,quad_num
                        !CONTINUES EXCLUDED IN THE NEXT PERIOD (NO BORROWING AND DEFAULT COST)
                        cdf = 0.5*(quad_x(i) + 1)*(cdf_sup - cdf_inf) + cdf_inf
                        y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
                        value_next = v1_fun((1+r)*b_LC, (1+r)*b_long, y_next, index_premium_next)
                        exp_v_excl = exp_v_excl + 0.5*(cdf_sup - cdf_inf)*quad_w(i) * value_next*trans_matrix(i_premium, index_premium_next)
                
                        !IF EXCLUSION ENDS IN THE NEXT PERIOD ==> LEAVE DEFAULT WITH Recovery DEBT
                        value_next = max(v0_fun(recovery*(1+r)*b_LC, recovery*(1+r)*b_long, y_next, index_premium_next), v1_fun(recovery*(1+r)*b_LC, recovery*(1+r)*b_long, y_next, index_premium_next))
                        exp_v_excl_ends = exp_v_excl_ends + 0.5*(cdf_sup - cdf_inf)*quad_w(i) * value_next*trans_matrix(i_premium, index_premium_next)
                  
                    end do
                    acum_excl(index_premium_next) = ((1d+0 - prob_excl_end) * exp_v_excl + prob_excl_end * exp_v_excl_ends)/(cdf_sup - cdf_inf)
                end do
                EV_excl_matrix(i_b_LC, i_b_long, i_y, i_premium) = SUM(acum_excl)
                q_menu_def_matrix(i_b_LC, i_b_long, i_y, i_premium) = MIN(MAX(SUM(prob_excl_end*(exp_m + exp_q) + (1d+0 - prob_excl_end)*exp_q_def), 0d+0), EXP(-r))
                q_LC_menu_def_matrix(i_b_LC, i_b_long, i_y, i_premium) = MIN(MAX(SUM(prob_excl_end*(exp_m_LC + exp_q_LC) + (1d+0 - prob_excl_end)*exp_q_def_LC), 0d+0), EXP(-r))
            end do
            !WRITE(nout, *) i_y, i_b_LC, EV_excl_matrix(i_b_LC, i_y, i_premium)
        end do
    end do
end do



end subroutine




subroutine iterate
USE param
DOUBLE PRECISION :: y_valor, b_valor, b_valor_def, v_valor, convergence, criteria, deviation, q_fun, b_next_LC,&
                    b_next_long, g, b0_next(3), b1_next(3), b, q_value, v_valor_exl, FVALUE, q_menu_fun, q_LC_menu_fun, &
                    v0_value, v1_value, q, acum_v, acum_excl, dev_q, dev_pi, dev_v, dev_debt, dev_LC_debt, objective_function_excl, &
                    FDATA(b_num_LC), DLEFT, DRIGHT, BREAK_GRID(b_num_LC), CSCOEF(4,b_num_LC), b_tirar, DCSVAL,&
                    FDATA_2D(b_num_LC,b_num_long), BSCOEF(b_num_LC, b_num_long), foc_vector(2), disut_default, q_vector(2), q_LC_menu_def_fun, q_menu_def_fun 

DOUBLE PRECISION, DIMENSION(b_num_LC, b_num_long,  y_num, theta_num) :: v0_matrix_new, v1_matrix_new, q_nodef_matrix_new, q_nodef_LC_matrix_new, &
                                                                   output0_matrix, output1_matrix, &
                                                                   c0_matrix, c1_matrix, g0_matrix_new, g1_matrix_new, q_def_matrix_new, q_LC_def_matrix_new

DOUBLE PRECISION, DIMENSION(b_num_LC, b_num_long,  y_num, theta_num) :: v_matrix_new, dev_matrix, q_paid_matrix_new, dev_matrix_q, &
                             b0_next_LC_matrix_new, b0_next_long_matrix_new, b1_next_LC_matrix_new, pi_matrix_new, b1_next_long_matrix_new
INTEGER :: d, i, i_b_LC, i_b_long, i_y, i_premium, i_def_opt, i_b_zero, i_b_optimal_LC0, i_b_optimal_long0, &
           i_b_optimal_LC1, i_b_optimal_long1, i_b_next_LC, i_b_next_long, ILEFT, IRIGHT,NINTV, N, &
           KXORD, KYORD, LDF, counter_final
parameter (KXORD =3, KYORD = 3)
EXTERNAL q_fun, DCSDEC, DCSVAL, objective_function_excl, q_menu_fun, DBS2IN, q_LC_menu_fun, q_LC_menu_def_fun, q_menu_def_fun



i_b_zero = b_num_long
criteria = 5d-5   !CRITERIA FOR CONVERGENCE
convergence = -1
indicator_global_search = 1
counter_final = 0
do WHILE(convergence<0)
    deviation = 0
    dev_q =0
    dev_v = 0
    dev_LC_debt = 0
    dev_debt = 0
    dev_pi = 0
    counter = counter + 1
    counter_final = counter_final + 1
    q_global =  kappa* coupon*(1d+0-((1d+0-delta)*exp(-r))**counter)/(r+delta)
    call compute_q_ev !COMPUTE MATRICES OF SPLINE COEFFICIENTS TO APPROXIMATE EXPECTED VALUE FUNCTIONS AND BOND PRICE FUNCTION
    do i_premium = 1, premium_num
        index_premium_initial = i_premium
        !g_initial = g_grid(i_premium)
        do i_y = 1,y_num
            i_y_global = i_y
            y_initial = y_grid(i_y)
            !WRITE(nout, *) i_premium, i_y
            do i_b_LC =  1, b_num_LC
                i_b_LC_global = i_b_LC
                b_LC_initial = b_grid_LC(i_b_LC)
                do i_b_long = 1, b_num_long
                    !write(nout, *) i_premium, i_y, i_b_LC, i_b_long
                    i_b_long_global = i_b_long
                    b_long_initial = b_grid_long(i_b_long)
                    i_default_global=2 !Country defaults
                    b1_next(1) = b_LC_initial!
                    b1_next(2) = b_long_initial
                    v1_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium) = objective_function_excl(b1_next)! - disut_default 
                    !write(nout, *) v_valor !v1_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium) 
    !                      write(nout, '(I4, X, I4, X, I4, X, I4, X, F20.8, X, F12.8, X, F12.8)') i_b_LC, i_b_long, i_y, i_premium, v1_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium) !,&
                    b1_next_LC_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium) = b1_next(1)
                    b1_next_long_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium) = b1_next(2)
                    c1_matrix(i_b_LC, i_b_long, i_y_global, i_premium) = c_global
                    output1_matrix(i_b_LC, i_b_long, i_y_global, i_premium) = output_global
                    indicator_tirar = 0

                    i_default_global=1   ! Country does not default and is not excluded for the next period
                    i_excl_global=1
                    call optimize(b0_next, v_valor)
                    indicator_tirar = 0
                    v0_matrix_new(i_b_LC, i_b_long, i_y, i_premium) = v_valor
                    b0_next_LC_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium) = b0_next(1)
                    b0_next_long_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium)  = b0_next(2)
                    v0_matrix_new(i_b_LC, i_b_long, i_y, i_premium) = v_valor
                    c0_matrix(i_b_LC, i_b_long, i_y_global, i_premium) = c_global
                    output0_matrix(i_b_LC, i_b_long, i_y_global, i_premium) = output_global
                    q_nodef_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium) = q_menu_fun(b0_next(1), b0_next(2), y_initial, i_premium)
                    q_nodef_LC_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium) = q_LC_menu_fun(b0_next(1), b0_next(2), y_initial, i_premium)
                    pi_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium) = b0_next(3) !pi_global
                    indicator_tirar =0
                    if (v1_matrix_new(i_b_LC, i_b_long, i_y, i_premium) <= v0_matrix_new(i_b_LC, i_b_long, i_y, i_premium)) then
                        default_decision(i_b_LC, i_b_long, i_y, i_premium) = 1
                        b_next_matrix(i_b_LC, i_b_long,  i_y, i_premium, 1) = b0_next(1)  !SAVINGS IF NO DEFAULT
                        b_next_matrix(i_b_LC, i_b_long,  i_y, i_premium, 2) = b0_next(2) !SAVINGS IF NO DEFAULT
                        q_paid_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium) =  q_nodef_LC_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium)
                    else
                        default_decision(i_b_LC, i_b_long, i_y, i_premium) = 2
                        b_next_matrix(i_b_LC, i_b_long,  i_y, i_premium, 1) = b1_next(1)  !SAVINGS IF NO DEFAULT
                        b_next_matrix(i_b_LC, i_b_long,  i_y, i_premium, 2) = b1_next(2)  !SAVINGS IF NO DEFAULT
                        q_paid_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium) =  0d+0
                    end if
                    v_matrix_new(i_b_LC, i_b_long, i_y, i_premium) = MAX(v1_matrix_new(i_b_LC, i_b_long, i_y, i_premium), v0_matrix_new(i_b_LC, i_b_long, i_y, i_premium))


                    dev_matrix(i_b_LC, i_b_long, i_y, i_premium) = ABS(v_matrix_new(i_b_LC, i_b_long, i_y, i_premium) - v_matrix(i_b_LC, i_b_long, i_y, i_premium))
                    !WRITE(nout, '(F12.8)') dev_matrix(i_b_LC, i_b_long, i_y, i_premium)
                    indicator_tirar=0
                    dev_matrix_q(i_b_LC, i_b_long, i_y, i_premium) = min(100d+0, ABS(q_paid_matrix_new(i_b_LC, i_b_long, i_y, i_premium) - q_paid_matrix(i_b_LC, i_b_long, i_y, i_premium)))
                    dev_q = MAX(dev_q, ABS(q_paid_matrix_new(i_b_LC, i_b_long, i_y, i_premium) - q_paid_matrix(i_b_LC, i_b_long, i_y, i_premium)))
                    dev_pi = MAX(dev_pi, ABS(pi_matrix_new(i_b_LC, i_b_long, i_y, i_premium) - pi_matrix(i_b_LC, i_b_long, i_y, i_premium)))
                    dev_v = MAX(dev_v, ABS(v_matrix_new(i_b_LC, i_b_long, i_y, i_premium) - v_matrix(i_b_LC, i_b_long, i_y, i_premium)))
                    dev_debt = MAX(dev_debt, ABS(b0_next_long_matrix_new(i_b_LC, i_b_long, i_y, i_premium) - b0_next_long_matrix(i_b_LC, i_b_long, i_y, i_premium)))
                    dev_LC_debt = MAX(dev_LC_debt, ABS(b0_next_LC_matrix_new(i_b_LC, i_b_long, i_y, i_premium) - b0_next_LC_matrix(i_b_LC, i_b_long, i_y, i_premium)))
                    deviation = MAX(deviation, MAX(dev_v, dev_q))

                end do
            end do
        end do
    end do

      
    WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)') deviation, dev_v, dev_q, dev_LC_debt, dev_debt, dev_pi


    open(100, FILE ='graphs\iteration.txt',POSITION ='append')
    write(100, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)') deviation, dev_q, dev_LC_debt, dev_debt, dev_pi, counter
    close(100)
    !3) SAVE RESULTS OF THE CURRENT ITERATION

    open (10, FILE='graphs\v.txt',STATUS='replace')
    open (11, FILE='graphs\default.txt',STATUS='replace')
    open (12, FILE='graphs\q.txt',STATUS='replace')
    open (13, FILE='graphs\b_next.txt',STATUS='replace')
    open (14, FILE='graphs\dev.txt',STATUS='replace')
    open (16, FILE='graphs\q_paid.txt',STATUS='replace')
    open (17, FILE='graphs\q_def.txt',STATUS='replace')
    open (18, FILE='graphs\b_guess.txt',STATUS='replace')
    open (19, FILE='graphs\c.txt',STATUS='replace')
    open (20, FILE='graphs\pi.txt',STATUS='replace')
    open (23, FILE='graphs\counter.txt',STATUS='replace')

    write(23, '(F16.8)') counter
    do i_b_LC = 1,b_num_LC
        do i_b_long = 1,b_num_long
            do i_y = 1,y_num
                do i_premium = 1,premium_num

                    i_y_global = i_y
                    y_initial = y_grid(i_y)
                    b_LC_initial = b_grid_LC(i_b_LC)
                    i_b_LC_global = i_b_LC

                    b_long_initial = b_grid_long(i_b_long)
                    i_b_long_global = i_b_long
                    indicator_tirar=0
                    WRITE(10, '(G15.6, X, G15.6, X, G15.6)') v_matrix_new(i_b_LC, i_b_long, i_y, i_premium), &
                    v0_matrix_new(i_b_LC, i_b_long, i_y, i_premium), v1_matrix_new(i_b_LC, i_b_long, i_y, i_premium)
                    WRITE(11, '(F6.2)') default_grid(default_decision(i_b_LC, i_b_long, i_y, i_premium))
                    d = default_decision(i_b_LC, i_b_long, i_y, i_premium)
                    WRITE(12, '(F15.8, X, F15.8)') q_LC_menu_fun(b_LC_initial, b_long_initial, y_initial, i_premium), q_menu_fun(b_LC_initial, b_long_initial, y_initial, i_premium)

                    WRITE(13, '(F15.11, X, F15.11)') b_next_matrix(i_b_LC, i_b_long, i_y, i_premium,1), b_next_matrix(i_b_LC, i_b_long, i_y, i_premium,2)
                    WRITE(14, '(F15.11, X, F15.11)') dev_matrix(i_b_LC, i_b_long, i_y, i_premium), dev_matrix_q(i_b_LC, i_b_long, i_y, i_premium)

                    q_def_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium) = q_menu_def_fun(b_LC_initial, b_long_initial, y_initial, i_premium)
                    q_LC_def_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium) = q_LC_menu_def_fun(b_LC_initial, b_long_initial, y_initial, i_premium)
                    WRITE(16, '(F15.8, X, F15.8, X, F15.8, X, F15.8, X, F15.8)') q_paid_matrix(i_b_LC, i_b_long, i_y, i_premium),&
                    q_nodef_matrix_new(i_b_LC, i_b_long, i_y, i_premium), q_nodef_LC_matrix_new(i_b_LC, i_b_long, i_y, i_premium)
                     
                    WRITE(17, '(F15.8, X, F15.8)') q_def_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium),&
                    q_LC_def_matrix_new(i_b_LC, i_b_long, i_y_global, i_premium)
                     
                    indicator_tirar=0

                    WRITE(18, '(G15.6, X, G15.6, X, G15.6)') b0_next_LC_matrix_new(i_b_LC, i_b_long, i_y, i_premium), &
                                        b0_next_long_matrix_new(i_b_LC, i_b_long, i_y, i_premium), &
                                        b1_next_LC_matrix_new(i_b_LC, i_b_long, i_y, i_premium)
 
                    WRITE(19, '(F15.8, X, F15.8)') c0_matrix(i_b_LC, i_b_long, i_y, i_premium), c1_matrix(i_b_LC, i_b_long, i_y, i_premium)
                    WRITE(20, '(G15.6)') pi_matrix_new(i_b_LC, i_b_long, i_y, i_premium)
                end do
            end do
        end do
    end do
    CLOSE(10)
    CLOSE(11)
    CLOSE(12)
    CLOSE(13)
    CLOSE(14)
    CLOSE(16)
    CLOSE(17)
    close(18)
    close(19)
    close(20)
    close(23)

    write(nout ,*) 'Iteration saved '
    !UPDATE VALUES OF MATRICES
    v0_matrix  = v0_matrix_new
    v1_matrix  = v1_matrix_new
    v_matrix   = v_matrix_new
    q_paid_matrix = q_paid_matrix_new
    q_nodef_matrix = q_nodef_matrix_new
    q_nodef_LC_matrix = q_nodef_LC_matrix_new
    q_def_matrix = q_def_matrix_new
    q_LC_def_matrix = q_LC_def_matrix_new
    b0_next_LC_matrix = b0_next_LC_matrix_new 
    b0_next_long_matrix  = b0_next_long_matrix_new
    b1_next_LC_matrix = b1_next_LC_matrix_new
    pi_matrix = pi_matrix_new
    b1_next_long_matrix  = b1_next_long_matrix_new
    !pause
    !PRINT*,'deviation =', deviation
    if (deviation < criteria .or. counter_final > 1000) then
        convergence =1
    end if
end do
!FINALLY, STORE VALUE FUNCTIONS, POLICY FUNCTIONS AND PRICES USING A FINER GRID.
!THIS HELPS TO VISUALIZE THE RESULTS.
end subroutine


subroutine simulate
USE param
integer :: period_num, i,j,k, gov_type, random_num, ivalue, sample_num, MAXFN, IOPT, num_samples_nodef,&
           index_observation, d_ss
parameter (sample_num = 250, period_num=501)
DOUBLE PRECISION :: random_matrix(period_num, sample_num, 2), random_vector(1:2*sample_num*period_num), &
                    z(period_num,sample_num), b_LC(period_num+1,sample_num), b_long(period_num+1,sample_num),&
                    q(period_num,sample_num), q_LC(period_num,sample_num), eps, v_def, v_no_def, b_next(3), q_fun, &
                    q_paid, b_zero, tb(period_num,sample_num), ind_payment_suspension, pi(period_num,sample_num),&
                    y(period_num,sample_num), STEP, BOUND, XACC, objective_function, b_next_guess,&
                    v_valor, v0_fun, v1_fun, DNORIN, q_risk_free, q_menu_fun, b_next_LC, &
                    acum_spread, acum_spread2, acum_debt, spread_value, acum_dev_spread2, acum_dev_spread,&
                    indicator_sample(sample_num), smoothing_parameter, labor_sim(period_num,sample_num), g_sim(period_num,sample_num),&
                    c_sim(period_num,sample_num), tax_sim(period_num,sample_num), recovery_rate, tfp, q_LC_menu_fun, q_LC_menu_def_fun, q_menu_def_fun
DOUBLE PRECISION, ALLOCATABLE :: mean_spread(:), mean_std_spread(:), mean_debt(:)
                    

INTEGER :: liquidity(period_num, sample_num), excl(period_num, sample_num), d(period_num, sample_num), i_premium_sim(period_num, sample_num)
EXTERNAL RNSET, DRNUN, DNORIN, q_fun, DUVMIF, objective_function, v0_fun, v1_fun, q_menu_fun, q_LC_menu_fun, q_LC_menu_def_fun, q_menu_def_fun

b_zero = 0d+0  !IF THE COUNTRY IS NOT EXCLUDED AFTER A DEFAULT EPISODE ==> BORROWS AS IF IT STARTS WITH ZERO DEBT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! NOTE !!!!!!
! z DENOTE UNDERLYING SHOCK TO THE ENDOWMENT
! y DENOTE REALIZED ENDOWMENT (EXP(z))
! CODE WAS WRITEN WITH y = SHOCK, SO mean_y ACTUALLY = E(z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



STEP = (b_grid_LC(2) - b_grid_LC(1))*0.5
BOUND = 10*y_sup
XACC = 1d-6
MAXFN = 1000

ivalue=139719
call RNSET(ivalue)       !INITIALIZE THE SEED SO THE SAME RANDOM NUMBERS ARE GENERATED IN EVERY REALIZATION
call DRNUN(2*period_num*sample_num, random_vector)
indicator_global_search = 1
counter = 1.0d+6
q_risk_free = coupon* (1d+0 - ((1d+0 - delta)/(1d+0 + r))**(counter-1d+0))/ (delta + r)

!open (UNIT=1, FILE="r_vector.txt", status = 'replace')
!   do i=1,period_num*sample_num*2
!      write(1,'(F12.8)') random_vector(i)
!      !READ
!      end do
!CLOSE(1)

!First column of random_matrix is used to generate transitory shocks.
!Second column of random_matrix is used to generate suddenstops.
do j=1,sample_num
    do i=1,period_num
       random_matrix(i,j,1)=random_vector((j-1)*period_num+i)
       random_matrix(i,j,2)=random_vector(period_num*sample_num + (j-1)*period_num+i)
    end do
end do


open (UNIT=21, FILE="graphs\data_sim.txt", status = 'replace')
open (UNIT=22, FILE="graphs\def_per.txt", status = 'replace')
open (UNIT=23, FILE="graphs\param.txt", status = 'replace')
open (UNIT=24, FILE="graphs\tirar.txt", status = 'replace')
open (UNIT=26, FILE="graphs\def_per_ss.txt", status = 'replace')

WRITE(23, '(I10)') period_num-1
WRITE(23, '(I10)') sample_num
CLOSE(23)



index_premium_initial = 1 !ALWAYS no premium TO MARKETS 
i_premium_sim = index_premium_initial 

do j=1,sample_num   !SOLVE FOR SAMPLE j
    WRITE(nout, *) j
    !Set initial values
    i = 1
    z(i,j) = mean_y
    b_LC(i,j) = b_inf_LC
    b_long(i,j) = b_inf
    d(i,j) = 1   !NO DEFAULT IN FIRST PERIOD
    i_default_global = 1
    excl(i,j) = 0
    b_LC_initial = b_LC(i,j)
    b_long_initial  = b_long(i,j)
    y_initial = z(i,j)
    call optimize(b_next, v_valor)
    b_LC(i+1,j) = b_next(1)
    b_long(i+1,j) = b_next(2)
    pi(i,j) = b_next(3)
    q(i,j) = q_menu_fun(b_next(1), b_next(2), z(i,j), index_premium_initial)
    q_LC(i,j) = q_LC_menu_fun(b_next(1), b_next(2), z(i,j), index_premium_initial)
    y(i,j) = EXP(z(i,j)) 
    !WHEN DOES NOT PAY COUPON OF COCOS BONDS, NUMBER OF COCOS BONDS INCREASE AT RATE exp(r)
    c_sim(i,j) = y(i,j) - coupon*b_long(i,j) + (b_long(i+1,j)- (1-delta)*b_long(i,j))*q(i,j) + &
        (b_LC(i+1,j)- b_LC(i,j)*(1-delta)) *q_LC(i,j) - coupon*b_LC(i,j)/(1d+0 + pi(i,j)) -absortion
    do i=2,period_num
        recovery_rate = recovery
        !epsilon = realization of standard gaussian * standard deviation
        eps = DNORIN(random_matrix(i,j,1)) * std_eps
        z(i,j) = rho*z(i-1,j) + (1-rho)*mean_y + eps   !current output = ex ante mean + epsilon
        y_initial = z(i,j)
        !trans_matrix(1,2) = min(prob_ss*exp(-slope_ss*y_initial - 0.5d+0*(std_y*slope_ss)**2d+0),1d+0) !PROBABILITY OF SWITCHING FROM LOW TO HIGH RISK PREMIUM IN THE NEXT PERIOD.
        !trans_matrix(1,1) = 1d+0 - trans_matrix(1,2)
        if (random_matrix(i,j,2) <=trans_matrix(i_premium_sim(i-1,j), 1) ) THEN  !TRANSIT TO STATE WITH LOW RISK PREMIUM
            i_premium_sim(i,j) = 1
        else !TRANSIT TO STATE WITH HIGH RISK PREMIUM
            i_premium_sim(i,j) = 2
        end if
        index_premium_initial = i_premium_sim(i,j)
        d_ss =0
       
        v_no_def = v0_fun(b_LC(i,j), b_long(i,j), y_initial, index_premium_initial)
        v_def    = v1_fun(b_LC(i,j), b_long(i,j), y_initial, index_premium_initial)
        if (excl(i-1,j) > 0) then !BORROWER WAS EXCLUDED IN THE PREVIOUS PERIOD
            if (random_matrix(i,j,2) <=prob_excl_end .and. v_no_def>v_def) THEN  !EXCLUSION ENDS and NO DEFAULT IS DECLARED
                excl(i,j) = 0
                d(i,j) = 1
                i_default_global = 1
                b_LC_initial = b_LC(i,j)
                b_long_initial = b_long(i,j)
                y_initial = z(i,j)
                call optimize(b_next, v_valor)
                y(i,j) = EXP(z(i,j)) 
                b_LC(i+1,j) = b_next(1)
                b_long(i+1,j) = b_next(2)
                q(i,j) = q_menu_fun(b_next(1), b_next(2), z(i,j), index_premium_initial)
                q_LC(i,j) = q_LC_menu_fun(b_next(1), b_next(2), z(i,j), index_premium_initial)
                c_sim(i,j) = y(i,j) - coupon*b_long(i,j) + (b_long(i+1,j)- (1-delta)*b_long(i,j))*q(i,j) + &
                        (b_LC(i+1,j)- b_LC(i,j)*(1-delta))*q_LC(i,j) - coupon*b_LC(i,j)/(1d+0 + pi(i,j))-absortion
                tb(i,j) = y(i,j) - c_sim(i,j)
                pi(i,j) = b_next(3)
            elseif (random_matrix(i,j,2) <=prob_excl_end .and. v_no_def<=v_def) THEN  !EXCLUSION ENDS and GOVERNMENT DEFAULTS AGAIN
                d(i,j) = 2
                excl(i,j) = 1
                y(i,j) = min(EXP(y_initial)*(1d+0 - d0) - d1*EXP(y_initial)**2d+0, exp(y_initial))
                !zero recovery of debt in default
                b_LC(i+1,j) = b_LC(i,j)*(1d+0 +r)*recovery
                b_long(i+1,j) = b_long(i,j)*(1d+0 +r)*recovery
                c_sim(i,j) = y(i,j) -absortion
                tb(i,j) = y(i,j) - c_sim(i,j) 
                q(i,j) = q_menu_def_fun(b_LC(i+1,j), b_long(i+1,j), z(i,j), index_premium_initial)
                q_LC(i,j) = q_LC_menu_def_fun(b_LC(i+1,j), b_long(i+1,j), z(i,j), index_premium_initial)
                pi(i,j) = zero
                if (i_premium_sim(i,j)==2 .AND. v0_fun(b_LC(i,j), b_long(i,j), y_initial, 1) > v1_fun(b_LC(i,j), b_long(i,j), y_initial, 1)) then
                   !WOULD NOT HAD DEFAULTED IF IT COULD HAVE RETAINED ACCESS TO MARKETS
                    WRITE(26, '(I7, X, I7)') i-1, j !DEFAULT CAUSED BY SUDDEN STOP ALONE
                    d_ss = 1 !DEFAULT IS CAUSED BY SS
                end if
            else !exclusion continues
                d(i,j) = 1
                excl(i,j) = 1
                y(i,j) = min(EXP(y_initial)*(1d+0 - d0) - d1*EXP(y_initial)**2d+0, exp(y_initial))
                !positive recovery of debt in default
                b_LC(i+1,j) = b_LC(i,j)*(1+r)
                b_long(i+1,j) = b_long(i,j)*(1+r)
                c_sim(i,j) = y(i,j) -absortion
                tb(i,j) = y(i,j) - c_sim(i,j)
                q(i,j) = q_menu_def_fun(b_LC(i+1,j), b_long(i+1,j), z(i,j), index_premium_initial)
                q_LC(i,j) = q_LC_menu_def_fun(b_LC(i+1,j), b_long(i+1,j), z(i,j), index_premium_initial)
                pi(i,j) = zero
            end if
        elseif(excl(i-1,j) ==0) THEN
            !WRITE(nout, '(I3, X, I3, F10.5, X, F10.5)') liquidity(i-1,j), i_premium_sim(i,j),trans_matrix(1,liquidity(i-1,j)),random_matrix(i,j,2)
            b_LC_initial = b_LC(i,j)
            b_long_initial = b_long(i,j)
            !WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8)') b_initial, y_initial, v_no_def, v_def
            if (v_def>v_no_def) then  !COUNTRY DEFAULTS
                d(i,j) = 2
                i_default_global = 2
                excl(i,j) = 1
                y(i,j) = min(EXP(y_initial)*(1d+0 - d0) - d1*EXP(y_initial)**2d+0, exp(y_initial))
                !positive recovery of debt in default
                b_LC(i+1,j) = recovery*b_LC(i,j)*(1+r)
                b_long(i+1,j) = recovery*b_long(i,j)*(1+r)
                c_sim(i,j) = y(i,j) -absortion
                tb(i,j) = y(i,j) - c_sim(i,j)
                q(i,j) = q_menu_def_fun(b_LC(i+1,j), b_long(i+1,j), z(i,j), index_premium_initial)
                q_LC(i,j) = q_LC_menu_def_fun(b_LC(i+1,j), b_long(i+1,j), z(i,j), index_premium_initial)
                pi(i,j) = zero
                WRITE(22, '(I7)') i-1 !NEED TO SUBSTRACT 1. REASON: files start saving data on period 2
                if (i_premium_sim(i,j)==2 .AND. v0_fun(b_LC(i,j), b_long(i,j), y_initial, 1) > v1_fun(b_LC(i,j), b_long(i,j), y_initial, 1)) then
                    !WOULD NOT HAD DEFAULTED IF IT COULD HAVE RETAINED ACCESS TO MARKETS
                    WRITE(26, '(I7, X, I7)') i-1, j !DEFAULT CAUSED BY SUDDEN STOP ALONE
                    d_ss = 1 !DEFAULT IS CAUSED BY SS
                end if
            else
                d(i,j) = 1   !COUNTRY DOES NOT DEFAULT
                excl(i,j) = 0
                i_default_global = 1
                call optimize(b_next, v_valor)
                y(i,j) = EXP(z(i,j))
                b_LC(i+1,j) = b_next(1)
                b_long(i+1,j) = b_next(2)
                q(i,j) = q_menu_fun(b_next(1), b_next(2), z(i,j), index_premium_initial)
                q_LC(i,j) = q_LC_menu_fun(b_next(1), b_next(2), z(i,j), index_premium_initial)
                !WHEN DOES NOT PAY COUPON OF COCOS BONDS, NUMBER OF COCOS BONDS INCREASE AT RATE exp(r)
                c_sim(i,j) = y(i,j) - coupon*b_long(i,j) + (b_long(i+1,j)- (1-delta)*b_long(i,j))*q(i,j) + &
                        (b_LC(i+1,j)- b_LC(i,j)*(1-delta))*q_LC(i,j) *coupon*b_LC(i,j)/(1d+0 + pi(i,j))-absortion
                tb(i,j) = y(i,j) - c_sim(i,j)
                pi(i,j) = b_next(3)
            end if
        end if
!                     1         2         3         4         5        6         7        8      9      10        11      12      13     14       15
        WRITE(21, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, I3, X, I3, X, I3, X, F10.6, X F10.6, X, I3, X F10.6, X F10.6)') &
 !           1           2            3            4          5              6              7       8           9            10            11            12      13           14
        LOG(y(i,j)), b_LC(i,j), b_long(i,j),  q(i,j),  q_LC(i,j), LOG(c_sim(i,j)), tb(i,j)/y(i,j), d(i,j), excl(i,j), i_premium_sim(i,j), z(i,j),  recovery_rate, d_ss, &
        trans_matrix(1,1), pi(i,j)

    end do
end do
CLOSE(21)
CLOSE(22)
CLOSE(24)
close(26)

end subroutine    



program main
include 'link_fnl_static.h'
USE param
DOUBLE PRECISION :: y_valor, b_valor, f_valor, q_fun, start_time, end_time, indicator_external, def, &
                    b_LC, b_long, consum, ind_payment_suspension, disut_default
INTEGER  i_b_LC, i_b_long, i_y, i_premium, i
EXTERNAL q_fun

trans_matrix(1,1) = 1d+0

premium_grid(1) = 0d+0
premium_grid(2) = premium_grid_value

write(6, '(A10, X, F17.4)') 'sigma = ', sigma
write(6, '(A10, X, F17.4)') 'beta =', beta
write(6, '(A10, X, F17.4)') 'd0 = ', d0
write(6, '(A10, X, F17.4)') 'd1 =', d1
write(6, '(A10, X, F17.4)') 'prob excl ends =', prob_excl_end
write(6, '(A10, X, F17.4)') 'delta =', delta
write(6, '(A10, X, F17.4)') 'hi premium =', premium_grid(2)
write(6, '(A10, X, F17.4)') 'Slope ss prob =', slope_ss

 
call cpu_time(start_time)
call compute_grid
call quadrature
call compute_knots

indicator_external = 0 !FROM EXTERNAL FILE
    
  if (indicator_external < 0.5) then
        counter = 0
        
         do i_premium = 1,premium_num
            index_premium_initial = i_premium
            ind_payment_suspension = default_grid(index_premium_initial) !INDICATOR THAT TAKES A VALUE OF ONE WHEN THERE IS PAYMENT SUSPENSION DUE TO HIGH g 
            do i_y = 1,y_num
                y_initial = y_grid(i_y)
                consum = min((EXP(y_initial)*(1d+0 - d0) - d1*EXP(y_initial)**2d+0), exp(y_initial))
                v1_matrix(:, :, i_y, i_premium) = (consum)**(1-sigma)/(1-sigma)! - disut_default
!                write(6, '(A12, X, F12.8, X, F12.8)') 'def ', output, (consum)**(1-sigma)/(1-sigma)
                do i_b_LC = 1,b_num_LC
                   b_LC_initial = b_grid_LC(i_b_LC)
                   b_LC = b_inf_LC
                   do i_b_long = 1,b_num_long
                      i_default_global = 1
                      b_long_initial = b_grid_long(i_b_long)
                      consum = EXP(y_initial) - coupon*b_long_initial - absortion
                      v0_matrix(i_b_LC, i_b_long, i_y, i_premium) = (consum)**(1-sigma)/(1-sigma)
                      v_matrix(i_b_LC, i_b_long, i_y, i_premium) = MAX(v1_matrix(i_b_LC, i_b_long, i_y, i_premium), v0_matrix(i_b_LC, i_b_long, i_y, i_premium))
                      pi_matrix(i_b_LC, i_b_long, i_y, i_premium) = 0d+0
                      q_nodef_LC_matrix(i_b_LC, i_b_long, i_y, i_premium) = 0d+0
                      q_nodef_matrix(i_b_LC, i_b_long, i_y, i_premium) = 0d+0
                      q_LC_def_matrix(i_b_LC, i_b_long, i_y, i_premium) = 0d+0
                      q_def_matrix(i_b_LC, i_b_long, i_y, i_premium) = 0d+0
                   end do
                end do
            end do
         end do
   else !READ DATA FROM EXTERNAL FILES

       open (10, FILE='graphs\v.txt')
       open (13, FILE='graphs\b_next.txt')
       open (16, FILE='graphs\q_paid.txt')
       open (17, FILE='graphs\q_def.txt')
       open (18, FILE='graphs\b_guess.txt')
       open (19, FILE='graphs\counter.txt')
       open (20, FILE='graphs\pi.txt')

       read (19, '(F16.8)') counter

     do i_b_LC = 1,b_num_LC
         do i_b_long = 1,b_num_long
             do i_y = 1,y_num
                do i_premium = 1,premium_num
                READ(10, '(G15.6, X, G15.6, X, G15.6)') v_matrix(i_b_LC, i_b_long, i_y, i_premium), &
                    v0_matrix(i_b_LC, i_b_long, i_y, i_premium), v1_matrix(i_b_LC, i_b_long, i_y, i_premium)
                if (v0_matrix(i_b_LC, i_b_long, i_y, i_premium)< v1_matrix(i_b_LC, i_b_long, i_y, i_premium)) then
                        default_decision(i_b_LC, i_b_long, i_y, i_premium) =2
                    else
                        default_decision(i_b_LC, i_b_long, i_y, i_premium) =1
                    end if
                READ(13, '(F15.11, X, F15.11)') b_next_matrix(i_b_LC, i_b_long, i_y, i_premium,1), b_next_matrix(i_b_LC, i_b_long, i_y, i_premium,2)
                READ(16, '(F15.8, X, F15.8, X, F15.8, X, F15.8, X, F15.8)')  q_paid_matrix(i_b_LC, i_b_long, i_y, i_premium),  q_nodef_matrix(i_b_LC, i_b_long, i_y, i_premium), &
                                                                   q_nodef_LC_matrix(i_b_LC, i_b_long, i_y, i_premium) 
                READ(17, '(F15.8, X, F15.8)') q_def_matrix(i_b_LC, i_b_long, i_y, i_premium),&
                     q_LC_def_matrix(i_b_LC, i_b_long, i_y, i_premium)
                
                READ(18, '(G15.6, X, G15.6, X, G15.6)') b0_next_LC_matrix(i_b_LC, i_b_long, i_y, i_premium), b0_next_long_matrix(i_b_LC, i_b_long, i_y, i_premium), &
                                        b1_next_LC_matrix(i_b_LC, i_b_long, i_y, i_premium)
                READ(20, '(G15.6)') pi_matrix(i_b_LC, i_b_long, i_y, i_premium)
                end do
              end do
         end do
     end do
     CLOSE(10)
     CLOSE(13)
     CLOSE(16)
     CLOSE(17)
     close(18)
     close(19)
     close(20)
end if

call iterate
call simulate
!call compute_initial_guess

call cpu_time(end_time)
WRITE(nout, '(A7, X, A7, X, A7)') 'Hours ', 'Minutes', 'Seconds'
WRITE(nout, '(I7, X, I7, X, I7)') INT((end_time - start_time) / 3600d+0), &
             INT((end_time-start_time)/60d+0 - INT((end_time - start_time) / 3600d+0)*60d+0),&
INT(end_time-start_time - INT((end_time - start_time) / 3600d+0)*3600d+0 - &
INT((end_time-start_time)/60d+0 - INT((end_time - start_time) / 3600d+0)*60d+0)*60d+0)
end program

SUBROUTINE POWELL(P,XI,N,NP,FTOL,ITER,FRET)
!-----------------------------------------------------------
! Minimization of a function  FUNC of N variables  (FUNC is
! not an argument, it is a fixed function name). Input con-
! sists of an initial starting point P  that is a vector of
! length N; an initial matrix XI  whose  logical dimensions
! are N by N, physical dimensions NP by NP, and whose columns
! contain the initial set of directions (usually the N unit
! vectors); and FTOL, the fractional tolerance in the func-
! tion value such that failure to decrease by more than this
! amount on one iteration signals doneness. On output, P is
! set to the best point found, XI is the then-current direc-
! tion set,  FRET is the returned function value at P,  and
! ITER is the number of iterations taken. The routine LINMIN
! is used.
!------------------------------------------------------------
  IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
  PARAMETER(NMAX=20,ITMAX=100)
  DIMENSION P(NP),XI(NP,NP),PT(NMAX),PTT(NMAX),XIT(NMAX)
  FRET=objective_function(P)
  DO J=1,N
    PT(J)=P(J)       !Save initial pont
  END DO
  ITER=0
1 ITER=ITER+1
  FP=FRET
  IBIG=0
  DEL=0.D0           !Will be the biggest function decrease.
  DO I=1,N           !In each iteration, loop over all directions in the set.
    DO J=1,N         !Copy the direction
      XIT(J)=XI(J,I)
    END DO
    FPTT=FRET
    CALL LINMIN(P,XIT,N,FRET)  !Minimize along it.
    IF (DABS(FPTT-FRET).GT.DEL) THEN
      DEL=DABS(FPTT-FRET)
      IBIG=I
    END IF
  END DO
  IF (2.D0*DABS(FP-FRET).LE.FTOL*(DABS(FP)+DABS(FRET))) RETURN !Termination criterion
  IF (ITER.EQ.ITMAX) Then
    Pause ' Powell exceeding maximum iterations.'
    return
  END IF		 
  DO J=1,N
    PTT(J)=2.D0*P(J)-PT(J)  !Construct the extrapolated point and the average
    XIT(J)=P(J)-PT(J)       !direction moved. Save the old starting point.
    PT(J)=P(J)
  END DO
  FPTT=objective_function(PTT)            !Function value at extrapolated point.
  IF (FPTT.GE.FP) GO TO 1   !One reason not to use new direction. 
  T=2.D0*(FP-2.D0*FRET+FPTT)*(FP-FRET-DEL)**2-DEL*(FP-FPTT)**2
  IF (T.GE.0.D0) GO TO 1    !Other reason not to use new direction.
  CALL LINMIN(P,XIT,N,FRET) !Move to the minimum of the new direction.
  DO J=1,N                  !and save the new direction
    XI(J,IBIG)=XIT(J)
  END DO
  GO TO 1
END

SUBROUTINE LINMIN(P,XI,N,FRET)
!----------------------------------------------------------
! Given an N dimensional point P and a N dimensional direc-
! tion XI, moves and resets P to where the function FUNC(P)
! takes on a minimum along the direction XI from P, and 
! replaces XI by the actual vector displacement that P was
! moved. Also returns as FRET the value of FUNC at the
! returned location P. This is actually all accomplished by
! calling the routines MNBRAK and BRENT.
!----------------------------------------------------------
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  PARAMETER(NMAX=50,TOL=1.D-4)
  DIMENSION P(N),XI(N)
  COMMON /F1COM/ PCOM(NMAX),XICOM(NMAX),NCOM
  NCOM=N
  DO J=1,N
    PCOM(J)=P(J)
    XICOM(J)=XI(J)
  END DO
  AX=0.D0
  XX=1.D0
  BX=2.D0
  CALL MNBRAK(AX,XX,BX,FA,FX,FB)
  FRET=BRENT(AX,XX,BX,TOL,XMIN)
  DO J=1,N
    XI(J)=XMIN*XI(J)
    P(J)=P(J)+XI(J)
  END DO
  RETURN
END

DOUBLE PRECISION FUNCTION F1DIM(X)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  PARAMETER(NMAX=50)
  COMMON /F1COM/ PCOM(NMAX),XICOM(NMAX),NCOM
  DIMENSION XT(NMAX)
  DO J=1, NCOM
    XT(J)=PCOM(J)+X*XICOM(J)
  END DO
  F1DIM = objective_function(XT)
  RETURN
END

SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC)
!----------------------------------------------------------------------
!Given a Function F1DIM(X), and given distinct initial points AX and
!BX, this routine searches in the downhill direction (defined by the
!F1DIMtion as evaluated at the initial points) and returns new points
!AX, BX, CX which bracket a minimum of the Function. Also returned
!are the Function values at the three points, FA, FB and FC.
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
PARAMETER(GOLD=1.618034,GLIMIT=100.,TINY=1.D-20)
!The first parameter is the default ratio by which successive intervals
!are magnified; the second is the maximum magnification allowed for
!a parabolic-fit step.
!----------------------------------------------------------------------
FA=F1DIM(AX)
FB=F1DIM(BX)
IF(FB.GT.FA) THEN
  DUM=AX
  AX=BX
  BX=DUM
  DUM=FB
  FB=FA
  FA=DUM
ENDIF
CX=BX+GOLD*(BX-AX)
FC=F1DIM(CX)
1 IF(FB.GE.FC) THEN
  R=(BX-AX)*(FB-FC)
  Q=(BX-CX)*(FB-FA)
  U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
  ULIM=BX+GLIMIT*(CX-BX)
  IF((BX-U)*(U-CX).GT.0) THEN
    FU=F1DIM(U)
    IF(FU.LT.FC) THEN
      AX=BX
      FA=FB
      BX=U
      FB=FU
      GOTO 1
    ELSE IF(FU.GT.FB) THEN
      CX=U
      FC=FU
      GOTO 1
    ENDIF
    U=CX+GOLD*(CX-BX)
    FU=F1DIM(U)
  ELSE IF((CX-U)*(U-ULIM).GT.0) THEN
    FU=F1DIM(U)
    IF(FU.LT.FC) THEN
      BX=CX
      CX=U
      U=CX+GOLD*(CX-BX)
      FB=FC
      FC=FU
      FU=F1DIM(U)
    ENDIF
  ELSE IF((U-ULIM)*(ULIM-CX).GE.0) THEN
    U=ULIM
    FU=F1DIM(U)
  ELSE
    U=CX+GOLD*(CX-BX)
    FU=F1DIM(U)
  ENDIF
  AX=BX
  BX=CX
  CX=U
  FA=FB
  FB=FC
  FC=FU
  GOTO 1
ENDIF
RETURN
END

DOUBLE PRECISION FUNCTION BRENT(AX,BX,CX,TOL,XMIN)
!-------------------------------------------------------------------
!Given a function F1DIM, and a bracketing triplet of abscissas
!AX,BX,CX (such that BX is between AX and CX, and F(BX) is less 
!than both F(AX) and F(CX)), this routine isolates the minimum 
!to a fractional precision of about TOL using Brent's method.
!The abscissa of the minimum is returned in XMIN, and the minimum
!function value is returned as BRENT, the returned function value.
!-------------------------------------------------------------------
PARAMETER(ITMAX=100,CGOLD=.3819660,ZEPS=1.D-10)
!Maximum allowed number of iterations; golden ration; and a small
!number which protects against trying to achieve fractional accuracy
!for a minimum that happens to be exactly zero.
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
A=MIN(AX,CX)
B=MAX(AX,CX)
V=BX
W=V
X=V
E=0.
FX=F1DIM(X)
FV=FX
FW=FX
DO 11 ITER=1,ITMAX	                                !main loop
  XM=0.5*(A+B)
  TOL1=TOL*ABS(X)+ZEPS
  TOL2=2.*TOL1
  IF (ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3  !Test for done here
  IF (ABS(E).GT.TOL1) THEN     !Construct a trial parabolic fit
    R=(X-W)*(FX-FV)
    Q=(X-V)*(FX-FW)
    P=(X-V)*Q-(X-W)*R
    Q=.2*(Q-R)
    IF (Q.GT.0)  P=-P
    Q=ABS(Q)
    ETEMP=E
    E=D
    IF (ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR.  &
	P.GE.Q*(B-X))  GOTO 1
!   The above conditions determine the acceptability of the 
!   parabolic fit. Here it is o.k.:
    D=P/Q
    U=X+D
    IF (U-A.LT.TOL2.OR.B-U.LT.TOL2)  D=SIGN(TOL1,XM-X)
    GOTO 2
  ENDIF
1 IF (X.GE.XM) THEN
    E=A-X
  ELSE
    E=B-X
  ENDIF
  D=CGOLD*E
2 IF (ABS(D).GE.TOL1) THEN
    U=X+D
  ELSE
    U=X+SIGN(TOL1,D)
  ENDIF
  FU=F1DIM(U)  !This the one function evaluation per iteration
  IF (FU.LE.FX) THEN
    IF (U.GE.X) THEN
      A=X
    ELSE
      B=X
    ENDIF
    V=W
    FV=FW
    W=X
    FW=FX
    X=U
    FX=FU
  ELSE
    IF (U.LT.X) THEN
      A=U
    ELSE
      B=U
    ENDIF
    IF (FU.LE.FW.OR.W.EQ.X) THEN
      V=W
      FV=FW
      W=U
      FW=FU
    ELSE IF (FU.LE.FV.OR.V.EQ.X.OR.V.EQ.W) THEN
      V=U
      FV=FU
    ENDIF
  ENDIF
11 CONTINUE
!Pause ' Brent exceed maximum iterations.'
3 XMIN=X   !exit section
  BRENT=FX
  RETURN
  END


!SUBROUTINE POWELL(FUNC,P,XI,N,NP,FTOL,ITER,FRET)
!!-----------------------------------------------------------
!! Minimization of a function  FUNC of N variables  (FUNC is
!! not an argument, it is a fixed function name). Input con-
!! sists of an initial starting point P  that is a vector of
!! length N; an initial matrix XI  whose  logical dimensions
!! are N by N, physical dimensions NP by NP, and whose columns26
!! contain the initial set of directions (usually the N unit
!! vectors); and FTOL, the fractional tolerance in the func-
!! tion value such that failure to decrease by more than this
!! amount on one iteration signals doneness. On output, P is
!! set to the best point found, XI is the then-current direc-
!! tion set,  FRET is the returned function value at P,  and
!! ITER is the number of iterations taken. The routine LINMIN
!! is used.
!!------------------------------------------------------------
!  !IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
!  INTEGER :: N, NP, NMAX, ITMAX, I, IBIG, J, ITER
!  PARAMETER(NMAX=2,ITMAX=100)
!  DOUBLE PRECISION :: P(NP),XI(NP,NP),PT(NMAX),PTT(NMAX),XIT(NMAX)
!  DOUBLE PRECISION :: FRET, FP, DEL, FPTT, FTOL, T
!  !EXTERNAL objective_function_nobind
!  
!  interface
!  DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
!    DOUBLE PRECISION, intent(in) :: ARG(2)
!  end function
!  end interface
!  
!  FRET=  FUNC(P) !objective_function_nobind(P) !
!  !write(nout, *) FRET, objective_function_nobind(P)
!  !pause
!  
!  DO J=1,N
!    PT(J)=P(J)       !Save initial pont
!  END DO
!  ITER=0
!1 ITER=ITER+1
!  FP=FRET
!  IBIG=0
!  DEL=0.0D+0           !Will be the biggest function decrease.
!  DO I=1,N           !In each iteration, loop over all directions in the set.
!    DO J=1,N         !Copy the direction
!      XIT(J)=XI(J,I)
!    END DO
!    FPTT=FRET
!    CALL LINMIN(FUNC,P,XIT,N,FRET)  !Minimize along it.
!    IF (DABS(FPTT-FRET).GT.DEL) THEN
!      DEL=DABS(FPTT-FRET)
!      IBIG=I
!    END IF
!  END DO
!  IF (2.D0*DABS(FP-FRET).LE.FTOL*(DABS(FP)+DABS(FRET))) RETURN !Termination criterion
!  IF (ITER.EQ.ITMAX) Then
!    Pause ' Powell exceeding maximum iterations.'
!    return
!  END IF		 
!  DO J=1,N
!    PTT(J)=2.D0*P(J)-PT(J)  !Construct the extrapolated point and the average
!    XIT(J)=P(J)-PT(J)       !direction moved. Save the old starting point.
!    PT(J)=P(J)
!  END DO
!  FPTT= FUNC(PTT) ! objective_function_nobind(PTT) !          !Function value at extrapolated point.
!!  write(nout, '(F12.8, X, F12.8, X, F15.6, X, F15.6)') PTT(1), PTT(2), FPTT, objective_function_nobind(PTT)
!  IF (FPTT.GE.FP) GO TO 1   !One reason not to use new direction. 
!  T=2.D0*(FP-2.D0*FRET+FPTT)*(FP-FRET-DEL)**2D+0-DEL*(FP-FPTT)**2D+0
!  IF (T.GE.0.D0) GO TO 1    !Other reason not to use new direction.
!  CALL LINMIN(FUNC,P,XIT,N,FRET) !Move to the minimum of the new direction.
!  DO J=1,N                  !and save the new direction
!    XI(J,IBIG)=XIT(J)
!  END DO
!  GO TO 1
!END
!
!SUBROUTINE LINMIN(FUNC,P,XI,N,FRET)
!!----------------------------------------------------------
!! Given an N dimensional point P and a N dimensional direc-
!! tion XI, moves and resets P to where the function FUNC(P)
!! takes on a minimum along the direction XI from P, and 
!! replaces XI by the actual vector displacement that P was
!! moved. Also returns as FRET the value of FUNC at the
!! returned location P. This is actually all accomplished by
!! calling the routines MNBRAK and BRENT.
!!----------------------------------------------------------
!  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!  DOUBLE PRECISION :: TOL
!  INTEGER :: N, NMAX, NCOM, J
!  PARAMETER(NMAX=50,TOL=1.D-4)
!  DOUBLE PRECISION :: P(N),XI(N) !, PCOM(N), XICOM(N)
!  DOUBLE PRECISION :: AX, XX, BX, FA, FX,FB, FRET, BRENT
!  COMMON /F1COM/ PCOM(NMAX),XICOM(NMAX),NCOM
!
!  interface
!  DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
!    DOUBLE PRECISION, intent(in) :: ARG(2)
!  end function
!  end interface
!  NCOM=N
!  DO J=1,N
!    PCOM(J)=P(J)
!    XICOM(J)=XI(J)
!  END DO
!  AX=0.D0
!  XX=1.D0
!  BX=2.D0
!  CALL MNBRAK(FUNC,AX,XX,BX,FA,FX,FB)
!  FRET=BRENT(FUNC,AX,XX,BX,TOL,XMIN)
!  DO J=1,N
!    XI(J)=XMIN*XI(J)
!    P(J)=P(J)+XI(J)
!  END DO
!  RETURN
!END
!
!DOUBLE PRECISION FUNCTION F1DIM(FUNC,X)
!  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!  INTEGER :: NMAX, J, NCOM
!  PARAMETER(NMAX=50)
!  COMMON /F1COM/ PCOM(NMAX),XICOM(NMAX),NCOM
!  DOUBLE PRECISION :: XT(NMAX), X
!  !EXTERNAL objective_function_nobind
! ! EXTERNAL FUNC
!
!  interface
!  DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
!    DOUBLE PRECISION, intent(in) :: ARG(2)
!  end function
!  end interface
!  
!  DO J=1, NCOM
!    XT(J)=PCOM(J)+X*XICOM(J)
!  END DO
!  F1DIM = FUNC(XT) !objective_function_nobind(XT) !
!  RETURN
!END
!
!SUBROUTINE MNBRAK(FUNC,AX,BX,CX,FA,FB,FC)
!!----------------------------------------------------------------------
!!Given a Function F1DIM(X), and given distinct initial points AX and
!!BX, this routine searches in the downhill direction (defined by the
!!F1DIMtion as evaluated at the initial points) and returns new points
!!AX, BX, CX which bracket a minimum of the Function. Also returned
!!are the Function values at the three points, FA, FB and FC.
!IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!DOUBLE PRECISION :: GOLD, GLIMIT, TINY, FA, FB,FC, AX,BX,CX, DUM,R,Q,U,ULIM,FU, F1DIM
!PARAMETER(GOLD=1.618034,GLIMIT=100.,TINY=1.D-20)
!  interface
!  DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
!    DOUBLE PRECISION, intent(in) :: ARG(2)
!  end function
!  end interface
!
!!The first parameter is the default ratio by which successive intervals
!!are magnified; the second is the maximum magnification allowed for
!!a parabolic-fit step.
!!----------------------------------------------------------------------
!FA=F1DIM(FUNC,AX)
!FB=F1DIM(FUNC,BX)
!IF(FB.GT.FA) THEN
!  DUM=AX
!  AX=BX
!  BX=DUM
!  DUM=FB
!  FB=FA
!  FA=DUM
!ENDIF
!CX=BX+GOLD*(BX-AX)
!FC=F1DIM(FUNC,CX)
!1 IF(FB.GE.FC) THEN
!  R=(BX-AX)*(FB-FC)
!  Q=(BX-CX)*(FB-FA)
!  U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
!  ULIM=BX+GLIMIT*(CX-BX)
!  IF((BX-U)*(U-CX).GT.0) THEN
!    FU=F1DIM(FUNC,U)
!    IF(FU.LT.FC) THEN
!      AX=BX
!      FA=FB
!      BX=U
!      FB=FU
!      GOTO 1
!    ELSE IF(FU.GT.FB) THEN
!      CX=U
!      FC=FU
!      GOTO 1
!    ENDIF
!    U=CX+GOLD*(CX-BX)
!    FU=F1DIM(FUNC,U)
!  ELSE IF((CX-U)*(U-ULIM).GT.0) THEN
!    FU=F1DIM(FUNC,U)
!    IF(FU.LT.FC) THEN
!      BX=CX
!      CX=U
!      U=CX+GOLD*(CX-BX)
!      FB=FC
!      FC=FU
!      FU=F1DIM(FUNC,U)
!    ENDIF
!  ELSE IF((U-ULIM)*(ULIM-CX).GE.0) THEN
!    U=ULIM
!    FU=F1DIM(FUNC,U)
!  ELSE
!    U=CX+GOLD*(CX-BX)
!    FU=F1DIM(FUNC,U)
!  ENDIF
!  AX=BX
!  BX=CX
!  CX=U
!  FA=FB
!  FB=FC
!  FC=FU
!  GOTO 1
!ENDIF
!RETURN
!END
!
!DOUBLE PRECISION FUNCTION BRENT(FUNC,AX,BX,CX,TOL,XMIN)
!!-------------------------------------------------------------------
!!Given a function F1DIM, and a bracketing triplet of abscissas
!!AX,BX,CX (such that BX is between AX and CX, and F(BX) is less 
!!than both F(AX) and F(CX)), this routine isolates the minimum 
!!to a fractional precision of about TOL using Brent's method.
!!The abscissa of the minimum is returned in XMIN, and the minimum
!!function value is returned as BRENT, the returned function value.
!!-------------------------------------------------------------------
!IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!INTEGER :: ITMAX,ITER
!DOUBLE PRECISION :: CGOLD, ZEPS
!PARAMETER(ITMAX=100,CGOLD=.3819660,ZEPS=1.D-10)
!!Maximum allowed number of iterations; golden ration; and a small
!!number which protects against trying to achieve fractional accuracy
!!for a minimum that happens to be exactly zero.
!DOUBLE PRECISION :: AX,BX, CX,TOL,XMIN,A,B,V,W,X,E,FX,FV,FU,FW,XM,TOL1,TOL2,R,P,Q,ETEMP,D,U
!interface
!  DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
!    DOUBLE PRECISION, intent(in) :: ARG(2)
!  end function
!  end interface
!
!
!A=MIN(AX,CX)
!B=MAX(AX,CX)
!V=BX
!W=V
!X=V
!E=0.
!FX=F1DIM(FUNC,X)
!FV=FX
!FW=FX
!DO 11 ITER=1,ITMAX	                                !main loop
!  XM=0.5*(A+B)
!  TOL1=TOL*ABS(X)+ZEPS
!  TOL2=2.*TOL1
!  IF (ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3  !Test for done here
!  IF (ABS(E).GT.TOL1) THEN     !Construct a trial parabolic fit
!    R=(X-W)*(FX-FV)
!    Q=(X-V)*(FX-FW)
!    P=(X-V)*Q-(X-W)*R
!    Q=.2*(Q-R)
!    IF (Q.GT.0)  P=-P
!    Q=ABS(Q)
!    ETEMP=E
!    E=D
!    IF (ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR.  &
!	P.GE.Q*(B-X))  GOTO 1
!!   The above conditions determine the acceptability of the 
!!   parabolic fit. Here it is o.k.:
!    D=P/Q
!    U=X+D
!    IF (U-A.LT.TOL2.OR.B-U.LT.TOL2)  D=SIGN(TOL1,XM-X)
!    GOTO 2
!  ENDIF
!1 IF (X.GE.XM) THEN
!    E=A-X
!  ELSE
!    E=B-X
!  ENDIF
!  D=CGOLD*E
!2 IF (ABS(D).GE.TOL1) THEN
!    U=X+D
!  ELSE
!    U=X+SIGN(TOL1,D)
!  ENDIF
!  FU=F1DIM(FUNC,U)  !This the one function evaluation per iteration
!  IF (FU.LE.FX) THEN
!    IF (U.GE.X) THEN
!      A=X
!    ELSE
!      B=X
!    ENDIF
!    V=W
!    FV=FW
!    W=X
!    FW=FX
!    X=U
!    FX=FU
!  ELSE
!    IF (U.LT.X) THEN
!      A=U
!    ELSE
!      B=U
!    ENDIF
!    IF (FU.LE.FW.OR.W.EQ.X) THEN
!      V=W
!      FV=FW
!      W=U
!      FW=FU
!    ELSE IF (FU.LE.FV.OR.V.EQ.X.OR.V.EQ.W) THEN
!      V=U
!      FV=FU
!    ENDIF
!  ENDIF
!11 CONTINUE
!!Pause ' Brent exceed maximum iterations.'
!3 XMIN=X   !exit section
!  BRENT=FX
!  RETURN
!  END
!!end of file tpowell.f90
!
!subroutine BISEC(func, x1, x2, ERRABS, ERRREL)
!implicit none
!DOUBLE PRECISION, INTENT(IN) ::  ERRABS, ERRREL
!DOUBLE PRECISION, INTENT(IN OUT) :: x2, x1
!DOUBLE PRECISION :: diff
!INTEGER :: nout
!
!interface
!        function func(x)
!        implicit none
!        DOUBLE PRECISION, INTENT(IN) :: x
!        DOUBLE PRECISION :: func
!        end function func
!END interface
!DOUBLE PRECISION :: fl, fh, f_mean, xl, xh, xmean, rel_error
!
!xl=x1
!xh=x2
!fl=func(xl)
!fh=func(xh)
!diff = fh-fl
!rel_error=1
!
!do WHILE( MIN((ABS(fh)-ERRABS), rel_error-ERRREL) > 0)
!
!xmean =  0.5d+0 * (xl+xh)
!!xl - fl / ((fh - fl)/(xh-xl) )  !Use linear approx to find the zero !(xl+xh)/2
!f_mean = func(xmean)
!
!if (fl*fh>0) then
!   !WRITE(nout,*) 'Wrong range'       !Need change of signs to compute the root
!   !call BREAK
!else
!   if (diff < 0) then
!      if (f_mean > 0) then
!         xl=xmean
!         fl=f_mean
!      else
!         xh=xmean
!         fh=f_mean
!      end if
!   else
!      if (f_mean > 0) then
!         xh=xmean
!         fh=f_mean
!      else
!         xl=xmean
!         fl=f_mean
!      end if
!    end if
!
!end if
!rel_error=ABS(xh-xl)
!
!!WRITE(6,*) rel_error, fh
!end do
!x1 = xl
!x2 = xh
!
!
!end subroutine    
    


    