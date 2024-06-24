
subroutine dthetaF1ast_f(dthetaF1ast_func, M, a, r, theta, l, s0, kappa, gamma)

    
    double precision, intent(in) :: M, a, r, theta, l, s0, kappa, gamma
    double precision, intent(out) :: dthetaF1ast_func


    dthetaF1ast_func = (s0*((a*M*r*(3*a**2 + r*(-4*M + 3*r))*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*&
&        Sin(theta)**2)/(2.*(a**2 - 2*M*r + r**2)*(r**2 + a**2*Cos(theta)**2)**3) - &
&     l*((M*r*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&           (-27*a**4 - 4*a**2*r**2 + 16*r**4 - &
&             4*a**2*(6*a**2 + 7*r**2)*Cos(2*theta) + 3*a**4*Cos(4*theta)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (24*a**2*M**2*r**2*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*Sin(theta)**2)/&
&         ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&           (a**2 + 2*r**2 + a**2*Cos(2*theta))**3)) - &
&     (8*M*r*(r**2 + a**2*Cos(theta)**2)*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*&
&        (2*a**4 + r**4 + a**2*r*(-2*M + 3*r) - &
&          a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**2*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) - &
&     (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&        ((12*a*M*r*(a**2 + r**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&             (3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*Sin(theta)**2)/&
&           ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&          (4*a*M**2*r**2*(-9*a**4 + 2*a**2*r**2 + 4*r**4 - &
&               2*a**2*(3*a**2 + 5*r**2)*Cos(2*theta) + 3*a**4*Cos(4*theta))*Sin(theta)**2)/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3)))/&
&      ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&        Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))&
&     )*((-2*a**2*(1/(Tan(theta))))/(a**2 + 2*r**2 + a**2*Cos(2*theta))**2 - &
&     (2*(r**2 + a**2*Cos(theta)**2)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)/&
&      (a**2 + 2*r**2 + a**2*Cos(2*theta))**2 + &
&     (4*a**2*(r**2 + a**2*Cos(theta)**2)*(1/(Sin(theta)))**2*Sin(2*theta))/&
&      (a**2 + 2*r**2 + a**2*Cos(2*theta))**3))/&
& (gamma*kappa*Sqrt(((r**2 + a**2*Cos(theta)**2)*(1/(Sin(theta)))**2)/&
&     (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)) + &
&(2*s0*Sqrt(((r**2 + a**2*Cos(theta)**2)*(1/(Sin(theta)))**2)/&
&     (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)*&
&   ((a*M*r*(3*a**2 + r*(-4*M + 3*r))*Cos(theta)*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*&
&        Sin(theta))/((a**2 - 2*M*r + r**2)*(r**2 + a**2*Cos(theta)**2)**3) + &
&     (3*a**3*M*r*(3*a**2 + r*(-4*M + 3*r))*Cos(theta)*&
&        (3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*Sin(theta)**3)/&
&      ((a**2 - 2*M*r + r**2)*(r**2 + a**2*Cos(theta)**2)**4) - &
&     (16*M*r*Cos(theta)*(r**2 + a**2*Cos(theta)**2)*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*&
&        (2*a**4 + r**4 + a**2*r*(-2*M + 3*r) - &
&          a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) + &
&     (16*a**2*M*r*Cos(theta)*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*&
&        (2*a**4 + r**4 + a**2*r*(-2*M + 3*r) - &
&          a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**3*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) - &
&     (8*M*r*(r**2 + a**2*Cos(theta)**2)*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*&
&        (2*a**4 + r**4 + a**2*r*(-2*M + 3*r) - &
&          a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**2*&
&        ((-4*a**2*l*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2)**2 + &
&          (4*a*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&          (4*a**3*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) - &
&     (l*((12*a*M*r*(a**2 + r**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&             (3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*Sin(theta)**2)/&
&           ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&          (4*a*M**2*r**2*(-9*a**4 + 2*a**2*r**2 + 4*r**4 - &
&               2*a**2*(3*a**2 + 5*r**2)*Cos(2*theta) + 3*a**4*Cos(4*theta))*Sin(theta)**2)/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3))*&
&        ((-4*a**2*l*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2)**2 + &
&          (4*a*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&          (4*a**3*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2))/&
&      ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&        Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))&
&       + (8*M*r*(r**2 + a**2*Cos(theta)**2)*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*&
&        (2*a**4 + r**4 + a**2*r*(-2*M + 3*r) - &
&          a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**2*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&        ((-4*a*l*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) - &
&          (4*a**3*l*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2 + &
&          2*Cos(theta)*Sin(theta)*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)) + &
&          Sin(theta)**2*((4*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&             (4*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2)))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&           Sin(theta)**2*(a**2 + r**2 + &
&              (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))**2) + &
&     (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&        ((12*a*M*r*(a**2 + r**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&             (3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*Sin(theta)**2)/&
&           ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&          (4*a*M**2*r**2*(-9*a**4 + 2*a**2*r**2 + 4*r**4 - &
&               2*a**2*(3*a**2 + 5*r**2)*Cos(2*theta) + 3*a**4*Cos(4*theta))*Sin(theta)**2)/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3))*&
&        ((-4*a*l*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) - &
&          (4*a**3*l*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2 + &
&          2*Cos(theta)*Sin(theta)*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)) + &
&          Sin(theta)**2*((4*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&             (4*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2)))/&
&      ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&         Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))&
&         )**2 - (3*a**3*M*r*(3*a**2 + r*(-4*M + 3*r))*Sin(theta)**2*Sin(2*theta))/&
&      ((a**2 - 2*M*r + r**2)*(r**2 + a**2*Cos(theta)**2)**3) - &
&     (16*a**2*M*r*(a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&        (3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*Sin(theta)**2*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*Sin(2*theta))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) + &
&     (48*a**2*M*r*(r**2 + a**2*Cos(theta)**2)*&
&        (2*a**4 + r**4 + a**2*r*(-2*M + 3*r) - &
&          a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**2*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*Sin(2*theta))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) - &
&     (64*a**2*M*r*(r**2 + a**2*Cos(theta)**2)*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*&
&        (2*a**4 + r**4 + a**2*r*(-2*M + 3*r) - &
&          a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**2*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*Sin(2*theta))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**5*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) - &
&     (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&        ((24*a*M*r*(a**2 + r**2)*Cos(theta)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&             (3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*Sin(theta))/&
&           ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&          (8*a*M**2*r**2*Cos(theta)*&
&             (-9*a**4 + 2*a**2*r**2 + 4*r**4 - &
&               2*a**2*(3*a**2 + 5*r**2)*Cos(2*theta) + 3*a**4*Cos(4*theta))*Sin(theta))/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&          (48*a**3*M**2*r**2*(a**2 + r**2)*Cos(theta)*&
&             (3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*Sin(theta)**3)/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**2*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&          (8*a**3*M**2*r**2*Cos(theta)*&
&             (-9*a**4 + 2*a**2*r**2 + 4*r**4 - &
&               2*a**2*(3*a**2 + 5*r**2)*Cos(2*theta) + 3*a**4*Cos(4*theta))*Sin(theta)**3)/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**2*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&          (72*a**3*M*r*(a**2 + r**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&             Sin(theta)**2*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&          (72*a**3*M*r*(a**2 + r**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&             (3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*Sin(theta)**2*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) + &
&          (24*a**3*M**2*r**2*&
&             (-9*a**4 + 2*a**2*r**2 + 4*r**4 - &
&               2*a**2*(3*a**2 + 5*r**2)*Cos(2*theta) + 3*a**4*Cos(4*theta))*Sin(theta)**2*&
&             Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**4) + &
&          (4*a*M**2*r**2*Sin(theta)**2*&
&             (4*a**2*(3*a**2 + 5*r**2)*Sin(2*theta) - 12*a**4*Sin(4*theta)))/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3)))/&
&      ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&        Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))&
&       - l*((48*a**2*M**2*r**2*Cos(theta)*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*&
&           Sin(theta))/&
&         ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&           (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (4*a**2*M**2*r**2*Cos(theta)*&
&           (-27*a**4 - 4*a**2*r**2 + 16*r**4 - &
&             4*a**2*(6*a**2 + 7*r**2)*Cos(2*theta) + 3*a**4*Cos(4*theta))*Sin(theta))/&
&         ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**2*&
&           (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (48*a**4*M**2*r**2*Cos(theta)*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*Sin(theta)**3)/&
&         ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**2*&
&           (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (6*a**2*M*r*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&           (-27*a**4 - 4*a**2*r**2 + 16*r**4 - &
&             4*a**2*(6*a**2 + 7*r**2)*Cos(2*theta) + 3*a**4*Cos(4*theta))*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) - &
&        (144*a**4*M**2*r**2*Sin(theta)**2*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&           (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (144*a**4*M**2*r**2*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*Sin(theta)**2*&
&           Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&           (a**2 + 2*r**2 + a**2*Cos(2*theta))**4) + &
&        (M*r*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&           (8*a**2*(6*a**2 + 7*r**2)*Sin(2*theta) - 12*a**4*Sin(4*theta)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3))))/(gamma*kappa)

end subroutine dthetaF1ast_f
