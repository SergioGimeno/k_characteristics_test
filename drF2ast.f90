
subroutine drF2ast_f(drF2ast_func, M, a, r, theta, l, s0, kappa, gamma)

    
    double precision, intent(in) :: M, a, r, theta, l, s0, kappa, gamma
    double precision, intent(out) :: drF2ast_func


    drF2ast_func = (2*s0*Sqrt(((r**2 + a**2*Cos(theta)**2)*(1/(Sin(theta)))**2)/&
&     (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)*&
&   ((-24*a**2*M*(a**2 + r**2)*(a**2 + r*(-2*M + r))*Cos(theta)*&
&        (r**2 + a**2*Cos(theta)**2)*(a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**3*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&        ((4*a*l*M*r**2*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)**2 - &
&          (2*a*l*M*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(2*r - (4*a**2*M*r**2*Sin(theta)**2)/&
&              (r**2 + a**2*Cos(theta)**2)**2 + &
&             (2*a**2*M*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&           Sin(theta)**2*(a**2 + r**2 + &
&              (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))**2) + &
&     (24*a**2*M*(a**2 + r**2)*(a**2 + r*(-2*M + r))*Cos(theta)*&
&        (r**2 + a**2*Cos(theta)**2)*(a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**3*&
&        (-(l*((-4*M*r**2)/(r**2 + a**2*Cos(theta)**2)**2 + &
&               (2*M)/(r**2 + a**2*Cos(theta)**2))) - &
&          (4*a*M*r**2*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)**2 + &
&          (2*a*M*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) - &
&     (384*a**2*M*r*(a**2 + r**2)*(a**2 + r*(-2*M + r))*Cos(theta)*&
&        (r**2 + a**2*Cos(theta)**2)*(a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**3*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**5*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) - &
&     (288*a**2*M*r*(a**2 + r**2)*(a**2 + r*(-2*M + r))*Cos(theta)*&
&        (r**2 + a**2*Cos(theta)**2)*Sin(theta)**3*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) + &
&     (48*a**2*M*r*(a**2 + r**2)*(a**2 + r*(-2*M + r))*Cos(theta)*&
&        (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**3*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) + &
&     (24*a**2*M*(-2*M + 2*r)*(a**2 + r**2)*Cos(theta)*(r**2 + a**2*Cos(theta)**2)*&
&        (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**3*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) - &
&     (24*a**2*M*(-2*M + 2*r)*(a**2 + r**2)*(a**2 + r*(-2*M + r))*Cos(theta)*&
&        (r**2 + a**2*Cos(theta)**2)*(a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**3*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((a**2 - 2*M*r + r**2)**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) + &
&     (48*a**2*M*r*(a**2 + r*(-2*M + r))*Cos(theta)*(r**2 + a**2*Cos(theta)**2)*&
&        (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**3*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) + &
&     (2*a*M*(a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&        (44*a**2*r + 48*r**3 - 28*a**2*r*Cos(2*theta))*Sin(2*theta))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) - &
&     (32*a*M*r*(a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&        (-3*a**4 + 22*a**2*r**2 + 12*r**4 - 2*a**2*(a**2 + 7*r**2)*Cos(2*theta) + &
&          a**4*Cos(4*theta))*Sin(2*theta))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**5) + &
&     (4*a*M*r*(a**2 + r*(-2*M + r))*&
&        (-3*a**4 + 22*a**2*r**2 + 12*r**4 - 2*a**2*(a**2 + 7*r**2)*Cos(2*theta) + &
&          a**4*Cos(4*theta))*Sin(2*theta))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) + &
&     (2*a*M*(-2*M + 2*r)*(r**2 + a**2*Cos(theta)**2)*&
&        (-3*a**4 + 22*a**2*r**2 + 12*r**4 - 2*a**2*(a**2 + 7*r**2)*Cos(2*theta) + &
&          a**4*Cos(4*theta))*Sin(2*theta))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) - &
&     (2*a*M*(-2*M + 2*r)*(a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&        (-3*a**4 + 22*a**2*r**2 + 12*r**4 - 2*a**2*(a**2 + 7*r**2)*Cos(2*theta) + &
&          a**4*Cos(4*theta))*Sin(2*theta))/&
&      ((a**2 - 2*M*r + r**2)**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) - &
&     l*((-96*a**2*M**2*r**2*Cos(theta)*(-2*a**2 + (2*M - r)*r + a**2*Cos(2*theta))*&
&           (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta))/&
&         ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&           (a**2 + 2*r**2 + a**2*Cos(2*theta))**4) - &
&        (96*a**2*M**2*r**2*Cos(theta)*(-2*a**2 + (2*M - r)*r + a**2*Cos(2*theta))*&
&           Sin(theta))/&
&         ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&           (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (8*a**2*M**2*(2*M - 2*r)*r*Cos(theta)*(a**2 - 6*r**2 + a**2*Cos(2*theta))*&
&           Sin(theta))/&
&         ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&           (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (16*a**2*M**2*r**2*Cos(theta)*(-2*a**2 + (2*M - r)*r + a**2*Cos(2*theta))*&
&           (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta))/&
&         ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**2*&
&           (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (8*a**2*M**2*r*(-2*M + 2*r)*Cos(theta)*&
&           (-2*a**2 + (2*M - r)*r + a**2*Cos(2*theta))*&
&           (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta))/&
&         ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)*&
&           (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (8*a**2*M**2*Cos(theta)*(-2*a**2 + (2*M - r)*r + a**2*Cos(2*theta))*&
&           (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta))/&
&         ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&           (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (24*a**2*M*r*(3*a**2 + r*(-2*M + 3*r))*&
&           (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&           (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) - &
&        (24*a**2*M*r*(3*a**2 + r*(-2*M + 3*r))*&
&           (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (2*a**2*M*(3*a**2 + r*(-2*M + 3*r))*&
&           ((-4*M*r**2)/(r**2 + a**2*Cos(theta)**2)**2 + &
&             (2*M)/(r**2 + a**2*Cos(theta)**2))*(a**2 - 6*r**2 + a**2*Cos(2*theta))*&
&           Sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3)&
&         + (2*a**2*M*(-2*M + 6*r)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&           (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (2*a**2*M*(-2*M + 2*r)*(3*a**2 + r*(-2*M + 3*r))*&
&           (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&           (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3)) - &
&     (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&        ((a*M*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&             (a**2 - 6*r**2 + a**2*Cos(2*theta))*&
&             (a**2*(2*M - 9*r) - 9*a**2*r - 16*r**3 + &
&               a**2*(-2*M + 2*r)*Cos(2*theta))*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&          (12*a*M*r*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&             (a**2 - 6*r**2 + a**2*Cos(2*theta))*&
&             (-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + &
&               a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) - &
&          (12*a*M*r*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&             (-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + &
&               a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&          (a*M*((-4*M*r**2)/(r**2 + a**2*Cos(theta)**2)**2 + &
&               (2*M)/(r**2 + a**2*Cos(theta)**2))*(a**2 - 6*r**2 + a**2*Cos(2*theta))*&
&             (-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + &
&               a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&          (a*M*(-2*M + 2*r)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&             (a**2 - 6*r**2 + a**2*Cos(2*theta))*&
&             (-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + &
&               a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&          (48*a**3*M**2*r**2*(3*a**2 + r*(-2*M + 3*r))*&
&             (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**2*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**4) - &
&          (48*a**3*M**2*r**2*(3*a**2 + r*(-2*M + 3*r))*Sin(theta)**2*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&          (8*a**3*M**2*r**2*(3*a**2 + r*(-2*M + 3*r))*&
&             (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**2*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**2*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&          (4*a**3*M**2*r*(-2*M + 6*r)*(a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**2*&
&             Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&          (4*a**3*M**2*r*(-2*M + 2*r)*(3*a**2 + r*(-2*M + 3*r))*&
&             (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**2*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&          (4*a**3*M**2*(3*a**2 + r*(-2*M + 3*r))*(a**2 - 6*r**2 + a**2*Cos(2*theta))*&
&             Sin(theta)**2*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3)))/&
&      ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&        Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))&
&       + (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&        ((4*a*l*M*r**2*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)**2 - &
&          (2*a*l*M*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(2*r - (4*a**2*M*r**2*Sin(theta)**2)/&
&              (r**2 + a**2*Cos(theta)**2)**2 + &
&             (2*a**2*M*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))*&
&        ((a*M*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&             (a**2 - 6*r**2 + a**2*Cos(2*theta))*&
&             (-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + &
&               a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&          (4*a**3*M**2*r*(3*a**2 + r*(-2*M + 3*r))*&
&             (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**2*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3)))/&
&      ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&         Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))&
&         )**2 - (l*(-(l*((-4*M*r**2)/(r**2 + a**2*Cos(theta)**2)**2 + &
&               (2*M)/(r**2 + a**2*Cos(theta)**2))) - &
&          (4*a*M*r**2*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)**2 + &
&          (2*a*M*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&        ((a*M*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&             (a**2 - 6*r**2 + a**2*Cos(2*theta))*&
&             (-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + &
&               a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&          (4*a**3*M**2*r*(3*a**2 + r*(-2*M + 3*r))*&
&             (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**2*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3)))/&
&      ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&        Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))&
&     ))/(gamma*kappa) + (s0*((-8*r*(r**2 + a**2*Cos(theta)**2)*(1/(Sin(theta)))**2)/&
&      (a**2 + 2*r**2 + a**2*Cos(2*theta))**3 + &
&     (2*r*(1/(Sin(theta)))**2)/(a**2 + 2*r**2 + a**2*Cos(2*theta))**2)*&
&   ((24*a**2*M*(a**2 + r**2)*(a**2 + r*(-2*M + r))*Cos(theta)*&
&        (r**2 + a**2*Cos(theta)**2)*(a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**3*&
&        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4*&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) + &
&     (2*a*M*(a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&        (-3*a**4 + 22*a**2*r**2 + 12*r**4 - 2*a**2*(a**2 + 7*r**2)*Cos(2*theta) + &
&          a**4*Cos(4*theta))*Sin(2*theta))/&
&      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) - &
&     l*((8*a**2*M**2*r*Cos(theta)*(-2*a**2 + (2*M - r)*r + a**2*Cos(2*theta))*&
&           (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta))/&
&         ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&           (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (2*a**2*M*(3*a**2 + r*(-2*M + 3*r))*&
&           (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&           (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3)) - &
&     (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&        ((a*M*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&             (a**2 - 6*r**2 + a**2*Cos(2*theta))*&
&             (-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + &
&               a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&          (4*a**3*M**2*r*(3*a**2 + r*(-2*M + 3*r))*&
&             (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**2*Sin(2*theta))/&
&           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3)))/&
&      ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&        Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))&
&     ))/(gamma*kappa*Sqrt(((r**2 + a**2*Cos(theta)**2)*(1/(Sin(theta)))**2)/&
&     (a**2 + 2*r**2 + a**2*Cos(2*theta))**2))

end subroutine drF2ast_f
