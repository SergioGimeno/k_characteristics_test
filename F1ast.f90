
subroutine F1ast_f(F1ast_func, M, a, r, theta, l, s0, kappa, gamma)

    
    double precision, intent(in) :: M, a, r, theta, l, s0, kappa, gamma
    double precision, intent(out) :: F1ast_func


    F1ast_func = (2*s0*Sqrt(((r**2 + a**2*Cos(theta)**2)*(1/Sin(theta))**2)/ &
    &     (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)* &
    &   ((a*M*r*(3*a**2 + r*(-4*M + 3*r))*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))* &
    &        Sin(theta)**2)/(2.*(a**2 - 2*M*r + r**2)*(r**2 + a**2*Cos(theta)**2)**3) -  &
    &     l*((M*r*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))* &
    &           (-27*a**4 - 4*a**2*r**2 + 16*r**4 -  &
    &             4*a**2*(6*a**2 + 7*r**2)*Cos(2*theta) + 3*a**4*Cos(4*theta)))/ &
    &         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) +  &
    &        (24*a**2*M**2*r**2*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*Sin(theta)**2)/ &
    &         ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)* &
    &           (a**2 + 2*r**2 + a**2*Cos(2*theta))**3)) -  &
    &     (8*M*r*(r**2 + a**2*Cos(theta)**2)*(3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))* &
    &        (2*a**4 + r**4 + a**2*r*(-2*M + 3*r) -  &
    &          a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**2* &
    &        (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) +  &
    &          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/ &
    &      ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4* &
    &        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) +  &
    &          Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)) &
    &          )) - (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) +  &
    &          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))* &
    &        ((12*a*M*r*(a**2 + r**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))* &
    &             (3*a**2 - 2*r**2 + 3*a**2*Cos(2*theta))*Sin(theta)**2)/ &
    &           ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) +  &
    &          (4*a*M**2*r**2*(-9*a**4 + 2*a**2*r**2 + 4*r**4 -  &
    &               2*a**2*(3*a**2 + 5*r**2)*Cos(2*theta) + 3*a**4*Cos(4*theta))*Sin(theta)**2)/ &
    &           ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)* &
    &             (a**2 + 2*r**2 + a**2*Cos(2*theta))**3)))/ &
    &      ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) +  &
    &        Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) &
    &   )/(gamma*kappa)

end subroutine F1ast_f