
subroutine F2ast_f(F2ast_func, M, a, r, theta, l, s0, kappa, gamma)

    
    double precision, intent(in) :: M, a, r, theta, l, s0, kappa, gamma
    double precision, intent(out) :: F2ast_func


    F2ast_func = (2*s0*Sqrt(((r**2 + a**2*Cos(theta)**2)*(1/(Sin(theta)))**2)/ &
    &      (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)* &
    &    ((24*a**2*M*(a**2 + r**2)*(a**2 + r*(-2*M + r))*Cos(theta)*(r**2 + a**2*Cos(theta)**2)* &
    &         (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**3* &
    &         (-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) +  &
    &           (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/ &
    &       ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4* &
    &         ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) +  &
    &           Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)) &
    &           )) + (2*a*M*(a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)* &
    &         (-3*a**4 + 22*a**2*r**2 + 12*r**4 - 2*a**2*(a**2 + 7*r**2)*Cos(2*theta) +  &
    &           a**4*Cos(4*theta))*Sin(2*theta))/ &
    &       ((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) -  &
    &      l*((8*a**2*M**2*r*Cos(theta)*(-2*a**2 + (2*M - r)*r + a**2*Cos(2*theta))* &
    &            (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta))/ &
    &          ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)* &
    &            (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) +  &
    &         (2*a**2*M*(3*a**2 + r*(-2*M + 3*r))* &
    &            (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*(a**2 - 6*r**2 + a**2*Cos(2*theta))* &
    &            Sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3)) &
    &       - (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) +  &
    &           (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))* &
    &         ((a*M*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))* &
    &            (a**2 - 6*r**2 + a**2*Cos(2*theta))* &
    &           (-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 +  &
    &           a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(2*theta))/ &
    &            ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) +  &
    &           (4*a**3*M**2*r*(3*a**2 + r*(-2*M + 3*r))* &
    &          (a**2 - 6*r**2 + a**2*Cos(2*theta))*Sin(theta)**2*Sin(2*theta))/ &
    &            ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)* &
    &         (a**2 + 2*r**2 + a**2*Cos(2*theta))**3)))/ &
    &       ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) +  &
    &         Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) &
    &    )/(gamma*kappa) 

end subroutine F2ast_f


