
subroutine dthetaGast_f(dthetaGast_func, M, a, r, theta, l, s0, kappa, gamma)
    
    double precision, intent(in) :: M, a, r, theta, l, s0, kappa, gamma
    double precision, intent(out) :: dthetaGast_func


    dthetaGast_func = (16*a**2*M*r*s0*Cos(theta)*Sqrt(((r**2 + a**2*Cos(theta)**2)*(1/(Sin(theta)))**2)/&
&     (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)*Sin(theta)*&
&   (1 - (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&        Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))&
&     )*(((-2*a*l*M*r*Sin(theta)**2*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&        (Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)))*&
&      ((4*Sqrt(2.)*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           (r**4*(-2*M + r) + a**4*r*Cos(theta)**4 - a**2*M*r**2*Sin(theta)**2 + &
&             Cos(theta)**2*(2*a**2*r**2*(-M + r) + a**4*M*Sin(theta)**2))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (2*Sqrt(2.)*a*M*(r**2 - a**2*Cos(theta)**2)*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        ((3*a**4 + 8*a**2*M*r + 8*a**2*r**2 + 8*r**4 + &
&             4*a**2*(a**2 + 2*r*(-M + r))*Cos(2*theta) + a**4*Cos(4*theta))*(1/(Tan(theta)))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         (Sqrt(2.)*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (4*Sqrt(2.)*a*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*(1/(Tan(theta)))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3)) + &
&     (l*Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&             (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&              (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&           (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&             (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&             Sin(theta)**2*(a**2 + r**2 + &
&                (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) - &
&        (2*a*M*r*Sin(theta)**2*Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/&
&                (r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))))*&
&      ((2*Sqrt(2.)*a*M*(a**4 - 3*a**2*r**2 - 6*r**4 + &
&             a**2*(a**2 - r**2)*Cos(2*theta))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**2*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (Sqrt(2.)*M*(a**2 + r**2)*(a**2 - 2*r**2 + a**2*Cos(2*theta))*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (8*Sqrt(2.)*a**3*M*r*Cos(theta)*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**3*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (2*Sqrt(2.)*a**2*M*r*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3))))/&
& ((a**2 - 2*M*r + r**2)*Gamma*kappa*(r**2 + a**2*Cos(theta)**2)) - &
&(8*a**2*s0*Cos(theta)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&   Sqrt(((r**2 + a**2*Cos(theta)**2)*(1/(Sin(theta)))**2)/(a**2 + 2*r**2 + a**2*Cos(2*theta))**2)*&
&   Sin(theta)*(1 - (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&        Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))&
&     )*(((-2*a*l*M*r*Sin(theta)**2*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&        (Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)))*&
&      ((4*Sqrt(2.)*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           (r**4*(-2*M + r) + a**4*r*Cos(theta)**4 - a**2*M*r**2*Sin(theta)**2 + &
&             Cos(theta)**2*(2*a**2*r**2*(-M + r) + a**4*M*Sin(theta)**2))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (2*Sqrt(2.)*a*M*(r**2 - a**2*Cos(theta)**2)*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        ((3*a**4 + 8*a**2*M*r + 8*a**2*r**2 + 8*r**4 + &
&             4*a**2*(a**2 + 2*r*(-M + r))*Cos(2*theta) + a**4*Cos(4*theta))*(1/(Tan(theta)))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         (Sqrt(2.)*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (4*Sqrt(2.)*a*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*(1/(Tan(theta)))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3)) + &
&     (l*Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&             (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&              (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&           (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&             (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&             Sin(theta)**2*(a**2 + r**2 + &
&                (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) - &
&        (2*a*M*r*Sin(theta)**2*Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/&
&                (r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))))*&
&      ((2*Sqrt(2.)*a*M*(a**4 - 3*a**2*r**2 - 6*r**4 + &
&             a**2*(a**2 - r**2)*Cos(2*theta))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**2*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (Sqrt(2.)*M*(a**2 + r**2)*(a**2 - 2*r**2 + a**2*Cos(2*theta))*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (8*Sqrt(2.)*a**3*M*r*Cos(theta)*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**3*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (2*Sqrt(2.)*a**2*M*r*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3))))/&
& ((a**2 - 2*M*r + r**2)*Gamma*kappa) + &
&(4*s0*(r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&   Sqrt(((r**2 + a**2*Cos(theta)**2)*(1/(Sin(theta)))**2)/(a**2 + 2*r**2 + a**2*Cos(2*theta))**2)*&
&   (-((l*((-4*a**2*l*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2)**2 + &
&            (4*a*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&            (4*a**3*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2))/&
&        ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) + &
&     (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&        ((-4*a*l*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) - &
&          (4*a**3*l*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2 + &
&          2*Cos(theta)*Sin(theta)*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)) + &
&          Sin(theta)**2*((4*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&             (4*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2)))/&
&      ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&         Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))&
&         )**2)*(((-2*a*l*M*r*Sin(theta)**2*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&        (Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)))*&
&      ((4*Sqrt(2.)*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           (r**4*(-2*M + r) + a**4*r*Cos(theta)**4 - a**2*M*r**2*Sin(theta)**2 + &
&             Cos(theta)**2*(2*a**2*r**2*(-M + r) + a**4*M*Sin(theta)**2))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (2*Sqrt(2.)*a*M*(r**2 - a**2*Cos(theta)**2)*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        ((3*a**4 + 8*a**2*M*r + 8*a**2*r**2 + 8*r**4 + &
&             4*a**2*(a**2 + 2*r*(-M + r))*Cos(2*theta) + a**4*Cos(4*theta))*(1/(Tan(theta)))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         (Sqrt(2.)*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (4*Sqrt(2.)*a*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*(1/(Tan(theta)))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3)) + &
&     (l*Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&             (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&              (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&           (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&             (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&             Sin(theta)**2*(a**2 + r**2 + &
&                (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) - &
&        (2*a*M*r*Sin(theta)**2*Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/&
&                (r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))))*&
&      ((2*Sqrt(2.)*a*M*(a**4 - 3*a**2*r**2 - 6*r**4 + &
&             a**2*(a**2 - r**2)*Cos(2*theta))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**2*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (Sqrt(2.)*M*(a**2 + r**2)*(a**2 - 2*r**2 + a**2*Cos(2*theta))*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (8*Sqrt(2.)*a**3*M*r*Cos(theta)*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**3*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (2*Sqrt(2.)*a**2*M*r*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3))))/&
& ((a**2 - 2*M*r + r**2)*Gamma*kappa) + &
&(2*s0*(r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&   (1 - (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&        Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))&
&     )*((-2*a**2*(1/(Tan(theta))))/(a**2 + 2*r**2 + a**2*Cos(2*theta))**2 - &
&     (2*(r**2 + a**2*Cos(theta)**2)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)/&
&      (a**2 + 2*r**2 + a**2*Cos(2*theta))**2 + &
&     (4*a**2*(r**2 + a**2*Cos(theta)**2)*(1/(Sin(theta)))**2*Sin(2*theta))/&
&      (a**2 + 2*r**2 + a**2*Cos(2*theta))**3)*&
&   (((-2*a*l*M*r*Sin(theta)**2*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&        (Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)))*&
&      ((4*Sqrt(2.)*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           (r**4*(-2*M + r) + a**4*r*Cos(theta)**4 - a**2*M*r**2*Sin(theta)**2 + &
&             Cos(theta)**2*(2*a**2*r**2*(-M + r) + a**4*M*Sin(theta)**2))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (2*Sqrt(2.)*a*M*(r**2 - a**2*Cos(theta)**2)*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        ((3*a**4 + 8*a**2*M*r + 8*a**2*r**2 + 8*r**4 + &
&             4*a**2*(a**2 + 2*r*(-M + r))*Cos(2*theta) + a**4*Cos(4*theta))*(1/(Tan(theta)))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         (Sqrt(2.)*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (4*Sqrt(2.)*a*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*(1/(Tan(theta)))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3)) + &
&     (l*Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&             (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&              (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&           (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&             (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&             Sin(theta)**2*(a**2 + r**2 + &
&                (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) - &
&        (2*a*M*r*Sin(theta)**2*Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/&
&                (r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))))*&
&      ((2*Sqrt(2.)*a*M*(a**4 - 3*a**2*r**2 - 6*r**4 + &
&             a**2*(a**2 - r**2)*Cos(2*theta))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**2*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (Sqrt(2.)*M*(a**2 + r**2)*(a**2 - 2*r**2 + a**2*Cos(2*theta))*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (8*Sqrt(2.)*a**3*M*r*Cos(theta)*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**3*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (2*Sqrt(2.)*a**2*M*r*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3))))/&
& ((a**2 - 2*M*r + r**2)*Gamma*kappa*&
&   Sqrt(((r**2 + a**2*Cos(theta)**2)*(1/(Sin(theta)))**2)/(a**2 + 2*r**2 + a**2*Cos(2*theta))**2))&
& + (4*s0*(r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&   Sqrt(((r**2 + a**2*Cos(theta)**2)*(1/(Sin(theta)))**2)/(a**2 + 2*r**2 + a**2*Cos(2*theta))**2)*&
&   (1 - (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&          (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&      ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&        Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))&
&     )*(((-4*a*l*M*r*Cos(theta)*Sin(theta)*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&        (8*a**3*l*M**2*r**2*Cos(theta)*Sin(theta)**3*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)**3*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))**2)&
&         - (4*a**3*l*M*r*Cos(theta)*Sin(theta)**3*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&        (2*Cos(theta)*Sin(theta)*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&        (4*a**2*M*r*Cos(theta)*Sin(theta)**3*&
&           (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))**2)&
&         + (Sin(theta)**2*((4*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&             (4*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2)*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&        (a*l*M*r*Sin(theta)**2*&
&           (-((((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&                    (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                     (a**2 + r**2 + &
&                       (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))*&
&                  ((4*a**2*l**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2)**2 - &
&                    (8*a*l*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) - &
&                    (8*a**3*l*M*r*Cos(theta)*Sin(theta)**3)/&
&                     (r**2 + a**2*Cos(theta)**2)**2 + &
&                    2*Cos(theta)*Sin(theta)*&
&                     (a**2 + r**2 + &
&                       (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)) + &
&                    Sin(theta)**2*&
&                     ((4*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&                       (4*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2))&
&                  )/&
&                (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&                   (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&                   Sin(theta)**2*&
&                    (a**2 + r**2 + &
&                      (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))**2) + &
&             ((16*a**2*M**2*r**2*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2 + &
&                (16*a**4*M**2*r**2*Cos(theta)*Sin(theta)**5)/&
&                 (r**2 + a**2*Cos(theta)**2)**3 - &
&                2*Cos(theta)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)*&
&                 (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))&
&                 - (4*a**2*M*r*Cos(theta)*Sin(theta)**3*&
&                   (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))&
&                   )/(r**2 + a**2*Cos(theta)**2)**2 - &
&                (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                 ((4*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&                   (4*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&                (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&                Sin(theta)**2*(a**2 + r**2 + &
&                   (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))))) + &
&        (Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&           (-((((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&                    (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                     (a**2 + r**2 + &
&                       (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))*&
&                  ((4*a**2*l**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2)**2 - &
&                    (8*a*l*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) - &
&                    (8*a**3*l*M*r*Cos(theta)*Sin(theta)**3)/&
&                     (r**2 + a**2*Cos(theta)**2)**2 + &
&                    2*Cos(theta)*Sin(theta)*&
&                     (a**2 + r**2 + &
&                       (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)) + &
&                    Sin(theta)**2*&
&                     ((4*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&                       (4*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2))&
&                  )/&
&                (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&                   (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&                   Sin(theta)**2*&
&                    (a**2 + r**2 + &
&                      (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))**2) + &
&             ((16*a**2*M**2*r**2*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2 + &
&                (16*a**4*M**2*r**2*Cos(theta)*Sin(theta)**5)/&
&                 (r**2 + a**2*Cos(theta)**2)**3 - &
&                2*Cos(theta)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)*&
&                 (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))&
&                 - (4*a**2*M*r*Cos(theta)*Sin(theta)**3*&
&                   (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))&
&                   )/(r**2 + a**2*Cos(theta)**2)**2 - &
&                (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                 ((4*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&                   (4*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&                (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&                Sin(theta)**2*(a**2 + r**2 + &
&                   (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         (2.*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))))))*&
&      ((4*Sqrt(2.)*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           (r**4*(-2*M + r) + a**4*r*Cos(theta)**4 - a**2*M*r**2*Sin(theta)**2 + &
&             Cos(theta)**2*(2*a**2*r**2*(-M + r) + a**4*M*Sin(theta)**2))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (2*Sqrt(2.)*a*M*(r**2 - a**2*Cos(theta)**2)*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        ((3*a**4 + 8*a**2*M*r + 8*a**2*r**2 + 8*r**4 + &
&             4*a**2*(a**2 + 2*r*(-M + r))*Cos(2*theta) + a**4*Cos(4*theta))*(1/(Tan(theta)))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         (Sqrt(2.)*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (4*Sqrt(2.)*a*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*(1/(Tan(theta)))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3)) + &
&     ((-4*a*M*r*Cos(theta)*Sin(theta)*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&        (8*a**3*M**2*r**2*Cos(theta)*Sin(theta)**3*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)**3*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))**2)&
&         - (4*a**3*M*r*Cos(theta)*Sin(theta)**3*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&        (l*(-((((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&                    (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                     (a**2 + r**2 + &
&                       (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))*&
&                  ((4*a**2*l**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2)**2 - &
&                    (8*a*l*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) - &
&                    (8*a**3*l*M*r*Cos(theta)*Sin(theta)**3)/&
&                     (r**2 + a**2*Cos(theta)**2)**2 + &
&                    2*Cos(theta)*Sin(theta)*&
&                     (a**2 + r**2 + &
&                       (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)) + &
&                    Sin(theta)**2*&
&                     ((4*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&                       (4*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2))&
&                  )/&
&                (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&                   (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&                   Sin(theta)**2*&
&                    (a**2 + r**2 + &
&                      (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))**2) + &
&             ((16*a**2*M**2*r**2*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2 + &
&                (16*a**4*M**2*r**2*Cos(theta)*Sin(theta)**5)/&
&                 (r**2 + a**2*Cos(theta)**2)**3 - &
&                2*Cos(theta)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)*&
&                 (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))&
&                 - (4*a**2*M*r*Cos(theta)*Sin(theta)**3*&
&                   (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))&
&                   )/(r**2 + a**2*Cos(theta)**2)**2 - &
&                (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                 ((4*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&                   (4*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&                (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&                Sin(theta)**2*(a**2 + r**2 + &
&                   (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         (2.*Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))))) - &
&        (a*M*r*Sin(theta)**2*(-((((4*a**2*M**2*r**2*Sin(theta)**4)/&
&                     (r**2 + a**2*Cos(theta)**2)**2 - &
&                    (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                     (a**2 + r**2 + &
&                       (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))*&
&                  ((4*a**2*l**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2)**2 - &
&                    (8*a*l*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) - &
&                    (8*a**3*l*M*r*Cos(theta)*Sin(theta)**3)/&
&                     (r**2 + a**2*Cos(theta)**2)**2 + &
&                    2*Cos(theta)*Sin(theta)*&
&                     (a**2 + r**2 + &
&                       (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)) + &
&                    Sin(theta)**2*&
&                     ((4*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&                       (4*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2))&
&                  )/&
&                (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&                   (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&                   Sin(theta)**2*&
&                    (a**2 + r**2 + &
&                      (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))**2) + &
&             ((16*a**2*M**2*r**2*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2 + &
&                (16*a**4*M**2*r**2*Cos(theta)*Sin(theta)**5)/&
&                 (r**2 + a**2*Cos(theta)**2)**3 - &
&                2*Cos(theta)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)*&
&                 (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))&
&                 - (4*a**2*M*r*Cos(theta)*Sin(theta)**3*&
&                   (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))&
&                   )/(r**2 + a**2*Cos(theta)**2)**2 - &
&                (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                 ((4*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&                   (4*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&                (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&                Sin(theta)**2*(a**2 + r**2 + &
&                   (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))))))*&
&      ((2*Sqrt(2.)*a*M*(a**4 - 3*a**2*r**2 - 6*r**4 + &
&             a**2*(a**2 - r**2)*Cos(2*theta))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**2*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (Sqrt(2.)*M*(a**2 + r**2)*(a**2 - 2*r**2 + a**2*Cos(2*theta))*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (8*Sqrt(2.)*a**3*M*r*Cos(theta)*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**3*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (2*Sqrt(2.)*a**2*M*r*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3)) + &
&     (l*Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&             (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&              (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&           (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&             (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&             Sin(theta)**2*(a**2 + r**2 + &
&                (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) - &
&        (2*a*M*r*Sin(theta)**2*Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/&
&                (r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))))*&
&      ((4*Sqrt(2.)*a*M*Cos(theta)*&
&           (a**4 - 3*a**2*r**2 - 6*r**4 + a**2*(a**2 - r**2)*Cos(2*theta))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (2*Sqrt(2.)*a*M*(a**4 - 3*a**2*r**2 - 6*r**4 + &
&             a**2*(a**2 - r**2)*Cos(2*theta))*&
&           (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - 2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**2*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (4*Sqrt(2.)*a**3*M*(a**2 - r**2)*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**2*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (12*Sqrt(2.)*a**3*M*(a**4 - 3*a**2*r**2 - 6*r**4 + &
&             a**2*(a**2 - r**2)*Cos(2*theta))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**2*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) + &
&        (2*Sqrt(2.)*a**2*M*(a**2 + r**2)*(a**2 - 2*r**2 + a**2*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (6*Sqrt(2.)*a**2*M*(a**2 + r**2)*(a**2 - 2*r**2 + a**2*Cos(2*theta))*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) + &
&        (2*Sqrt(2.)*a**2*M*(a**2 + r**2)*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (Sqrt(2.)*a*M*(a**4 - 3*a**2*r**2 - 6*r**4 + a**2*(a**2 - r**2)*Cos(2*theta))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*Sin(theta)**2*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))*&
&           ((-2*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*(1/(Tan(theta)))*&
&                (1/(Sin(theta)))**2)/&
&              (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                r*(1/(Sin(theta)))**2*&
&                 (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)) - &
&             ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*(1/(Sin(theta)))**2*&
&                (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                  2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                  2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                  2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))/&
&              (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                 r*(1/(Sin(theta)))**2*&
&                  (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**2 - &
&             (2*a**2*(a**2 + r*(-2*M + r))*(1/(Sin(theta)))**2*Sin(2*theta))/&
&              (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                r*(1/(Sin(theta)))**2*&
&                 (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))) + &
&        (24*Sqrt(2.)*a**3*M*r*Cos(theta)**2*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**2*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (8*Sqrt(2.)*a**3*M*r*Cos(theta)*&
&           (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - 2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**3*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (8*Sqrt(2.)*a**3*M*r*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**4*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (4*Sqrt(2.)*a**2*M*r*Cos(2*theta)*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (48*Sqrt(2.)*a**5*M*r*Cos(theta)*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**3*Sin(2*theta)*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) + &
&        (4*Sqrt(2.)*a**4*M*r*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)**2*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         (a**2 + 2*r**2 + a**2*Cos(2*theta))**3 - &
&        (12*Sqrt(2.)*a**4*M*r*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)**2*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) + &
&        (4*Sqrt(2.)*a**3*M*r*Cos(theta)*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*Sin(theta)**3*&
&           ((-2*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*(1/(Tan(theta)))*&
&                (1/(Sin(theta)))**2)/&
&              (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                r*(1/(Sin(theta)))**2*&
&                 (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)) - &
&             ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*(1/(Sin(theta)))**2*&
&                (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                  2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                  2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                  2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))/&
&              (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                 r*(1/(Sin(theta)))**2*&
&                  (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**2 - &
&             (2*a**2*(a**2 + r*(-2*M + r))*(1/(Sin(theta)))**2*Sin(2*theta))/&
&              (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                r*(1/(Sin(theta)))**2*&
&                 (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))) - &
&        (M*(a**2 + r**2)*(a**2 - 2*r**2 + a**2*Cos(2*theta))*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))*&
&           ((2*(a**2 + r*(-2*M + r))*Cos(theta)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                Sin(theta))/&
&              (l**2*(2*M - r)*r + &
&                r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)) - &
&             ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2*&
&                (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                  4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                  2*r*Cos(theta)*Sin(theta)*&
&                   (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                  2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&              (l**2*(2*M - r)*r + &
&                 r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                 a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2 - &
&             (2*a**2*(a**2 + r*(-2*M + r))*Sin(theta)**2*Sin(2*theta))/&
&              (l**2*(2*M - r)*r + &
&                r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&         (Sqrt(2.)*(a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))) - &
&        (Sqrt(2.)*a**2*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(2*theta)*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2))*&
&           ((2*(a**2 + r*(-2*M + r))*Cos(theta)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                Sin(theta))/&
&              (l**2*(2*M - r)*r + &
&                r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)) - &
&             ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2*&
&                (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                  4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                  2*r*Cos(theta)*Sin(theta)*&
&                   (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                  2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&              (l**2*(2*M - r)*r + &
&                 r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                 a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2 - &
&             (2*a**2*(a**2 + r*(-2*M + r))*Sin(theta)**2*Sin(2*theta))/&
&              (l**2*(2*M - r)*r + &
&                r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))) + &
&        (2*Sqrt(2.)*a*M*(a**4 - 3*a**2*r**2 - 6*r**4 + &
&             a**2*(a**2 - r**2)*Cos(2*theta))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**2*&
&           ((2*a**2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*(1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                  l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             ((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                  2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                  2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                  2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**2)&
&              + ((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-4*r*Cos(theta)*Sin(theta) + &
&                  (8*a**2*M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**3)/&
&                   (r**2 + a**2*Cos(theta)**2)**2 + &
&                  (8*a**4*M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**5)/&
&                   (r**2 + a**2*Cos(theta)**2)**3 + &
&                  (4*a**4*M*Cos(theta)*Sin(theta)**5)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (4*a*M*Cos(theta)*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (12*a**3*M*Cos(theta)*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**3*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**4*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) - &
&             (4*a**3*M*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**3*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (3*a**2*M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((r**2 + a**2*Cos(theta)**2)**4*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2*&
&                (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                  4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                  2*r*Cos(theta)*Sin(theta)*&
&                   (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                  2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2) + &
&             (2*a**2*(1/(Sin(theta)))**2*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2)*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (4*a**3*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))**2)&
&              + (8*a**3*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&             (2*a**2*M*(a**2 + r*(-2*M + r))*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**2*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&                ((-2*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Tan(theta)))*(1/(Sin(theta)))**2)/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Sin(theta)))**2*&
&                     (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                       2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                       2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                       2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                        (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))&
&                    /(2*a**2*M*r + &
&                      a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                      r*(1/(Sin(theta)))**2*&
&                       (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**&
&                    2 - (2*a**2*(a**2 + r*(-2*M + r))*(1/(Sin(theta)))**2*Sin(2*theta))/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))&
&              + (a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*((2*(a**2 + r*(-2*M + r))*Cos(theta)*&
&                     (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     Sin(theta)**2*&
&                     (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                       4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                       2*r*Cos(theta)*Sin(theta)*&
&                        (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                       2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&                   (l**2*(2*M - r)*r + &
&                      r*Sin(theta)**2*&
&                       (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                      a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2 - &
&                  (2*a**2*(a**2 + r*(-2*M + r))*Sin(theta)**2*Sin(2*theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (Sqrt(2.)*M*(a**2 + r**2)*(a**2 - 2*r**2 + a**2*Cos(2*theta))*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           ((2*a**2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*(1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                  l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             ((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                  2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                  2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                  2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**2)&
&              + ((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-4*r*Cos(theta)*Sin(theta) + &
&                  (8*a**2*M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**3)/&
&                   (r**2 + a**2*Cos(theta)**2)**2 + &
&                  (8*a**4*M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**5)/&
&                   (r**2 + a**2*Cos(theta)**2)**3 + &
&                  (4*a**4*M*Cos(theta)*Sin(theta)**5)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (4*a*M*Cos(theta)*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (12*a**3*M*Cos(theta)*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**3*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**4*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) - &
&             (4*a**3*M*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**3*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (3*a**2*M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((r**2 + a**2*Cos(theta)**2)**4*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2*&
&                (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                  4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                  2*r*Cos(theta)*Sin(theta)*&
&                   (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                  2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2) + &
&             (2*a**2*(1/(Sin(theta)))**2*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2)*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (4*a**3*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))**2)&
&              + (8*a**3*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&             (2*a**2*M*(a**2 + r*(-2*M + r))*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**2*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&                ((-2*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Tan(theta)))*(1/(Sin(theta)))**2)/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Sin(theta)))**2*&
&                     (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                       2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                       2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                       2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                        (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))&
&                    /(2*a**2*M*r + &
&                      a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                      r*(1/(Sin(theta)))**2*&
&                       (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**&
&                    2 - (2*a**2*(a**2 + r*(-2*M + r))*(1/(Sin(theta)))**2*Sin(2*theta))/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))&
&              + (a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*((2*(a**2 + r*(-2*M + r))*Cos(theta)*&
&                     (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     Sin(theta)**2*&
&                     (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                       4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                       2*r*Cos(theta)*Sin(theta)*&
&                        (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                       2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&                   (l**2*(2*M - r)*r + &
&                      r*Sin(theta)**2*&
&                       (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                      a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2 - &
&                  (2*a**2*(a**2 + r*(-2*M + r))*Sin(theta)**2*Sin(2*theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (8*Sqrt(2.)*a**3*M*r*Cos(theta)*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(theta)**3*&
&           ((4*a**2*Cos(theta)**2*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**2*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (2*(1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (4*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                  2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                  2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                  2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**2)&
&              + (2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                ((-8*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) - &
&                  (16*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2 - &
&                  (8*a**6*M*r*Cos(theta)*Sin(theta)**5)/(r**2 + a**2*Cos(theta)**2)**3))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (4*a*M*r*(a**2 + r**2)*Cos(2*theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) - &
&             (3*a**2*M*r*Cos(theta)**2*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (6*a**4*M*r*Cos(theta)**2*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**4)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**4*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**4)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3*&
&                (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                  4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                  2*r*Cos(theta)*Sin(theta)*&
&                   (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                  2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2) + &
&             (4*a**2*(1/(Tan(theta)))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2)*Sin(2*theta))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (12*a**3*M*r*(a**2 + r**2)*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**4*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (4*a**4*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**3*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (2*a**4*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3*Sin(2*theta))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (4*a**3*M*r*(a**2 + r**2)*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&                Sin(2*theta)**2)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (8*a**3*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&                Sin(2*theta)**2)/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&             (a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)*&
&                ((-2*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Tan(theta)))*(1/(Sin(theta)))**2)/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Sin(theta)))**2*&
&                     (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                       2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                       2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                       2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                        (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))&
&                    /(2*a**2*M*r + &
&                      a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                      r*(1/(Sin(theta)))**2*&
&                       (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**&
&                    2 - (2*a**2*(a**2 + r*(-2*M + r))*(1/(Sin(theta)))**2*Sin(2*theta))/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))&
&              + (a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(2*theta)*((2*(a**2 + r*(-2*M + r))*Cos(theta)*&
&                     (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     Sin(theta)**2*&
&                     (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                       4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                       2*r*Cos(theta)*Sin(theta)*&
&                        (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                       2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&                   (l**2*(2*M - r)*r + &
&                      r*Sin(theta)**2*&
&                       (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                      a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2 - &
&                  (2*a**2*(a**2 + r*(-2*M + r))*Sin(theta)**2*Sin(2*theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (2*Sqrt(2.)*a**2*M*r*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)*&
&           ((4*a**2*Cos(theta)**2*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**2*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (2*(1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (4*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                  2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                  2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                  2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**2)&
&              + (2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                ((-8*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) - &
&                  (16*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2 - &
&                  (8*a**6*M*r*Cos(theta)*Sin(theta)**5)/(r**2 + a**2*Cos(theta)**2)**3))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (4*a*M*r*(a**2 + r**2)*Cos(2*theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) - &
&             (3*a**2*M*r*Cos(theta)**2*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (6*a**4*M*r*Cos(theta)**2*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**4)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**4*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**4)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3*&
&                (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                  4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                  2*r*Cos(theta)*Sin(theta)*&
&                   (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                  2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2) + &
&             (4*a**2*(1/(Tan(theta)))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2)*Sin(2*theta))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (12*a**3*M*r*(a**2 + r**2)*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**4*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (4*a**4*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**3*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (2*a**4*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3*Sin(2*theta))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (4*a**3*M*r*(a**2 + r**2)*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&                Sin(2*theta)**2)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (8*a**3*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&                Sin(2*theta)**2)/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&             (a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)*&
&                ((-2*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Tan(theta)))*(1/(Sin(theta)))**2)/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Sin(theta)))**2*&
&                     (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                       2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                       2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                       2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                        (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))&
&                    /(2*a**2*M*r + &
&                      a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                      r*(1/(Sin(theta)))**2*&
&                       (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**&
&                    2 - (2*a**2*(a**2 + r*(-2*M + r))*(1/(Sin(theta)))**2*Sin(2*theta))/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))&
&              + (a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(2*theta)*((2*(a**2 + r*(-2*M + r))*Cos(theta)*&
&                     (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     Sin(theta)**2*&
&                     (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                       4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                       2*r*Cos(theta)*Sin(theta)*&
&                        (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                       2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&                   (l**2*(2*M - r)*r + &
&                      r*Sin(theta)**2*&
&                       (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                      a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2 - &
&                  (2*a**2*(a**2 + r*(-2*M + r))*Sin(theta)**2*Sin(2*theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3)) + &
&     ((-2*a*l*M*r*Sin(theta)**2*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&        (Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&           Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&               (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&                (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&             (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&               (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&               Sin(theta)**2*(a**2 + r**2 + &
&                  (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&         (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)))*&
&      ((4*Sqrt(2.)*(-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&             2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           (r**4*(-2*M + r) + a**4*r*Cos(theta)**4 - a**2*M*r**2*Sin(theta)**2 + &
&             Cos(theta)**2*(2*a**2*r**2*(-M + r) + a**4*M*Sin(theta)**2))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (4*Sqrt(2.)*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           (-2*a**2*M*r**2*Cos(theta)*Sin(theta) + 2*a**4*M*Cos(theta)**3*Sin(theta) - &
&             4*a**4*r*Cos(theta)**3*Sin(theta) - &
&             2*Cos(theta)*Sin(theta)*(2*a**2*r**2*(-M + r) + a**4*M*Sin(theta)**2))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (4*Sqrt(2.)*a**3*M*Cos(theta)*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (24*Sqrt(2.)*a**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&             l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           (r**4*(-2*M + r) + a**4*r*Cos(theta)**4 - a**2*M*r**2*Sin(theta)**2 + &
&             Cos(theta)**2*(2*a**2*r**2*(-M + r) + a**4*M*Sin(theta)**2))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) - &
&        (4*Sqrt(2.)*a**3*M*(r**2 - a**2*Cos(theta)**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (12*Sqrt(2.)*a**3*M*(r**2 - a**2*Cos(theta)**2)*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))*Sin(2*theta))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) + &
&        (2*Sqrt(2.)*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           (r**4*(-2*M + r) + a**4*r*Cos(theta)**4 - a**2*M*r**2*Sin(theta)**2 + &
&             Cos(theta)**2*(2*a**2*r**2*(-M + r) + a**4*M*Sin(theta)**2))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))*&
&           ((-2*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*(1/(Tan(theta)))*&
&                (1/(Sin(theta)))**2)/&
&              (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                r*(1/(Sin(theta)))**2*&
&                 (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)) - &
&             ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*(1/(Sin(theta)))**2*&
&                (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                  2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                  2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                  2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))/&
&              (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                 r*(1/(Sin(theta)))**2*&
&                  (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**2 - &
&             (2*a**2*(a**2 + r*(-2*M + r))*(1/(Sin(theta)))**2*Sin(2*theta))/&
&              (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                r*(1/(Sin(theta)))**2*&
&                 (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))) - &
&        ((3*a**4 + 8*a**2*M*r + 8*a**2*r**2 + 8*r**4 + &
&             4*a**2*(a**2 + 2*r*(-M + r))*Cos(2*theta) + a**4*Cos(4*theta))*(1/(Sin(theta)))**2*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         (Sqrt(2.)*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        ((3*a**4 + 8*a**2*M*r + 8*a**2*r**2 + 8*r**4 + &
&             4*a**2*(a**2 + 2*r*(-M + r))*Cos(2*theta) + a**4*Cos(4*theta))*(1/(Tan(theta)))*&
&           (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - 2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         (Sqrt(2.)*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (4*Sqrt(2.)*a*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*(1/(Sin(theta)))**2*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (3*Sqrt(2.)*a**2*(3*a**4 + 8*a**2*M*r + 8*a**2*r**2 + 8*r**4 + &
&             4*a**2*(a**2 + 2*r*(-M + r))*Cos(2*theta) + a**4*Cos(4*theta))*(1/(Tan(theta)))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*Sin(2*theta)*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) + &
&        (8*Sqrt(2.)*a**3*M*r*(1/(Tan(theta)))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         (a**2 + 2*r**2 + a**2*Cos(2*theta))**3 - &
&        (24*Sqrt(2.)*a**3*M*r*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*(1/(Tan(theta)))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**4) + &
&        ((3*a**4 + 8*a**2*M*r + 8*a**2*r**2 + 8*r**4 + &
&             4*a**2*(a**2 + 2*r*(-M + r))*Cos(2*theta) + a**4*Cos(4*theta))*(1/(Tan(theta)))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           ((-2*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*(1/(Tan(theta)))*&
&                (1/(Sin(theta)))**2)/&
&              (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                r*(1/(Sin(theta)))**2*&
&                 (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)) - &
&             ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*(1/(Sin(theta)))**2*&
&                (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                  2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                  2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                  2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))/&
&              (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                 r*(1/(Sin(theta)))**2*&
&                  (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**2 - &
&             (2*a**2*(a**2 + r*(-2*M + r))*(1/(Sin(theta)))**2*Sin(2*theta))/&
&              (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                r*(1/(Sin(theta)))**2*&
&                 (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&         (2.*Sqrt(2.)*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))) + &
&        (Sqrt(2.)*a*M*(r**2 - a**2*Cos(theta)**2)*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))*&
&           ((2*(a**2 + r*(-2*M + r))*Cos(theta)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                Sin(theta))/&
&              (l**2*(2*M - r)*r + &
&                r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)) - &
&             ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2*&
&                (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                  4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                  2*r*Cos(theta)*Sin(theta)*&
&                   (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                  2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&              (l**2*(2*M - r)*r + &
&                 r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                 a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2 - &
&             (2*a**2*(a**2 + r*(-2*M + r))*Sin(theta)**2*Sin(2*theta))/&
&              (l**2*(2*M - r)*r + &
&                r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))) - &
&        (2*Sqrt(2.)*a*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*(1/(Tan(theta)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2))*&
&           ((2*(a**2 + r*(-2*M + r))*Cos(theta)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                Sin(theta))/&
&              (l**2*(2*M - r)*r + &
&                r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)) - &
&             ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2*&
&                (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                  4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                  2*r*Cos(theta)*Sin(theta)*&
&                   (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                  2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&              (l**2*(2*M - r)*r + &
&                 r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                 a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2 - &
&             (2*a**2*(a**2 + r*(-2*M + r))*Sin(theta)**2*Sin(2*theta))/&
&              (l**2*(2*M - r)*r + &
&                r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))) + &
&        (4*Sqrt(2.)*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           (r**4*(-2*M + r) + a**4*r*Cos(theta)**4 - a**2*M*r**2*Sin(theta)**2 + &
&             Cos(theta)**2*(2*a**2*r**2*(-M + r) + a**4*M*Sin(theta)**2))*&
&           ((2*a**2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*(1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                  l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             ((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                  2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                  2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                  2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**2)&
&              + ((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-4*r*Cos(theta)*Sin(theta) + &
&                  (8*a**2*M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**3)/&
&                   (r**2 + a**2*Cos(theta)**2)**2 + &
&                  (8*a**4*M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**5)/&
&                   (r**2 + a**2*Cos(theta)**2)**3 + &
&                  (4*a**4*M*Cos(theta)*Sin(theta)**5)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (4*a*M*Cos(theta)*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (12*a**3*M*Cos(theta)*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**3*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**4*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) - &
&             (4*a**3*M*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**3*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (3*a**2*M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((r**2 + a**2*Cos(theta)**2)**4*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2*&
&                (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                  4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                  2*r*Cos(theta)*Sin(theta)*&
&                   (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                  2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2) + &
&             (2*a**2*(1/(Sin(theta)))**2*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2)*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (4*a**3*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))**2)&
&              + (8*a**3*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&             (2*a**2*M*(a**2 + r*(-2*M + r))*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**2*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&                ((-2*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Tan(theta)))*(1/(Sin(theta)))**2)/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Sin(theta)))**2*&
&                     (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                       2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                       2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                       2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                        (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))&
&                    /(2*a**2*M*r + &
&                      a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                      r*(1/(Sin(theta)))**2*&
&                       (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**&
&                    2 - (2*a**2*(a**2 + r*(-2*M + r))*(1/(Sin(theta)))**2*Sin(2*theta))/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))&
&              + (a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*((2*(a**2 + r*(-2*M + r))*Cos(theta)*&
&                     (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     Sin(theta)**2*&
&                     (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                       4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                       2*r*Cos(theta)*Sin(theta)*&
&                        (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                       2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&                   (l**2*(2*M - r)*r + &
&                      r*Sin(theta)**2*&
&                       (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                      a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2 - &
&                  (2*a**2*(a**2 + r*(-2*M + r))*Sin(theta)**2*Sin(2*theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        (2*Sqrt(2.)*a*M*(r**2 - a**2*Cos(theta)**2)*&
&           (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           ((2*a**2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*(1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                  l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             ((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                  2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                  2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                  2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**2)&
&              + ((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-4*r*Cos(theta)*Sin(theta) + &
&                  (8*a**2*M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**3)/&
&                   (r**2 + a**2*Cos(theta)**2)**2 + &
&                  (8*a**4*M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**5)/&
&                   (r**2 + a**2*Cos(theta)**2)**3 + &
&                  (4*a**4*M*Cos(theta)*Sin(theta)**5)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (4*a*M*Cos(theta)*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (12*a**3*M*Cos(theta)*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**3*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**4*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) - &
&             (4*a**3*M*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**3*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (3*a**2*M*Cos(theta)*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((r**2 + a**2*Cos(theta)**2)**4*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2*&
&                (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                  4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                  2*r*Cos(theta)*Sin(theta)*&
&                   (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                  2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&              (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2) + &
&             (2*a**2*(1/(Sin(theta)))**2*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-2*r*Sin(theta)**2 + &
&                  (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                   (r**2 + a**2*Cos(theta)**2)**2)*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (4*a**3*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))**2)&
&              + (8*a**3*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&             (2*a**2*M*(a**2 + r*(-2*M + r))*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**2*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*(r**2 - a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&                ((-2*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Tan(theta)))*(1/(Sin(theta)))**2)/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Sin(theta)))**2*&
&                     (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                       2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                       2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                       2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                        (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))&
&                    /(2*a**2*M*r + &
&                      a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                      r*(1/(Sin(theta)))**2*&
&                       (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**&
&                    2 - (2*a**2*(a**2 + r*(-2*M + r))*(1/(Sin(theta)))**2*Sin(2*theta))/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))&
&              + (a*M*(-r**2 + a**2*Cos(theta)**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)**2*((2*(a**2 + r*(-2*M + r))*Cos(theta)*&
&                     (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     Sin(theta)**2*&
&                     (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                       4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                       2*r*Cos(theta)*Sin(theta)*&
&                        (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                       2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&                   (l**2*(2*M - r)*r + &
&                      r*Sin(theta)**2*&
&                       (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                      a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2 - &
&                  (2*a**2*(a**2 + r*(-2*M + r))*Sin(theta)**2*Sin(2*theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))))/&
&         ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        ((3*a**4 + 8*a**2*M*r + 8*a**2*r**2 + 8*r**4 + &
&             4*a**2*(a**2 + 2*r*(-M + r))*Cos(2*theta) + a**4*Cos(4*theta))*(1/(Tan(theta)))*&
&           (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           ((4*a**2*Cos(theta)**2*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**2*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (2*(1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (4*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                  2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                  2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                  2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**2)&
&              + (2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                ((-8*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) - &
&                  (16*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2 - &
&                  (8*a**6*M*r*Cos(theta)*Sin(theta)**5)/(r**2 + a**2*Cos(theta)**2)**3))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (4*a*M*r*(a**2 + r**2)*Cos(2*theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) - &
&             (3*a**2*M*r*Cos(theta)**2*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (6*a**4*M*r*Cos(theta)**2*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**4)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**4*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**4)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3*&
&                (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                  4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                  2*r*Cos(theta)*Sin(theta)*&
&                   (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                  2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2) + &
&             (4*a**2*(1/(Tan(theta)))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2)*Sin(2*theta))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (12*a**3*M*r*(a**2 + r**2)*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**4*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (4*a**4*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**3*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (2*a**4*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3*Sin(2*theta))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (4*a**3*M*r*(a**2 + r**2)*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&                Sin(2*theta)**2)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (8*a**3*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&                Sin(2*theta)**2)/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&             (a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)*&
&                ((-2*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Tan(theta)))*(1/(Sin(theta)))**2)/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Sin(theta)))**2*&
&                     (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                       2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                       2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                       2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                        (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))&
&                    /(2*a**2*M*r + &
&                      a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                      r*(1/(Sin(theta)))**2*&
&                       (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**&
&                    2 - (2*a**2*(a**2 + r*(-2*M + r))*(1/(Sin(theta)))**2*Sin(2*theta))/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))&
&              + (a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(2*theta)*((2*(a**2 + r*(-2*M + r))*Cos(theta)*&
&                     (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     Sin(theta)**2*&
&                     (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                       4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                       2*r*Cos(theta)*Sin(theta)*&
&                        (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                       2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&                   (l**2*(2*M - r)*r + &
&                      r*Sin(theta)**2*&
&                       (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                      a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2 - &
&                  (2*a**2*(a**2 + r*(-2*M + r))*Sin(theta)**2*Sin(2*theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))))/&
&         (Sqrt(2.)*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&        (4*Sqrt(2.)*a*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&             a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*(1/(Tan(theta)))*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               Sin(theta)**2)/&
&             (l**2*(2*M - r)*r + &
&               r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&               a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&           ((4*a**2*Cos(theta)**2*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**2*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (2*(1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (4*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                  2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                  2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                  2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**2)&
&              + (2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + &
&                   l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                ((-8*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) - &
&                  (16*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2 - &
&                  (8*a**6*M*r*Cos(theta)*Sin(theta)**5)/(r**2 + a**2*Cos(theta)**2)**3))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (4*a*M*r*(a**2 + r**2)*Cos(2*theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) - &
&             (3*a**2*M*r*Cos(theta)**2*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (6*a**4*M*r*Cos(theta)**2*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**4)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**4*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**4)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3*&
&                (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                  4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                  2*r*Cos(theta)*Sin(theta)*&
&                   (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                  2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2) + &
&             (4*a**2*(1/(Tan(theta)))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&                (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2)*Sin(2*theta))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (-2*a**2*l*(1/(Tan(theta)))*(1/(Sin(theta)))**2 - &
&                  2*l*r*(-2*M + r)*(1/(Tan(theta)))*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (12*a**3*M*r*(a**2 + r**2)*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(theta)*Sqrt(((a**2 + r*(-2*M + r))*&
&                    (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**4*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (4*a**4*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*Sin(theta)**3*Sin(2*theta))/&
&              ((r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (2*a**4*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3*Sin(2*theta))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) - &
&             (4*a**3*M*r*(a**2 + r**2)*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&                Sin(2*theta)**2)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&             (8*a**3*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&                Sin(2*theta)**2)/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&             (a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)*&
&                ((-2*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Tan(theta)))*(1/(Sin(theta)))**2)/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     (1/(Sin(theta)))**2*&
&                     (2*a**2*l**2*(1/(Tan(theta)))**3*(1/(Sin(theta)))**2 - &
&                       2*l**2*(2*M - r)*r*(1/(Tan(theta)))*(1/(Sin(theta)))**4 - &
&                       2*a**2*(1/(Tan(theta)))*(1/(Sin(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) - &
&                       2*r*(1/(Tan(theta)))*(1/(Sin(theta)))**2*&
&                        (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))&
&                    /(2*a**2*M*r + &
&                      a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                      r*(1/(Sin(theta)))**2*&
&                       (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))**&
&                    2 - (2*a**2*(a**2 + r*(-2*M + r))*(1/(Sin(theta)))**2*Sin(2*theta))/&
&                   (2*a**2*M*r + &
&                     a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                     r*(1/(Sin(theta)))**2*&
&                      (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))))&
&              + (a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sin(2*theta)*((2*(a**2 + r*(-2*M + r))*Cos(theta)*&
&                     (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)) - &
&                  ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                     Sin(theta)**2*&
&                     (2*a**2*(a**2 + r**2)*Cos(theta)**3*Sin(theta) + &
&                       4*a**2*M*r*Cos(theta)*Sin(theta)**3 + &
&                       2*r*Cos(theta)*Sin(theta)*&
&                        (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) - &
&                       2*a**2*Cos(theta)*Sin(theta)*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))/&
&                   (l**2*(2*M - r)*r + &
&                      r*Sin(theta)**2*&
&                       (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                      a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))**2 - &
&                  (2*a**2*(a**2 + r*(-2*M + r))*Sin(theta)**2*Sin(2*theta))/&
&                   (l**2*(2*M - r)*r + &
&                     r*Sin(theta)**2*&
&                      (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                     a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&              ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))))/&
&         ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&        ((1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&           Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (1/(Sin(theta)))**2)/&
&             (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&               r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + &
&                  l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&           ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                 2*(-a**2 - r**2 - &
&                  (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                  (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                  r*(1/(Sin(theta)))**2*&
&                   (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&             (a**2*M*r*Cos(theta)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                   a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&              ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                (l**2*(2*M - r)*r + &
&                  r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                  a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&             (2*a*M*r*(a**2 + r**2)*&
&                (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&                (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    (1/(Sin(theta)))**2)/&
&                  (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                    r*(1/(Sin(theta)))**2*&
&                     (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&                Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                    Sin(theta)**2)/&
&                  (l**2*(2*M - r)*r + &
&                    r*Sin(theta)**2*&
&                     (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                    a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))&
&               /((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&                (a**2 + 2*r**2 + a**2*Cos(2*theta))**2))*&
&           (-8*a**2*(a**2 + 2*r*(-M + r))*Sin(2*theta) - 4*a**4*Sin(4*theta)))/&
&         (Sqrt(2.)*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3))))/&
& ((a**2 - 2*M*r + r**2)*Gamma*kappa)

end subroutine dthetaGast_f
