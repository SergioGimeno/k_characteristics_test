
subroutine Gast_f(Gast_func, M, a, r, theta, l, s0, kappa, gamma)
    
    double precision, intent(in) :: M, a, r, theta, l, s0, kappa, gamma
    double precision, intent(out) :: Gast_func


    Gast_func = (4*s0*(r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*&
&  Sqrt(((r**2 + a**2*Cos(theta)**2)*(1/(Sin(theta)))**2)/(a**2 + 2*r**2 + a**2*Cos(2*theta))**2)*&
&  (1 - (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&         (2*a*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&     ((-2*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&       Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))))&
&   *(((-2*a*l*M*r*Sin(theta)**2*&
&          Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&              (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&               (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&            (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&              (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&              Sin(theta)**2*(a**2 + r**2 + &
&     (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&        ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))) + &
&       (Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))*&
&          Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&              (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&               (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&            (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&              (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&              Sin(theta)**2*(a**2 + r**2 + &
&                 (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&        (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)))*&
&     ((4*Sqrt(2.)*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&          Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*(1/(Sin(theta)))**2)/&
&            (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&              r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)&
&              ))*(r**4*(-2*M + r) + a**4*r*Cos(theta)**4 - a**2*M*r**2*Sin(theta)**2 + &
&            Cos(theta)**2*(2*a**2*r**2*(-M + r) + a**4*M*Sin(theta)**2))*&
&          (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                2*(-2*r*Sin(theta)**2 + &
&                 (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                  (r**2 + a**2*Cos(theta)**2)**2))/&
&             ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                 r*(1/(Sin(theta)))**2*&
&                  (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&            (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                 a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&               (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&               Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                   (1/(Sin(theta)))**2)/&
&                 (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&               Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                   (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                 (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&             ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&            (M*(r**2 - a**2*Cos(theta)**2)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&             (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (l**2*(2*M - r)*r + &
&                 r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                 a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&        ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&       (2*Sqrt(2.)*a*M*(r**2 - a**2*Cos(theta)**2)*&
&          (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&            a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&          Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&            (l**2*(2*M - r)*r + &
&              r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&              a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&          (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                2*(-2*r*Sin(theta)**2 + &
&                 (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                  (r**2 + a**2*Cos(theta)**2)**2))/&
&             ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                 r*(1/(Sin(theta)))**2*&
&                  (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&            (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                 a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&               (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&               Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                   (1/(Sin(theta)))**2)/&
&                 (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&               Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                   (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                 (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&             ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&            (M*(r**2 - a**2*Cos(theta)**2)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&             (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (l**2*(2*M - r)*r + &
&                 r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                 a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&        ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&       ((3*a**4 + 8*a**2*M*r + 8*a**2*r**2 + 8*r**4 + &
&            4*a**2*(a**2 + 2*r*(-M + r))*Cos(2*theta) + a**4*Cos(4*theta))*(1/(Tan(theta)))*&
&          (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&          Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*(1/(Sin(theta)))**2)/&
&            (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&              r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)&
&              ))*((2*(1/(Tan(theta)))*&
&               (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&               (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                 (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&             ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                 r*(1/(Sin(theta)))**2*&
&                  (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&            (a**2*M*r*Cos(theta)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&             ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (l**2*(2*M - r)*r + &
&                 r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                 a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&            (2*a*M*r*(a**2 + r**2)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                 a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&               (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&               Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                   (1/(Sin(theta)))**2)/&
&                 (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&               Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                   Sin(theta)**2)/&
&                 (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))/&
&             ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&        (Sqrt(2.)*(a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&       (4*Sqrt(2.)*a*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&            a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*(1/(Tan(theta)))*&
&          Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&            (l**2*(2*M - r)*r + &
&              r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&              a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&          ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                2*(-a**2 - r**2 - &
&                 (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                 (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&             ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                 r*(1/(Sin(theta)))**2*&
&                  (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&            (a**2*M*r*Cos(theta)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&             ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (l**2*(2*M - r)*r + &
&                 r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                 a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&            (2*a*M*r*(a**2 + r**2)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                 a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&               (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&               Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                   (1/(Sin(theta)))**2)/&
&                 (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&               Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                   Sin(theta)**2)/&
&                 (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))/&
&             ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&        ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3)) + &
&    (l*Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&            (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&             (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&          (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&            (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&            Sin(theta)**2*(a**2 + r**2 + &
&               (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))) - &
&       (2*a*M*r*Sin(theta)**2*Sqrt(((4*a**2*M**2*r**2*Sin(theta)**4)/&
&               (r**2 + a**2*Cos(theta)**2)**2 - &
&              (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&               (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&            (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&              (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&              Sin(theta)**2*(a**2 + r**2 + &
&                 (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))))/&
&        ((r**2 + a**2*Cos(theta)**2)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))))*&
&     ((2*Sqrt(2.)*a*M*(a**4 - 3*a**2*r**2 - 6*r**4 + a**2*(a**2 - r**2)*Cos(2*theta))*&
&          (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&          Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*(1/(Sin(theta)))**2)/&
&            (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&              r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)&
&              ))*Sin(theta)**2*(((1/(Sin(theta)))**2*&
&               (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&               (-2*r*Sin(theta)**2 + &
&                 (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                  (r**2 + a**2*Cos(theta)**2)**2))/&
&             ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                 r*(1/(Sin(theta)))**2*&
&                  (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&            (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                 a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&               (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&               Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                   (1/(Sin(theta)))**2)/&
&                 (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&               Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                   (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                 (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&             ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&            (M*(r**2 - a**2*Cos(theta)**2)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&             (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (l**2*(2*M - r)*r + &
&                 r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                 a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&        ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&       (Sqrt(2.)*M*(a**2 + r**2)*(a**2 - 2*r**2 + a**2*Cos(2*theta))*&
&          (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&            a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&          Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&            (l**2*(2*M - r)*r + &
&              r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&              a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*&
&          (((1/(Sin(theta)))**2*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                2*(-2*r*Sin(theta)**2 + &
&                 (2*a**2*M*(r**2 - a**2*Cos(theta)**2)*Sin(theta)**4)/&
&                  (r**2 + a**2*Cos(theta)**2)**2))/&
&             ((r**2 + a**2*Cos(theta)**2)*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                 r*(1/(Sin(theta)))**2*&
&                  (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) + &
&            (2*a*M*(-r**2 + a**2*Cos(theta)**2)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                 a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&               (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&               Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                   (1/(Sin(theta)))**2)/&
&                 (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&               Sin(theta)**2*Sqrt(((a**2 + r*(-2*M + r))*&
&                   (a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&                 (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))))/&
&             ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))**2) + &
&            (M*(r**2 - a**2*Cos(theta)**2)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**2)/&
&             (2.*(r**2 + a**2*Cos(theta)**2)**3*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (l**2*(2*M - r)*r + &
&                 r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                 a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))))/&
&        ((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) + &
&       (8*Sqrt(2.)*a**3*M*r*Cos(theta)*&
&          (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&          Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*(1/(Sin(theta)))**2)/&
&            (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&              r*(1/(Sin(theta)))**2*(-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)&
&              ))*Sin(theta)**3*((2*(1/(Tan(theta)))*&
&               (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**2*&
&               (-a**2 - r**2 - (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                 (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&             ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                 r*(1/(Sin(theta)))**2*&
&                  (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&            (a**2*M*r*Cos(theta)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&             ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (l**2*(2*M - r)*r + &
&                 r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                 a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&            (2*a*M*r*(a**2 + r**2)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                 a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&               (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&               Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                   (1/(Sin(theta)))**2)/&
&                 (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&               Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                   Sin(theta)**2)/&
&                 (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))/&
&             ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&        ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3) - &
&       (2*Sqrt(2.)*a**2*M*r*(a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&            a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&          Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*Sin(theta)**2)/&
&            (l**2*(2*M - r)*r + &
&              r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&              a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta)*&
&          ((2*(1/(Tan(theta)))*(2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)**&
&                2*(-a**2 - r**2 - &
&                 (4*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) - &
&                 (2*a**4*M*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2))/&
&             ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                 r*(1/(Sin(theta)))**2*&
&                  (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2))) - &
&            (a**2*M*r*Cos(theta)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                  a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))**2*Sin(theta)**3)/&
&             ((a**2 + r*(-2*M + r))*(r**2 + a**2*Cos(theta)**2)**3*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&               (l**2*(2*M - r)*r + &
&                 r*Sin(theta)**2*(-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                 a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2))) + &
&            (2*a*M*r*(a**2 + r**2)*&
&               (a**4 - 4*a*l*M*r + 2*r**4 + a**2*r*(2*M + 3*r) + &
&                 a**2*(a**2 + r*(-2*M + r))*Cos(2*theta))*&
&               (2*a*M*r + a**2*l*(1/(Tan(theta)))**2 + l*r*(-2*M + r)*(1/(Sin(theta)))**2)*&
&               Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                   (1/(Sin(theta)))**2)/&
&                 (2*a**2*M*r + a**2*(1/(Tan(theta)))**2*(a**2 + r**2 - l**2*(1/(Sin(theta)))**2) + &
&                   r*(1/(Sin(theta)))**2*&
&                    (-4*a*l*M + a**2*r + r**3 + l**2*(2*M - r)*(1/(Sin(theta)))**2)))*&
&               Sqrt(((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))*&
&                   Sin(theta)**2)/&
&                 (l**2*(2*M - r)*r + &
&                   r*Sin(theta)**2*&
&                    (-4*a*l*M + a**2*r + r**3 + 2*a**2*M*Sin(theta)**2) + &
&                   a**2*Cos(theta)**2*(-l**2 + (a**2 + r**2)*Sin(theta)**2)))*Sin(2*theta))/&
&             ((a**2 + r*(-2*M + r))**2*(r**2 + a**2*Cos(theta)**2)**3*&
&               (a**2 + 2*r**2 + a**2*Cos(2*theta))**2)))/&
&        ((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*Cos(2*theta))**3))))/&
&((a**2 - 2*M*r + r**2)*gamma*kappa)

end subroutine Gast_f




