
subroutine drlogut_f(drlogut_func, M, a, r, theta, l)
    
    double precision, intent(in) :: M, a, r, theta, l
    double precision, intent(out) :: drlogut_func


    drlogut_func = ((l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&    (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&    Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))*&
&  (((-16*a**2*M**2*r**3*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**3 + &
&       (8*a**2*M**2*r*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&       (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&        (2*r - (4*a**2*M*r**2*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)**2 + &
&          (2*a**2*M*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)) - &
&       ((-4*M*r**2)/(r**2 + a**2*Cos(theta)**2)**2 + (2*M)/(r**2 + a**2*Cos(theta)**2))*&
&        Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))&
&      /(l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&       (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&       Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))&
&     - ((l**2*((-4*M*r**2)/(r**2 + a**2*Cos(theta)**2)**2 + &
&            (2*M)/(r**2 + a**2*Cos(theta)**2)) + &
&         (8*a*l*M*r**2*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)**2 - &
&         (4*a*l*M*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&         Sin(theta)**2*(2*r - (4*a**2*M*r**2*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)**2 + &
&            (2*a**2*M*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))*&
&       ((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&         (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&          (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))))/&
&     (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&        (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&        Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))&
&       **2))/&
&(2.*((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&    (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&     (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))))

end subroutine drlogut_f