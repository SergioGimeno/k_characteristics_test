
subroutine dthetalogut_f(dthetalogut_func, M, a, r, theta, l)
    
    double precision, intent(in) :: M, a, r, theta, l
    double precision, intent(out) :: dthetalogut_func


    dthetalogut_func = ((l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&    (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&    Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))*&
&  (-((((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&           (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&            (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))*&
&         ((4*a**2*l**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2)**2 - &
&           (8*a*l*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) - &
&           (8*a**3*l*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2 + &
&           2*Cos(theta)*Sin(theta)*(a**2 + r**2 + &
&              (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)) + &
&           Sin(theta)**2*((4*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&              (4*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2)))/&
&       (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&          (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&          Sin(theta)**2*(a**2 + r**2 + &
&             (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))**2) + &
&    ((16*a**2*M**2*r**2*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2 + &
&       (16*a**4*M**2*r**2*Cos(theta)*Sin(theta)**5)/(r**2 + a**2*Cos(theta)**2)**3 - &
&       2*Cos(theta)*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)*&
&        (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)) - &
&       (4*a**2*M*r*Cos(theta)*Sin(theta)**3*&
&          (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2)))/&
&        (r**2 + a**2*Cos(theta)**2)**2 - &
&       (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&        ((4*a**2*M*r*Cos(theta)*Sin(theta))/(r**2 + a**2*Cos(theta)**2) + &
&          (4*a**4*M*r*Cos(theta)*Sin(theta)**3)/(r**2 + a**2*Cos(theta)**2)**2))/&
&     (l**2*(-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2)) - &
&       (4*a*l*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2) + &
&       Sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))))&
&  )/(2.*((4*a**2*M**2*r**2*Sin(theta)**4)/(r**2 + a**2*Cos(theta)**2)**2 - &
&    (-1 + (2*M*r)/(r**2 + a**2*Cos(theta)**2))*Sin(theta)**2*&
&     (a**2 + r**2 + (2*a**2*M*r*Sin(theta)**2)/(r**2 + a**2*Cos(theta)**2))))

end subroutine dthetalogut_f