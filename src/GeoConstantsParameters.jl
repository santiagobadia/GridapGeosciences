const rₑ  = 6371220.0 # Radius of the Earth (meters)
const Ωₑ  = 7.292e-5  # Rotational frequency of the Earth (1/second)
const g   = 9.80616   # Gravitational constant (meters/second^2)
function f(x)         # Coriolis term
  2.0*Ωₑ*x[3]/rₑ
end
