# PIMC.jl - Path Integral Monte Carlo 

"""
Perform path integral Monte Carlo simulation to estimate the 
current-current correlation function for a Holstein polaron system.

# Arguments
- `G::Real=0.5`: Holstein Electron-phonon coupling strength
- `ω0::Real=1`:  Phonon frequency 
- `T::Real=2`:   Temperature
- `m::Int=10`:   Number of imaginary time slices
- `Q::Int=5`:    Number of real time slices

# Returns
- `C_jj::Float64`: Estimate of current-current correlation function
"""
function PIMC(;G::Real=0.5, ω0::Real=1, T::Real=2, m::Int=10, Q::Int=5)
    β = 1/T    

    C_jj=0

    return C_jj
end

