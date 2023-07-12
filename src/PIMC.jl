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
- `N::Int=10`:   Number of (electron) sites

# Returns
- `C_jj::Float64`: Estimate of current-current correlation function
"""
function PIMC(;G::Real=0.5, ω0::Real=1, T::Real=2, m::Int=10, Q::Int=5, N::Int=10)
    β = 1/T    

    r= rand(1:N, m+Q+m+1) 
    # random electron positions for Im,Re,Im time slices 

    w=real_time_weight(r,m,Q, N=N)
    println("Real time weight: $(w)")

    C_jj=0
    return C_jj
end

"""
    real_time_weight(r, m, Q; kwargs...)

Compute the weight function w(r) for the real-time path integral.

Equation (26) in Miladić and Vukmirović 2023.

# Arguments
- `r::Vector{Int}`: Electron position state labels
- `m::Int`: Number of imaginary time slices
- `Q::Int`: Number of real time slices
- `J::Real=1`: Hopping integral
- `β::Real=1`: Inverse temperature
- `Δt::Real=0.1`: Time step

# Returns
- `w::Float64`: Weight function value w(r)
"""
function real_time_weight(r, m, Q; J=1, N=10, β=1, Δt=0.1)

    τ = β/m
    t = Q*Δt

    w = 1
    # Imaginary time slices before?
    @inbounds for j = 1:m-1
        w *= real_time_I(τ, r[j+1]-r[j]; J=J)
    end

    # I'm really confused about the indices in in the second and third product in (26)
    # Real time slices?
    @inbounds for p = m:m+Q-1
        w *= abs(real_time_I(-im*Δt, r[p+1]-r[p]; J=J))
    end

    # Imaginary time slices afterwards?! El-Ph coupling?
    @inbounds for q = m+Q+1:m+2Q
        w *= abs(real_time_I(im*Δt, r[q+1]-r[q]; J=J))
    end

    return w
end

"""
    real_time_I(z, Δr; J=1, C=1000)

Compute the Fourier transformed electron propagator I(z; Δr)
for the real-time path integral.

Equation (21) in Miladić and Vukmirović 2023.
(Naive implementation.)

# Arguments
- `z::Complex`: Complex time
- `Δr::Int`: Distance between electron positions
- `J::Real=1`: Hopping integral
- `C::Int=1000`: Number of Fourier components

# Returns
- `I::Complex`: Value of the propagator I(z; Δr)
"""
function real_time_I(z, Δr; J=1, N=1000)

    I = zero(ComplexF64)
    for n = 0:N-1
        I += cis(2π/N * n * Δr) *
             exp(2z*J*cos(2π/N*n))
    end
    I /= N

    return I
end



