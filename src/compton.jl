export ComptonKleinNishina, sample
struct ComptonKleinNishina end

function cross_section_per_atom(::ComptonKleinNishina, Z, E)

end

function dσ_dθ(::ComptonKleinNishina, Z, E)

end

function differential_cross_section_per_atom(::ComptonKleinNishina,
                                         Z, E, energy_ratio)
end

"""
    compton_formula(E0::UF.Energy, theta)

Compute energy of scattered photon.

* `E0`: Energy of primary photon.
* `theta`: scattering angle
"""
function compton_formula(E0::UF.Energy, theta)
    mec2 = electron_rest_energy
    E0 * mec2 / (mec2 + E0 * (1-cos(theta)))
end

function rejection_function(eps, sintheta2)
    ret = 1 - (eps / (1 + eps^2))*sintheta2
end

function sample_buggy(rng::AbstractRNG, ::ComptonKleinNishina, E0::UF.Energy)
    # composition rejection with (alpha1*f1 + alpha2*f2) * g
    E_elec = electron_rest_energy
    eps0 = uconvert(NoUnits, E_elec / (E_elec + 2*E0)) # backward scattering energy_ratio

    alpha1 = log(1/eps0)
    alpha2 = (1-eps0^2) / 2

    while true
        r1 = rand(rng)
        r2 = rand(rng)
        if r1 < alpha1 / (alpha1 + alpha2)
            eps = exp(-r2*alpha1)
        else
            eps2 = eps0^2 + r2*(1 - eps0)^2
            eps = sqrt(eps2)
        end

        # t = 1 - cos(theta)
        t = uconvert(NoUnits, E_elec*(1-eps) / (E0*eps))
        sintheta2 = t*(2-t)

        if rejection_function(eps, sintheta2) > rand(rng)
            return (energy_ratio=eps, costheta = 1-t)
        end
    end
end

function sample(rng::AbstractRNG, ::ComptonKleinNishina, E0::UF.Energy)
    E0_m = E0/electron_rest_energy
    eps0       = 1/(1 + 2*E0_m)
    alpha1     = -log(eps0)
    alpha2     = alpha1 + (1- eps0^2)/2
    while true
        if alpha1 > alpha2*rand(rng)
            eps   = exp(-alpha1*rand(rng))
        else
            eps2  = eps0^2 + (1 - eps0^2)*rand(rng)
            eps   = sqrt(eps2)
        end
  
        one_costheta    = uconvert(NoUnits, (1 - eps)/(eps*E0_m))
        sintheta2   = one_costheta*(2-one_costheta)
        if rand(rng) < rejection_function(eps, sintheta2)
            costheta = 1-one_costheta
            ret = (energy_ratio=eps, costheta=costheta)
            return ret
        end
    end
end
