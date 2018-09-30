using CrossSections
using Test

using Unitful: MeV
using Random
rng = MersenneTwister()
proc = ComptonKleinNishina()

for E in [0.1MeV, 1MeV, 10MeV, 100MeV]
    for _ in 1:10
        ret = sample(rng, proc, E)
        @show E
        @show ret
    end
end
