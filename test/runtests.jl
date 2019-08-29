using Rambo
using Test

const me = 0.510998928; # electron mass
const mμ = 105.6583715; # muon mass
const αem = 1.0 / 137.04; # EM fine structure
const qe = sqrt(4π * αem); # electric charge
const GF = 1.1663787e-11; # Fermi-constant
const mw = 80.385e3;  # W-mass MeV
const mz = 91.1876e3;  # Z-Mass MeV
const sw = sqrt(0.2223); # sine of weak-mixing angle
const cw = sqrt(1.0 - sw^2)  # cosine of weak-mixing angle

@testset "Rambo.jl" begin
    @test begin
        function msqrd_ee_to_mumu(momenta::Array{FourMomentum{Float64}, 1})
            p3::FourMomentum{Float64} = momenta[1]
            p4::FourMomentum{Float64} = momenta[2]

            P::FourMomentum{Float64} = sum(momenta)
            Q::Float64 = P.e
            pmag::Float64 = sqrt(Q^2/4 - me^2)

            p1::FourMomentum{Float64} = FourMomentum(Q/2, 0.0, 0.0, pmag)
            p2::FourMomentum{Float64} = FourMomentum(Q/2, 0.0, 0.0, -pmag)

            s::Float64 = scalar_product(p1+p2, p1+p2)
            t::Float64 = scalar_product(p1-p3, p1-p3)

            return ((32π^2 * (2me^4 + 2mμ^4 + s^2 + 4me^2 * (mμ^2 - t) -
                    4mμ^2*t + 2s * t + 2t^2) * αem^2)/s^2)
        end

        isp_masses = [me, me]
        fsp_masses = [mμ, mμ]
        cme = 1000.0

        rambo = cross_section(cme, isp_masses, fsp_masses; nevents=10000, msqrd=msqrd_ee_to_mumu)
        analytic = ((4π * sqrt(1 - 4mμ^2 / cme^2) * (2me^2 + cme^2) * (2mμ^2 + cme^2) * αem^2) /
                    (3 * sqrt(1 - 4me^2 / cme^2) * cme^6))
        # println((rambo[1] - analytic) / analytic)
        # println(rambo[2] / rambo[1] * 100)
        abs((rambo[1] - analytic) / analytic) < 1e-2
    end

    @test begin
        function msqrd_μ_to_eνν(momenta)
            pe = momenta[1]
            peν = momenta[2]
            pμν = momenta[3]
            pμ = sum(momenta)
            return 64GF^2 * scalar_product(pe, pμν) * scalar_product(pμ, peν)
        end

        fsp_masses = [me, 0.0, 0.0]
        rambo = decay_width(mμ, fsp_masses; nevents=100000, msqrd=msqrd_μ_to_eνν)
        r = me^2 / mμ^2
        corr_fac = 1.0 - 8r + 8r^3 - r^4 - 12r^2 * log(r)
        analytic = GF^2 * mμ^5 / 192π^3 * corr_fac
        # println((rambo[1] - analytic) / analytic)
        # println(rambo[2] / rambo[1] * 100)
        abs((rambo[1] - analytic) / analytic) < 1e-2
    end

    @test begin
        sw::Float64 = sqrt(0.2223)
        cw::Float64 = sqrt(1.0 - sw^2)
        fsp_masses = [me, me]

        function msqd_Z_to_ee(momenta):
            p1 = momenta[0]
            p2 = momenta[1]

            return (qe^2* (
                    2 * (1 - 4sw^2 + 8sw^4) * scalar_product(p1, p2)^2 +
                    2 * (1 - 4sw^2 + 8sw^4) * me^4 +
                    12 * sw^2 * (-1 + 2 * sw^2) * me^2 * mz^2 +
                    (1 - 4 * sw^2 + 8 * sw^4) * scalar_product(p1, p2) *
                    (4 * me^2 + mz^2))) /
                    (6.0 * cw^2 * sw^2 * mz^2)

            rambo = compute_decay_width(
                fsp_masses, mz, num_ps_pts=10, mat_elem_sqrd=msqd_Z_to_ee
            )

            num = qe^2 * (8.0sw^4 - 4sw^2 + 1) * mz
            den = 96π * cw^2 * sw^2
            analytic = num / den

            assert np.isclose(rambo[0], analytic, atol=0.0, rtol=5e-2)
    end

end
