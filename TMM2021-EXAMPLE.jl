#************************************************************************#
#                           TMM2021-EXAMPLE.jl                           #
#************************************************************************#
#                                                                        #
# Inverse dynamics simulation for a family of planar mechanisms via MMM. #
#                                                                        #
# From: ORSINO, R.M.M; HESS-COELHO, T.A.; MALVEZZI, F. (2021) Exporing   #
#    the common modularity among parallel manipulators through the       #
#    Modular Modeling Methodology (MMM). Submitted to: XIII Int. Conf.   #
#    on the Theory of Machines and Mechanisms.                           #
#                                                                        #
# Version: May, 5th 2021                                                 #
#                                                                        #
# Author: Dr. Renato Maia Matarazzo Orsino                               #
#                                                                        #
#************************************************************************#

using LinearAlgebra
using Plots


## -- TOPOLOGY -- ##

n_chains = 3
scenario = 1

if scenario == 1
    PRR_chains = [1, 2, 3]
    RPR_chains = []
    label_ = "3PRR"
elseif scenario == 2
    PRR_chains = [2, 3]
    RPR_chains = [1]
    label_ = "2PRR+RPR"
end


## -- PARAMETERS -- ##

# - acceleration of gravity [m/s^2] - #
g = 9.8

# - lengths [m] - #
l = [0.25, 0.25, 0.25]
b = [0.5, 0.4, 0.5]
l_E = 0.5
σ = [-1, 0, 1]
h = [0.0, 0.5, 1.0]

# - masses [kg] - #
m_A = 4.0 * l
m_B = 1.0 * b
m_E = 2.0

# - moments of inertia [kg*m^2] - #
I_A = @. m_A * l^2/3
I_B = @. m_B * b^2/3
I_E = m_E * l_E^2/3


## -- KINEMATICS -- ##

# - initialization - #

global nf = 3
global nx = nf
global nξ = Dict(
    1 => 2,
    2 => 2,
    3 => 2
    )
global x = zeros(nx)             # [x, y, ϕ]
global ξ = Dict(
    1 => zeros(nξ[1]),           # [q1, θ1]
    2 => zeros(nξ[2]),           # [q2, θ2]
    3 => zeros(nξ[3])            # [q3, θ3]
    )
global vx = zeros(nx)
global vξ = Dict(
    1 => zeros(nξ[1]),
    2 => zeros(nξ[2]),
    3 => zeros(nξ[3])
    )
global ax = zeros(nx)
global aξ = Dict(
    1 => zeros(nξ[1]),
    2 => zeros(nξ[2]),
    3 => zeros(nξ[3])
    )
global Jξ = Dict(
    1 => zeros(nξ[1], nξ[1]),
    2 => zeros(nξ[2], nξ[2]),
    3 => zeros(nξ[3], nξ[3])
    )
global Jx = Dict(
    1 => zeros(nξ[1], nx),
    2 => zeros(nξ[2], nx),
    3 => zeros(nξ[3], nx)
    )
global C = Dict(
    1 => zeros(nξ[1], nx),
    2 => zeros(nξ[2], nx),
    3 => zeros(nξ[3], nx)
    )
global d = Dict(
    1 => zeros(nξ[1]),
    2 => zeros(nξ[2]),
    3 => zeros(nξ[3])
    )
global V = Dict(
    (1,1) => zeros(2, nξ[1]),
    (1,2) => zeros(2, nξ[1]),
    (2,1) => zeros(2, nξ[2]),
    (2,2) => zeros(2, nξ[2]),
    (3,1) => zeros(2, nξ[3]),
    (3,2) => zeros(2, nξ[3]),
    "E" => zeros(2, nx)
    )
global Ω = Dict(
    (1,1) => zeros(1, nξ[1]),
    (1,2) => zeros(1, nξ[1]),
    (2,1) => zeros(1, nξ[2]),
    (2,2) => zeros(1, nξ[2]),
    (3,1) => zeros(1, nξ[3]),
    (3,2) => zeros(1, nξ[3]),
    "E" => zeros(1, nx)
    )
global e = Dict(
    (1,1) => zeros(2),
    (1,2) => zeros(2),
    (2,1) => zeros(2),
    (2,2) => zeros(2),
    (3,1) => zeros(2),
    (3,2) => zeros(2),
    "E" => zeros(2),
    )
global ϵ = Dict(
    (1,1) => zeros(1),
    (1,2) => zeros(1),
    (2,1) => zeros(1),
    (2,2) => zeros(1),
    (3,1) => zeros(1),
    (3,2) => zeros(1),
    "E" => zeros(1),
    )


# - inverse kinematics - #

function inv_kinematics!()

    global x, ξ, vx, vξ, ax, aξ, Jx, Jξ, C, d, V, Ω, e, ϵ

    # - position - #

    for i in PRR_chains

        rq0 = x[1] + σ[i] * l_E * cos(x[3])
        rΔ = sqrt((2 * b[i]  - x[2] - σ[i] * l_E * sin(x[3])) * (2 * b[i] + x[2] + σ[i] * l_E * sin(x[3])))
        rq1 = rq0 - rΔ
        rq2 = rq0 + rΔ

        if abs(ξ[i][1] - rq1) < abs(ξ[i][1] - rq2)
            q = rq1
        else
            q = rq2
        end

        θ = atan(x[2] + σ[i] * l_E * sin(x[3]), x[1] + σ[i] * l_E * cos(x[3]) - q)

        Jξ[i] = [1 (-2 * b[i] * sin(θ)); 0 (2 * b[i] * cos(θ))]
        Jx[i] = [1 0 (-σ[i] * l_E * sin(x[3])); 0 1 (σ[i] * l_E * cos(x[3]))]

        C[i] = Jξ[i] \ Jx[i]
        vq, vθ = C[i] * vx

        d[i] = [(2 * b[i] * vθ^2 - σ[i] * l_E * vx[3]^2 * cos(θ - x[3])); (vθ^2 * sin(θ) - σ[i] * l_E * vx[3]^2 * sin(x[3])/(2 * b[i]))] / cos(θ)
        aq, aθ = C[i] * ax + d[i]

        V[i,1] = [1.0 0; 0 0]
        V[i,2] = [cos(θ) 0; -sin(θ) b[i]]
        Ω[i,2] = [0 1]

        e[i,2] = [-b[i] * vθ^2; 0]

        ξ[i][1] = q
        ξ[i][2] = θ
        vξ[i][1] = vq
        vξ[i][2] = vθ
        aξ[i][1] = aq
        aξ[i][2] = aθ

    end

    for i in RPR_chains

        rq0 = - (l[i] + 2 * b[i])
        rΔ = sqrt((x[1] + σ[i] * l_E * cos(x[3]) - h[i])^2 + (x[2] + σ[i] * l_E * sin(x[3]))^2)
        rq1 = rq0 - rΔ
        rq2 = rq0 + rΔ

        if abs(ξ[i][1] - rq1) < abs(ξ[i][1] - rq2)
            q = rq1
        else
            q = rq2
        end

        θ = atan(x[2] + σ[i] * l_E * sin(x[3]), x[1] + σ[i] * l_E * cos(x[3]) - h[i])

        lP = q + l[i] + 2 * b[i]
        Jξ[i] = [cos(θ) (- lP * sin(θ)); sin(θ) (lP * cos(θ))]
        Jx[i] = [1 0 (-σ[i] * l_E * sin(x[3])); 0 1 (σ[i] * l_E * cos(x[3]))]

        C[i] = Jξ[i] \ Jx[i]
        vq, vθ = C[i] * vx

        d[i] = [(lP * vθ^2 - σ[i] * l_E * vx[3]^2 * cos(θ - x[3])); (-2 * vq * vθ + σ[i] * l_E * vx[3]^2 * sin(θ - x[3]))/lP]
        aq, aθ = C[i] * ax + d[i]

        lG = q + l[i] + b[i]
        V[i,2] = [1 0; 0 lG]
        Ω[i,1] = [0 1]
        Ω[i,2] = [0 1]

        e[i,2] = [-lG * vθ^2; 2 * vq * vθ]

        ξ[i][1] = q
        ξ[i][2] = θ
        vξ[i][1] = vq
        vξ[i][2] = vθ
        aξ[i][1] = aq
        aξ[i][2] = aθ

    end

    V["E"] = [1 0 0; 0 1 0]
    Ω["E"] = [0 0 1]

end


## -- EQUATIONS OF MOTION -- ##

global Mξ = Dict(
    1 => zeros(nξ[1], nξ[1]),
    2 => zeros(nξ[2], nξ[2]),
    3 => zeros(nξ[3], nξ[3]),
    )
global hξ = Dict(
    1 => zeros(nξ[1]),
    2 => zeros(nξ[2]),
    3 => zeros(nξ[3]),
    )
global gξ = Dict(
    1 => zeros(nξ[1]),
    2 => zeros(nξ[2]),
    3 => zeros(nξ[3]),
    )
global Mx = zeros(nx, nx)
global hx = zeros(nx)
global gx = zeros(nx)


function inv_dynamics!(t)

    global x, ξ, vx, vξ, ax, aξ, Jx, Jξ, C, d, V, Ω, e, ϵ
    global Mξ, hξ, gξ, τξ, Mx, hx, gx, τx

    x = r(t)
    vx = v(t)
    ax = a(t)

    inv_kinematics!()

    for i in PRR_chains
        Mξ[i] = m_A[i] * V[i,1]' * V[i,1] + m_B[i] * V[i,2]' * V[i,2] + Ω[i,2]' * I_B[i] * Ω[i,2]
        hξ[i] = m_A[i] * V[i,1]' * e[i,1] + m_B[i] * V[i,2]' * e[i,2] + Ω[i,2]' * I_B[i] * ϵ[i,2]
        gξ[i] = - m_A[i] * V[i,1]' * g * [0; -1] - m_B[i] * V[i,2]' * g * [-sin(ξ[i][2]); -cos(ξ[i][2])]
    end

    for i in RPR_chains
        Mξ[i] = Ω[i,1]' * I_A[i] * Ω[i,1] + m_B[i] * V[i,2]' * V[i,2] + Ω[i,2]' * I_B[i] * Ω[i,2]
        hξ[i] = Ω[i,1]' * I_A[i] * ϵ[i,1] + m_B[i] * V[i,2]' * e[i,2] + Ω[i,2]' * I_B[i] * ϵ[i,2]
        gξ[i] = - (m_A[i] * V[i,1]' + m_B[i] * V[i,2]') * g * [-sin(ξ[i][2]), -cos(ξ[i][2])]
    end

    Mx = (m_E * V["E"]' * V["E"] + Ω["E"]' * I_E * Ω["E"]) + sum(C[i]' * Mξ[i] * C[i] for i in 1:n_chains)
    hx = (m_E * V["E"]' * e["E"] + Ω["E"]' * I_E * ϵ["E"]) + sum(C[i]' * (Mξ[i] * d[i] + hξ[i]) for i in 1:n_chains)
    gx = - m_E * V["E"]' * g * [0; -1] + sum(C[i]' * gξ[i] for i in 1:n_chains)

    Cx = hcat(C[1][1,:], C[2][1,:], C[3][1,:])

    Cx \ (Mx * ax + hx + gx)
end


## -- PRESCRIBED MOTION -- ##

x0 = 2.0
y0 = 0.3
Ax = 0.1
ω = 4 * π
ϕ0 = 0
Aϕ = π/9
ϖ = 8 * π
r(t) = [x0 + Ax * cos(ω * t); y0 + Ax * sin(ω * t); ϕ0 + Aϕ * sin(ϖ * t)]
v(t) = [- Ax * ω * sin(ω * t); Ax * ω * cos(ω * t); Aϕ * ϖ * cos(ϖ * t)]
a(t) = [- Ax * ω * ω * cos(ω * t); - Ax * ω * ω * sin(ω * t); - Aϕ * ϖ * ϖ * cos(ϖ * t)]

x = r(0)
vx = v(0)
ax = a(0)
ξ[1] = [x0/4; π/3]
ξ[2] = [x0/4; π/2]
ξ[3] = [x0; 2 * π/3]

inv_kinematics!()
inv_dynamics!(0)


## -- INVERSE DYNAMICS SIMULATION -- ##

n_steps = 1001
τ = zeros(n_steps, nf)
q_ = zeros(n_steps, nf)
x_ = Dict(
    1 => zeros(n_steps, 2),
    2 => zeros(n_steps, 2),
    3 => zeros(n_steps, 2)
    )
t = LinRange(0.0, 1.0, n_steps)

for k in 1:n_steps
    τ[k,:] = inv_dynamics!(t[k])
    q_[k,:] = [ξ[i][1] for i in 1:n_chains]
    for i in 1:n_chains
        x_[i][k,:] = [x[1] + σ[i] * l_E * cos(x[3]), x[2] + σ[i] * l_E * sin(x[3])]
    end
end


## -- GRAFICOS -- ##

begin
    Px = plot(x_[1][:,1], x_[1][:,2], label="P_1", color=:green)
    plot!(x_[2][:,1], x_[2][:,2], label="P_2", color=:red)
    plot!(x_[3][:,1], x_[3][:,2], label="P_3", color=:blue)
    plot!(
        xlabel="x [m]",
        ylabel="y [m]",
        xlims=[1.35, 2.65],
        ylims=[0.0, 0.55],
        aspect_ratio=1,
        size=(350,167),
        foreground_color_legend = nothing,
        background_color_legend = RGBA(0.64,0.64,0.55,0.4)
    )
end

Pτ = plot(
    t, τ,
    title=label_*" mechanism",
    label=["τ_1" "τ_2" "τ_3"],
    color=[:green :red :blue],
    xlabel="time [s]",
    ylabel="actuator force [N]",
    xlims=[0.0, 1.0],
    ylims=[-1000, 1500],
    size=(250,250),
    foreground_color_legend = nothing,
    background_color_legend = RGBA(0.64,0.64,0.55,0.4)
    )
