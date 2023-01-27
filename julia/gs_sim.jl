using DifferentialEquations, Sundials, NetCDF

const dx = 1/102 # grid spacing
const dt = 1. # timestep
const N = 256 # number of 
const dom = collect(0:1:N-1)*dx;

# Choice of diffusion constants is from 
# https://github.com/pmneila/jsexp/blob/master/grayscott/index.html

const Du, Dv = 0.2097, 0.1050

const tf = parse(Float32,ARGS[1])

print(tf)
println()
const F = parse(Float32,ARGS[2])
print(F)
println()
const k = parse(Float32,ARGS[3])
print(k)
println()

function ic_pearson(dom)
    """
    This is the choice of initial pattern from Pearson's 1993 paper. Here we assume a
    256 x 256 grid is used and he specified the middle box of 20 grid points length and 
    width be modified to "seed" patterns.
    """
    k = 0.01 # noisiness of the random
#    k=0
    
    N = length(dom)
    u = zeros(N, N, 2)
    u[:,:,1] = 1 .+ k*rand(N,N) # trivial steady state u = 1
    u[:,:,2] = k*rand(N,N) # trivial steady state v = 0
    u[119:138,119:138,1] = 1/2 .+ k*rand(20,20) # Pearson's initial square
    u[119:138,119:138,2] = 1/4 .+ k*rand(20,20) # Pearson's initial square

    return u
end

# From SciML tutorial, for convenience defining the matrix resulting from periodic BC
limit(a, N) = a == N+1 ? 1 : a == 0 ? N : a 

#= The name gs2d! with ! is a convention to remind we do in-place modifications, and not
returning anything. =#

function gs2d!(du, u, p, t)
    F, k, Du, Dv = p
    # @inbounds asks that index bound checks be skipped for speed.
    for I in CartesianIndices((N, N)) 
        i, j = Tuple(I)
        x, y = dom[I[1]], dom[I[2]]
        ip1, im1, jp1, jm1 = limit(i+1, N), limit(i-1, N), limit(j+1, N), limit(j-1, N)
        
        du[i,j,1] = Du*(u[im1,j,1] + u[ip1,j,1] + u[i,jp1,1] + u[i,jm1,1] - 4u[i,j,1]) -
                    u[i,j,1]*u[i,j,2]^2 + F*(1-u[i,j,1])
        du[i,j,2] = Dv*(u[im1,j,2] + u[ip1,j,2] + u[i,jp1,2] + u[i,jm1,2] - 4u[i,j,2]) +
                    u[i,j,1]*u[i,j,2]^2 - (F+k)*u[i,j,2]
    end
end

function simulate(tf, F, k)
    p = (F, k, Du, Dv)
    u0 = ic_pearson(dom)
    gs2d_prob = ODEProblem(gs2d!,u0,(0,tf),p)
    sol = solve(gs2d_prob,CVODE_BDF(linear_solver=:GMRES),saveat=tf/1800);
    return sol.t, sol.u
end

t, u = simulate(tf,F,k)

filename="gs_F=(F)_k=(k)_tf=(tf).nc"
varname = "u"
attribs = Dict("data_min" => 0.0, "data_max" => 1.0)
nccreate(filename,varname,"x", dom, "y", dom, "t", t, atts=attribs)
ncwrite(u[:,:,1],filename,varname)


# simulate(18000,0.014,0.054) # Moving Spots

# simulate(18000,0.037,0.06) # Default

# simulate(18000,0.03,0.062) # Solitons

# simulate(18000,0.025,0.06) # Pulsating Solitons

# simulate(18000,0.078,0.061) # Worms

# simulate(18000,0.029,0.057) # Mazes

# simulate(18000,0.039,0.058) # Holes

# simulate(18000,0.026,0.051) # Chaos

# simulate(18000,0.034,0.056) # Chaos and Holes

# simulate(18000,0.018,0.051) # Spots and Loops

# simulate(18000,0.014,0.045) # Waves

# simulate(180000,0.062,0.06093) # The U-Skate World
