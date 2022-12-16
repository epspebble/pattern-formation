"""
    Solution generator for the Gray-Scott pattern-forming system.

    Julia packages required: DifferentialEquations, Sundials, NetCDF

    Command line usage of this script:

    Syntax #1:

    >>> julia gs_gen.il [tf] [F] [k]

        where [tf] = the time value to which we evolve the Gray-Scott system to
              [F] = the "fill" rate, e.g. F = 0.0620
              [k] = the "kill" rate, e.g. k = 0.0630
    
    A NetCDF file will be saved to contain the simulation results.

    Example:

    >>> julia gs_gen.jl 10000 0.0620 0.0630

    generates the output in "gs_F=0.123_k=0.2_tf=10000.0.nc" located in the current folder.

    For a gallery of possible pattern with their F, k values, see http://www.mrob.com/pub/comp/xmorphia/.


    To specify also the diffusion constants, use the syntax:

    Syntax #2:

    >>> julia gs_gen.jl [tf] [F] [k] [Du] [Dv]

        where [Du] and [Dv] are the diffusivities for the U and V components.

    The built-in values are Du = 0.2 and Dv = 0.1 which come from this project, which is close
    to that found in the following project (Du = 0.2097, Dv = 0.1050):
    
    https://github.com/pmneila/jsexp/blob/master/grayscott/index.html

    However, Robert Munafo's project, as well as the Pearson's paper on Science both use 
    much smaller Du and Dv values: Du = 2e-5, Dv = 1e-5 according the scientific summary 
    in this paper: http://www.mrob.com/sci/papers/munafo2014-1218.pdf

"""

const tf = parse(Float64, ARGS[1])
const F = parse(Float64, ARGS[2])
const k = parse(Float64, ARGS[3])

if length(ARGS) > 3
    @assert length(ARGS) == 5
    const Du = parse(Float64, ARGS[4])
    const Dv = parse(Float64, ARGS[5])
    custom_D = true # only used to determine filename for output
else
    #const Du, Dv = 0.2097, 0.1050
    const Du, Dv = 0.2, 0.1
    custom_D = false
end

println("Setting up Gray-Scott system with F = $(F), k = $(k), Du = $(Du), Dv = $(Dv)...")

@time begin
    using DifferentialEquations, Sundials, NetCDF

    const dx = 1/102 # grid spacing
    const dt = 1. # timestep
    const N = 256 # number of 
    const dom = collect(0:1:N-1)*dx;

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

    limit(a, N) = a == N+1 ? 1 : a == 0 ? N : a
    function gs2d!(du, u, p, t)
        F, k, Du, Dv = p
        # @inbounds asks that index bound checks be skipped for speed.
        @inbounds for I in CartesianIndices((N, N)) 
            i, j = Tuple(I)
            x, y = dom[I[1]], dom[I[2]]
            ip1, im1, jp1, jm1 = limit(i+1, N), limit(i-1, N), limit(j+1, N), limit(j-1, N)
            
            du[i,j,1] = Du*(u[im1,j,1] + u[ip1,j,1] + u[i,jp1,1] + u[i,jm1,1] - 4u[i,j,1]) -
                        u[i,j,1]*u[i,j,2]^2 + F*(1-u[i,j,1])
            du[i,j,2] = Dv*(u[im1,j,2] + u[ip1,j,2] + u[i,jp1,2] + u[i,jm1,2] - 4u[i,j,2]) +
                        u[i,j,1]*u[i,j,2]^2 - (F+k)*u[i,j,2]
        end
    end
    u0 = ic_pearson(dom)
    tspan = (0.0,tf)
    p = (F, k, Du, Dv)
    prob = ODEProblem(gs2d!,u0,tspan,p)
end # Done setting up.


println("Solving for t in $tspan...")
@time sol = solve(prob, CVODE_BDF(linear_solver=:GMRES),saveat=0:100:tf)

if custom_D
    filename="gs_F=$(F)_k=$(k)_tf=$(tf)_Du=$(Du)_Dv=$(Dv).nc"
else
    filename="gs_F=$(F)_k=$(k)_tf=$(tf).nc"
end 

println("Saving data to $filename ...")
@time begin
    # Copy
    t = sol.t
    M = length(sol.t)
    u = zeros(M,N,N)
    #v = zeros(M,N,N)
    for i in range(1,M)
        u[i,:,:] = sol.u[i][:,:,1]
	#v[i,:,:] = sol.u[i][:,:,2]
    end
        
    attribs = Dict("data_min" => 0.0, "data_max" => 1.0)
    nccreate(filename,"u","t", t, "x", dom, "y", dom, atts=attribs)
    #nccreate(filename,"v","t", t, "x", dom, "y", dom, atts=attribs)    
    ncwrite(t,filename,"t")
    ncwrite(u,filename,"u") # could annotate with PDE parameters here.
    #ncwrite(v,filename,"v") # could annotate with PDE parameters here.

    # Get ready for garbage collection.
    sol = Nothing
    u = Nothing
    t = Nothing 
end