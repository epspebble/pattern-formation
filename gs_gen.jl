const tf = parse(Float64, ARGS[1])
const F = parse(Float64, ARGS[2])
const k = parse(Float64, ARGS[3])
const Du, Dv = 0.2097, 0.1050
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
@time sol = solve(prob, CVODE_BDF(linear_solver=:GMRES))

filename="gs_F=$(F)_k=$(k)_tf=$(tf).nc"
varname = "u"

println("Saving $varname to $filename ...")
@time begin
    # Copy
    t = sol.t
    M = length(sol.t)
    u = zeros(M,N,N)
    for i in range(1,M)
        u[i,:,:] = sol.u[i][:,:,1]
    end
        
    attribs = Dict("data_min" => 0.0, "data_max" => 1.0)
    nccreate(filename,varname,"t", t, "x", dom, "y", dom, atts=attribs)
    ncwrite(u,filename,varname)

    # Get ready for garbage collection.
    sol = Nothing
    u = Nothing
    t = Nothing 
end