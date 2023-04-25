# using DSP
using MAT
using Dates
using GLPK
if isdefined(Base, :Grisu)
    import Base.Grisu
else
    import Grisu
end
using JuMP
# using Plots
using Dierckx
# include("conv_test.jl")

const MATLAB_EPOCH = Dates.DateTime(-1,12,31)

date2num(d::Dates.DateTime) = Dates.value(d-MATLAB_EPOCH)/(1000*60*60*24)
num2date(n::Number) =  MATLAB_EPOCH + Dates.Millisecond(round(Int64, n*1000*60*60*24))

interp1(x,v,xq) = Spline1D(x, v; k=1)(xq)
# nanmean(x) = mean(filter(!isnan,x))
# nanstd(x) = std(filter(!isnan,x))

function fall_velocity(D, Tw)
    # D = Grain size [m]
    # Tw = Temp in degrees [C]
    # w returned in m/s
    D=D*100

    ROWs=2.75	# Density of sand (Mark 2.75, Kristen, 2.65?)
    g=981		# Gravity cm/s^2

    T   =[5, 10, 15, 20, 25]
    v   =[0.0157, 0.0135, 0.0119, 0.0105, 0.0095]
    ROW =[1.028, 1.027, 1.026, 1.025, 1.024] 

    vw=interp1(T,v,Tw)
    ROWw=interp1(T,ROW,Tw)

    A = ((ROWs-ROWw)*g*(D.^3))./(ROWw*(vw.^2))

    if A < 39
        w=((ROWs-ROWw)*g*(D.^2))./(18*ROWw*vw)
    else
        if A < 10^4   
            w=((((ROWs-ROWw)*g./ROWw).^0.7)*(D.^1.1))./(6*(vw.^0.4))
        else
            w=sqrt(((ROWs-ROWw)*g*D)./(0.91*ROWw))
        end
    end

    w=w./100 # convert to SI (m/s)
    return w
end

function calcCg(h,T)
    T = T[:, 1]
    g = 9.81;
    y=4.03*h./(T.^2);
    kd2 = y.^2 .+ y ./ (1 .+ (0.666 .* y)+(0.355 .* y.^2)+(0.161 .* y.^3)+(0.0632 .* y.^4)+(0.0218 .* y.^5)+(0.00564 .* y.^6));
    kh=sqrt.(kd2);
    Cg = g .* T ./ (2*pi) .* (tanh.(kh)) .* (0.5.*(1 .+ 2 .* kh ./ sinh.(2 .* kh)));  
    return Cg
end

function calcHoShoal(H,T,h1)
    #
    #  function Ho = calcHoShoal(H,T,h)
    # function to reverse shoal wave height to deep water equivalent.
    # 
    # kristen, 10
    #
    g = 9.81;
    Cgo = 1/4*g.*T./pi;
    Cg=calcCg(h1,T);
    Ks = sqrt.(Cgo./Cg);
    Ho = H./(Ks);   
    return Ho 
end

function calcPb(H)  # ,T)
    #function to calculate wave power at breaking
    #P = ECn, where E = wave Energy, Cn=Cg = group velocity

    g=9.81;
    rho = 1025;
    gamma=0.78;
    E = 1 ./ 16 .* rho .* g .* H.^2;
    Cg=sqrt.(g.*H./gamma);
    P=E.*Cg;
    return P
end

function WS85FilterConv_TB(omega, D, phi, dt)
    # code to calculate the omegaMean value based on Wright and Short 1985 paper.
    # inputs 
    # ~~~~~~
    # omega = time series of dimensionless fall velocity
    # dt = omega time step in seconds
    # D = number of days used in back filtering - set this to 2*phi
    # phi = number of days when beach memory is 10%

    # The method utilises convolution theorem to apply the filter. It is much faster than using loops! 
    # However, the methodology does not give good results for the last phi data points which must be 
    # calculated using the slow looping method
    # 
    # Mark davidson 11/7/12
    # 

    # println("Computing equilibrium omega value ... WS85 convolution")

    dt = dt ./ (3600 .* 24)
    D = round(D ./ dt)
    phi = round.(phi ./ dt)
    meanOmega = nanmean(omega)
    omega = omega .- meanOmega

    # Define filter back-to-front for convolution
    ii = collect(0:D-1)
    padding = zeros(trunc(Int, D-1), 1)
    filterCoeff = 10 .^ (-abs.(ii) ./ phi)
    filterCoeff = vcat(padding, filterCoeff)
    window = filterCoeff ./ sum(filterCoeff)
    # Nw = length(window)

    #### perform convolution
    k = size(window)[1]
    padding_size = trunc(Int, ceil((k - 1) / 2))
    input_size = size(omega)[1]
    omegaFiltered = conv(omega, window)
    omegaFiltered = omegaFiltered[padding_size + 1:padding_size+input_size]
    omegaFiltered = round.(omegaFiltered, digits=4)
    ####
    # omegaFiltered = mconv(omega,window)
    ####
    
    # Finally add on mean
    omegaFiltered = omegaFiltered .+ meanOmega
    return omegaFiltered
end

function calc_omega_eq(omega, D, phi, dt)
    dt = dt ./ (3600 .* 24)
    D = round(D ./ dt)
    phi = round.(phi ./ dt)
    meanOmega = nanmean(omega)
    omega = omega .- meanOmega

    y2mn = omega.*0
    for i in D+1:length(omega)
        if i % ceil(length(omega) / 10) == 0 || i==length(omega)
            println(i, "/", length(omega))
        end
        # weighted_omega = omega[Int(i-D)+1:Int(i)] .* 10.0 .^(-reverse(collect(1:D))./phi)
        weights = 10 .^(-reverse(collect(1:D))./phi)
        weighted_omega = omega[Int(i-D)+1:Int(i)] .* weights
        o = sum(weighted_omega)
        w = sum(weights)
        y2mn[Int(i)] = o/w
    end
    y2mn .+ meanOmega
end

function cumtrapz(X, Y)
    # Check matching vector length
    @assert length(X) == length(Y)
    # Initialize Output
    out = similar(X)
    out[1] = 0.0
    # Iterate over arrays
    for i in 2:length(X)
        out[i] = out[i-1] + 0.5 * (X[i] - X[i-1]) * (Y[i] + Y[i-1])
    end
    # Return output
    return out
end

function mcumtrapz(vx, a)
    vals = [0.0]
    for i in collect(2:length(vx))
        # global last_val
        v = vals[i-1] + 0.5 * (vx[i] - vx[i-1]) * (a[i] + a[i-1])
        push!(vals, v)
    end
    vals
end

function calc_r(Fa, Fe)
    vx = collect(1:length(Fa))
    num = NumericalIntegration.integrate(vx, Fa)
    denum = NumericalIntegration.integrate(vx, Fe)
    return abs(num / denum)
end

function mconv(v, u)
    kernel = copy(length(u) < length(v) ? u : v)
    vector = copy(length(u) > length(v) ? u : v)

    vmean = mean(vector)
    vector .-= vmean
    out = vector .* 0
    l = length(kernel)
    for i in l+1:length(vector)
        wv = vector[Int(i-D)+1:Int(i)] .* kernel # (kernel ./ sum(kernel))
        out[Int(i)] = sum(wv)
    end

    return out .+ vmean
end
