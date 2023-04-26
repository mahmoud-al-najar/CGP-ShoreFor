using CSV
using MAT
using DataFrames
using GLPK
using JuMP
using Dates
using Dierckx
# using Statistics
using NaNStatistics

const MATLAB_EPOCH = Dates.DateTime(-1,12,31)
const MONTHLY_DATADIR = "/home/mn/Documents/data"
date2num(d::Dates.DateTime) = Dates.value(d-MATLAB_EPOCH)/(1000*60*60*24)
num2date(n::Number) =  MATLAB_EPOCH + Dates.Millisecond(round(Int64, n*1000*60*60*24))
const INDEX_t = 4

interp1(x,v,xq) = Spline1D(x, v; k=1)(xq)
# nanmean(x) = mean(filter(!isnan,x))
# nanmaximum(x) = maximum(filter(!isnan,x))
# nanminimum(x) = minimum(filter(!isnan,x))
# nanstd(x) = std(filter(!isnan,x))
# nansum(x) = std(filter(!isnan,x))

function remove_monthly_means(t, x)
    monthly_vals = Dict(1=>Array{Float64}(undef, 0), 2=>Array{Float64}(undef, 0), 3=>Array{Float64}(undef, 0), 
        4=>Array{Float64}(undef, 0), 5=>Array{Float64}(undef, 0), 6=>Array{Float64}(undef, 0), 
        7=>Array{Float64}(undef, 0), 8=>Array{Float64}(undef, 0), 9=>Array{Float64}(undef, 0), 
        10=>Array{Float64}(undef, 0), 11=>Array{Float64}(undef, 0), 12=>Array{Float64}(undef, 0)
    )
    months = Dates.month.(num2date.(t))
    
    for i in eachindex(months)
        m = months[i]
        push!(monthly_vals[m], x[i])
    end
    
    new_x = x .* 0
    for i in eachindex(months)
        m = months[i]
        new_x[i] = x[i] - nanmean(monthly_vals[m])
    end

    return new_x, monthly_vals
end

function minmax_normalize(x::Array{Float64})
    return (x .- nanminimum(x)) ./ (nanmaximum(x) - nanminimum(x))
end

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

function _load_monthly_data(ds; d50=0.25, Tw=15, h=8, remove_means=false, normalize=false)
    t_waves = convert(Vector{Float64}, ds["time"][:])
    Tp = convert(Vector{Float64}, ds["tp"][:])
    Hs = convert(Vector{Float64}, ds["hs"][:])
    Dir = convert(Vector{Float64}, ds["dir"][:])
    E = convert(Vector{Float64}, ds["ewave"][:])
    Sla = convert(Vector{Float64}, ds["sla"][:])
    rivdis = convert(Vector{Float64}, ds["rivdis"][:])
    
    t_shore = t_waves
    X = convert(Vector{Float64}, ds["X"][:])

    hWaveH = 8                               # water depth
    Ho = calcHoShoal(Hs, Tp, hWaveH)  #calculate deep water equivalent
    g = 9.81
    Hsb = 0.39.*g.^(1/5).*(Tp.*Ho.^2).^(2/5)  # eq.10

    w = fall_velocity(d50/1000, Tw)
    P = calcPb(Hsb) .^ 0.5
    omega = Hsb ./ (Tp * w)
    dt_waves = Int(round(nanmean(diff(t_waves))))
    t = collect(0:length(t_waves)-1) .* dt_waves

    if remove_means
        X, monthly_vals = remove_monthly_means(t_shore, X)
    end
    if normalize
        Tp = minmax_normalize(Tp)
        Hsb = minmax_normalize(Hsb)
        Dir = minmax_normalize(Dir)
        X = minmax_normalize(X)
        E = minmax_normalize(E)
        Sla = minmax_normalize(Sla)
        rivdis = minmax_normalize(rivdis)
        omega = minmax_normalize(omega)
        P = minmax_normalize(P)
    end

    data = [Tp, Hsb, Dir, t_waves, X, t_shore, E, Sla, rivdis, omega, P, t]
    # @assert count(isnan, data[i]) == 0
    return data
end

function load_monthly_dataset_GRANDPOPO(; calibration=false, remove_means=false, normalize=false)  # NO NAN's
    # BEST::
    phi=30.0
    # corr=0.81
    # rmse=3.18
    # ccc=0.79
    # mielke_skill=0.8

    ds = matread("$MONTHLY_DATADIR/Sites_X_sla_ewave_rivdis.mat")["GrandPopo"]
    data = _load_monthly_data(ds; remove_means=remove_means, normalize=normalize)
    data = [data..., phi]
    if !calibration
        return data
    else
        t = data[INDEX_t]
        dates = num2date.(t)
        calibration_date = DateTime(2015, 1, 1)
        i_end = findlast(x->x<=calibration_date, dates)[1]

        for i in eachindex(data)
            if typeof(data[i]) == Array{Float64,1}
                data[i] = data[i][1:i_end]
            end
        end
        return data
    end
end

function load_monthly_dataset_NARRABEEN(;calibration=false, remove_means=false, normalize=false)  # last index == NAN
    phi=78.0
    # corr=0.56
    # rmse=2.14
    # ccc=0.48
    # mielke_skill=0.48
    ds = matread("$MONTHLY_DATADIR/Sites_X_sla_ewave_rivdis.mat")["Narrabeen"]
    data = _load_monthly_data(ds; remove_means=remove_means, normalize=normalize)
    for i in eachindex(data)
        data[i] = data[i][1:end-1]
    end
    
    data = [data..., phi]
    
    if !calibration
        return data
    else
        t = data[INDEX_t]
        dates = num2date.(t)
        calibration_date = DateTime(2010, 1, 1)
        i_end = findlast(x->x<=calibration_date, dates)[1]

        for i in eachindex(data)
            if typeof(data[i]) == Array{Float64,1}
                data[i] = data[i][1:i_end]
            end
        end
        return data
    end

end

function load_monthly_dataset_DUCK(;calibration=false, remove_means=false, normalize=false)  # 70 first indices == NAN
    phi=351.0
    # corr=0.24
    # rmse=7.66
    # ccc=0.11
    # mielke_skill=0.12
    ds = matread("$MONTHLY_DATADIR/Sites_X_sla_ewave_rivdis.mat")["Duck"]
    data = _load_monthly_data(ds; remove_means=remove_means, normalize=normalize)
    for i in eachindex(data)
        data[i] = data[i][71:end]
    end
    data = [data..., phi]
    
    if !calibration
        return data
    else
        t = data[INDEX_t]
        dates = num2date.(t)
        calibration_date = DateTime(2004, 1, 1)
        i_end = findlast(x->x<=calibration_date, dates)[1]

        for i in eachindex(data)
            if typeof(data[i]) == Array{Float64,1}
                data[i] = data[i][1:i_end]
            end
        end
        return data
    end
end

function load_monthly_dataset_TORREYPINES(;calibration=false, remove_means=false, normalize=false)  # NO NAN's
    phi=34.0
    # corr=0.49
    # rmse=10.88
    # ccc=0.38
    # mielke_skill=0.39
    ds = matread("$MONTHLY_DATADIR/Sites_X_sla_ewave_rivdis.mat")["TorreyPines"]
    data = _load_monthly_data(ds; remove_means=remove_means, normalize=normalize)
    data = [data..., phi]
    if !calibration
        return data
    else
        t = data[INDEX_t]
        dates = num2date.(t)
        calibration_date = DateTime(2012, 1, 1)
        i_end = findlast(x->x<=calibration_date, dates)[1]

        for i in eachindex(data)
            if typeof(data[i]) == Array{Float64,1}
                data[i] = data[i][1:i_end]
            end
        end
        return data
    end
end

function load_monthly_dataset_TRUCVERT(;calibration=false, remove_means=false, normalize=false)  # NO NAN's
    phi=1000.0
    # corr=0.25
    # rmse=29.47
    # ccc=0.12
    # mielke_skill=0.13
    ds = matread("$MONTHLY_DATADIR/Sites_X_sla_ewave_rivdis.mat")["TrucVert"]
    data = _load_monthly_data(ds; remove_means=remove_means, normalize=normalize)
    data = [data..., phi]
    if !calibration
        return data
    else
        t = data[INDEX_t]
        dates = num2date.(t)
        calibration_date = DateTime(2011, 1, 1)
        i_end = findlast(x->x<=calibration_date, dates)[1]

        for i in eachindex(data)
            if typeof(data[i]) == Array{Float64,1}
                data[i] = data[i][1:i_end]
            end
        end
        return data
    end
end
