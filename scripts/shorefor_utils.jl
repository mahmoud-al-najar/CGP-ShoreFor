function moving_sum(x::AbstractArray)
    if length(x) > 1
        x = deepcopy(x)
        for i in collect(1:length(x)-1)
            x[i] = x[i] + x[i+1]
        end
        x[1:end-1]
    else
        x
    end
end
function run_shorefor(omega, P, t)
    v = 0.0
    xs = []
    for i in collect(D+1:length(omega))
        input_omega = omega[Int(i-D)+1:Int(i)]
        input_P = P[Int(i-D)+1:Int(i)]
        current_time = t[Int(i)]

        indices = collect(1:D)
        flipped_indices = reverse(indices)
        negative_flipped_indices = -flipped_indices
        negative_flipped_indices_divby_phi = negative_flipped_indices ./ phi
        weights = 10 .^ negative_flipped_indices_divby_phi

        weighted_omega = input_omega .* weights

        sum_weights = sum(weights)
        sum_weighted_omega = sum(weighted_omega)
        omega_eq = sum_weighted_omega / sum_weights
        
        d_omega = omega_eq .- input_omega
        norm_d_omega = d_omega ./ std(d_omega)
        
        P_pow_half = input_P .^ 0.5
        F = P_pow_half .* norm_d_omega
        
        erosion_bools = omega_eq .<= input_omega
        accretion_bools = omega_eq .> input_omega
        Fe = F .* erosion_bools
        Fa = F .* accretion_bools
        
        s1 = r .* Fe
        s2 = s1 .+ Fa
        # v = last_val + 0.5 * (i - (i-1)) * (s2[end] + s2[end-1]) * dt_waves  # ORIGINAL (i - (i-1)) == 1
        # v = last_val + 0.5 * (s2[end] + s2[end-1]) * dt_waves  # SIMPLIFIED
        # SIMPLIFIED version, broken down:
        ms_s2 = moving_sum(s2)
        end_ms_s2 = ms_s2[end]
        half_end_ms_s2 = 0.5 * end_ms_s2
        dtw_half_end_ms_s2 = half_end_ms_s2 * dt_waves
        v = v + dtw_half_end_ms_s2

        # x = A[1] + (A[2] * t[Int(i)]) + (A[3] * v)  # ORIGINAL
        # broken down:
        sx1 = coeff_time * current_time
        sx2 = coeff_v * v
        sx3 = sx1 + sx2
        # x = intercept + sx3 
        x = sx3

        push!(xs, x)
    end
    xs
end


function intersect_lists(shoreline_dnums, wave_dnums)
    shoreline_dnums = vec(shoreline_dnums)
    wave_dnums = vec(wave_dnums)
    i_shorelines = findall(x -> x in wave_dnums, shoreline_dnums)
    i_waves = findall(x -> x in shoreline_dnums, wave_dnums)
    return i_shorelines, i_waves
end

function eval(t_pred, pred, t_target, target)
    i_pred, i_target = intersect_lists(t_pred, t_target)
    eval_pred = pred[i_pred]
    eval_target = target[i_target]
    corr_str = round(cor(eval_target, eval_pred), digits=2)
    rmse_str = round(sqrt(mean((eval_target - eval_pred).^2)), digits=2)

    return "corr=$corr_str \t rmse=$rmse_str"
end
