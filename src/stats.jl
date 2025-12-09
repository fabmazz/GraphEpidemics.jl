import StatsBase: counts

export moving_average, compute_R_eff

function moving_average(data, window::Integer)
    n = length(data)
    result = similar(data)
    
    for i in 1:n
        start_idx = max(1, i - window)
        end_idx = min(n, i + window)
        result[i] = mean(data[start_idx:end_idx])
    end
    
    return result
end


function compute_R_eff(simdata, Tmax, nseeds;method::Symbol = :daily, smoothing::Integer=2, window::Integer=7)

    raw_R = fill(NaN, Tmax+1)

    if method == :daily
        iact = get_idx_infect(simdata)
        
        if length(iact) == 0
            0 ## no infecteds
        else
            getTR(ii) = get_t_R(simdata, ii)
            getTI(ii) = get_t_I(simdata, ii)

            date_rec =@.( convert(Int32,getTR(iact) +1) )
            date_I = @.( convert(Int32,getTI(iact) +1) )
            mean_infect_time = mean(date_rec) - mean(date_I)

            date_infected = get_t_infected(simdata, iact)
            
            # this is an array now
            new_infections = StatsBase.counts(@.(convert(Int32,date_infected +1)), Tmax+1)
            
            ### adjusted the arrays above, this is an array starting from 0
            N_I = cumsum(StatsBase.counts(date_I, Tmax+1)) .- cumsum(StatsBase.counts(date_rec, Tmax+1))
            raw_R[:] = @.(new_infections /N_I )* mean_infect_time
            raw_R[@. N_I == 0] .= 0
            #raw_R[@.(isnan(raw_R) || isinf(raw_R))] .= 0.

        end

        #return raw_R#moving_average(raw_R, window)

    elseif (method==:infectious) || (method==:outcome)
        sources = zeros(Tmax+1)
        targets = zeros(Tmax+1)

        iact = get_idx_infect(simdata)
        if length(iact) == 0
            return raw_R
        end

        ### Find the people becoming infectious or recovered in each time step
        tims =convert.(Int32, if (method==:infectious)
            get_t_I(simdata, iact)
        else
            get_t_R(simdata, iact)
        end
        )
        allts = StatsBase.counts(tims)
        mm = min(length(allts), length(sources))
        sources[1:mm] .= allts[1:mm]

        t_infected = convert.(Int32,get_t_infected(simdata, iact))
        targets[2:end] = StatsBase.counts(t_infected, Tmax) ## this automatically searches for k=1:Tmax
        
        ma = sources .> 0
        fill!(raw_R, 0.0)
        raw_R[ma] .= (targets./sources)[ma]

        num = cumsum(raw_R .* sources)
        num[window+1:end] = num[window+1:end] .- num[1:end-window]

        den = cumsum(sources)
        den[window+1:end] = den[window+1:end] .- den[1:end-window]

        values = fill(NaN, Tmax+1)
        mask = den.>0
        values[mask] = num[mask] ./ den[mask]
        return values 
    else
        throw(error("Method must be :daily, :infectious or :outcome"))
    end

    return raw_R
end