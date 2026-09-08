# Independent branch-limit validation (issue #386)
#
# Deliberately independent of the OPF constraint/stamping code. It consumes only
# solved terminal phasors and reported branch currents and recomputes endpoint
# limits from BMOPF network data.

struct IndependentBranchFinding
    severity::Symbol
    code::String
    line::String
    endpoint::Symbol
    terminal::String
    message::String
    detail::Dict{String,Any}
end

@inline function _ib_wrap_pi(x::Real)
    y = mod(Float64(x) + π, 2π) - π
    return y == -π ? π : y
end

@inline function _ib_rating(v, k::Int)
    v isa Number && return Float64(v)
    v isa AbstractVector && k <= length(v) && return Float64(v[k])
    return nothing
end

@inline function _ib_finite_complex(vr, vi)
    vr isa Real && vi isa Real && isfinite(vr) && isfinite(vi)
end

function _ib_endpoint_current(vals, side::Symbol)
    vals isa AbstractDict || return nothing
    suffix = side === :from ? "fr" : "to"
    x = get(vals, "cm_$suffix", nothing)
    x isa Real && isfinite(x) && return abs(Float64(x))
    rk = side === :from ? "cr_fr" : "cr_to"
    ik = side === :from ? "ci_fr" : "ci_to"
    cr = get(vals, rk, nothing); ci = get(vals, ik, nothing)
    if _ib_finite_complex(cr, ci)
        return hypot(Float64(cr), Float64(ci))
    end
    cur = get(vals, side === :from ? "current_from" : "current_to", nothing)
    if cur isa AbstractDict
        cr = get(cur, "cr", get(cur, "real", nothing))
        ci = get(cur, "ci", get(cur, "imag", nothing))
        _ib_finite_complex(cr, ci) && return hypot(Float64(cr), Float64(ci))
    end
    nothing
end

function _ib_voltage(bus_results, bus::String, terminal::String)
    vals = get(bus_results, bus, nothing)
    vals isa AbstractDict || return nothing
    tv = get(vals, terminal, nothing)
    tv isa AbstractDict || return nothing
    vr = get(tv, "vr", nothing); vi = get(tv, "vi", nothing)
    _ib_finite_complex(vr, vi) || return nothing
    return complex(Float64(vr), Float64(vi))
end

"""
    independent_branch_limit_check(net, result) -> Vector{IndependentBranchFinding}

Independently verify both terminals of each line carrying `i_max` or `s_max`,
and signed branch-angle windows. No OPF model, constraint builder, or line
stamping routine is called.

Receiving-end apparent power is recomputed as `|V_to|*|I_to|`, so a receiving-
end violation cannot be hidden by a sending-end proxy. Reverse flow is naturally
supported because thermal limits use current magnitude. A zero-voltage endpoint
makes an angle undefined and is reported as informational rather than as a
physical angle violation.
"""
function independent_branch_limit_check(net::AbstractDict, result::AbstractDict)
    findings = IndependentBranchFinding[]
    buses = get(result, "bus", Dict())
    linecodes = get(net, "linecode", Dict())
    line_results = get(result, "line", Dict())

    for (lid_raw, line) in get(net, "line", Dict())
        line isa AbstractDict || continue
        lid = string(lid_raw)
        lr = get(line_results, lid_raw, get(line_results, lid, nothing))
        lr isa AbstractDict || continue
        lcid = get(line, "linecode", nothing)
        lc = lcid isa AbstractString ? get(linecodes, lcid, Dict()) : Dict()
        tm_fr = string.(get(line, "terminal_map_from", String[]))
        tm_to = string.(get(line, "terminal_map_to", String[]))
        length(tm_fr) == length(tm_to) || continue
        bf = string(get(line, "bus_from", "")); bt = string(get(line, "bus_to", ""))
        (isempty(bf) || isempty(bt)) && continue

        i_fr = get(line, "i_max", get(lc, "i_max", nothing))
        s_fr = get(line, "s_max", get(lc, "s_max", nothing))
        i_to = get(line, "i_max_to", get(lc, "i_max_to", i_fr))
        s_to = get(line, "s_max_to", get(lc, "s_max_to", s_fr))
        amin = get(line, "va_diff_min", get(lc, "va_diff_min", nothing))
        amax = get(line, "va_diff_max", get(lc, "va_diff_max", nothing))

        for k in eachindex(tm_fr)
            tf, tt = tm_fr[k], tm_to[k]
            row = get(lr, tf, nothing)
            row isa AbstractDict || (row = get(lr, tt, nothing))
            row isa AbstractDict || continue

            for (endpoint, term, bus, irating, srating) in
                    ((:from, tf, bf, i_fr, s_fr), (:to, tt, bt, i_to, s_to))
                i_lim = _ib_rating(irating, k)
                s_lim = _ib_rating(srating, k)
                i_mag = _ib_endpoint_current(row, endpoint)
                if i_mag !== nothing && i_lim !== nothing && i_mag > i_lim * (1 + 1e-6)
                    code = endpoint === :to ? "E.SOL.THERMAL_RECEIVING_CURRENT" :
                                               "E.SOL.THERMAL_SENDING_CURRENT"
                    push!(findings, IndependentBranchFinding(
                        :error, code, lid, endpoint, term,
                        "Line '$lid' $endpoint-end current exceeds i_max.",
                        Dict("line"=>lid, "endpoint"=>String(endpoint),
                             "terminal"=>term, "current_A"=>i_mag,
                             "i_max_A"=>i_lim)))
                end

                if s_lim !== nothing && i_mag !== nothing
                    v = _ib_voltage(buses, bus, term)
                    if v !== nothing
                        smag = abs(v) * i_mag
                        if smag > s_lim * (1 + 1e-6)
                            code = endpoint === :to ? "E.SOL.THERMAL_RECEIVING_APPARENT_POWER" :
                                                       "E.SOL.THERMAL_SENDING_APPARENT_POWER"
                            push!(findings, IndependentBranchFinding(
                                :error, code, lid, endpoint, term,
                                "Line '$lid' $endpoint-end apparent power exceeds s_max.",
                                Dict("line"=>lid, "endpoint"=>String(endpoint),
                                     "terminal"=>term, "voltage_V"=>abs(v),
                                     "current_A"=>i_mag, "s_VA"=>smag,
                                     "s_max_VA"=>s_lim)))
                        end
                    end
                end
            end

            if amin !== nothing || amax !== nothing
                vf = _ib_voltage(buses, bf, tf); vt = _ib_voltage(buses, bt, tt)
                if vf === nothing || vt === nothing || abs(vf) <= 1e-9 || abs(vt) <= 1e-9
                    push!(findings, IndependentBranchFinding(
                        :info, "I.IBRANCH.ANGLE_UNDEFINED", lid, :branch, tf,
                        "Line '$lid' angle check is indeterminate because an endpoint voltage is zero or missing.",
                        Dict("line"=>lid, "from_terminal"=>tf, "to_terminal"=>tt)))
                else
                    dθ = _ib_wrap_pi(angle(vf) - angle(vt))
                    lo = _ib_rating(amin, k); hi = _ib_rating(amax, k)
                    violated = (lo !== nothing && dθ < lo - 1e-9) ||
                               (hi !== nothing && dθ > hi + 1e-9)
                    if violated
                        push!(findings, IndependentBranchFinding(
                            :error, "E.SOL.BRANCH_ANGLE_VIOLATION", lid, :branch, tf,
                            "Line '$lid' signed endpoint angle difference violates its angle window.",
                            Dict("line"=>lid, "from_terminal"=>tf, "to_terminal"=>tt,
                                 "angle_from_minus_to_rad"=>dθ,
                                 "va_diff_min"=>lo, "va_diff_max"=>hi)))
                    end
                end
            end
        end
    end
    return findings
end
