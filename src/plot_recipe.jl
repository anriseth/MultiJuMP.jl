@recipe function f{Tx,To}(md::MultiData{Tx,To})
    seriestype --> :scatter
    xguide --> L"$f_1$"
    yguide --> L"$f_2$"
    zguide --> L"$f_3$"
    label --> ""
    numobjectives = length(md.objectives)
    if numobjectives > 3 || numobjectives < 2
        Base.error("Only plotting 2d and 3d fronts")
    end

    f1arr = convert(Array{Float64},
                    [val[1] for val in md.paretofront])
    f2arr = convert(Array{Float64},
                    [val[2] for val in md.paretofront])

    title --> "Pareto front with $(length(f1arr)) points"

    if numobjectives == 2
        xyz = (f1arr, f2arr)
    else
        f3arr = convert(Array{Float64},
                        [val[3] for val in md.paretofront])
        xyz = (f1arr, f2arr, f3arr)
    end
    (xyz...,)
end

@recipe function f(m::JuMP.Model)
    plot(get_multidata(m))
end
