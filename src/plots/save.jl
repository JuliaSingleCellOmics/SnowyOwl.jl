save(plt, path::String, ::Nothing) = nothing
save(plt, path::String, ext) = Plots.savefig(plt, "$(path).$(ext)")

function save(plt, path::String, exts::AbstractVector)
    for ext in exts
        Plots.savefig(plt, "$(path).$(ext)")
    end
end
