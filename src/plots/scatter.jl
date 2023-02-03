scatterplot(x, y, z; kwargs...) =
    Plots.scatter(x, y, z; markersize=1, markerstrokewidth=0, kwargs...)

scatterplot(x, y; figsize=(1000, 600), dpi=300, kwargs...) =
    Plots.scatter(x, y; legend=:outerright, markersize=2, markerstrokewidth=0,
    thickness_scaling=2, widen=false, size=figsize, dpi=dpi, kwargs...)

pca(x, y, z; kwargs...) =
    scatterplot(x, y, z; xlabel="PC1", ylabel="PC2", zlabel="PC3", kwargs...)
pca(x, y; kwargs...) = scatterplot(x, y; xlabel="PC1", ylabel="PC2", kwargs...)

umap(x, y, z; kwargs...) =
    scatterplot(x, y, z; xlabel="UMAP1", ylabel="UMAP2", zlabel="UMAP3", kwargs...)
umap(x, y; kwargs...) = scatterplot(x, y; xlabel="UMAP1", ylabel="UMAP2", kwargs...)
