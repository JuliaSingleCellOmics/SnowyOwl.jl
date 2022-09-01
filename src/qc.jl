"""
Calculate quality control metrics
"""
function quality_control_metrics!(p::Profile; qc_vars=["mt"], percent_top=nothing, log1p=false)
    describe_obs!(p.obs, p.data, qc_vars, log1p; percent_top=percent_top)
    describe_var!(p.var, p.data, log1p)
    return p
end


function quality_control_metrics(p::Profile, qc_vars=["mt"], percent_top=nothing, log1p=false)
    obs_metrics = describe_obs(
        p.obs;
        qc_vars=qc_vars,
        percent_top=percent_top,
        X = p.data,
        log1p=log1p,
    )
    var_metrics = describe_var(
        p.var;
        X = p.data,
        log1p=log1p,
    )
    return obs_metrics, var_metrics
end


function describe_obs!(obs::DataFrame, X::AbstractMatrix, qc_vars::AbstractVector{String}, log1p::Bool;
    percent_top=nothing)

    #obs_metrics = p.obs
    
    # TODO: implement log1p

    # if percent_top:
    #     percent_top = sorted(percent_top)
    #     proportions = top_segment_proportions(X, percent_top)
    #     for i, n in enumerate(percent_top):
    #         obs_metrics[f"pct_{expr_type}_in_top_{n}_{var_type}"] = (
    #             proportions[:, i] * 100
    #         )
    # for qc_var in qc_vars:
    #     obs_metrics[f"total_{expr_type}_{qc_var}"] = np.ravel(
    #         X[:, adata.var[qc_var].values].sum(axis=1)
    #     )
    #     if log1p:
    #         obs_metrics[f"log1p_total_{expr_type}_{qc_var}"] = np.log1p(
    #             obs_metrics[f"total_{expr_type}_{qc_var}"]
    #         )
    #     obs_metrics[f"pct_{expr_type}_{qc_var}"] = (
    #         obs_metrics[f"total_{expr_type}_{qc_var}"]
    #         / obs_metrics[f"total_{expr_type}"]
    #         * 100
    #     )
    
    #return obs_metrics

end

function describe_var!(var::DataFrame, X::AbstractMatrix, log1p::Bool)

end