function unspliced(τ, u₀, α, β)
    expᵘ = exp(-β * τ)
    u₀ * expᵘ + α / β * (1 - expᵘ)
end

function spliced(τ, s₀, u₀, α, β, γ)
    c = (α - u₀ * β) / (γ - β)
    expᵘ = exp(-β * τ)
    expˢ = exp(-γ * τ)
    s₀ * expˢ + α / γ * (1 - expˢ) + c * (expˢ - expᵘ)
end

function mRNA(τ, u₀, s₀, α, β, γ)
    expᵘ, expˢ = exp(-β * τ), exp(-γ * τ)
    expᵘˢ = (α - u₀ * β) * inv(γ - β) * (expˢ - expᵘ)
    u = u₀ * expᵘ + α / β * (1 - expᵘ)
    s = s₀ * expˢ + α / γ * (1 - expˢ) + expᵘˢ
    u, s
end
