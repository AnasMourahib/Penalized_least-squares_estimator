using Statistics

function penalized_loss(theta_pred, theta_true; lambda=1, which=1, mode=:l1)
    base_loss = Statistics.mean(abs.(theta_pred .- theta_true))

    if ndims(theta_pred) == 2 && size(theta_pred,1) >= which
        pen_vals = abs.(theta_pred[which, :])
    else
        pen_vals = abs.(theta_pred)
    end

    pen = (mode == :l1) ? Statistics.mean(pen_vals) : Statistics.mean(pen_vals .^ 2)

    return base_loss + lambda * pen
end


custom_loss(θ_pred, θ_true) = penalized_loss(θ_pred, θ_true)

