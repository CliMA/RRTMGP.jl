
function setup_gray_as_pr_grid!(::ClimaComms.CUDADevice, ncol, args...)
    max_threads = 256
    tx = min(ncol, max_threads)
    bx = cld(ncol, tx)
    @cuda always_inline = true threads = (tx) blocks = (bx) _setup_gray_as_pr_grid_kernel!(ncol, args...)
end

function _setup_gray_as_pr_grid_kernel!(ncol, args...)
    gcol = threadIdx().x + (blockIdx().x - 1) * blockDim().x # global id
    if gcol ≤ ncol
        setup_gray_as_pr_grid_kernel!(args..., gcol)
    end
    return nothing
end
