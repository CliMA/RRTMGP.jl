module Device

using CUDA

export array_type

array_type() = CUDA.has_cuda_gpu() ? CuArray : Array

end
