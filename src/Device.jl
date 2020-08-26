module Device
using KernelAbstractions
using CUDA

export array_type, array_device

array_type() = CUDA.has_cuda_gpu() ? CuArray : Array
array_device(::Union{Array}) = CPU()
array_device(::CuArray) = CUDADevice()

end
