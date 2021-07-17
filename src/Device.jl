module Device
using CUDA

export array_type, array_device, CPU, CUDADevice

struct CPU end
struct CUDADevice end

array_type() = CUDA.has_cuda_gpu() ? CuArray : Array
array_device(::Union{Array}) = CPU()
array_device(::CuArray) = CUDADevice()

end
