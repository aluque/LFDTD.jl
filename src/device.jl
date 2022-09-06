const CPU = Val{:CPU}
const GPU = Val{:GPU}

array_type(::CPU, T) = Array{T, 2}
array_type(::GPU, T) = CuArray{T, 2}
vector_type(::CPU, T) = Array{T, 1}
vector_type(::GPU, T) = CuArray{T, 1}

Base.zeros(::CPU, args...) = zeros(args...)
Base.zeros(::GPU, args...) = CUDA.zeros(args...)

