#pragma once

#include "../core/types.hpp"
#include <immintrin.h>
#include <cstring>

namespace math::optimization {

template<typename T>
struct SIMDTraits {};

template<>
struct SIMDTraits<float> {
    using VectorType = __m256;
    static constexpr core::Size Width = 8;
    static constexpr core::Size Alignment = 32;
    
    static VectorType load(const float* ptr) { return _mm256_load_ps(ptr); }
    static VectorType loadUnaligned(const float* ptr) { return _mm256_loadu_ps(ptr); }
    static void store(float* ptr, VectorType vec) { _mm256_store_ps(ptr, vec); }
    static void storeUnaligned(float* ptr, VectorType vec) { _mm256_storeu_ps(ptr, vec); }
    
    static VectorType add(VectorType a, VectorType b) { return _mm256_add_ps(a, b); }
    static VectorType sub(VectorType a, VectorType b) { return _mm256_sub_ps(a, b); }
    static VectorType mul(VectorType a, VectorType b) { return _mm256_mul_ps(a, b); }
    static VectorType div(VectorType a, VectorType b) { return _mm256_div_ps(a, b); }
    static VectorType fma(VectorType a, VectorType b, VectorType c) { return _mm256_fmadd_ps(a, b, c); }
    
    static VectorType sqrt(VectorType a) { return _mm256_sqrt_ps(a); }
    static VectorType exp(VectorType a) { return _mm256_exp_ps(a); }
    static VectorType log(VectorType a) { return _mm256_log_ps(a); }
    
    static VectorType set1(float value) { return _mm256_set1_ps(value); }
    static VectorType setzero() { return _mm256_setzero_ps(); }
};

template<>
struct SIMDTraits<double> {
    using VectorType = __m256d;
    static constexpr core::Size Width = 4;
    static constexpr core::Size Alignment = 32;
    
    static VectorType load(const double* ptr) { return _mm256_load_pd(ptr); }
    static VectorType loadUnaligned(const double* ptr) { return _mm256_loadu_pd(ptr); }
    static void store(double* ptr, VectorType vec) { _mm256_store_pd(ptr, vec); }
    static void storeUnaligned(double* ptr, VectorType vec) { _mm256_storeu_pd(ptr, vec); }
    
    static VectorType add(VectorType a, VectorType b) { return _mm256_add_pd(a, b); }
    static VectorType sub(VectorType a, VectorType b) { return _mm256_sub_pd(a, b); }
    static VectorType mul(VectorType a, VectorType b) { return _mm256_mul_pd(a, b); }
    static VectorType div(VectorType a, VectorType b) { return _mm256_div_pd(a, b); }
    static VectorType fma(VectorType a, VectorType b, VectorType c) { return _mm256_fmadd_pd(a, b, c); }
    
    static VectorType sqrt(VectorType a) { return _mm256_sqrt_pd(a); }
    static VectorType exp(VectorType a) { return _mm256_exp_pd(a); }
    static VectorType log(VectorType a) { return _mm256_log_pd(a); }
    
    static VectorType set1(double value) { return _mm256_set1_pd(value); }
    static VectorType setzero() { return _mm256_setzero_pd(); }
};

#ifdef __AVX512F__
template<>
struct SIMDTraits<double> {
    using VectorType = __m512d;
    static constexpr core::Size Width = 8;
    static constexpr core::Size Alignment = 64;
    
    static VectorType load(const double* ptr) { return _mm512_load_pd(ptr); }
    static VectorType loadUnaligned(const double* ptr) { return _mm512_loadu_pd(ptr); }
    static void store(double* ptr, VectorType vec) { _mm512_store_pd(ptr, vec); }
    static void storeUnaligned(double* ptr, VectorType vec) { _mm512_storeu_pd(ptr, vec); }
    
    static VectorType add(VectorType a, VectorType b) { return _mm512_add_pd(a, b); }
    static VectorType sub(VectorType a, VectorType b) { return _mm512_sub_pd(a, b); }
    static VectorType mul(VectorType a, VectorType b) { return _mm512_mul_pd(a, b); }
    static VectorType div(VectorType a, VectorType b) { return _mm512_div_pd(a, b); }
    static VectorType fma(VectorType a, VectorType b, VectorType c) { return _mm512_fmadd_pd(a, b, c); }
    
    static VectorType sqrt(VectorType a) { return _mm512_sqrt_pd(a); }
    static VectorType exp(VectorType a) { return _mm512_exp_pd(a); }
    static VectorType log(VectorType a) { return _mm512_log_pd(a); }
    
    static VectorType set1(double value) { return _mm512_set1_pd(value); }
    static VectorType setzero() { return _mm512_setzero_pd(); }
};
#endif

template<typename T>
class SIMDVector {
private:
    using Traits = SIMDTraits<T>;
    using VectorType = typename Traits::VectorType;
    
    alignas(Traits::Alignment) T* data_;
    core::Size size_;
    core::Size capacity_;
    
public:
    explicit SIMDVector(core::Size capacity = 0) : data_(nullptr), size_(0), capacity_(0) {
        if (capacity > 0) {
            reserve(capacity);
        }
    }
    
    ~SIMDVector() {
        if (data_) {
            std::free(data_);
        }
    }
    
    SIMDVector(const SIMDVector&) = delete;
    SIMDVector& operator=(const SIMDVector&) = delete;
    
    SIMDVector(SIMDVector&& other) noexcept
        : data_(other.data_), size_(other.size_), capacity_(other.capacity_) {
        other.data_ = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
    }
    
    SIMDVector& operator=(SIMDVector&& other) noexcept {
        if (this != &other) {
            if (data_) std::free(data_);
            data_ = other.data_;
            size_ = other.size_;
            capacity_ = other.capacity_;
            other.data_ = nullptr;
            other.size_ = 0;
            other.capacity_ = 0;
        }
        return *this;
    }
    
    void reserve(core::Size newCapacity) {
        if (newCapacity > capacity_) {
            const core::Size alignedCapacity = ((newCapacity + Traits::Width - 1) / Traits::Width) * Traits::Width;
            T* newData = static_cast<T*>(std::aligned_alloc(Traits::Alignment, alignedCapacity * sizeof(T)));
            
            if (!newData) throw std::bad_alloc();
            
            if (data_) {
                std::memcpy(newData, data_, size_ * sizeof(T));
                std::free(data_);
            }
            
            data_ = newData;
            capacity_ = alignedCapacity;
        }
    }
    
    void resize(core::Size newSize) {
        if (newSize > capacity_) {
            reserve(newSize * 2);
        }
        size_ = newSize;
    }
    
    void push_back(const T& value) {
        if (size_ >= capacity_) {
            reserve(capacity_ == 0 ? Traits::Width : capacity_ * 2);
        }
        data_[size_++] = value;
    }
    
    T& operator[](core::Size index) { return data_[index]; }
    const T& operator[](core::Size index) const { return data_[index]; }
    
    T* data() { return data_; }
    const T* data() const { return data_; }
    
    core::Size size() const { return size_; }
    core::Size capacity() const { return capacity_; }
    bool empty() const { return size_ == 0; }
    
    void add(const SIMDVector& other) {
        const core::Size minSize = std::min(size_, other.size_);
        const core::Size vectorizedSize = (minSize / Traits::Width) * Traits::Width;
        
        for (core::Size i = 0; i < vectorizedSize; i += Traits::Width) {
            const VectorType a = Traits::load(&data_[i]);
            const VectorType b = Traits::load(&other.data_[i]);
            const VectorType result = Traits::add(a, b);
            Traits::store(&data_[i], result);
        }
        
        for (core::Size i = vectorizedSize; i < minSize; ++i) {
            data_[i] += other.data_[i];
        }
    }
    
    void multiply(const SIMDVector& other) {
        const core::Size minSize = std::min(size_, other.size_);
        const core::Size vectorizedSize = (minSize / Traits::Width) * Traits::Width;
        
        for (core::Size i = 0; i < vectorizedSize; i += Traits::Width) {
            const VectorType a = Traits::load(&data_[i]);
            const VectorType b = Traits::load(&other.data_[i]);
            const VectorType result = Traits::mul(a, b);
            Traits::store(&data_[i], result);
        }
        
        for (core::Size i = vectorizedSize; i < minSize; ++i) {
            data_[i] *= other.data_[i];
        }
    }
    
    void scale(T scalar) {
        const core::Size vectorizedSize = (size_ / Traits::Width) * Traits::Width;
        const VectorType scalarVec = Traits::set1(scalar);
        
        for (core::Size i = 0; i < vectorizedSize; i += Traits::Width) {
            const VectorType a = Traits::load(&data_[i]);
            const VectorType result = Traits::mul(a, scalarVec);
            Traits::store(&data_[i], result);
        }
        
        for (core::Size i = vectorizedSize; i < size_; ++i) {
            data_[i] *= scalar;
        }
    }
    
    T sum() const {
        const core::Size vectorizedSize = (size_ / Traits::Width) * Traits::Width;
        VectorType accum = Traits::setzero();
        
        for (core::Size i = 0; i < vectorizedSize; i += Traits::Width) {
            const VectorType a = Traits::load(&data_[i]);
            accum = Traits::add(accum, a);
        }
        
        alignas(Traits::Alignment) T temp[Traits::Width];
        Traits::store(temp, accum);
        
        T result = T{0};
        for (core::Size i = 0; i < Traits::Width; ++i) {
            result += temp[i];
        }
        
        for (core::Size i = vectorizedSize; i < size_; ++i) {
            result += data_[i];
        }
        
        return result;
    }
    
    T dot(const SIMDVector& other) const {
        const core::Size minSize = std::min(size_, other.size_);
        const core::Size vectorizedSize = (minSize / Traits::Width) * Traits::Width;
        VectorType accum = Traits::setzero();
        
        for (core::Size i = 0; i < vectorizedSize; i += Traits::Width) {
            const VectorType a = Traits::load(&data_[i]);
            const VectorType b = Traits::load(&other.data_[i]);
            const VectorType product = Traits::mul(a, b);
            accum = Traits::add(accum, product);
        }
        
        alignas(Traits::Alignment) T temp[Traits::Width];
        Traits::store(temp, accum);
        
        T result = T{0};
        for (core::Size i = 0; i < Traits::Width; ++i) {
            result += temp[i];
        }
        
        for (core::Size i = vectorizedSize; i < minSize; ++i) {
            result += data_[i] * other.data_[i];
        }
        
        return result;
    }
};

class VectorizedMath {
public:
    template<typename T>
    static void add(const T* a, const T* b, T* result, core::Size count) {
        using Traits = SIMDTraits<T>;
        const core::Size vectorizedCount = (count / Traits::Width) * Traits::Width;
        
        for (core::Size i = 0; i < vectorizedCount; i += Traits::Width) {
            const auto vecA = Traits::loadUnaligned(&a[i]);
            const auto vecB = Traits::loadUnaligned(&b[i]);
            const auto vecResult = Traits::add(vecA, vecB);
            Traits::storeUnaligned(&result[i], vecResult);
        }
        
        for (core::Size i = vectorizedCount; i < count; ++i) {
            result[i] = a[i] + b[i];
        }
    }
    
    template<typename T>
    static void multiply(const T* a, const T* b, T* result, core::Size count) {
        using Traits = SIMDTraits<T>;
        const core::Size vectorizedCount = (count / Traits::Width) * Traits::Width;
        
        for (core::Size i = 0; i < vectorizedCount; i += Traits::Width) {
            const auto vecA = Traits::loadUnaligned(&a[i]);
            const auto vecB = Traits::loadUnaligned(&b[i]);
            const auto vecResult = Traits::mul(vecA, vecB);
            Traits::storeUnaligned(&result[i], vecResult);
        }
        
        for (core::Size i = vectorizedCount; i < count; ++i) {
            result[i] = a[i] * b[i];
        }
    }
    
    template<typename T>
    static void fma(const T* a, const T* b, const T* c, T* result, core::Size count) {
        using Traits = SIMDTraits<T>;
        const core::Size vectorizedCount = (count / Traits::Width) * Traits::Width;
        
        for (core::Size i = 0; i < vectorizedCount; i += Traits::Width) {
            const auto vecA = Traits::loadUnaligned(&a[i]);
            const auto vecB = Traits::loadUnaligned(&b[i]);
            const auto vecC = Traits::loadUnaligned(&c[i]);
            const auto vecResult = Traits::fma(vecA, vecB, vecC);
            Traits::storeUnaligned(&result[i], vecResult);
        }
        
        for (core::Size i = vectorizedCount; i < count; ++i) {
            result[i] = a[i] * b[i] + c[i];
        }
    }
    
    template<typename T>
    static void scale(const T* input, T scalar, T* result, core::Size count) {
        using Traits = SIMDTraits<T>;
        const core::Size vectorizedCount = (count / Traits::Width) * Traits::Width;
        const auto scalarVec = Traits::set1(scalar);
        
        for (core::Size i = 0; i < vectorizedCount; i += Traits::Width) {
            const auto vecInput = Traits::loadUnaligned(&input[i]);
            const auto vecResult = Traits::mul(vecInput, scalarVec);
            Traits::storeUnaligned(&result[i], vecResult);
        }
        
        for (core::Size i = vectorizedCount; i < count; ++i) {
            result[i] = input[i] * scalar;
        }
    }
    
    template<typename T>
    static T sum(const T* input, core::Size count) {
        using Traits = SIMDTraits<T>;
        const core::Size vectorizedCount = (count / Traits::Width) * Traits::Width;
        auto accum = Traits::setzero();
        
        for (core::Size i = 0; i < vectorizedCount; i += Traits::Width) {
            const auto vecInput = Traits::loadUnaligned(&input[i]);
            accum = Traits::add(accum, vecInput);
        }
        
        alignas(Traits::Alignment) T temp[Traits::Width];
        Traits::store(temp, accum);
        
        T result = T{0};
        for (core::Size i = 0; i < Traits::Width; ++i) {
            result += temp[i];
        }
        
        for (core::Size i = vectorizedCount; i < count; ++i) {
            result += input[i];
        }
        
        return result;
    }
    
    template<typename T>
    static T dot(const T* a, const T* b, core::Size count) {
        using Traits = SIMDTraits<T>;
        const core::Size vectorizedCount = (count / Traits::Width) * Traits::Width;
        auto accum = Traits::setzero();
        
        for (core::Size i = 0; i < vectorizedCount; i += Traits::Width) {
            const auto vecA = Traits::loadUnaligned(&a[i]);
            const auto vecB = Traits::loadUnaligned(&b[i]);
            const auto product = Traits::mul(vecA, vecB);
            accum = Traits::add(accum, product);
        }
        
        alignas(Traits::Alignment) T temp[Traits::Width];
        Traits::store(temp, accum);
        
        T result = T{0};
        for (core::Size i = 0; i < Traits::Width; ++i) {
            result += temp[i];
        }
        
        for (core::Size i = vectorizedCount; i < count; ++i) {
            result += a[i] * b[i];
        }
        
        return result;
    }
    
    template<typename T>
    static void exp(const T* input, T* result, core::Size count) {
        using Traits = SIMDTraits<T>;
        const core::Size vectorizedCount = (count / Traits::Width) * Traits::Width;
        
        for (core::Size i = 0; i < vectorizedCount; i += Traits::Width) {
            const auto vecInput = Traits::loadUnaligned(&input[i]);
            const auto vecResult = Traits::exp(vecInput);
            Traits::storeUnaligned(&result[i], vecResult);
        }
        
        for (core::Size i = vectorizedCount; i < count; ++i) {
            result[i] = std::exp(input[i]);
        }
    }
    
    template<typename T>
    static void sqrt(const T* input, T* result, core::Size count) {
        using Traits = SIMDTraits<T>;
        const core::Size vectorizedCount = (count / Traits::Width) * Traits::Width;
        
        for (core::Size i = 0; i < vectorizedCount; i += Traits::Width) {
            const auto vecInput = Traits::loadUnaligned(&input[i]);
            const auto vecResult = Traits::sqrt(vecInput);
            Traits::storeUnaligned(&result[i], vecResult);
        }
        
        for (core::Size i = vectorizedCount; i < count; ++i) {
            result[i] = std::sqrt(input[i]);
        }
    }
};

class MatrixOperations {
public:
    template<typename T>
    static void matrixMultiply(const T* A, const T* B, T* C, 
                              core::Size M, core::Size N, core::Size K) {
        using Traits = SIMDTraits<T>;
        const core::Size vectorWidth = Traits::Width;
        
        for (core::Size i = 0; i < M; ++i) {
            for (core::Size j = 0; j < N; j += vectorWidth) {
                auto sum = Traits::setzero();
                
                for (core::Size k = 0; k < K; ++k) {
                    const auto a_ik = Traits::set1(A[i * K + k]);
                    const auto b_kj = Traits::loadUnaligned(&B[k * N + j]);
                    sum = Traits::fma(a_ik, b_kj, sum);
                }
                
                Traits::storeUnaligned(&C[i * N + j], sum);
            }
            
            for (core::Size j = (N / vectorWidth) * vectorWidth; j < N; ++j) {
                T sum = T{0};
                for (core::Size k = 0; k < K; ++k) {
                    sum += A[i * K + k] * B[k * N + j];
                }
                C[i * N + j] = sum;
            }
        }
    }
    
    template<typename T>
    static void matrixTranspose(const T* input, T* output, core::Size rows, core::Size cols) {
        using Traits = SIMDTraits<T>;
        const core::Size tileSize = 8;
        
        for (core::Size i = 0; i < rows; i += tileSize) {
            for (core::Size j = 0; j < cols; j += tileSize) {
                const core::Size maxI = std::min(i + tileSize, rows);
                const core::Size maxJ = std::min(j + tileSize, cols);
                
                for (core::Size ii = i; ii < maxI; ++ii) {
                    for (core::Size jj = j; jj < maxJ; ++jj) {
                        output[jj * rows + ii] = input[ii * cols + jj];
                    }
                }
            }
        }
    }
    
    template<typename T>
    static void matrixAdd(const T* A, const T* B, T* C, core::Size rows, core::Size cols) {
        const core::Size totalElements = rows * cols;
        VectorizedMath::add(A, B, C, totalElements);
    }
    
    template<typename T>
    static void matrixScale(const T* input, T scalar, T* output, core::Size rows, core::Size cols) {
        const core::Size totalElements = rows * cols;
        VectorizedMath::scale(input, scalar, output, totalElements);
    }
};

using SIMDVectorF = SIMDVector<float>;
using SIMDVectorD = SIMDVector<double>;

}