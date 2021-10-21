
#pragma once
#ifdef _MSC_VER
#include <intrin.h>
#include <stdexcept>

#define  __builtin_popcount(t) __popcnt(t)
#else
#include <x86intrin.h>
#endif
#define USE_AVX
#if defined(__GNUC__)
#define PORTABLE_ALIGN32 __attribute__((aligned(32)))
#else
#define PORTABLE_ALIGN32 __declspec(align(32))
#endif

#include <cmath>
#include "hnswlib.h"
namespace hnswlib {
	using namespace std;
	static float
		InnerProduct(const void *pVect1, const void *pVect2, const void *qty_ptr)
	{
    const float* a = (const float*)pVect1;
    const float* b = (const float*)pVect2;
    unsigned size = *((unsigned*)qty_ptr);
    float result = 0;
    #pragma omp simd reduction(+:result)
    for (int i = 0; i< size; i++)
        result += a[i]*b[i];
    /*
    #ifdef __GNUC__
    #ifdef __AVX__
          #define AVX_DOT(addr1, addr2, dest, tmp1, tmp2) \
              tmp1 = _mm256_loadu_ps(addr1);\
              tmp2 = _mm256_loadu_ps(addr2);\
              tmp1 = _mm256_mul_ps(tmp1, tmp2); \
              dest = _mm256_add_ps(dest, tmp1);

        __m256 sum;
          __m256 l0, l1;
          __m256 r0, r1;
           unsigned D = (size + 7) & ~7U;
          unsigned DR = D % 16;
          unsigned DD = D - DR;
          const float *l = a;
          const float *r = b;
          const float *e_l = l + DD;
          const float *e_r = r + DD;
          float unpack[8] __attribute__ ((aligned (32))) = {0, 0, 0, 0, 0, 0, 0, 0};

          sum = _mm256_loadu_ps(unpack);
          if(DR){AVX_DOT(e_l, e_r, sum, l0, r0);}

          for (unsigned i = 0; i < DD; i += 16, l += 16, r += 16) {
          AVX_DOT(l, r, sum, l0, r0);
          AVX_DOT(l + 8, r + 8, sum, l1, r1);
          }
          _mm256_storeu_ps(unpack, sum);
          result = unpack[0] + unpack[1] + unpack[2] + unpack[3] + unpack[4] + unpack[5] + unpack[6] + unpack[7];

    #else
    #ifdef __SSE2__
          #define SSE_DOT(addr1, addr2, dest, tmp1, tmp2) \
              tmp1 = _mm128_loadu_ps(addr1);\
              tmp2 = _mm128_loadu_ps(addr2);\
              tmp1 = _mm128_mul_ps(tmp1, tmp2); \
              dest = _mm128_add_ps(dest, tmp1);
          __m128 sum;
          __m128 l0, l1, l2, l3;
          __m128 r0, r1, r2, r3;
          unsigned D = (size + 3) & ~3U;
          unsigned DR = D % 16;
          unsigned DD = D - DR;
          const float *l = a;
          const float *r = b;
          const float *e_l = l + DD;
          const float *e_r = r + DD;
          float unpack[4] __attribute__ ((aligned (16))) = {0, 0, 0, 0};

          sum = _mm_load_ps(unpack);
          switch (DR) {
              case 12:
              SSE_DOT(e_l+8, e_r+8, sum, l2, r2);
              case 8:
              SSE_DOT(e_l+4, e_r+4, sum, l1, r1);
              case 4:
              SSE_DOT(e_l, e_r, sum, l0, r0);
            default:
              break;
          }
          for (unsigned i = 0; i < DD; i += 16, l += 16, r += 16) {
              SSE_DOT(l, r, sum, l0, r0);
              SSE_DOT(l + 4, r + 4, sum, l1, r1);
              SSE_DOT(l + 8, r + 8, sum, l2, r2);
              SSE_DOT(l + 12, r + 12, sum, l3, r3);
          }
          _mm_storeu_ps(unpack, sum);
          result += unpack[0] + unpack[1] + unpack[2] + unpack[3];
    #else

          float dot0, dot1, dot2, dot3;
          const float* last = a + size;
          const float* unroll_group = last - 3;

          // Process 4 items with each loop for efficiency. 
          while (a < unroll_group) {
              dot0 = a[0] * b[0];
              dot1 = a[1] * b[1];
              dot2 = a[2] * b[2];
              dot3 = a[3] * b[3];
              result += dot0 + dot1 + dot2 + dot3;
              a += 4;
              b += 4;
          }
          // Process last 0-3 pixels.  Not needed for standard vector lengths.
          while (a < last) {
              result += *a++ * *b++;
          }
    #endif
    #endif
    #endif
          */
    return -result;
	}
	
 //  static float
	// 	CosineSimilarity(const void *pVect1, const void *pVect2, const void *qty_ptr)
	// {
 //    size_t qty = *((size_t *)qty_ptr);
	// 	float norm1 = 0;
 //    float norm2 = 0;
	// 	for (int i = 0; i < qty; i++) {
	// 		norm1 += ((float*)pVect1)[i] * ((float*)pVect1)[i];
	// 		norm2 += ((float*)pVect2)[i] * ((float*)pVect2)[i];
	// 	}
 //    if(norm1 == 0.0f) norm1 = 1.0f;
 //    if(norm2 == 0.0f) norm2 = 1.0f;

	// 	return InnerProduct(pVect1, pVect2, qty_ptr) / (sqrt(norm1) *  sqrt(norm2));
	// }
	
  class L2Space : public SpaceInterface<float> {
		
		DISTFUNC<float> fstdistfunc_;
		size_t data_size_;
		size_t dim_;
	public:
		L2Space(size_t dim) {
			fstdistfunc_ = InnerProduct;
			//fstdistfunc_ = CosineSimilarity;
			dim_ = dim;
			data_size_ = dim * sizeof(float);
		}

		size_t get_data_size() {
			return data_size_;
		}
		DISTFUNC<float> get_dist_func() {
			return fstdistfunc_;
		}
		void *get_dist_func_param() {
			return &dim_;
		}

    };

}
