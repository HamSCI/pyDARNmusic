/* 
*Author: Francis Tholley
*Date: 12/03/2022
*Functionality: For doing intrinsic addition,subtraction, and multiplication to imporve performance
*          in python pyDarnmusic code
*/
#include <stdlib.h>
#include <complex.h>
#include <stdio.h>
#include <emmintrin.h>
#include <immintrin.h>
#include <xmmintrin.h>
#include <tmmintrin.h>
#include <time.h>
//intrinsic addition for integer values
__m128i add4(__m128i a, __m128i b) {
  return _mm_add_epi32(a, b);
}
//intrinsic for loading integer values
__m128i loadu_si128(void *m) {
  return _mm_load_si128((__m128i*) m);
}
//intrinsic for storeing integer values
void storeu_si128(void *m, __m128i res){
    return _mm_store_si128((__m128i*) m, res); 
}
//intrinsic addition for float values
__m128 add4ps(__m128 a, __m128 b) {
  return _mm_add_ps(a, b);
}
//intrinsic multiplication for float values
__m128 mul4ps(__m128 a, __m128 b) {
  return _mm_mul_ps(a, b);
}
//intrinsic subtraction for float values
__m128 sub4ps(__m128 a, __m128 b) {
  return _mm_sub_ps(a, b);
}
//intrinsic for storing float values
void storeu_ps(float *m, __m128 res){
    return _mm_storeu_ps(m, res); 
}
//intrinsic for loading float values
__m128 loadu_ps(float *m) {
  return _mm_loadu_ps(m);
}

int sum_of_array(int *A,int *B, int len) {
    const __m128i vk0 = _mm_set1_epi8(0);       // constant vector of all 0s for use with _mm_unpacklo_epi8/_mm_unpackhi_epi8
    const __m128i vk1 = _mm_set1_epi16(1);      // constant vector of all 1s for use with _mm_madd_epi16
    __m128i vsum = _mm_set1_epi32(0);           // initialise vector of four partial 32 bit sums
    unsigned int sum;
    int i;

    for (i = 0; i < len; i += 4)
    {
        __m128i v = loadu_si128(&A[i]);      // load vector of 8 bit values
        __m128i vl = _mm_unpacklo_epi8(v, vk0); // unpack to two vectors of 16 bit values
        __m128i vh = _mm_unpackhi_epi8(v, vk0);
        vsum = _mm_add_epi32(vsum, _mm_madd_epi16(vl, vk1));
        vsum = _mm_add_epi32(vsum, _mm_madd_epi16(vh, vk1));
                                                // unpack and accumulate 16 bit values to
                                                // 32 bit partial sum vector

    }
    // horizontal add of four 32 bit partial sums and return result
    vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 8));
    vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 4));
    sum = _mm_cvtsi128_si32(vsum);
    return sum;
}
void real(float* real1,float* real2,float* realres,char arithtype, int n){
    /* Function can add,subtract, and multiply real values in two arrays 
    * it takes in three float pointers with two holding
    * the data and one for storing the final result. In addition it takes in a character 
    * for the type of arithmitic and the length of the array
    * ensure that the lenghth of all arrays is the same
    *
    */
    if(arithtype=='a'){
        //ensures length of array is a multiples of 4
        int len = n;
        int end = (len % 4 == 4) ? len : len - 4;
        int i = 0;
        for (; i < end; i += 4) {
            __m128 v1 = loadu_ps(&real1[i]);// loads four values from array 1
            __m128 v2 = loadu_ps(&real2[i]);// loads four values from array 2
            __m128 res = add4ps(v1, v2);  // adds values using vectorization intrinsc
            storeu_ps(&realres[i], res); //stores result in result array
        }
        // Last 4 iterations and we may repeat some.
        int index = len - 4;
        __m128 v1 = loadu_ps(&real1[index]);
        __m128 v2 = loadu_ps(&real2[index]);
        __m128 res = add4ps(v1, v2);
        storeu_ps(&realres[index], res);
    }
    if(arithtype == 's'){
        //ensures length of array is a multiples of 4
        int len = n;
        int end = (len % 4 == 4) ? len : len - 4;
        int i = 0;
        for (; i < end; i += 4) {
            __m128 v1 = loadu_ps(&real1[i]); // loads four values from array 1
            __m128 v2 = loadu_ps(&real2[i]); // loads four values from array 2
            __m128 res = sub4ps(v1, v2); // subtracts values using vectorization intrinsc
            storeu_ps(&realres[i], res); //stores result in result array
        }
        // Last 4 iterations and we may repeat some.
        int index = len - 4;
        __m128 v1 = loadu_ps(&real1[index]);
        __m128 v2 = loadu_ps(&real2[index]);
        __m128 res = sub4ps(v1, v2);
        storeu_ps(&realres[index], res);
    }
    if(arithtype=='m'){
        //ensures length of array is a multiples of 4
        int len = n;
        int end = (len % 4 == 4) ? len : len - 4;
        int i = 0;
        for (; i < end; i += 4) {
            __m128 v1 = loadu_ps(&real1[i]);// loads four values from array 1
            __m128 v2 = loadu_ps(&real2[i]);// loads four values from array 2
            __m128 res = mul4ps(v1, v2); // multiplys values using vectorization intrinsc
            storeu_ps(&realres[i], res); //stores result in result array
        }
        // Last 4 iterations and we may repeat some.
        int index = len - 4;
        __m128 v1 = loadu_ps(&real1[index]);
        __m128 v2 = loadu_ps(&real2[index]);
        __m128 res = mul4ps(v1, v2);
        storeu_ps(&realres[index], res);
    }
 

}
void img(float* img1,float* img2, float* imgres, char arithtype,int n){
    /* Function can add,subtract, and multiply imaginary values in two arrays 
    * it takes in three float pointers with two holding
    * the data and one for storing the final result. In addition it takes in a character 
    * for the type of arithmitic and the length of the array
    * ensure that the lenghth of all arrays is the same
    *
    */
    if(arithtype == 'a'){
        // For imaginary
        //ensures length of array is a multiples of 4
        int len = n;
        int end = (len % 4 == 4) ? len : len - 4;
        int i = 0;
        
        for (; i < end; i += 4) {
            __m128 v3 = loadu_ps(&img1[i]); // loads four values from array 1
            __m128 v4 = loadu_ps(&img2[i]); // loads four values from array 2
            __m128 ires = add4ps(v3, v4); // adds values using vectorization intrinsc
            storeu_ps(&imgres[i], ires); //stores result in result array
        }
        // Last 4 iterations and we may repeat some.
        int index = len - 4;
        __m128 v3 = loadu_ps(&img1[index]);
        __m128 v4 = loadu_ps(&img2[index]);
        __m128 ires = add4ps(v3, v4);
        storeu_ps(&imgres[index], ires);
    }
    if(arithtype == 's'){
        // For imaginary
        int len = n;
        int end = (len % 4 == 4) ? len : len - 4;
        int i = 0;
        
        for (; i < end; i += 4) {
            __m128 v3 = loadu_ps(&img1[i]);
            __m128 v4 = loadu_ps(&img2[i]);
            __m128 ires = sub4ps(v3, v4);
            storeu_ps(&imgres[i], ires);
        }
        // Last 4 iterations and we may repeat some.
        int index = len - 4;
        __m128 v3 = loadu_ps(&img1[index]);
        __m128 v4 = loadu_ps(&img2[index]);
        __m128 ires = sub4ps(v3, v4);
        storeu_ps(&imgres[index], ires);
    }

    if(arithtype == 'm'){
        // For imaginary
        int len = n;
        int end = (len % 4 == 4) ? len : len - 4;
        int i = 0;
        
        for (; i < end; i += 4) {
            __m128 v3 = loadu_ps(&img1[i]);
            __m128 v4 = loadu_ps(&img2[i]);
            __m128 ires = mul4ps(v3, v4);
            storeu_ps(&imgres[i], ires);
        }
        // Last 4 iterations and we may repeat some.
        int index = len - 4;
        __m128 v3 = loadu_ps(&img1[index]);
        __m128 v4 = loadu_ps(&img2[index]);
        __m128 ires = mul4ps(v3, v4);
        storeu_ps(&imgres[index], ires);
    }
 
}
int c_sum(_Complex float* matrix1, _Complex float* matrix2,_Complex float* finalresult, int n){
    /* Function add complex values in two arrays 
    * it takes in three complex type of float pointers with two holding
    * the data and one for storing the final result. In addition it takes in the length of the array
    * ensure that the lenghth of all arrays is the same
    *
    */
    float * real1 = (float *)malloc(sizeof(float) * n);
    float * img1 = (float *)malloc(sizeof(float) * n);
    float * real2 = (float *)malloc(sizeof(float) * n);
    float * img2 = (float *)malloc(sizeof(float) * n);
    float * realres = (float *)malloc(sizeof(float) * n);
    float * imgres = (float *)malloc(sizeof(float) * n);

    // separates real parts of complex values in both data arrays
    // separates imaginary parts of complex values in both data arrays
    for(int i = 0; i < n; i++){
        real1[i] = creal(matrix1[i]);
        img1[i] = cimag(matrix1[i]);
        real2[i] = creal(matrix1[i]);
        img2[i] = cimag(matrix1[i]);        
    }
    char arithtype= 'a'; // used to specify which type of arithematic operation to do (a --> addition)
    real(real1,real2,realres,arithtype,n);// add real values from both arrays
    img(img1,img2,imgres,arithtype,n); // add imaginary values from both arrays

    //combines the real part and the imaginary part and stores values in final result
    for(int i = 0; i < n; i++){
        finalresult[i] = realres[i] + imgres[i]*_Complex_I;
    }

    for(int i = 0; i < n; i++){
        printf("%.9f + %.8fj\n",creal(finalresult[i]),cimag(finalresult[i]));
    }
    //releases unneeded memory
    free(real1);
    free(img1);
    free(real2);
    free(img2);
    return 0;
}
int c_sub(_Complex float* matrix1, _Complex float* matrix2,_Complex float* finalresult, int n){
    /* Function subtracts complex values in two arrays 
    * it takes in three complex type of float pointers with two holding
    * the data and one for storing the final result. In addition it takes in the length of the array
    * ensure that the lenghth of all arrays is the same
    *
    */
    float * real1 = (float *)malloc(sizeof(float) * n);
    float * img1 = (float *)malloc(sizeof(float) * n);
    float * real2 = (float *)malloc(sizeof(float) * n);
    float * img2 = (float *)malloc(sizeof(float) * n);
    float * realres = (float *)malloc(sizeof(float) * n);
    float * imgres = (float *)malloc(sizeof(float) * n);

    // separates real parts of complex values in both data arrays
    // separates imaginary parts of complex values in both data arrays
    for(int i = 0; i < n; i++){
        real1[i] = creal(matrix1[i]);
        img1[i] = cimag(matrix1[i]);
        real2[i] = creal(matrix1[i]);
        img2[i] = cimag(matrix1[i]);        
    }
    char arithtype= 's'; // used to specify which type of arithematic operation to do (s --> subtraction)
    real(real1,real2,realres,arithtype,n); // subtract real values from both arrays
    img(img1,img2,imgres,arithtype,n);  // subtract imaginary values from both arrays

    //combines the real part and the imaginary part and stores values in final result
    for(int i = 0; i < n; i++){
        finalresult[i] = realres[i]+ imgres[i]*_Complex_I;
    }

    for(int i = 0; i < n; i++){
        printf("%.9f + %.8fj\n",creal(finalresult[i]),cimag(finalresult[i]));
    }
    //releases unneeded memory
    free(real1);
    free(img1);
    free(real2);
    free(img2);
    free(realres);
    free(imgres);
    return 0;
}
void real_times_imag(float* real1,float* img2, float* real1imag2res,int n){
    /* Function can multiply real values with imaginary values in two arrays 
    * it takes in three float pointers with two holding
    * the data and one for storing the final result. In addition it takes in a character 
    * for the type of arithmitic and the length of the array
    * ensure that the lenghth of all arrays is the same
    *
    */
        //ensures length of array is a multiples of 4
        int len = n;
        int end = (len % 4 == 4) ? len : len - 4;
        int i = 0;
        
        for (; i < end; i += 4) {
            __m128 v3 = loadu_ps(&real1[i]); // loads four values from array 1
            __m128 v4 = loadu_ps(&img2[i]); // loads four values from array 2
            __m128 ires = mul4ps(v3, v4); // multply values using vectorization intrinsc
            storeu_ps(&real1imag2res[i], ires); //stores result in result array
        }
        // Last 4 iterations and we may repeat some.
        int index = len - 4;
        __m128 v3 = loadu_ps(&real1[index]);
        __m128 v4 = loadu_ps(&img2[index]);
        __m128 ires = mul4ps(v3, v4);
        storeu_ps(&real1imag2res[index], ires);
}

int c_mul(_Complex float* matrix1, _Complex float* matrix2,_Complex float* finalresult, int n){
    /* Function takes in three complex type of float pointers with two holding
    * the data and one for storing the final result. In addition it takes in the length of the array
    * ensure that the lenghth of all arrays is the same
    *
    */
    float * real1 = (float *)malloc(sizeof(float) * n);
    float * img1 = (float *)malloc(sizeof(float) * n);
    float * real2 = (float *)malloc(sizeof(float) * n);
    float * img2 = (float *)malloc(sizeof(float) * n);
    float * realres = (float *)malloc(sizeof(float) * n);
    float * imgres = (float *)malloc(sizeof(float) * n);
    float* real1imag2res = (float *)malloc(sizeof(float) * n);
    float* real2imag1res = (float *)malloc(sizeof(float) * n);
    float* finalimg1res = (float *)malloc(sizeof(float) * n);
    float* finalimg2res = (float *)malloc(sizeof(float) * n);
    float* realpart = (float *)malloc(sizeof(float) * n);
    float* imgpart = (float *)malloc(sizeof(float) * n);

    // separates real parts of complex values in both data arrays
    // separates imaginary parts of complex values in both data arrays
    for(int i = 0; i < n; i++){
        real1[i] = creal(matrix1[i]);
        img1[i] = cimag(matrix1[i]);
        real2[i] = creal(matrix1[i]);
        img2[i] = cimag(matrix1[i]);        
    }
    char arithtype= 'm'; // used to specify which type of arithematic operation to do (m --> multiplication)
    // (x+yi)(u+vi) = (xu - yv)(xv + yu)i
    real(real1,real2,realres,arithtype,n); //xu
    img(img1,img2,imgres,arithtype,n);   //yv
    real_times_imag(real1,img2,real1imag2res,n);  //xvi
    real_times_imag(real2,img1,real2imag1res,n);  //yui

    arithtype ='s'; // (s --> subtraction)
    real(realres,imgres,realpart,arithtype,n);

    arithtype='a'; // (a --> subtraction)
    img(real1imag2res,real2imag1res,imgpart,arithtype,n);

    //combines the real part and the imaginary part and stores values in final result
    for(int i = 0; i < n; i++){
        finalresult[i] = realpart[i] + imgpart[i]*_Complex_I;
    }

    for(int i = 0; i < n; i++){
        printf("%.9f + %.8fj\n",creal(finalresult[i]),cimag(finalresult[i]));
    }
    // releases unneeded pointers in memory
    free(real1);
    free(img1);
    free(real2);
    free(img2);
    return 0;
}
