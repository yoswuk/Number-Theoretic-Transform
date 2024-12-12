/******************************************************************************
 * Integrating the improved Plantard arithmetic into NTTRU.
 *
 * Efficient Plantard arithmetic enables a faster NTTRU implementation with the
 * same stack usage.
 *
 * See the paper at https://eprint.iacr.org/2022/956.pdf for more details.
 *
 * @author   Junhao Huang, BNU-HKBU United International College, Zhuhai, China
 *           jhhuang_nuaa@126.com
 *
 * @date     September 2022
 ******************************************************************************/

#ifndef MACROS_NTT_I
#define MACROS_NTT_I

#include "macros.i"
#include "params.h"

### 7 instructions
.macro doublebutterfly_plant a0, a1, twiddle, tmp, q, qa
	smulwb \tmp, \twiddle, \a1 
	smulwt \a1, \twiddle, \a1
	smlabt \tmp, \tmp, \q, \qa
	smlabt \a1, \a1, \q, \qa
	pkhtb \tmp, \a1, \tmp, asr#16
	usub16 \a1, \a0, \tmp
	uadd16 \a0, \a0, \tmp
.endm

.macro singlebutterfly_plant a0, a1, twiddle, tmp, q, qa
	smulwb \tmp, \twiddle, \a1
	smlabt \tmp, \tmp, \q, \qa
	asr \tmp, \tmp, #16
	usub16 \a1, \a0, \tmp
	uadd16 \a0, \a0, \tmp
.endm

.macro two_doublebutterfly_plant a0, a1, a2, a3, twiddle0, twiddle1, tmp, q, qa
	doublebutterfly_plant \a0, \a1, \twiddle0, \tmp, \q, \qa
	doublebutterfly_plant \a2, \a3, \twiddle1, \tmp, \q, \qa
.endm


.macro two_singlebutterfly_plant a0, a1, a2, a3, twiddle0, twiddle1, tmp, q, qa
	singlebutterfly_plant \a0, \a1, \twiddle0, \tmp, \q, \qa
	singlebutterfly_plant \a2, \a3, \twiddle1, \tmp, \q, \qa
.endm

//a1=a0+a1-0.5q; a0=a0+0.5q; a1:[-1q,1q];a0:[-1q,1q]
.macro double_radix2_plant_first_layer a0, a1, twiddle, tmp, tmp1, q, qa
	smulwb \tmp, \twiddle, \a1
	smulwt \tmp1, \twiddle, \a1
	smlabt \tmp, \tmp, \q, \qa
	smlabt \tmp1, \tmp1, \q, \qa
	pkhtb  \tmp, \tmp1, \tmp, asr#16
	uadd16 \a1, \a0, \a1
	usub16 \a1, \a1, \tmp
	uadd16 \a0, \a0, \tmp
.endm

.macro double_radix3_plant a0, a1, a2, twiddle1, twiddle2, tmp, tmp1, q, qa
	smulwb \tmp, \twiddle1, \a1
	smulwt \tmp1, \twiddle1, \a1
	smlabt \tmp, \tmp, \q, \qa
	smlabt \tmp1, \tmp1, \q, \qa
	pkhtb  \tmp, \tmp1, \tmp, asr#16

	smulwb \tmp, \twiddle2, \a2
	smulwt \tmp1, \twiddle2, \a2
	smlabt \tmp, \tmp, \q, \qa
	smlabt \tmp1, \tmp1, \q, \qa
	pkhtb  \tmp1, \tmp1, \tmp, asr#16	
	
	usub16 \tmp2, \tmp, \tmp1
	smulwb \tmp, \twiddle2, \a2
	smulwt \tmp1, \twiddle2, \a2
	smlabt \tmp, \tmp, \q, \qa
	smlabt \tmp1, \tmp1, \q, \qa
	pkhtb  \tmp1, \tmp1, \tmp, asr#16
	
	usub16 \a2, \a2, \t0
	usub16 \a2, \a2, \t2	
	usub16 \a1, \a1, \t0
	uadd16 \a1, \a2, \t2	
	uadd16 \a0, \a0, \t0
	uadd16 \a0, \a0, \t1
.endm

.macro two_doublebutterfly_plant_first_layer a0, a1, a2, a3, twiddle0, tmp, tmp1, q, qa
	double_radix2_plant_first_layer \a0, \a1, \twiddle0, \tmp, \tmp1, \q, \qa
	double_radix2_plant_first_layer \a2, \a3, \twiddle0, \tmp, \tmp1, \q, \qa
.endm

.macro three_double_radix2_plant_first_layer a0, a1, a2, a3, a4, a5, twiddle0, tmp, tmp1, q, qa
	double_radix2_plant_first_layer \a0, \a1, \twiddle0, \tmp, \tmp1, \q, \qa
	double_radix2_plant_first_layer \a2, \a3, \twiddle0, \tmp, \tmp1, \q, \qa
	double_radix2_plant_first_layer \a4, \a5, \twiddle0, \tmp, \tmp1, \q, \qa
.endm

#endif /* MACROS_NTT_I */
