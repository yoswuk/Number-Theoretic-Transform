	.section	__TEXT,__text,regular,pure_instructions
	.build_version macos, 15, 0	sdk_version 15, 1
	.globl	_plantard                       ; -- Begin function plantard
	.p2align	2
_plantard:                              ; @plantard
	.cfi_startproc
; %bb.0:
	sub	sp, sp, #16
	.cfi_def_cfa_offset 16
	strh	w0, [sp, #14]
	strh	w1, [sp, #12]
	ldrsh	w8, [sp, #14]
	ldrsh	w9, [sp, #12]
	mul	w8, w8, w9
	mov	w9, #12929                      ; =0x3281
	movk	w9, #29782, lsl #16
	mul	w8, w8, w9
	str	w8, [sp, #8]gcc
	ldr	w8, [sp, #8]
	asr	w8, w8, #16
	str	w8, [sp, #8]
	ldr	w8, [sp, #8]
	add	w8, w8, #8
	str	w8, [sp, #8]
	ldr	w8, [sp, #8]
	mov	w9, #3457                       ; =0xd81
	mul	w8, w8, w9
	asr	w8, w8, #16
	str	w8, [sp, #8]
	ldr	w8, [sp, #8]
	sxth	w0, w8
	add	sp, sp, #16
	ret
	.cfi_endproc
                                        ; -- End function
.subsections_via_symbols
