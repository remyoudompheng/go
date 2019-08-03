// Copyright 2017 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "textflag.h"

// An encryption round
// VLD1.P 32(R2), {Q1-Q2} # unaligned 32 bytes
// AESE  Q1, Q0
// AESMC Q0, Q0
// AESE  Q2, Q0
// AESMC Q0, Q0
#define ROUND2E \
	WORD	$0xf422220d \
	WORD	$0xf3b00302 \
	WORD	$0xf3b00380 \
	WORD	$0xf3b00304 \
	WORD	$0xf3b00380

// An decryption round
// VLD1.P 32(R2), {Q1-Q2} # unaligned 32 bytes
// AESD  Q1, Q0
// AESIMC Q0, Q0
// AESD  Q2, Q0
// AESIMC Q0, Q0
#define ROUND2D \
	WORD	$0xf422220d \
	WORD	$0xf3b00342 \
	WORD	$0xf3b003c0 \
	WORD	$0xf3b00344 \
	WORD	$0xf3b003c0

// func encryptBlockAsm(nr int, xk *uint32, dst, src *byte)
TEXT ·encryptBlockAsm(SB),NOSPLIT,$0
	MOVW	nr+0(FP), R1
	MOVW	xk+4(FP), R2
	MOVW	dst+8(FP), R3
	MOVW	src+12(FP), R4

	// Load a 16-byte block
	WORD	$0xf4240a0f // VLD1 (R4), Q0

	CMP	$12, R1
	BLT	enc128
	BEQ	enc192
	// 14, 12, 10 times AESE, AESMC (except the last time) then VEOR
enc256:
	ROUND2E
enc192:
	ROUND2E
enc128:
	ROUND2E
	ROUND2E
	ROUND2E
	ROUND2E
	// last rounds
	WORD	$0xf422220d // VLD1.P 32(R2), {Q1-Q2}
	WORD	$0xf3b00302 // AESE Q1, Q0
	WORD	$0xf3b00380 // AESMC Q0, Q0
	WORD	$0xf3b00304 // AESE Q2, Q0
	// xor with the last key
	WORD	$0xf4222a0d // VLD1.P 16(R2), Q1
	WORD	$0xf3020150 // VEOR Q1, Q0
	// store result
	WORD	$0xf4030a0f // VSTR1 Q0, (R3)
	RET

// func decryptBlockAsm(nr int, xk *uint32, dst, src *byte)
TEXT ·decryptBlockAsm(SB),NOSPLIT,$0
	MOVW	nr+0(FP), R1
	MOVW	xk+4(FP), R2
	MOVW	dst+8(FP), R3
	MOVW	src+12(FP), R4

	// Load a 16-byte block
	WORD	$0xf4240a0f // VLD1 (R4), Q0

	CMP	$12, R1
	BLT	dec128
	BEQ	dec192
	// 14, 12, 10 times AESD, AESIMC (except the last time) then VEOR
dec256:
	ROUND2D
dec192:
	ROUND2D
dec128:
	ROUND2D
	ROUND2D
	ROUND2D
	ROUND2D
	// last rounds
	WORD	$0xf422220d // VLD1.P 32(R2), {Q1-Q2}
	WORD	$0xf3b00342 // AESD Q1, Q0
	WORD	$0xf3b003c0 // AESIMC Q0, Q0
	WORD	$0xf3b00344 // AESD Q2, Q0
	// xor with the last key
	WORD	$0xf4222a0d // VLD1.P 16(R2), Q1
	WORD	$0xf3020150 // VEOR Q1, Q0
	// store result
	WORD	$0xf4030a0f // VSTR1 Q0, (R3)
	RET

// func expandKeyAsm(nr int, key *byte, enc, dec *uint32)
// unsupported
TEXT ·expandKeyAsm(SB),NOSPLIT,$0
	B	runtime·sigpanic(SB)

