// Copyright 2017 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "textflag.h"

// The encoding of CRC instructions is as follows (ARMv8 manual F7.1.41)
// Castagnoli: 0xe1SYZ24X
// IEEE:       0xe1SYZ04X
// where
// - X is the input register
// - S is the size (4 for word, 2 for halfword, 0 for byte)
// - Y is the input accumulator reg
// - Z is the output accumulator reg

// castagnoliUpdate updates the non-inverted crc with the given data.

// func castagnoliUpdate(crc uint32, p []byte) uint32
TEXT ·castagnoliUpdate(SB),NOSPLIT,$0-20
	MOVW	crc+0(FP), R1   // CRC value
	MOVW	p+4(FP), R2     // data pointer
	MOVW	p_len+8(FP), R4 // len(p)

	CMP	$4, R4
	BLT	less_than_4

update:
	MOVW.P	4(R2), R5
	WORD	$0xe1411245 // CRC32CW R5, R1
	SUB	$4, R4

	CMP	$4, R4
	BLT	less_than_4

	JMP	update

less_than_4:
	CMP	$2, R4
	BLT	less_than_2

	MOVHU.P	2(R2), R5
	WORD	$0xe1211245 // CRC32CH R5, R1
	SUB	$2, R4

less_than_2:
	CMP	$1, R4
	BLT	done

	MOVBU	(R2), R5
	WORD	$0xe1011245 // CRC32CB R5, R1

done:
	MOVW	R1, ret+16(FP)
	RET

// ieeeUpdate updates the non-inverted crc with the given data.

// func ieeeUpdate(crc uint32, p []byte) uint32
TEXT ·ieeeUpdate(SB),NOSPLIT,$0-20
	MOVW	crc+0(FP), R1   // CRC value
	MOVW	p+4(FP), R2     // data pointer
	MOVW	p_len+8(FP), R4 // len(p)

	CMP	$4, R4
	BLT	less_than_4

update:
	MOVW.P	4(R2), R5
	WORD	$0xe1411045 // CRC32W R5, R1
	SUB	$4, R4

	CMP	$4, R4
	BLT	less_than_4

	JMP	update

less_than_4:
	CMP	$2, R4
	BLT	less_than_2

	MOVHU.P	2(R2), R5
	WORD	$0xe1211045 // CRC32W R5, R1
	SUB	$2, R4

less_than_2:
	CMP	$1, R4
	BLT	done

	MOVBU	(R2), R5
	WORD	$0xe1011045 // CRC32B R5, R1

done:
	MOVW	R1, ret+16(FP)
	RET
