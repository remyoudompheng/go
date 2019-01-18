// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package strconv

import (
	"math/bits"
)

/*
This file implements conversion to decimals using the techniques
from "Ryū: Fast Float-to-String Conversion" (doi:10.1145/3192366.3192369)
by Ulf Adams.

The reference implementation is C code licensed under the Apache 2.0
license, found at https://github.com/ulfjack/ryu

The conversion problem is as follows:

1. Let x be a floating-point number to convert.

2. Compute bounds l < x < u such that any number
   in the (l, u) interval rounds back to x. The bounds may
   be inclusive.

3. Write
   l = (A + ε_A) × 10^q
   x = (B + ε_B) × 10^q
   u = (C + ε_C) × 10^q

   Exponent q is chosen such that C-A >= 2.
   Then the knowledge of integers A, B, C, and of the zeroness of the ε_i
   allow to determine the final decimal digits.

The implementation is as follows:
- determine q and precision k such that fixed-precision multiplication
  (and right-shift) gives correct values for A, B, C. This is a hard
  precomputation.
  It turns out that a precision of 128 bits is enough for IEEE-754 64-bit floats.
- compute exactly whether ε_A, ε_B, ε_C are zero
  If q > 0, this implies that the fixed precision value for 10^q was exact.
  If q < 0, this implies that the mantissa was actually divisible by 5^q.
  Looking at the lower bits of the extended precision product gives the final answer.
- if C is not allowed (that is, u does not round to x, and ε_C == 0),
  replace C by C-1
- produce digits using A, B, C
  This is done by choosing the most "rounded" version of C: 10^k * (C / 10^k)
  which is still larger than A
*/

func RyuShortest(d *decimalSlice, mant uint64, exp int) {
	if mant == 0 {
		d.nd, d.dp = 0, 0
		return
	}
	ml, mc, mu, e2 := computeBounds(mant, exp)
	if e2 == 0 {
		ryuDigits(d, ml, mc, mu, true, true, false)
		return
	}
	var q int
	var pow *extfloat128 // a representation of 10^q
	if e2 < 0 {
		// Find 10^q *larger* than 2^-e2
		e := uint(-e2)
		q = int(Exp2toExponent10(e) + 1)
		pow = &RyuPowersOfTen[q]
	} else {
		// Divide by a power of 10 strictly less than 2^e2
		q = int(Exp2toExponent10(uint(e2)) - 1)
		if q < 0 {
			q = 0
		}
		pow = &RyuInvPowersOfTen[q]
		q = -q
	}
	// We are going to multiply by 10^q using 128-bit arithmetic.
	// Is it an exact computation?
	lexact, cexact, uexact := false, false, false
	switch {
	case q > 55:
		// large positive powers of ten are not exact
	case 54 >= q && q >= 0:
		lexact, cexact, uexact = true, true, true
	case 0 > q && q > -25:
		// division by a power of ten might be exact
		// if mantissas are multiples of 5. This is because
		// the inverse powers were correctly rounded up,
		// and because of the fundamental property that
		// no extra carry may happen.
		// 5^25 is a 59-bit number so is out of range of ml, mc, mu.
		if divisibleByPower5(ml, -q) {
			lexact = true
		}
		if divisibleByPower5(mc, -q) {
			cexact = true
		}
		if divisibleByPower5(mu, -q) {
			uexact = true
		}
		// In this particular case, the *binary* mantissa
		// of m/10^q will have fewer bits. As a consequence,
		// lower bits of the computed m/10^q *must* be ignored
		// (see below).
	default:
		// division by 10^q (q >= 25) cannot be exact.
	}

	// Compute the decimal mantissas (Floor((l, c, u)*10^q)
	dl, dl0 := ryuMultiply(ml, pow.Hi, pow.Lo)
	dc, dc0 := ryuMultiply(mc, pow.Hi, pow.Lo)
	du, du0 := ryuMultiply(mu, pow.Hi, pow.Lo)
	// If computation was an exact division, lower bits must be ignored.
	if q < 0 {
		if lexact {
			dl0 = true
		}
		if cexact {
			dc0 = true
		}
		if uexact {
			du0 = true
		}
	}
	// The 64 upper bits of the product are more than needed for
	// floor((l,c,u)*10^q). The floating-point exponent of
	// the product e2+pow.Exp+64+55.
	extra := uint(-(e2 + pow.Exp + 64 + 55))
	extramask := uint64(1<<extra - 1)
	// Now compute the floored, integral base 10 mantissas.
	dl, fracl := dl>>extra, dl&extramask
	dc, fracc := dc>>extra, dc&extramask
	du, fracu := du>>extra, du&extramask
	// Is it allowed to use 'du' as a result?
	// It is always allowed when it is truncated, but also
	// if it is exact and the original binary mantissa is even
	// When disallowed, we can substract 1.
	uok := !uexact || !du0 || fracu > 0
	if uexact && du0 && fracu == 0 {
		uok = mant&1 == 0
	}
	if !uok {
		du--
	}
	// Is 'dc' the correctly rounded base 10 mantissa?
	// The correct rounding might be dc+1
	cup := false // don't round up.
	if cexact {
		// If we computed an exact product, the half integer
		// should round to next integer if 'dc' is odd.
		cup = fracc > 1<<(extra-1) ||
			(fracc == 1<<(extra-1) && !dc0) ||
			(fracc == 1<<(extra-1) && dc0 && dc&1 == 1)
	} else {
		// otherwise, the result is a lower truncation of the ideal
		// result.
		cup = fracc>>(extra-1) == 1
	}
	// Is 'dl' an allowed representation?
	// Only if it is an exact value, and if the original binary mantissa
	// was even.
	lok := lexact && dl0 && fracl == 0 && (mant&1 == 0)
	// We need to remember whether the trimmed digits of 'dc' are zero.
	c0 := cexact && dc0 && fracc == 0
	// render digits
	ryuDigits(d, dl, dc, du, lok, c0, cup)
	d.dp -= q
	return
}

// RyuFixed is a variation of the original Ryu algorithm for fixed precision
// output. It can output 16 digits reliably.
func RyuFixed(d *decimalSlice, mant uint64, exp int, prec int, flt *floatInfo) {
	if prec > 16 {
		panic("RyuFixed called with prec > 16")
	}
	// Fixed precision output proceeds as shortest output
	// except that we don't compute bounds, and that denormals
	// must be processed as normalized floats (and require higher
	// multipliers).
	//
	// In shortest mode, denormals are processed with their native
	// exponent and lose precision.
	if mant == 0 {
		d.nd, d.dp = 0, 0
		return
	}
	// substract mantbits to interpret mantissa as integer
	e2 := exp - int(flt.mantbits)
	// renormalize denormals to a 53-bit mantissa.
	if b := bits.Len64(mant); b < 53 {
		mant = mant << uint(53-b)
		e2 += int(b) - 53
	}
	// multiply by a power of 10. It is required to know
	// whether the computation is exact.
	var q int
	var pow *extfloat128 // a representation of 10^q
	if e2 < 0 {
		// Find 10^q *larger* than 2^-exp
		e := uint(-e2)
		q = int(Exp2toExponent10(e) + 1)
		pow = &RyuPowersOfTen[q]
	} else {
		// Divide by a power of 10 strictly less than 2^exp
		q = int(Exp2toExponent10(uint(e2)) - 1)
		if q < 0 {
			q = 0
		}
		pow = &RyuInvPowersOfTen[q]
		q = -q
	}
	// Is it an exact computation?
	exact := false
	switch {
	case q > 55:
		// large positive powers of ten are not exact
	case 54 >= q && q >= 0:
		exact = true
	case 0 > q && q >= -22:
		// division by a power of ten might be exact
		// if mantissas are multiples of 5
		if divisibleByPower5(mant, -q) {
			exact = true
		}
	default:
		// division by 10^23 cannot be exact
		// as 5^23 has 54 bits.
	}

	// Compute Floor(x*10^q)
	di, d0 := ryuMultiply(mant, pow.Hi, pow.Lo)
	dexp2 := e2 + pow.Exp + 64 + 55
	// If computation was an exact division, lower bits must be ignored.
	if q < 0 && exact {
		d0 = true
	}
	// Remove extra lower bits and keep rounding info.
	extra := uint(-dexp2)
	extramask := uint64(1<<extra - 1)

	di, dfrac := di>>extra, di&extramask
	roundUp := false
	if exact {
		// If we computed an exact product, d + 1/2
		// should round to d+1 if 'd' is odd.
		roundUp = dfrac > 1<<(extra-1) ||
			(dfrac == 1<<(extra-1) && !d0) ||
			(dfrac == 1<<(extra-1) && d0 && di&1 == 1)
	} else {
		// otherwise, d+1/2 always rounds up because
		// we truncated below.
		roundUp = dfrac>>(extra-1) == 1
	}
	if dfrac != 0 {
		d0 = false
	}
	// Proceed to the requested number of digits
	size := 16
	max := uint64pow10[prec]
	if di > max {
		// We reduce to exactly prec digits.
		size = prec
	} else {
		// d is larger than the original 53-bit mantissa
		// so has at least 16 digits. But prec <= 16,
		// so d has 16 digits exactly.
		size = 16
	}
	trimmed := 0
	for di >= max {
		a, b := di/10, di%10
		di = a
		trimmed++
		if b > 5 {
			roundUp = true
		} else if b < 5 {
			roundUp = false
		} else { // b == 5
			roundUp = !d0 || (d0 && di&1 == 1)
		}
		// update d0 for next iteration
		d0 = d0 && b == 0
	}
	if roundUp {
		di++
	}
	if di >= max {
		// Happens if di was originally 99999....xx
		di = di / 10
		trimmed++
	}
	// render digits
	n := uint(size)
	d.nd = int(size)
	for v := di; v > 0; {
		v1, v2 := v/10, v%10
		v = v1
		if int(n) == d.nd && v2 == 0 {
			d.nd-- // trim trailing zeros
			trimmed++
		}
		n--
		d.d[n&31] = byte(v2 + '0')
	}
	d.dp = d.nd + trimmed - q
	return
}

func RyuFromDecimal(mant uint64, exp int, flt *floatInfo) (fbits uint64, ovf, ok bool) {
	// Conversion from decimal to binary floating-point
	// can be achieved by reusing the same building blocks
	// as the Ryū algorithm.
	//
	// Given a decimal mantissa, we can multiply by the requested
	// power of ten using the same routines. The 64 bit result
	// is guaranteed to be correctly truncated (floored), when
	// the decimal mantissa fits in 55 bits.
	//
	// This covers 16-digit mantissas, and a few 17-digits values.

	const maxMantBitlen = 55 // could be larger?
	bitlen := bits.Len64(mant)
	var e2 int
	switch {
	case bitlen < maxMantBitlen:
		mant <<= uint(maxMantBitlen - bitlen)
		e2 = bitlen - maxMantBitlen
	case bitlen == maxMantBitlen:
		e2 = 0
	case bitlen > maxMantBitlen:
		zeros := bits.TrailingZeros64(mant)
		if zeros >= bitlen-maxMantBitlen {
			// Actually only maxMantBitlen significant bits
			mant = mant >> uint(bitlen-maxMantBitlen)
			e2 = bitlen - maxMantBitlen
		} else {
			return 0, false, false // cannot handle values that large.
		}
	}

	// multiply by a power of 10. It is required to know
	// whether the computation is exact.
	var pow *extfloat128 // a representation of 10^q
	switch {
	case exp > 309:
		return 0x7ff << 52, true, true
	case exp < -342:
		return 0, false, true
	case exp > 0:
		pow = &RyuPowersOfTen[exp]
	case exp == 0:
		// no multiply
	case exp < 0:
		pow = &RyuInvPowersOfTen[-exp]
	}
	// Is it an exact computation?
	exact := false
	switch {
	case exp > 55:
		// large positive powers of ten are not exact
	case 54 >= exp && exp >= 0:
		exact = true
	case 0 > exp && exp >= -25:
		// division by a power of ten might be exact
		// if mantissas are multiples of 5
		if divisibleByPower5(mant, -exp) {
			exact = true
		}
	default:
		// division by 10^25 cannot be exact
		// as 5^25 has 59 bits.
	}

	// Compute Floor(x*10^q)
	var di uint64
	var d0 bool
	if exp == 0 {
		di, d0 = mant, true
		exact = true
	} else {
		di, d0 = ryuMultiply(mant, pow.Hi, pow.Lo)
		e2 += pow.Exp + 64 + 55
	}
	// If computation was an exact division, lower bits must be ignored.
	if exp < 0 && exact {
		d0 = true
	}
	// Is exponent too low? Shrink mantissa for denormals.
	blen := bits.Len64(di)
	e2 += blen - 1
	extra := uint(blen - 53) // number of lower bits to remove
	if e2 < flt.bias+1 {
		extra += uint(flt.bias + 1 - e2)
		e2 = flt.bias + 1
	}
	if extra > uint(blen) {
		return 0.0, false, true
	}
	// Compute correct rounding.
	extramask := uint64(1<<extra - 1)
	di, dfrac := di>>extra, di&extramask
	roundUp := false
	if exact {
		// If we computed an exact product, d + 1/2
		// should round to d+1 if 'd' is odd.
		roundUp = dfrac > 1<<(extra-1) ||
			(dfrac == 1<<(extra-1) && !d0) ||
			(dfrac == 1<<(extra-1) && d0 && di&1 == 1)
	} else {
		// otherwise, d+1/2 always rounds up because
		// we truncated below.
		roundUp = dfrac>>(extra-1) == 1
	}
	if dfrac != 0 {
		d0 = false
	}
	if roundUp {
		di++
	}

	// Rounding might have added a bit; shift down.
	if di == 2<<flt.mantbits {
		di >>= 1
		e2++
	}

	// Infinities.
	if e2-flt.bias >= 1<<flt.expbits-1 {
		// ±Inf
		di = 0
		e2 = 1<<flt.expbits - 1 + flt.bias
		ovf = true
	} else if di&(1<<flt.mantbits) == 0 {
		// Denormalized?
		e2 = flt.bias
	}
	// Assemble bits.
	fbits = di & (uint64(1)<<flt.mantbits - 1)
	fbits |= uint64((e2-flt.bias)&(1<<flt.expbits-1)) << flt.mantbits
	return fbits, ovf, true
}

func divisibleByPower5(m uint64, k int) bool {
	for i := 0; i < k; i++ {
		a, b := m/5, m%5
		if b != 0 {
			return false
		}
		m = a
	}
	return true
}

func ryuDigits(d *decimalSlice, lower, central, upper uint64,
	lok, c0, cup bool) {
	lhi, llo := uint32(lower/1e9), uint32(lower%1e9)
	chi, clo := uint32(central/1e9), uint32(central%1e9)
	uhi, ulo := uint32(upper/1e9), uint32(upper%1e9)
	if uhi == 0 {
		// only low digits (for denormals)
		ryuDigits32(d, llo, clo, ulo,
			lok, c0, cup, 0)
	} else if lhi < uhi {
		// truncate 9 digits at once.
		lok = lok && llo == 0
		c0 = c0 && clo == 0
		cup = (clo > 5e8) || (clo == 5e8 && cup)
		d.nd = 0
		ryuDigits32(d, lhi, uint32(central/1e9), uhi,
			lok, c0, cup, 0)
		d.dp += 9
	} else {
		d.nd = 0
		// emit high part
		n := uint(16)
		for v := chi; v > 0; {
			v1, v2 := v/10, v%10
			v = v1
			n--
			d.d[n&31] = byte(v2 + '0')
		}
		copy(d.d[0:], d.d[n:16])
		d.nd = int(16 - n)
		// emit low part
		ryuDigits32(d, llo, clo, ulo,
			lok, c0, cup, d.nd+8)
	}
}

// ryuDigits32 emits decimal digits for a number less than 1e9.
func ryuDigits32(d *decimalSlice, lower, central, upper uint32,
	lok, c0, cup bool, endindex int) {
	trimmed := 0
	// Remember last trimmed digit to check for round-up.
	// c0 will be used to remember zeroness of following digits.
	cNextDigit := 0
	for {
		// Trim digits as long as it is possible.
		// Note that is it forbidden to go below 'lower'.
		l, ldigit := lower/10, lower%10
		c, cdigit := central/10, central%10
		u, _ := upper/10, upper%10
		lok = lok && ldigit == 0
		if !lok && l == u {
			break
		}
		// Check that we didn't cross the lower boundary.
		// The case where l == c < u is extremely rare,
		// and means that 'central' is very close but less than
		// an integer ending with many zeros, and usually
		// the "round-up" logic hides the problem.
		if l == c && !lok && c < u {
			c++
			cdigit = 0
			cup = false
		}
		trimmed++
		// Remember trimmed digits of c
		c0 = c0 && cNextDigit == 0
		cNextDigit = int(cdigit)
		lower, central, upper = l, c, u
	}
	// should we round up?
	if trimmed > 0 {
		cup = cNextDigit > 5 ||
			(cNextDigit == 5 && !c0) ||
			(cNextDigit == 5 && c0 && central&1 == 1)
	}
	if central < upper && cup {
		central++
	}
	if endindex > 0 {
		// We know where the number ends, fill directly
		endindex -= trimmed
		v := central
		for n := endindex; n >= d.nd; n-- {
			v1, v2 := v/10, v%10
			d.d[n] = byte(v2 + '0')
			v = v1
		}
		d.nd = endindex + 1
		d.dp = d.nd + trimmed
		return
	}
	// stupid
	n := uint(32)
	for v := central; v > 0; {
		v1, v2 := v/10, v%10
		n--
		d.d[n&31] = byte(v2 + '0')
		v = v1
	}
	copy(d.d[:], d.d[n:32])
	d.nd = int(32 - n)
	d.dp = d.nd + trimmed
}

// computeBounds returns a floating-point vector (l, c, u)×2^e2
// where the mantissas are 55-bit integers, describing the interval
// represented by the input float64.
func computeBounds(mant uint64, exp int) (lower, central, upper uint64, e2 int) {
	// substract mantbits to interpret mantissa as integer
	exp = exp - int(float64info.mantbits)
	expBiased := exp - float64info.bias

	if mant != 1<<float64info.mantbits || expBiased == 1 {
		// regular case (or denormals)
		lower, central, upper = 2*mant-1, 2*mant, 2*mant+1
		e2 = exp - 1
		return
	} else {
		// border of an exponent
		lower, central, upper = 4*mant-1, 4*mant, 4*mant+2
		e2 = exp - 2
		return
	}
}

// Exp2toExponent10 returns q = math.Floor(e * log10(2))
func Exp2toExponent10(e uint) uint {
	if e > 1600 {
		panic("out of approx range")
	}
	// log10(2) = 0.3010299956639812 = 78913.207... / 2**18
	return (e * 78913) >> 18
}

// ryuMultiply returns the 64 highest bits of the product:
// mant * (hi<<64|lo), where mant is a 55-bit integer.
// Also the boolean is set to true if the result is "exact",
// in the sense that all lower bits were zero.
func ryuMultiply(mant uint64, hi, lo uint64) (uint64, bool) {
	// long multiplication
	mant <<= 9
	l1, l0 := bits.Mul64(mant, lo)
	h1, h0 := bits.Mul64(mant, hi)
	mid, carry := bits.Add64(l1, h0, 0)
	return h1 + carry, mid == 0 && l0 == 0
}

type extfloat128 struct {
	Hi  uint64
	Lo  uint64
	Exp int
}

// ryuPowersOfTen[q] stores floating-point representations of 10^q,
// with 128-bit mantissas. The mantissa is always rounded down.
var RyuPowersOfTen = [...]extfloat128{
	{Hi: 0x8000000000000000, Lo: 0x0000000000000000, Exp: -127},
	{Hi: 0xa000000000000000, Lo: 0x0000000000000000, Exp: -124},
	{Hi: 0xc800000000000000, Lo: 0x0000000000000000, Exp: -121},
	{Hi: 0xfa00000000000000, Lo: 0x0000000000000000, Exp: -118},
	{Hi: 0x9c40000000000000, Lo: 0x0000000000000000, Exp: -114},
	{Hi: 0xc350000000000000, Lo: 0x0000000000000000, Exp: -111},
	{Hi: 0xf424000000000000, Lo: 0x0000000000000000, Exp: -108},
	{Hi: 0x9896800000000000, Lo: 0x0000000000000000, Exp: -104},
	{Hi: 0xbebc200000000000, Lo: 0x0000000000000000, Exp: -101},
	{Hi: 0xee6b280000000000, Lo: 0x0000000000000000, Exp: -98},
	{Hi: 0x9502f90000000000, Lo: 0x0000000000000000, Exp: -94},
	{Hi: 0xba43b74000000000, Lo: 0x0000000000000000, Exp: -91},
	{Hi: 0xe8d4a51000000000, Lo: 0x0000000000000000, Exp: -88},
	{Hi: 0x9184e72a00000000, Lo: 0x0000000000000000, Exp: -84},
	{Hi: 0xb5e620f480000000, Lo: 0x0000000000000000, Exp: -81},
	{Hi: 0xe35fa931a0000000, Lo: 0x0000000000000000, Exp: -78},
	{Hi: 0x8e1bc9bf04000000, Lo: 0x0000000000000000, Exp: -74},
	{Hi: 0xb1a2bc2ec5000000, Lo: 0x0000000000000000, Exp: -71},
	{Hi: 0xde0b6b3a76400000, Lo: 0x0000000000000000, Exp: -68},
	{Hi: 0x8ac7230489e80000, Lo: 0x0000000000000000, Exp: -64},
	{Hi: 0xad78ebc5ac620000, Lo: 0x0000000000000000, Exp: -61},
	{Hi: 0xd8d726b7177a8000, Lo: 0x0000000000000000, Exp: -58},
	{Hi: 0x878678326eac9000, Lo: 0x0000000000000000, Exp: -54},
	{Hi: 0xa968163f0a57b400, Lo: 0x0000000000000000, Exp: -51},
	{Hi: 0xd3c21bcecceda100, Lo: 0x0000000000000000, Exp: -48},
	{Hi: 0x84595161401484a0, Lo: 0x0000000000000000, Exp: -44},
	{Hi: 0xa56fa5b99019a5c8, Lo: 0x0000000000000000, Exp: -41},
	{Hi: 0xcecb8f27f4200f3a, Lo: 0x0000000000000000, Exp: -38},
	{Hi: 0x813f3978f8940984, Lo: 0x4000000000000000, Exp: -34},
	{Hi: 0xa18f07d736b90be5, Lo: 0x5000000000000000, Exp: -31},
	{Hi: 0xc9f2c9cd04674ede, Lo: 0xa400000000000000, Exp: -28},
	{Hi: 0xfc6f7c4045812296, Lo: 0x4d00000000000000, Exp: -25},
	{Hi: 0x9dc5ada82b70b59d, Lo: 0xf020000000000000, Exp: -21},
	{Hi: 0xc5371912364ce305, Lo: 0x6c28000000000000, Exp: -18},
	{Hi: 0xf684df56c3e01bc6, Lo: 0xc732000000000000, Exp: -15},
	{Hi: 0x9a130b963a6c115c, Lo: 0x3c7f400000000000, Exp: -11},
	{Hi: 0xc097ce7bc90715b3, Lo: 0x4b9f100000000000, Exp: -8},
	{Hi: 0xf0bdc21abb48db20, Lo: 0x1e86d40000000000, Exp: -5},
	{Hi: 0x96769950b50d88f4, Lo: 0x1314448000000000, Exp: -1},
	{Hi: 0xbc143fa4e250eb31, Lo: 0x17d955a000000000, Exp: 2},
	{Hi: 0xeb194f8e1ae525fd, Lo: 0x5dcfab0800000000, Exp: 5},
	{Hi: 0x92efd1b8d0cf37be, Lo: 0x5aa1cae500000000, Exp: 9},
	{Hi: 0xb7abc627050305ad, Lo: 0xf14a3d9e40000000, Exp: 12},
	{Hi: 0xe596b7b0c643c719, Lo: 0x6d9ccd05d0000000, Exp: 15},
	{Hi: 0x8f7e32ce7bea5c6f, Lo: 0xe4820023a2000000, Exp: 19},
	{Hi: 0xb35dbf821ae4f38b, Lo: 0xdda2802c8a800000, Exp: 22},
	{Hi: 0xe0352f62a19e306e, Lo: 0xd50b2037ad200000, Exp: 25},
	{Hi: 0x8c213d9da502de45, Lo: 0x4526f422cc340000, Exp: 29},
	{Hi: 0xaf298d050e4395d6, Lo: 0x9670b12b7f410000, Exp: 32},
	{Hi: 0xdaf3f04651d47b4c, Lo: 0x3c0cdd765f114000, Exp: 35},
	{Hi: 0x88d8762bf324cd0f, Lo: 0xa5880a69fb6ac800, Exp: 39},
	{Hi: 0xab0e93b6efee0053, Lo: 0x8eea0d047a457a00, Exp: 42},
	{Hi: 0xd5d238a4abe98068, Lo: 0x72a4904598d6d880, Exp: 45},
	{Hi: 0x85a36366eb71f041, Lo: 0x47a6da2b7f864750, Exp: 49},
	{Hi: 0xa70c3c40a64e6c51, Lo: 0x999090b65f67d924, Exp: 52},
	{Hi: 0xd0cf4b50cfe20765, Lo: 0xfff4b4e3f741cf6d, Exp: 55},
	{Hi: 0x82818f1281ed449f, Lo: 0xbff8f10e7a8921a4, Exp: 59},
	{Hi: 0xa321f2d7226895c7, Lo: 0xaff72d52192b6a0d, Exp: 62},
	{Hi: 0xcbea6f8ceb02bb39, Lo: 0x9bf4f8a69f764490, Exp: 65},
	{Hi: 0xfee50b7025c36a08, Lo: 0x02f236d04753d5b4, Exp: 68},
	{Hi: 0x9f4f2726179a2245, Lo: 0x01d762422c946590, Exp: 72},
	{Hi: 0xc722f0ef9d80aad6, Lo: 0x424d3ad2b7b97ef5, Exp: 75},
	{Hi: 0xf8ebad2b84e0d58b, Lo: 0xd2e0898765a7deb2, Exp: 78},
	{Hi: 0x9b934c3b330c8577, Lo: 0x63cc55f49f88eb2f, Exp: 82},
	{Hi: 0xc2781f49ffcfa6d5, Lo: 0x3cbf6b71c76b25fb, Exp: 85},
	{Hi: 0xf316271c7fc3908a, Lo: 0x8bef464e3945ef7a, Exp: 88},
	{Hi: 0x97edd871cfda3a56, Lo: 0x97758bf0e3cbb5ac, Exp: 92},
	{Hi: 0xbde94e8e43d0c8ec, Lo: 0x3d52eeed1cbea317, Exp: 95},
	{Hi: 0xed63a231d4c4fb27, Lo: 0x4ca7aaa863ee4bdd, Exp: 98},
	{Hi: 0x945e455f24fb1cf8, Lo: 0x8fe8caa93e74ef6a, Exp: 102},
	{Hi: 0xb975d6b6ee39e436, Lo: 0xb3e2fd538e122b44, Exp: 105},
	{Hi: 0xe7d34c64a9c85d44, Lo: 0x60dbbca87196b616, Exp: 108},
	{Hi: 0x90e40fbeea1d3a4a, Lo: 0xbc8955e946fe31cd, Exp: 112},
	{Hi: 0xb51d13aea4a488dd, Lo: 0x6babab6398bdbe41, Exp: 115},
	{Hi: 0xe264589a4dcdab14, Lo: 0xc696963c7eed2dd1, Exp: 118},
	{Hi: 0x8d7eb76070a08aec, Lo: 0xfc1e1de5cf543ca2, Exp: 122},
	{Hi: 0xb0de65388cc8ada8, Lo: 0x3b25a55f43294bcb, Exp: 125},
	{Hi: 0xdd15fe86affad912, Lo: 0x49ef0eb713f39ebe, Exp: 128},
	{Hi: 0x8a2dbf142dfcc7ab, Lo: 0x6e3569326c784337, Exp: 132},
	{Hi: 0xacb92ed9397bf996, Lo: 0x49c2c37f07965404, Exp: 135},
	{Hi: 0xd7e77a8f87daf7fb, Lo: 0xdc33745ec97be906, Exp: 138},
	{Hi: 0x86f0ac99b4e8dafd, Lo: 0x69a028bb3ded71a3, Exp: 142},
	{Hi: 0xa8acd7c0222311bc, Lo: 0xc40832ea0d68ce0c, Exp: 145},
	{Hi: 0xd2d80db02aabd62b, Lo: 0xf50a3fa490c30190, Exp: 148},
	{Hi: 0x83c7088e1aab65db, Lo: 0x792667c6da79e0fa, Exp: 152},
	{Hi: 0xa4b8cab1a1563f52, Lo: 0x577001b891185938, Exp: 155},
	{Hi: 0xcde6fd5e09abcf26, Lo: 0xed4c0226b55e6f86, Exp: 158},
	{Hi: 0x80b05e5ac60b6178, Lo: 0x544f8158315b05b4, Exp: 162},
	{Hi: 0xa0dc75f1778e39d6, Lo: 0x696361ae3db1c721, Exp: 165},
	{Hi: 0xc913936dd571c84c, Lo: 0x03bc3a19cd1e38e9, Exp: 168},
	{Hi: 0xfb5878494ace3a5f, Lo: 0x04ab48a04065c723, Exp: 171},
	{Hi: 0x9d174b2dcec0e47b, Lo: 0x62eb0d64283f9c76, Exp: 175},
	{Hi: 0xc45d1df942711d9a, Lo: 0x3ba5d0bd324f8394, Exp: 178},
	{Hi: 0xf5746577930d6500, Lo: 0xca8f44ec7ee36479, Exp: 181},
	{Hi: 0x9968bf6abbe85f20, Lo: 0x7e998b13cf4e1ecb, Exp: 185},
	{Hi: 0xbfc2ef456ae276e8, Lo: 0x9e3fedd8c321a67e, Exp: 188},
	{Hi: 0xefb3ab16c59b14a2, Lo: 0xc5cfe94ef3ea101e, Exp: 191},
	{Hi: 0x95d04aee3b80ece5, Lo: 0xbba1f1d158724a12, Exp: 195},
	{Hi: 0xbb445da9ca61281f, Lo: 0x2a8a6e45ae8edc97, Exp: 198},
	{Hi: 0xea1575143cf97226, Lo: 0xf52d09d71a3293bd, Exp: 201},
	{Hi: 0x924d692ca61be758, Lo: 0x593c2626705f9c56, Exp: 205},
	{Hi: 0xb6e0c377cfa2e12e, Lo: 0x6f8b2fb00c77836c, Exp: 208},
	{Hi: 0xe498f455c38b997a, Lo: 0x0b6dfb9c0f956447, Exp: 211},
	{Hi: 0x8edf98b59a373fec, Lo: 0x4724bd4189bd5eac, Exp: 215},
	{Hi: 0xb2977ee300c50fe7, Lo: 0x58edec91ec2cb657, Exp: 218},
	{Hi: 0xdf3d5e9bc0f653e1, Lo: 0x2f2967b66737e3ed, Exp: 221},
	{Hi: 0x8b865b215899f46c, Lo: 0xbd79e0d20082ee74, Exp: 225},
	{Hi: 0xae67f1e9aec07187, Lo: 0xecd8590680a3aa11, Exp: 228},
	{Hi: 0xda01ee641a708de9, Lo: 0xe80e6f4820cc9495, Exp: 231},
	{Hi: 0x884134fe908658b2, Lo: 0x3109058d147fdcdd, Exp: 235},
	{Hi: 0xaa51823e34a7eede, Lo: 0xbd4b46f0599fd415, Exp: 238},
	{Hi: 0xd4e5e2cdc1d1ea96, Lo: 0x6c9e18ac7007c91a, Exp: 241},
	{Hi: 0x850fadc09923329e, Lo: 0x03e2cf6bc604ddb0, Exp: 245},
	{Hi: 0xa6539930bf6bff45, Lo: 0x84db8346b786151c, Exp: 248},
	{Hi: 0xcfe87f7cef46ff16, Lo: 0xe612641865679a63, Exp: 251},
	{Hi: 0x81f14fae158c5f6e, Lo: 0x4fcb7e8f3f60c07e, Exp: 255},
	{Hi: 0xa26da3999aef7749, Lo: 0xe3be5e330f38f09d, Exp: 258},
	{Hi: 0xcb090c8001ab551c, Lo: 0x5cadf5bfd3072cc5, Exp: 261},
	{Hi: 0xfdcb4fa002162a63, Lo: 0x73d9732fc7c8f7f6, Exp: 264},
	{Hi: 0x9e9f11c4014dda7e, Lo: 0x2867e7fddcdd9afa, Exp: 268},
	{Hi: 0xc646d63501a1511d, Lo: 0xb281e1fd541501b8, Exp: 271},
	{Hi: 0xf7d88bc24209a565, Lo: 0x1f225a7ca91a4226, Exp: 274},
	{Hi: 0x9ae757596946075f, Lo: 0x3375788de9b06958, Exp: 278},
	{Hi: 0xc1a12d2fc3978937, Lo: 0x0052d6b1641c83ae, Exp: 281},
	{Hi: 0xf209787bb47d6b84, Lo: 0xc0678c5dbd23a49a, Exp: 284},
	{Hi: 0x9745eb4d50ce6332, Lo: 0xf840b7ba963646e0, Exp: 288},
	{Hi: 0xbd176620a501fbff, Lo: 0xb650e5a93bc3d898, Exp: 291},
	{Hi: 0xec5d3fa8ce427aff, Lo: 0xa3e51f138ab4cebe, Exp: 294},
	{Hi: 0x93ba47c980e98cdf, Lo: 0xc66f336c36b10137, Exp: 298},
	{Hi: 0xb8a8d9bbe123f017, Lo: 0xb80b0047445d4184, Exp: 301},
	{Hi: 0xe6d3102ad96cec1d, Lo: 0xa60dc059157491e5, Exp: 304},
	{Hi: 0x9043ea1ac7e41392, Lo: 0x87c89837ad68db2f, Exp: 308},
	{Hi: 0xb454e4a179dd1877, Lo: 0x29babe4598c311fb, Exp: 311},
	{Hi: 0xe16a1dc9d8545e94, Lo: 0xf4296dd6fef3d67a, Exp: 314},
	{Hi: 0x8ce2529e2734bb1d, Lo: 0x1899e4a65f58660c, Exp: 318},
	{Hi: 0xb01ae745b101e9e4, Lo: 0x5ec05dcff72e7f8f, Exp: 321},
	{Hi: 0xdc21a1171d42645d, Lo: 0x76707543f4fa1f73, Exp: 324},
	{Hi: 0x899504ae72497eba, Lo: 0x6a06494a791c53a8, Exp: 328},
	{Hi: 0xabfa45da0edbde69, Lo: 0x0487db9d17636892, Exp: 331},
	{Hi: 0xd6f8d7509292d603, Lo: 0x45a9d2845d3c42b6, Exp: 334},
	{Hi: 0x865b86925b9bc5c2, Lo: 0x0b8a2392ba45a9b2, Exp: 338},
	{Hi: 0xa7f26836f282b732, Lo: 0x8e6cac7768d7141e, Exp: 341},
	{Hi: 0xd1ef0244af2364ff, Lo: 0x3207d795430cd926, Exp: 344},
	{Hi: 0x8335616aed761f1f, Lo: 0x7f44e6bd49e807b8, Exp: 348},
	{Hi: 0xa402b9c5a8d3a6e7, Lo: 0x5f16206c9c6209a6, Exp: 351},
	{Hi: 0xcd036837130890a1, Lo: 0x36dba887c37a8c0f, Exp: 354},
	{Hi: 0x802221226be55a64, Lo: 0xc2494954da2c9789, Exp: 358},
	{Hi: 0xa02aa96b06deb0fd, Lo: 0xf2db9baa10b7bd6c, Exp: 361},
	{Hi: 0xc83553c5c8965d3d, Lo: 0x6f92829494e5acc7, Exp: 364},
	{Hi: 0xfa42a8b73abbf48c, Lo: 0xcb772339ba1f17f9, Exp: 367},
	{Hi: 0x9c69a97284b578d7, Lo: 0xff2a760414536efb, Exp: 371},
	{Hi: 0xc38413cf25e2d70d, Lo: 0xfef5138519684aba, Exp: 374},
	{Hi: 0xf46518c2ef5b8cd1, Lo: 0x7eb258665fc25d69, Exp: 377},
	{Hi: 0x98bf2f79d5993802, Lo: 0xef2f773ffbd97a61, Exp: 381},
	{Hi: 0xbeeefb584aff8603, Lo: 0xaafb550ffacfd8fa, Exp: 384},
	{Hi: 0xeeaaba2e5dbf6784, Lo: 0x95ba2a53f983cf38, Exp: 387},
	{Hi: 0x952ab45cfa97a0b2, Lo: 0xdd945a747bf26183, Exp: 391},
	{Hi: 0xba756174393d88df, Lo: 0x94f971119aeef9e4, Exp: 394},
	{Hi: 0xe912b9d1478ceb17, Lo: 0x7a37cd5601aab85d, Exp: 397},
	{Hi: 0x91abb422ccb812ee, Lo: 0xac62e055c10ab33a, Exp: 401},
	{Hi: 0xb616a12b7fe617aa, Lo: 0x577b986b314d6009, Exp: 404},
	{Hi: 0xe39c49765fdf9d94, Lo: 0xed5a7e85fda0b80b, Exp: 407},
	{Hi: 0x8e41ade9fbebc27d, Lo: 0x14588f13be847307, Exp: 411},
	{Hi: 0xb1d219647ae6b31c, Lo: 0x596eb2d8ae258fc8, Exp: 414},
	{Hi: 0xde469fbd99a05fe3, Lo: 0x6fca5f8ed9aef3bb, Exp: 417},
	{Hi: 0x8aec23d680043bee, Lo: 0x25de7bb9480d5854, Exp: 421},
	{Hi: 0xada72ccc20054ae9, Lo: 0xaf561aa79a10ae6a, Exp: 424},
	{Hi: 0xd910f7ff28069da4, Lo: 0x1b2ba1518094da04, Exp: 427},
	{Hi: 0x87aa9aff79042286, Lo: 0x90fb44d2f05d0842, Exp: 431},
	{Hi: 0xa99541bf57452b28, Lo: 0x353a1607ac744a53, Exp: 434},
	{Hi: 0xd3fa922f2d1675f2, Lo: 0x42889b8997915ce8, Exp: 437},
	{Hi: 0x847c9b5d7c2e09b7, Lo: 0x69956135febada11, Exp: 441},
	{Hi: 0xa59bc234db398c25, Lo: 0x43fab9837e699095, Exp: 444},
	{Hi: 0xcf02b2c21207ef2e, Lo: 0x94f967e45e03f4bb, Exp: 447},
	{Hi: 0x8161afb94b44f57d, Lo: 0x1d1be0eebac278f5, Exp: 451},
	{Hi: 0xa1ba1ba79e1632dc, Lo: 0x6462d92a69731732, Exp: 454},
	{Hi: 0xca28a291859bbf93, Lo: 0x7d7b8f7503cfdcfe, Exp: 457},
	{Hi: 0xfcb2cb35e702af78, Lo: 0x5cda735244c3d43e, Exp: 460},
	{Hi: 0x9defbf01b061adab, Lo: 0x3a0888136afa64a7, Exp: 464},
	{Hi: 0xc56baec21c7a1916, Lo: 0x088aaa1845b8fdd0, Exp: 467},
	{Hi: 0xf6c69a72a3989f5b, Lo: 0x8aad549e57273d45, Exp: 470},
	{Hi: 0x9a3c2087a63f6399, Lo: 0x36ac54e2f678864b, Exp: 474},
	{Hi: 0xc0cb28a98fcf3c7f, Lo: 0x84576a1bb416a7dd, Exp: 477},
	{Hi: 0xf0fdf2d3f3c30b9f, Lo: 0x656d44a2a11c51d5, Exp: 480},
	{Hi: 0x969eb7c47859e743, Lo: 0x9f644ae5a4b1b325, Exp: 484},
	{Hi: 0xbc4665b596706114, Lo: 0x873d5d9f0dde1fee, Exp: 487},
	{Hi: 0xeb57ff22fc0c7959, Lo: 0xa90cb506d155a7ea, Exp: 490},
	{Hi: 0x9316ff75dd87cbd8, Lo: 0x09a7f12442d588f2, Exp: 494},
	{Hi: 0xb7dcbf5354e9bece, Lo: 0x0c11ed6d538aeb2f, Exp: 497},
	{Hi: 0xe5d3ef282a242e81, Lo: 0x8f1668c8a86da5fa, Exp: 500},
	{Hi: 0x8fa475791a569d10, Lo: 0xf96e017d694487bc, Exp: 504},
	{Hi: 0xb38d92d760ec4455, Lo: 0x37c981dcc395a9ac, Exp: 507},
	{Hi: 0xe070f78d3927556a, Lo: 0x85bbe253f47b1417, Exp: 510},
	{Hi: 0x8c469ab843b89562, Lo: 0x93956d7478ccec8e, Exp: 514},
	{Hi: 0xaf58416654a6babb, Lo: 0x387ac8d1970027b2, Exp: 517},
	{Hi: 0xdb2e51bfe9d0696a, Lo: 0x06997b05fcc0319e, Exp: 520},
	{Hi: 0x88fcf317f22241e2, Lo: 0x441fece3bdf81f03, Exp: 524},
	{Hi: 0xab3c2fddeeaad25a, Lo: 0xd527e81cad7626c3, Exp: 527},
	{Hi: 0xd60b3bd56a5586f1, Lo: 0x8a71e223d8d3b074, Exp: 530},
	{Hi: 0x85c7056562757456, Lo: 0xf6872d5667844e49, Exp: 534},
	{Hi: 0xa738c6bebb12d16c, Lo: 0xb428f8ac016561db, Exp: 537},
	{Hi: 0xd106f86e69d785c7, Lo: 0xe13336d701beba52, Exp: 540},
	{Hi: 0x82a45b450226b39c, Lo: 0xecc0024661173473, Exp: 544},
	{Hi: 0xa34d721642b06084, Lo: 0x27f002d7f95d0190, Exp: 547},
	{Hi: 0xcc20ce9bd35c78a5, Lo: 0x31ec038df7b441f4, Exp: 550},
	{Hi: 0xff290242c83396ce, Lo: 0x7e67047175a15271, Exp: 553},
	{Hi: 0x9f79a169bd203e41, Lo: 0x0f0062c6e984d386, Exp: 557},
	{Hi: 0xc75809c42c684dd1, Lo: 0x52c07b78a3e60868, Exp: 560},
	{Hi: 0xf92e0c3537826145, Lo: 0xa7709a56ccdf8a82, Exp: 563},
	{Hi: 0x9bbcc7a142b17ccb, Lo: 0x88a66076400bb691, Exp: 567},
	{Hi: 0xc2abf989935ddbfe, Lo: 0x6acff893d00ea435, Exp: 570},
	{Hi: 0xf356f7ebf83552fe, Lo: 0x0583f6b8c4124d43, Exp: 573},
	{Hi: 0x98165af37b2153de, Lo: 0xc3727a337a8b704a, Exp: 577},
	{Hi: 0xbe1bf1b059e9a8d6, Lo: 0x744f18c0592e4c5c, Exp: 580},
	{Hi: 0xeda2ee1c7064130c, Lo: 0x1162def06f79df73, Exp: 583},
	{Hi: 0x9485d4d1c63e8be7, Lo: 0x8addcb5645ac2ba8, Exp: 587},
	{Hi: 0xb9a74a0637ce2ee1, Lo: 0x6d953e2bd7173692, Exp: 590},
	{Hi: 0xe8111c87c5c1ba99, Lo: 0xc8fa8db6ccdd0437, Exp: 593},
	{Hi: 0x910ab1d4db9914a0, Lo: 0x1d9c9892400a22a2, Exp: 597},
	{Hi: 0xb54d5e4a127f59c8, Lo: 0x2503beb6d00cab4b, Exp: 600},
	{Hi: 0xe2a0b5dc971f303a, Lo: 0x2e44ae64840fd61d, Exp: 603},
	{Hi: 0x8da471a9de737e24, Lo: 0x5ceaecfed289e5d2, Exp: 607},
	{Hi: 0xb10d8e1456105dad, Lo: 0x7425a83e872c5f47, Exp: 610},
	{Hi: 0xdd50f1996b947518, Lo: 0xd12f124e28f77719, Exp: 613},
	{Hi: 0x8a5296ffe33cc92f, Lo: 0x82bd6b70d99aaa6f, Exp: 617},
	{Hi: 0xace73cbfdc0bfb7b, Lo: 0x636cc64d1001550b, Exp: 620},
	{Hi: 0xd8210befd30efa5a, Lo: 0x3c47f7e05401aa4e, Exp: 623},
	{Hi: 0x8714a775e3e95c78, Lo: 0x65acfaec34810a71, Exp: 627},
	{Hi: 0xa8d9d1535ce3b396, Lo: 0x7f1839a741a14d0d, Exp: 630},
	{Hi: 0xd31045a8341ca07c, Lo: 0x1ede48111209a050, Exp: 633},
	{Hi: 0x83ea2b892091e44d, Lo: 0x934aed0aab460432, Exp: 637},
	{Hi: 0xa4e4b66b68b65d60, Lo: 0xf81da84d5617853f, Exp: 640},
	{Hi: 0xce1de40642e3f4b9, Lo: 0x36251260ab9d668e, Exp: 643},
	{Hi: 0x80d2ae83e9ce78f3, Lo: 0xc1d72b7c6b426019, Exp: 647},
	{Hi: 0xa1075a24e4421730, Lo: 0xb24cf65b8612f81f, Exp: 650},
	{Hi: 0xc94930ae1d529cfc, Lo: 0xdee033f26797b627, Exp: 653},
	{Hi: 0xfb9b7cd9a4a7443c, Lo: 0x169840ef017da3b1, Exp: 656},
	{Hi: 0x9d412e0806e88aa5, Lo: 0x8e1f289560ee864e, Exp: 660},
	{Hi: 0xc491798a08a2ad4e, Lo: 0xf1a6f2bab92a27e2, Exp: 663},
	{Hi: 0xf5b5d7ec8acb58a2, Lo: 0xae10af696774b1db, Exp: 666},
	{Hi: 0x9991a6f3d6bf1765, Lo: 0xacca6da1e0a8ef29, Exp: 670},
	{Hi: 0xbff610b0cc6edd3f, Lo: 0x17fd090a58d32af3, Exp: 673},
	{Hi: 0xeff394dcff8a948e, Lo: 0xddfc4b4cef07f5b0, Exp: 676},
	{Hi: 0x95f83d0a1fb69cd9, Lo: 0x4abdaf101564f98e, Exp: 680},
	{Hi: 0xbb764c4ca7a4440f, Lo: 0x9d6d1ad41abe37f1, Exp: 683},
	{Hi: 0xea53df5fd18d5513, Lo: 0x84c86189216dc5ed, Exp: 686},
	{Hi: 0x92746b9be2f8552c, Lo: 0x32fd3cf5b4e49bb4, Exp: 690},
	{Hi: 0xb7118682dbb66a77, Lo: 0x3fbc8c33221dc2a1, Exp: 693},
	{Hi: 0xe4d5e82392a40515, Lo: 0x0fabaf3feaa5334a, Exp: 696},
	{Hi: 0x8f05b1163ba6832d, Lo: 0x29cb4d87f2a7400e, Exp: 700},
	{Hi: 0xb2c71d5bca9023f8, Lo: 0x743e20e9ef511012, Exp: 703},
	{Hi: 0xdf78e4b2bd342cf6, Lo: 0x914da9246b255416, Exp: 706},
	{Hi: 0x8bab8eefb6409c1a, Lo: 0x1ad089b6c2f7548e, Exp: 710},
	{Hi: 0xae9672aba3d0c320, Lo: 0xa184ac2473b529b1, Exp: 713},
	{Hi: 0xda3c0f568cc4f3e8, Lo: 0xc9e5d72d90a2741e, Exp: 716},
	{Hi: 0x8865899617fb1871, Lo: 0x7e2fa67c7a658892, Exp: 720},
	{Hi: 0xaa7eebfb9df9de8d, Lo: 0xddbb901b98feeab7, Exp: 723},
	{Hi: 0xd51ea6fa85785631, Lo: 0x552a74227f3ea565, Exp: 726},
	{Hi: 0x8533285c936b35de, Lo: 0xd53a88958f87275f, Exp: 730},
	{Hi: 0xa67ff273b8460356, Lo: 0x8a892abaf368f137, Exp: 733},
	{Hi: 0xd01fef10a657842c, Lo: 0x2d2b7569b0432d85, Exp: 736},
	{Hi: 0x8213f56a67f6b29b, Lo: 0x9c3b29620e29fc73, Exp: 740},
	{Hi: 0xa298f2c501f45f42, Lo: 0x8349f3ba91b47b8f, Exp: 743},
	{Hi: 0xcb3f2f7642717713, Lo: 0x241c70a936219a73, Exp: 746},
	{Hi: 0xfe0efb53d30dd4d7, Lo: 0xed238cd383aa0110, Exp: 749},
	{Hi: 0x9ec95d1463e8a506, Lo: 0xf4363804324a40aa, Exp: 753},
	{Hi: 0xc67bb4597ce2ce48, Lo: 0xb143c6053edcd0d5, Exp: 756},
	{Hi: 0xf81aa16fdc1b81da, Lo: 0xdd94b7868e94050a, Exp: 759},
	{Hi: 0x9b10a4e5e9913128, Lo: 0xca7cf2b4191c8326, Exp: 763},
	{Hi: 0xc1d4ce1f63f57d72, Lo: 0xfd1c2f611f63a3f0, Exp: 766},
	{Hi: 0xf24a01a73cf2dccf, Lo: 0xbc633b39673c8cec, Exp: 769},
	{Hi: 0x976e41088617ca01, Lo: 0xd5be0503e085d813, Exp: 773},
	{Hi: 0xbd49d14aa79dbc82, Lo: 0x4b2d8644d8a74e18, Exp: 776},
	{Hi: 0xec9c459d51852ba2, Lo: 0xddf8e7d60ed1219e, Exp: 779},
	{Hi: 0x93e1ab8252f33b45, Lo: 0xcabb90e5c942b503, Exp: 783},
	{Hi: 0xb8da1662e7b00a17, Lo: 0x3d6a751f3b936243, Exp: 786},
	{Hi: 0xe7109bfba19c0c9d, Lo: 0x0cc512670a783ad4, Exp: 789},
	{Hi: 0x906a617d450187e2, Lo: 0x27fb2b80668b24c5, Exp: 793},
	{Hi: 0xb484f9dc9641e9da, Lo: 0xb1f9f660802dedf6, Exp: 796},
	{Hi: 0xe1a63853bbd26451, Lo: 0x5e7873f8a0396973, Exp: 799},
	{Hi: 0x8d07e33455637eb2, Lo: 0xdb0b487b6423e1e8, Exp: 803},
	{Hi: 0xb049dc016abc5e5f, Lo: 0x91ce1a9a3d2cda62, Exp: 806},
	{Hi: 0xdc5c5301c56b75f7, Lo: 0x7641a140cc7810fb, Exp: 809},
	{Hi: 0x89b9b3e11b6329ba, Lo: 0xa9e904c87fcb0a9d, Exp: 813},
	{Hi: 0xac2820d9623bf429, Lo: 0x546345fa9fbdcd44, Exp: 816},
	{Hi: 0xd732290fbacaf133, Lo: 0xa97c177947ad4095, Exp: 819},
	{Hi: 0x867f59a9d4bed6c0, Lo: 0x49ed8eabcccc485d, Exp: 823},
	{Hi: 0xa81f301449ee8c70, Lo: 0x5c68f256bfff5a74, Exp: 826},
	{Hi: 0xd226fc195c6a2f8c, Lo: 0x73832eec6fff3111, Exp: 829},
	{Hi: 0x83585d8fd9c25db7, Lo: 0xc831fd53c5ff7eab, Exp: 833},
	{Hi: 0xa42e74f3d032f525, Lo: 0xba3e7ca8b77f5e55, Exp: 836},
	{Hi: 0xcd3a1230c43fb26f, Lo: 0x28ce1bd2e55f35eb, Exp: 839},
	{Hi: 0x80444b5e7aa7cf85, Lo: 0x7980d163cf5b81b3, Exp: 843},
	{Hi: 0xa0555e361951c366, Lo: 0xd7e105bcc332621f, Exp: 846},
	{Hi: 0xc86ab5c39fa63440, Lo: 0x8dd9472bf3fefaa7, Exp: 849},
	{Hi: 0xfa856334878fc150, Lo: 0xb14f98f6f0feb951, Exp: 852},
	{Hi: 0x9c935e00d4b9d8d2, Lo: 0x6ed1bf9a569f33d3, Exp: 856},
	{Hi: 0xc3b8358109e84f07, Lo: 0x0a862f80ec4700c8, Exp: 859},
	{Hi: 0xf4a642e14c6262c8, Lo: 0xcd27bb612758c0fa, Exp: 862},
	{Hi: 0x98e7e9cccfbd7dbd, Lo: 0x8038d51cb897789c, Exp: 866},
	{Hi: 0xbf21e44003acdd2c, Lo: 0xe0470a63e6bd56c3, Exp: 869},
	{Hi: 0xeeea5d5004981478, Lo: 0x1858ccfce06cac74, Exp: 872},
	{Hi: 0x95527a5202df0ccb, Lo: 0x0f37801e0c43ebc8, Exp: 876},
	{Hi: 0xbaa718e68396cffd, Lo: 0xd30560258f54e6ba, Exp: 879},
	{Hi: 0xe950df20247c83fd, Lo: 0x47c6b82ef32a2069, Exp: 882},
	{Hi: 0x91d28b7416cdd27e, Lo: 0x4cdc331d57fa5441, Exp: 886},
	{Hi: 0xb6472e511c81471d, Lo: 0xe0133fe4adf8e952, Exp: 889},
	{Hi: 0xe3d8f9e563a198e5, Lo: 0x58180fddd97723a6, Exp: 892},
	{Hi: 0x8e679c2f5e44ff8f, Lo: 0x570f09eaa7ea7648, Exp: 896},
	{Hi: 0xb201833b35d63f73, Lo: 0x2cd2cc6551e513da, Exp: 899},
	{Hi: 0xde81e40a034bcf4f, Lo: 0xf8077f7ea65e58d1, Exp: 902},
	{Hi: 0x8b112e86420f6191, Lo: 0xfb04afaf27faf782, Exp: 906},
	{Hi: 0xadd57a27d29339f6, Lo: 0x79c5db9af1f9b563, Exp: 909},
	{Hi: 0xd94ad8b1c7380874, Lo: 0x18375281ae7822bc, Exp: 912},
	{Hi: 0x87cec76f1c830548, Lo: 0x8f2293910d0b15b5, Exp: 916},
	{Hi: 0xa9c2794ae3a3c69a, Lo: 0xb2eb3875504ddb22, Exp: 919},
	{Hi: 0xd433179d9c8cb841, Lo: 0x5fa60692a46151eb, Exp: 922},
	{Hi: 0x849feec281d7f328, Lo: 0xdbc7c41ba6bcd333, Exp: 926},
	{Hi: 0xa5c7ea73224deff3, Lo: 0x12b9b522906c0800, Exp: 929},
	{Hi: 0xcf39e50feae16bef, Lo: 0xd768226b34870a00, Exp: 932},
	{Hi: 0x81842f29f2cce375, Lo: 0xe6a1158300d46640, Exp: 936},
	{Hi: 0xa1e53af46f801c53, Lo: 0x60495ae3c1097fd0, Exp: 939},
	{Hi: 0xca5e89b18b602368, Lo: 0x385bb19cb14bdfc4, Exp: 942},
	{Hi: 0xfcf62c1dee382c42, Lo: 0x46729e03dd9ed7b5, Exp: 945},
	{Hi: 0x9e19db92b4e31ba9, Lo: 0x6c07a2c26a8346d1, Exp: 949},
	{Hi: 0xc5a05277621be293, Lo: 0xc7098b7305241885, Exp: 952},
	{Hi: 0xf70867153aa2db38, Lo: 0xb8cbee4fc66d1ea7, Exp: 955},
	{Hi: 0x9a65406d44a5c903, Lo: 0x737f74f1dc043328, Exp: 959},
	{Hi: 0xc0fe908895cf3b44, Lo: 0x505f522e53053ff2, Exp: 962},
	{Hi: 0xf13e34aabb430a15, Lo: 0x647726b9e7c68fef, Exp: 965},
	{Hi: 0x96c6e0eab509e64d, Lo: 0x5eca783430dc19f5, Exp: 969},
	{Hi: 0xbc789925624c5fe0, Lo: 0xb67d16413d132072, Exp: 972},
	{Hi: 0xeb96bf6ebadf77d8, Lo: 0xe41c5bd18c57e88f, Exp: 975},
	{Hi: 0x933e37a534cbaae7, Lo: 0x8e91b962f7b6f159, Exp: 979},
	{Hi: 0xb80dc58e81fe95a1, Lo: 0x723627bbb5a4adb0, Exp: 982},
	{Hi: 0xe61136f2227e3b09, Lo: 0xcec3b1aaa30dd91c, Exp: 985},
	{Hi: 0x8fcac257558ee4e6, Lo: 0x213a4f0aa5e8a7b1, Exp: 989},
	{Hi: 0xb3bd72ed2af29e1f, Lo: 0xa988e2cd4f62d19d, Exp: 992},
	{Hi: 0xe0accfa875af45a7, Lo: 0x93eb1b80a33b8605, Exp: 995},
	{Hi: 0x8c6c01c9498d8b88, Lo: 0xbc72f130660533c3, Exp: 999},
}

// ryuInvPowersOfTen[q] stores floating-point representations of 1/10^q,
// with 128-bit mantissas. The mantissa is always rounded up (ceil(m)),
// as in the traditional "divide by multiply high and shift right" method.
// This allows obtaining correct results when computing 10^q * (1/10^q).
var RyuInvPowersOfTen = [...]extfloat128{
	{Hi: 0x8000000000000000, Lo: 0x0000000000000000, Exp: -127},
	{Hi: 0xcccccccccccccccc, Lo: 0xcccccccccccccccd, Exp: -131},
	{Hi: 0xa3d70a3d70a3d70a, Lo: 0x3d70a3d70a3d70a4, Exp: -134},
	{Hi: 0x83126e978d4fdf3b, Lo: 0x645a1cac083126ea, Exp: -137},
	{Hi: 0xd1b71758e219652b, Lo: 0xd3c36113404ea4a9, Exp: -141},
	{Hi: 0xa7c5ac471b478423, Lo: 0x0fcf80dc33721d54, Exp: -144},
	{Hi: 0x8637bd05af6c69b5, Lo: 0xa63f9a49c2c1b110, Exp: -147},
	{Hi: 0xd6bf94d5e57a42bc, Lo: 0x3d32907604691b4d, Exp: -151},
	{Hi: 0xabcc77118461cefc, Lo: 0xfdc20d2b36ba7c3e, Exp: -154},
	{Hi: 0x89705f4136b4a597, Lo: 0x31680a88f8953031, Exp: -157},
	{Hi: 0xdbe6fecebdedd5be, Lo: 0xb573440e5a884d1c, Exp: -161},
	{Hi: 0xafebff0bcb24aafe, Lo: 0xf78f69a51539d749, Exp: -164},
	{Hi: 0x8cbccc096f5088cb, Lo: 0xf93f87b7442e45d4, Exp: -167},
	{Hi: 0xe12e13424bb40e13, Lo: 0x2865a5f206b06fba, Exp: -171},
	{Hi: 0xb424dc35095cd80f, Lo: 0x538484c19ef38c95, Exp: -174},
	{Hi: 0x901d7cf73ab0acd9, Lo: 0x0f9d37014bf60a11, Exp: -177},
	{Hi: 0xe69594bec44de15b, Lo: 0x4c2ebe687989a9b4, Exp: -181},
	{Hi: 0xb877aa3236a4b449, Lo: 0x09befeb9fad487c3, Exp: -184},
	{Hi: 0x9392ee8e921d5d07, Lo: 0x3aff322e62439fd0, Exp: -187},
	{Hi: 0xec1e4a7db69561a5, Lo: 0x2b31e9e3d06c32e6, Exp: -191},
	{Hi: 0xbce5086492111aea, Lo: 0x88f4bb1ca6bcf585, Exp: -194},
	{Hi: 0x971da05074da7bee, Lo: 0xd3f6fc16ebca5e04, Exp: -197},
	{Hi: 0xf1c90080baf72cb1, Lo: 0x5324c68b12dd6339, Exp: -201},
	{Hi: 0xc16d9a0095928a27, Lo: 0x75b7053c0f178294, Exp: -204},
	{Hi: 0x9abe14cd44753b52, Lo: 0xc4926a9672793543, Exp: -207},
	{Hi: 0xf79687aed3eec551, Lo: 0x3a83ddbd83f52205, Exp: -211},
	{Hi: 0xc612062576589dda, Lo: 0x95364afe032a819e, Exp: -214},
	{Hi: 0x9e74d1b791e07e48, Lo: 0x775ea264cf55347e, Exp: -217},
	{Hi: 0xfd87b5f28300ca0d, Lo: 0x8bca9d6e188853fd, Exp: -221},
	{Hi: 0xcad2f7f5359a3b3e, Lo: 0x096ee45813a04331, Exp: -224},
	{Hi: 0xa2425ff75e14fc31, Lo: 0xa1258379a94d028e, Exp: -227},
	{Hi: 0x81ceb32c4b43fcf4, Lo: 0x80eacf948770ced8, Exp: -230},
	{Hi: 0xcfb11ead453994ba, Lo: 0x67de18eda5814af3, Exp: -234},
	{Hi: 0xa6274bbdd0fadd61, Lo: 0xecb1ad8aeacdd58f, Exp: -237},
	{Hi: 0x84ec3c97da624ab4, Lo: 0xbd5af13bef0b113f, Exp: -240},
	{Hi: 0xd4ad2dbfc3d07787, Lo: 0x955e4ec64b44e865, Exp: -244},
	{Hi: 0xaa242499697392d2, Lo: 0xdde50bd1d5d0b9ea, Exp: -247},
	{Hi: 0x881cea14545c7575, Lo: 0x7e50d64177da2e55, Exp: -250},
	{Hi: 0xd9c7dced53c72255, Lo: 0x96e7bd358c904a22, Exp: -254},
	{Hi: 0xae397d8aa96c1b77, Lo: 0xabec975e0a0d081b, Exp: -257},
	{Hi: 0x8b61313bbabce2c6, Lo: 0x2323ac4b3b3da016, Exp: -260},
	{Hi: 0xdf01e85f912e37a3, Lo: 0x6b6c46dec52f6689, Exp: -264},
	{Hi: 0xb267ed1940f1c61c, Lo: 0x55f038b237591ed4, Exp: -267},
	{Hi: 0x8eb98a7a9a5b04e3, Lo: 0x77f3608e92adb243, Exp: -270},
	{Hi: 0xe45c10c42a2b3b05, Lo: 0x8cb89a7db77c506b, Exp: -274},
	{Hi: 0xb6b00d69bb55c8d1, Lo: 0x3d607b97c5fd0d23, Exp: -277},
	{Hi: 0x9226712162ab070d, Lo: 0xcab3961304ca70e9, Exp: -280},
	{Hi: 0xe9d71b689dde71af, Lo: 0xaab8f01e6e10b4a7, Exp: -284},
	{Hi: 0xbb127c53b17ec159, Lo: 0x5560c018580d5d53, Exp: -287},
	{Hi: 0x95a8637627989aad, Lo: 0xdde7001379a44aa9, Exp: -290},
	{Hi: 0xef73d256a5c0f77c, Lo: 0x963e66858f6d4441, Exp: -294},
	{Hi: 0xbf8fdb78849a5f96, Lo: 0xde98520472bdd034, Exp: -297},
	{Hi: 0x993fe2c6d07b7fab, Lo: 0xe546a8038efe402a, Exp: -300},
	{Hi: 0xf53304714d9265df, Lo: 0xd53dd99f4b3066a9, Exp: -304},
	{Hi: 0xc428d05aa4751e4c, Lo: 0xaa97e14c3c26b887, Exp: -307},
	{Hi: 0x9ced737bb6c4183d, Lo: 0x55464dd69685606c, Exp: -310},
	{Hi: 0xfb158592be068d2e, Lo: 0xeed6e2f0f0d56713, Exp: -314},
	{Hi: 0xc8de047564d20a8b, Lo: 0xf245825a5a445276, Exp: -317},
	{Hi: 0xa0b19d2ab70e6ed6, Lo: 0x5b6aceaeae9d0ec5, Exp: -320},
	{Hi: 0x808e17555f3ebf11, Lo: 0xe2bbd88bbee40bd1, Exp: -323},
	{Hi: 0xcdb02555653131b6, Lo: 0x3792f412cb06794e, Exp: -327},
	{Hi: 0xa48ceaaab75a8e2b, Lo: 0x5fa8c3423c052dd8, Exp: -330},
	{Hi: 0x83a3eeeef9153e89, Lo: 0x1953cf68300424ad, Exp: -333},
	{Hi: 0xd29fe4b18e88640e, Lo: 0x8eec7f0d19a03aae, Exp: -337},
	{Hi: 0xa87fea27a539e9a5, Lo: 0x3f2398d747b36225, Exp: -340},
	{Hi: 0x86ccbb52ea94baea, Lo: 0x98e947129fc2b4ea, Exp: -343},
	{Hi: 0xd7adf884aa879177, Lo: 0x5b0ed81dcc6abb10, Exp: -347},
	{Hi: 0xac8b2d36eed2dac5, Lo: 0xe272467e3d222f40, Exp: -350},
	{Hi: 0x8a08f0f8bf0f156b, Lo: 0x1b8e9ecb641b5900, Exp: -353},
	{Hi: 0xdcdb1b2798182244, Lo: 0xf8e431456cf88e66, Exp: -357},
	{Hi: 0xb0af48ec79ace837, Lo: 0x2d835a9df0c6d852, Exp: -360},
	{Hi: 0x8d590723948a535f, Lo: 0x579c487e5a38ad0f, Exp: -363},
	{Hi: 0xe2280b6c20dd5232, Lo: 0x25c6da63c38de1b1, Exp: -367},
	{Hi: 0xb4ecd5f01a4aa828, Lo: 0x1e38aeb6360b1af4, Exp: -370},
	{Hi: 0x90bd77f3483bb9b9, Lo: 0xb1c6f22b5e6f48c3, Exp: -373},
	{Hi: 0xe7958cb87392c2c2, Lo: 0xb60b1d1230b20e05, Exp: -377},
	{Hi: 0xb94470938fa89bce, Lo: 0xf808e40e8d5b3e6a, Exp: -380},
	{Hi: 0x9436c0760c86e30b, Lo: 0xf9a0b6720aaf6522, Exp: -383},
	{Hi: 0xed246723473e3813, Lo: 0x290123e9aab23b69, Exp: -387},
	{Hi: 0xbdb6b8e905cb600f, Lo: 0x5400e987bbc1c921, Exp: -390},
	{Hi: 0x97c560ba6b0919a5, Lo: 0xdccd879fc967d41b, Exp: -393},
	{Hi: 0xf2d56790ab41c2a2, Lo: 0xfae27299423fb9c4, Exp: -397},
	{Hi: 0xc24452da229b021b, Lo: 0xfbe85badce996169, Exp: -400},
	{Hi: 0x9b69dbe1b548ce7c, Lo: 0xc986afbe3ee11abb, Exp: -403},
	{Hi: 0xf8a95fcf88747d94, Lo: 0x75a44c6397ce912b, Exp: -407},
	{Hi: 0xc6ede63fa05d3143, Lo: 0x91503d1c79720dbc, Exp: -410},
	{Hi: 0x9f24b832e6b0f436, Lo: 0x0dd9ca7d2df4d7ca, Exp: -413},
	{Hi: 0xfea126b7d78186bc, Lo: 0xe2f610c84987bfa9, Exp: -417},
	{Hi: 0xcbb41ef979346bca, Lo: 0x4f2b40a03ad2ffba, Exp: -420},
	{Hi: 0xa2f67f2dfa90563b, Lo: 0x728900802f0f32fb, Exp: -423},
	{Hi: 0x825ecc24c873782f, Lo: 0x8ed400668c0c28c9, Exp: -426},
	{Hi: 0xd097ad07a71f26b2, Lo: 0x7e2000a41346a7a8, Exp: -430},
	{Hi: 0xa6dfbd9fb8e5b88e, Lo: 0xcb4ccd500f6bb953, Exp: -433},
	{Hi: 0x857fcae62d8493a5, Lo: 0x6f70a4400c562ddc, Exp: -436},
	{Hi: 0xd59944a37c0752a2, Lo: 0x4be76d3346f04960, Exp: -440},
	{Hi: 0xaae103b5fcd2a881, Lo: 0xd652bdc29f26a11a, Exp: -443},
	{Hi: 0x88b402f7fd75539b, Lo: 0x11dbcb0218ebb415, Exp: -446},
	{Hi: 0xdab99e59958885c4, Lo: 0xe95fab368e45ecee, Exp: -450},
	{Hi: 0xaefae51477a06b03, Lo: 0xede622920b6b23f2, Exp: -453},
	{Hi: 0x8bfbea76c619ef36, Lo: 0x57eb4edb3c55b65b, Exp: -456},
	{Hi: 0xdff9772470297ebd, Lo: 0x59787e2b93bc56f8, Exp: -460},
	{Hi: 0xb32df8e9f3546564, Lo: 0x47939822dc96abfa, Exp: -463},
	{Hi: 0x8f57fa54c2a9eab6, Lo: 0x9fa946824a12232e, Exp: -466},
	{Hi: 0xe55990879ddcaabd, Lo: 0xcc420a6a101d0516, Exp: -470},
	{Hi: 0xb77ada0617e3bbcb, Lo: 0x09ce6ebb40173745, Exp: -473},
	{Hi: 0x92c8ae6b464fc96f, Lo: 0x3b0b8bc90012929e, Exp: -476},
	{Hi: 0xeadab0aba3b2dbe5, Lo: 0x2b45ac74ccea842f, Exp: -480},
	{Hi: 0xbbe226efb628afea, Lo: 0x890489f70a55368c, Exp: -483},
	{Hi: 0x964e858c91ba2655, Lo: 0x3a6a07f8d510f870, Exp: -486},
	{Hi: 0xf07da27a82c37088, Lo: 0x5d767327bb4e5a4d, Exp: -490},
	{Hi: 0xc06481fb9bcf8d39, Lo: 0xe45ec2862f71e1d7, Exp: -493},
	{Hi: 0x99ea0196163fa42e, Lo: 0x504bced1bf8e4e46, Exp: -496},
	{Hi: 0xf64335bcf065d37d, Lo: 0x4d4617b5ff4a16d6, Exp: -500},
	{Hi: 0xc5029163f384a931, Lo: 0x0a9e795e65d4df12, Exp: -503},
	{Hi: 0x9d9ba7832936edc0, Lo: 0xd54b944b84aa4c0e, Exp: -506},
	{Hi: 0xfc2c3f3841f17c67, Lo: 0xbbac2078d443ace3, Exp: -510},
	{Hi: 0xc9bcff6034c13052, Lo: 0xfc89b393dd02f0b6, Exp: -513},
	{Hi: 0xa163ff802a3426a8, Lo: 0xca07c2dcb0cf26f8, Exp: -516},
	{Hi: 0x811ccc668829b887, Lo: 0x0806357d5a3f5260, Exp: -519},
	{Hi: 0xce947a3da6a9273e, Lo: 0x733d226229feea33, Exp: -523},
	{Hi: 0xa54394fe1eedb8fe, Lo: 0xc2974eb4ee658829, Exp: -526},
	{Hi: 0x843610cb4bf160cb, Lo: 0xcedf722a585139bb, Exp: -529},
	{Hi: 0xd389b47879823479, Lo: 0x4aff1d108d4ec2c4, Exp: -533},
	{Hi: 0xa93af6c6c79b5d2d, Lo: 0xd598e40d3dd89bd0, Exp: -536},
	{Hi: 0x87625f056c7c4a8b, Lo: 0x11471cd764ad4973, Exp: -539},
	{Hi: 0xd89d64d57a607744, Lo: 0xe871c7bf077ba8b8, Exp: -543},
	{Hi: 0xad4ab7112eb3929d, Lo: 0x86c16c98d2c953c7, Exp: -546},
	{Hi: 0x8aa22c0dbef60ee4, Lo: 0x6bcdf07a423aa96c, Exp: -549},
	{Hi: 0xddd0467c64bce4a0, Lo: 0xac7cb3f6d05ddbdf, Exp: -553},
	{Hi: 0xb1736b96b6fd83b3, Lo: 0xbd308ff8a6b17cb3, Exp: -556},
	{Hi: 0x8df5efabc5979c8f, Lo: 0xca8d3ffa1ef463c2, Exp: -559},
	{Hi: 0xe3231912d5bf60e6, Lo: 0x10e1fff697ed6c6a, Exp: -563},
	{Hi: 0xb5b5ada8aaff80b8, Lo: 0x0d819992132456bb, Exp: -566},
	{Hi: 0x915e2486ef32cd60, Lo: 0x0ace1474dc1d122f, Exp: -569},
	{Hi: 0xe896a0d7e51e1566, Lo: 0x77b020baf9c81d18, Exp: -573},
	{Hi: 0xba121a4650e4ddeb, Lo: 0x92f34d62616ce414, Exp: -576},
	{Hi: 0x94db483840b717ef, Lo: 0xa8c2a44eb4571cdd, Exp: -579},
	{Hi: 0xee2ba6c0678b597f, Lo: 0x746aa07ded582e2d, Exp: -583},
	{Hi: 0xbe89523386091465, Lo: 0xf6bbb397f1135824, Exp: -586},
	{Hi: 0x986ddb5c6b3a76b7, Lo: 0xf89629465a75e01d, Exp: -589},
	{Hi: 0xf3e2f893dec3f126, Lo: 0x5a89dba3c3efccfb, Exp: -593},
	{Hi: 0xc31bfa0fe5698db8, Lo: 0x486e494fcff30a63, Exp: -596},
	{Hi: 0x9c1661a651213e2d, Lo: 0x06bea10ca65c084f, Exp: -599},
	{Hi: 0xf9bd690a1b68637b, Lo: 0x3dfdce7aa3c673b1, Exp: -603},
	{Hi: 0xc7caba6e7c5382c8, Lo: 0xfe64a52ee96b8fc1, Exp: -606},
	{Hi: 0x9fd561f1fd0f9bd3, Lo: 0xfeb6ea8bedefa634, Exp: -609},
	{Hi: 0xffbbcfe994e5c61f, Lo: 0xfdf17746497f7053, Exp: -613},
	{Hi: 0xcc963fee10b7d1b3, Lo: 0x318df905079926a9, Exp: -616},
	{Hi: 0xa3ab66580d5fdaf5, Lo: 0xc13e60d0d2e0ebbb, Exp: -619},
	{Hi: 0x82ef85133de648c4, Lo: 0x9a984d73dbe722fc, Exp: -622},
	{Hi: 0xd17f3b51fca3a7a0, Lo: 0xf75a15862ca504c6, Exp: -626},
	{Hi: 0xa798fc4196e952e7, Lo: 0x2c48113823b73705, Exp: -629},
	{Hi: 0x8613fd0145877585, Lo: 0xbd06742ce95f5f37, Exp: -632},
	{Hi: 0xd686619ba27255a2, Lo: 0xc80a537b0efefebe, Exp: -636},
	{Hi: 0xab9eb47c81f5114f, Lo: 0x066ea92f3f326565, Exp: -639},
	{Hi: 0x894bc396ce5da772, Lo: 0x6b8bba8c328eb784, Exp: -642},
	{Hi: 0xdbac6c247d62a583, Lo: 0xdf45f746b74abf3a, Exp: -646},
	{Hi: 0xafbd2350644eeacf, Lo: 0xe5d1929ef90898fb, Exp: -649},
	{Hi: 0x8c974f7383725573, Lo: 0x1e414218c73a13fc, Exp: -652},
	{Hi: 0xe0f218b8d25088b8, Lo: 0x306869c13ec3532d, Exp: -656},
	{Hi: 0xb3f4e093db73a093, Lo: 0x59ed216765690f57, Exp: -659},
	{Hi: 0x8ff71a0fe2c2e6dc, Lo: 0x47f0e785eaba72ac, Exp: -662},
	{Hi: 0xe65829b3046b0afa, Lo: 0x0cb4a5a3112a5113, Exp: -666},
	{Hi: 0xb84687c269ef3bfb, Lo: 0x3d5d514f40eea743, Exp: -669},
	{Hi: 0x936b9fcebb25c995, Lo: 0xcab10dd900beec35, Exp: -672},
	{Hi: 0xebdf661791d60f56, Lo: 0x111b495b3464ad22, Exp: -676},
	{Hi: 0xbcb2b812db11a5de, Lo: 0x7415d448f6b6f0e8, Exp: -679},
	{Hi: 0x96f5600f15a7b7e5, Lo: 0x29ab103a5ef8c0ba, Exp: -682},
	{Hi: 0xf18899b1bc3f8ca1, Lo: 0xdc44e6c3cb279ac2, Exp: -686},
	{Hi: 0xc13a148e3032d6e7, Lo: 0xe36a52363c1faf02, Exp: -689},
	{Hi: 0x9a94dd3e8cf578b9, Lo: 0x82bb74f8301958cf, Exp: -692},
	{Hi: 0xf7549530e188c128, Lo: 0xd12bee59e68ef47d, Exp: -696},
	{Hi: 0xc5dd44271ad3cdba, Lo: 0x40eff1e1853f29fe, Exp: -699},
	{Hi: 0x9e4a9cec15763e2e, Lo: 0x9a598e4e043287ff, Exp: -702},
	{Hi: 0xfd442e4688bd304a, Lo: 0x908f4a166d1da664, Exp: -706},
	{Hi: 0xca9cf1d206fdc03b, Lo: 0xa6d90811f0e4851d, Exp: -709},
	{Hi: 0xa21727db38cb002f, Lo: 0xb8ada00e5a506a7d, Exp: -712},
	{Hi: 0x81ac1fe293d599bf, Lo: 0xc6f14cd848405531, Exp: -715},
	{Hi: 0xcf79cc9db955c2cc, Lo: 0x7182148d4066eeb5, Exp: -719},
	{Hi: 0xa5fb0a17c777cf09, Lo: 0xf468107100525891, Exp: -722},
	{Hi: 0x84c8d4dfd2c63f3b, Lo: 0x29ecd9f40041e074, Exp: -725},
	{Hi: 0xd47487cc8470652b, Lo: 0x7647c32000696720, Exp: -729},
	{Hi: 0xa9f6d30a038d1dbc, Lo: 0x5e9fcf4ccd211f4d, Exp: -732},
	{Hi: 0x87f8a8d4cfa417c9, Lo: 0xe54ca5d70a80e5d7, Exp: -735},
	{Hi: 0xd98ddaee19068c76, Lo: 0x3badd624dd9b0958, Exp: -739},
	{Hi: 0xae0b158b4738705e, Lo: 0x9624ab50b148d446, Exp: -742},
	{Hi: 0x8b3c113c38f9f37e, Lo: 0xde83bc408dd3dd05, Exp: -745},
	{Hi: 0xdec681f9f4c31f31, Lo: 0x6405fa00e2ec94d5, Exp: -749},
	{Hi: 0xb23867fb2a35b28d, Lo: 0xe99e619a4f23aa44, Exp: -752},
	{Hi: 0x8e938662882af53e, Lo: 0x547eb47b7282ee9d, Exp: -755},
	{Hi: 0xe41f3d6a7377eeca, Lo: 0x20caba5f1d9e4a94, Exp: -759},
	{Hi: 0xb67f6455292cbf08, Lo: 0x1a3bc84c17b1d543, Exp: -762},
	{Hi: 0x91ff83775423cc06, Lo: 0x7b6306a34627ddd0, Exp: -765},
	{Hi: 0xe998d258869facd7, Lo: 0x2bd1a438703fc94c, Exp: -769},
	{Hi: 0xbae0a846d2195712, Lo: 0x8974836059cca10a, Exp: -772},
	{Hi: 0x9580869f0e7aac0e, Lo: 0xd45d35e6ae3d4da1, Exp: -775},
	{Hi: 0xef340a98172aace4, Lo: 0x86fb897116c87c35, Exp: -779},
	{Hi: 0xbf5cd54678eef0b6, Lo: 0xd262d45a78a0635e, Exp: -782},
	{Hi: 0x991711052d8bf3c5, Lo: 0x751bdd152d4d1c4b, Exp: -785},
	{Hi: 0xf4f1b4d515acb93b, Lo: 0xee92fb5515482d45, Exp: -789},
	{Hi: 0xc3f490aa77bd60fc, Lo: 0xbedbfc4411068a9d, Exp: -792},
	{Hi: 0x9cc3a6eec6311a63, Lo: 0xcbe3303674053bb1, Exp: -795},
	{Hi: 0xfad2a4b13d1b5d6c, Lo: 0x796b805720085f82, Exp: -799},
	{Hi: 0xc8a883c0fdaf7df0, Lo: 0x6122cd128006b2ce, Exp: -802},
	{Hi: 0xa086cfcd97bf97f3, Lo: 0x80e8a40eccd228a5, Exp: -805},
	{Hi: 0x806bd9714632dff6, Lo: 0x00ba1cd8a3db53b7, Exp: -808},
	{Hi: 0xcd795be870516656, Lo: 0x67902e276c921f8c, Exp: -812},
	{Hi: 0xa46116538d0deb78, Lo: 0x52d9be85f074e609, Exp: -815},
	{Hi: 0x8380dea93da4bc60, Lo: 0x4247cb9e59f71e6e, Exp: -818},
	{Hi: 0xd267caa862a12d66, Lo: 0xd072df63c324fd7c, Exp: -822},
	{Hi: 0xa8530886b54dbdeb, Lo: 0xd9f57f830283fdfd, Exp: -825},
	{Hi: 0x86a8d39ef77164bc, Lo: 0xae5dff9c02033198, Exp: -828},
	{Hi: 0xd77485cb25823ac7, Lo: 0x7d633293366b828c, Exp: -832},
	{Hi: 0xac5d37d5b79b6239, Lo: 0x311c2875c522ced6, Exp: -835},
	{Hi: 0x89e42caaf9491b60, Lo: 0xf41686c49db57245, Exp: -838},
	{Hi: 0xdca04777f541c567, Lo: 0xecf0d7a0fc5583a1, Exp: -842},
	{Hi: 0xb080392cc4349dec, Lo: 0xbd8d794d96aacfb4, Exp: -845},
	{Hi: 0x8d3360f09cf6e4bd, Lo: 0x64712dd7abbbd95d, Exp: -848},
	{Hi: 0xe1ebce4dc7f16dfb, Lo: 0xd3e8495912c62895, Exp: -852},
	{Hi: 0xb4bca50b065abe63, Lo: 0x0fed077a756b53aa, Exp: -855},
	{Hi: 0x9096ea6f3848984f, Lo: 0x3ff0d2c85def7622, Exp: -858},
	{Hi: 0xe757dd7ec07426e5, Lo: 0x331aeada2fe589d0, Exp: -862},
	{Hi: 0xb913179899f68584, Lo: 0x28e2557b59846e40, Exp: -865},
	{Hi: 0x940f4613ae5ed136, Lo: 0x871b7795e136be9a, Exp: -868},
	{Hi: 0xece53cec4a314ebd, Lo: 0xa4f8bf5635246429, Exp: -872},
	{Hi: 0xbd8430bd08277231, Lo: 0x50c6ff782a838354, Exp: -875},
	{Hi: 0x979cf3ca6cec5b5a, Lo: 0xa705992ceecf9c43, Exp: -878},
	{Hi: 0xf294b943e17a2bc4, Lo: 0x3e6f5b7b17b2939e, Exp: -882},
	{Hi: 0xc21094364dfb5636, Lo: 0x985915fc12f542e5, Exp: -885},
	{Hi: 0x9b407691d7fc44f8, Lo: 0x79e0de63425dcf1e, Exp: -888},
	{Hi: 0xf867241c8cc6d4c0, Lo: 0xc30163d203c94b63, Exp: -892},
	{Hi: 0xc6b8e9b0709f109a, Lo: 0x359ab6419ca1091c, Exp: -895},
	{Hi: 0x9efa548d26e5a6e1, Lo: 0xc47bc5014a1a6db0, Exp: -898},
	{Hi: 0xfe5d54150b090b02, Lo: 0xd3f93b35435d7c4d, Exp: -902},
	{Hi: 0xcb7ddcdda26da268, Lo: 0xa9942f5dcf7dfd0a, Exp: -905},
	{Hi: 0xa2cb1717b52481ed, Lo: 0x54768c4b0c64ca6f, Exp: -908},
	{Hi: 0x823c12795db6ce57, Lo: 0x76c53d08d6b70859, Exp: -911},
	{Hi: 0xd0601d8efc57b08b, Lo: 0xf13b94daf124da27, Exp: -915},
	{Hi: 0xa6b34ad8c9dfc06f, Lo: 0xf42faa48c0ea481f, Exp: -918},
	{Hi: 0x855c3be0a17fcd26, Lo: 0x5cf2eea09a550680, Exp: -921},
	{Hi: 0xd5605fcdcf32e1d6, Lo: 0xfb1e4a9a90880a65, Exp: -925},
	{Hi: 0xaab37fd7d8f58178, Lo: 0xc8e5087ba6d33b84, Exp: -928},
	{Hi: 0x888f99797a5e012d, Lo: 0x6d8406c952429604, Exp: -931},
	{Hi: 0xda7f5bf590966848, Lo: 0xaf39a475506a899f, Exp: -935},
	{Hi: 0xaecc49914078536d, Lo: 0x58fae9f773886e19, Exp: -938},
	{Hi: 0x8bd6a141006042bd, Lo: 0xe0c8bb2c5c6d24e1, Exp: -941},
	{Hi: 0xdfbdcece67006ac9, Lo: 0x67a791e093e1d49b, Exp: -945},
	{Hi: 0xb2fe3f0b8599ef07, Lo: 0x861fa7e6dcb4aa16, Exp: -948},
	{Hi: 0x8f31cc0937ae58d2, Lo: 0xd1b2ecb8b0908811, Exp: -951},
	{Hi: 0xe51c79a85916f484, Lo: 0x82b7e12780e7401b, Exp: -955},
	{Hi: 0xb749faed14125d36, Lo: 0xcef980ec671f667c, Exp: -958},
	{Hi: 0x92a1958a7675175f, Lo: 0x0bfacd89ec191eca, Exp: -961},
	{Hi: 0xea9c227723ee8bcb, Lo: 0x465e15a979c1cadd, Exp: -965},
	{Hi: 0xbbb01b9283253ca2, Lo: 0x9eb1aaedfb016f17, Exp: -968},
	{Hi: 0x96267c7535b763b5, Lo: 0x4bc1558b2f3458df, Exp: -971},
	{Hi: 0xf03d93eebc589f88, Lo: 0x793555ab7eba27cb, Exp: -975},
	{Hi: 0xc0314325637a1939, Lo: 0xfa911155fefb5309, Exp: -978},
	{Hi: 0x99c102844f94e0fb, Lo: 0x2eda7444cbfc426e, Exp: -981},
	{Hi: 0xf6019da07f549b2b, Lo: 0x7e2a53a146606a49, Exp: -985},
	{Hi: 0xc4ce17b399107c22, Lo: 0xcb550fb4384d21d4, Exp: -988},
	{Hi: 0x9d71ac8fada6c9b5, Lo: 0x6f773fc3603db4aa, Exp: -991},
	{Hi: 0xfbe9141915d7a922, Lo: 0x4bf1ff9f0062baa9, Exp: -995},
	{Hi: 0xc987434744ac874e, Lo: 0xa327ffb266b56221, Exp: -998},
	{Hi: 0xa139029f6a239f72, Lo: 0x1c1fffc1ebc44e81, Exp: -1001},
	{Hi: 0x80fa687f881c7f8e, Lo: 0x7ce66634bc9d0b9a, Exp: -1004},
	{Hi: 0xce5d73ff402d98e3, Lo: 0xfb0a3d212dc81290, Exp: -1008},
	{Hi: 0xa5178fff668ae0b6, Lo: 0x626e974dbe39a873, Exp: -1011},
	{Hi: 0x8412d9991ed58091, Lo: 0xe858790afe9486c3, Exp: -1014},
	{Hi: 0xd3515c2831559a83, Lo: 0x0d5a5b44ca873e04, Exp: -1018},
	{Hi: 0xa90de3535aaae202, Lo: 0x711515d0a205cb37, Exp: -1021},
	{Hi: 0x873e4f75e2224e68, Lo: 0x5a7744a6e804a292, Exp: -1024},
	{Hi: 0xd863b256369d4a40, Lo: 0x90bed43e40076a83, Exp: -1028},
	{Hi: 0xad1c8eab5ee43b66, Lo: 0xda3243650005eed0, Exp: -1031},
	{Hi: 0x8a7d3eef7f1cfc52, Lo: 0x482835ea666b2573, Exp: -1034},
	{Hi: 0xdd95317f31c7fa1d, Lo: 0x40405643d711d584, Exp: -1038},
	{Hi: 0xb1442798f49ffb4a, Lo: 0x99cd11cfdf41779d, Exp: -1041},
	{Hi: 0x8dd01fad907ffc3b, Lo: 0xae3da7d97f6792e4, Exp: -1044},
	{Hi: 0xe2e69915b3fff9f9, Lo: 0x16c90c8f323f516d, Exp: -1048},
	{Hi: 0xb58547448ffffb2d, Lo: 0xabd40a0c2832a78b, Exp: -1051},
	{Hi: 0x91376c36d99995be, Lo: 0x23100809b9c21fa2, Exp: -1054},
	{Hi: 0xe858ad248f5c22c9, Lo: 0xd1b3400f8f9cff69, Exp: -1058},
	{Hi: 0xb9e08a83a5e34f07, Lo: 0xdaf5ccd93fb0cc54, Exp: -1061},
	{Hi: 0x94b3a202eb1c3f39, Lo: 0x7bf7d71432f3d6aa, Exp: -1064},
	{Hi: 0xedec366b11c6cb8f, Lo: 0x2cbfbe86b7ec8aa9, Exp: -1068},
	{Hi: 0xbe5691ef416bd60c, Lo: 0x23cc986bc656d554, Exp: -1071},
	{Hi: 0x9845418c345644d6, Lo: 0x830a13896b78aaaa, Exp: -1074},
	{Hi: 0xf3a20279ed56d48a, Lo: 0x6b43527578c11110, Exp: -1078},
	{Hi: 0xc2e801fb244576d5, Lo: 0x229c41f793cda740, Exp: -1081},
	{Hi: 0x9becce62836ac577, Lo: 0x4ee367f9430aec33, Exp: -1084},
	{Hi: 0xf97ae3d0d2446f25, Lo: 0x4b0573286b44ad1e, Exp: -1088},
	{Hi: 0xc795830d75038c1d, Lo: 0xd59df5b9ef6a2418, Exp: -1091},
	{Hi: 0x9faacf3df73609b1, Lo: 0x77b191618c54e9ad, Exp: -1094},
	{Hi: 0xff77b1fcbebcdc4f, Lo: 0x25e8e89c13bb0f7b, Exp: -1098},
	{Hi: 0xcc5fc196fefd7d0c, Lo: 0x1e53ed49a96272c9, Exp: -1101},
	{Hi: 0xa37fce126597973c, Lo: 0xe50ff107bab528a1, Exp: -1104},
	{Hi: 0x82cca4db847945ca, Lo: 0x50d98d9fc890ed4e, Exp: -1107},
	{Hi: 0xd1476e2c07286faa, Lo: 0x1af5af660db4aee2, Exp: -1111},
	{Hi: 0xa76c582338ed2621, Lo: 0xaf2af2b80af6f24f, Exp: -1114},
	{Hi: 0x85f0468293f0eb4e, Lo: 0x25bbf56008c58ea6, Exp: -1117},
	{Hi: 0xd64d3d9db981787d, Lo: 0x092cbbccdad5b109, Exp: -1121},
	{Hi: 0xab70fe17c79ac6ca, Lo: 0x6dbd630a48aaf407, Exp: -1124},
	{Hi: 0x892731ac9faf056e, Lo: 0xbe311c083a225cd3, Exp: -1127},
	{Hi: 0xdb71e91432b1a24a, Lo: 0xc9e82cd9f69d6151, Exp: -1131},
	{Hi: 0xaf8e5410288e1b6f, Lo: 0x07ecf0ae5ee44dda, Exp: -1134},
	{Hi: 0x8c71dcd9ba0b4925, Lo: 0x9ff0c08b7f1d0b15, Exp: -1137},
	{Hi: 0xe0b62e2929aba83c, Lo: 0x331acdabfe94de88, Exp: -1141},
	{Hi: 0xb3c4f1ba87bc8696, Lo: 0x8f48a4899877186d, Exp: -1144},
	{Hi: 0x8fd0c16206306bab, Lo: 0xa5d3b6d479f8e057, Exp: -1147},
	{Hi: 0xe61acf033d1a45df, Lo: 0x6fb92487298e33be, Exp: -1151},
	{Hi: 0xb8157268fdae9e4c, Lo: 0x5960ea05bad82965, Exp: -1154},
	{Hi: 0x93445b8731587ea3, Lo: 0x7ab3ee6afbe0211e, Exp: -1157},
	{Hi: 0xeba09271e88d976b, Lo: 0xf7864a44c633682f, Exp: -1161},
	{Hi: 0xbc807527ed3e12bc, Lo: 0xc605083704f5ecf3, Exp: -1164},
	{Hi: 0x96cd2a865764dbca, Lo: 0x380406926a5e5729, Exp: -1167},
	{Hi: 0xf148440a256e2c76, Lo: 0xc00670ea43ca250e, Exp: -1171},
	{Hi: 0xc1069cd4eabe89f8, Lo: 0x999ec0bb696e840b, Exp: -1174},
	{Hi: 0x9a6bb0aa55653b2d, Lo: 0x47b233c92125366f, Exp: -1177},
	{Hi: 0xf712b443bbd52b7b, Lo: 0xa5e9ec7501d523e5, Exp: -1181},
	{Hi: 0xc5a890362fddbc62, Lo: 0xeb2189f734aa831e, Exp: -1184},
	{Hi: 0x9e20735e8cb16382, Lo: 0x55b46e5f5d5535b1, Exp: -1187},
	{Hi: 0xfd00b897478238d0, Lo: 0x8920b098955522b5, Exp: -1191},
	{Hi: 0xca66fa129f9b60a6, Lo: 0xd41a26e077774ef7, Exp: -1194},
	{Hi: 0xa1ebfb4219491a1f, Lo: 0x1014ebe6c5f90bf9, Exp: -1197},
	{Hi: 0x818995ce7aa0e1b2, Lo: 0x7343efebd1940994, Exp: -1200},
	{Hi: 0xcf42894a5dce35ea, Lo: 0x52064cac828675ba, Exp: -1204},
	{Hi: 0xa5ced43b7e3e9188, Lo: 0x419ea3bd35385e2e, Exp: -1207},
	{Hi: 0x84a57695fe98746d, Lo: 0x014bb630f7604b58, Exp: -1210},
	{Hi: 0xd43bf0effdc0ba48, Lo: 0x0212bd1b2566def3, Exp: -1214},
	{Hi: 0xa9c98d8ccb009506, Lo: 0x680efdaf511f18c3, Exp: -1217},
	{Hi: 0x87d4713d6f33aa6b, Lo: 0x8672648c40e5ad69, Exp: -1220},
	{Hi: 0xd953e8624b85dd78, Lo: 0xd71d6dad34a2af0e, Exp: -1224},
	{Hi: 0xaddcb9e83c6b1793, Lo: 0xdf4abe242a1bbf3e, Exp: -1227},
	{Hi: 0x8b16fb203055ac76, Lo: 0x4c3bcb5021afcc32, Exp: -1230},
	{Hi: 0xde8b2b66b3bc4723, Lo: 0xad2c788035e61383, Exp: -1234},
	{Hi: 0xb208ef855c969f4f, Lo: 0xbdbd2d335e51a936, Exp: -1237},
	{Hi: 0x8e6d8c6ab0787f72, Lo: 0xfe30f0f5e50e20f8, Exp: -1240},
	{Hi: 0xe3e27a444d8d98b7, Lo: 0xfd1b1b2308169b26, Exp: -1244},
	{Hi: 0xb64ec836a47146f9, Lo: 0x9748e2826cdee285, Exp: -1247},
	{Hi: 0x91d8a02bb6c10594, Lo: 0x79071b9b8a4be86a, Exp: -1250},
	{Hi: 0xe95a99df8ace6f53, Lo: 0xf4d82c2c107973dd, Exp: -1254},
	{Hi: 0xbaaee17fa23ebf76, Lo: 0x5d79bcf00d2df64a, Exp: -1257},
	{Hi: 0x9558b4661b6565f8, Lo: 0x4ac7ca59a424c508, Exp: -1260},
	{Hi: 0xeef453d6923bd65a, Lo: 0x113faa2906a13b40, Exp: -1264},
	{Hi: 0xbf29dcaba82fdeae, Lo: 0x7432ee873880fc34, Exp: -1267},
}
