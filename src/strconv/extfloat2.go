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

func ryuShortest(d *decimalSlice, mant uint64, exp int) {
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
	var pow [2]uint64 // a representation of 10^q
	var powexp int
	if e2 < 0 {
		// Find 10^q *larger* than 2^-e2
		q = int(exp2toExponent10(-e2) + 1)
	} else {
		// Divide by a power of 10 strictly less than 2^e2
		q = int(exp2toExponent10(e2) - 1)
		if q < 0 {
			q = 0
		}
		q = -q
	}
	if q >= 0 {
		pow = pow10wide[q]
	} else {
		pow = invpow10wide[-q]
	}
	powexp = exp10toExponent2(q) - 127
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
	dl, dl0 := ryuMultiply(ml, pow[0], pow[1])
	dc, dc0 := ryuMultiply(mc, pow[0], pow[1])
	du, du0 := ryuMultiply(mu, pow[0], pow[1])
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
	extra := uint(-(e2 + powexp + ryuMultiplyShift))
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

// ryuFixed is a variation of the original Ryu algorithm for fixed precision
// output. It can output 18 digits reliably.
func ryuFixed(d *decimalSlice, mant uint64, exp int, prec int, flt *floatInfo) {
	if prec > 18 {
		panic("ryuFixed called with prec > 18")
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
	// renormalize denormals (and float32s) to a 53-bit mantissa.
	if b := bits.Len64(mant); b < 53 {
		mant = mant << uint(53-b)
		e2 += int(b) - 53
	}
	// multiply by a power of 10, such that mant*2^e2*10^q
	// is guaranteed to be larger than 10^(prec-1), i.e
	//     2^(mantbits+e2) >= 10^(-q+prec-1)
	// or q = -exp2toExponent10(mantbits+e2) + prec - 1
	//
	// The maximal required exponent is exp2toExponent10(1074)+18 = 342
	q := -exp2toExponent10(e2+52) + prec - 1
	var pow [2]uint64 // a representation of 10^q
	var powexp int
	if q >= 0 {
		pow = pow10wide[q]
	} else {
		pow = invpow10wide[-q]
	}
	powexp = exp10toExponent2(q) - 127
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

	// Compute Floor(x*10^q). Optimal use of ryuMultiply
	// is for 55-bit inputs, so we shift mant again.
	mant <<= 2
	e2 -= 2
	// Since 10^(prec-1) <= x*10^q < 2*10^prec, it may be required
	// to divide by 10. The result of ryuMultiply is
	// 63 or 64 bits large, so we can guarantee only 18 digits.
	di, d0 := ryuMultiply(mant, pow[0], pow[1])
	dexp2 := e2 + powexp + ryuMultiplyShift
	if dexp2 >= 0 {
		println(mant, "x", pow[0], pow[1], "=", di, dexp2)
		panic("not enough significant bits after ryuMultiply")
	}
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
	max := uint64pow10[prec]
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
	n := uint(prec)
	d.nd = int(prec)
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

func ryuFromDecimal(mant uint64, exp int, flt *floatInfo) (fbits uint64, ovf bool) {
	// Conversion from decimal to binary floating-point
	// can be achieved by reusing the same building blocks
	// as the Ryū algorithm.
	//
	// Given a decimal mantissa, we can multiply by the requested
	// power of ten using the same routines. Then for any 64-bit
	// multiplier, the result is correctly truncated by
	// right-shift by 135 bits.
	// (TestRyuMultiplyCarry shows no error for 62-bit multipliers,
	// and no warning for 64-bit multipliers).

	bitlen := bits.Len64(mant)
	e2 := 0
	if bitlen < 64 {
		mant <<= uint(64 - bitlen)
		e2 = bitlen - 64
	}

	// multiply by a power of 10. It is required to know
	// whether the computation is exact.
	var pow [2]uint64 // a representation of 10^q
	var powexp int
	switch {
	case exp > 309:
		return 0x7ff << 52, true
	case exp < -342:
		return 0, false
	case exp > 0:
		pow = pow10wide[exp]
		powexp = exp10toExponent2(exp) - 127
	case exp == 0:
		// no multiply
	case exp < 0:
		pow = invpow10wide[-exp]
		powexp = exp10toExponent2(exp) - 127
	}
	// Is it an exact computation?
	exact := false
	switch {
	case exp > 55:
		// large positive powers of ten are not exact
	case 54 >= exp && exp >= 0:
		exact = true
	case 0 > exp && exp >= -27:
		// division by a power of ten might be exact
		// if mantissas are multiples of 5
		if divisibleByPower5(mant, -exp) {
			exact = true
		}
	default:
		// division by 10^28 cannot be exact
		// as 5^28 has 66 bits.
	}

	// Compute Floor(x*10^q). ryuMultiply2 returns
	// 54 bits, enough to find a float64 mantissa.
	var di uint64
	var d0 bool
	if exp == 0 {
		di, d0 = mant, true
		exact = true
	} else {
		di, d0 = ryuMultiply2(mant, pow[0], pow[1])
		e2 += powexp + ryuMultiply2Shift
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
		return 0.0, false
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
	return fbits, ovf
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
	// trim trailing zeros
	for d.nd > 0 && d.d[d.nd-1] == '0' {
		d.nd--
	}
}

// ryuDigits32 emits decimal digits for a number less than 1e9.
func ryuDigits32(d *decimalSlice, lower, central, upper uint32,
	lok, c0, cup bool, endindex int) {
	if upper == 0 {
		if endindex > 0 {
			d.dp = endindex
		} else {
			d.nd, d.dp = 0, 0
		}
		return
	}
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
			// don't trim the last digit as it is forbidden to go below l
			// other, trim and exit now.
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
		if l == u {
			break
		}
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

// exp2toExponent10 returns q = math.Floor(e * log10(2))
func exp2toExponent10(e int) int {
	if e > 1600 || e < -1600 {
		panic("out of approx range " + Itoa(int(e)))
	}
	// log10(2) = 0.3010299956639812 = 78913.207... / 2**18
	return (e * 78913) >> 18 // even for negative e
}

// ryuMultiply returns the 63 or 64 highest bits of the product:
// mant * (hi<<64|lo), where mant is a 55-bit integer,
// using a right shift by 119 bits.
// Also the boolean is set to true if the result is "exact",
// in the sense that all lower bits were zero.
func ryuMultiply(mant uint64, hi, lo uint64) (uint64, bool) {
	// long multiplication
	l1, l0 := bits.Mul64(mant, lo)
	h1, h0 := bits.Mul64(mant, hi)
	mid, carry := bits.Add64(l1, h0, 0)
	h1 += carry
	return h1<<9 | mid>>55, mid<<9 == 0 && l0 == 0
}

const ryuMultiplyShift = 119

// ryuMultiply2 is like ryuMultiply but takes a 64-bit multiplier
// and returns at least 54 bits from the full product which is
// 191 or 192 bit wide.
// It computes mant*(hi<<64|lo) >> 135.
func ryuMultiply2(mant uint64, hi, lo uint64) (uint64, bool) {
	// long multiplication
	l1, l0 := bits.Mul64(mant, lo)
	h1, h0 := bits.Mul64(mant, hi)
	mid, carry := bits.Add64(l1, h0, 0)
	h1 += carry
	return h1 >> 7, h1&(1<<7-1) == 0 && mid == 0 && l0 == 0
}

const ryuMultiply2Shift = 135

type extfloat128 struct {
	Hi  uint64
	Lo  uint64
	Exp int
}
