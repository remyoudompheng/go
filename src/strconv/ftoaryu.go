// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package strconv

import (
	"math/bits"
)

// binary to decimal conversion using the Ryū algorithm.
//
// See Ulf Adams, "Ryū: Fast Float-to-String Conversion" (doi:10.1145/3192366.3192369)
//
// Fixed precision formatting is a variant of the original paper's
// algorithm, where a single multiplication by 10^k is required,
// sharing the same rounding guarantees.

func ryuFtoaFixed32(d *decimalSlice, mant uint32, exp int, prec int) {
	if prec > 9 {
		panic("ryuFtoaFixed64 called with prec > 9")
	}
	// Zero input.
	if mant == 0 {
		d.nd, d.dp = 0, 0
		return
	}
	// Renormalize to a 25-bit mantissa.
	e2 := exp
	if b := bits.Len32(mant); b < 25 {
		mant = mant << uint(25-b)
		e2 += int(b) - 25
	}
	// Choose an exponent such that mant*2^e2*10^q has at least prec digits,
	// i.e
	//     2^(e2+54) >= 10^(-q+prec-1)
	// or q = -log10pow2(e2+54) + prec - 1
	//
	// The maximal required exponent is log10pow2(1074)+18 = 342
	q := -log10pow2(e2+24) + prec - 1

	// Now compute mant*2^e2*10^q.
	// Is it an exact computation?
	// Only small positive powers of 10 are exact (5^28 has 66 bits).
	exact := q <= 27 && q >= 0

	di, dexp2, d0 := mult64bitPow10(mant, e2, q)
	if dexp2 >= 0 {
		panic("not enough significant bits after mult64bitPow10")
	}
	// As a special case, computation might still be exact, if exponent
	// was negative and if it amounts to computing an exact division.
	// In that case, we ignore all lower bits.
	// Note that division by 10^11 cannot be exact as 5^11 has 26 bits.
	if q < 0 && q >= -10 && divisibleByPower5(uint64(mant), -q) {
		exact = true
		d0 = true
	}
	// Remove extra lower bits and keep rounding info.
	extra := uint(-dexp2)
	extramask := uint32(1<<extra - 1)

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
	max := uint32(uint64pow10[prec])
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
	// render digits (similar to formatBits)
	n := uint(prec)
	d.nd = int(prec)
	v := di
	for v >= 100 {
		v1, v2 := v/100, v%100
		n -= 2
		d.d[n+1] = smallsString[2*v2+1]
		d.d[n+0] = smallsString[2*v2+0]
		v = v1
	}
	if v > 0 {
		n--
		d.d[n] = smallsString[2*v+1]
	}
	if v >= 10 {
		n--
		d.d[n] = smallsString[2*v]
	}
	for d.d[d.nd-1] == '0' {
		d.nd--
		trimmed++
	}
	d.dp = d.nd + trimmed - q
	return
}

// ryuFtoaFixed64 formats mant*2^exp with prec decimal digits.
func ryuFtoaFixed64(d *decimalSlice, mant uint64, exp int, prec int) {
	if prec > 18 {
		panic("ryuFtoaFixed64 called with prec > 18")
	}
	// Zero input.
	if mant == 0 {
		d.nd, d.dp = 0, 0
		return
	}
	// Renormalize to a 55-bit mantissa.
	e2 := exp
	if b := bits.Len64(mant); b < 55 {
		mant = mant << uint(55-b)
		e2 += int(b) - 55
	}
	// Choose an exponent such that mant*2^e2*10^q has at least prec digits,
	// i.e
	//     2^(e2+54) >= 10^(-q+prec-1)
	// or q = -log10pow2(e2+54) + prec - 1
	//
	// The maximal required exponent is log10pow2(1074)+18 = 342
	q := -log10pow2(e2+54) + prec - 1

	// Now compute mant*2^e2*10^q.
	// Is it an exact computation?
	// Only small positive powers of 10 are exact (5^55 has 128 bits).
	exact := q <= 55 && q >= 0

	di, dexp2, d0 := mult128bitPow10(mant, e2, q)
	if dexp2 >= 0 {
		panic("not enough significant bits after mult128bitPow10")
	}
	// As a special case, computation might still be exact, if exponent
	// was negative and if it amounts to computing an exact division.
	// In that case, we ignore all lower bits.
	// Note that division by 10^23 cannot be exact as 5^23 has 54 bits.
	if q < 0 && q >= -22 && divisibleByPower5(mant, -q) {
		exact = true
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
	// render digits (similar to formatBits)
	n := uint(prec)
	d.nd = int(prec)
	v := di
	for v >= 100 {
		v1, v2 := v/100, v%100
		n -= 2
		d.d[n+1] = smallsString[2*v2+1]
		d.d[n+0] = smallsString[2*v2+0]
		v = v1
	}
	if v > 0 {
		n--
		d.d[n] = smallsString[2*v+1]
	}
	if v >= 10 {
		n--
		d.d[n] = smallsString[2*v]
	}
	for d.d[d.nd-1] == '0' {
		d.nd--
		trimmed++
	}
	d.dp = d.nd + trimmed - q
	return
}

// log10pow2 returns q = math.Floor(e * log10(2)).
// It expects an input exponent such that |e| < 1600.
func log10pow2(e int) int {
	// log10(2) = 0.3010299956639812 = 78913.207... / 2**18
	return (e * 78913) >> 18 // even for negative e
}

// mult64bitPow10 takes a floating-point input with a 25-bit
// mantissa and multiplies it with 10**e10. The resulting mantissa
// is m*P >> 57 where P is a 64-bit element of the pow10wide tables.
// It which is typically 31 or 32-bit wide.
// The returned boolean is true is all trimmed bits were zero.
func mult64bitPow10(m uint32, e2 int, e10 int) (uint32, int, bool) {
	if e10 == 0 {
		return m << 7, e2 - 7, true
	}
	var pow uint64
	if e10 >= 0 {
		pow = pow10wide[e10][0]
	} else {
		pow = invpow10wide[-e10][0] + 1
	}
	hi, lo := bits.Mul64(uint64(m), pow)
	e2 += pow10exp(e10) - 63 + 57
	return uint32(hi<<7 | lo>>57), e2, lo<<7 == 0
}

// mult128bitPow10 takes a floating-point input with a 55-bit
// mantissa and multiplies it with 10**e10. The resulting mantissa
// is m*P >> 119 where P is a 128-bit element of the pow10wide tables.
// It is typically 63 or 64-bit wide.
// The returned boolean is true is all trimmed bits were zero.
func mult128bitPow10(m uint64, e2 int, e10 int) (uint64, int, bool) {
	if e10 == 0 {
		return m << 9, e2 - 9, true
	}
	var pow [2]uint64
	if e10 >= 0 {
		pow = pow10wide[e10]
	} else {
		pow = invpow10wide[-e10]
	}
	e2 += pow10exp(e10) - 127 + 119

	// long multiplication
	l1, l0 := bits.Mul64(m, pow[1])
	h1, h0 := bits.Mul64(m, pow[0])
	mid, carry := bits.Add64(l1, h0, 0)
	h1 += carry
	return h1<<9 | mid>>55, e2, mid<<9 == 0 && l0 == 0
}
