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
		a, b := divmod10(di)
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
		di, _ = divmod10(di)
		trimmed++
	}
	// render digits (similar to formatBits)
	n := uint(prec)
	d.nd = int(prec)
	v := di
	for v >= 100 {
		v1, v2 := divmod100(v)
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

// ryuFtoaShortest formats mant*2^exp with prec decimal digits.
func ryuFtoaShortest(d *decimalSlice, mant uint64, exp int, flt *floatInfo) {
	if mant == 0 {
		d.nd, d.dp = 0, 0
		return
	}
	// If input is an exact integer with fewer bits than the mantissa,
	// the previous and next integer are not admissible representations.
	if exp <= 0 && bits.TrailingZeros64(mant) >= -exp {
		mant >>= uint(-exp)
		ryuDigits(d, mant, mant, mant, true, false)
		return
	}
	ml, mc, mu, e2 := computeBounds(mant, exp, flt)
	if e2 == 0 {
		ryuDigits(d, ml, mc, mu, true, false)
		return
	}
	// Find 10^q *larger* than 2^-e2
	q := log10pow2(-e2) + 1

	// We are going to multiply by 10^q using 128-bit arithmetic.
	// The exponent is the same for all 3 numbers.
	var dl, dc, du uint64
	var dl0, dc0, du0 bool
	if flt == &float32info {
		var dl32, dc32, du32 uint32
		dl32, _, dl0 = mult64bitPow10(uint32(ml), e2, q)
		dc32, _, dc0 = mult64bitPow10(uint32(mc), e2, q)
		du32, e2, du0 = mult64bitPow10(uint32(mu), e2, q)
		dl, dc, du = uint64(dl32), uint64(dc32), uint64(du32)
	} else {
		dl, _, dl0 = mult128bitPow10(ml, e2, q)
		dc, _, dc0 = mult128bitPow10(mc, e2, q)
		du, e2, du0 = mult128bitPow10(mu, e2, q)
	}
	if e2 >= 0 {
		panic("not enough significant bits after mult128bitPow10")
	}
	// Is it an exact computation?
	if q > 55 {
		// Large positive powers of ten are not exact
		dl0, dc0, du0 = false, false, false
	}
	if q < 0 && q >= -24 {
		// Division by a power of ten may be exact.
		// (note that 5^25 is a 59-bit number so division by 5^25 is never exact).
		if divisibleByPower5(ml, -q) {
			dl0 = true
		}
		if divisibleByPower5(mc, -q) {
			dc0 = true
		}
		if divisibleByPower5(mu, -q) {
			du0 = true
		}
	}
	// Express the results (dl, dc, du)*2^e2 as integers.
	// Extra bits must be removed and rounding hints computed.
	extra := uint(-e2)
	extramask := uint64(1<<extra - 1)
	// Now compute the floored, integral base 10 mantissas.
	dl, fracl := dl>>extra, dl&extramask
	dc, fracc := dc>>extra, dc&extramask
	du, fracu := du>>extra, du&extramask
	// Is it allowed to use 'du' as a result?
	// It is always allowed when it is truncated, but also
	// if it is exact and the original binary mantissa is even
	// When disallowed, we can substract 1.
	uok := !du0 || fracu > 0
	if du0 && fracu == 0 {
		uok = mant&1 == 0
	}
	if !uok {
		du--
	}
	// Is 'dc' the correctly rounded base 10 mantissa?
	// The correct rounding might be dc+1
	cup := false // don't round up.
	if dc0 {
		// If we computed an exact product, the half integer
		// should round to next (even) integer if 'dc' is odd.
		cup = fracc > 1<<(extra-1) ||
			(fracc == 1<<(extra-1) && dc&1 == 1)
	} else {
		// otherwise, the result is a lower truncation of the ideal
		// result.
		cup = fracc>>(extra-1) == 1
	}
	// Is 'dl' an allowed representation?
	// Only if it is an exact value, and if the original binary mantissa
	// was even.
	lok := dl0 && fracl == 0 && (mant&1 == 0)
	if !lok {
		dl++
	}
	// We need to remember whether the trimmed digits of 'dc' are zero.
	c0 := dc0 && fracc == 0
	// render digits
	ryuDigits(d, dl, dc, du, c0, cup)
	d.dp -= q
	return
}

// log10pow2 returns q = math.Floor(e * log10(2)).
// It expects an input exponent such that |e| < 1600.
func log10pow2(e int) int {
	// log10(2) = 0.3010299956639812 = 78913.207... / 2**18
	return (e * 78913) >> 18 // even for negative e
}

// computeBounds returns a floating-point vector (l, c, u)×2^e2
// where the mantissas are 55-bit (or 25-bit) integers, describing the interval
// represented by the input float64 or float32.
func computeBounds(mant uint64, exp int, flt *floatInfo) (lower, central, upper uint64, e2 int) {
	if mant != 1<<flt.mantbits || exp-flt.bias == 1 {
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

func ryuDigits(d *decimalSlice, lower, central, upper uint64,
	c0, cup bool) {
	lhi, llo := divmod1e9(lower)
	chi, clo := divmod1e9(central)
	uhi, ulo := divmod1e9(upper)
	if uhi == 0 {
		// only low digits (for denormals)
		ryuDigits32(d, llo, clo, ulo,
			c0, cup, 8)
	} else if lhi < uhi {
		// truncate 9 digits at once.
		if llo != 0 {
			lhi++
		}
		c0 = c0 && clo == 0
		cup = (clo > 5e8) || (clo == 5e8 && cup)
		ryuDigits32(d, lhi, chi, uhi,
			c0, cup, 8)
		d.dp += 9
	} else {
		d.nd = 0
		// emit high part
		n := uint(9)
		for v := chi; v > 0; {
			v1, v2 := v/10, v%10
			v = v1
			n--
			d.d[n] = byte(v2 + '0')
		}
		d.d = d.d[n:]
		d.nd = int(9 - n)
		// emit low part
		ryuDigits32(d, llo, clo, ulo,
			c0, cup, d.nd+8)
	}
	// trim trailing zeros
	for d.nd > 0 && d.d[d.nd-1] == '0' {
		d.nd--
	}
	// trim initial zeros
	for d.nd > 0 && d.d[0] == '0' {
		d.nd--
		d.dp--
		d.d = d.d[1:]
	}
}

// ryuDigits32 emits decimal digits for a number less than 1e9.
func ryuDigits32(d *decimalSlice, lower, central, upper uint32,
	c0, cup bool, endindex int) {
	if upper == 0 {
		d.dp = endindex + 1
		return
	}
	trimmed := 0
	// Remember last trimmed digit to check for round-up.
	// c0 will be used to remember zeroness of following digits.
	cNextDigit := 0
	for upper > 0 {
		// Repeatedly compute:
		// l = Ceil(lower / 10^k)
		// c = Round(central / 10^k)
		// u = Floor(upper / 10^k)
		// and stop when c goes out of the (l, u) interval.
		l := (lower + 9) / 10
		c, cdigit := central/10, central%10
		u := upper / 10
		if l > u {
			// don't trim the last digit as it is forbidden to go below l
			// other, trim and exit now.
			break
		}
		// Check that we didn't cross the lower boundary.
		// The case where l < u but c == l-1 is essentially impossible,
		// but may happen if:
		//    lower   = ..11
		//    central = ..19
		//    upper   = ..31
		// and means that 'central' is very close but less than
		// an integer ending with many zeros, and usually
		// the "round-up" logic hides the problem.
		if l == c+1 && c < u {
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
	// We know where the number ends, fill directly
	endindex -= trimmed
	v := central
	n := endindex
	for n > d.nd {
		v1, v2 := v/100, v%100
		d.d[n] = smallsString[2*v2+1]
		d.d[n-1] = smallsString[2*v2+0]
		n -= 2
		v = v1
	}
	if n == d.nd {
		d.d[n] = byte(v + '0')
	}
	d.nd = endindex + 1
	d.dp = d.nd + trimmed
	return
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

// divmod10 computes quotient and remainder of division by 10,
// avoiding runtime uint64 division on 32-bit platforms.
func divmod10(x uint64) (uint64, uint64) {
	if !host32bit {
		return x / 10, x % 10
	}
	hi, _ := bits.Mul64(x, 0xcccccccccccccccd) // invpow10wide[1][0]+1
	q := hi >> 3
	return q, x - q*10
}

// divmod10 computes quotient and remainder of division by 100,
// avoiding runtime uint64 division on 32-bit platforms.
func divmod100(x uint64) (uint64, uint64) {
	if !host32bit {
		return x / 100, x % 100
	}
	hi, _ := bits.Mul64(x, 0xa3d70a3d70a3d70b) // invpow10wide[2][0]+1
	q := hi >> 6
	return q, x - q*100
}

// divmod1e9 computes quotient and remainder of division by 1e9,
// avoiding runtime uint64 division on 32-bit platforms.
func divmod1e9(x uint64) (uint32, uint32) {
	if !host32bit {
		return uint32(x / 1e9), uint32(x % 1e9)
	}
	// Use the same sequence of operations as the amd64 compiler.
	hi, _ := bits.Mul64(x>>1, 0x89705f4136b4a598) // invpow10wide[9][0]+1
	q := hi >> 28
	return uint32(q), uint32(x - q*1e9)
}
