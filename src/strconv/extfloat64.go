package strconv

import (
	"math/bits"
)

func ryuFromDecimal32(mant uint64, exp int) (fbits uint32, ovf bool) {
	const mantbits = 23
	const expbits = 8
	const bias = -127
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
	// A representation of 10^q: powmant × 2^powexp
	var powmant uint64
	var powexp int
	switch {
	case exp >= 39:
		// float32 cannot be larger than 1e39
		return 0xff << 23, true
	case exp <= -65:
		// float32 cannot be smaller than (1<<64) * 1e-65
		return 0, false
	case exp > 0:
		powmant = pow10wide[exp][0]
		powexp = exp10toExponent2(exp) - 63
	case exp == 0:
		// no multiply
	case exp < 0:
		powmant = invpow10wide[-exp][0] + 1 // round up
		powexp = exp10toExponent2(exp) - 63
	}
	// Is it an exact computation?
	exact := false
	switch {
	case 27 >= exp && exp >= 0:
		exact = true
	case 0 > exp && exp >= -27:
		// division by a power of ten might be exact
		// if mantissas are multiples of 5
		if divisibleByPower5(mant, -exp) {
			exact = true
		}
	default:
		// positive powers of ten (>= 28) are not exact

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
		hi, lo := bits.Mul64(mant, powmant)
		di, d0 = hi, lo == 0
		e2 += powexp + 64
	}
	// If computation was an exact division, lower bits must be ignored.
	if exp < 0 && exact {
		d0 = true
	}
	// Is exponent too low? Shrink mantissa for denormals.
	blen := bits.Len64(di)
	e2 += blen - 1
	extra := uint(blen - mantbits - 1) // number of lower bits to remove
	if e2 < bias+1 {
		extra += uint(bias + 1 - e2)
		e2 = bias + 1
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
	if di == 2<<mantbits {
		di >>= 1
		e2++
	}

	// Infinities.
	if e2-bias >= 1<<expbits-1 {
		// ±Inf
		di = 0
		e2 = 1<<expbits - 1 + bias
		ovf = true
	} else if di&(1<<mantbits) == 0 {
		// Denormalized?
		e2 = bias
	}
	// Assemble bits.
	fbits = uint32(di) & (uint32(1)<<mantbits - 1)
	fbits |= uint32((e2-bias)&(1<<expbits-1)) << mantbits
	return fbits, ovf
}

// exp10toExponent2 returns math.Floor(q * log2(10))
// It is only guaranteed to work when -500 <= q <= 500.
func exp10toExponent2(q int) int {
	// log2(10) = 3.3219280948873626 ~= 108853 / 2**15
	return (q * 108853) >> 15 // even for negative e
}
