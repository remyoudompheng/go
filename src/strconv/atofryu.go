// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package strconv

import (
	"math/bits"
)

// decimal to binary conversion following the principles of the Ryū algorithm.
//
// See Ulf Adams, "Ryū: Fast Float-to-String Conversion" (doi:10.1145/3192366.3192369)
//
// The steps are:
// 1) Ensure the decimal input fits in 31/64 bits (for float32, float64)
// 2) Multiply decimal by the requested power 10^e, truncating it down.
//    This can be done using a fixed precision representation of 10^e
//    tabulated in pow10wide and invpow10wide tables.
// 3) Remove extra bits and apply rounding rules to obtain mantissa.

func ryuAtof32(mant uint32, exp int) (fbits uint64, ovf bool) {
	// 32-bit decimal mantissas cannot be accepted because
	// of the (only) edge case 3465080571e-48 where we would
	// end up off by one.
	bitlen := bits.Len32(mant)
	if bitlen >= 32 {
		panic("ryuAtof32 called with decimal mantissa larger than 31 bits")
	}
	// Normalize input to be 31-bit wide.
	e2 := 0
	if bitlen < 31 {
		mant <<= uint(31 - bitlen)
		e2 = bitlen - 31
	}

	// multiply by a power of 10. It is required to know
	// whether the computation is exact.
	//
	// Using a 64-bit representation (P+ε)*2^E of 10^q, the following
	// property can be proved:
	//     (mant × P) >> 69
	// which is usually a 25-bit integer, is the correct lower
	// truncation of (mant × (P+ε)) >> 69.
	const atof32Shift = 69

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
		powexp = pow10exp(exp) - 63
	case exp == 0:
		// no multiply
	case exp < 0:
		powmant = invpow10wide[-exp][0] + 1 // round up
		powexp = pow10exp(exp) - 63
	}
	// Is it an exact computation?
	// Only small positive powers of 10 are exact (5^28 has 66 bits).
	exact := 27 >= exp && exp >= 0

	var di uint64
	if exp == 0 {
		di = uint64(mant)
	} else {
		// Multiply and shift right by 69 bits.
		hi, lo := bits.Mul64(uint64(mant), powmant)
		di = hi >> (atof32Shift - 64)
		d0 := hi&(1<<(atof32Shift-64)-1) == 0 && lo == 0
		if !d0 {
			exact = false
		}
		e2 += powexp + atof32Shift
	}
	// As a special case, computation might still be exact, if exponent
	// was negative and if it amounts to computing an exact division.
	// In that case, we ignore all lower bits.
	// Note that division by 10^14 cannot be exact as 5^14 has 33 bits.
	if exp < 0 && exp >= -13 && divisibleByPower5(uint64(mant), -exp) {
		exact = true
	}
	return ryuFloatBits(di, e2, exact, &float32info)
}

func ryuAtof64(mant uint64, exp int, flt *floatInfo) (fbits uint64, ovf bool) {
	bitlen := bits.Len64(mant)
	e2 := 0
	if bitlen < 64 {
		mant <<= uint(64 - bitlen)
		e2 = bitlen - 64
	}

	// multiply by a power of 10. It is required to know
	// whether the computation is exact.
	//
	// Using a 128-bit representation (P+ε)*2^E of 10^q, the following
	// property can be proved:
	//     (mant × P) >> 137
	// which is usually a 55-bit integer, is the correct lower
	// truncation of (mant × (P+ε)) >> 137.
	const atof64Shift = 137

	var pow [2]uint64 // a representation of 10^q
	var powexp int
	switch {
	case exp >= 309:
		// float64 cannot be larger than 1e309
		return 0x7ff << 52, true
	case exp < -342:
		// Even (1<<64-1) * 1e-343 converts to zero.
		return 0, false
	case exp > 0:
		pow = pow10wide[exp]
		powexp = pow10exp(exp) - 127
	case exp == 0:
		// no multiply
	case exp < 0:
		pow = invpow10wide[-exp]
		powexp = pow10exp(exp) - 127
	}
	// Is it an exact computation?
	// Only small positive powers of 10 are exact (5^55 has 128 bits).
	exact := exp <= 55 && exp >= 0

	// Compute the truncated down value of (x*10^q).
	var di uint64
	if exp == 0 {
		di = mant
	} else {
		// Compute (mant * pow) >> atof64Shift
		l1, l0 := bits.Mul64(mant, pow[1])
		h1, h0 := bits.Mul64(mant, pow[0])
		mid, carry := bits.Add64(l1, h0, 0)

		di = h1 + carry
		d0 := (di&(1<<(atof64Shift-128)-1) == 0) &&
			mid == 0 && l0 == 0
		if !d0 {
			exact = false
		}
		di >>= (atof64Shift - 128)
		e2 += powexp + atof64Shift
	}
	// As a special case, computation might still be exact, if exponent
	// was negative and if it amounts to computing an exact division.
	// In that case, we ignore all lower bits.
	// Note that division by 10^28 cannot be exact as 5^28 has 66 bits.
	if exp < 0 && exp >= -27 && divisibleByPower5(mant, -exp) {
		exact = true
	}
	return ryuFloatBits(di, e2, exact, flt)
}

// ryuFloatBits returns the correctly rounded floating point
// representation of a number x = (di + ε) * 2^e2,
// where ε is an unknown number < 1, which is known to be zero
// if 'exact' is true.
func ryuFloatBits(di uint64, e2 int, exact bool, flt *floatInfo) (fbits uint64, ovf bool) {
	mantbits := flt.mantbits
	bias := flt.bias
	expbits := flt.expbits
	// If the input mantissa is too large, truncate it.
	blen := bits.Len64(di)
	e2 += blen - 1
	extra := uint(blen) - mantbits - 1 // number of lower bits to remove
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
			(dfrac == 1<<(extra-1) && !exact) ||
			(dfrac == 1<<(extra-1) && exact && di&1 == 1)
	} else {
		// otherwise, d+1/2 always rounds up because
		// we truncated below.
		roundUp = dfrac>>(extra-1) == 1
	}
	if dfrac != 0 {
		exact = false
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
	fbits = di & (1<<mantbits - 1)
	fbits |= (uint64(e2-bias) & (1<<expbits - 1)) << mantbits
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
