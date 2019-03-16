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
		powmant = ryuPowersOfTen32[exp]
		powexp = exp10toExponent2(exp) - 63
	case exp == 0:
		// no multiply
	case exp < 0:
		powmant = ryuInvPowersOfTen32[-exp]
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

// ryuPowersOfTen32[q] is the 64-bit mantissa of 10^q.
// ryuPowersOfTen32[q] = 10^q << 63 >> floor(q*log2(10))
//                     = ryuPowersOfTen[q].Hi
var ryuPowersOfTen32 = [...]uint64{
	0x8000000000000000,
	0xa000000000000000,
	0xc800000000000000,
	0xfa00000000000000,
	0x9c40000000000000,
	0xc350000000000000,
	0xf424000000000000,
	0x9896800000000000,
	0xbebc200000000000,
	0xee6b280000000000,
	0x9502f90000000000,
	0xba43b74000000000,
	0xe8d4a51000000000,
	0x9184e72a00000000,
	0xb5e620f480000000,
	0xe35fa931a0000000,
	0x8e1bc9bf04000000,
	0xb1a2bc2ec5000000,
	0xde0b6b3a76400000,
	0x8ac7230489e80000,
	0xad78ebc5ac620000,
	0xd8d726b7177a8000,
	0x878678326eac9000,
	0xa968163f0a57b400,
	0xd3c21bcecceda100,
	0x84595161401484a0,
	0xa56fa5b99019a5c8,
	0xcecb8f27f4200f3a,
	0x813f3978f8940984,
	0xa18f07d736b90be5,
	0xc9f2c9cd04674ede,
	0xfc6f7c4045812296,
	0x9dc5ada82b70b59d,
	0xc5371912364ce305,
	0xf684df56c3e01bc6,
	0x9a130b963a6c115c,
	0xc097ce7bc90715b3,
	0xf0bdc21abb48db20,
	0x96769950b50d88f4,
	0xbc143fa4e250eb31,
	0xeb194f8e1ae525fd,
	0x92efd1b8d0cf37be,
	0xb7abc627050305ad,
	0xe596b7b0c643c719,
	0x8f7e32ce7bea5c6f,
	0xb35dbf821ae4f38b,
	0xe0352f62a19e306e,
	0x8c213d9da502de45,
	0xaf298d050e4395d6,
	0xdaf3f04651d47b4c,
}

// ryuInvPowersOfTen32[q] is the 64-bit mantissa (rounded up) of 10^-q.
// ryuInvPowersOfTen32[q]
// = ryuInvPowersOfTen[q].Hi + 1    (when q > 0)
// = Ceil(2^(64+q*log2(10)) / 10^q) (when q > 0)
var ryuInvPowersOfTen32 = [...]uint64{
	0x8000000000000000,
	0xcccccccccccccccd,
	0xa3d70a3d70a3d70b,
	0x83126e978d4fdf3c,
	0xd1b71758e219652c,
	0xa7c5ac471b478424,
	0x8637bd05af6c69b6,
	0xd6bf94d5e57a42bd,
	0xabcc77118461cefd,
	0x89705f4136b4a598,
	0xdbe6fecebdedd5bf,
	0xafebff0bcb24aaff,
	0x8cbccc096f5088cc,
	0xe12e13424bb40e14,
	0xb424dc35095cd810,
	0x901d7cf73ab0acda,
	0xe69594bec44de15c,
	0xb877aa3236a4b44a,
	0x9392ee8e921d5d08,
	0xec1e4a7db69561a6,
	0xbce5086492111aeb,
	0x971da05074da7bef,
	0xf1c90080baf72cb2,
	0xc16d9a0095928a28,
	0x9abe14cd44753b53,
	0xf79687aed3eec552,
	0xc612062576589ddb,
	0x9e74d1b791e07e49,
	0xfd87b5f28300ca0e,
	0xcad2f7f5359a3b3f,
	0xa2425ff75e14fc32,
	0x81ceb32c4b43fcf5,
	0xcfb11ead453994bb,
	0xa6274bbdd0fadd62,
	0x84ec3c97da624ab5,
	0xd4ad2dbfc3d07788,
	0xaa242499697392d3,
	0x881cea14545c7576,
	0xd9c7dced53c72256,
	0xae397d8aa96c1b78,
	0x8b61313bbabce2c7,
	0xdf01e85f912e37a4,
	0xb267ed1940f1c61d,
	0x8eb98a7a9a5b04e4,
	0xe45c10c42a2b3b06,
	0xb6b00d69bb55c8d2,
	0x9226712162ab070e,
	0xe9d71b689dde71b0,
	0xbb127c53b17ec15a,
	0x95a8637627989aae,
}
