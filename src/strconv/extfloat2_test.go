// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package strconv_test

import (
	"fmt"
	"math"
	"math/big"
	. "strconv"
	"testing"
)

func TestRyuShortest(t *testing.T) {
	tests := []float64{
		// Basic cases
		0,
		12345,
		123451234512345,
		12345671234567890, // bounds are integers
		123456123456123456,
		// power of two: the bounds are x-1<<31 and x+1<<32
		// taking incorrect bounds leads to shortening.
		1 << 85,
		// negative power of 2,
		math.Ldexp(1, -1006),
		math.Ldexp(1, -1003),

		622666234635.3213e-320,
		3.0702409010742164e-151,
		// Mantissa is a (multiple of) power of 5
		math.Ldexp(476837158203125, 25),
		// 8450086975892335p+333 = 1.4785967089999999e+116
		// This is the case where we need to round up after
		// trimming more than 9 digits.
		1.478596709e+116,
		// 5320979720846975p+83 = 51461358161421999999999999999...
		5.1461358161422e+40,

		// Almost decimal midpoint: the decimal mantissa
		// is close to the midpoint between 2 integers.
		// The correct rounding direction requires precise computation.

		// 5549842152709732p-314 = 16628841458313728.501e-95 (rounds up)
		1.6628841458313729e-79,
		// 5693508526008355p-31 = 2651246.50951863965019 (rounds up)
		2.6512465095186397e6,

		// Decimal midpoints: the number is an exact decimal,
		// and the shortest representation removes a final '5'
		// Exactly representable as decimal 8436130790825957p-2
		2109032697706489.2,
		// 7868288860232000p-11 = 3841937920035.15625 exactly
		3.8419379200351562e12,
		// 4933269278992698p-3, must round (down) to even: 6.166586598740872e14
		616658659874087.25,
		// 4747277867737230p-3, must round (up) to even: 5.934097334671538e14
		593409733467153.75,

		// Binary midpoint is shorter and the round-to-even convention applies.

		// Actually 61861299179594376, upper midpoint 61861299179594380
		// is excluded by round-to-even convention
		6.1861299179594376e16,
		// 8750786544792037p+5, upper midpoint 280025169433345200 is
		// shorter but forbidden by the base 2 round-to-even convention.
		280025169433345184,
		// 5960464477539062p+71, upper midpoint is allowed: 140737488355328 × 10^23
		1.40737488355328e37,
		// 5960464477539063p+71, lower midpoint is forbidden
		1.4073748835532801e+37,
		// 8875518846412187p+13, upper midpoint is 72708250389808640000
		7.2708250389808636e19,
		// https://github.com/golang/go/issues/29491
		// lower midpoint 123382031554022200 is allowed
		1.233820315540222e+17,
	}
	for _, x := range tests {
		mant, exp := mantExp(x)

		dold := NewShortDecimal()
		dnew := NewShortDecimal()
		// slow algo
		oldShortest(&dold, mant, exp)
		// new algo
		RyuShortest(&dnew, mant, exp)
		if !dold.Equals(&dnew) {
			t.Error("ERROR:")
		}
		t.Logf("%b %v %v", x, ShowDecimal(&dold), ShowDecimal(&dnew))
	}
}

func oldShortest(d *ShortDecimal, mant uint64, exp int) {
	const neg = false
	const mantbits = 52

	f := new(ExtFloat)
	lower, upper := f.AssignComputeBounds(mant, exp, neg, &Float64info)
	ok := f.ShortestDecimal(d, &lower, &upper)
	if !ok {
		dec := new(Decimal)
		dec.Assign(mant)
		dec.Shift(exp - int(mantbits))
		RoundShortest(dec, mant, exp, &Float64info)
		*d = ToShort(dec)
	}
}

func TestRyuFixed(t *testing.T) {
	tests := []fixedTest{
		{1, 5}, // 1.0000

		{0.9, 1},    // 0.9
		{0.09, 1},   // 0.1
		{0.0999, 1}, // 0.1
		{0.05, 1},   // 0.1
		{0.5, 1},    // 0.5
		{1.5, 1},    // 2
		// Rounding final '5' in exact decimal
		{616658659874087.25, 16}, // 6.166586598740872e14
		{593409733467153.75, 16}, // 5.934097334671538e14
		{4428582144510765, 15},   // 4.42858214451076e14
		// Smallest denormal
		{5e-324, 15}, // 4.94065645841247e-324
	}
	for _, test := range tests {
		mant, exp := mantExp(test.x)

		dold := NewShortDecimal()
		dnew := NewShortDecimal()
		// slow algo
		oldFixed(&dold, mant, exp, test.prec)
		// new algo
		RyuFixed(&dnew, mant, exp, test.prec, &Float64info)
		if !dold.Equals(&dnew) {
			t.Error("ERROR:")
		}
		t.Logf("%b %.17g %v %v", test.x, test.x,
			ShowDecimal(&dold), ShowDecimal(&dnew))
	}
}

func oldFixed(d *ShortDecimal, mant uint64, exp int, prec int) {
	const mantbits = 52

	if prec <= 15 {
		// try fast algorithm when the number of digits is reasonable.
		f := new(ExtFloat)
		_, _ = f.AssignComputeBounds(mant, exp, false, &Float64info)
		ok := f.FixedDecimal(d, prec)
		if ok {
			return
		}
	}
	dec := new(Decimal)
	dec.Assign(mant)
	dec.Shift(exp - int(mantbits))
	dec.Round(prec)
	*d = ToShort(dec)
}

func TestRyuFtoa(t *testing.T) {
	// A standard desktop machine can check a few million numbers
	// per second.
	N := int(1e7)
	if testing.Short() {
		N = 1e6
	}
	ok, ko := 0, 0
	t.Logf("testing %d random numbers with fast and slow FormatFloat", N)

	dold := NewShortDecimal()
	dnew := NewShortDecimal()
	for i := 0; i < N; i++ {
		bits := uint64(i) * 0xdeadbeefdeadbeef
		bits = (bits << 1) >> 1 // clear sign bit
		if bits>>52 == 2047 {
			// only finite numbers
			bits ^= (1 << 60)
		}
		x := math.Float64frombits(bits)

		mant, exp := mantExp(x)

		// slow algo
		oldShortest(&dold, mant, exp)
		// new algo
		RyuShortest(&dnew, mant, exp)

		// compare
		if !dold.Equals(&dnew) {
			t.Logf("%b old=%s new=%s", x, ShowDecimal(&dold), ShowDecimal(&dnew))
			ko++
		} else {
			ok++
		}
	}
	t.Logf("%d ok, %d ko", ok, ko)
}

func TestRyuFtoaHard(t *testing.T) {
	const neg = false
	const mantbits = 52

	// test difficult cases. All these cases are rejected by Grisu3.
	hardFloats := GenerateHardFloat64s()
	t.Logf("testing %d float64 corner cases with Ryū and slow FormatFloat",
		len(hardFloats))

	for _, f := range hardFloats {
		mant, exp := mantExp(f)

		dold := NewShortDecimal()
		dnew := NewShortDecimal()
		// slow algo
		d := new(Decimal)
		d.Assign(mant)
		d.Shift(exp - int(mantbits))
		RoundShortest(d, mant, exp, &Float64info)
		dold = ToShort(d)
		// new algo
		RyuShortest(&dnew, mant, exp)

		// compare
		if !dold.Equals(&dnew) {
			t.Errorf("%b old=%s new=%s", f,
				ShowDecimal(&dold), ShowDecimal(&dnew))
		}
	}
}

func TestRyuFtoaFixed(t *testing.T) {
	// A standard desktop machine can check a few million numbers
	// per second.
	N := int(1e7)
	if testing.Short() {
		N = 1e6
	}
	ok, ko := 0, 0
	t.Logf("testing %d random numbers with fast and slow FormatFloat", N)

	dold := NewShortDecimal()
	dnew := NewShortDecimal()
	for i := 0; i < N; i++ {
		bits := uint64(i) * 0xdeadbeefdeadbeef
		bits = (bits << 1) >> 1 // clear sign bit
		if bits>>52 == 2047 {
			// only finite numbers
			bits ^= (1 << 60)
		}
		x := math.Float64frombits(bits)
		prec := (int(i)*0xc0dedead)%16 + 1

		mant, exp := mantExp(x)

		// slow algo
		oldFixed(&dold, mant, exp, prec)
		// new algo
		RyuFixed(&dnew, mant, exp, prec, &Float64info)

		// compare
		if !dold.Equals(&dnew) {
			t.Logf("%b old=%s new=%s", x, ShowDecimal(&dold), ShowDecimal(&dnew))
			ko++
		} else {
			ok++
		}
	}
	t.Logf("%d ok, %d ko", ok, ko)
}

func mantExp(x float64) (mant uint64, exp int) {
	const (
		mantbits = 52
		expbits  = 11
		bias     = -1023
	)
	bits := math.Float64bits(x)
	exp = int(bits>>mantbits) & (1<<expbits - 1)
	mant = bits & (uint64(1)<<mantbits - 1)

	switch exp {
	case 1<<expbits - 1:
		return 0, 0
	case 0:
		// denormalized
		exp++
	default:
		// add implicit top bit
		mant |= uint64(1) << mantbits
	}
	exp += bias
	return mant, exp
}

func TestRyuPowersOfTen(t *testing.T) {
	for q := int64(0); q <= 339; q++ {
		// positive exponents (324 to 339) are needed for fixed
		// precision handling of denormals
		pow := big.NewInt(10)
		pow = pow.Exp(pow, big.NewInt(q), nil)
		sz := pow.BitLen()
		if sz < 128 {
			pow = pow.Lsh(pow, uint(128-sz))
		} else {
			pow = pow.Rsh(pow, uint(sz-128))
		}
		digits := pow.Bytes()
		hi := new(big.Int).SetBytes(digits[:8]).Uint64()
		lo := new(big.Int).SetBytes(digits[8:]).Uint64()
		exp := sz - 128
		//t.Logf(`{Hi: 0x%016x, Lo: 0x%016x, Exp: %d},`, hi, lo, exp)
		expect := ExtFloat128{Hi: hi, Lo: lo, Exp: exp}
		if int(q) >= len(RyuPowersOfTen) || RyuPowersOfTen[q] != expect {
			t.Errorf("wrong entry, wants %#v", expect)
		}
	}

	for q := 0; q < 320; q++ {
		// negative exponents
		// Let's compute 2^128 * 8^q / 5^q, which is never an integer
		pow := big.NewInt(5)
		pow = pow.Exp(pow, big.NewInt(int64(q)), nil)
		p := big.NewInt(1)
		p = p.Lsh(p, 3*uint(q)+128)

		bits := p.Div(p, pow)
		sz := bits.BitLen()
		if sz < 128 {
			bits = bits.Lsh(bits, uint(128-sz))
		} else {
			bits = bits.Rsh(bits, uint(sz-128))
		}
		digits := bits.Bytes()
		hi := new(big.Int).SetBytes(digits[:8]).Uint64()
		lo := new(big.Int).SetBytes(digits[8:]).Uint64()
		if q > 0 {
			lo++ // round up
		}
		exp := sz - 256 - 4*q
		//t.Logf(`{Hi: 0x%016x, Lo: 0x%016x, Exp: %d},`, hi, lo, exp)
		expect := ExtFloat128{Hi: hi, Lo: lo, Exp: exp}
		if RyuInvPowersOfTen[q] != expect {
			t.Errorf("wrong entry")
		}
	}
}

func TestRyuExp2toExp10(t *testing.T) {
	for i := 1; i < 1600; i++ {
		// Is it really math.Floor(i * log10(2))
		exact := math.Ln2 / math.Ln10 * float64(i)
		approx := Exp2toExponent10(uint(i))

		if exact < float64(approx)+0.0001 {
			t.Fatalf("%d*log10(2): approx=%d, exact=%v",
				i, approx, exact)
		}
		if exact > float64(approx)+0.9999 {
			t.Fatalf("%d*log10(2): approx=%d, exact=%v",
				i, approx, exact)
		}
	}
}

/*
// TestRyuNoCarry checks the fundamental assumption of the Ryu algorithm.
// Let x = mant × 2**exp be a floating-point number used as upper/lower bound (mant <= 2**54)
// Let 10^q = (M + ε) × 2**E be the corresponding entry in the ryuPowersOfTen table.
// then the upper 64 bits of x × 10^q can always be computed by omitting ε,
// i.e. the contribution of ε cannot generate enough carries in the multiplication.
func TestRyuNoCarry(t *testing.T) {
	t.Skip("skip")
	// Let B = 54 be the number of bits in mant
	// Then mant × 10^q = mant × M + mant × ε
	//                  = H + L + mant × ε
	// where H are the upper 64 bits
	//       L are the lower 64+B bits
	//       mant×ε is less than 2**B
	//
	// Note than L = mant × M mod 2**(64+B)
	// so it is enough to prove that L < 2**(64+B) - 2**B always
	// to prevent the carry from propagating to H.
	// This typically happens if M is pseudo-random enough.
	//
	// In Grisu3, where the precision is 64 bits:
	//     H are the upper 64 bits
	//     L are the lower B bits
	//     mant×ε can be up to B bits
	// so H is usually uncertain by 1.
	testNoCarry(t, 0xbd2b1a2d54a54a58, 0xd5c4b8f4a5c4d7ef)
}

// testNoCarry checks that if m = hi<<64|lo
// then k*m mod W = 2**(64+B) is always less than 2**(64+B) - 2**B.
func testNoCarry(t *testing.T, hi, lo uint64) {
	const B = 54

	mask := big.NewInt(1)
	mask = mask.Lsh(mask, 64+B)
	mask = mask.Sub(mask, big.NewInt(1))
	// The bound that should not be passed
	bound := big.NewInt(1)
	bound = bound.Lsh(bound, 64+B)
	bound = bound.Sub(bound, big.NewInt(1<<B))
	// Prepare m as big integer, truncated to 128-B bits.
	m := new(big.Int).SetUint64(hi)
	m = m.Lsh(m, 64)
	m = m.Add(m, new(big.Int).SetUint64(lo))
	m = m.And(m, mask)

	// Look for a small stride such that (stride*m) mod W is very small.
	// We want a stride of order 2**(B/2) to reduce the complexity to sqrt(W).
	x := big.NewInt(0)
	var stride uint64
	step := ^big.Word(0)
	for s := uint64(1); s < 1<<(B/2); s++ {
		x = x.Add(x, m)
		x = x.And(x, mask)
		if x.BitLen() > 64 && x.Bits()[1] < step {
			step = x.Bits()[1]
			stride = s
			println(stride, step)
		}
	}
	if stride == 0 {
		t.Fatal("no small stride")
	}
	t.Logf("m=%s, found stride %d, stride*m~=%d<<64",
		m, stride, step)
	println("stride", stride, "step", step)

	multiply := func(n *big.Int, k uint64) {
		n.SetUint64(k)
		n.Mul(n, m)
		n.And(n, mask)
		if n.Cmp(bound) >= 0 {
			t.Fatalf("found %d * %s = %x, bound=%x", k, m, n, bound)
		}
	}

	// Now look for bad values.
	// All possible multipliers can be written as:
	// k = b + a*stride
	// It is not necessary to test all a's if we can skip.
	skip := (1 << (B - 2)) / step
	t.Logf("using skip=%d", skip)
	println("skip", skip, "max strides", (1<<B)/stride)
	for b := uint64(1); b <= stride; b++ {
		if b%1000 == 0 {
			println("b =", b)
		}
		k := b
		x := new(big.Int)
		for k < 1<<B {
			multiply(x, k)

			topWord := big.Word(0)
			if x.BitLen() > 64 {
				topWord = x.Bits()[1]
			}
			if topWord < 3<<(B-2) {
				//println("testing", k, "product", x.String(), "limit", bound.String(), "skip", skipWord*skip)
				k += uint64(skip) * stride
				println("now k =", k)
			} else {
				// be careful
				k += stride
			}
		}
	}
}
*/

func BenchmarkRyuShortest(b *testing.B) {
	d := NewShortDecimal()
	for _, c := range ftoaBenches {
		b.Run(c.name, func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				mant, exp := mantExp(c.float)
				RyuShortest(&d, mant, exp)
			}
		})
	}
}

func BenchmarkRyuShortestHard(b *testing.B) {
	d := NewShortDecimal()
	for _, c := range hardFloatSamples {
		b.Run(fmt.Sprintf("%b", c), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				mant, exp := mantExp(c)
				RyuShortest(&d, mant, exp)
			}
		})
	}
}

type fixedTest struct {
	x    float64
	prec int
}

var fixedBenches = []fixedTest{
	{0.9, 1}, // 0.9
	{0.5, 1}, // 0.5
	{1.5, 1}, // 2
	{123456, 3},
	{123.456, 3},
	{1.23456e+78, 3},
	{1.23456e-78, 3},
	{1.23456e-78, 15},
	{4428582144510765, 15}, // 4.42858214451076e14
	{1.2345678912345654e+249, 15},
	{5e-324, 15}, // 4.94065645841247e-324
	{1.23456e-78, 16},
	{4428582144510765, 16}, // 4.42858214451076e14
	{1.2345678912345654e+249, 16},
	{5e-324, 16}, // 4.94065645841247e-324
}

func BenchmarkRyuFixed(b *testing.B) {
	d := NewShortDecimal()
	for _, c := range fixedBenches {
		b.Run(fmt.Sprintf("%v/%d", c.x, c.prec), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				mant, exp := mantExp(c.x)
				RyuFixed(&d, mant, exp, c.prec, &Float64info)
			}
		})
	}
}

func BenchmarkOldFixed(b *testing.B) {
	d := NewShortDecimal()
	for _, c := range fixedBenches {
		b.Run(fmt.Sprintf("%v/%d", c.x, c.prec), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				mant, exp := mantExp(c.x)
				oldFixed(&d, mant, exp, c.prec)
			}
		})
	}
}
