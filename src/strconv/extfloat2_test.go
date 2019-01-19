// Copyright 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package strconv_test

import (
	"fmt"
	"math"
	"math/big"
	"math/bits"
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
		// trailing zeros
		1230000,
		12300000000,
		1.23e22, // exact integer with few digits
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

func TestRyuFtoaRandom(t *testing.T) {
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

func TestRyuFtoaFixedRandom(t *testing.T) {
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
		prec := int((uint64(i)*0xc0dedead)%16 + 1)

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

func TestRyuAtofRandom(t *testing.T) {
	// A standard desktop machine can check about 10e6 numbers/second.
	N := int(1e7)
	if testing.Short() {
		N = 1e6
	}
	ok, ko := 0, 0
	t.Logf("testing %d random numbers with fast and slow ParseFloat", N)

	for i := 0; i < N; i++ {
		mant := uint64(i) * 0xdeadbeefdeadbeef
		mant &= (1<<55 - 1)
		exp := (i % 640) - 342 // range (-342, 298)

		fold := OldAtof(mant, exp)
		bnew, _, ryuOk := RyuFromDecimal(mant, exp, &Float64info)
		fnew := math.Float64frombits(bnew)

		if !ryuOk || fold != fnew {
			t.Logf("%de%d old=%v new=%v", mant, exp, fold, fnew)
			ko++
		} else {
			ok++
		}
	}
	t.Logf("%d ok, %d ko", ok, ko)
}

// TestRyuAtofCoverage tests coverage of fast algorithms for the shortest
// representation of float64s. Grisu3 covers 99.5% of values, and naïve Ryū
// covers 92.9% (when rejecting 17-digit mantissas over 55 bits).
func TestRyuAtofCoverage(t *testing.T) {
	N := int(1e7)
	if testing.Short() {
		N = 1e6
	}
	t.Logf("testing ParseFloat of the shortest representation of %d random float64s", N)

	oldOk, ryuOk := 0, 0
	for i := 0; i < N; i++ {
		bits := uint64(i) * 0xdeadbeefdeadbeef
		bits = (bits << 1) >> 1 // clear sign bit
		if bits>>52 == 2047 {
			// only finite numbers
			bits ^= (1 << 60)
		}
		x := math.Float64frombits(bits)
		s := FormatFloat(x, 'g', -1, 64)

		// Compute decimal mantissa, exponent.
		mant, exp := ReadFloat(s)
		f1, ok1 := FastAtof(mant, exp)
		b2, _, ok2 := RyuFromDecimal(mant, exp, &Float64info)
		f2 := math.Float64frombits(b2)
		if ok1 {
			oldOk++
		}
		if ok2 {
			ryuOk++
		}
		if ok1 && ok2 && f1 != f2 {
			t.Fatalf("inconsistent results %s => %v %v", s, f1, f2)
		}
	}
	t.Logf("%d successes with old fast paths (%.2f%% failures)",
		oldOk, 100*float64(N-oldOk)/float64(N))
	t.Logf("%d successes with Ryū (%.2f%% failures)",
		ryuOk, 100*float64(N-ryuOk)/float64(N))
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
		expect := ExtFloat128{Hi: hi, Lo: lo, Exp: exp}
		if int(q) >= len(RyuPowersOfTen) || RyuPowersOfTen[q] != expect {
			t.Errorf("wrong entry, wants %#v", expect)
		}
	}

	for q := 0; q < 344; q++ {
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
		expect := ExtFloat128{Hi: hi, Lo: lo, Exp: exp}
		if int(q) >= len(RyuInvPowersOfTen) || RyuInvPowersOfTen[q] != expect {
			t.Errorf("wrong entry, wants %#v", expect)
		}
	}
}

// TestRyuExp2toExp10 checks that Exp2toExponent10(i) == math.Floor(i * log10(2))
func TestRyuExp2toExp10(t *testing.T) {
	for i := 1; i < 1600; i++ {
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

// TestRyuMultiplyNoCarry checks the fundamental assumption of the Ryu algorithm,
// which is that "ryuMultiply" upper 64 bits are always the exact truncation
// of m * 10^q even if the tables contain rounded floating-point versions of 10^q.
func TestRyuMultiplyCarry(t *testing.T) {
	// The tables contain (HI, LO, EXP) such that:
	// 10^q = (HI << 64 + LO ± ε) * 2**EXP
	// where ε is a residue < 1 and the sign of ±ε is the same as q.
	//
	// We are using ryuMultiply with input 55-bit wide (or more):
	//     k × 10^q = (K << 128 + H << 64 + L ± kε) * 2 ** EXP
	// so to return correct 64 bits, the upper 9 bits of H must be proved correct.
	// In other words, the hidden residue kε (less than k) must never
	// overflow when added to (k*P) mod (2^119), which would happen if:
	//     k × P = m × 2^119 + r
	//     r < k or r > 2^119-k (depending on the sign of ±ε above)
	// where P = (HI<<64|LO) mod 2^119
	//
	// Equivalently, we are looking for rational numbers m/k such that:
	//     0 < abs(m/k - P/2^119) < 1/2^119
	// with sign dependent on the sign of q, which is what this test is doing.

	check := func(inBits, outBits uint, pow ExtFloat128, e10 int) {
		// tests that when multiplying k (size inBits)
		// we can extract reliably the outBits top bits.
		Z := big.NewInt(1)
		Z = Z.Lsh(Z, 128+inBits-outBits)
		P := new(big.Int).SetUint64(pow.Hi)
		P = P.Lsh(P, 64)
		P = P.Or(P, new(big.Int).SetUint64(pow.Lo))

		Pm1 := new(big.Int).Set(P)
		if e10 > 0 {
			Pm1 = Pm1.Add(Pm1, big.NewInt(1))
		} else {
			Pm1 = Pm1.Sub(Pm1, big.NewInt(1))
		}

		_, q, err := findRational(
			new(big.Rat).SetFrac(Pm1, Z),
			new(big.Rat).SetFrac(P, Z))
		if err != nil {
			// overflow ensures that any 64-bit multiplier will be safe.
			return
		}
		product := new(big.Int).SetUint64(q)
		product = product.Mul(product, P)
		// If q has more than inBits, multiplier of that size are safe.
		if bits.Len64(q) <= int(inBits) {
			t.Errorf("10^%d -> 0x%x%016xp%d worst %d (%d bits, product=%x)",
				e10, pow.Hi, pow.Lo, pow.Exp, q, bits.Len64(q), product)
		} else if bits.Len64(q) <= int(inBits+2) {
			// no carry issue, but that was close.
			t.Logf("edge case (in=%db, out=%db): 10^%d -> 0x%x%016xp%d worst %d (%d bits, product=%x)",
				inBits, outBits, e10, pow.Hi, pow.Lo, pow.Exp, q, bits.Len64(q), product)
		}
	}

	for exp, pow := range RyuInvPowersOfTen {
		if exp <= 27 {
			// Multiplier 5^exp is the worst case (exact division),
			// we already know that.
			continue
		}
		check(55, 64, pow, -exp) // for ftoa
		check(60, 54, pow, -exp) // for atof
	}
	for exp, pow := range RyuPowersOfTen {
		if exp < 50 {
			continue // exact power of ten, no problem.
		}
		check(55, 64, pow, exp) // for ftoa
		check(60, 54, pow, exp) // for atof
	}
}

// findRational looks for a small rational p/q such that a < p/q < b.
// It panics if no such result can be found using uint64s.
func findRational(a, b *big.Rat) (p, q uint64, err error) {
	// The search will be based on continued fractions and the Stern-Brocot tree.
	// Assuming the continued fraction expansions of a and b are:
	//     [x1, x2, ...,xn, y, ...]
	//     [x1, x2, ...,xn, z, ...]
	// where y < z, then
	//   p/q = [x1, x2, ...,xn, y+1]
	// is the fraction with smallest denominator between a and b.
	// See https://en.wikipedia.org/wiki/Stern%E2%80%93Brocot_tree
	//     https://en.wikipedia.org/wiki/Continued_fraction#Best_rational_within_an_interval

	if a.Cmp(b) == 0 {
		return 0, 0, fmt.Errorf("a == b")
	}
	pa, qa := a.Num(), a.Denom()
	pb, qb := b.Num(), b.Denom()
	var pp, qq uint64
	p, pp = 1, 0
	q, qq = 0, 1
	overflow := func(x, y uint64) bool {
		return y > 0 && x >= (1<<64-1)/y
	}
	for qa.BitLen() > 0 && qb.BitLen() > 0 {
		// construct the xi as Euclidean quotients
		q1, r1 := new(big.Int), new(big.Int)
		q1, r1 = q1.DivMod(pa, qa, r1)
		q2, r2 := new(big.Int), new(big.Int)
		q2, r2 = q2.DivMod(pb, qb, r2)
		if !q1.IsUint64() || !q2.IsUint64() {
			return 0, 0, fmt.Errorf("overflow %s, %s", q1, q2)
		}
		x1 := q1.Uint64()
		x2 := q2.Uint64()
		// ...and accumulate in p, q.
		// Actually we are not interested in p so overflow will be ignored.
		switch {
		case x1 < x2:
			if overflow(x1+1, q) {
				return 0, 0, fmt.Errorf("overflow %d*%d+%d", q, x1+1, qq)
			}
			return p*(x1+1) + pp, q*(x1+1) + qq, nil
		case x1 > x2:
			if overflow(x2+1, q) {
				return 0, 0, fmt.Errorf("overflow %d*%d+%d", q, x2+1, qq)
			}
			return p*(x2+1) + pp, q*(x2+1) + qq, nil
		default:
			if overflow(x1, q) {
				return 0, 0, fmt.Errorf("overflow %d*%d+%d", q, x1, qq)
			}
			p, pp = p*x1+pp, p
			q, qq = q*x1+qq, q
		}
		pa, qa = qa, r1
		pb, qb = qb, r2
	}
	// If we reached zero remainder, it means that either
	// a == b, or that either a or b is actually the shortest rational.
	return p, q, nil
}

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

var ryuAtofBenches = []struct {
	name string
	mant uint64
	exp  int
}{
	{"64Decimal", 33909, 0},
	{"64Float", 3397784, -4},
	{"64FloatExp", 509, 73},
	{"64Denormal", 6226662346353213, -324},
	// Almost halfway, with less than 1e-16 ulp difference
	// with only 16 decimal digits.
	{"64HalfwayHard1", 6808957268280643, 116},  // from ftoahard
	{"64HalfwayHard2", 4334126125515466, -225}, // from ftoahard
	// Only 3e-13*ulp larger than halfway between denormals,
	{"64HalfwayDenormal", 168514038588815, -323},
	// Few digits, but 9.11691642378e-312 = 0x1ada385d67b.7fffffff5d9...p-1074
	// so naive, rounded 64-bit arithmetic is not enough to round it correctly.
	{"64HalfwayDenormalShort", 911691642378, -323},
	// 1.62420278e-315 = 0x1398359e.7fffe022p-1074,
	// should parsable using 64-bit arithmetic.
	{"64HalfwayDenormalVeryShort", 162420278, -323},
	// https://www.exploringbinary.com/php-hangs-on-numeric-value-2-2250738585072011e-308/
	{"64Denormal", 22250738585072011, -324},
}

func BenchmarkRyuAtof(b *testing.B) {
	for _, c := range ryuAtofBenches {
		b.Run(fmt.Sprintf("%de%d", c.mant, c.exp), func(b *testing.B) {
			var v uint64
			for i := 0; i < b.N; i++ {
				f, ovf, ok := RyuFromDecimal(c.mant, c.exp, &Float64info)
				if ovf || !ok {
					b.Fatal("could not parse")
				}
				v = f
			}
			b.Logf("parsed to %v", math.Float64frombits(v))
		})
	}
}

func BenchmarkOldAtof(b *testing.B) {
	for _, c := range ryuAtofBenches {
		b.Run(fmt.Sprintf("%de%d", c.mant, c.exp), func(b *testing.B) {
			var f float64
			for i := 0; i < b.N; i++ {
				f = OldAtof(c.mant, c.exp)
			}
			b.Logf("parsed to %v", f)
		})
	}
}
