package strconv_test

import (
	"fmt"
	"math"
	"math/big"
	"math/bits"
	. "strconv"
	"testing"
)

/*
This test generates floating point numbers which are hard
for the "shortest decimal" problem.

Let f  = m × 2**e (where m and e are integers)
    f+ = (2m+1) × 2**(e-1)
    f- = (2m-1) × 2**(e-1)

Let q be the smallest exponent such that 10^q × 2^e > 1
i.e. q = Floor(log10(2) * -e) + 1
then [f- × 10^q, f+ × 10^q] contains at least one integer,
and the shortest decimal for f is n × 10^-q where n belongs
to that interval (and is divisible by the largest power of 10).

A floating point number f is "hard" if f± × 10^q is very
close to an integer. Sample mantissas for these corner cases
can be found by computing continued fractions.

These values are typically rejected by the Grisu3 algorithm.
*/
func GenerateHardFloat64s() []float64 {
	var hards []float64
	for e := -1022 - 52; e <= 1023-52; e++ {
		if -10 <= e && e <= 10 {
			continue // nothing interesting here
		}

		q := int(math.Floor(math.Ln2/math.Ln10*float64(-e))) + 1
		//step := math.Ldexp(math.Pow(10, float64(q)), (e - 1))
		//t.Logf("exponent %d, power 10^%d, step=%.6f", e, q, step)

		// We are looking for a fraction x/y very close
		// to wd = 10^q × 2^(e-1), where y is a 54-bit odd integer.
		// Also, x should be a multiple of 10 to be a candidate
		// for shortest decimal.
		var x, y uint64
		var prec float64
		if q >= 0 { // e < 0
			a := big.NewInt(10)
			a = a.Exp(a, big.NewInt(int64(q)), nil)
			b := big.NewInt(1)
			b = b.Lsh(b, uint(-(e - 1)))
			x, y, prec = findFrac(a, b, 54)
		} else {
			a := big.NewInt(10)
			a = a.Exp(a, big.NewInt(int64(-q)), nil)
			b := big.NewInt(1)
			b = b.Lsh(b, uint(e-1))
			x, y, prec = findFrac(b, a, 54)
		}

		if bits.Len64(y) == 54 {
			f := math.Ldexp(float64(y>>1), e)
			_ = fmt.Sprintf("f=ldexp(%d,%d)=%v, f+=(%d+%.3e)e%d\n",
				y>>1, e, f, x, prec, -q)
			hards = append(hards, f)
		}

		if e == -1074 {
			for bitlen := 30; bitlen < 54; bitlen++ {
				// also find hard denormals
				a := big.NewInt(10)
				a = a.Exp(a, big.NewInt(int64(q)), nil)
				b := big.NewInt(1)
				b = b.Lsh(b, uint(-(e - 1)))
				x, y, prec = findFrac(a, b, bitlen)

				f := math.Ldexp(float64(y>>1), e)
				_ = fmt.Sprintf("f=ldexp(%d,%d)=%v, f+=(%d+%.3e)e%d\n",
					y>>1, e, f, x, prec, -q)
				hards = append(hards, f)
			}
		}
	}
	return hards
}

// findFrac returns a fraction x/y very close to u/v,
// such that y*(u/v) = x+prec
func findFrac(u, v *big.Int, bitlen int) (x, y uint64, prec float64) {
	q := new(big.Rat).SetFrac(u, v)
	for seed := uint64(1); seed < 90; seed += 3 {
		x, y = contFrac(u, v, seed, 1<<uint(bitlen-1))
		if bits.Len64(y) == bitlen && y%2 == 1 && x%10 == 0 {
			break
		}
	}
	q = q.Mul(q, big.NewRat(int64(y), 1))
	q = q.Sub(q, big.NewRat(int64(x), 1))
	prec, _ = q.Float64()
	return x, y, prec
}

func contFrac(u, v *big.Int, seed uint64, max uint64) (x, y uint64) {
	var a, b uint64 = 1, 0
	var c, d uint64 = 0, seed
	for c < max {
		if v.BitLen() == 0 {
			break
		}
		q, r := new(big.Int), new(big.Int)
		q, r = q.DivMod(u, v, r)
		if !q.IsUint64() {
			panic("!q.IsUint64")
		}
		quo := q.Uint64()
		a, b = quo*a+b, a
		c, d = quo*c+d, c
		u, v = v, r
	}
	return a * seed, c
}

var hardFloatSamples = []float64{
	// Denormals
	math.Ldexp(328742302, -1074),
	math.Ldexp(1845284427387, -1074),
	math.Ldexp(341076211242912, -1074),
	// Difficulty < 1e-15
	math.Ldexp(6417092537094053, -748),
	math.Ldexp(7675932596762664, -653),
	math.Ldexp(6419534400875886, -426),
	math.Ldexp(4566633709189828, -328),
	math.Ldexp(8640368759831959, 385),
	math.Ldexp(6503767923869541, 602),
	math.Ldexp(5662764645683412, 635),
	math.Ldexp(7953761449385755, 828),
	math.Ldexp(7953761449385755, 831),
	math.Ldexp(6018986745823044, 858),
	math.Ldexp(6018986745823044, 861),
	math.Ldexp(6018986745823044, 862),
	math.Ldexp(4787903260141515, 897),
	math.Ldexp(5349337776366262, 949),
	math.Ldexp(6073849323345086, 962),
	// Difficulty < 1e-14
	math.Ldexp(5969291480317302, -146),
	math.Ldexp(5130627738529412, -134),
	math.Ldexp(5130627738529412, -133),
	math.Ldexp(6931776026129216, -131),
	math.Ldexp(6146622122784629, -99),
	math.Ldexp(4528599518205136, -81),
	math.Ldexp(5660749397756420, -78),
	math.Ldexp(8040837212722187, -75),
	math.Ldexp(4576042559928398, 81),
	math.Ldexp(4576042559928398, 82),
	math.Ldexp(5853077692931672, 84),
	math.Ldexp(4800294408018791, 89),
	math.Ldexp(5240375412144155, 104),
	math.Ldexp(6319502805243561, 114),
	math.Ldexp(7869598596808504, 127),
	math.Ldexp(5889671799622512, 138),
	math.Ldexp(5889671799622512, 139),
	math.Ldexp(5353445750064544, 148),
}

func TestFtoaHard(t *testing.T) {
	hards := GenerateHardFloat64s()
	t.Logf("testing %d float64 corner cases with fast and slow FormatFloat",
		len(hards))
	for _, x := range hards {
		shortFast := FormatFloat(x, 'g', -1, 64)
		SetOptimize(false)
		shortSlow := FormatFloat(x, 'g', -1, 64)
		SetOptimize(true)
		if shortSlow != shortFast {
			t.Errorf("%b printed as %s, want %s", x, shortFast, shortSlow)
		}

		for prec := 5; prec < 17; prec++ {
			shortFast = FormatFloat(x, 'e', prec, 64)
			SetOptimize(false)
			shortSlow = FormatFloat(x, 'e', prec, 64)
			SetOptimize(true)
			if shortSlow != shortFast {
				t.Errorf("%b printed as %s, want %s", x, shortFast, shortSlow)
			}
		}
	}
}

func BenchmarkAppendFloatHard(b *testing.B) {
	dst := make([]byte, 30)
	for _, c := range hardFloatSamples {
		b.Run(fmt.Sprintf("%b", c), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				AppendFloat(dst[:0], c, 'g', -1, 64)
			}
		})
	}
}
