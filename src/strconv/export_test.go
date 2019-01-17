// Copyright 2017 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package strconv

import (
	"math"
)

var (
	BitSizeError = bitSizeError
	BaseError    = baseError
)

type ExtFloat = extFloat
type ExtFloat128 = extfloat128
type ShortDecimal = decimalSlice
type Decimal = decimal

var Float64info = float64info

var BigFtoa = bigFtoa
var RoundShortest = roundShortest

func NewShortDecimal() decimalSlice {
	var buf [32]byte
	var d decimalSlice
	d.d = buf[:]
	return d
}

func ToShort(d *decimal) decimalSlice {
	sd := NewShortDecimal()
	copy(sd.d[:], d.d[:])
	sd.nd, sd.dp = d.nd, d.dp
	return sd
}

func (d1 *ShortDecimal) Equals(d2 *ShortDecimal) bool {
	if d1.nd != d2.nd && d1.dp != d2.dp {
		return false
	}
	for i := 0; i < d1.nd; i++ {
		if d1.d[i] != d2.d[i] {
			return false
		}
	}
	return true
}

func ShowDecimal(d *decimalSlice) string {
	if d.nd == 0 {
		return "0"
	}
	exp := d.dp - 1
	return string(d.d[0]) + "." + string(d.d[1:d.nd]) + "e" + Itoa(exp)
}

func OldAtof(mant uint64, exp int) float64 {
	if f, ok := atof64exact(mant, exp, false); ok {
		return f
	}

	// Try another fast path.
	ext := new(extFloat)
	if ok := ext.AssignDecimal(mant, exp, false, false, &float64info); ok {
		b, _ := ext.floatBits(&float64info)
		return math.Float64frombits(b)
	}

	var d decimal
	d.Assign(mant)
	d.dp += exp
	b, _ := d.floatBits(&float64info)
	return math.Float64frombits(b)
}

// FastAtof optimistically performs Atof, and only tries the float64 and Grisu fast paths.
func FastAtof(mant uint64, exp int) (float64, bool) {
	if f, ok := atof64exact(mant, exp, false); ok {
		return f, true
	}

	// Try another fast path.
	ext := new(extFloat)
	if ok := ext.AssignDecimal(mant, exp, false, false, &float64info); ok {
		b, _ := ext.floatBits(&float64info)
		return math.Float64frombits(b), true
	}

	return 0, false
}

func ReadFloat(s string) (mant uint64, exp int) {
	mantissa, exp, _, trunc, ok := readFloat(s)
	if !ok || trunc {
		panic("readFloat failure")
	}
	return mantissa, exp
}
