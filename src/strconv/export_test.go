// Copyright 2017 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package strconv

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
