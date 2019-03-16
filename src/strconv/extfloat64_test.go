package strconv_test

import (
	"fmt"
	. "strconv"
	"strings"
	"testing"
)

func TestRyuFromDecimal32(t *testing.T) {
	for exp := -2000; exp < 2000; exp++ {
		tests := []string{
			fmt.Sprintf("1e%d", exp),
			fmt.Sprintf("1000000000000000000e%d", exp),
			fmt.Sprintf("999999999999999999e%d", exp),
		}
		for _, s := range tests {
			f64, err := ParseFloat(s, 64)
			if err != nil && strings.HasSuffix(err.Error(), ErrRange.Error()) {
				err = nil
			}
			if err != nil {
				t.Fatal(err)
			}
			f32, err := ParseFloat(s, 32)
			if err != nil && strings.HasSuffix(err.Error(), ErrRange.Error()) {
				err = nil
			}
			if err != nil {
				t.Fatal(err)
			}
			if float32(f64) != float32(f32) {
				t.Errorf("%s => %v %v", s, f64, float32(f32))
			}
			t.Logf("%s => %v %v", s, f64, float32(f32))
		}
	}
}
