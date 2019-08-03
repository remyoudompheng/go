// Copyright 2017 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package cpu

const CacheLinePadSize = 32

// arm doesn't have a 'cpuid' equivalent, so we rely on HWCAP/HWCAP2.
// These are linknamed in runtime/os_(linux|freebsd)_arm.go and are
// initialized by archauxv().
// These should not be changed after they are initialized.
var HWCap uint
var HWCap2 uint

// HWCAP/HWCAP2 bits. These are exposed by Linux and FreeBSD.
const (
	hwcap_VFPv4 = 1 << 16
	hwcap_IDIVA = 1 << 17

	hwcap2_AES   = 1 << 0
	hwcap2_PMULL = 1 << 1
	hwcap2_SHA1  = 1 << 2
	hwcap2_SHA2  = 1 << 3
	hwcap2_CRC32 = 1 << 4
)

func doinit() {
	options = []option{
		{Name: "vfpv4", Feature: &ARM.HasVFPv4},
		{Name: "idiva", Feature: &ARM.HasIDIVA},

		{Name: "aes", Feature: &ARM.HasAES},
		{Name: "pmull", Feature: &ARM.HasPMULL},
		{Name: "sha1", Feature: &ARM.HasSHA1},
		{Name: "sha2", Feature: &ARM.HasSHA2},
		{Name: "crc32", Feature: &ARM.HasCRC32},
	}

	// HWCAP feature bits
	ARM.HasVFPv4 = isSet(HWCap, hwcap_VFPv4)
	ARM.HasIDIVA = isSet(HWCap, hwcap_IDIVA)

	// HWCAP2 feature bits
	ARM.HasAES = isSet(HWCap2, hwcap2_AES)
	ARM.HasPMULL = isSet(HWCap2, hwcap2_PMULL)
	ARM.HasSHA1 = isSet(HWCap2, hwcap2_SHA1)
	ARM.HasSHA2 = isSet(HWCap2, hwcap2_SHA2)
	ARM.HasCRC32 = isSet(HWCap2, hwcap2_CRC32)
}

func isSet(hwc uint, value uint) bool {
	return hwc&value != 0
}
