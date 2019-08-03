package aes

// aesCipherGCM implements crypto/cipher.gcmAble so that crypto/cipher.NewGCM
// will use the optimised implementation in this file when possible. Instances
// of this type only exist when hasGCMAsm returns true.
type aesCipherGCM struct {
	aesCipherAsm
}
