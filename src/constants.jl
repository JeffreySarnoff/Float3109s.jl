const IntFormat = Int16

const ByteDatum = UInt8
const WordDatum = UInt16
const UIntDatum = Union{ByteDatum, WordDatum}

# short string forms for types
const UnsignedFormatChar = '𝙐'        # \bisansU
const SignedFormatChar = '𝙎'          # \bisansS
const UnsignedStr = string(UnsignedFormatChar)
const SignedStr = string(SignedFormatChar)

const FiniteFormatChar = '𝙁'          # \bisansF
const ExtendedFormatChar = '𝙀'        # \bisansE
const FiniteStr = string(FiniteFormatChar)
const ExtendedStr = string(ExtendedFormatChar)

# short string forms for instances
const UnsignedInstanceChar = '𝘜'  # \isansU
const SignedInstanceChar = '𝘚'    # \isansS
const UnsignedFormatStr = string(UnsignedInstanceChar)
const SignedFormatStr = string(SignedInstanceChar)

const FiniteInstanceChar = '𝘍'    # \isansF
const ExtendedInstanceChar = '𝘌'  # \isansE
const FiniteFormatStr = string(FiniteInstanceChar)
const ExtendedFormatStr = string(ExtendedInstanceChar)


