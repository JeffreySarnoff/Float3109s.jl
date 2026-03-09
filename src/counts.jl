"""Number of NaN code points (always 1)."""
nNaNsOf(@nospecialize fmt::Format) = 1

"""Number of zero code points (always 1)."""
nZerosOf(@nospecialize fmt::Format) = 1

"""Total number of infinity code points: 0 (finite), 1 (unsigned extended), 2 (signed extended)."""
nInfsOf(@nospecialize fmt::Format{T,is_finite}) where T = 0
nInfsOf(@nospecialize fmt::Format{is_unsigned,is_extended}) = 1
nInfsOf(@nospecialize fmt::Format{is_signed,is_extended}) = 2

"""Number of positive infinity code points: 0 or 1."""
nPosInfsOf(@nospecialize fmt::Format{T,is_finite}) where T = 0
nPosInfsOf(@nospecialize fmt::Format{is_unsigned,is_extended}) = 1
nPosInfsOf(@nospecialize fmt::Format{is_signed,is_extended}) = 1

"""Number of negative infinity code points: 0 or 1."""
nNegInfsOf(@nospecialize fmt::Format{T,is_finite}) where T = 0
nNegInfsOf(@nospecialize fmt::Format{is_unsigned,is_extended}) = 0
nNegInfsOf(@nospecialize fmt::Format{is_signed,is_extended}) = 1

"""Total number of code points: `2^K`."""
nValuesOf(@nospecialize fmt::Format) = twopow(BitwidthOf(fmt))

"""Number of numerical (non-NaN) code points: `2^K - 1`."""
nNumericalValuesOf(@nospecialize fmt::Format) = nValuesOf(fmt) - 1

"""Number of non-zero numerical code points."""
nNonZeroNumericalValuesOf(@nospecialize fmt::Format) = nNumericalValuesOf(fmt) - nZerosOf(fmt)

"""Number of finite values (numerical minus infinities), including zero."""
nFiniteValuesOf(@nospecialize fmt::Format) = nNumericalValuesOf(fmt) - nInfsOf(fmt)

"""Number of non-zero finite values."""
nNonZeroFiniteValuesOf(@nospecialize fmt::Format) = nFiniteValuesOf(fmt) - nZerosOf(fmt)

"""Number of negative non-zero values (including -Inf if present). Always 0 for unsigned."""
nNegValuesOf(@nospecialize fmt::Format{is_unsigned,T}) where T = 0
nNegValuesOf(@nospecialize fmt::Format{is_signed,T}) where T = nNonZeroNumericalValuesOf(fmt) ÷ 2

"""Number of positive non-zero values (including +Inf if present)."""
nPosValuesOf(@nospecialize fmt::Format{is_unsigned,T}) where T = nNonZeroNumericalValuesOf(fmt)
nPosValuesOf(@nospecialize fmt::Format{is_signed,T}) where T = nNonZeroNumericalValuesOf(fmt) ÷ 2

"""Number of non-negative values (positive values plus zero)."""
nNonNegValuesOf(@nospecialize fmt::Format) = nPosValuesOf(fmt) + nZerosOf(fmt)

"""Number of negative finite values."""
nNegFiniteValuesOf(@nospecialize fmt::Format) = nNegValuesOf(fmt) - nNegInfsOf(fmt)

"""Number of positive finite values."""
nPosFiniteValuesOf(@nospecialize fmt::Format) = nPosValuesOf(fmt) - nPosInfsOf(fmt)

"""Number of non-negative finite values (positive finite plus zero)."""
nNonNegFiniteValuesOf(@nospecialize fmt::Format) = nPosFiniteValuesOf(fmt) + nZerosOf(fmt)

#

"""Number of non-negative prenormal code points (zero + positive subnormals): `2^(P-1)`, or 1 when `P=1`."""
nNonNegPrenormalsOf(@nospecialize fmt::Format) = PrecisionOf(fmt) == 1 ? 1 : twopow(PrecisionOf(fmt) - 1)

"""Total number of prenormal code points (zero + all subnormals)."""
nPrenormalsOf(@nospecialize fmt::Format{is_unsigned,T}) where T = PrecisionOf(fmt) == 1 ? 1 : nNonNegPrenormalsOf(fmt)
nPrenormalsOf(@nospecialize fmt::Format{is_signed,T}) where T = PrecisionOf(fmt) == 1 ? 1 : nNonNegPrenormalsOf(fmt) + (nNonNegPrenormalsOf(fmt) - 1)

"""Number of positive subnormal values: `2^(P-1) - 1`, or 0 when `P=1`."""
nPosSubnormalsOf(@nospecialize fmt::Format) = PrecisionOf(fmt) == 1 ? 0 : nNonNegPrenormalsOf(fmt) - 1

"""Number of negative subnormal values. Always 0 for unsigned."""
nNegSubnormalsOf(@nospecialize fmt::Format{is_unsigned,T}) where T = 0
nNegSubnormalsOf(@nospecialize fmt::Format{is_signed,T}) where T = PrecisionOf(fmt) == 1 ? 0 : nPosSubnormalsOf(fmt)

"""Total number of subnormal values (positive + negative)."""
nSubnormalsOf(@nospecialize fmt::Format) = PrecisionOf(fmt) == 1 ? 0 : nPosSubnormalsOf(fmt) + nNegSubnormalsOf(fmt)

#

"""Total number of normal values (finite minus prenormals)."""
nNormalsOf(@nospecialize fmt::Format) = nFiniteValuesOf(fmt) - nPrenormalsOf(fmt)

"""Number of positive normal values."""
nPosNormalsOf(@nospecialize fmt::Format{is_unsigned,T}) where T = nNormalsOf(fmt)
nPosNormalsOf(@nospecialize fmt::Format{is_signed,T}) where T = nNormalsOf(fmt) ÷ 2

"""Number of negative normal values. Always 0 for unsigned."""
nNegNormalsOf(@nospecialize fmt::Format{is_unsigned,T}) where T = 0
nNegNormalsOf(@nospecialize fmt::Format{is_signed,T}) where T = nPosNormalsOf(fmt)

#

"""Number of normal binades: `2B - 1` where `B = ExponentBiasOf(fmt)`."""
nBinadesOf(@nospecialize fmt::Format) = 2 * ExponentBiasOf(fmt) - 1