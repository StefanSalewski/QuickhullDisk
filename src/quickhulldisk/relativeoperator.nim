import constset

proc ABS*(x: float): float =
  return if x > 0: x else: -x

proc EQ*(r1: float; r2: float; res: float = MATH_RES): bool =
  return ABS(r1 - r2) <= res

proc GT*(r1: float; r2: float; res: float = MATH_RES): bool =
  return r1 > r2 + res

proc LT*(r1: float; r2: float; res: float = MATH_RES): bool =
  return r1 < r2 - res

proc GE*(r1: float; r2: float; res: float = MATH_RES): bool =
  return not LT(r1, r2, res)

proc LE*(r1: float; r2: float; res: float = MATH_RES): bool =
  return not GT(r1, r2, res)

proc NE*(r1: float; r2: float; res: float = MATH_RES): bool =
  return not EQ(r1, r2, res)

proc ZERO*(r: float; res: float = MATH_RES): bool =
  return ABS(r) <= res

proc POS*(r: float; res: float = MATH_RES): bool =
  return r > res

proc NEG*(r: float; res: float = MATH_RES): bool =
  return r < -res

proc NZERO*(r: float; res: float = MATH_RES): bool =
  return not ZERO(r, res)

proc NPOS*(r: float; res: float = MATH_RES): bool =
  return not POS(r, res)

proc NNEG*(r: float; res: float = MATH_RES): bool =
  return not NEG(r, res)

#       Between Two Ordered Reals: (the real values are ordered, i.e., left <= right)
proc BTOR*(left: float; value: float; right: float; res: float = MATH_RES): bool =
  return LE(left, value, res) and LE(value, right, res)

proc BTORexclusive*(left: float; value: float; right: float; res: float = MATH_RES): bool =
  return LT(left, value, res) and LT(value, right, res)

#       Between Two Reals: (the real values are not ordered.)
proc BTR*(r1: float; r: float; r2: float; res: float = MATH_RES): bool =
  return (LE(r1, r, res) and LE(r, r2, res)) or (LE(r2, r, res) and LE(r, r1, res))

# ----- if r1 > r2 then 1 is returned
# ----- if r1 = r2 then 0 is returned
# ----- if r1 < r2 then -1is  returned
#
proc COMPARE*(r1, r2, res: float = MATH_RES): int =
  if ( ABS(r1 - r2) < res ):
    return 0
  elif ( r1 > r2 ):
    return 1
  else:
    return -1
