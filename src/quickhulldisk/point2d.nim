from math import TAU, hypot, arctan2
from relativeoperator import EQ

type
  Point2D* = object
    m_X: float
    m_Y: float

const
  Point2DNil* = Point2D(m_X: system.INF, m_Y: system.INF)

proc initPoint2D*(): Point2D =
  Point2D()

proc initPoint2D*(point: Point2D): Point2D =
  result.m_x = point.m_x
  result.m_y = point.m_y

proc initPoint2D*(px, py: float): Point2D =
  result.m_x = px
  result.m_y = py

proc get_x*(p: Point2d): float = p.m_X

proc get_y*(p: Point2d): float = p.m_Y

proc set_x*(p: var Point2D; px: float) =
  p.m_X = px

proc set_y*(p: var Point2D; py: float) =
  p.m_Y = py

proc set_point*(p: var Point2D; point: Point2D) =
  p.m_X = point.m_X
  p.m_Y = point.m_Y

proc set_point*(p: var Point2D; px, py: float) =
  p.m_X = px
  p.m_Y = py

proc magnitude*(p: Point2D): float =
  return hypot(p.m_X, p.m_Y)

proc magnitude_square*(p: Point2D): float =
  return p.m_X * p.m_X + p.m_Y * p.m_Y

proc distance*(p: Point2D; point: Point2D): float =
  return hypot(p.m_X - point.m_X, p.m_Y - point.m_Y)

#[ not allowed in Nim ?
proc `=`*(p: var Point2D; point: Point2D): Point2D =
  p.m_X = point.m_X
  p.m_Y = point.m_Y
  return p
]#

proc `==`*(p: Point2D; point: Point2D): bool =
  EQ(p.m_X, point.m_X) and EQ(p.m_Y, point.m_Y)

# we get that for free in Nim
#proc `!=`*(p: Point2D; point: Point2D): bool =
#  return not (p == point)

proc `+`*(p: Point2D; point: Point2D): Point2D =
  result.m_X = p.m_X + point.m_X
  result.m_Y = p.m_Y + point.m_Y
  #initPoint2D(p.m_X + point.m_X, p.m_Y + point.m_Y)

proc `+=`*(p: var Point2D; point: Point2D): Point2D =
  p.m_X += point.m_X
  p.m_Y += point.m_Y
  return p

proc `-`*(p: Point2D; point: Point2D): Point2D =
  return initPoint2D(p.m_X - point.m_X, p.m_Y - point.m_Y)

proc `-=`*(p: var Point2D; point: Point2D): Point2D =
  p.m_X -= point.m_X
  p.m_Y -= point.m_Y
  return p

proc `*`*(p: Point2D; point: Point2D): float = # cross product
  return p.m_X * point.m_Y - p.m_Y * point.m_X

proc `*`*(p: Point2D; n: float): Point2D =
  return initPoint2D(n * p.m_X, n * p.m_Y)

proc `/`*(p: Point2D; n: float): Point2D =
  return initPoint2D(p.m_X / n, p.m_Y / n)

proc `=====`*(p: var Point2D; n: float): Point2D = # plain `=` not allowed!
  p.m_X = n
  p.m_Y = n
  return p

proc `%`*(p: Point2D; point: Point2D): float = # dot product
  return p.m_X * point.get_x() + p.m_Y * point.get_y()

proc `*`*(n: float; point: Point2D): Point2D =
  return point * n

proc `-`*(point: Point2D): Point2D =
  return initPoint2D(-point.m_X, -point.m_Y)

proc get_unit_vector*(p: Point2D): Point2D =
  return p / magnitude(p)

# should be arctan2()
proc no_angle_from_vec1_to_vec2*(vector1: Point2D; vector2: Point2D): float =
  var vec1: Point2D = vector1.get_unit_vector()
  var vec2: Point2D = vector2.get_unit_vector()
  var pt2pt: Point2D = vec1 - vec2
  var length: float = pt2pt.magnitude()
  var cosine: float = (2.0 - length * length) / (2.0)
  if cosine > 1.0:
    cosine = 1.0
  elif cosine < -1.0:
    cosine = -1.0
  if vec1 * vec2 > 0.0:
    return math.arccos(cosine)
  else:
    return math.TAU - math.arccos(cosine)

proc angle_from_vec1_to_vec2*(p1: Point2D; p2: Point2D): float =
  result = arctan2(p2.m_Y, p2.m_x) - arctan2(p1.m_Y, p1.m_X)
  if result < 0:
    result += TAU
  assert (result - no_angle_from_vec1_to_vec2(p1, p2)).abs < 1e-13

when isMainModule:
  import random
  for i in 0 .. 30:
    var p1 = initPoint2D(rand(99.0), rand(50.0))
    var p2 = initPoint2D(rand(99.0), rand(50.0))
    echo angle_from_vec1_to_vec2(p1, p2)
    echo no_angle_from_vec1_to_vec2(p1, p2)

