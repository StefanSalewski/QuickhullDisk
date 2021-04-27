import relativeoperator, point2d, line2d
from math import PI, TAU, sqrt, hypot

type
  Circle2D* = object of RootObj
    m_CenterPt: Point2D
    m_Radius: float

# https://forum.nim-lang.org/t/7708
proc initCircle2D*(center: Point2D, r: float; T: typedesc[Circle2d] = Circle2D): T =
  result.m_CenterPt = center
  result.m_Radius   = r

proc initCircle2D*(pointX, pointY, r: float; T: typedesc[Circle2d] = Circle2D): T =
  result.m_CenterPt.set_point(pointX, pointY)
  result.m_Radius   = r

proc initCircle2D*(circle: Circle2D; T: typedesc[Circle2d] = Circle2D): T =
  result.m_CenterPt = circle.get_center_pt()
  result.m_Radius   = circle.get_radius()

proc get_center_pt*(c: Circle2D): Point2D =
  c.m_CenterPt

proc get_radius*(c: Circle2D): float =
  c.m_Radius

proc get_circle*(c: Circle2D): Circle2D =
  return c

proc set_center_pt*(c: var Circle2D; center: Point2D) =
  c.m_CenterPt = center

proc set_radius*(c: var Circle2D; r: float) =
  c.m_Radius   = r

proc set_circle*(c: var Circle2D;  circle: Circle2D) =
  c.m_CenterPt = circle.m_CenterPt
  c.m_Radius = circle.m_Radius

proc set_circle*(c: var Circle2D; center: Point2D; r: float) =
  c.m_CenterPt = center
  c.m_Radius = r

proc compute_area*(c: Circle2D): float =
  math.PI * c.m_Radius * c.m_Radius

proc compute_perimeter*(c: Circle2D): float =
  math.TAU * c.m_Radius

proc has_intersection_with*(c: Circle2D; lineSegment: Line2D): bool =
  return if LT(lineSegment.get_distance(c.m_CenterPt), c.m_Radius) == true: true else: false

proc is_intersected_with*(c: Circle2D; circle: Circle2D): bool =
  var center2center = initPoint2D(circle.m_CenterPt - c.m_CenterPt)
  var distance: float = center2center.magnitude()
  return not LT(c.m_Radius + circle.m_Radius, distance)

proc is_this_circle_boundary_intersected_with_another_circle_boundary*(c: Circle2D; circle: Circle2D): bool =
  var center2center = initPoint2D(circle.m_CenterPt - c.m_CenterPt)
  var distance = center2center.magnitude()
  return (GT(distance, ABS(c.m_Radius - circle.m_Radius)) and LT(distance, c.m_Radius + circle.m_Radius))

proc is_tangent_to*(c: Circle2D; circle: Circle2D): bool =
  var center2center = initPoint2D(circle.m_CenterPt - c.m_CenterPt)
  var distance: float = center2center.magnitude()
  return EQ(c.m_Radius + circle.m_Radius, distance) or EQ(ABS(c.m_Radius - circle.m_Radius), distance)

proc is_included_in*(c: Circle2D; circle: Circle2D): bool =
  var center2center = initPoint2D(circle.m_CenterPt - c.m_CenterPt)
  var distance: float = center2center.magnitude()
  return GT((circle.m_Radius - c.m_Radius), distance)

proc does_contain*(c: Circle2D; point: Point2D): bool =
  var distance: float = c.m_CenterPt.distance(point)
  return GT(c.m_Radius, distance)

proc contain*(c: Circle2D; circle: Circle2D; tolerance: float): bool =
  var center2center = initPoint2D(circle.m_CenterPt - c.m_CenterPt)
  var distance: float = center2center.magnitude()
  return  c.m_Radius >= (circle.m_Radius + distance) - tolerance

proc get_intersection_point*(cc: Circle2D; circle: Circle2D): (Point2D, Point2D) =
  #          *  <-- intersect point
  #         /|\
  #          / |  \
  #         /   |h   \
  #        /_c_|______\
  #     center  foot   circle.center
  #
  var center2center = initPoint2D(circle.m_CenterPt - cc.m_CenterPt)
  var distance: float = center2center.magnitude()
  var foot: Point2D
  if LT(cc.m_Radius + circle.m_Radius, distance):
    return (Point2DNil, Point2DNil)
  if GT(ABS(cc.m_Radius - circle.m_Radius), distance):
    return (Point2DNil, Point2DNil)
  if EQ(cc.m_Radius, circle.m_Radius) and EQ(distance, 0.0):
    return (Point2DNil, Point2DNil)
  var
    c: float = 0
    h: float = 0
  c = (distance * distance + cc.m_Radius * cc.m_Radius - circle.m_Radius * circle.m_Radius) / (2 * distance)
  if cc.m_Radius * cc.m_Radius - c * c > 0:
    h = math.sqrt(cc.m_Radius * cc.m_Radius - c * c)
  var temp: Point2D = initPoint2D(-1 * center2center.get_y(), center2center.get_x())
  foot = cc.m_CenterPt + c * (center2center.get_unit_vector())
  result[0] = foot + temp.get_unit_vector() * h
  result[1] = foot - temp.get_unit_vector() * h

#[
proc `=`*(c: Circle2D; circle: Circle2D): Circle2D =
  if c == addr(circle):
    return this[]
  c.m_CenterPt = circle.m_CenterPt
  c.m_Radius = circle.m_Radius
  return this[]
]#

proc `==`*(c: Circle2D; circle: Circle2D): bool =
  return c.m_CenterPt == circle.m_CenterPt and c.m_Radius == circle.m_Radius

proc distance*(c: Circle2D; point: Point2D): float =
  return c.m_CenterPt.distance(point) - c. m_Radius

proc distance*(c: Circle2D; circle: Circle2D): float =
  var dist: float = 0
  var distC2C: float = c.m_CenterPt.distance(circle.m_CenterPt)
  if GT(distC2C, (c.m_Radius + circle.m_Radius)):
    dist = distC2C - (c.m_Radius + circle.m_Radius)
  else:
    dist = (c.m_Radius + circle.m_Radius) - distC2C
  return dist
