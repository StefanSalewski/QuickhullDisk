import constset, point2d, circle2d, circulardisk2d
from math import TAU, sin, cos

type
  CircularArc2D* = object of CircularDisk2D
    m_StartPoint: Point2D
    m_EndPoint: Point2D

# https://forum.nim-lang.org/t/7708
proc initCircularArc2D*(T:typedesc[CircularArc2D] = CircularArc2D): T =
  result = initCircularDisk2D(T)
  result.m_StartPoint = initPoint2D(DBL_MAX, DBL_MAX)
  result.m_EndPoint = initPoint2D(DBL_MAX, DBL_MAX)

proc initCircularArc2D*(disk: CircularDisk2D; T:typedesc[CircularArc2D] = CircularArc2D): T =
  result = initCircularDisk2D(disk, T)
  result.m_StartPoint = initPoint2D(DBL_MAX, DBL_MAX)
  result.m_EndPoint = initPoint2D(DBL_MAX, DBL_MAX)

proc initCircularArc2D*(disk: CircularDisk2D; startPoint, endPoint: Point2D; T:typedesc[CircularArc2D] = CircularArc2D): T =
  result = initCircularDisk2D(disk, T)
  result.m_StartPoint = initPoint2D(startPoint)
  result.m_EndPoint = initPoint2D(endPoint)

proc initCircularArc2D*(arc: CircularArc2D; T:typedesc[CircularArc2D] = CircularArc2D): T =
  result = initCircularDisk2D(arc, T)
  result.m_StartPoint = arc.m_StartPoint
  result.m_EndPoint = arc.m_EndPoint

proc get_start_point*(ca: CircularArc2D): Point2D =
  ca.m_StartPoint

proc get_end_point*(ca: CircularArc2D): Point2D =
  ca.m_EndPoint

proc set_start_point*(ca: var CircularArc2D; startPoint: Point2D) =
  ca.m_StartPoint = startPoint

proc set_end_point*(ca: var CircularArc2D; endPoint: Point2D) =
  ca.m_EndPoint = endPoint

proc `==`*(ca, arc: CircularArc2D): bool =
  return CircularDisk2D(ca) == CircularDisk2D(arc) and ca.m_StartPoint == arc.m_StartPoint and ca.m_EndPoint == arc.m_EndPoint

proc angle_btw_start_and_end_pts*(ca: CircularArc2D): float =
  if ca.m_StartPoint.get_x() == DBL_MAX or ca.m_StartPoint.get_y() == DBL_MAX:
    return math.TAU
  elif ca.m_StartPoint == ca.m_EndPoint:
    return 0.0
  else:
    var center: Point2D = ca.get_center_pt()
    var vecCenterToStart: Point2D = ca.m_StartPoint - center
    var vecCenterToEnd: Point2D = ca.m_EndPoint - center
    var angleStartToEnd: float = angle_from_vec1_to_vec2(vecCenterToStart, vecCenterToEnd)
    return angleStartToEnd

#[
inline double CircularArc2D::caculate_arc_length() const
{
    return angle_btw_start_and_end_pts() * get_radius();
}
]#

proc caculate_arc_length*(ca: CircularArc2D): float =
  return ca.angle_btw_start_and_end_pts() * ca.get_radius()

proc evaluate_points_on_arc_given_resolution*(ca: CircularArc2D; circleResolution: int; pointsOnArc: var seq[Point2D]) =
  var center: Point2D = ca.get_center_pt()
  var radius: float = ca.get_radius()
  if ca.m_StartPoint == ca.m_EndPoint:
    var i: int = 0
    while i < circleResolution:
      var  angle: float = math.TAU * i.float / circleResolution.float
      var point: Point2D = initPoint2D(center.get_x() + radius * cos(angle), center.get_y() + radius * sin(angle))
      pointsOnArc.add(point)
      inc(i)
    pointsOnArc.add(pointsOnArc[0])
  else:
    var startVec: Point2D = ca.m_StartPoint - center
    var endVec: Point2D = ca.m_EndPoint - center
    var uVector: Point2D = startVec.get_unit_vector()
    var vVector: Point2D = initPoint2D(-uVector.get_y(), uVector.get_x())
    var arcAngle: float = angle_from_vec1_to_vec2(startVec, endVec)
    var angle: float = 0.0
    var i: int = 0
    while true:
      var point: Point2D = center + radius * (cos(angle) * uVector + sin(angle) * vVector)
      pointsOnArc.add(point)
      inc(i)
      angle = math.TAU * i.float / circleResolution.float
      if not (angle < arcAngle):
        break
    pointsOnArc.add(ca.m_EndPoint)
