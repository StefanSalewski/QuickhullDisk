import point2d, relativeoperator

type
  Line2D* = object
    m_StartPt: Point2D
    m_EndPt: Point2D

proc initLine2D*(): Line2D =
  Line2D()

proc initLine2D*(sPt, ePt: Point2D): Line2D =
  result.m_StartPt = sPt
  result.m_EndPt = ePt

proc initLine2D*(line: Line2D): Line2D =
  result.m_StartPt = line.m_StartPt
  result.m_EndPt = line.m_StartPt

proc get_start_point*(l: Line2D): Point2D =
  return l.m_StartPt

proc get_end_point*(l: Line2D): Point2D =
  return l.m_EndPt

proc get_reversed_line2D*(l: Line2D): Line2D =
  return initLine2D(l.m_EndPt, l.m_StartPt)

proc set_start_pt*(l: var Line2D; pt: Point2D) =
  l.m_StartPt = pt

proc set_end_pt*(l: var Line2D; pt: Point2D) =
  l.m_EndPt = pt

proc set_line2d*(l: var Line2D; sPt, ePt: Point2D) =
  l.m_StartPt = sPt
  l.m_EndPt = ePt

proc set_line2d*(l: var Line2D; line: Line2D) =
  l.m_StartPt = line.m_StartPt
  l.m_EndPt = line.m_EndPt

proc get_length*(l: Line2D): float =
  var  length = initPoint2D(l.m_EndPt - l.m_StartPt)
  return length.magnitude()

proc get_distance*(l: Line2D; point: Point2D): float =
  var spPtVector: Point2D = initPoint2D(point - l.m_StartPt)
  var epPtVector: Point2D = initPoint2D(point - l.m_EndPt)
  var lineVector: Point2D = initPoint2D(l.m_EndPt - l.m_StartPt)
  var spPtLen: float = spPtVector.magnitude()
  var epPtLen: float = epPtVector.magnitude()
  var lineLen: float = lineVector.magnitude()
  if GE(epPtLen * epPtLen, (spPtLen * spPtLen) + (lineLen * lineLen)):
    return spPtLen
  elif GE((spPtLen * spPtLen), (epPtLen * epPtLen) + (lineLen * lineLen)):
    return epPtLen
  else:
    var output: float = ABS(lineVector * spPtVector) / lineLen
    return output

proc signed_distance*(l: Line2D; point: Point2D): float =
  var spPtVector: Point2D = initPoint2D(point - l.m_StartPt)
  var lineVector: Point2D = initPoint2D(l.m_EndPt - l.m_StartPt)
  var output: float = lineVector * spPtVector / lineVector.magnitude()
  return output

proc does_contain*(l: Line2D; point: Point2D): bool =
  return if (ZERO(l.signed_distance(point))): true else: false

proc evaluate_vector*(l: Line2D; ): Point2D =
  var output: Point2D = l.m_EndPt - l.m_StartPt
  return output

proc get_normal_vector*(l: Line2D; ): Point2D =
  var dirVec: Point2D = evaluate_vector(l)
  return initPoint2D(-dirVec.get_y(), dirVec.get_x())

# moved to module geometricfunction2d
#[
proc compute_intersection_with_line*(l: Line2D; lineSegment: Line2D; bTwoLinesAreParallel: bool): Point2D =
  return compute_intersection_between_two_lines(l.m_StartPt,
      l.m_EndPt, lineSegment.m_StartPt, lineSegment.m_EndPt, bTwoLinesAreParallel)
]#

proc make_perpendicular_line*(l: Line2D; passingPoint: Point2D): Line2D =
  var dirVec: Point2D = l.get_normal_vector()
  return initLine2D(passingPoint, passingPoint + dirVec)

proc `==`*(l: Line2D; line: Line2D): bool =
  (l.m_StartPt == line.m_StartPt and l.m_EndPt == line.m_EndPt) or (l.m_StartPt == line.m_EndPt and l.m_EndPt == line.m_StartPt)

#[
proc `=`*(l: Line2D; line: Line2D): Line2D =
  if this == addr(line):
    return this[]
  m_StartPt = line.m_StartPt
  m_EndPt = line.m_EndPt
  return this[]
]#
