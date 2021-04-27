import constset, relativeoperator, point2d, line2d, circle2d
from math import sqrt, hypot, pow

# The line implicit equation are as follows.
# coefficient[0] * x + coefficient[1] * y + coefficient[2] = 0
type
  LineImplicitEquation = object
    coefficient: array[3, float]

proc evaluateImpEquation(lie: LineImplicitEquation; valueForX, valueForY: float): float =
  return lie.coefficient[0] * valueForX + lie.coefficient[1] * valueForY + lie.coefficient[2]

type
  GeometricFunction2D* = object
    dummy: int

proc compute_intersection_between_two_lines*(gf: GeometricFunction2D; startPointOfOrientedLineSegment1: Point2D; endPointOfOrientedLineSegment1: Point2D;
    startPointOfOrientedLineSegment2: Point2D; endPointOfOrientedLineSegment2: Point2D; bTwoLinesAreParallel: var bool): Point2D =
  var parameterOfIntersectionForLineSeg1: float = DBL_MAX
  var lineSeg1Vec: Point2D = endPointOfOrientedLineSegment1 - startPointOfOrientedLineSegment1
  var lineSeg2Vec: Point2D = endPointOfOrientedLineSegment2 - startPointOfOrientedLineSegment2
  if ZERO(lineSeg1Vec * (lineSeg2Vec)):
    bTwoLinesAreParallel = true
    return initPoint2D(DBL_MAX, DBL_MAX)
  else:
    bTwoLinesAreParallel = false
    var startX1: float = startPointOfOrientedLineSegment1.get_x()
    var startY1: float = startPointOfOrientedLineSegment1.get_y()
    var endX1: float = endPointOfOrientedLineSegment1.get_x()
    var endY1: float = endPointOfOrientedLineSegment1.get_y()
    var startX2: float = startPointOfOrientedLineSegment2.get_x()
    var startY2: float = startPointOfOrientedLineSegment2.get_y()
    var endX2: float = endPointOfOrientedLineSegment2.get_x()
    var endY2: float = endPointOfOrientedLineSegment2.get_y()
    parameterOfIntersectionForLineSeg1 = (endX2 * startY1 - endX2 * startY2 - endY2 * startX1 + endY2 * startX2 + startX1 * startY2 - startX2 * startY1) /
        (endX1 * endY2 - endX1 * startY2 - endX2 * endY1 + endX2 * startY1 + endY1 * startX2 - endY2 * startX1 + startX1 * startY2 - startX2 * startY1)
    return (1.0 - parameterOfIntersectionForLineSeg1) * startPointOfOrientedLineSegment1 + parameterOfIntersectionForLineSeg1 * endPointOfOrientedLineSegment1

proc compute_intersection_between_two_circles*(gf: GeometricFunction2D; circle1: Circle2D; circle2: Circle2D; intersectionPt: var (Point2D, Point2D)): int =
  var radius1: float = circle1.get_radius()
  var radius2: float = circle2.get_radius()
  var center1: Point2D = circle1.get_center_pt()
  var center2: Point2D = circle2.get_center_pt()
  var vecCenter12: Point2D = center2 - center1
  var vecOrthogonalToVecCenter12: Point2D = initPoint2D(vecCenter12.get_y(), -vecCenter12.get_x())
  var radiu1_plus_radius2: float = radius1 + radius2
  var radiu1_minus_radius2: float = radius1 - radius2
  var normOfVecCenter12_squared: float = vecCenter12.magnitude_square()
  var parameterOfVecCenter12: float = (radiu1_plus_radius2 * radiu1_minus_radius2 / normOfVecCenter12_squared + 1.0) * 0.5
  var parameterOfVecOrthogonalToVecCenter12_squared: float = (-normOfVecCenter12_squared + radiu1_plus_radius2 * radiu1_plus_radius2) *
      (normOfVecCenter12_squared - radiu1_minus_radius2 * radiu1_minus_radius2) / (4.0 * normOfVecCenter12_squared * normOfVecCenter12_squared)
  if parameterOfVecOrthogonalToVecCenter12_squared < 0.0:
    return 0
  else:
    var parameterOfVecOrthogonalToVecCenter12: array[2, float] = [sqrt(parameterOfVecOrthogonalToVecCenter12_squared),
        -sqrt(parameterOfVecOrthogonalToVecCenter12_squared)]
    intersectionPt[0] = center1 + parameterOfVecCenter12 * vecCenter12 + parameterOfVecOrthogonalToVecCenter12[0] * vecOrthogonalToVecCenter12
    intersectionPt[1] = center1 + parameterOfVecCenter12 * vecCenter12 + parameterOfVecOrthogonalToVecCenter12[1] * vecOrthogonalToVecCenter12
    return 2

# this function returns tangent lines
proc make_exterior_tangent_lines_of_two_circles(gf: GeometricFunction2D; circle1, circle2: Circle2D; result1, result2: var LineImplicitEquation) =
  var
    w1: Circle2D
    w2: Circle2D
  # the radius of w2 is less than that of w1
  if circle1.get_radius() < circle2.get_radius():
    w1 = circle2
    w2 = circle1
  else:
    w1 = circle1
    w2 = circle2
  var c2cVector: Point2D = w1.get_center_pt() - w2.get_center_pt()
  var r: float = w1.get_radius() - w2.get_radius()
  var length: float = c2cVector.magnitude()
  var sine: float = r / length
  var cosine: float = sqrt(length * length - r * r) / length
  # rotate theta  /  -theta
  var normal1: Point2D = initPoint2D(c2cVector.get_x() * cosine - c2cVector.get_y() * sine , c2cVector.get_x() * sine + c2cVector.get_y() * cosine)
  var normal2: Point2D = initPoint2D(c2cVector.get_x() * cosine + c2cVector.get_y() * sine, c2cVector.get_y() * cosine - c2cVector.get_x() * sine)
  normal1 = normal1.get_unit_vector()
  normal2 = normal2.get_unit_vector()
  # rotate -PI/2  /  PI/2
  normal1.set_point(normal1.get_y(), -1 * normal1.get_x())
  normal2.set_point(-1 * normal2.get_y(), normal2.get_x())
  # set results
  result1.coefficient[0] = normal1.get_x()
  result1.coefficient[1] = normal1.get_y()
  result1.coefficient[2] = -(1 * normal1.get_x() * w2.get_center_pt.get_x()) - normal1.get_y() * w2.get_center_pt.get_y() + w2.get_radius()
  result2.coefficient[0] = normal2.get_x()
  result2.coefficient[1] = normal2.get_y()
  result2.coefficient[2] = -(1 * normal2.get_x() * w2.get_center_pt.get_x()) - normal2.get_y() * w2.get_center_pt.get_y() + w2.get_radius()

proc compute_tangent_point_between_line_and_circle(gf: GeometricFunction2D; circle: Circle2D; lineImplicitEq: LineImplicitEquation): Point2D =
  var tangentPt: Point2D
  var x: float = circle.get_center_pt.get_x()
  var y: float = circle.get_center_pt.get_y()
  var a: float = lineImplicitEq.coefficient[0]
  var b: float = lineImplicitEq.coefficient[1]
  # double c = lineImplicitEq.coefficient[2];
  var a2_b2: float = pow(a, 2) + pow(b, 2)
  if EQ(a2_b2, 1.0):
    tangentPt.set_x(-(a * lineImplicitEq.evaluateImpEquation(x, y)) + x)
    tangentPt.set_y(-(b * lineImplicitEq.evaluateImpEquation(x, y)) + y)
  else:
    tangentPt.set_x(-(a * lineImplicitEq.evaluateImpEquation(x, y) / a2_b2) + x)
    tangentPt.set_y(-(b * lineImplicitEq.evaluateImpEquation(x, y) / a2_b2) + y)
  return tangentPt

proc compute_two_oriented_tangent_line_segments_from_circle_c1_to_c2*(gf: GeometricFunction2D; circle1, circle2: Circle2D; twoTangentOrientedLineSegments: var (Line2D, Line2D)) =
  var implicitEqOfTangentLine1: LineImplicitEquation
  var implicitEqOfTangentLine2: LineImplicitEquation
  make_exterior_tangent_lines_of_two_circles(gf, circle1, circle2, implicitEqOfTangentLine1, implicitEqOfTangentLine2)
  var tangentPtOnDisk1: array[2, Point2D]
  var tangentPtOnDisk2: array[2, Point2D]
  if circle1.get_radius() != 0.0:
    tangentPtOnDisk1[0] = compute_tangent_point_between_line_and_circle(gf, circle1, implicitEqOfTangentLine1)
    tangentPtOnDisk1[1] = compute_tangent_point_between_line_and_circle(gf, circle1, implicitEqOfTangentLine2)
  else:
    tangentPtOnDisk1[0] = circle1.get_center_pt()
    tangentPtOnDisk1[1] = circle1.get_center_pt()
  if circle2.get_radius() != 0.0:
    tangentPtOnDisk2[0] = compute_tangent_point_between_line_and_circle(gf, circle2, implicitEqOfTangentLine1)
    tangentPtOnDisk2[1] = compute_tangent_point_between_line_and_circle(gf, circle2, implicitEqOfTangentLine2)
  else:
    tangentPtOnDisk2[0] = circle2.get_center_pt()
    tangentPtOnDisk2[1] = circle2.get_center_pt()
  twoTangentOrientedLineSegments[0] = initLine2D(tangentPtOnDisk1[0], tangentPtOnDisk2[0])
  twoTangentOrientedLineSegments[1] = initLine2D(tangentPtOnDisk1[1], tangentPtOnDisk2[1])

proc compute_CCW_oriented_tangent_line_segment_from_circle_c1_to_c2*(gf: GeometricFunction2D; circle1, circle2: Circle2D): Line2D =
  if circle1.get_radius() == 0.0 and circle2.get_radius() == 0.0:
    return initLine2D(circle1.get_center_pt(), circle2.get_center_pt())
  else:
    var twoTangentOrientedLineSegments: (Line2D, Line2D)
    compute_two_oriented_tangent_line_segments_from_circle_c1_to_c2(gf, circle1, circle2, twoTangentOrientedLineSegments)
    var sumOfSignedDistance: float = twoTangentOrientedLineSegments[0].signed_distance(circle1.get_center_pt()) +
        twoTangentOrientedLineSegments[0].signed_distance(circle2.get_center_pt())
    if sumOfSignedDistance > 0.0:
      return twoTangentOrientedLineSegments[0]
    else:
      return twoTangentOrientedLineSegments[1]

proc compute_CW_oriented_tangent_line_segment_from_circle_c1_to_c2*(gf: GeometricFunction2D; circle1, circle2: Circle2D): Line2D =
  if circle1.get_radius() == 0.0 and circle2.get_radius() == 0.0:
    return initLine2D(circle1.get_center_pt(), circle2.get_center_pt())
  else:
    var twoTangentOrientedLineSegments: (Line2D, Line2D)
    compute_two_oriented_tangent_line_segments_from_circle_c1_to_c2(gf, circle1, circle2, twoTangentOrientedLineSegments)
    var sumOfSignedDistance: float = twoTangentOrientedLineSegments[0].signed_distance(circle1.get_center_pt()) +
        twoTangentOrientedLineSegments[0].signed_distance(circle2.get_center_pt())
    if sumOfSignedDistance < 0.0:
      return twoTangentOrientedLineSegments[0]
    else:
      return twoTangentOrientedLineSegments[1]

proc compute_intersection_between_circle_and_line*(gf: GeometricFunction2D; circle: Circle2D; lineSegment: Line2D; intersection: var (Point2D, Point2D)): int =
  var startPoint: Point2D = lineSegment.get_start_point()
  var endPoint: Point2D = lineSegment.get_end_point()
  var center: Point2D = circle.get_center_pt()
  var radius: float = circle.get_radius()
  var startX: float = startPoint.get_x()
  var startY: float = startPoint.get_y()
  var endX: float = endPoint.get_x()
  var endY: float = endPoint.get_y()
  var centerX: float = center.get_x()
  var centerY: float = center.get_y()
  var coefficient: array[3, float] = [0.0, 0.0, 0.0]
  coefficient[0] = pow(startX - centerX, 2.0) + pow(startY - centerY, 2.0) - radius * radius
  coefficient[1] = 2.0 * (startX - centerX) * (-startX + endX) + 2.0 * (startY - centerY) * (-startY + endY)
  coefficient[2] = pow(-startX + endX, 2.0) + pow(-startY + endY, 2.0)
  var realSolution: array[2, float] = [DBL_MAX, DBL_MAX]
  if not ZERO(coefficient[2]): #  coefficient[ 2 ] == 0.0
    var p: float = coefficient[1] / (2.0 * coefficient[2])
    var q: float = coefficient[0] / coefficient[2]
    var determinant: float = p * p - q
    #  determinant >= 0
    if NNEG(determinant):
      realSolution[0] = -p + sqrt(determinant)
      realSolution[1] = -p - sqrt(determinant)
  else:
    realSolution[0] = -(coefficient[0] / coefficient[1])
  var numberOfIntersections: int = 0
  var i: int = 0
  var hhh: array[2, Point2D]
  while i < 2:
    if realSolution[i] != DBL_MAX:
      var r: float = realSolution[i]
      hhh[numberOfIntersections] = (1.0 - r) * startPoint + r * endPoint
      inc(numberOfIntersections)
    inc(i)
  intersection = (hhh[0], hhh[1])
  return numberOfIntersections
