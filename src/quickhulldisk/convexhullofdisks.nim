import constset, relativeoperator, point2d, line2d, circle2d, circulardisk2d, circulararc2d, geometricfunction2d
from math import TAU
from algorithm import sort

type
  ListCircularDisk2D* = seq[CircularDisk2D]

type
  ConvexHullOfDisks* = object of RootObj
    m_InputDisks*: ListCircularDisk2D # field export needed for quickhulldisk
    m_HullDisks*: ListCircularDisk2D
    
#[
struct compare_two_disk_for_unique_existence_in_set
{
    bool operator()(const CircularDisk2D& disk1, const CircularDisk2D& disk2) const
    {
        const int tolerance = _resNeg6;
        if (disk1.get_radius() < disk2.get_radius() - tolerance)
        {
            if (disk1.get_center_pt().get_x() < disk2.get_center_pt().get_x() - tolerance)
            {
                if (disk1.get_center_pt().get_y() < disk2.get_center_pt().get_y() - tolerance)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
        else
        {
            if (disk1.get_center_pt().get_x() < disk2.get_center_pt().get_x() - tolerance)
            {
                if (disk1.get_center_pt().get_y() < disk2.get_center_pt().get_y() - tolerance)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
    }
};
]#

#proc get_input_disks(chd: ConvexHullOfDisks): ListCircularDisk2D =
#  chd.m_InputDisks

#proc get_hull_disks(chd: ConvexHullOfDisks): ListCircularDisk2D =
#  chd.m_HullDisks

#  \brief Find extreme disks. Should be implemented in inherited class following the below parameter rules.
#
#  \param[in]  inputDisks    The input disks of convex hull problem.
#  \param[out] hullDisks     The output hull disks stored in CCW direction starting from a disk having high left extreme point.
#                            The same disk may occur many times but not consecutively;
#                             ex. Possible result:(d_1, d_2, d_1, d_3, d_4, d_1, ...),
#                                 impossible result: (d_1, d_2, d_2, d_2, d_4, d_1, ...) because d_2 occurs consecutively.
#                            The first disk must occurs at the last; ex. (d_1, ..., d_1).
#
#     virtual void find_hull_disks(const listCircularDisk2D& inputDisks, listCircularDisk2D& hullDisks) = 0;
#     virtual void clear_convexhull();
#

#[
void ConvexHullOfDisks::find_out_unique_disks_using_binary_search_tree(const list<CircularDisk2D>& inputDisks, list<CircularDisk2D>& uniqueDisks)
{
    set<CircularDisk2D, compare_two_disk_for_unique_existence_in_set> setOfDisks;
    setOfDisks.insert(inputDisks.begin(), inputDisks.end());

    uniqueDisks.insert(uniqueDisks.end(), setOfDisks.begin(), setOfDisks.end());
    uniqueDisks.sort(disk1_id_is_smaller_than_disk2_id);
}
]#

#[
inline bool ConvexHullOfDisks::disk1_id_is_smaller_than_disk2_id(CircularDisk2D& disk1, CircularDisk2D& disk2)
{
    return disk1.get_ID() < disk2.get_ID();
}
]#
proc disk1_id_is_smaller_than_disk2_id*(chd: ConvexHullOfDisks; disk1, disk2: CircularDisk2D): bool =
  return disk1.get_ID() < disk2.get_ID()

proc this_disk_is_a_member_of_expanded_non_negative_set_wrt_line*(chd: ConvexHullOfDisks;
    candidateDisk: CircularDisk2D; orientedLineSegment: Line2D; includingOnNegative: bool = false): bool  =
  if includingOnNegative:
    return GE(orientedLineSegment.signed_distance(candidateDisk.get_center_pt()), -candidateDisk.get_radius())
  else:
    return GT(orientedLineSegment.signed_distance(candidateDisk.get_center_pt()), -candidateDisk.get_radius())

proc this_disk_is_a_member_of_expanded_non_positive_set_wrt_line*(chd: ConvexHullOfDisks;
    candidateDisk: CircularDisk2D; orientedLineSegment: Line2D; includingOnPositive: bool = false): bool =
  if includingOnPositive:
    return LE(orientedLineSegment.signed_distance(candidateDisk.get_center_pt()), candidateDisk.get_radius())
  else:
    return LT(orientedLineSegment.signed_distance(candidateDisk.get_center_pt()), candidateDisk.get_radius())

proc this_disk_is_a_member_of_expanded_non_positive_set_wrt_line_segment*(chd: ConvexHullOfDisks;
    candidateDisk: CircularDisk2D; orientedLineSegment: Line2D; includingOnPositive: bool = false): bool  =
  var orientedLine: Line2D = orientedLineSegment
  var orthogonalLineAtStartPt: Line2D = orientedLine.make_perpendicular_line(orientedLineSegment.get_start_point())
  var orthogonalLineAtEndPt: Line2D = orientedLine.make_perpendicular_line(orientedLineSegment.get_end_point())
  var signedDistanceFromOriendtedLineToDiskCenterPt: float = orientedLine.signed_distance(candidateDisk.get_center_pt())
  # negative or on-negative
  if LE(signedDistanceFromOriendtedLineToDiskCenterPt, -candidateDisk.get_radius()):
    #  a candidate disk is fully contained between two lines
    return LT(orthogonalLineAtStartPt.signed_distance(candidateDisk.get_center_pt()), -candidateDisk.get_radius()) and
        GT(orthogonalLineAtEndPt.signed_distance(candidateDisk.get_center_pt()), candidateDisk.get_radius())
  # crossing 
  elif GT(signedDistanceFromOriendtedLineToDiskCenterPt, -candidateDisk.get_radius()) and
      LT(signedDistanceFromOriendtedLineToDiskCenterPt, candidateDisk.get_radius()):
    #  a center of disk is contained betweeen two lines
    return LT(orthogonalLineAtStartPt.signed_distance(candidateDisk.get_center_pt()), 0.0) and
        GT(orthogonalLineAtEndPt.signed_distance(candidateDisk.get_center_pt()), 0.0)
  # on-positive
  elif includingOnPositive and EQ(signedDistanceFromOriendtedLineToDiskCenterPt, candidateDisk.get_radius()):
    #  a center of disk is contained betweeen two lines
    return LT(orthogonalLineAtStartPt.signed_distance(candidateDisk.get_center_pt()), 0.0) and
        GT(orthogonalLineAtEndPt.signed_distance(candidateDisk.get_center_pt()), 0.0)
  else: # positive
    return false

proc disk1_is_smaller_than_disk2*(chd: ConvexHullOfDisks; disk1: CircularDisk2D; disk2: CircularDisk2D): bool  =
  return disk1.get_radius() < disk2.get_radius()

proc disk1_is_larger_than_disk2*(disk1: CircularDisk2D; disk2: CircularDisk2D): int =
  let res = disk1.get_radius() - disk2.get_radius()
  if res < 0:
    -1
  elif res > 0:
    1
  else:
    0

proc get_pointers_of_disks*(chd: ConvexHullOfDisks; disks: ListCircularDisk2D; diskPointers: var ListCircularDisk2D) =
  diskPointers = disks

proc get_input_disks*(chd: var ConvexHullOfDisks; inputDisks: var ListCircularDisk2D) =
  get_pointers_of_disks(chd, chd.m_InputDisks, inputDisks)

proc copy_disks_from_its_pointers*(chd: ConvexHullOfDisks;  diskPointers: ListCircularDisk2D; disks: var ListCircularDisk2D) =
  disks = diskPointers

proc get_hull_disks*(chd: ConvexHullOfDisks; hullDisks: var ListCircularDisk2D) =
  copy_disks_from_its_pointers(chd, chd.m_HullDisks, hullDisks)

proc compute_a_CCW_oriented_tangent_line_from_disk_d1_to_d2*(chd: ConvexHullOfDisks; disk1: CircularDisk2D; disk2: CircularDisk2D): Line2D =
  var dummy: GeometricFunction2D
  return compute_CCW_oriented_tangent_line_segment_from_circle_c1_to_c2(dummy, disk1, disk2)

proc clear_convexhull*(chd: var ConvexHullOfDisks) =
  chd.m_InputDisks.setLen(0)
  chd.m_HullDisks.setLen(0)

#  \brief Extract convex hull boundary from the input hull disks.
#
#  \details It assumes the input 'hullDisks' are hull disks. If not, you cannot get a correct solution.
#
#  \param[in]  hullDisks           The input hull disks of convex hull. ex. ( d_1, d_2, ..., d_k ).
#  \param[out] convexHullBoundary  The output convex hull boundary consist of pairs of arc(a_i) and line(t_i) so that the boundary is stored as ( (a_1, tl_1), (a_2, tl_2), ..., (a_k, tl_k) )
#                                   where arc is on boundary of disk (d_i) and line is tangent to two consecutive hull disks (d_i, d_i+1).
#                                  The start point of line 'tl_i' equals to the end point of arc 'a_i' and the end point of line 'tl_i' equals to the start point of arc 'a_i+1'.
#                                  Notice that the end point of 'tl_k' is equals to the start point of 'a_1'.
#                                  If the size of hullDisks is 1, it returns an circle arc and unmeaningful line having two end points as (x, y) = (DBL_MAX, DBL_MAX).
#  \param[return] length of convex hull boundary
#
proc extract_convexhull_boundary*(chd: ConvexHullOfDisks; 
    hullDisks: ListCircularDisk2D; boundaryArcNLines: var seq[(CircularArc2D, Line2D)]): float =
  if hullDisks.len == 0:
    return 0.0
  elif hullDisks.len == 1:
    var arc: CircularArc2D = initCircularArc2D(hullDisks[0])
    var tangentLine: Line2D = initLine2D(initPoint2D(DBL_MAX, DBL_MAX), initPoint2D(DBL_MAX, DBL_MAX))
    boundaryArcNLines.add((arc, tangentLine))
    return math.TAU * arc.get_radius()
  var tangentLineSegments: seq[Line2D]
  var geometryOfHullDisks: seq[CircularDisk2D]
  let numOfTangentLines: int = hullDisks.len - 1
  tangentLineSegments.setLen(numOfTangentLines)
  geometryOfHullDisks.setLen(numOfTangentLines)
  for idx in 0 .. hullDisks.high:
    var it_Disk = hullDisks[idx]
    var currDisk: CircularDisk2D = it_Disk
    geometryOfHullDisks[idx] = currDisk
    var nextDisk: CircularDisk2D = hullDisks[idx + 1]
    var orientedTangentLineSegment: Line2D = compute_a_CCW_oriented_tangent_line_from_disk_d1_to_d2(chd, currDisk, nextDisk)
    tangentLineSegments[idx] = orientedTangentLineSegment
  # adjusting
  var i: int = 0
  while i < numOfTangentLines:
    var currLineSegment: Line2D = tangentLineSegments[i mod numOfTangentLines]
    var nextLineSegment: Line2D = tangentLineSegments[(i + 1) mod numOfTangentLines]
    if EQ(currLineSegment.get_end_point().get_x(), nextLineSegment.get_start_point().get_x(), resNeg5) and
        EQ(currLineSegment.get_end_point().get_y(), nextLineSegment.get_start_point().get_y(), resNeg5):
      nextLineSegment.set_start_pt(currLineSegment.get_end_point())
    inc(i)
  i = 0
  while i < numOfTangentLines:
    var startPtOfArc: Point2D = tangentLineSegments[(i - 1 + numOfTangentLines) mod numOfTangentLines].get_end_point()
    var endPtOfArc: Point2D = tangentLineSegments[i].get_start_point()
    var currArc: CircularArc2D = initCircularArc2D(geometryOfHullDisks[i], startPtOfArc, endPtOfArc)
    var tangentLineSegmentNextToCurrArc: Line2D = tangentLineSegments[i]
    boundaryArcNLines.add((currArc, tangentLineSegmentNextToCurrArc))
    inc(i)
  var lengthOfHullBoundary: float = 0.0
  for  it_ArcNLine in boundaryArcNLines:
    let arcNLine: (CircularArc2D, Line2D) = it_ArcNLine
    var arc: CircularArc2D = arcNLine[0]
    var line: Line2D = arcNLine[1]
    # lengthOfHullBoundary += arc.angle_btw_start_and_end_pts() * arc.get_radius()
    lengthOfHullBoundary += arc.caculate_arc_length()
    lengthOfHullBoundary += line.get_length()
  return lengthOfHullBoundary

proc sort_disks_by_their_radius_in_decreasing_order*(chd: ConvexHullOfDisks; 
    inputDisks: ListCircularDisk2D; sortedDisks: var ListCircularDisk2D) =
  sortedDisks = inputDisks
  sortedDisks.sort(disk1_is_larger_than_disk2)
  if sortedDisks.len > 1:
    assert sortedDisks[0].get_radius() > sortedDisks[1].get_radius()

proc find_high_left_extreme_point_N_its_disk*(chd: ConvexHullOfDisks; inputDisks: ListCircularDisk2D; highLeftExtremePt: var Point2D;
    diskHavingHighLeftExtremePt: var CircularDisk2D) =
  diskHavingHighLeftExtremePt = inputDisks[0]
  var leftMostX: float = diskHavingHighLeftExtremePt.get_center_pt().get_x() - diskHavingHighLeftExtremePt.get_radius()
  for it_Disk in inputDisks[1 .. ^1]: # caution: start at index 1 ???
    var currDisk: CircularDisk2D = it_Disk
    var leftMostXOfCurrDisk: float = currDisk.get_center_pt().get_x() - currDisk.get_radius()
    if leftMostXOfCurrDisk < leftMostX:
      leftMostX = leftMostXOfCurrDisk
      diskHavingHighLeftExtremePt = currDisk
    elif leftMostXOfCurrDisk == leftMostX:
      if currDisk.get_center_pt().get_y() > diskHavingHighLeftExtremePt.get_center_pt().get_y():
        diskHavingHighLeftExtremePt = currDisk
      elif currDisk.get_center_pt().get_y() == diskHavingHighLeftExtremePt.get_center_pt().get_y():
        if currDisk.get_radius() > diskHavingHighLeftExtremePt.get_radius():
          diskHavingHighLeftExtremePt = currDisk
  highLeftExtremePt = initPoint2D(leftMostX, diskHavingHighLeftExtremePt.get_center_pt().get_y())

proc rearrange_hull_disks_for_output_format_by_begining_with_disk_having_high_left_extreme_point*(chd: ConvexHullOfDisks; hullDisks: var ListCircularDisk2D) =
  var diskHavingHighLeftExtremePt: CircularDisk2D 
  var highLeftExtremePt: Point2D
  find_high_left_extreme_point_N_its_disk(chd, hullDisks, highLeftExtremePt, diskHavingHighLeftExtremePt)
  if hullDisks[0] != diskHavingHighLeftExtremePt: # is this ever used?
    hullDisks.setLen(hullDisks.high) # Because the first disk and the last is same for output format.
    while true: # really slow for a seq :-(
      var currDisk: CircularDisk2D = hullDisks[0]
      hullDisks.delete(0)
      hullDisks.add(currDisk)
      if not (hullDisks[0] != diskHavingHighLeftExtremePt):
        break
    hullDisks.add(hullDisks[0]) # Because the first disk and the last is same for output format.

proc find_low_right_extreme_point_N_its_disk*(chd: ConvexHullOfDisks;  inputDisks: ListCircularDisk2D; lowRightExtremePt: var Point2D;
    diskHavingLowRightExtremePt: var CircularDisk2D) =
  diskHavingLowRightExtremePt = inputDisks[0]
  var rightMostX: float = diskHavingLowRightExtremePt.get_center_pt().get_x() + diskHavingLowRightExtremePt.get_radius()
  for it_Disk in inputDisks[1 .. ^1]:
    var currDisk: CircularDisk2D = it_Disk
    var rightMostXOfCurrDisk: float = currDisk.get_center_pt().get_x() + currDisk.get_radius()
    if rightMostXOfCurrDisk > rightMostX:
      rightMostX = rightMostXOfCurrDisk
      diskHavingLowRightExtremePt = currDisk
    elif rightMostXOfCurrDisk == rightMostX:
      if currDisk.get_center_pt().get_y() < diskHavingLowRightExtremePt.get_center_pt().get_y():
        diskHavingLowRightExtremePt = currDisk
      elif currDisk.get_center_pt().get_y() == diskHavingLowRightExtremePt.get_center_pt().get_y():
        if currDisk.get_radius() > diskHavingLowRightExtremePt.get_radius():
          diskHavingLowRightExtremePt = currDisk
  lowRightExtremePt = initPoint2D(rightMostX, diskHavingLowRightExtremePt.get_center_pt().get_y())

# 265 lines

