import constset, relativeoperator, convexhullofdisks,  circulardisk2d, point2d, circle2d, line2d
from times import getTime, nanosecond
from sequtils import deduplicate
import random

# dividing option when disks are cotangent to one common line
type
  CotangentDisksDividingOption* = enum
    RANDOM_DIVIDE, BISECTION_DIVIDE

# configuration type of sliver triangle filter
type
  SliverConfiguration* = enum
    SLIVER_CASE_A, SLIVER_CASE_B, SLIVER_CASE_C1, SLIVER_CASE_C2

# checking computation time of each step
type
  QHullDiskAlgorithmStep* = enum
    QH_TIME_QUICKHULLDISK, QH_TIME_TOTAL, QH_TIME_SIZE

# representation for input of FindHull at each iteration
type
  InputForFindHull* = object
    disks*: ListCircularDisk2D
    preApexDisk*: CircularDisk2D
    postApexDisk*: CircularDisk2D
    hullPointOfPreApexDisk*: Point2D
    hullPointOfPostApexDisk*: Point2D

type
  StackInputForFindHull = seq[InputForFindHull]
  
type
  QuickhullDisk* = object of ConvexHullOfDisks
    m_CotangentDisksDividingOption: CotangentDisksDividingOption

proc find_the_fartest_point_of_disk_from_this_line(disk: Circle2D; orientedLine: Line2D): Point2D =
  var negativeDirection: Point2D = -orientedLine.get_normal_vector()
  var farthestPoint: Point2D = negativeDirection.get_unit_vector() * disk.get_radius() + disk.get_center_pt()
  return farthestPoint

proc prepare_and_insert_input_data_for_finding_hull_to_stack(
    disks: ListCircularDisk2D; preApexDisk: CircularDisk2D;
    postApexDisk: CircularDisk2D; hullPointOfPreApexDisk: Point2D;
    hullPointOfPostApexDisk: Point2D; stackForFindingHull: var StackInputForFindHull) =
  var inputForFindHull = InputForFindHull(disks: disks, preApexDisk: preApexDisk, postApexDisk: postApexDisk, hullPointOfPreApexDisk: hullPointOfPreApexDisk,
    hullPointOfPostApexDisk: hullPointOfPostApexDisk)
  stackForFindingHull.add(inputForFindHull)

proc triangle_filter_is_sliver(numOfInputDisks: int;
    numOfDisksOnFrontEdge: int; numOfDisksOnBackEdge: int;
    preApexDisk_dp: CircularDisk2D; postApexDisk_dq: CircularDisk2D): bool =
  return (preApexDisk_dp != postApexDisk_dq) and ((numOfDisksOnFrontEdge == numOfInputDisks and
    numOfDisksOnBackEdge == 1) or (numOfDisksOnFrontEdge == 1 and numOfDisksOnBackEdge == numOfInputDisks))

proc find_high_left_N_low_right_extreme_points_and_their_disks_to_divide_disk_set(qhd: QuickhullDisk; 
    inputDisks_D: ListCircularDisk2D; highLeftExtremePoint_p: var Point2D;
    lowRightExtremePoint_q: var Point2D;
    diskHavingHighLeftExtremePoint_p: var CircularDisk2D;
    diskHavingLowRightExtremePoint_q: var CircularDisk2D) =
  find_high_left_extreme_point_N_its_disk(qhd, inputDisks_D, highLeftExtremePoint_p, diskHavingHighLeftExtremePoint_p)
  find_low_right_extreme_point_N_its_disk(qhd, inputDisks_D, lowRightExtremePoint_q, diskHavingLowRightExtremePoint_q)

proc find_the_highest_triangle_apex_and_the_apex_disk_wrt_this_oriented_line_segment*(qhd: QuickhullDisk; 
    disks_D: ListCircularDisk2D; orientedLineSegment_pq: Line2D;
    apexNDiskPairs: var seq[(Point2D, CircularDisk2D)]; orientedLineIsNonNegativeSupportOf_dp_N_dq: bool = false;
    preApexDisk_dp: CircularDisk2D = CircularDisk2DNil;
    postApexDisk_dq: CircularDisk2D = CircularDisk2DNil) =
  var largestMaxPerpendicularDistanceFromLineToBoundaryOfDiskAmongDisks: float = 0.0
  for it_Disk in disks_D:
    var currDisk: CircularDisk2D = it_Disk
    if orientedLineIsNonNegativeSupportOf_dp_N_dq and (currDisk == preApexDisk_dp or currDisk == postApexDisk_dq):
      #  avoid numerical decision for preApexDisk_dp and preApexDisk_dp
      continue
    var maxPerpendicularDistanceFromLineToBoundaryOfCurrDisk: float = currDisk.get_radius() -
      orientedLineSegment_pq.signed_distance(currDisk.get_center_pt())
    if GT(maxPerpendicularDistanceFromLineToBoundaryOfCurrDisk, largestMaxPerpendicularDistanceFromLineToBoundaryOfDiskAmongDisks):
      largestMaxPerpendicularDistanceFromLineToBoundaryOfDiskAmongDisks = maxPerpendicularDistanceFromLineToBoundaryOfCurrDisk
      apexNDiskPairs.setLen(0)
      var farthestPointOnCurrDiskTouchingTangentLine: Point2D = find_the_fartest_point_of_disk_from_this_line(currDisk, orientedLineSegment_pq)
      apexNDiskPairs.add((farthestPointOnCurrDiskTouchingTangentLine, currDisk))
    elif EQ(maxPerpendicularDistanceFromLineToBoundaryOfCurrDisk, largestMaxPerpendicularDistanceFromLineToBoundaryOfDiskAmongDisks):
      var farthestPointOnCurrDiskTouchingTangentLine: Point2D = find_the_fartest_point_of_disk_from_this_line(currDisk, orientedLineSegment_pq)
      apexNDiskPairs.add((farthestPointOnCurrDiskTouchingTangentLine,currDisk))
  if orientedLineIsNonNegativeSupportOf_dp_N_dq and ZERO(largestMaxPerpendicularDistanceFromLineToBoundaryOfDiskAmongDisks):
    if preApexDisk_dp != CircularDisk2DNil:
      apexNDiskPairs.add((orientedLineSegment_pq.get_start_point(), preApexDisk_dp))
    if postApexDisk_dq != CircularDisk2DNil:
      apexNDiskPairs.add((orientedLineSegment_pq.get_end_point(), postApexDisk_dq))

proc pick_one_as_triangle_apex_and_apex_disk(qhd: QuickhullDisk; 
    apexNDiskPairs: seq[(Point2D, CircularDisk2D)]; triangleApex_x: var Point2D;
    apexDisk_dx: var CircularDisk2D;
    cotangentDisksDividingOption: CotangentDisksDividingOption) =
  if apexNDiskPairs.len >= 2:
    var generatedNumber: int = 0
    case cotangentDisksDividingOption
    of RANDOM_DIVIDE:
      generatedNumber = rand(apexNDiskPairs.len)
    of BISECTION_DIVIDE:
      generatedNumber = apexNDiskPairs.len div 2
    var it_ApexNItsDisk: (Point2D, CircularDisk2D)
    for i, el in apexNDiskPairs:
      it_ApexNItsDisk = el
      if i == generatedNumber:
        break
    triangleApex_x = it_ApexNItsDisk[0]
    apexDisk_dx = it_ApexNItsDisk[1]
  else: #  (apexNDiskPairs.size() == 1)
    triangleApex_x = apexNDiskPairs[0][0]
    apexDisk_dx = apexNDiskPairs[0][1]

proc find_largest_apex_disk_containing_this_apex_disk_selected_from_candidates(qhd: QuickhullDisk; 
    candidateTriangleApexNDiskPairs: seq[(Point2D, CircularDisk2D)];
    selectedApexDisk: CircularDisk2D;
    containdDisksInOthers: var ListCircularDisk2D): CircularDisk2D =
  var largestApexDiskContainingCandidateApexDisk: CircularDisk2D = selectedApexDisk
  for it_ApexNItsDisk in candidateTriangleApexNDiskPairs:
    var currApexDisk: CircularDisk2D = it_ApexNItsDisk[1]
    if currApexDisk == largestApexDiskContainingCandidateApexDisk:
      continue
    if currApexDisk.contain(largestApexDiskContainingCandidateApexDisk, resNeg6):
      containdDisksInOthers.add(largestApexDiskContainingCandidateApexDisk)
      largestApexDiskContainingCandidateApexDisk = currApexDisk
    elif largestApexDiskContainingCandidateApexDisk.contain(currApexDisk, resNeg6):
      containdDisksInOthers.add(currApexDisk)
  return largestApexDiskContainingCandidateApexDisk

proc remove_contained_disks_from_input_disks*(qhd: QuickhullDisk; 
    containdDisksInOthers: ListCircularDisk2D; disks_D: var ListCircularDisk2D) =
  for it_ContainedDisk in containdDisksInOthers:
    let h = disks_D.find(it_ContainedDisk)
    if h >= 0:
      disks_D.delete(h)

proc pick_one_as_triangle_apex_and_apex_disk_among_disks_with_identical_height_and_remove_disks_contained_in_others_from_input_disks_if_exist*(qhd: QuickhullDisk; 
    apexNDiskPairs: seq[(Point2D, CircularDisk2D)]; triangleApex_x: var Point2D;
    apexDisk_dx: var CircularDisk2D; disks_D: var seq[CircularDisk2D];
    cotangentDisksDividingOption: CotangentDisksDividingOption) =
  var candidateOfTriangleApex: Point2D
  var candidateOfApexDisk: CircularDisk2D = CircularDisk2DNil
  pick_one_as_triangle_apex_and_apex_disk(qhd, apexNDiskPairs, candidateOfTriangleApex, candidateOfApexDisk, cotangentDisksDividingOption)
  var containdDisksInOthers: ListCircularDisk2D
  triangleApex_x = candidateOfTriangleApex
  apexDisk_dx = find_largest_apex_disk_containing_this_apex_disk_selected_from_candidates(qhd, apexNDiskPairs, candidateOfApexDisk, containdDisksInOthers)
  if containdDisksInOthers.len > 0:
    remove_contained_disks_from_input_disks(qhd, containdDisksInOthers, disks_D)

proc find_expanded_non_positive_disks_wrt_oriented_line_segment*(qhd: QuickhullDisk; 
    disks: ListCircularDisk2D; orientedLineSegmentOfTwoPtsOn_d1_N_d2: Line2D;
    disk_d1: CircularDisk2D; disk_d2: CircularDisk2D;
    outputDisks: var ListCircularDisk2D; includingOnPositive: bool = false) =
  if orientedLineSegmentOfTwoPtsOn_d1_N_d2.get_start_point() != orientedLineSegmentOfTwoPtsOn_d1_N_d2.get_end_point():
    for it_Disk in disks:
      var currDisk: CircularDisk2D = it_Disk
      if currDisk == disk_d1 or currDisk == disk_d2:
        #  do not want numerical decision for d1 and d2, so pass those and insert to list at the end of this function.
        continue
      else:
        if this_disk_is_a_member_of_expanded_non_positive_set_wrt_line_segment(qhd, currDisk, orientedLineSegmentOfTwoPtsOn_d1_N_d2, includingOnPositive):
          outputDisks.add(currDisk)
  if disk_d1 == disk_d2:
    outputDisks.add(disk_d1)
  else:
    outputDisks.add(disk_d1)
    outputDisks.add(disk_d2)

proc remove_pre_and_post_apex_disks_and_contained_disks_in_one_of_the_two_from_candidate_apex_disks_if_exist*(qhd: QuickhullDisk; 
    preApexDisk_dp: CircularDisk2D; postApexDisk_dq: CircularDisk2D;
    candidateTriangleApexNDiskPairs: var seq[(Point2D, CircularDisk2D)];
    containedDisksInPreOrPostApexDisk: var ListCircularDisk2D) =
  var i = 0
  var j = 0
  while i < candidateTriangleApexNDiskPairs.len:
    var candidateApexDisk: CircularDisk2D = candidateTriangleApexNDiskPairs[i][1]
    if candidateApexDisk == preApexDisk_dp or candidateApexDisk == postApexDisk_dq:
      inc(i)
      continue
    if preApexDisk_dp.contain(candidateApexDisk, resNeg6) or postApexDisk_dq.contain(candidateApexDisk, resNeg6):
      containedDisksInPreOrPostApexDisk.add(candidateApexDisk)
      inc(i)
      continue
    candidateTriangleApexNDiskPairs[j] = candidateTriangleApexNDiskPairs[i]
    inc(i)
    inc(j)
  candidateTriangleApexNDiskPairs.setLen(candidateTriangleApexNDiskPairs.len - (i - j))

proc find_sliver_triangle_configuration_and_apex_disks*(qhd: QuickhullDisk; 
    disks_D: var ListCircularDisk2D;
    nonNegativeSupportingTangentLineSegmentFrom_dp_to_dq: Line2D;
    preApexDisk_dp: CircularDisk2D; postApexDisk_dq: CircularDisk2D;
    triangleApexNDiskPairsOfExpandedNonPositive: var seq[(Point2D, CircularDisk2D)];
    triangleApexNDiskPairsOfOnPositiveExcept_dp_N_dq: var seq[(Point2D, CircularDisk2D)]): SliverConfiguration =
  var sliverConfiguration: SliverConfiguration
  var triangleApexNDiskPairs: seq[(Point2D, CircularDisk2D)]
  find_the_highest_triangle_apex_and_the_apex_disk_wrt_this_oriented_line_segment(qhd,
      disks_D, nonNegativeSupportingTangentLineSegmentFrom_dp_to_dq, triangleApexNDiskPairs, true, preApexDisk_dp, postApexDisk_dq)
  var candidateApexDisk: CircularDisk2D = triangleApexNDiskPairs[0][1]
  var heightOfTriangleApex: float = candidateApexDisk.get_radius() -
      nonNegativeSupportingTangentLineSegmentFrom_dp_to_dq.signed_distance(candidateApexDisk.get_center_pt())
  if ZERO(heightOfTriangleApex):
    var containedDisksInPreOrPostApexDisk: ListCircularDisk2D
    remove_pre_and_post_apex_disks_and_contained_disks_in_one_of_the_two_from_candidate_apex_disks_if_exist(qhd,
        preApexDisk_dp, postApexDisk_dq, triangleApexNDiskPairs,containedDisksInPreOrPostApexDisk)
    if containedDisksInPreOrPostApexDisk.len > 0:
      remove_contained_disks_from_input_disks(qhd, containedDisksInPreOrPostApexDisk, disks_D)
    if triangleApexNDiskPairs.len == 0:
      sliverConfiguration = SLIVER_CASE_A
    else: # (triangleApexNDiskPairs.size() >= 1)
      sliverConfiguration = SLIVER_CASE_B
      triangleApexNDiskPairsOfOnPositiveExcept_dp_N_dq = triangleApexNDiskPairs
  else:
    if triangleApexNDiskPairs.len == 1:
      sliverConfiguration = SLIVER_CASE_C1
      triangleApexNDiskPairsOfExpandedNonPositive.add(triangleApexNDiskPairs[0])
    else: #  (triangleApexNDiskPairs.size() >= 2)
      sliverConfiguration = SLIVER_CASE_C2
      triangleApexNDiskPairsOfExpandedNonPositive = triangleApexNDiskPairs
  return sliverConfiguration

proc regularize_sliver_triangle_and_repivot_disks(qhd: QuickhullDisk; 
    inputForFindHull: InputForFindHull; disks_D_FrontEdge: var ListCircularDisk2D;
    disks_D_BackEdge: var ListCircularDisk2D; pivotDisk_dx: var CircularDisk2D; pivotPoint_x: var Point2D) =
  var disks_D: ListCircularDisk2D = inputForFindHull.disks
  var preApexDisk_dp: CircularDisk2D = inputForFindHull.preApexDisk
  var postApexDisk_dq: CircularDisk2D = inputForFindHull.postApexDisk
  var hullPoint_p: Point2D = inputForFindHull.hullPointOfPreApexDisk
  var hullPoint_q: Point2D = inputForFindHull.hullPointOfPostApexDisk
  var nonNegativeSupportingTangentLineSegmentFrom_dp_to_dq: Line2D = compute_a_CCW_oriented_tangent_line_from_disk_d1_to_d2(qhd,
      preApexDisk_dp, postApexDisk_dq)
  var triangleApexNDiskPairsOfExpandedNonPositive: seq[(Point2D, CircularDisk2D)]
  var triangleApexNDiskPairsOfOnPositiveExcept_dp_N_dq: seq[(Point2D, CircularDisk2D)]
  var sliverConfiguration: SliverConfiguration = find_sliver_triangle_configuration_and_apex_disks(qhd, disks_D, nonNegativeSupportingTangentLineSegmentFrom_dp_to_dq,
      preApexDisk_dp, postApexDisk_dq, triangleApexNDiskPairsOfExpandedNonPositive, triangleApexNDiskPairsOfOnPositiveExcept_dp_N_dq)
  case sliverConfiguration
  of SLIVER_CASE_A:
    var triangleApexNDiskPairsOf_dp_N_dq: seq[(Point2D, CircularDisk2D)]
    triangleApexNDiskPairsOf_dp_N_dq.add((nonNegativeSupportingTangentLineSegmentFrom_dp_to_dq.get_start_point(), preApexDisk_dp))
    triangleApexNDiskPairsOf_dp_N_dq.add((nonNegativeSupportingTangentLineSegmentFrom_dp_to_dq.get_end_point(),postApexDisk_dq))
    pick_one_as_triangle_apex_and_apex_disk(qhd, triangleApexNDiskPairsOf_dp_N_dq, pivotPoint_x, pivotDisk_dx, qhd.m_CotangentDisksDividingOption)
  of SLIVER_CASE_B:
    if triangleApexNDiskPairsOfOnPositiveExcept_dp_N_dq.len == 1:
      pivotPoint_x = triangleApexNDiskPairsOfOnPositiveExcept_dp_N_dq[0][0]
      pivotDisk_dx = triangleApexNDiskPairsOfOnPositiveExcept_dp_N_dq[0][1]
    else: #  (triangleApexNDiskPairsOfOnPositiveExcept_dp_N_dq.size() >= 2)
      pick_one_as_triangle_apex_and_apex_disk_among_disks_with_identical_height_and_remove_disks_contained_in_others_from_input_disks_if_exist(qhd,
          triangleApexNDiskPairsOfOnPositiveExcept_dp_N_dq, pivotPoint_x,pivotDisk_dx, disks_D, qhd.m_CotangentDisksDividingOption)
  of SLIVER_CASE_C1:
    pivotPoint_x = triangleApexNDiskPairsOfExpandedNonPositive[0][0]
    pivotDisk_dx = triangleApexNDiskPairsOfExpandedNonPositive[0][1]
  of SLIVER_CASE_C2:
    pick_one_as_triangle_apex_and_apex_disk_among_disks_with_identical_height_and_remove_disks_contained_in_others_from_input_disks_if_exist(qhd,
        triangleApexNDiskPairsOfExpandedNonPositive, pivotPoint_x, pivotDisk_dx,disks_D, qhd.m_CotangentDisksDividingOption)
  var orientedLineSegment_FrontEdge_px: Line2D = initLine2D(hullPoint_p, pivotPoint_x)
  find_expanded_non_positive_disks_wrt_oriented_line_segment(qhd, disks_D, orientedLineSegment_FrontEdge_px, preApexDisk_dp, pivotDisk_dx,disks_D_FrontEdge, true)
  var orientedLineSegment_BackEdge_xq: Line2D = initLine2D(pivotPoint_x, hullPoint_q)
  find_expanded_non_positive_disks_wrt_oriented_line_segment(qhd, disks_D,
      orientedLineSegment_BackEdge_xq, pivotDisk_dx, postApexDisk_dq, disks_D_BackEdge, true)

proc find_hull_disks_by_iteration*(qhd: QuickhullDisk; 
    stackForFindingHull: var StackInputForFindHull; hullDisks: var ListCircularDisk2D) =
  while stackForFindingHull.len > 0:
    var inputOfCurrentStep: InputForFindHull = stackForFindingHull[^1]
    stackForFindingHull.setLen(stackForFindingHull.high)
    var disks_D: ListCircularDisk2D = inputOfCurrentStep.disks
    var preApexDisk_dp: CircularDisk2D = inputOfCurrentStep.preApexDisk
    var postApexDisk_dq: CircularDisk2D = inputOfCurrentStep.postApexDisk
    var hullPoint_p: Point2D = inputOfCurrentStep.hullPointOfPreApexDisk
    var hullPoint_q: Point2D = inputOfCurrentStep.hullPointOfPostApexDisk
    var numOfDisksInD: int = disks_D.len
    if numOfDisksInD == 1:
      hullDisks.add(preApexDisk_dp) # Output
    elif numOfDisksInD == 2 and (preApexDisk_dp != postApexDisk_dq):
      hullDisks.add(preApexDisk_dp) # Output
      hullDisks.add(postApexDisk_dq) # Output
    else:
      # *. Make two sets for next interation.
      var
        disks_D_FrontEdge: ListCircularDisk2D
        disks_D_BackEdge: ListCircularDisk2D
      var triangleApex_x: Point2D
      var apexDisk_dx: CircularDisk2D# = _NULL
      var orientedLineSegment_BaseEdge_pq: Line2D = initLine2D(hullPoint_p, hullPoint_q)
      var candidateApexNDiskPairs: seq[(Point2D, CircularDisk2D)]
      find_the_highest_triangle_apex_and_the_apex_disk_wrt_this_oriented_line_segment(qhd,
          disks_D, orientedLineSegment_BaseEdge_pq, candidateApexNDiskPairs)
      if candidateApexNDiskPairs.len == 1:
        triangleApex_x = candidateApexNDiskPairs[0][0]
        apexDisk_dx = candidateApexNDiskPairs[0][1]
      else: # (apexPointNItsDisks.size() > 1)
        pick_one_as_triangle_apex_and_apex_disk_among_disks_with_identical_height_and_remove_disks_contained_in_others_from_input_disks_if_exist(qhd,
            candidateApexNDiskPairs, triangleApex_x, apexDisk_dx, disks_D, qhd.m_CotangentDisksDividingOption)
      var orientedLineSegment_FrontEdge_px: Line2D = initLine2D(hullPoint_p, triangleApex_x)
      var orientedLineSegment_BackEdge_xq: Line2D = initLine2D(triangleApex_x, hullPoint_q)
      find_expanded_non_positive_disks_wrt_oriented_line_segment(qhd, disks_D,
          orientedLineSegment_FrontEdge_px, preApexDisk_dp, apexDisk_dx, disks_D_FrontEdge, true)
      find_expanded_non_positive_disks_wrt_oriented_line_segment(qhd, disks_D,
          orientedLineSegment_BackEdge_xq, apexDisk_dx, postApexDisk_dq, disks_D_BackEdge, true)
      # * Intriguing case
      if triangle_filter_is_sliver(disks_D.len, disks_D_FrontEdge.len, disks_D_BackEdge.len, preApexDisk_dp, postApexDisk_dq):
        disks_D_FrontEdge.setLen(0)
        disks_D_BackEdge.setLen(0)
        regularize_sliver_triangle_and_repivot_disks(qhd, inputOfCurrentStep, disks_D_FrontEdge, disks_D_BackEdge, apexDisk_dx, triangleApex_x)
      prepare_and_insert_input_data_for_finding_hull_to_stack(disks_D_BackEdge, apexDisk_dx, postApexDisk_dq, triangleApex_x, hullPoint_q,
          stackForFindingHull)
      prepare_and_insert_input_data_for_finding_hull_to_stack(disks_D_FrontEdge,
          preApexDisk_dp, apexDisk_dx, hullPoint_p, triangleApex_x, stackForFindingHull)
  # *. Current 'hullDisks' has same disks consecutively, so delete the repeated disks and make one.
  hulldisks = hulldisks.deduplicate(isSorted = true)

proc divide_input_disks_into_two_initial_subsets*(qhd: QuickhullDisk; 
    inputDisks_D: ListCircularDisk2D; orientedLineSegment_pq: Line2D;
    initialDisks_D_right: var ListCircularDisk2D;
    initialDisks_D_left: var ListCircularDisk2D) =
  for it_Disk in inputDisks_D:
    var currDisk: CircularDisk2D = it_Disk
    if this_disk_is_a_member_of_expanded_non_positive_set_wrt_line(qhd, currDisk, orientedLineSegment_pq, true): # last arg true to fix for zero radius
      initialDisks_D_right.add(currDisk)
    if this_disk_is_a_member_of_expanded_non_negative_set_wrt_line(qhd, currDisk, orientedLineSegment_pq, true):
      initialDisks_D_left.add(currDisk)

proc find_hull_disks_by_QHULLDISK*(qhd: QuickhullDisk;  inputDisks_D: ListCircularDisk2D; hullDisks: var ListCircularDisk2D) =
  case inputDisks_D.len
  of 0:
    discard
  of 1:
    hullDisks.add(inputDisks_D[0])
  else:                       #  inputDisks_D.size() >= 2
    # *. Initial step
    var
      highLeftExtremePoint_p: Point2D
      lowRightExtremePoint_q: Point2D
      diskHavingHighLeftPoint_p: CircularDisk2D
      diskHavingLowRightPoint_q: CircularDisk2D
    find_high_left_N_low_right_extreme_points_and_their_disks_to_divide_disk_set(qhd,
        inputDisks_D, highLeftExtremePoint_p, lowRightExtremePoint_q,
        diskHavingHighLeftPoint_p, diskHavingLowRightPoint_q)
    var
      initialExpandedNonPositiveDisks_D_right: ListCircularDisk2D
      initialExpandedNonNegativeDisks_D_left: ListCircularDisk2D
      orientedLineSegment_BaseLine_pq: Line2D = initLine2D(highLeftExtremePoint_p, lowRightExtremePoint_q)
    divide_input_disks_into_two_initial_subsets(qhd, inputDisks_D,
        orientedLineSegment_BaseLine_pq, initialExpandedNonPositiveDisks_D_right,
        initialExpandedNonNegativeDisks_D_left)
    var stackForFindingHull: StackInputForFindHull
    prepare_and_insert_input_data_for_finding_hull_to_stack(
        initialExpandedNonNegativeDisks_D_left, diskHavingLowRightPoint_q,
        diskHavingHighLeftPoint_p, lowRightExtremePoint_q, highLeftExtremePoint_p,
        stackForFindingHull)
    prepare_and_insert_input_data_for_finding_hull_to_stack(
        initialExpandedNonPositiveDisks_D_right, diskHavingHighLeftPoint_p,
        diskHavingLowRightPoint_q, highLeftExtremePoint_p, lowRightExtremePoint_q,
        stackForFindingHull)
    # *. Main step
    find_hull_disks_by_iteration(qhd, stackForFindingHull, hullDisks)

#  \brief Find hull disks by QuickhullDisk; This function calls core function 'find_hull_disks_by_QHULLDISK'.
#
#  \param[in]  inputDisks                     The input disks of convex hull problem.
#  \param[in]  cotangentDisksDividingOption   The dividing option for cotangent disks.
#
#  \param[out] hullDisks     The output hull disks stored in CCW direction starting from a disk having high left extreme point.
#                            The same disk may occur many times but not consecutively;
#                             ex. Possible result:(d_1, d_2, d_1, d_3, d_4, d_1, ...),
#                                 impossible result: (d_1, d_2, d_2, d_2, d_4, d_1, ...) because d_2 occurs consecutively.
#                            The first disk must occurs at the last; ex. (d_1, ..., d_1).
#
proc find_hull_disks*(qhd: var QuickhullDisk; inputDisks: ListCircularDisk2D; hullDisks: var ListCircularDisk2D;
  cotangentDisksDividingOption: CotangentDisksDividingOption = RANDOM_DIVIDE) =
  # Erase input disks, output extreme disks, and time statistics used before.
  clear_convexhull(qhd)
  # Set memory of the computation times of each step(#steps: QH_TIME_SIZE)
  # Set cotangent disks' dividing option between RANDOM_DIVIDE and BISECTION_DIVIDE.
  qhd.m_CotangentDisksDividingOption = cotangentDisksDividingOption
  # *. For random selection of disks tangent to a common line, we initialize the seed number.
  if qhd.m_CotangentDisksDividingOption == RANDOM_DIVIDE:
    let seed = times.getTime().nanoSecond
    echo "Init random with seed: ", seed
    random.randomize(seed)
  # Set the input disks
  qhd.m_InputDisks = inputDisks
  # *. Extract pointers of input disks, which is a input format for 'find_hull_disks_by_QHULLDISK'.
  var inputDisksForQHullDisk: ListCircularDisk2D
  get_pointers_of_disks(qhd, qhd.m_InputDisks, inputDisksForQHullDisk)
  # *. Run QuickhullDisk [core function].
  find_hull_disks_by_QHULLDISK(qhd, inputDisksForQHullDisk, qhd.m_HullDisks)
  # *. Rearrange for output format
  rearrange_hull_disks_for_output_format_by_begining_with_disk_having_high_left_extreme_point(qhd, qhd.m_HullDisks)
  # *. Copy extreme disks from disk pointers, which is a output format for 'find_hull_disks'.
  copy_disks_from_its_pointers(qhd, qhd.m_HullDisks, hullDisks)

# 392 lines

when isMainModule:

  var qhd: QuickhullDisk
  var d1: CircularDisk2D = initCircularDisk2D(1.0, 2.0, 1e-7, 0)
  var s: ListCircularDisk2D
  var res: ListCircularDisk2D
  s.add(d1)
  s.add(initCircularDisk2D(-1.0, 2.0, 1e-7, 1))
  s.add(initCircularDisk2D(-1.2, 2.2, 1e-7, 5))
  s.add(initCircularDisk2D(3.0, 5.0, 1e-7, 2))
  echo s
  find_hull_disks(qhd, s, res)
  echo "==="
  for el in res:
    echo el
