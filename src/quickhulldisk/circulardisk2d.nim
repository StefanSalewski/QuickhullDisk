import circle2d, point2d

type
  CircularDisk2D* = object of Circle2D
    m_ID: int

proc initCircularDisk2D*(): CircularDisk2D =
  CircularDisk2D(m_ID: -1)

# https://forum.nim-lang.org/t/7708
proc initCircularDisk2D*(center: Point2D; r: float; id: int; T:typedesc[CircularDisk2D] = CircularDisk2D): T =
  result = initCircle2D(center, r, T)
  result.m_ID = id

proc initCircularDisk2D*(pointX, pointY, r: float; id: int; T:typedesc[CircularDisk2D] = CircularDisk2D): T =
  result = initCircle2D(pointX, pointY, r, T)
  result.m_ID = id

let
  CircularDisk2DNil* = initCircularDisk2D(system.INF, system.INF, system.INF, int.low)

# CAUTION: maybe first par should be of type Circle2d ? Typo in C++ code?
proc initCircularDisk2D*(circle: CircularDisk2D; id: int; T:typedesc[CircularDisk2D] = CircularDisk2D): T =
  result = initCircle2D(circle.get_circle(), T) # get_circle() not really neded
  result.m_ID = id

proc initCircularDisk2D*(disk: CircularDisk2D; T:typedesc[CircularDisk2D] = CircularDisk2D): T =
  result = initCircle2D(disk, T)
  result.m_ID = disk.get_ID()

proc get_ID*(cd: CircularDisk2D): int = cd.m_ID

proc set_ID*(cd: var CircularDisk2D; id: int) = cd.m_ID = id

proc `<`*(cd, disk: CircularDisk2D): bool =
  return cd.get_radius() < disk.get_radius()

proc `>`*(cd, disk: CircularDisk2D): bool =
  return cd.get_radius() > disk.get_radius()

proc `==`*(cd, disk: CircularDisk2D): bool =
  return disk.m_ID == cd.m_ID and Circle2D(disk) == Circle2D(cd)
