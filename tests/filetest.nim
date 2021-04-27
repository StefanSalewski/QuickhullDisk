import quickhulldisk/[quickhulldisk, convexhullofdisks, circulardisk2d]

import parsecsv, strutils
from os import paramStr
from streams import newFileStream

proc main =
  var qhd: QuickhullDisk
  var s: ListCircularDisk2D
  var res: ListCircularDisk2D
  var stream1 = newFileStream(paramStr(1), fmRead)
  if stream1 == nil:
    quit("cannot open the file" & paramStr(1))
  var parser: CsvParser
  open(parser, stream1, paramStr(1), separator = '\t')
  parser.readHeaderRow()
  while readRow(parser):
    let id = parseInt(parser.row[0])
    let x = parseFloat(parser.row[1])
    let y = parseFloat(parser.row[2])
    let r = parseFloat(parser.row[3])
    s.add(initCircularDisk2D(x, y, r, id))
  close(parser)
  find_hull_disks(qhd, s, res)
  echo "Number of hull disk: ", res.len
  stream1 = newFileStream(paramStr(2), fmRead)
  if stream1 == nil:
    quit("cannot open the file" & paramStr(2))
  open(parser, stream1, paramStr(2), separator = '\t')
  parser.readHeaderRow()
  var line: int
  while readRow(parser):
    if not res[line].get_id == parseInt(parser.row[0]):
      echo "Computed result differ from result file at line ", line
      assert false
  echo "OK"

main()


