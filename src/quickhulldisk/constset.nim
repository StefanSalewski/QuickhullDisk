from fenv import maximumPositiveValue

const
  DBL_MAX* = float.maximumPositiveValue

  # DegreeToRad* = 57.2957795130823
  RadToDegree = 57.2957795130823 # By Joonghyun on April 18, 2019 rename from _DegreeToRad 

# TOLERANCES
const
  resNeg20* = 1.0e-20 
  resNeg19* = 1.0e-19 
  resNeg18* = 1.0e-18 
  resNeg17* = 1.0e-17 
  resNeg16* = 1.0e-16 
  resNeg15* = 1.0e-15 
  resNeg14* = 1.0e-14 
  resNeg13* = 1.0e-13 
  resNeg12* = 1.0e-12 
  resNeg11* = 1.0e-11 
  resNeg10* = 1.0e-10 
  resNeg9* = 1.0e-9 
  resNeg8* = 1.0e-8 
  resNeg7* = 1.0e-7 
  resNeg6* = 1.0e-6 
  resNeg5* = 1.0e-5 
  resNeg4* = 1.0e-4 
  resNeg3* = 1.0e-3 
  resNeg2* = 1.0e-2 
  resNeg1* = 1.0e-1 

  SYSTEM_RES* = 1.0e-15 #system resolution 
  MATH_RES*   = 1.0e-6 # mathematical resolution 
  MAX_REAL*   = 1.0e38 #the largest possible real value 

