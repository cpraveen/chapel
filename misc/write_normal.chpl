extern
{
#include <float.h>
}

proc normal(x) do return if abs(x) < DBL_MIN then 0.0 else x;

proc main()
{
   var x = 1.0e-200;
   var y = 1.0e-311;

   writeln("x, normal(x) = ", x, " ", normal(x));
   writeln("y, normal(y) = ", y, " ", normal(y));
}
