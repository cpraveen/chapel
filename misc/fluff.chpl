use IO;
use StencilDist;
use Math;

config const nx = 5,   // number of cells in x direction
             ny = 5;   // number of cells in y direction

// main function
proc main()
{
  writeln("Grid size is ",nx," x ",ny);

  const D  = {1..nx, 1..ny};
  param rank = D.rank;

  // problem space
  var halo: rank*int = (2,2);
  const PSpace = stencilDist.createDomain(D, fluff=halo, periodic=true);

  var u : [PSpace] real;

  // Set initial condition
  forall (i,j) in D
  {
    u[i,j] = i + ny * (j-1);
  }
  u.updateFluff();

  for j in -1..ny+2
  {
     for i in -1..nx+2
     {
        writef("%4.0r", u[i,j]);
     }
     writeln();
  }


}
