module Problem
{

use Constants;

const periodic = false;
const xmin = 0.0, xmax = 1.0;
const ymin = 0.0, ymax = 1.0;
const Uwall = 1.0;

const west  = BC.noslip;
const east  = BC.noslip;
const south = BC.noslip;
const north = BC.noslip;

proc initial(ref flag, ref u, ref v, ref p)
{
}

proc boundary(ref u, ref v, ref p)
{
   const ny = p.domain.dim(1);

   // Top lid velocity
   u[..,ny+1] = 2.0*Uwall - u[..,ny];
}

}

