use StencilDist;
use Problem;
use Constants;

config const nx  = 100,
             ny  = 100,
             Tf  = 1.0,
             omg = 1.5,
             nu  = 1.0,
             itmax = 1000,
             tol = 1.0e-3;

const dx   = (xmax - xmin)/nx;
const dy   = (xmax - xmin)/ny;
const idx  = 1.0/dx;
const idy  = 1.0/dy;
const idx2 = 1.0/dx**2;
const idy2 = 1.0/dy**2;
const beta = 0.5/(idx2 + idy2);

const D = {1..nx,   1..ny};
const Dx= {1..nx+1, 1..ny};
const Dy= {1..nx,   1..ny+1};

const cell  = stencilDist.createDomain(D,  fluff=(1,1), periodic=periodic);
const xface = stencilDist.createDomain(Dx, fluff=(1,1), periodic=periodic);
const yface = stencilDist.createDomain(Dy, fluff=(1,1), periodic=periodic);

const uinner= {2..nx, 1..ny};
const vinner= {1..nx, 2..ny};

//-----------------------------------------------------------------------------
proc initial_condition(ref flag, ref u, ref v, ref p)
{
   // set all real cells to fluid
   forall ij in cell do
      flag[ij] = C_F;

   // set external ghost cells to obstacle
   forall i in 0..nx+1
   {
      flag[i,0   ] = C_B;
      flag[i,ny+1] = C_B;
   }
   forall j in 0..ny+1
   {
      flag[0,   j] = C_B;
      flag[nx+1,j] = C_B;
   }

   u = 0.0;
   v = 0.0;
   p = 0.0;

   // Set problem specfic flag and ic
   initial(flag, u, v, p);
}

//-----------------------------------------------------------------------------
proc boundary_condition(ref u, ref v, ref p)
{
   select west
   {
      when BC.slip   
      { 
         u[1,..] = 0.0; 
         v[0,..] = v[1,..]; 
      }
      when BC.noslip 
      { 
         u[1,..] =  0.0; 
         v[0,..] = -v[1,..]; 
      }
      when BC.outflow
      { 
         u[1,..] = u[2,..]; 
         v[0,..] = v[1,..]; 
      }
   }

   select east
   {
      when BC.slip
      {
         u[nx+1,..] = 0.0;
         v[nx+1,..] = v[nx,..];
      }
      when BC.noslip
      {
         u[nx+1,..] =  0.0;
         v[nx+1,..] = -v[nx,..];
      }
      when BC.outflow
      {
         u[nx+1,..] = u[nx,..];
         v[nx+1,..] = v[nx,..];
      }
   }

   // problem specific bc
   boundary(u, v, p);

   // this does periodic also
   u.updateFluff();
   v.updateFluff();
   p.updateFluff();
}

//-----------------------------------------------------------------------------
proc time_step(u, v)
{
   const dtx = min reduce [ij in xface] dx/(abs(u[ij]) + 1.0e-12);
   const dty = min reduce [ij in yface] dy/(abs(v[ij]) + 1.0e-12);

   return min(dtx,dty);
}

//-----------------------------------------------------------------------------
proc predict(flag, u, v, dt, ref U, ref V)
{
   forall (i,j) in uinner do
   if (flag[i,j]   & C_F != 0) && 
      (flag[i,j]   < C_E     ) && 
      (flag[i-1,j] & C_F != 0) &&
      (flag[i-1,j] < C_E     )
   {
      const uw = 0.5*(u[i,j  ] + u[i-1,j  ]);
      const ue = 0.5*(u[i,j  ] + u[i+1,j  ]);
      const us = 0.5*(u[i,j  ] + u[i  ,j-1]);
      const un = 0.5*(u[i,j  ] + u[i  ,j+1]);
      const vs = 0.5*(v[i,j  ] + v[i-1,j  ]);
      const vn = 0.5*(v[i,j+1] + v[i-1,j+1]);
      const C  = idx*(ue**2 - uw**2) + idy*(un * vn - us * vs)/dy;
      const uxx= idx2*(u[i-1,j] - 2.0*u[i,j] + u[i+1,j]);
      const uyy= idy2*(u[i,j-1] - 2.0*u[i,j] + u[i,j+1]);
      U[i,j] = u[i,j] + dt * (nu*(uxx + uyy) - C);
   }
   else
   {
      U[i,j] = u[i,j];
   }

   forall (i,j) in vinner do
   if (flag[i,j]   & C_F != 0) && 
      (flag[i,j]   < C_E     ) && 
      (flag[i,j-1] & C_F != 0) &&
      (flag[i,j-1] < C_E     )
   {
      const uw = 0.5*(u[i  ,j] + u[i  ,j-1]);
      const ue = 0.5*(u[i+1,j] + u[i+1,j-1]);
      const vw = 0.5*(v[i  ,j] + v[i-1,j]);
      const ve = 0.5*(v[i  ,j] + v[i+1,j]);
      const vs = 0.5*(v[i  ,j] + v[i  ,j-1]);
      const vn = 0.5*(v[i  ,j] + v[i  ,j+1]);
      const C  = idx*(ue * ve - uw * vw) + idy*(vn**2 - vs**2)/dy;
      const vxx= idx2*(v[i-1,j] - 2.0*v[i,j] + v[i+1,j]);
      const vyy= idy2*(v[i,j-1] - 2.0*v[i,j] + v[i,j+1]);
      V[i,j] = v[i,j] + dt * (nu*(vxx + vyy) - C);
   }
   else
   {
      V[i,j] = v[i,j];
   }

   forall i in 1..nx+1
   {
      U[i,0   ] = u[i,0   ]; // bottom
      U[i,ny+1] = u[i,ny+1]; // top
   }

   forall j in 1..ny+1
   {
      V[0   ,j] = v[0   ,j]; // left
      V[nx+1,j] = v[nx+1,j]; // right
   }
}

//-----------------------------------------------------------------------------
proc pressure(flag, U, V, dt, ref p)
{
   // rhs of poisson equation
   var rhs : [cell] real;
   forall (i,j) in cell do
   if (flag[i,j] & C_F != 0) && (flag[i,j] < 0x0100) then
      rhs[i,j] = (idx*(U[i+1,j] - U[i,j]) + idy*(V[i,j+1] - V[i,j]))/dt;

   // pressure scale to decide convergence
   var p0 = 0.0, ifull = 0;
   forall ij in cell with (+ reduce ifull, + reduce p0) do
   if flag[ij] & C_F
   {
      p0    += p[ij]**2;
      ifull += 1;
   }

   p0 = sqrt(p0/ifull);
   if p0 < 0.0001 then
      p0 = 1.0;

   var it = 0, res = 1.0e20;
   while it < itmax && res > tol
   {
      // Neumann bc on pressure
      forall i in 1..nx
      {
         p[i,0]    = p[i,1];  // bottom
         p[i,ny+1] = p[i,ny]; // top
      }

      forall j in 1..ny
      {
         p[0,j]    = p[1,j];  // left
         p[nx+1,j] = p[nx,j]; // right
      }

      // Boundaries in the interior for embedded case
      forall (i,j) in cell do
      if (flag[i,j] >= B_N) && (flag[i,j] <= B_S) then
      select flag[i,j]
      {
         when B_N  { p[i,j] = p[i  ,j+1]; }
         when B_O  { p[i,j] = p[i+1,j  ]; }
         when B_S  { p[i,j] = p[i  ,j-1]; }
         when B_W  { p[i,j] = p[i-1,j  ]; }
         when B_NO { p[i,j] = 0.5*(p[i,j+1] + p[i+1,j]); }
         when B_SO { p[i,j] = 0.5*(p[i,j-1] + p[i+1,j]); }
         when B_SW { p[i,j] = 0.5*(p[i,j-1] + p[i-1,j]); }
         when B_NW { p[i,j] = 0.5*(p[i,j+1] + p[i-1,j]); }
      }

      // Red-black Gauss-Seidel
      // Do Red
      forall i in cell.dim(0) do
      forall j in cell.dim(1) by 2 align i do
      if (flag[i,j] & C_F != 0) && (flag[i,j] < 0x0100)
      {
         const tmp = beta*(  idx2*(p[i+1,j] + p[i-1,j])
                           + idy2*(p[i,j+1] + p[i,j-1]) - rhs[i,j]);
         p[i,j] = (1.0-omg)*p[i,j] + omg*tmp;
      }

      // Do Black
      forall i in cell.dim(0) do
      forall j in cell.dim(1) by 2 align i+1 do
      if (flag[i,j] & C_F != 0) && (flag[i,j] < 0x0100)
      {
         const tmp = beta*(  idx2*(p[i+1,j] + p[i-1,j])
                           + idy2*(p[i,j+1] + p[i,j-1]) - rhs[i,j]);
         p[i,j] = (1.0-omg)*p[i,j] + omg*tmp;
      }

      // Compute residual
      res = 0.0;
      forall (i,j) in cell with (+ reduce res) do
      if (flag[i,j] & C_F != 0) && (flag[i,j] < 0x0100) then
      {
         const tmp =   idx2*(p[i-1,j] - 2.0*p[i,j] + p[i+1,j])
                     + idy2*(p[i,j-1] - 2.0*p[i,j] + p[i,j+1])
                     - rhs[i,j];
         res += tmp**2;
      }
      res = sqrt(res/ifull) / p0;
      it += 1;
   }
}

//-----------------------------------------------------------------------------
proc correct(flag, U, V, p, dt, ref u, ref v)
{
   forall (i,j) in uinner do    
   if (flag[i,j]   & C_F != 0) && 
      (flag[i,j]   < C_E     ) && 
      (flag[i-1,j] & C_F != 0) &&
      (flag[i-1,j] < C_E     ) // only if both adjacent cells are fluid cells
   {
      u[i,j] = U[i,j] - dt*idx*(p[i,j] - p[i-1,j]);
   }

   forall (i,j) in vinner do
   if (flag[i,j]   & C_F != 0) && 
      (flag[i,j]   < C_E     ) && 
      (flag[i,j-1] & C_F != 0) &&
      (flag[i,j-1] < C_E     )  // only if both adjacent cells are fluid cells
   {
      v[i,j] = V[i,j] - dt*idy*(p[i,j] - p[i,j-1]);
   }
}

//-----------------------------------------------------------------------------
proc main()
{
   var u, U : [xface] real;
   var v, V : [yface] real;
   var p    : [cell]  real;
   var flag : [cell]  int;

   initial_condition(flag, u, v, p);
   var dt = time_step(u, v);

   var t = 0.0, it = 0;
   while t < Tf
   {
      boundary_condition(u,v,p);
      predict(flag, u, v, dt, U, V);
      pressure(flag, U, V, dt, p);
      correct(flag, U, V, p, dt, u, v);
      t += dt; it += 1;
   }
}
