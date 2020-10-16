#
# Very simple flip without level set
# and without any particle resampling
# 
from manta import *

# solver params
dim = 3
particleNumber = 2
res = 64
gs = vec3(res,res,res)

if (dim==2):
	gs.z=1
	particleNumber = 3 # use more particles in 2d

s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep = 0.5

# prepare grids and particles
flags    = s.create(FlagGrid)
vel      = s.create(MACGrid)
velOld   = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)
pp       = s.create(BasicParticleSystem) 
# add velocity data to particles
pVel     = pp.create(PdataVec3)
mesh     = s.create(Mesh) 

# scene setup
flags.initDomain(boundaryWidth=0) 

# The scale is 1 = 10 micrometer (1 x 10E-7 meter)
# Agar
agar = Box(parent=s, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.1,1.0))
phi = agar.computeLevelset()

# Bacillus
center = vec3(0.5,0.2,0.5)
radius = 0.05
zpoint = vec3(0, 0, 0.1)
# The mantaflow cylinder has a capsule-like geometry
bacillus = Cylinder(parent=s, center=gs*center, radius=res*radius, z=gs*zpoint)
phiBacillus = bacillus.computeLevelset()

# Fluid Cylinder
centerTwo = vec3(0.3,0.2,0.5)
cylinder = Cylinder(parent=s, center=gs*centerTwo, radius=res*radius, z=gs*zpoint)
phiCylinder = cylinder.computeLevelset()

# phi.join(phiBacillus)

# Set the flags to solid for the bacillus 
setSolidFlags(flags=flags, phiSolid=phiBacillus)
# setSolidFlags(flags=flags, phiSolid=phi)

# Join the agar, bacillus, and fluid cylinder levelsets
phi.join(phiBacillus)
phi.join(phiCylinder)

# Set the flags to fluid for the agar and cylinder
flags.updateFromLevelset(phi)

# There's no resamplig here, so we need _LOTS_ of particles...
sampleFlagsWithParticles(flags=flags, parts=pp, discretization=particleNumber, randomness=0.2)

indices = []

for index in range(pp.pySize()): 
	inSolid = isParticleInSolid(index=index, particles=pp, flags=flags)

	if inSolid:
		indices.append(index)
		
# print(indices)

if (GUI):
	gui = Gui()
	gui.show()
	gui.pause()

# if (dim==3):
# 	phi.createMesh(mesh)
	
#main loop
for t in range(2500):
	# mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))

	# Algorithm 1 in Gao, 2018
	
	# 1: Advect velocities of particles 
	# 2: Enforce external forces (gravity)
	# 3: Verify fluid and solid particle flags via Fsolid, Ffluid
	# 4: Map all the particles to grid ug ← up, xg ← xp
	# 6: Project up ← ug,xp ← xg0 + ugΔt
	# 7: if particle ∈ Fsolid then
	# 8: 	Project shape matching constraint
	# 9: 	Compute target position g
	# 10: 	Update x∗+=α 􏰉g − x 􏰊
	# 11: 	Update u p∗ ← ( u p , ( x p∗ − x p0 ) / Δ t )
	# 12: else continue
	# 13: Correct boundary conditions
	# 14: Conduct the mechanism of penetration prevention
	# 15: Update the velocities and positions of all particles
	# 16: Update particles’ flags
	
	# FLIP
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False ) 
	mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
	extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 ) 
	markFluidCells(parts=pp, flags=flags)
	clearSolidFlags(flags = flags)

	for index in indices:
		markCellSolid(index=index, particles=pp, flags=flags)
	
	addGravity(flags=flags, vel=vel, gravity=(0,-0.002,0))

	# pressure solve
	setWallBcs(flags=flags, vel=vel)    
	solvePressure(flags=flags, vel=vel, pressure=pressure)
	setWallBcs(flags=flags, vel=vel)

	# we dont have any levelset, ie no extrapolation, so make sure the velocities are valid
	extrapolateMACSimple(flags=flags, vel=vel)
	
	# FLIP velocity update
	flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )

	s.step()