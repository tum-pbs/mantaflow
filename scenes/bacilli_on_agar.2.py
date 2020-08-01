#
# Bacilli on Agar
# 

from manta import *

# solver params
dim = 3
res = 64
# res = 128
gs = vec3(res,res,res)

if (dim==2):
	gs.z=1

solver = Solver(name='main', gridSize = gs, dim=dim)
solver.timestep = 0.8
minParticles = pow(2,dim)

# save particles for separate surface generation pass?
saveParts = False

# size of particles 
radiusFactor = 1.0

# prepare grids and particles
flags    = solver.create(FlagGrid)
phi 	 = solver.create(LevelsetGrid)
vel      = solver.create(MACGrid)
velOld   = solver.create(MACGrid)
pressure = solver.create(RealGrid)
tmpVec3  = solver.create(VecGrid)
force    = solver.create(VecGrid)
pp       = solver.create(BasicParticleSystem)  
pVel     = pp.create(PdataVec3) 
mesh     = solver.create(Mesh)

# acceleration data for particle nbs
pindex = solver.create(ParticleIndexSystem) 
gpi    = solver.create(IntGrid)
bWidth=1
flags.initDomain(boundaryWidth=bWidth)
fluidVel = 0
fluidSetVel = 0

# The scale is 1 = 10 micrometer (1 x 10E-7 meter)
# Agar
agar = Box(parent=solver, p0=gs*vec3(0,0,0), p1=gs*vec3(1.0,0.1,1.0))
phi = agar.computeLevelset()

# Bacillus
center = vec3(0.5,0.1125,0.5)
radius = 0.0125
zpoint = vec3(0, 0, 0.1)
# The mantaflow cylinder has a capsule-like geometry
bacillus = Cylinder(parent=solver, center=gs*center, radius=res*radius, z=gs*zpoint)
phiBacillus = bacillus.computeLevelset()
setObstacleFlags(flags=flags, phiObs=phiBacillus) 

# Join the agar and bacillus levelsets
phi.join(phiBacillus)

flags.updateFromLevelset(phi)
sampleLevelsetWithParticles(phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.05)

# save reference any grid, to automatically determine grid size
if saveParts:
	pressure.save('ref_flipParts_0000.uni');

if 1 and (GUI):
	gui = Gui()
	gui.show()
	gui.pause()
   

#main loop
for t in range(2500):
	# mantaMsg('\nFrame %i, simulation time %f' % (solver.frame, solver.timeTotal))
	
	# Algorithm 1 in Gao, 2018
	# 1: Advect velocities of particles 
	pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False )
	
	# make sure we have velocities throughout liquid region
	mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 ) 
	extrapolateMACFromWeight(vel=vel , distance=2, weight=tmpVec3)  # note, tmpVec3 could be free'd now...
	markFluidCells(parts=pp, flags=flags)

	# create approximate surface level set, resample particles
	gridParticleIndex(parts=pp , flags=flags, indexSys=pindex, index=gpi)
	unionParticleLevelset(pp, pindex, flags, gpi, phi , radiusFactor) 
	resetOutflow(flags=flags,parts=pp,index=gpi,indexSys=pindex) 
	# extend levelset somewhat, needed by particle resampling in adjustNumber
	extrapolateLsSimple(phi=phi, distance=4, inside=True); 

	setWallBcs(flags=flags, vel=vel)	
	solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)

	# set source grids for resampling, used in adjustNumber!
	pVel.setSource( vel, isMAC=True )
	adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor ) 

	# make sure we have proper velocities
	extrapolateMACSimple( flags=flags, vel=vel )

	# 2: Enforce external forces (gravity)
	addGravity(flags=flags, vel=vel, gravity=(0,-0.001,0))
	
	if t == 3:
		force.setConst(vec3(0.0, 0.0, 0.5))
		addForceField(flags=flags, vel=vel, force=force)

	# 3: Verify fluid and solid particle flags via Fsolid, Ffluid
	# 4: Map all the particles to grid ug ← up, xg ← xp
	# 5: Compute level set Φ and velocity on grid ug
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
	flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )

	# 16: Update particles’ flags

	if (dim==3):
		phi.createMesh(mesh)
		mesh.fromShape(bacillus, append=True)

	#solver.printMemInfo()
	solver.step()

	# generate data for flip03_gen.py surface generation scene
	# if saveParts:
	# 	pp.save( 'flipParts_%04d.uni' % t );



