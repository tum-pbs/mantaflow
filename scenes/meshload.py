#
# simple example to rasterize an obj file and create an SDF for an obstacle
#
from manta import *
import os

# mesh to load
# meshfile = '../resources/simpletorus.obj'
meshfile = './resources/simpletorus.obj'
mantaMsg("Loading %s. Note: relative path by default, assumes this scene is called from the 'scenes' directory.")

# resolution for level set / output mesh
res = 50 
gs = vec3(res,res,res)
s = Solver(name='main', gridSize = vec3(res,res,res) , dim=3)

flags    = s.create(FlagGrid)
density  = s.create(RealGrid)
vel      = s.create(MACGrid)
pressure = s.create(RealGrid)
phiObs   = s.create(LevelsetGrid)
mesh     = s.create(Mesh)


# === load mesh, and turn into SDF ===
mesh.load( meshfile )
mesh.scale( vec3(res/3.0) )
mesh.offset( gs* (Vec3(0.5) + Vec3(0.1,0.05,0)) ) # center + slight offset
mesh.computeLevelset(phiObs, 2.)

flags.initDomain()
setObstacleFlags(flags=flags, phiObs=phiObs) #, fractions=fractions)
flags.fillGrid()


# run simple smoke sim 

if 1 and (GUI):
	gui = Gui()
	gui.show()

source = s.create(Cylinder, center=gs*vec3(0.35,0.2,0.5), radius=res*0.15, z=gs*vec3(0, 0.05, 0))

for t in range(250):
	mantaMsg('\nFrame %i' % (s.frame))
	source.applyToGrid(grid=density, value=1.)
	#densityInflow(flags=flags, density=density, shape=source, scale=1, sigma=0.5)

	advectSemiLagrange(flags=flags, vel=vel, grid=density, order=2)    
	advectSemiLagrange(flags=flags, vel=vel, grid=vel    , order=2, strength=1.0)
	
	setWallBcs(flags=flags, vel=vel)    
	addBuoyancy(density=density, vel=vel, gravity=vec3(0,-1e-3,0), flags=flags)
	
	solvePressure( flags=flags, vel=vel, pressure=pressure )
	s.step()

