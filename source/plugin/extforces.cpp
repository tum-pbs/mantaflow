/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Set boundary conditions, gravity
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "grid.h"
#include "commonkernels.h"
#include "particle.h"

using namespace std;

namespace Manta { 

//! add constant force between fl/fl and fl/em cells
KERNEL(bnd=1) void KnApplyForceField(const FlagGrid& flags, MACGrid& vel, const Grid<Vec3>& force, const Grid<Real>* include, bool additive, bool isMAC) {
	bool curFluid = flags.isFluid(i,j,k);
	bool curEmpty = flags.isEmpty(i,j,k);
	if (!curFluid && !curEmpty) return;
	if (include && ((*include)(i,j,k) > 0.)) return;

	Real forceX = (isMAC) ? force(i,j,k).x : 0.5*(force(i-1,j,k).x + force(i,j,k).x);
	Real forceY = (isMAC) ? force(i,j,k).y : 0.5*(force(i,j-1,k).y + force(i,j,k).y);

	Real forceZ = 0.;
	if(vel.is3D())  
		forceZ =  (isMAC) ? force(i,j,k).z : 0.5*(force(i,j,k-1).z + force(i,j,k).z);

	if (flags.isFluid(i-1,j,k) || (curFluid && flags.isEmpty(i-1,j,k))) 
		vel(i,j,k).x = (additive) ? vel(i,j,k).x+forceX : forceX;
	if (flags.isFluid(i,j-1,k) || (curFluid && flags.isEmpty(i,j-1,k))) 
		vel(i,j,k).y = (additive) ? vel(i,j,k).y+forceY : forceY;
	if (vel.is3D() && (flags.isFluid(i,j,k-1) || (curFluid && flags.isEmpty(i,j,k-1))))
		vel(i,j,k).z = (additive) ? vel(i,j,k).z+forceZ : forceZ;
}

// ! add constant force between fl/fl and fl/em cells
// KERNEL(bnd=1) void KnApplyForce(const FlagGrid& flags, MACGrid& vel, Vec3 force, const Grid<Real>* exclude, bool additive) {
// 	bool curFluid = flags.isFluid(i,j,k);
// 	bool curEmpty = flags.isEmpty(i,j,k);
// 	if (!curFluid && !curEmpty) return;
// 	if (exclude && ((*exclude)(i,j,k) < 0.)) return;

// 	if (flags.isFluid(i-1,j,k) || (curFluid && flags.isEmpty(i-1,j,k))) 
// 		vel(i,j,k).x = (additive) ? vel(i,j,k).x+force.x : force.x;
// 	if (flags.isFluid(i,j-1,k) || (curFluid && flags.isEmpty(i,j-1,k))) 
// 		vel(i,j,k).y = (additive) ? vel(i,j,k).y+force.y : force.y;
// 	if (vel.is3D() && (flags.isFluid(i,j,k-1) || (curFluid && flags.isEmpty(i,j,k-1))))
// 		vel(i,j,k).z = (additive) ? vel(i,j,k).z+force.z : force.z;
// }

KERNEL(bnd=1) void KnApplyForce(const FlagGrid& flags, MACGrid& vel, Vec3 force, const Grid<Real>* exclude, bool additive) {
	bool curFluid = flags.isFluid(i,j,k);
	bool curSolid = flags.isSolid(i,j,k);
	bool curEmpty = flags.isEmpty(i,j,k);
	
	if (!curFluid && !curEmpty && !curSolid) return;
	if (exclude && ((*exclude)(i,j,k) < 0.)) return;

	// fl/fl and fl/em
	if (flags.isFluid(i-1,j,k) || (curFluid && flags.isEmpty(i-1,j,k))) 
		vel(i,j,k).x = (additive) ? vel(i,j,k).x+force.x : force.x;
	
	if (flags.isFluid(i,j-1,k) || (curFluid && flags.isEmpty(i,j-1,k))) 
		vel(i,j,k).y = (additive) ? vel(i,j,k).y+force.y : force.y;
	
	if (vel.is3D() && (flags.isFluid(i,j,k-1) || (curFluid && flags.isEmpty(i,j,k-1))))
		vel(i,j,k).z = (additive) ? vel(i,j,k).z+force.z : force.z;
	
	// solid/solid and solid/em
	
	if (flags.isSolid(i-1,j,k) || (curSolid && flags.isEmpty(i-1,j,k))) 
		vel(i,j,k).x = (additive) ? vel(i,j,k).x+force.x : force.x;
	
	if (flags.isSolid(i,j-1,k) || (curSolid && flags.isEmpty(i,j-1,k))) 
		vel(i,j,k).y = (additive) ? vel(i,j,k).y+force.y : force.y;
	
	if (vel.is3D() && (flags.isSolid(i,j,k-1) || (curSolid && flags.isEmpty(i,j,k-1))))
		vel(i,j,k).z = (additive) ? vel(i,j,k).z+force.z : force.z;

	// fluid/solid	
}

//! add gravity forces to all fluid cells, optionally  adapts to different grid sizes automatically
PYTHON() void addGravity(const FlagGrid& flags, MACGrid& vel, Vec3 gravity, const Grid<Real>* exclude=NULL, bool scale=true) {
	float gridScale = (scale) ? flags.getDx() : 1;
	Vec3 f = gravity * flags.getParent()->getDt() / gridScale;
	KnApplyForce(flags, vel, f, exclude, true);
}
//! Deprecated: use addGravity(scale=false) instead
PYTHON() void addGravityNoScale(const FlagGrid& flags, MACGrid& vel, const Vec3& gravity, const Grid<Real>* exclude=NULL) {
	addGravity(flags, vel, gravity, exclude, false);
}

//! kernel to add Buoyancy force 
KERNEL(bnd=1) void KnAddBuoyancy(const FlagGrid& flags, const Grid<Real>& factor, MACGrid& vel, Vec3 strength) {
	if (!flags.isFluid(i,j,k)) return;
	if (flags.isFluid(i-1,j,k))
		vel(i,j,k).x += (0.5 * strength.x) * (factor(i,j,k)+factor(i-1,j,k));
	if (flags.isFluid(i,j-1,k))
		vel(i,j,k).y += (0.5 * strength.y) * (factor(i,j,k)+factor(i,j-1,k));
	if (vel.is3D() && flags.isFluid(i,j,k-1))
		vel(i,j,k).z += (0.5 * strength.z) * (factor(i,j,k)+factor(i,j,k-1));
}

//! add Buoyancy force based on factor (e.g. smoke density), optionally adapts to different grid sizes automatically
PYTHON() void addBuoyancy(const FlagGrid& flags, const Grid<Real>& density, MACGrid& vel, Vec3 gravity, Real coefficient=1., bool scale=true) {
	float gridScale = (scale) ? flags.getDx() : 1;
	Vec3 f = -gravity * flags.getParent()->getDt() / gridScale * coefficient;
	KnAddBuoyancy(flags,density, vel, f);
}

// inflow / outflow boundaries

//! helper to parse openbounds string [xXyYzZ] , convert to vec3 
inline void convertDescToVec(const string& desc, Vector3D<bool>& lo, Vector3D<bool>& up) {
	for (size_t i = 0; i<desc.size(); i++) {
		if (desc[i] == 'x') lo.x = true;
		else if (desc[i] == 'y') lo.y = true;
		else if (desc[i] == 'z') lo.z = true;
		else if (desc[i] == 'X') up.x = true;
		else if (desc[i] == 'Y') up.y = true;
		else if (desc[i] == 'Z') up.z = true;
		else errMsg("invalid character in boundary description string. Only [xyzXYZ] allowed.");
	}
}

//! add empty and outflow flag to cells of open boundaries 
PYTHON() void setOpenBound(FlagGrid& flags, int bWidth, string openBound = "", int type = FlagGrid::TypeOutflow | FlagGrid::TypeEmpty){
	if (openBound == "") return;
	Vector3D<bool> lo, up;
	convertDescToVec(openBound, lo, up);

	FOR_IJK(flags){
		bool loX = lo.x && i <= bWidth; // a cell which belongs to the lower x open bound
		bool loY = lo.y && j <= bWidth; 
		bool upX = up.x && i >= flags.getSizeX() - bWidth - 1; // a cell which belongs to the upper x open bound
		bool upY = up.y && j >= flags.getSizeY() - bWidth - 1; 
		bool innerI = i>bWidth && i<flags.getSizeX() - bWidth - 1; // a cell which does not belong to the lower or upper x bound
		bool innerJ = j>bWidth && j<flags.getSizeY() - bWidth - 1; 

		// when setting boundaries to open: don't set shared part of wall to empty if neighboring wall is not open
		if ( (!flags.is3D()) && (loX||upX||loY||upY)){
			if ((loX || upX || innerI) && (loY || upY || innerJ) && flags.isObstacle(i, j, k)) flags(i, j, k) = type;
		} else {
			bool loZ = lo.z && k <= bWidth; // a cell which belongs to the lower z open bound
			bool upZ = up.z && k >= flags.getSizeZ() - bWidth - 1; // a cell which belongs to the upper z open bound
			bool innerK = k>bWidth && k<flags.getSizeZ() - bWidth - 1; // a cell which does not belong to the lower or upper z bound
			if (loX || upX || loY || upY || loZ || upZ) {
				if ((loX || upX || innerI) && (loY || upY || innerJ) && (loZ || upZ || innerK) && flags.isObstacle(i, j, k)) flags(i, j, k) = type;
			}
		}
	}
}

//! delete fluid and ensure empty flag in outflow cells, delete particles and density and set phi to 0.5
PYTHON() void resetOutflow(FlagGrid& flags, Grid<Real>* phi = 0, BasicParticleSystem* parts = 0, Grid<Real>* real = 0, Grid<int>* index = 0, ParticleIndexSystem* indexSys = 0){
	// check if phi and parts -> pindex and gpi already created -> access particles from cell index, avoid extra looping over particles
	if (parts && (!index || !indexSys)){
		if (phi) debMsg("resetOpenBound for phi and particles, but missing index and indexSys for enhanced particle access!",1);
		for (int idx = 0; idx < (int)parts->size(); idx++) 
			if (parts->isActive(idx) && flags.isInBounds(parts->getPos(idx)) && flags.isOutflow(parts->getPos(idx))) parts->kill(idx);
	}
	FOR_IJK(flags){
		if (flags.isOutflow(i,j,k)){
			flags(i, j, k) = (flags(i, j, k) | FlagGrid::TypeEmpty) & ~FlagGrid::TypeFluid; // make sure there is not fluid flag set and to reset the empty flag
			// the particles in a cell i,j,k are particles[index(i,j,k)] to particles[index(i+1,j,k)-1]
			if (parts && index && indexSys){
				int isysIdxS = index->index(i, j, k);
				int pStart = (*index)(isysIdxS), pEnd = 0;
				if (flags.isInBounds(isysIdxS + 1)) pEnd = (*index)(isysIdxS + 1);
				else								pEnd = indexSys->size();
				// now loop over particles in cell
				for (int p = pStart; p<pEnd; ++p) {
					int psrc = (*indexSys)[p].sourceIndex;
					if (parts->isActive(psrc) && flags.isInBounds(parts->getPos(psrc))) parts->kill(psrc);
				}
			}
			if (phi) (*phi)(i, j, k) = 0.5;
			if (real) (*real)(i, j, k) = 0;
		}
	}
	if (parts) parts->doCompress();
}

//! enforce a constant inflow/outflow at the grid boundaries
KERNEL() void KnSetInflow(MACGrid& vel, int dim, int p0, const Vec3& val) {
	Vec3i p(i,j,k);
	if (p[dim] == p0 || p[dim] == p0+1)
		vel(i,j,k) = val;
}

//! enforce a constant inflow/outflow at the grid boundaries
PYTHON() void setInflowBcs(MACGrid& vel, string dir, Vec3 value) {
	for(size_t i=0; i<dir.size(); i++) {
		if (dir[i] >= 'x' && dir[i] <= 'z') { 
			int dim = dir[i]-'x';
			KnSetInflow(vel,dim,0,value);
		} else if (dir[i] >= 'X' && dir[i] <= 'Z') {
			int dim = dir[i]-'X';
			KnSetInflow(vel,dim,vel.getSize()[dim]-1,value);
		} else 
			errMsg("invalid character in direction string. Only [xyzXYZ] allowed.");
	}
}

// set obstacle boundary conditions

//! set no-stick wall boundary condition between ob/fl and ob/ob cells
KERNEL() void KnSetWallBcs(const FlagGrid& flags, MACGrid& vel, const MACGrid* obvel) {

	bool curFluid = flags.isFluid(i,j,k);
	bool curSolid = flags.isSolid(i,j,k);
	bool curObs   = flags.isObstacle(i,j,k);
	Vec3 bcsVel(0.,0.,0.);
	
	if (!curFluid && !curSolid && !curObs) return;

	if (obvel) {
		bcsVel.x = (*obvel)(i,j,k).x;
		bcsVel.y = (*obvel)(i,j,k).y;
		if((*obvel).is3D()) bcsVel.z = (*obvel)(i,j,k).z;
	}

	// we use i>0 instead of bnd=1 to check outer wall
	if (i > 0 && flags.isObstacle(i-1,j,k))			vel(i,j,k).x = bcsVel.x;
	if (i > 0 && curObs && flags.isFluid(i-1,j,k))	vel(i,j,k).x = bcsVel.x;
	if (i > 0 && curObs && flags.isSolid(i-1,j,k))	vel(i,j,k).x = bcsVel.x;
	if (j > 0 && flags.isObstacle(i,j-1,k))			vel(i,j,k).y = bcsVel.y;
	if (j > 0 && curObs && flags.isFluid(i,j-1,k))	vel(i,j,k).y = bcsVel.y;
	if (j > 0 && curObs && flags.isSolid(i,j-1,k))	vel(i,j,k).y = bcsVel.y;

	if(!vel.is3D()) {vel(i,j,k).z = 0;} 
	else {
		if (k > 0 && flags.isObstacle(i,j,k-1))		 	vel(i,j,k).z = bcsVel.z;
		if (k > 0 && curObs && flags.isFluid(i,j,k-1)) 	vel(i,j,k).z = bcsVel.z;
		if (k > 0 && curObs && flags.isSolid(i,j,k-1)) 	vel(i,j,k).z = bcsVel.z;
	}
	
	if (curFluid || curSolid) {
		if ((i>0 && flags.isStick(i-1,j,k)) || (i<flags.getSizeX()-1 && flags.isStick(i+1,j,k)))
			vel(i,j,k).y = vel(i,j,k).z = 0;
		if ((j>0 && flags.isStick(i,j-1,k)) || (j<flags.getSizeY()-1 && flags.isStick(i,j+1,k)))
			vel(i,j,k).x = vel(i,j,k).z = 0;
		if (vel.is3D() && ((k>0 && flags.isStick(i,j,k-1)) || (k<flags.getSizeZ()-1 && flags.isStick(i,j,k+1))))
			vel(i,j,k).x = vel(i,j,k).y = 0;
	}
}

// KERNEL() void KnSetWallBcs(const FlagGrid& flags, MACGrid& vel, const MACGrid* obvel) {

// 	bool curFluid = flags.isFluid(i,j,k);
// 	bool curObs   = flags.isObstacle(i,j,k);
// 	Vec3 bcsVel(0.,0.,0.);
// 	if (!curFluid && !curObs) return;

// 	if (obvel) {
// 		bcsVel.x = (*obvel)(i,j,k).x;
// 		bcsVel.y = (*obvel)(i,j,k).y;
// 		if((*obvel).is3D()) bcsVel.z = (*obvel)(i,j,k).z;
// 	}

// 	// we use i>0 instead of bnd=1 to check outer wall
// 	if (i>0 && flags.isObstacle(i-1,j,k))						 vel(i,j,k).x = bcsVel.x;
// 	if (i>0 && curObs && flags.isFluid(i-1,j,k))				 vel(i,j,k).x = bcsVel.x;
// 	if (j>0 && flags.isObstacle(i,j-1,k))						 vel(i,j,k).y = bcsVel.y;
// 	if (j>0 && curObs && flags.isFluid(i,j-1,k))				 vel(i,j,k).y = bcsVel.y;

// 	if(!vel.is3D()) {                            				vel(i,j,k).z = 0; } else {
// 	if (k>0 && flags.isObstacle(i,j,k-1))		 				vel(i,j,k).z = bcsVel.z;
// 	if (k>0 && curObs && flags.isFluid(i,j,k-1)) 				vel(i,j,k).z = bcsVel.z; }
	
// 	if (curFluid) {
// 		if ((i>0 && flags.isStick(i-1,j,k)) || (i<flags.getSizeX()-1 && flags.isStick(i+1,j,k)))
// 			vel(i,j,k).y = vel(i,j,k).z = 0;
// 		if ((j>0 && flags.isStick(i,j-1,k)) || (j<flags.getSizeY()-1 && flags.isStick(i,j+1,k)))
// 			vel(i,j,k).x = vel(i,j,k).z = 0;
// 		if (vel.is3D() && ((k>0 && flags.isStick(i,j,k-1)) || (k<flags.getSizeZ()-1 && flags.isStick(i,j,k+1))))
// 			vel(i,j,k).x = vel(i,j,k).y = 0;
// 	}
// }

//! set wall BCs for fill fraction mode, note - only needs obstacle SDF
KERNEL() void KnSetWallBcsFrac(const FlagGrid& flags, const MACGrid& vel, MACGrid& velTarget, const MACGrid* obvel,
							const Grid<Real>* phiObs, const int &boundaryWidth=0) 
{ 
	bool curFluid = flags.isFluid(i,j,k);
	bool curObs   = flags.isObstacle(i,j,k);
	velTarget(i,j,k) = vel(i,j,k);
	if (!curFluid && !curObs) return;

	// zero normal component in all obstacle regions
	if(flags.isInBounds(Vec3i(i,j,k),1)) {

	if( curObs | flags.isObstacle(i-1,j,k) )  { 
		Vec3 dphi(0.,0.,0.);
		const Real tmp1 = (phiObs->get(i,j,k)+phiObs->get(i-1,j,k))*.5;
		Real tmp2 = (phiObs->get(i,j+1,k)+phiObs->get(i-1,j+1,k))*.5;
		Real phi1 = (tmp1+tmp2)*.5;
		tmp2 = (phiObs->get(i,j-1,k)+phiObs->get(i-1,j-1,k))*.5;
		Real phi2 = (tmp1+tmp2)*.5;
		
		dphi.x = phiObs->get(i,j,k)-phiObs->get(i-1,j,k);
		dphi.y = phi1-phi2;

		if(phiObs->is3D()) {
			tmp2 = (phiObs->get(i,j,k+1)+phiObs->get(i-1,j,k+1))*.5;
			phi1 = (tmp1+tmp2)*.5;
			tmp2 = (phiObs->get(i,j,k-1)+phiObs->get(i-1,j,k-1))*.5;
			phi2 = (tmp1+tmp2)*.5;
			dphi.z = phi1-phi2;
		}

		normalize(dphi); 
		Vec3 velMAC = vel.getAtMACX(i,j,k);
		velTarget(i,j,k).x = velMAC.x - dot(dphi, velMAC) * dphi.x;
	}

	if( curObs | flags.isObstacle(i,j-1,k) )  { 
		Vec3 dphi(0.,0.,0.);
		const Real tmp1 = (phiObs->get(i,j,k)+phiObs->get(i,j-1,k))*.5;
		Real tmp2 = (phiObs->get(i+1,j,k)+phiObs->get(i+1,j-1,k))*.5;
		Real phi1 = (tmp1+tmp2)*.5;
		tmp2 = (phiObs->get(i-1,j,k)+phiObs->get(i-1,j-1,k))*.5;
		Real phi2 = (tmp1+tmp2)*.5;

		dphi.x = phi1-phi2;
		dphi.y = phiObs->get(i,j,k)-phiObs->get(i,j-1,k);
		if(phiObs->is3D()) {
			tmp2 = (phiObs->get(i,j,k+1)+phiObs->get(i,j-1,k+1))*.5;
			phi1 = (tmp1+tmp2)*.5;
			tmp2 = (phiObs->get(i,j,k-1)+phiObs->get(i,j-1,k-1))*.5;
			phi2 = (tmp1+tmp2)*.5;
			dphi.z = phi1-phi2;
		}

		normalize(dphi); 
		Vec3 velMAC = vel.getAtMACY(i,j,k);
		velTarget(i,j,k).y = velMAC.y - dot(dphi, velMAC) * dphi.y;
	}

	if( phiObs->is3D() && (curObs | flags.isObstacle(i,j,k-1)) )  {
		Vec3 dphi(0.,0.,0.); 
		const Real tmp1 = (phiObs->get(i,j,k)+phiObs->get(i,j,k-1))*.5;

		Real tmp2;
		tmp2      = (phiObs->get(i+1,j,k)+phiObs->get(i+1,j,k-1))*.5;
		Real phi1 = (tmp1+tmp2)*.5;
		tmp2      = (phiObs->get(i-1,j,k)+phiObs->get(i-1,j,k-1))*.5;
		Real phi2 = (tmp1+tmp2)*.5; 
		dphi.x    = phi1-phi2;

		tmp2      = (phiObs->get(i,j+1,k)+phiObs->get(i,j+1,k-1))*.5;
		phi1      = (tmp1+tmp2)*.5;
		tmp2      = (phiObs->get(i,j-1,k)+phiObs->get(i,j-1,k-1))*.5;
		phi2      = (tmp1+tmp2)*.5; 
		dphi.y    = phi1-phi2;

		dphi.z = phiObs->get(i,j,k) - phiObs->get(i,j,k-1);

		normalize(dphi); 
		Vec3 velMAC = vel.getAtMACZ(i,j,k);
		velTarget(i,j,k).z = velMAC.z - dot(dphi, velMAC) * dphi.z;
	}
	} // not at boundary

}

//! set zero normal velocity boundary condition on walls
// (optionally with second order accuracy using the obstacle SDF , fractions grid currentlyl not needed)
PYTHON() void setWallBcs(const FlagGrid& flags, MACGrid& vel, const MACGrid* obvel = 0, const MACGrid* fractions = 0, const Grid<Real>* phiObs = 0, int boundaryWidth=0) {
	if(!phiObs || !fractions) {
		KnSetWallBcs(flags, vel, obvel);
	} else {
		MACGrid tmpvel(vel.getParent());
		KnSetWallBcsFrac(flags, vel, tmpvel, obvel, phiObs, boundaryWidth);
		vel.swap(tmpvel);
	}
}

//! add Forces between fl/fl and fl/em cells (interpolate cell centered forces to MAC grid)
KERNEL(bnd=1) void KnAddForceIfLower(const FlagGrid& flags, MACGrid& vel, const Grid<Vec3>& force) {
	bool curFluid = flags.isFluid(i,j,k);
	bool curEmpty = flags.isEmpty(i,j,k);
	if (!curFluid && !curEmpty) return;

	if (flags.isFluid(i-1,j,k) || (curFluid && flags.isEmpty(i-1,j,k))) {
		Real forceMACX = 0.5*(force(i-1,j,k).x + force(i,j,k).x);
		Real min = std::min(vel(i,j,k).x, forceMACX);
		Real max = std::max(vel(i,j,k).x, forceMACX);
		Real sum = vel(i,j,k).x + forceMACX;
		vel(i,j,k).x = (forceMACX > 0) ? std::min(sum, max) : std::max(sum, min);
	}
	if (flags.isFluid(i,j-1,k) || (curFluid && flags.isEmpty(i,j-1,k))) {
		Real forceMACY = 0.5*(force(i,j-1,k).y + force(i,j,k).y);
		Real min = std::min(vel(i,j,k).y, forceMACY);
		Real max = std::max(vel(i,j,k).y, forceMACY);
		Real sum = vel(i,j,k).y + forceMACY;
		vel(i,j,k).y = (forceMACY > 0) ? std::min(sum, max) : std::max(sum, min);
	}
	if (vel.is3D() && (flags.isFluid(i,j,k-1) || (curFluid && flags.isEmpty(i,j,k-1)))) {
		Real forceMACZ = 0.5*(force(i,j,k-1).z + force(i,j,k).z);
		Real min = std::min(vel(i,j,k).z, forceMACZ);
		Real max = std::max(vel(i,j,k).z, forceMACZ);
		Real sum = vel(i,j,k).z + forceMACZ;
		vel(i,j,k).z = (forceMACZ > 0) ? std::min(sum, max) : std::max(sum, min);
	}
}

// Initial velocity for smoke
PYTHON() void setInitialVelocity(const FlagGrid& flags, MACGrid& vel, const Grid<Vec3>& invel) {
	KnAddForceIfLower(flags, vel, invel);
}

//! Kernel: gradient norm operator
KERNEL(bnd=1) void KnConfForce(Grid<Vec3>& force, const Grid<Real>& grid, const Grid<Vec3>& curl, Real str, const Grid<Real>* strGrid) {
	Vec3 grad = 0.5 * Vec3(        grid(i+1,j,k)-grid(i-1,j,k), 
								   grid(i,j+1,k)-grid(i,j-1,k), 0.);
	if(grid.is3D()) grad[2]= 0.5*( grid(i,j,k+1)-grid(i,j,k-1) );
	normalize(grad);
	if (strGrid) str += (*strGrid)(i,j,k);
	force(i,j,k) = str * cross(grad, curl(i,j,k));
}

PYTHON() void vorticityConfinement(MACGrid& vel, const FlagGrid& flags, Real strength=0, const Grid<Real>* strengthCell=NULL) {
	Grid<Vec3> velCenter(flags.getParent()), curl(flags.getParent()), force(flags.getParent());
	Grid<Real> norm(flags.getParent());
	
	GetCentered(velCenter, vel);
	CurlOp(velCenter, curl);
	GridNorm(norm, curl);
	KnConfForce(force, norm, curl, strength, strengthCell);
	KnApplyForceField(flags, vel, force, NULL, true, false);
}

PYTHON() void addForceField(const FlagGrid& flags, MACGrid& vel, const Grid<Vec3>& force, const Grid<Real>* region=NULL, bool isMAC=false) {
	KnApplyForceField(flags, vel, force, region, true, isMAC);
}

PYTHON() void setForceField(const FlagGrid& flags, MACGrid& vel, const Grid<Vec3>& force, const Grid<Real>* region=NULL, bool isMAC=false) {
	KnApplyForceField(flags, vel, force, region, false, isMAC);
}

KERNEL() void KnDissolveSmoke(const FlagGrid& flags, Grid<Real>& density, Grid<Real>* heat,
	Grid<Real>* red, Grid<Real>* green, Grid<Real>* blue, int speed, bool logFalloff, float dydx, float fac) {

	bool curFluid = flags.isFluid(i,j,k);
	if (!curFluid) return;

	if (logFalloff) {
		density(i,j,k) *= fac;
		if (heat) {
			(*heat)(i,j,k) *= fac;
		}
		if (red) {
			(*red)(i,j,k) *= fac;
			(*green)(i,j,k) *= fac;
			(*blue)(i,j,k) *= fac;
		}
	}
	else { // linear falloff
		float d = density(i,j,k);
		density(i,j,k) -= dydx;
		if (density(i,j,k) < 0.0f)
			density(i,j,k) = 0.0f;
		if (heat) {
			if      (fabs((*heat)(i,j,k)) < dydx) (*heat)(i,j,k) = 0.0f;
			else if ((*heat)(i,j,k) > 0.0f) (*heat)(i,j,k) -= dydx;
			else if ((*heat)(i,j,k) < 0.0f) (*heat)(i,j,k) += dydx;
		}
		if (red && notZero(d) ) {
			(*red)(i,j,k)   *= (density(i,j,k)/d);
			(*green)(i,j,k) *= (density(i,j,k)/d);
			(*blue)(i,j,k)  *= (density(i,j,k)/d);
		}
	}
}

PYTHON() void dissolveSmoke(const FlagGrid& flags, Grid<Real>& density, Grid<Real>* heat=NULL,
	Grid<Real>* red=NULL, Grid<Real>* green=NULL, Grid<Real>* blue=NULL, int speed=5, bool logFalloff=true)
{
	float dydx = 1.0f / (float)speed; // max density/speed = dydx
	float fac = 1.0f - dydx;
	KnDissolveSmoke(flags, density, heat, red, green, blue, speed, logFalloff, dydx, fac);
}

} // namespace
