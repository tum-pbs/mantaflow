/******************************************************************************
 *
 * Plugins for Position-Based Dynamics (PDB)
 * Copyright 2020 Cesar Rodriguez, James Mossell 
 * 
 * The implementation is based on the publication An efficient FLIP and shape
 * matching coupled method for fluid-solid and two-phase fluid simulations by
 * Yang Gao, Shuai Li, Hong Qin, Yinghao Xu, and Aimin Hao
 * https://doi.org/10.1007/s00371-018-1569-8
 *
 * This source code is free and distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 *
 ******************************************************************************/

#include "vectorbase.h"
#include "grid.h"
#include "kernel.h"
#include <limits>

using namespace std;

namespace Manta { 

	//! add constant force between fl/fl and fl/em cells
	KERNEL(bnd=1) void KnApplyForce(const FlagGrid& flags, MACGrid& vel, Vec3 force, const Grid<Real>* exclude, bool additive) {
		bool curFluid = flags.isFluid(i,j,k);
		bool curEmpty = flags.isEmpty(i,j,k);
		if (!curFluid && !curEmpty) return;
		if (exclude && ((*exclude)(i,j,k) < 0.)) return;

		if (flags.isFluid(i-1,j,k) || (curFluid && flags.isEmpty(i-1,j,k))) 
			vel(i,j,k).x = (additive) ? vel(i,j,k).x+force.x : force.x;
		if (flags.isFluid(i,j-1,k) || (curFluid && flags.isEmpty(i,j-1,k))) 
			vel(i,j,k).y = (additive) ? vel(i,j,k).y+force.y : force.y;
		if (vel.is3D() && (flags.isFluid(i,j,k-1) || (curFluid && flags.isEmpty(i,j,k-1))))
			vel(i,j,k).z = (additive) ? vel(i,j,k).z+force.z : force.z;
	}

	//! add gravity forces to all fluid cells, optionally  adapts to different grid sizes automatically
	PYTHON() void addGravity(const FlagGrid& flags, MACGrid& vel, Vec3 gravity, const Grid<Real>* exclude=NULL, bool scale=true) {
		float gridScale = (scale) ? flags.getDx() : 1;
		Vec3 f = gravity * flags.getParent()->getDt() / gridScale;
		KnApplyForce(flags, vel, f, exclude, true);
	}


	//! Semi-Lagrange interpolation kernel
	KERNEL(bnd=1) template<class T> 
	void SemiLagrange (const FlagGrid& flags, const MACGrid& vel, Grid<T>& dst, const Grid<T>& src, Real dt, bool isLevelset, int orderSpace, int orderTrace)
	{
		if (orderTrace == 1) {
			// traceback position
			Vec3 pos = Vec3(i+0.5f,j+0.5f,k+0.5f) - vel.getCentered(i,j,k) * dt;
			dst(i,j,k) = src.getInterpolatedHi(pos, orderSpace);
		} else if (orderTrace == 2) {
			// backtracing using explicit midpoint
			Vec3 p0 = Vec3(i+0.5f,j+0.5f,k+0.5f);
			Vec3 p1 = p0 - vel.getCentered(i,j,k)*dt*0.5;
			Vec3 p2 = p0 - vel.getInterpolated(p1)*dt;
			dst(i,j,k) = src.getInterpolatedHi(p2, orderSpace);
		} else {
			assertMsg(false, "Unknown backtracing order "<<orderTrace);
		}
	}

	//! Perform semi-lagrangian advection of target Real- or Vec3 grid
	//! Open boundary handling needs information about width of border
	//! Clamping modes: 1 regular clamp leading to more overshoot and sharper results, 2 revert to 1st order slightly smoother less overshoot (enable when 1 gives artifacts)
	PYTHON() void advectSemiLagrange (const FlagGrid* flags, const MACGrid* vel, GridBase* grid,
									int order = 1, Real strength = 1.0, int orderSpace = 1, bool openBounds = false, int boundaryWidth = -1, int clampMode = 2, int orderTrace = 1)
	{    
		assertMsg(order==1 || order==2, "AdvectSemiLagrange: Only order 1 (regular SL) and 2 (MacCormack) supported");
		if((boundaryWidth!=-1)||(openBounds)) { debMsg("Warning: boundaryWidth and openBounds parameters in AdvectSemiLagrange plugin are deprecated (and have no more effect), please remove.", 0); } 
		
		// determine type of grid    
		if (grid->getType() & GridBase::TypeReal) {
			fnAdvectSemiLagrange< Grid<Real> >(flags->getParent(), *flags, *vel, *((Grid<Real>*) grid), order, strength, orderSpace, clampMode, orderTrace);
		}
		else if (grid->getType() & GridBase::TypeMAC) {    
			fnAdvectSemiLagrange< MACGrid >(flags->getParent(), *flags, *vel, *((MACGrid*) grid), order, strength, orderSpace, clampMode, orderTrace);
		}
		else if (grid->getType() & GridBase::TypeVec3) {    
			fnAdvectSemiLagrange< Grid<Vec3> >(flags->getParent(), *flags, *vel, *((Grid<Vec3>*) grid), order, strength, orderSpace, clampMode, orderTrace);
		}
		else
			errMsg("AdvectSemiLagrange: Grid Type is not supported (only Real, Vec3, MAC, Levelset)");    
	}

} // end namespace DDF 

