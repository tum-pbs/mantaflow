/******************************************************************************
 *
 * Plugins for Position-Based Dynamics (PBD)
 * Copyright 2020 Cesar Rodriguez, James Mossell 
 * 
 * The implementation is based on the publication: 
 * 
 * Gao, Y., Li, S., Qin, H. et al. An efficient FLIP and shape matching 
 * coupled method for fluid–solid and two-phase fluid simulations. 
 * Vis Comput 35, 1741–1753 (2019). 
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
	//--------------------- Kernel Plugins ---------------------------

	// KERNEL (bnd=boundaryWidth)
	// void knUpdateFlagsSolid(FlagGrid& flags, const Grid<Real>& phiSolid, int boundaryWidth) {
	// 	bool isSolid;
		
	// 	if(phiSolid(i,j,k) < 0.) isSolid = true;

	// 	if (isSolid) flags(i,j,k) = FlagGrid::TypeSolid;

	// }

	KERNEL (bnd=boundaryWidth)
	void KnUpdateFlagsSolid(FlagGrid& flags, const MACGrid* fractions, const Grid<Real>& phiSolid, const Grid<Real>* phiOut, const Grid<Real>* phiIn, int boundaryWidth) {

		bool isSolid = false;
		
		if(fractions) {
			Real f = 0.;
			f += fractions->get(i  ,j,k).x;
			f += fractions->get(i+1,j,k).x;
			f += fractions->get(i,j  ,k).y;
			f += fractions->get(i,j+1,k).y;
			if (flags.is3D()) {
			f += fractions->get(i,j,k  ).z;
			f += fractions->get(i,j,k+1).z; }
			if(f==0.) isSolid = true;
		} else {
			if(phiSolid(i,j,k) < 0.) isSolid = true;
		}

		bool isOutflow = false;
		bool isInflow = false;
		if (phiOut && (*phiOut)(i,j,k) < 0.) isOutflow = true;
		if (phiIn && (*phiIn)(i,j,k) < 0.) isInflow = true;

		if (isSolid)          flags(i,j,k) = FlagGrid::TypeSolid;
		else if (isInflow)  flags(i,j,k) = (FlagGrid::TypeFluid | FlagGrid::TypeInflow);
		else if (isOutflow) flags(i,j,k) = (FlagGrid::TypeEmpty | FlagGrid::TypeOutflow);
		else                flags(i,j,k) = FlagGrid::TypeEmpty;
	}

	//! update obstacle and outflow flags from levelsets
	//! optionally uses fill fractions for obstacle
	PYTHON() void setSolidFlags(FlagGrid& flags, const Grid<Real>& phiSolid, const MACGrid* fractions=NULL, const Grid<Real>* phiOut=NULL, const Grid<Real>* phiIn=NULL, int boundaryWidth=1){
		KnUpdateFlagsSolid(flags, fractions, phiSolid, phiOut, phiIn, boundaryWidth);
	}

	// Update solid flags from levelsets
	// PYTHON() void setSolidFlags(FlagGrid& flags, const Grid<Real>& phiSolid, int boundaryWidth=1) {
	// 	knUpdateFlagsSolid(flags, phiSolid, boundaryWidth);
	// }

	// PYTHON() bool initShapeMatchingConstraint(
	// const Vector3r x0[], const Real invMasses[], int numPoints,
	// Vector3r &restCm, Matrix3r &invRestMat)
	// {
	// 	invRestMat.setIdentity();

	// 	// center of mass
	// 	restCm.setZero();
	// 	Real wsum = 0.0;
	// 	for (int i = 0; i < numPoints; i++) {
	// 		Real wi = static_cast<Real>(1.0) / (invMasses[i] + eps);
	// 		restCm += x0[i] * wi;
	// 		wsum += wi;
	// 	}
	// 	if (wsum == 0.0)
	// 		return false;
	// 	restCm /= wsum;

	// 	// A
	// 	Matrix3r A;
	// 	A.setZero();
	// 	for (int i = 0; i < numPoints; i++) {
	// 		const Vector3r qi = x0[i] - restCm;
	// 		Real wi = static_cast<Real>(1.0) / (invMasses[i] + eps);
	// 		Real x2 = wi * qi[0] * qi[0];
	// 		Real y2 = wi * qi[1] * qi[1];
	// 		Real z2 = wi * qi[2] * qi[2];
	// 		Real xy = wi * qi[0] * qi[1];
	// 		Real xz = wi * qi[0] * qi[2];
	// 		Real yz = wi * qi[1] * qi[2];
	// 		A(0, 0) += x2; A(0, 1) += xy; A(0, 2) += xz;
	// 		A(1, 0) += xy; A(1, 1) += y2; A(1, 2) += yz;
	// 		A(2, 0) += xz; A(2, 1) += yz; A(2, 2) += z2;
	// 	}
	// 	Real det = A.determinant();
	// 	if (fabs(det) > eps)
	// 	{
	// 		invRestMat = A.inverse();
	// 		return true;
	// 	}
	// 	return false;
	// }

	// PYTHON() bool solveShapeMatchingConstraint(
	// 	const Vector3r x0[], const Vector3r x[], const Real invMasses[], int numPoints,
	// 	const Vector3r &restCm, 
	// 	const Matrix3r &invRestMat,
	// 	const Real stiffness,
	// 	const bool allowStretch,
	// 	Vector3r corr[], Matrix3r *rot)
	// {
	// 	for (int i = 0; i < numPoints; i++)
	// 		corr[i].setZero();

	// 	// center of mass
	// 	Vector3r cm(0.0, 0.0, 0.0);
	// 	Real wsum = 0.0;
	// 	for (int i = 0; i < numPoints; i++) 
	// 	{
	// 		Real wi = static_cast<Real>(1.0) / (invMasses[i] + eps);
	// 		cm += x[i] * wi;
	// 		wsum += wi;
	// 	}
	// 	if (wsum == 0.0)
	// 		return false;
	// 	cm /= wsum;

	// 	// A
	// 	Matrix3r mat;
	// 	mat.setZero();
	// 	for (int i = 0; i < numPoints; i++) {
	// 		Vector3r q = x0[i] - restCm;
	// 		Vector3r p = x[i] - cm;

	// 		Real w = static_cast<Real>(1.0) / (invMasses[i] + eps);
	// 		p *= w;

	// 		mat(0, 0) += p[0] * q[0]; mat(0, 1) += p[0] * q[1]; mat(0, 2) += p[0] * q[2];
	// 		mat(1, 0) += p[1] * q[0]; mat(1, 1) += p[1] * q[1]; mat(1, 2) += p[1] * q[2];
	// 		mat(2, 0) += p[2] * q[0]; mat(2, 1) += p[2] * q[1]; mat(2, 2) += p[2] * q[2];
	// 	}

	// 	mat = mat * invRestMat;

	// 	Matrix3r R, U, D;
	// 	R = mat;
	// 	if (allowStretch)
	// 		R = mat;
	// 	else
	// 		//MathFunctions::polarDecomposition(mat, R, U, D);
	// 		MathFunctions::polarDecompositionStable(mat, eps, R);

	// 	for (int i = 0; i < numPoints; i++) {
	// 		Vector3r goal = cm + R * (x0[i] - restCm);
	// 		corr[i] = (goal - x[i]) * stiffness;
	// 	}

	// 	if (rot)
	// 		*rot = R;

	// 	return true;
	// }

} // end namespace DDF 

