/******************************************************************************
*
* MantaFlow fluid solver framework
* Copyright 2017 Georg Kohl, Nils Thuerey
*
* This program is free software, distributed under the terms of the
* GNU General Public License (GPL)
* http://www.gnu.org/licenses
*
* Secondary particle plugin for FLIP simulations
*
******************************************************************************/

#include "particle.h"
#include "commonkernels.h"

namespace Manta {

#pragma region Secondary Particles for FLIP
//----------------------------------------------------------------------------------------------------------------------------------------------------
// Secondary Particles for FLIP
//----------------------------------------------------------------------------------------------------------------------------------------------------

// helper function that clamps the value in potential to the interval [tauMin, tauMax] and normalizes it to [0, 1] afterwards
Real clampPotential(Real potential, Real tauMin, Real tauMax) {
	return (std::min(potential, tauMax) - std::min(potential, tauMin)) / (tauMax - tauMin);
}

// computes all three potentials(trapped air, wave crest, kinetic energy) and the neighbor ratio for every fluid cell and stores it in the respective grid.
// Is less readable but significantly faster than using separate potential computation
KERNEL(bnd = radius)
void knFlipComputeSecondaryParticlePotentials(
	Grid<Real> &potTA, Grid<Real> &potWC, Grid<Real> &potKE, Grid<Real> &neighborRatio, const FlagGrid &flags, const MACGrid &v, const Grid<Vec3> &normal,
	const int radius, const Real tauMinTA, const Real tauMaxTA, const Real tauMinWC, const Real tauMaxWC, const Real tauMinKE, const Real tauMaxKE,
	const Real scaleFromManta, const int itype = FlagGrid::TypeFluid, const int jtype = FlagGrid::TypeObstacle | FlagGrid::TypeOutflow | FlagGrid::TypeInflow) {

	if (!(flags(i, j, k) & itype)) return;

	//compute trapped air potential + wave crest potential + neighbor ratio at once
	const Vec3 &xi = scaleFromManta * Vec3(i, j, k);	//scale to unit cube
	const Vec3 &vi = scaleFromManta * v.getCentered(i, j, k);
	const Vec3 &ni = getNormalized(normal(i, j, k));
	Real vdiff = 0;			//for trapped air
	Real kappa = 0;			//for wave crests
	int countFluid = 0;		//for neighbor ratio
	int countMaxFluid = 0;	//for neighbor ratio

	//iterate over neighboring cells within radius
	for (IndexInt x = i - radius; x <= i + radius; x++) {
		for (IndexInt y = j - radius; y <= j + radius; y++) {
			for (IndexInt z = k - radius; z <= k + radius; z++) {
				//ensure that xyz is in bounds: use bnd=1 to ensure that vel.getCentered() always has a neighbor cell
				if ((x == i && y == j && z == k) || !flags.isInBounds(Vec3i(x, y, z), 1) || (flags(x, y, z) & jtype)) continue;

				if (flags(x, y, z) & itype) {
					countFluid++;
					countMaxFluid++;
				}
				else {
					countMaxFluid++;
				}

				const Vec3 &xj = scaleFromManta * Vec3(x, y, z); //scale to unit cube
				const Vec3 &vj = scaleFromManta * v.getCentered(x, y, z);
				const Vec3 &nj = getNormalized(normal(x, y, z));
				const Vec3 xij = xi - xj;
				const Vec3 vij = vi - vj;
				Real h = !potTA.is3D() ? 1.414*radius : 1.732*radius; //estimate sqrt(2)*radius resp. sqrt(3)*radius for h, due to squared resp. cubic neighbor area
				vdiff += norm(vij) * (1 - dot(getNormalized(vij), getNormalized(xij))) * (1 - norm(xij) / h);

				if (dot(getNormalized(xij), ni) < 0) {	//identifies wave crests
					kappa += (1 - dot(ni, nj)) * (1 - norm(xij) / h);
				}
			}
		}
	}

	neighborRatio(i, j, k) = float(countFluid) / float(countMaxFluid);

	potTA(i, j, k) = clampPotential(vdiff, tauMinTA, tauMaxTA);
	if (dot(getNormalized(vi), ni) >= 0.6) {	//avoid to mark boarders of the scene as wave crest
		potWC(i, j, k) = clampPotential(kappa, tauMinWC, tauMaxWC);
	}
	else {
		potWC(i, j, k) = Real(0);
	}

	//compute kinetic energy potential
	Real ek = Real(0.5) * 125 * normSquare(vi);	//use arbitrary constant for mass, potential adjusts with thresholds anyways
	potKE(i, j, k) = clampPotential(ek, tauMinKE, tauMaxKE);
}
PYTHON()
void flipComputeSecondaryParticlePotentials(
	Grid<Real> &potTA, Grid<Real> &potWC, Grid<Real> &potKE, Grid<Real> &neighborRatio, const FlagGrid &flags, const MACGrid &v, Grid<Vec3>& normal, const Grid<Real>& phi,
	const int radius, const Real tauMinTA, const Real tauMaxTA, const Real tauMinWC, const Real tauMaxWC, const Real tauMinKE, const Real tauMaxKE,
	const Real scaleFromManta, const int itype = FlagGrid::TypeFluid, const int jtype = FlagGrid::TypeObstacle | FlagGrid::TypeOutflow | FlagGrid::TypeInflow) {
	potTA.clear();
	potWC.clear();
	potKE.clear();
	neighborRatio.clear();
	GradientOp(normal, phi);
	knFlipComputeSecondaryParticlePotentials(potTA, potWC, potKE, neighborRatio, flags, v, normal, radius, tauMinTA, tauMaxTA, tauMinWC, tauMaxWC, tauMinKE, tauMaxKE, scaleFromManta, itype, jtype);
}

// adds secondary particles to &pts_sec for every fluid cell in &flags according to the potential grids &potTA, &potWC and &potKE
// secondary particles are uniformly sampled in every fluid cell in a randomly offset cylinder in fluid movement direction
// In contrast to flipSampleSecondaryParticles this uses more cylinders per cell and interpolates velocity and potentials.
// To control number of cylinders in each dimension adjust radius(0.25=>2 cyl, 0.1666=>3 cyl, 0.125=>3cyl etc.).
KERNEL(single)
void knFlipSampleSecondaryParticlesMoreCylinders(
	const FlagGrid &flags, const MACGrid &v, BasicParticleSystem &pts_sec,
	ParticleDataImpl<Vec3> &v_sec, ParticleDataImpl<Real> &l_sec, const Real lMin, const Real lMax,
	const Grid<Real> &potTA, const Grid<Real> &potWC, const Grid<Real> &potKE, const Grid<Real> &neighborRatio,
	const Real c_s, const Real c_b, const Real k_ta, const Real k_wc, const Real dt, const int itype = FlagGrid::TypeFluid) {

	if (!(flags(i, j, k) & itype)) return;

	static RandomStream mRand(9832);
	Real radius = 0.25;	//diameter=0.5 => sampling with two cylinders in each dimension since cell size=1
	for (Real x = i - radius; x <= i + radius; x += 2 * radius) {
		for (Real y = j - radius; y <= j + radius; y += 2 * radius) {
			for (Real z = k - radius; z <= k + radius; z += 2 * radius) {

				Vec3 xi = Vec3(x, y, z);
				Real KE = potKE.getInterpolated(xi);
				Real TA = potTA.getInterpolated(xi);
				Real WC = potWC.getInterpolated(xi);

				const int n = KE * (k_ta*TA + k_wc*WC) * dt;		//number of secondary particles
				if (n == 0) continue;
				Vec3 vi = v.getInterpolated(xi);
				Vec3 dir = dt*vi;									//direction of movement of current particle
				Vec3 e1 = getNormalized(Vec3(dir.z, 0, -dir.x));	//perpendicular to dir
				Vec3 e2 = getNormalized(cross(e1, dir));			//perpendicular to dir and e1, so e1 and e1 create reference plane

				for (int di = 0; di < n; di++) {
					const Real r = radius * sqrt(mRand.getReal());			//distance to cylinder axis
					const Real theta = mRand.getReal() * Real(2) * M_PI;	//azimuth
					const Real h = mRand.getReal() * norm(dt*vi);			//distance to reference plane
					Vec3 xd = xi + r*cos(theta)*e1 + r*sin(theta)*e2 + h*getNormalized(vi);
					if (!flags.is3D()) xd.z = 0;
					pts_sec.add(xd);

					v_sec[v_sec.size() - 1] = r*cos(theta)*e1 + r*sin(theta)*e2 + vi;	//init velocity of new particle
					Real temp = (KE + TA + WC) / 3;
					l_sec[l_sec.size() - 1] = ((lMax - lMin) * temp) + lMin + mRand.getReal()*0.1;	//init lifetime of new particle

					//init type of new particle
					if (neighborRatio(i, j, k) < c_s) { pts_sec[pts_sec.size() - 1].flag = ParticleBase::PSPRAY; }
					else if (neighborRatio(i, j, k) > c_b) { pts_sec[pts_sec.size() - 1].flag = ParticleBase::PBUBBLE; }
					else { pts_sec[pts_sec.size() - 1].flag = ParticleBase::PFOAM; }
				}
			}
		}
	}
}

// adds secondary particles to &pts_sec for every fluid cell in &flags according to the potential grids &potTA, &potWC and &potKE
// secondary particles are uniformly sampled in every fluid cell in a randomly offset cylinder in fluid movement direction
KERNEL(single)
void knFlipSampleSecondaryParticles(
	const FlagGrid &flags, const MACGrid &v, BasicParticleSystem &pts_sec,
	ParticleDataImpl<Vec3> &v_sec, ParticleDataImpl<Real> &l_sec, const Real lMin, const Real lMax,
	const Grid<Real> &potTA, const Grid<Real> &potWC, const Grid<Real> &potKE, const Grid<Real> &neighborRatio,
	const Real c_s, const Real c_b, const Real k_ta, const Real k_wc, const Real dt, const int itype = FlagGrid::TypeFluid) {

	if (!(flags(i, j, k) & itype)) return;

	Real KE = potKE(i, j, k);
	Real TA = potTA(i, j, k);
	Real WC = potWC(i, j, k);

	const int n = KE * (k_ta*TA + k_wc*WC) * dt;		//number of secondary particles
	if (n == 0) return;
	static RandomStream mRand(9832);

	Vec3 xi = Vec3(i, j, k) + mRand.getVec3(); //randomized offset uniform in cell
	Vec3 vi = v.getInterpolated(xi);
	Vec3 dir = dt*vi;									//direction of movement of current particle
	Vec3 e1 = getNormalized(Vec3(dir.z, 0, -dir.x));	//perpendicular to dir
	Vec3 e2 = getNormalized(cross(e1, dir));			//perpendicular to dir and e1, so e1 and e1 create reference plane

	for (int di = 0; di < n; di++) {
		const Real r = Real(0.5) * sqrt(mRand.getReal());		//distance to cylinder axis
		const Real theta = mRand.getReal() * Real(2) * M_PI;	//azimuth
		const Real h = mRand.getReal() * norm(dt*vi);			//distance to reference plane
		Vec3 xd = xi + r*cos(theta)*e1 + r*sin(theta)*e2 + h*getNormalized(vi);
		if (!flags.is3D()) xd.z = 0;
		pts_sec.add(xd);

		v_sec[v_sec.size() - 1] = r*cos(theta)*e1 + r*sin(theta)*e2 + vi;	//init velocity of new particle
		Real temp = (KE + TA + WC) / 3;
		l_sec[l_sec.size() - 1] = ((lMax - lMin) * temp) + lMin + mRand.getReal()*0.1;	//init lifetime of new particle

		//init type of new particle
		if (neighborRatio(i, j, k) < c_s) { pts_sec[pts_sec.size() - 1].flag = ParticleBase::PSPRAY; }
		else if (neighborRatio(i, j, k) > c_b) { pts_sec[pts_sec.size() - 1].flag = ParticleBase::PBUBBLE; }
		else { pts_sec[pts_sec.size() - 1].flag = ParticleBase::PFOAM; }
	}
}
PYTHON()
void flipSampleSecondaryParticles(
	const std::string mode, const FlagGrid &flags, const MACGrid &v, BasicParticleSystem &pts_sec,
	ParticleDataImpl<Vec3> &v_sec, ParticleDataImpl<Real> &l_sec, const Real lMin, const Real lMax,
	const Grid<Real> &potTA, const Grid<Real> &potWC, const Grid<Real> &potKE, const Grid<Real> &neighborRatio,
	const Real c_s, const Real c_b, const Real k_ta, const Real k_wc, const Real dt=0, const int itype = FlagGrid::TypeFluid) {

	float timestep = dt;
	if (dt <= 0) timestep = flags.getParent()->getDt();

	if (mode == "single") {
		knFlipSampleSecondaryParticles(flags, v, pts_sec, v_sec, l_sec, lMin, lMax, potTA, potWC, potKE, neighborRatio, c_s, c_b, k_ta, k_wc, timestep, itype);
	}
	else if (mode == "multiple") {
		knFlipSampleSecondaryParticlesMoreCylinders(flags, v, pts_sec, v_sec, l_sec, lMin, lMax, potTA, potWC, potKE, neighborRatio, c_s, c_b, k_ta, k_wc, timestep, itype);
	}
	else {
		throw std::invalid_argument("Unknown mode: use \"single\" or \"multiple\" instead!");
	}
}




// evaluates cubic spline with radius h and distance l in dim dimensions
Real cubicSpline(const Real h, const Real l, const int dim) {
	const Real h2 = square(h), h3 = h2*h;
	const Real c[] = { Real(2e0 / (3e0*h)), Real(10e0 / (7e0*M_PI*h2)), Real(1e0 / (M_PI*h3)) };
	const Real q = l / h;
	if (q<1e0) return c[dim - 1] * (1e0 - 1.5*square(q) + 0.75*cubed(q));
	else if (q<2e0) return c[dim - 1] * (0.25*cubed(2e0 - q));
	return 0;
}


// updates position &pts_sec.pos and velocity &v_sec of secondary particles according to the particle type determined by the neighbor ratio with linear interpolation
KERNEL(pts)
void knFlipUpdateSecondaryParticlesLinear(
	BasicParticleSystem &pts_sec, ParticleDataImpl<Vec3> &v_sec, ParticleDataImpl<Real> &l_sec, const ParticleDataImpl<Vec3> &f_sec,
	const FlagGrid &flags, const MACGrid &v, const Grid<Real> &neighborRatio, const Vec3 gravity, const Real k_b, const Real k_d,
	const Real c_s, const Real c_b, const Real dt, const int exclude, const int antitunneling) {

	if (!pts_sec.isActive(idx) || pts_sec[idx].flag & exclude) return;
	if (!flags.isInBounds(pts_sec[idx].pos)) {
		pts_sec.kill(idx);
		return;
	}

	Vec3i gridpos = toVec3i(pts_sec[idx].pos);

	//spray particle
	if (neighborRatio(gridpos) < c_s) {
		pts_sec[idx].flag |= ParticleBase::PSPRAY;
		pts_sec[idx].flag &= ~(ParticleBase::PBUBBLE | ParticleBase::PFOAM);
		v_sec[idx] += dt * ((f_sec[idx] / 1) + gravity);	//TODO: if forces are added (e.g. fluid guiding), add parameter for mass instead of 1

		//anti tunneling for small obstacles
		for (int ct = 1; ct < antitunneling; ct++) {
			Vec3i tempPos = toVec3i(pts_sec[idx].pos + ct * (1 / Real(antitunneling)) * dt * v_sec[idx]);
			if (!flags.isInBounds(tempPos) || flags(tempPos) & FlagGrid::TypeObstacle) {
				pts_sec.kill(idx);
				return;
			}
		}
		pts_sec[idx].pos += dt * v_sec[idx];
	}

	//air bubble particle
	else if (neighborRatio(gridpos) > c_b) {
		pts_sec[idx].flag |= ParticleBase::PBUBBLE;
		pts_sec[idx].flag &= ~(ParticleBase::PSPRAY | ParticleBase::PFOAM);
		
		const Vec3 vj = (v.getInterpolated(pts_sec[idx].pos) - v_sec[idx]) / dt;
		v_sec[idx] += dt * (k_b * -gravity + k_d * vj);

		//anti tunneling for small obstacles
		for (int ct = 1; ct < antitunneling; ct++) {
			Vec3i tempPos = toVec3i(pts_sec[idx].pos + ct * (1 / Real(antitunneling)) * dt * v_sec[idx]);
			if (!flags.isInBounds(tempPos) || flags(tempPos) & FlagGrid::TypeObstacle) {
				pts_sec.kill(idx);
				return;
			}
		}
		pts_sec[idx].pos += dt * v_sec[idx];
	}

	//foam particle
	else {
		pts_sec[idx].flag |= ParticleBase::PFOAM;
		pts_sec[idx].flag &= ~(ParticleBase::PBUBBLE | ParticleBase::PSPRAY);

		const Vec3 vj = v.getInterpolated(pts_sec[idx].pos);
		//anti tunneling for small obstacles
		for (int ct = 1; ct < antitunneling; ct++) {
			Vec3i tempPos = toVec3i(pts_sec[idx].pos + ct * (1 / Real(antitunneling)) * dt * vj);
			if (!flags.isInBounds(tempPos) || flags(tempPos) & FlagGrid::TypeObstacle) {
				pts_sec.kill(idx);
				return;
			}
		}
		pts_sec[idx].pos += dt * v.getInterpolated(pts_sec[idx].pos);
	}

	//lifetime
	l_sec[idx] -= dt;
	if (l_sec[idx] <= Real(0)) {
		pts_sec.kill(idx);
	}
}
// updates position &pts_sec.pos and velocity &v_sec of secondary particles according to the particle type determined by the neighbor ratio with cubic spline interpolation
KERNEL(pts)
void knFlipUpdateSecondaryParticlesCubic(
	BasicParticleSystem &pts_sec, ParticleDataImpl<Vec3> &v_sec, ParticleDataImpl<Real> &l_sec, const ParticleDataImpl<Vec3> &f_sec,
	const FlagGrid &flags, const MACGrid &v, const Grid<Real> &neighborRatio,
	const int radius, const Vec3 gravity, const Real k_b, const Real k_d,
	const Real c_s, const Real c_b, const Real dt, const int exclude, const int antitunneling, const int itype) {

	if (!pts_sec.isActive(idx) || pts_sec[idx].flag & exclude) return;
	if (!flags.isInBounds(pts_sec[idx].pos)) {
		pts_sec.kill(idx);
		return;
	}

	Vec3i gridpos = toVec3i(pts_sec[idx].pos);
	int i = gridpos.x;
	int j = gridpos.y;
	int k = gridpos.z;

	//spray particle
	if (neighborRatio(gridpos) < c_s) {
		pts_sec[idx].flag |= ParticleBase::PSPRAY;
		pts_sec[idx].flag &= ~(ParticleBase::PBUBBLE | ParticleBase::PFOAM);
		v_sec[idx] += dt * ((f_sec[idx] / 1) + gravity);	//TODO: if forces are added (e.g. fluid guiding), add parameter for mass instead of 1

		//anti tunneling for small obstacles
		for (int ct = 1; ct < antitunneling; ct++) {
			Vec3i tempPos = toVec3i(pts_sec[idx].pos + ct * (1 / Real(antitunneling)) * dt * v_sec[idx]);
			if (!flags.isInBounds(tempPos) || flags(tempPos) & FlagGrid::TypeObstacle) {
				pts_sec.kill(idx);
				return;
			}
		}
		pts_sec[idx].pos += dt * v_sec[idx];
	}

	//air bubble particle
	else if (neighborRatio(gridpos) > c_b) {
		pts_sec[idx].flag |= ParticleBase::PBUBBLE;
		pts_sec[idx].flag &= ~(ParticleBase::PSPRAY | ParticleBase::PFOAM);
		const Vec3 &xi = pts_sec[idx].pos;
		Vec3 sumNumerator = Vec3(0, 0, 0);
		Real sumDenominator = 0;
		for (IndexInt x = i - radius; x <= i + radius; x++) {
			for (IndexInt y = j - radius; y <= j + radius; y++) {
				for (IndexInt z = k - radius; z <= k + radius; z++) {
					Vec3i xj = Vec3i(x, y, z);
					if ((x == i && y == j && z == k) || !flags.isInBounds(xj)) continue;
					if (!(flags(xj) & itype)) continue;
					const Real len_xij = norm(xi - Vec3(x, y, z));

					int dim = flags.is3D() ? 3 : 2;
					Real dist = flags.is3D() ? 1.732 : 1.414;
					Real weight = cubicSpline(radius * dist, len_xij, dim);
					sumNumerator += v.getCentered(xj) * weight;	//estimate next position by current velocity
					sumDenominator += weight;
				}
			}
		}
		const Vec3 temp = ((sumNumerator / sumDenominator) - v_sec[idx]) / dt;
		v_sec[idx] += dt * (k_b * -gravity + k_d * temp);

		//anti tunneling for small obstacles
		for (int ct = 1; ct < antitunneling; ct++) {
			Vec3i tempPos = toVec3i(pts_sec[idx].pos + ct * (1 / Real(antitunneling)) * dt * v_sec[idx]);
			if (!flags.isInBounds(tempPos) || flags(tempPos) & FlagGrid::TypeObstacle) {
				pts_sec.kill(idx);
				return;
			}
		}
		pts_sec[idx].pos += dt * v_sec[idx];
	}

	//foam particle
	else {
		pts_sec[idx].flag |= ParticleBase::PFOAM;
		pts_sec[idx].flag &= ~(ParticleBase::PBUBBLE | ParticleBase::PSPRAY);
		const Vec3 &xi = pts_sec[idx].pos;
		Vec3 sumNumerator = Vec3(0, 0, 0);
		Real sumDenominator = 0;
		for (IndexInt x = i - radius; x <= i + radius; x++) {
			for (IndexInt y = j - radius; y <= j + radius; y++) {
				for (IndexInt z = k - radius; z <= k + radius; z++) {
					Vec3i xj = Vec3i(x, y, z);
					if ((x == i && y == j && z == k) || !flags.isInBounds(xj)) continue;
					if (!(flags(xj) & itype)) continue;
					const Real len_xij = norm(xi - Vec3(x, y, z));

					int dim = flags.is3D() ? 3 : 2;
					Real dist = flags.is3D() ? 1.732 : 1.414;
					Real weight = cubicSpline(radius * dist, len_xij, dim);
					sumNumerator += v.getCentered(xj) * weight;	//estimate next position by current velocity
					sumDenominator += weight;
				}
			}
		}

		//anti tunneling for small obstacles
		for (int ct = 1; ct < antitunneling; ct++) {
			Vec3i tempPos = toVec3i(pts_sec[idx].pos + ct * (1 / Real(antitunneling)) * dt * (sumNumerator / sumDenominator));
			if (!flags.isInBounds(tempPos) || flags(tempPos) & FlagGrid::TypeObstacle) {
				pts_sec.kill(idx);
				return;
			}
		}
		pts_sec[idx].pos += dt * (sumNumerator / sumDenominator);
	}

	//lifetime
	l_sec[idx] -= dt;
	if (l_sec[idx] <= Real(0)) {
		pts_sec.kill(idx);
	}
}
PYTHON()
void flipUpdateSecondaryParticles(
	const std::string mode, BasicParticleSystem &pts_sec, ParticleDataImpl<Vec3> &v_sec, ParticleDataImpl<Real> &l_sec, const ParticleDataImpl<Vec3> &f_sec,
	FlagGrid &flags, const MACGrid &v, const Grid<Real> &neighborRatio,
	const int radius, const Vec3 gravity, const Real k_b, const Real k_d,
	const Real c_s, const Real c_b, const Real dt=0, bool scale=true, const int exclude = ParticleBase::PTRACER, const int antitunneling=0, const int itype = FlagGrid::TypeFluid) {

	float gridScale = (scale) ? flags.getParent()->getDx() : 1;
	Vec3 g = gravity / gridScale;

	float timestep = dt;
	if (dt <= 0) timestep = flags.getParent()->getDt();

	if (mode == "linear") {
		knFlipUpdateSecondaryParticlesLinear(pts_sec, v_sec, l_sec, f_sec, flags, v, neighborRatio, g, k_b, k_d, c_s, c_b, timestep, exclude, antitunneling);
	}
	else if (mode == "cubic") {
		knFlipUpdateSecondaryParticlesCubic(pts_sec, v_sec, l_sec, f_sec, flags, v, neighborRatio, radius, g, k_b, k_d, c_s, c_b, timestep, exclude, antitunneling, itype);
	}
	else {
		throw std::invalid_argument("Unknown mode: use \"linear\" or \"cubic\" instead!");
	}
	pts_sec.doCompress();
}


// removes secondary particles in &pts_sec that are inside boundaries (cells that are marked as obstacle/outflow in &flags)
KERNEL(pts)
void knFlipDeleteParticlesInObstacle(
	BasicParticleSystem &pts, const FlagGrid &flags) {

	if (!pts.isActive(idx)) return;

	const Vec3 &xi = pts[idx].pos;
	const Vec3i xidx = toVec3i(xi);
	//remove particles that completely left the bounds
	if (!flags.isInBounds(xidx)) {
		pts.kill(idx);
		return;
	}
	int gridIndex = flags.index(xidx);
	//remove particles that penetrate obstacles
	if (flags.isObstacle(gridIndex) || flags.isOutflow(gridIndex)) {
		pts.kill(idx);
	}
}
PYTHON()
void flipDeleteParticlesInObstacle(
	BasicParticleSystem &pts, const FlagGrid &flags) {

	knFlipDeleteParticlesInObstacle(pts, flags);
	pts.doCompress();
}

//helper method to debug statistical data from grid
PYTHON()
void debugGridInfo(
	const FlagGrid &flags, Grid<Real> &grid, std::string name, const int itype = FlagGrid::TypeFluid) {
	FluidSolver* s = flags.getParent();
	int countFluid = 0;
	int countLargerZero = 0;
	Real avg = 0;
	Real max = 0;
	Real sum = 0;
	Real avgLargerZero = 0;
	FOR_IJK_BND(grid, 1) {
		if (!(flags(i, j, k) & itype)) continue;
		countFluid++;
		if (grid(i, j, k) > 0) countLargerZero++;
		sum += grid(i, j, k);
		if (grid(i, j, k) > max) max = grid(i, j, k);
	}
	avg = sum / std::max(Real(countFluid), Real(1));
	avgLargerZero = sum / std::max(Real(countLargerZero), Real(1));

	debMsg("Step: " << s->mFrame  << " - Grid " << name <<
		"\n\tcountFluid \t\t" << countFluid <<
		"\n\tcountLargerZero \t" << countLargerZero <<
		"\n\tsum \t\t\t" << sum <<
		"\n\tavg \t\t\t" << avg <<
		"\n\tavgLargerZero \t\t" << avgLargerZero <<
		"\n\tmax \t\t\t" << max, 1);
}



// The following methods are helper functions to recreate the velocity and flag grid from the underlying FLIP simulation.
// They cannot simply be loaded because of the upres to a higher resolution, instead a levelset is used.
KERNEL(idx)
void knSetFlagsFromLevelset(
	FlagGrid &flags, const Grid<Real> &phi, const int exclude = FlagGrid::TypeObstacle, const int itype = FlagGrid::TypeFluid) {
	if (phi(idx) < 0 && !(flags(idx) & exclude)) flags(idx) = itype;

}
PYTHON()
void setFlagsFromLevelset(
	FlagGrid &flags, const Grid<Real> &phi, const int exclude = FlagGrid::TypeObstacle, const int itype = FlagGrid::TypeFluid) {
	knSetFlagsFromLevelset(flags, phi, exclude, itype);
}

KERNEL()
void knSetMACFromLevelset(
	MACGrid &v, const Grid<Real> &phi, const Vec3 c) {
	if (phi.getInterpolated(Vec3(i, j, k)) > 0) v(i, j, k) = c;
}
PYTHON()
void setMACFromLevelset(
	MACGrid &v, const Grid<Real> &phi, const Vec3 c) {
	knSetMACFromLevelset(v, phi, c);
}



//----------------------------------------------------------------------------------------------------------------------------------------------------
// END Secondary Particles for FLIP
//----------------------------------------------------------------------------------------------------------------------------------------------------
#pragma endregion



#pragma region Legacy Methods (still useful for debugging)
//-----------------------------------------------------------------------------------------------------------------------------------	-----------------
// Legacy Methods (still useful for debugging)
//----------------------------------------------------------------------------------------------------------------------------------------------------

// LEGACY METHOD! Use flipComputeSecondaryParticlePotentials instead!
// computes trapped air potential for all fluid cells in &flags and saves it in &pot
KERNEL(bnd = 1)
void knFlipComputePotentialTrappedAir(
	Grid<Real> &pot, const FlagGrid &flags, const MACGrid &v,
	const int radius, const Real tauMin, const Real tauMax,
	const Real scaleFromManta, const int itype = FlagGrid::TypeFluid,
	const int jtype = FlagGrid::TypeFluid) {

	if (!(flags(i, j, k) & itype)) return;

	const Vec3 &xi = scaleFromManta * Vec3(i, j, k);	//scale to unit cube
	const Vec3 &vi = scaleFromManta * v.getCentered(i, j, k);
	Real vdiff = 0;
	for (IndexInt x = i - radius; x <= i + radius; x++) {
		for (IndexInt y = j - radius; y <= j + radius; y++) {
			for (IndexInt z = k - radius; z <= k + radius; z++) {
				if ((x == i && y == j && z == k) || !(flags(x, y, z) & jtype)) continue;

				const Vec3 &xj = scaleFromManta * Vec3(x, y, z); //scale to unit cube
				const Vec3 &vj = scaleFromManta * v.getCentered(x, y, z);
				const Vec3 xij = xi - xj;
				const Vec3 vij = vi - vj;
				Real h = !pot.is3D() ? 1.414*radius : 1.732*radius; //estimate sqrt(2)*radius resp. sqrt(3)*radius for h, due to squared resp. cubic neighbor area
				vdiff += norm(vij) * (1 - dot(getNormalized(vij), getNormalized(xij))) * (1 - norm(xij) / h);
			}
		}
	}
	pot(i, j, k) = (std::min(vdiff, tauMax) - std::min(vdiff, tauMin)) / (tauMax - tauMin);
}
PYTHON()
void flipComputePotentialTrappedAir(
	Grid<Real> &pot, const FlagGrid &flags, const MACGrid &v,
	const int radius, const Real tauMin, const Real tauMax,
	const Real scaleFromManta, const int itype = FlagGrid::TypeFluid,
	const int jtype = FlagGrid::TypeFluid) {
	pot.clear();
	knFlipComputePotentialTrappedAir(pot, flags, v, radius, tauMin, tauMax, scaleFromManta, itype, jtype);
}


// LEGACY METHOD! Use flipComputeSecondaryParticlePotentials instead!
// computes kinetic energy potential for all fluid cells in &flags and saves it in &pot
KERNEL()
void knFlipComputePotentialKineticEnergy(
	Grid<Real> &pot, const FlagGrid &flags, const MACGrid &v,
	const Real tauMin, const Real tauMax, const Real scaleFromManta,
	const int itype = FlagGrid::TypeFluid) {

	if (!(flags(i, j, k) & itype)) return;

	const Vec3 &vi = scaleFromManta * v.getCentered(i, j, k); //scale to unit cube
	Real ek = Real(0.5) * 125 * normSquare(vi);	//use arbitrary constant for mass, potential adjusts with thresholds anyways
	pot(i, j, k) = (std::min(ek, tauMax) - std::min(ek, tauMin)) / (tauMax - tauMin);
}
PYTHON()
void flipComputePotentialKineticEnergy(
	Grid<Real> &pot, const FlagGrid &flags, const MACGrid &v,
	const Real tauMin, const Real tauMax, const Real scaleFromManta,
	const int itype = FlagGrid::TypeFluid) {
	pot.clear();
	knFlipComputePotentialKineticEnergy(pot, flags, v, tauMin, tauMax, scaleFromManta, itype);
}


// LEGACY METHOD! Use flipComputeSecondaryParticlePotentials instead!
// computes wave crest potential for all fluid cells in &flags and saves it in &pot
KERNEL(bnd = 1)
void knFlipComputePotentialWaveCrest(
	Grid<Real> &pot, const FlagGrid &flags, const MACGrid &v, const int radius,
	Grid<Vec3> &normal, const Real tauMin, const Real tauMax, const Real scaleFromManta,
	const int itype = FlagGrid::TypeFluid, const int jtype = FlagGrid::TypeFluid) {

	if (!(flags(i, j, k) & itype)) return;

	const Vec3 &xi = scaleFromManta * Vec3(i, j, k);	//scale to unit cube
	const Vec3 &vi = scaleFromManta * v.getCentered(i, j, k);
	const Vec3 &ni = normal(i, j, k);
	Real kappa = 0;
	for (IndexInt x = i - radius; x <= i + radius; x++) {
		for (IndexInt y = j - radius; y <= j + radius; y++) {
			for (IndexInt z = k - radius; z <= k + radius; z++) {
				if ((x == i && y == j && z == k) || !(flags(x, y, z) & jtype)) continue;
				const Vec3 &xj = scaleFromManta * Vec3(x, y, z); //scale to unit cube
				const Vec3 &nj = normal(x, y, z);
				const Vec3 xij = xi - xj;
				if (dot(getNormalized(xij), ni) < 0) {	//identifies wave crests
					Real h = !pot.is3D() ? 1.414*radius : 1.732*radius; //estimate sqrt(2)*radius resp. sqrt(3)*radius for h, due to squared resp. cubic neighbor area
					kappa += (1 - dot(ni, nj)) * (1 - norm(xij) / h);
				}
			}
		}
	}

	if (dot(getNormalized(vi), ni) >= 0.6) {	//avoid to mark boarders of the scene as wave crest
		pot(i, j, k) = (std::min(kappa, tauMax) - std::min(kappa, tauMin)) / (tauMax - tauMin);
	}
	else {
		pot(i, j, k) = Real(0);
	}
}
PYTHON()
void flipComputePotentialWaveCrest(
	Grid<Real> &pot, const FlagGrid &flags, const MACGrid &v, const int radius,
	Grid<Vec3> &normal, const Real tauMin, const Real tauMax, const Real scaleFromManta,
	const int itype = FlagGrid::TypeFluid, const int jtype = FlagGrid::TypeFluid) {

	pot.clear();
	knFlipComputePotentialWaveCrest(pot, flags, v, radius, normal, tauMin, tauMax, scaleFromManta, itype, jtype);
}

// LEGACY METHOD! Use flipComputeSecondaryParticlePotentials instead!
// computes normal grid &normal as gradient of levelset &phi and normalizes it
KERNEL(idx)
void knFlipComputeSurfaceNormals(Grid<Vec3>& normal, const Grid<Real>& phi) {
	normal[idx] = getNormalized(normal[idx]);
}
PYTHON()
void flipComputeSurfaceNormals(Grid<Vec3>& normal, const Grid<Real>& phi) {
	GradientOp(normal, phi);
	knFlipComputeSurfaceNormals(normal, phi);
}

// LEGACY METHOD! Use flipComputeSecondaryParticlePotentials instead!
// computes the neighbor ratio for every fluid cell in &flags as the number of fluid neighbors over the maximum possible number of fluid neighbors
KERNEL(bnd = 1)
void knFlipUpdateNeighborRatio(
	const FlagGrid &flags, Grid<Real> &neighborRatio, const int radius,
	const int itype = FlagGrid::TypeFluid, const int jtype = FlagGrid::TypeObstacle) {

	if (!(flags(i, j, k) & itype)) return;

	int countFluid = 0;
	int countMaxFluid = 0;
	for (IndexInt x = i - radius; x <= i + radius; x++) {
		for (IndexInt y = j - radius; y <= j + radius; y++) {
			for (IndexInt z = k - radius; z <= k + radius; z++) {
				if ((x == i && y == j && z == k) || (flags(x, y, z) & jtype)) continue;
				if (flags(x, y, z) & itype) {
					countFluid++;
					countMaxFluid++;
				}
				else {
					countMaxFluid++;
				}
			}
		}
	}
	neighborRatio(i, j, k) = float(countFluid) / float(countMaxFluid);
}
PYTHON()
void flipUpdateNeighborRatio(
	const FlagGrid &flags, Grid<Real> &neighborRatio, const int radius,
	const int itype = FlagGrid::TypeFluid, const int jtype = FlagGrid::TypeObstacle) {

	neighborRatio.clear();
	knFlipUpdateNeighborRatio(flags, neighborRatio, radius, itype, jtype);
}


//----------------------------------------------------------------------------------------------------------------------------------------------------
// Legacy Methods (still useful for debugging)
//----------------------------------------------------------------------------------------------------------------------------------------------------
#pragma endregion



} //namespace