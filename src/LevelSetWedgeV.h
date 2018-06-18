#pragma once

#ifndef OPENVDB_TOOLS_LEVELSETWEDGEV_HAS_BEEN_INCLUDED
#define OPENVDB_TOOLS_LEVELSETWEDGEV_HAS_BEEN_INCLUDED

#include <openvdb/Grid.h>
#include <openvdb/Types.h>
#include <openvdb/math/Math.h>
#include <openvdb/util/NullInterrupter.h>
#include <boost/utility.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <openvdb/tools/SignedFloodFill.h>
#include <openvdb/math/Proximity.h>

namespace openvdb {
	OPENVDB_USE_VERSION_NAMESPACE
	namespace OPENVDB_VERSION_NAME {
		namespace tools {

			/// @brief Return a grid of type @c GridType containing a narrow-band level set
			/// representation of a sphere.
			///
			/// @param radius       radius of the sphere in world units
			/// @param center       center of the sphere in world units
			/// @param voxelSize    voxel size in world units
			/// @param halfWidth    half the width of the narrow band, in voxel units
			/// @param interrupt    a pointer adhering to the util::NullInterrupter interface
			///
			/// @note @c GridType::ValueType must be a floating-point scalar.
			/// @note The leapfrog algorithm employed in this method is best suited
			/// for a single large sphere.  For multiple small spheres consider
			/// using the faster algorithm in ParticlesToLevelSet.h
			template<typename GridType, typename InterruptT>
			typename GridType::Ptr
				createLevelSetWedgeV(const openvdb::Vec3f& p1, const openvdb::Vec3f& p2, const openvdb::Vec3f& p3, float t1, float t2, float t3, float voxelSize,
					float halfWidth = float(LEVEL_SET_HALF_WIDTH), InterruptT* interrupt = NULL);

			/// @brief Return a grid of type @c GridType containing a narrow-band level set
			/// representation of a Capsule.
			///
			/// @param radius       radius of the Capsule in world units
			/// @param start       center of the sphere in world units
			/// @param voxelSize    voxel size in world units
			/// @param halfWidth    half the width of the narrow band, in voxel units
			///
			/// @note @c GridType::ValueType must be a floating-point scalar.
			/// @note The leapfrog algorithm employed in this method is best suited
			/// for a single large sphere.  For multiple small spheres consider
			/// using the faster algorithm in ParticlesToLevelSet.h
			template<typename GridType>
			typename GridType::Ptr
				createLevelSetWedgeV(const openvdb::Vec3f& p1, const openvdb::Vec3f& p2, const openvdb::Vec3f& p3, float t1, float t2, float t3, float voxelSize,
					float halfWidth = float(LEVEL_SET_HALF_WIDTH))
			{
				return createLevelSetWedgeV<GridType, util::NullInterrupter>(radius, start, end, voxelSize, halfWidth);
			}


			////////////////////////////////////////


			/// @brief Generates a signed distance field (or narrow band level
			/// set) to a single sphere.
			///
			/// @note The leapfrog algorithm employed in this class is best
			/// suited for a single large sphere. For multiple small spheres consider
			/// using the faster algorithm in tools/ParticlesToLevelSet.h
			template<typename GridT, typename InterruptT = util::NullInterrupter>
			class LevelSetWedgeV
			{
			public:
				typedef typename GridT::ValueType   ValueT;
				typedef typename math::Vec3<ValueT> Vec3T;
				BOOST_STATIC_ASSERT(boost::is_floating_point<ValueT>::value);

				/// @brief Constructor
				///
				/// @param radius radius of the sphere in world units
				/// @param center center of the sphere in world units
				/// @param interrupt pointer to optional interrupter. Use template
				/// argument util::NullInterrupter if no interruption is desired.
				///
				/// @note If the radius of the sphere is smaller than
				/// 1.5*voxelSize, i.e. the sphere is smaller than the Nyquist
				/// frequency of the grid, it is ignored!
				LevelSetWedgeV(const Vec3T &p1, const Vec3T &p2, const Vec3T &p3, float t1, float t2, float t3, InterruptT* interrupt = NULL)
					: mT1(t1), mT2(t2), mT3(t3), mP1(p1), mP2(p2), mP3(p3), mInterrupt(interrupt)
				{
					//if (mThickness <= 0) OPENVDB_THROW(ValueError, "radius must be positive");
					//check for negative thickness
				}

				/// @return a narrow-band level set of the sphere
				///
				/// @param voxelSize  Size of voxels in world units
				/// @param halfWidth  Half-width of narrow-band in voxel units
				typename GridT::Ptr getLevelSet(ValueT voxelSize, ValueT halfWidth)
				{
					mGrid = createLevelSet<GridT>(voxelSize, halfWidth);
					this->rasterWedge(voxelSize, halfWidth);
					mGrid->setGridClass(GRID_LEVEL_SET);
					return mGrid;
				}

				void rasterWedge(ValueT dx, ValueT w)
				{
					if (!(dx>0.0f)) OPENVDB_THROW(ValueError, "voxel size must be positive");
					if (!(w>1)) OPENVDB_THROW(ValueError, "half-width must be larger than one");

					// Define radius of sphere and narrow-band in voxel units
					const ValueT r1 = mT1 / dx;
					const ValueT r2 = mT2 / dx;
					const ValueT r3 = mT3 / dx;
					ValueT rMax = max(r1, r2, r3);

					// Radius below the Nyquist frequency
					if (r0 < 1.5f)  return;

					// Define dir
					const Vec3T a(mP1[0] / dx, mP1[1] / dx, mP1[2] / dx);
					const Vec3T b(mP2[0] / dx, mP2[1] / dx, mP2[2] / dx);
					const Vec3T c(mP3[0] / dx, mP3[1] / dx, mP3[2] / dx);

					Vec3T dir = a - b;
					const ValueT len = dir.length();
					dir /= len;
					// Define index coordinates and their respective bounds
					openvdb::Coord ijk;
					int &i = ijk[0], &j = ijk[1], &k = ijk[2], m = 1;
					const int imin = min(min(math::Floor(a[0] - rmax), math::Floor(b[0] - rmax)), math::Floor(c[0] - rmax));
					const int imax = max(max(math::Ceil(a[0] + rmax), math::Ceil(b[0] + rmax)), math::Ceil(c[0] + rmax));
					const int jmin = min(min(math::Floor(a[1] - rmax), math::Floor(b[1] - rmax)), math::Floor(c[1] - rmax));
					const int jmax = max(max(math::Ceil(a[1] + rmax), math::Ceil(b[1] + rmax)), math::Ceil(c[1] + rmax));
					const int kmin = min(min(math::Floor(a[2] - rmax), math::Floor(b[2] - rmax)), math::Floor(c[2] - rmax));
					const int kmax = max(max(math::Ceil(a[2] + rmax), math::Ceil(b[2] + rmax)), math::Ceil(c[2] + rmax));

					// Allocate a ValueAccessor for accelerated random access
					typename GridT::Accessor accessor = mGrid->getAccessor();

					ValueT pval;
					const ValueT inside = -mGrid->background();
					Vec3d uvw;
					Vec3d closest, vijk;
					if (mInterrupt) mInterrupt->start("Generating level set of wedge");
					// Compute signed distances to sphere using leapfrogging in k
					for (i = imin; i <= imax; ++i) {
						if (util::wasInterrupted(mInterrupt)) return;
						//const float x2 = math::Pow2(i - c[0]);
						vijk[0] = i;

						for (j = jmin; j <= jmax; ++j) {
							vijk[1] = j;
							//const float x2y2 = math::Pow2(j - c[1]) + x2;
							for (k = kmin; k <= kmax; k += m) {
								vijk[2]= k;
								m = 1;
								closest = (vijk - math::closestPointOnTriangleToPoint(a, b, c, vijk, uvw));
								ValueT v = closest.length()-r0,
									d = math::Abs(v);

								if (v > w || (!accessor.probeValue(ijk, pval) && pval < 0)) {
									m += math::Floor(d - w);//leapfrog
								}
								else if (v < -w) {
									accessor.setValueOff(ijk, inside);
								}
								else {
									if (dx*v < pval) accessor.setValue(ijk, dx*v);
								}
								//if (d < w) { // inside narrow band
								//	accessor.setValue(ijk, dx*v);// distance in world units
								//}
								//else {// outside narrow band
								//	m += math::Floor(d - w);// leapfrog
								//}
							}//end leapfrog over k
						}//end loop over j
					}//end loop over i

					 // Define consistent signed distances outside the narrow-band
					 //tools::signedFloodFill(mGrid->tree());

					if (mInterrupt) mInterrupt->end();
				}

				const ValueT        mT1, mT2, mT3;
				const Vec3T         mP1, mP2, mP3;
				InterruptT*         mInterrupt;
				typename GridT::Ptr mGrid;
			};// LevelSetSphere

			  ////////////////////////////////////////


			  /// @brief Generates a signed distance field (or narrow band level
			  /// set) to a single sphere.
			  ///
			  /// @note The leapfrog algorithm employed in this class is best
			  /// suited for a single large sphere. For multiple small spheres consider
			  /// using the faster algorithm in tools/ParticlesToLevelSet.h
			template<typename GridT, typename InterruptT = util::NullInterrupter>
			class LevelSetWedgeSpecial
			{
			public:
				typedef typename GridT::ValueType   ValueT;
				typedef typename math::Vec3<ValueT> Vec3T;
				BOOST_STATIC_ASSERT(boost::is_floating_point<ValueT>::value);

				/// @brief Constructor
				///
				/// @param radius radius of the sphere in world units
				/// @param center center of the sphere in world units
				/// @param interrupt pointer to optional interrupter. Use template
				/// argument util::NullInterrupter if no interruption is desired.
				///
				/// @note If the radius of the sphere is smaller than
				/// 1.5*voxelSize, i.e. the sphere is smaller than the Nyquist
				/// frequency of the grid, it is ignored!
				LevelSetWedgeSpecial(const Vec3T &p1, const Vec3T &p2, const Vec3T &p3, float thickness, InterruptT* interrupt = NULL)
					: mThickness(thickness), mP1(p1), mP2(p2), mP3(p3), mInterrupt(interrupt)
				{
					if (mThickness <= 0) OPENVDB_THROW(ValueError, "radius must be positive");
				}

				/// @return a narrow-band level set of the sphere
				///
				/// @param voxelSize  Size of voxels in world units
				/// @param halfWidth  Half-width of narrow-band in voxel units
				typename GridT::Ptr getLevelSet(ValueT voxelSize, ValueT halfWidth)
				{
					mGrid = createLevelSet<GridT>(voxelSize, halfWidth);
					this->rasterWedge(voxelSize, halfWidth);
					mGrid->setGridClass(GRID_LEVEL_SET);
					return mGrid;
				}

				void rasterWedge(ValueT dx, ValueT w)
				{
					if (!(dx>0.0f)) OPENVDB_THROW(ValueError, "voxel size must be positive");
					if (!(w>1)) OPENVDB_THROW(ValueError, "half-width must be larger than one");

					// Define radius of sphere and narrow-band in voxel units
					const ValueT r0 = mThickness / dx, rmax = r0 + w;

					// Radius below the Nyquist frequency
					if (r0 < 1.5f)  return;

					// Define dir
					const Vec3T a(mP1[0] / dx, mP1[1] / dx, mP1[2] / dx);
					const Vec3T b(mP2[0] / dx, mP2[1] / dx, mP2[2] / dx);
					const Vec3T c(mP3[0] / dx, mP3[1] / dx, mP3[2] / dx);

					Vec3d dir = a - b;
					Vec3d dir2 = c - b;
					Vec3d normD = dir.cross(dir2);
					normD.normalize();
					const ValueT len = dir.length();
					dir /= len;
					// Define index coordinates and their respective bounds
					openvdb::Coord ijk;
					int &i = ijk[0], &j = ijk[1], &k = ijk[2], m = 1;
					const int imin = min(min(math::Floor(a[0] - rmax), math::Floor(b[0] - rmax)), math::Floor(c[0] - rmax));
					const int imax = max(max(math::Ceil(a[0] + rmax), math::Ceil(b[0] + rmax)), math::Ceil(c[0] + rmax));
					const int jmin = min(min(math::Floor(a[1] - rmax), math::Floor(b[1] - rmax)), math::Floor(c[1] - rmax));
					const int jmax = max(max(math::Ceil(a[1] + rmax), math::Ceil(b[1] + rmax)), math::Ceil(c[1] + rmax));
					const int kmin = min(min(math::Floor(a[2] - rmax), math::Floor(b[2] - rmax)), math::Floor(c[2] - rmax));
					const int kmax = max(max(math::Ceil(a[2] + rmax), math::Ceil(b[2] + rmax)), math::Ceil(c[2] + rmax));

					// Allocate a ValueAccessor for accelerated random access
					typename GridT::Accessor accessor = mGrid->getAccessor();

					ValueT pval;
					const ValueT inside = -mGrid->background();
					Vec3d uvw;
					Vec3d closest, vijk;
					if (mInterrupt) mInterrupt->start("Generating level set of wedge");
					// Compute signed distances to sphere using leapfrogging in k
					for (i = imin; i <= imax; ++i) {
						if (util::wasInterrupted(mInterrupt)) return;
						//const float x2 = math::Pow2(i - c[0]);
						vijk[0] = i;

						for (j = jmin; j <= jmax; ++j) {
							vijk[1] = j;
							//const float x2y2 = math::Pow2(j - c[1]) + x2;
							for (k = kmin; k <= kmax; k += m) {
								vijk[2] = k;
								m = 1;
								closest = (vijk - math::closestPointOnTriangleToPoint(a, b, c, vijk, uvw));
								if (uvw[0] > 0 && uvw[1] > 0 && uvw[2] > 0) {
									ValueT dist = closest.length();
									ValueT v = dist - r0,
										d = math::Abs(v);

									if (v > w) {// || (!accessor.probeValue(ijk, pval) && pval < 0)) {
										m += math::Floor(d - w);//leapfrog
									}
									else {
										ValueT theta = theta1*uvw[0] + theta2*uvw[1] + theta3*uvw[2];
										ValueT beta = beta1*uvw[0] + beta2*uvw[1] + beta3*uvw[2];
										ValueT lLog = 0.;
										if (uvw.z() <= uvw.x() && uvw.z() <= uvw.y()) lLog = (M_PI / 3.) * (1. + (uvw.y() - uvw.x()) / (1. - 3.*uvw.z()));
										else if (uvw.x() <= uvw.y() && uvw.x() <= uvw.z()) lLog = (M_PI / 3.) * (3. + (uvw.z() - uvw.y()) / (1. - 3.*uvw.x()));
										else                                lLog = (M_PI / 3.) * (5. + (uvw.x() - uvw.z()) / (1. - 3.*uvw.y()));

										//ValueT nDist = math::Sign(closest.dot(normD))*dist;
										theta += lLog*nu;
										beta += lLog*nv;
										/*
										if (cos(theta * 2) > .8 || cos(beta * 2) > .8) {
											v = dist - r0*.3;
											if (!accessor.probeValue(ijk, pval) && pval < 0) {

											}
											else if (v <= w) {
												if (dx*v < pval) accessor.setValue(ijk, dx*v);
											}
										}
										*/
										//get closest pt on ribbon
										ValueT cTheta = theta;
										float cosT = cos(theta);
										float cosB = cos(beta);
										//if (cos(theta * 2) > .75) {
											ValueT h = r0*0.6*cosT*cosB;
											closest -= normD*h;
											dist = closest.length();

											ValueT tMod = math::Clamp((cos(theta * 2) - .5) / .3, -1.0, 1.0);
											v = dist - r0*tMod*0.4;
											if (!accessor.probeValue(ijk, pval) && pval < 0) {

											}if (v <= w) {
												if (dx*v < pval) accessor.setValue(ijk, dx*v);
											}

											closest += 2*normD*h;
											dist = closest.length();

											tMod = math::Clamp((cos(beta * 2) - .5) / .3, -1.0, 1.0);
											v = dist - r0*tMod*0.4;
											if (!accessor.probeValue(ijk, pval) && pval < 0) {

											}if (v <= w) {
												if (dx*v < pval) accessor.setValue(ijk, dx*v);
											}
										//}
										/*
										if (cos(beta * 2) > .75) {
											ValueT h = -r0*.6*cos(theta)*cos(beta);
											v = math::Abs(nDist - h) - r0*.3;
											if (!accessor.probeValue(ijk, pval) && pval < 0) {

											}if (v <= w) {
												if (dx*v < pval) accessor.setValue(ijk, dx*v);
											}
										}
										else {

										}
										*/
										//ValueT h = max(r0*cos(theta),ValueT(0.));
										//v = dist - h;
										//if (v > w || (!accessor.probeValue(ijk, pval) && pval < 0)) {
										//} else if (v < -w) {
										//	accessor.setValueOff(ijk, inside);
										//}
										//else {
										//	if (dx*v < pval) accessor.setValue(ijk, dx*v);
										//}
									}
									//if (d < w) { // inside narrow band
									//	accessor.setValue(ijk, dx*v);// distance in world units
									//}
									//else {// outside narrow band
									//	m += math::Floor(d - w);// leapfrog
									//}
								}
							}//end leapfrog over k
						}//end loop over j
					}//end loop over i

					 // Define consistent signed distances outside the narrow-band
					 //tools::signedFloodFill(mGrid->tree());

					if (mInterrupt) mInterrupt->end();
				}

				const ValueT        mThickness;
				ValueT				nu, nv, theta1,theta2,theta3, beta1, beta2, beta3;
				const Vec3T         mP1, mP2, mP3;
				InterruptT*         mInterrupt;
				typename GridT::Ptr mGrid;
			};// LevelSetSphere
			  ////////////////////////////////////////


			template<typename GridType, typename InterruptT>
			typename GridType::Ptr
				createLevelSetWedge(const openvdb::Vec3f& p1, const openvdb::Vec3f& p2, const openvdb::Vec3f& p3, float thickness, float voxelSize,
					float halfWidth, InterruptT* interrupt)
			{
				// GridType::ValueType is required to be a floating-point scalar.
				BOOST_STATIC_ASSERT(boost::is_floating_point<typename GridType::ValueType>::value);

				typedef typename GridType::ValueType ValueT;
				LevelSetWedge<GridType, InterruptT> factory(ValueT(thickness), p1,p2,p3, interrupt);
				return factory.getLevelSet(ValueT(voxelSize), ValueT(halfWidth));
			}

		} // namespace tools
	} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

#endif // OPENVDB_TOOLS_LEVELSETWEDGE_HAS_BEEN_INCLUDED
