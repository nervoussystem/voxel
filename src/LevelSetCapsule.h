///////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2012-2015 DreamWorks Animation LLC
//
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
//
// Redistributions of source code must retain the above copyright
// and license notice and the following restrictions and disclaimer.
//
// *     Neither the name of DreamWorks Animation nor the names of
// its contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// IN NO EVENT SHALL THE COPYRIGHT HOLDERS' AND CONTRIBUTORS' AGGREGATE
// LIABILITY FOR ALL CLAIMS REGARDLESS OF THEIR BASIS EXCEED US$250.00.
//
///////////////////////////////////////////////////////////////////////////
///
/// @file LevelSetSphere.h
///
/// @brief Generate a narrow-band level set of sphere.
///
/// @note By definition a level set has a fixed narrow band width
/// (the half width is defined by LEVEL_SET_HALF_WIDTH in Types.h),
/// whereas an SDF can have a variable narrow band width.

#ifndef OPENVDB_TOOLS_LEVELSETCAPSULE_HAS_BEEN_INCLUDED
#define OPENVDB_TOOLS_LEVELSETCAPSULE_HAS_BEEN_INCLUDED

#include <openvdb/Grid.h>
#include <openvdb/Types.h>
#include <openvdb/math/Math.h>
#include <openvdb/util/NullInterrupter.h>
#include <boost/utility.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <openvdb/tools/SignedFloodFill.h>

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
				createLevelSetCapsule(float radius, const openvdb::Vec3f& start, const openvdb::Vec3f& end, float voxelSize,
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
				createLevelSetCapsule(float radius, const openvdb::Vec3f& start, const openvdb::Vec3f& end, float voxelSize,
					float halfWidth = float(LEVEL_SET_HALF_WIDTH))
			{
				return createLevelSetCapsule<GridType, util::NullInterrupter>(radius, start, end, voxelSize, halfWidth);
			}


			////////////////////////////////////////


			/// @brief Generates a signed distance field (or narrow band level
			/// set) to a single sphere.
			///
			/// @note The leapfrog algorithm employed in this class is best
			/// suited for a single large sphere. For multiple small spheres consider
			/// using the faster algorithm in tools/ParticlesToLevelSet.h
			template<typename GridT, typename InterruptT = util::NullInterrupter>
			class LevelSetCapsule
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
				LevelSetCapsule(ValueT radius, const Vec3T &start, const Vec3T &end, InterruptT* interrupt = NULL)
					: mRadius(radius), mStart(start), mEnd(end), mInterrupt(interrupt)
				{
					if (mRadius <= 0) OPENVDB_THROW(ValueError, "radius must be positive");
				}

				/// @return a narrow-band level set of the sphere
				///
				/// @param voxelSize  Size of voxels in world units
				/// @param halfWidth  Half-width of narrow-band in voxel units
				typename GridT::Ptr getLevelSet(ValueT voxelSize, ValueT halfWidth)
				{
					mGrid = createLevelSet<GridT>(voxelSize, halfWidth);
					this->rasterCapsule(voxelSize, halfWidth);
					mGrid->setGridClass(GRID_LEVEL_SET);
					return mGrid;
				}

				void rasterCapsule(ValueT dx, ValueT w)
				{
					if (!(dx>0.0f)) OPENVDB_THROW(ValueError, "voxel size must be positive");
					if (!(w>1)) OPENVDB_THROW(ValueError, "half-width must be larger than one");

					// Define radius of sphere and narrow-band in voxel units
					const ValueT r0 = mRadius / dx, rmax = r0 + w;

					// Radius below the Nyquist frequency
					if (r0 < 1.5f)  return;

					// Define dir
					const Vec3T s(mStart[0] / dx, mStart[1] / dx, mStart[2] / dx);
					const Vec3T e(mEnd[0] / dx, mEnd[1] / dx, mEnd[2] / dx);

					Vec3T dir = e - s;
					const ValueT len = dir.length();
					dir /= len;
					// Define index coordinates and their respective bounds
					openvdb::Coord ijk;
					int &i = ijk[0], &j = ijk[1], &k = ijk[2], m = 1;
					const int imin = min(math::Floor(s[0] - rmax), math::Floor(e[0] - rmax));
					const int imax = max(math::Ceil(s[0] + rmax), math::Ceil(e[0] + rmax));
					const int jmin = min(math::Floor(s[1] - rmax), math::Floor(e[1] - rmax));
					const int jmax = max(math::Ceil(s[1] + rmax), math::Ceil(e[1] + rmax));
					const int kmin = min(math::Floor(s[2] - rmax), math::Floor(e[2] - rmax));
					const int kmax = max(math::Ceil(s[2] + rmax), math::Ceil(e[2] + rmax));

					// Allocate a ValueAccessor for accelerated random access
					typename GridT::Accessor accessor = mGrid->getAccessor();

					ValueT pval;
					const ValueT inside = -mGrid->background();

					if (mInterrupt) mInterrupt->start("Generating level set of capsule");
					// Compute signed distances to sphere using leapfrogging in k
					for (i = imin; i <= imax; ++i) {
						if (util::wasInterrupted(mInterrupt)) return;
						//const float x2 = math::Pow2(i - c[0]);
						const float x = i - s[0];
						const float xdot = x*dir[0];

						for (j = jmin; j <= jmax; ++j) {
							//const float x2y2 = math::Pow2(j - c[1]) + x2;
							const float y = j - s[1];
							const float xydot = y*dir[1]+xdot;
							for (k = kmin; k <= kmax; k += m) {
								m = 1;
								const float z = k - s[2];
								float dot = z*dir[2]+xydot;
								dot = min(max(dot, 0.0f), len);
								const Vec3T projPt = dot*dir - Vec3T(x,y,z);
								/// Distance in voxel units to capsule
								//const float v = math::Sqrt(x2y2 + math::Pow2(k - c[2])) - r0,
								const float v = projPt.length()-r0,
									d = math::Abs(v);
								if (v > w || (!accessor.probeValue(ijk, pval) && pval < 0)) {

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

				const ValueT        mRadius;
				const Vec3T         mStart;
				const Vec3T         mEnd;
				InterruptT*         mInterrupt;
				typename GridT::Ptr mGrid;
			};// LevelSetSphere


			  ////////////////////////////////////////


			template<typename GridType, typename InterruptT>
			typename GridType::Ptr
				createLevelSetCapsule(float radius, const openvdb::Vec3f& start, const openvdb::Vec3f& end, float voxelSize,
					float halfWidth, InterruptT* interrupt)
			{
				// GridType::ValueType is required to be a floating-point scalar.
				BOOST_STATIC_ASSERT(boost::is_floating_point<typename GridType::ValueType>::value);

				typedef typename GridType::ValueType ValueT;
				LevelSetCapsule<GridType, InterruptT> factory(ValueT(radius), start, end, interrupt);
				return factory.getLevelSet(ValueT(voxelSize), ValueT(halfWidth));
			}

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
				createLevelSetCapsuleV(float radius1, float radius2, const openvdb::Vec3f& start, const openvdb::Vec3f& end, float voxelSize,
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
				createLevelSetCapsuleV(float radius1, float radius2, const openvdb::Vec3f& start, const openvdb::Vec3f& end, float voxelSize,
					float halfWidth = float(LEVEL_SET_HALF_WIDTH))
			{
				return createLevelSetCapsuleV<GridType, util::NullInterrupter>(radius1, radius2, start, end, voxelSize, halfWidth);
			}

			////////////////////////////////////////


			/// @brief Generates a signed distance field (or narrow band level
			/// set) to a single sphere.
			///
			/// @note The leapfrog algorithm employed in this class is best
			/// suited for a single large sphere. For multiple small spheres consider
			/// using the faster algorithm in tools/ParticlesToLevelSet.h
			template<typename GridT, typename InterruptT = util::NullInterrupter>
			class LevelSetCapsuleV
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
				LevelSetCapsuleV(ValueT radius1, ValueT radius2, const Vec3T &start, const Vec3T &end, InterruptT* interrupt = NULL)
					: mRadius1(radius1), mRadius2(radius2), mStart(start), mEnd(end), mInterrupt(interrupt)
				{
					if (mRadius1 <= 0 || mRadius2 <= 0) OPENVDB_THROW(ValueError, "radius must be positive");
				}

				/// @return a narrow-band level set of the sphere
				///
				/// @param voxelSize  Size of voxels in world units
				/// @param halfWidth  Half-width of narrow-band in voxel units
				typename GridT::Ptr getLevelSet(ValueT voxelSize, ValueT halfWidth)
				{
					mGrid = createLevelSet<GridT>(voxelSize, halfWidth);
					this->rasterCapsule(voxelSize, halfWidth);
					mGrid->setGridClass(GRID_LEVEL_SET);
					return mGrid;
				}

				void rasterCapsule(ValueT dx, ValueT w)
				{
					if (!(dx>0.0f)) OPENVDB_THROW(ValueError, "voxel size must be positive");
					if (!(w>1)) OPENVDB_THROW(ValueError, "half-width must be larger than one");

					// Define radius of sphere and narrow-band in voxel units
					const ValueT r1 = mRadius1 / dx;
					const ValueT r2 = mRadius2 / dx;

					ValueT rmax = max(r1, r2) + w;
					ValueT r0 = min(r1, r2);
					// Radius below the Nyquist frequency
					if (r0 < 1.5f)  return;

					// Define dir
					const Vec3T s(mStart[0] / dx, mStart[1] / dx, mStart[2] / dx);
					const Vec3T e(mEnd[0] / dx, mEnd[1] / dx, mEnd[2] / dx);

					Vec3T dir = e - s;
					const ValueT len = dir.length();
					dir /= len;
					// Define index coordinates and their respective bounds
					openvdb::Coord ijk;
					int &i = ijk[0], &j = ijk[1], &k = ijk[2], m = 1;
					const int imin = min(math::Floor(s[0] - rmax), math::Floor(e[0] - rmax));
					const int imax = max(math::Ceil(s[0] + rmax), math::Ceil(e[0] + rmax));
					const int jmin = min(math::Floor(s[1] - rmax), math::Floor(e[1] - rmax));
					const int jmax = max(math::Ceil(s[1] + rmax), math::Ceil(e[1] + rmax));
					const int kmin = min(math::Floor(s[2] - rmax), math::Floor(e[2] - rmax));
					const int kmax = max(math::Ceil(s[2] + rmax), math::Ceil(e[2] + rmax));

					// Allocate a ValueAccessor for accelerated random access
					typename GridT::Accessor accessor = mGrid->getAccessor();

					ValueT pval, currR, t;
					const ValueT inside = -mGrid->background();

					if (mInterrupt) mInterrupt->start("Generating level set of capsule");
					// Compute signed distances to sphere using leapfrogging in k
					for (i = imin; i <= imax; ++i) {
						if (util::wasInterrupted(mInterrupt)) return;
						//const float x2 = math::Pow2(i - c[0]);
						const float x = i - s[0];
						const float xdot = x*dir[0];

						for (j = jmin; j <= jmax; ++j) {
							//const float x2y2 = math::Pow2(j - c[1]) + x2;
							const float y = j - s[1];
							const float xydot = y*dir[1] + xdot;
							for (k = kmin; k <= kmax; k += m) {
								m = 1;
								const float z = k - s[2];
								float dot = z*dir[2] + xydot;
								dot = min(max(dot, 0.0f), len);
								t = dot / len;
								const Vec3T projPt = dot*dir - Vec3T(x, y, z);
								currR = t*r2 + (1 - t)*r1;
								/// Distance in voxel units to capsule
								//const float v = math::Sqrt(x2y2 + math::Pow2(k - c[2])) - r0,
								const float v = projPt.length() - currR,
									d = math::Abs(v);
								if (v > w || (!accessor.probeValue(ijk, pval) && pval < 0)) {

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

				const ValueT        mRadius1;
				const ValueT        mRadius2;
				const Vec3T         mStart;
				const Vec3T         mEnd;
				InterruptT*         mInterrupt;
				typename GridT::Ptr mGrid;
			};// LevelSetSphere


			  ////////////////////////////////////////


			template<typename GridType, typename InterruptT>
			typename GridType::Ptr
				createLevelSetCapsuleV(float radius1, float radius2, const openvdb::Vec3f& start, const openvdb::Vec3f& end, float voxelSize,
					float halfWidth, InterruptT* interrupt)
			{
				// GridType::ValueType is required to be a floating-point scalar.
				BOOST_STATIC_ASSERT(boost::is_floating_point<typename GridType::ValueType>::value);

				typedef typename GridType::ValueType ValueT;
				LevelSetCapsuleV<GridType, InterruptT> factory(ValueT(radius1), ValueT(radius2), start, end, interrupt);
				return factory.getLevelSet(ValueT(voxelSize), ValueT(halfWidth));
			}

		} // namespace tools
	} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

#endif // OPENVDB_TOOLS_LEVELSETSPHERE_HAS_BEEN_INCLUDED

  // Copyright (c) 2012-2015 DreamWorks Animation LLC
  // All rights reserved. This software is distributed under the
  // Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
