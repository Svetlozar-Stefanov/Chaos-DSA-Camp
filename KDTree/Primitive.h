#pragma once

#include "MemUtils.h"
#include "Utils.hpp"
#include "Material.h"

#include <vector>
#include <memory>

/// Data for an intersection between a ray and scene primitive
struct Intersection {
	float t = -1.f; ///< Position of the intersection along the ray
	vec3 p; ///< The intersection point
	vec3 normal; ///< The normal at the intersection
	Material *material = nullptr; ///< Material of the intersected primitive
};

/// Interface for anything that can be intersected
struct Intersectable {
	/// @brief Intersect a ray with the primitive allowing intersection in (tMin, tMax) along the ray
	/// @param ray - the ray
	/// @param tMin - near clip distance
	/// @param tMax - far clip distance
	/// @param intersection [out] - data for intersection if one is found
	/// @return true when intersection is found, false otherwise
	virtual bool intersect(const Ray &ray, float tMin, float tMax, Intersection &intersection) = 0;

	/// @brief Test intersection of the primitive with a box, used by IntersectionAccelerator
	/// @param box - bounding box to test against
	virtual bool boxIntersect(const BBox &box) = 0;

	/// @brief Expand given bounding box with all points part of the Intersectable
	/// @param box [out] - the box to expand
	virtual void expandBox(BBox &box) = 0;

	//virtual const BBox& getBBox() const = 0;

	virtual ~Intersectable() = default;
};

/// Base class for scene object
struct Primitive : Intersectable {
	BBox box;

	/// @brief Called after scene is fully created and before rendering starts
	///	       Used to build acceleration structures
	virtual void onBeforeRender() {}

	/// @brief Default implementation intersecting the bbox of the primitive, overriden if possible more efficiently
	bool boxIntersect(const BBox& other) override {
		return !box.boxIntersection(other).isEmpty();
	}

	/// @brief Default implementation adding the whole bbox, overriden for special cases
	void expandBox(BBox &other) override {
		other.add(box);
	}

	~Primitive() override = default;
};

typedef std::unique_ptr<Primitive> PrimPtr;
typedef std::shared_ptr<Primitive> SharedPrimPtr;

/// Interface for an acceleration structure for any intersectable primitives
struct IntersectionAccelerator {
	enum class Purpose {
		Generic,
		Mesh,
		Instances
	};

	/// @brief Add the primitive to the accelerated list
	/// @param prim - non owning pointer
	virtual void addPrimitive(Intersectable *prim) = 0;

	/// @brief Clear all data allocated by the accelerator
	virtual void clear() = 0;

	/// @brief Build all the internal data for the accelerator
	///	@param purpose - the purpose of the tree, implementation can use it as hint for internal parameters
	virtual void build(Purpose purpose = Purpose::Generic) = 0;

	/// @brief Check if the accelerator is built
	virtual bool isBuilt() const = 0;

	/// @brief Implement intersect from Intersectable but don't inherit the Interface
	virtual bool intersect(const Ray &ray, float tMin, float tMax, Intersection &intersection) = 0;

	virtual ~IntersectionAccelerator() = default;
};

typedef std::unique_ptr<IntersectionAccelerator> AcceleratorPtr;
AcceleratorPtr makeDefaultAccelerator();

/// Simple smooth sphere primitive
struct SpherePrim : Primitive {
	vec3 center;
	float radius;
	std::unique_ptr<Material> material;

	SpherePrim(vec3 center, float radius, MaterialPtr material);

	bool intersect(const Ray &ray, float tMin, float tMax, Intersection &intersection) override;
};

/// Primitive that contains a list of other primitives along with offset and scale for each one
///	Each primitive is tested on intersect call and intersected with its offset and scale
struct Instancer : Primitive {
private:
	struct Instance : Intersectable {
		SharedPrimPtr primitive;
		vec3 offset;
		float scale;
		SharedMaterialPtr material;

		bool intersect(const Ray &ray, float tMin, float tMax, Intersection &intersection) override;
		bool boxIntersect(const BBox &other) override;
		void expandBox(BBox &other) override;
	};
	std::vector<Instance> instances;

	AcceleratorPtr accelerator;
public:
	void onBeforeRender() override;

	void addInstance(SharedPrimPtr prim, const vec3 &offset = vec3(0.f), float scale = 1.f, SharedMaterialPtr material = nullptr);

	bool intersect(const Ray &ray, float tMin, float tMax, Intersection &intersection) override;
};
