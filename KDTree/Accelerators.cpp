#include "Primitive.h"
#include "threading.hpp"
#include <iostream>

/// TODO: Implement one/both or any other acceleration structure and change makeDefaultAccelerator to create it
struct KDTree : IntersectionAccelerator {
	struct Node //8 byte to fit into 64-byte cashe line
	{
		union
		{
			float split; //split point in space along a given axis
			int onePrimitive; //index to primitive in allPrimitives if only one is overlapping
			int primitiveIndicesOffset; //offset in primitiveIndices if more than one primitive is overlapping
		};
		union
		{
			int flags; //last two bits hold the split axis (0-x, 1-y, 2-z) when node is interior and value 3 when node is leaf
			int nOverlappingPrim; //number of overlapping primitives, stored in the upper 30 bits
			int farChild; //offset in Nodes array, other child immediately after parent 
		};

		void InitLeaf(int* primNums, int nOvpPrims, std::vector<int>* primitiveIndices)
		{
			flags = 3;
			nOverlappingPrim |= (nOvpPrims << 2);
			if (nOvpPrims == 0)
			{
				onePrimitive = 0;
			}
			else if (nOvpPrims == 1)
			{
				onePrimitive = primNums[0];
			}
			else
			{
				primitiveIndicesOffset = primitiveIndices->size();
				for (int i = 0; i < nOvpPrims; i++)
				{
					primitiveIndices->push_back(primNums[i]);
				}
			}
		}

		void InitInterior(int axis, int fC, float splitPoint)
		{
			split = splitPoint;
			flags = axis;
			flags |= (fC << 2);
		}

		float SplitPos() const { return split; }
		int NumberOfPrim() const { return nOverlappingPrim >> 2; }
		int Axis() const { return flags & 3; }
		bool IsLeaf() const { return (flags & 3) == 3; }
		int FarChild() const { return farChild >> 2; }
	};

	//Struct for edge comparison
	enum class EdgeType { Start, End };
	struct BBoxEdge
	{
		float t;
		int primNum;
		EdgeType type;
		BBoxEdge() {};
		BBoxEdge(float t, int primNum, bool starting)
			: t(t), primNum(primNum) {
			type = starting ? EdgeType::Start : EdgeType::End;
		}
	};

	//Struct for node ordering
	struct KdToDo
	{
		const Node* node;
		float tMin, tMax;
	};

	std::vector<Intersectable*> allPrimitives;
	std::vector<int> primitiveIndices;

	int MAX_DEPTH;
	int MAX_PRIMITIVES;
	int interCost;
	int travCost;
	float emptyBonus;

	Node* nodes = nullptr;
	int nextFreeNode;
	int nAllocNodes;

	BBox bounds;

	void addPrimitive(Intersectable* prim) override {
		allPrimitives.push_back(prim);
	}

	void clear() override {
		_aligned_free(nodes);
		nodes = nullptr;
		allPrimitives.clear();
	}

	void buildHelper(int idx, BBox& nodeBounds, const std::vector<BBox>& allPrimBounds, int* primNums, int nPrimitives, int depth,
		const std::unique_ptr<BBoxEdge[]> edges[3], int* prims0, int* prims1, int badRefines = 0) {

		//Get next free node idx
		if (nextFreeNode == nAllocNodes)
		{
			int nAllocNodesNew = std::max(2 * nAllocNodes, 512);
			Node* temp = (Node*)_aligned_malloc(nAllocNodesNew * sizeof(Node), 64);
			if (nAllocNodes > 0)
			{
				memcpy(temp, nodes, nAllocNodes * sizeof(Node));
				_aligned_free(nodes);
			}
			nodes = temp;
			nAllocNodes = nAllocNodesNew;
		}
		nextFreeNode++;

		//Recursion bottom/ Leaf creation
		if (nPrimitives <= MAX_PRIMITIVES || depth == 0)
		{
			nodes[idx].InitLeaf(primNums, nPrimitives, &primitiveIndices);
			return;
		}

		//Choosing split axis position for interior node using SAH
		int bestAxis = -1;
		int bestOffset = -1;
		float bestCost = std::numeric_limits<float>::infinity();
		float oldCost = interCost * (float)nPrimitives;
		float totalSurfaceArea = nodeBounds.SurfaceArea();
		float invTotalSA = 1 / totalSurfaceArea;
		vec3 diag = nodeBounds.Diagonal();
		int axis = nodeBounds.MaximumExtent();
		int retries = 0;

		while (bestAxis == -1 && retries <= 2)
		{
			//Sort all edges at current axis
			for (int i = 0; i < nPrimitives; i++)
			{
				int pn = primNums[i];
				const BBox& bounds = allPrimBounds[pn];
				edges[axis][2 * i] = BBoxEdge(bounds.min[axis], pn, true);
				edges[axis][2 * i + 1] = BBoxEdge(bounds.max[axis], pn, false);
			}
			std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
				[](const BBoxEdge& e1, const BBoxEdge& e2)->bool
				{
					if (e1.t == e2.t)
					{
						return (int)e1.type < (int)e2.type;
					}
					else
					{
						return e1.t < e2.t;
					}
				});

			//Calculate cost
			int belowAxis = 0; int aboveAxis = nPrimitives;
			for (int i = 0; i < 2 * nPrimitives; i++)
			{
				if (edges[axis][i].type == EdgeType::End)
				{
					aboveAxis--;
				}
				float t = edges[axis][i].t;
				if (t > nodeBounds.min[axis] && t < nodeBounds.max[axis])
				{
					int otherAxis0 = (axis + 1) % 3;
					int otherAxis1 = (axis + 2) % 3;
					float belowSurface = 2 * (diag[otherAxis0] * diag[otherAxis1] + (t - nodeBounds.min[axis]) * (diag[otherAxis0] + diag[otherAxis1]));
					float aboveSurface = 2 * (diag[otherAxis0] * diag[otherAxis1] + (nodeBounds.max[axis] - t) * (diag[otherAxis0] + diag[otherAxis1]));

					float pB = belowSurface * invTotalSA;
					float pA = aboveSurface * invTotalSA;
					float eb = (aboveAxis == 0 || belowAxis == 0) ? emptyBonus : 0;
					float cost = travCost + interCost * (1 - eb) * (pB * belowAxis + pA * aboveAxis);

					//Get best cost
					if (cost < bestCost)
					{
						bestCost = cost;
						bestAxis = axis;
						bestOffset = i;
					}
				}
				if (edges[axis][i].type == EdgeType::Start)
				{
					belowAxis++;
				}
			}


			//If current axis unsatisfactory try other 2
			retries++;
			axis = (axis + 1) % 3;
		}

		//Validation for when no good splits can be made
		if (bestCost > oldCost)
		{
			badRefines++;
		}
		if ((bestCost > 4 * oldCost && nPrimitives < MAX_PRIMITIVES + 5)
			|| bestAxis == -1 || badRefines == 3)
		{
			nodes[idx].InitLeaf(primNums, nPrimitives, &primitiveIndices);
			return;
		}

		//Set primitives according to split
		int n0 = 0, n1 = 0;
		for (int i = 0; i < bestOffset; i++)
		{
			if (edges[bestAxis][i].type == EdgeType::Start)
			{
				prims0[n0++] = edges[bestAxis][i].primNum;
			}
		}
		for (int i = bestOffset + 1; i < 2 * nPrimitives; i++)
		{
			if (edges[bestAxis][i].type == EdgeType::End)
			{
				prims1[n1++] = edges[bestAxis][i].primNum;
			}
		}

		//Init children nodes
		float split = edges[bestAxis][bestOffset].t;
		BBox bounds0 = nodeBounds;
		BBox bounds1 = nodeBounds;
		bounds0.max[bestAxis] = split;
		bounds1.min[bestAxis] = split;
		buildHelper(idx + 1, bounds0, allPrimBounds, prims0, n0, depth - 1, edges, prims0, prims1 + nPrimitives, badRefines);
		int farChild = nextFreeNode;
		nodes[idx].InitInterior(bestAxis, farChild, split);
		buildHelper(farChild, bounds1, allPrimBounds, prims1, n1, depth - 1, edges, prims0, prims1 + nPrimitives, badRefines);
	}

	void build(Purpose purpose = Purpose::Generic) override {
		const char* treePurpose = "";
		if (purpose == Purpose::Instances) {
			MAX_DEPTH = std::round(8 + 1.3f * log2f(allPrimitives.size()));
			MAX_PRIMITIVES = 10;
			treePurpose = " instances";
		}
		else if (purpose == Purpose::Mesh) {
			MAX_DEPTH = std::round(8 + 1.3f * log2f(allPrimitives.size()));
			MAX_PRIMITIVES = 32;
			treePurpose = " mesh";
		}
		else if (purpose == Purpose::Generic)
		{
			MAX_DEPTH = std::round(8 + 1.3f * log2f(allPrimitives.size()));
			MAX_PRIMITIVES = 32;
			treePurpose = " generic";
		}
		interCost = 80;
		travCost = 1;
		emptyBonus = 0.5;
		nextFreeNode = nAllocNodes = 0;


		//Saving all BBoxes to avoid calling expandBox in buildHelper
		std::vector<BBox> allPrimBounds;
		for (Intersectable*& prim : allPrimitives)
		{
			BBox temp;
			prim->expandBox(bounds);
			prim->expandBox(temp);
			allPrimBounds.push_back(temp);
		}

		std::unique_ptr<int[]> primNums(new int[allPrimitives.size()]);
		for (int i = 0; i < allPrimitives.size(); i++)
		{
			primNums[i] = i;
		}

		std::unique_ptr<BBoxEdge[]> edges[3];
		for (int i = 0; i < 3; i++)
		{
			edges[i].reset(new BBoxEdge[2 * allPrimitives.size()]);
		}

		std::unique_ptr<int[]> prims0(new int[allPrimitives.size()]);
		std::unique_ptr<int[]> prims1(new int[(MAX_DEPTH + 1) * allPrimitives.size()]);

		Timer timer;
		buildHelper(0, bounds, allPrimBounds, primNums.get(), allPrimitives.size(), MAX_DEPTH, edges, prims0.get(), prims1.get());
		printf(" done in %lldms, nodes %d\n", timer.toMs(timer.elapsedNs()), nAllocNodes);
	}

	bool isBuilt() const override { return nodes != nullptr; }

	bool intersect(const Ray& ray, float lMin, float lMax, Intersection& intersection) override {
		//Check if ray hits the scene
		float tMin = 0, tMax = 0;
		if (!bounds.intersect(ray, lMin, std::numeric_limits<float>::infinity(), &tMin, &tMax))
		{
			return false;
		}

		float rayMax = std::numeric_limits<float>::infinity();

		vec3 invertedDir(1.0f / ray.dir.x, 1.0f / ray.dir.y, 1.0f / ray.dir.z);
		constexpr int maxTodo = 64;
		KdToDo todo[maxTodo];
		int todoPos = 0;

		bool hit = false;
		const Node* node = &nodes[0];
		while (node != nullptr)
		{
			if (rayMax < tMin)
			{
				break;
			}

			//Interior intersection
			if (!node->IsLeaf())
			{
				int axis = node->Axis();
				float tPlane = (node->SplitPos() - ray.origin[axis]) * invertedDir[axis];
				const Node* firstChild, * secondChild;
				int belowFirst =
					(ray.origin[axis] < node->SplitPos())
					|| (ray.origin[axis] == node->SplitPos() && ray.dir[axis] <= 0);

				//Get children
				if (belowFirst)
				{
					firstChild = node + 1;
					secondChild = &nodes[node->FarChild()];
				}
				else
				{
					firstChild = &nodes[node->FarChild()];
					secondChild = node + 1;
				}

				//Prepare next nodes to be processed
				if (tPlane > tMax || tPlane <= 0) //Ray does not cross second child
				{
					node = firstChild;
				} 
				else if (tPlane < tMin) //Ray does not cross first child
				{
					node = secondChild;
				}
				else //Both children are crossed
				{
					todo[todoPos].node = secondChild;
					todo[todoPos].tMin = tPlane;
					todo[todoPos].tMax = tMax;
					todoPos++;

					node = firstChild;
					tMax = tPlane;
				}
			}
			else //Leaf primitives intersection
			{
				int nPrimitives = node->NumberOfPrim();
				if (nPrimitives == 1)
				{
					Intersection data;
					if (allPrimitives[node->onePrimitive]->intersect(ray, tMin, rayMax, data))
					{
						if (rayMax > data.t)
						{
							rayMax = data.t;
							hit = true;
							intersection = data;
						}
					}
				}
				else
				{
					Intersection data;
					for (int i = 0; i < nPrimitives; i++)
					{
						int index = primitiveIndices[node->primitiveIndicesOffset + i];
						if (allPrimitives[index]->intersect(ray, tMin, rayMax, data))
						{
							if (rayMax > data.t)
							{
								rayMax = data.t;
								hit = true;
								intersection = data;
							}
						}
					} 
				}

				//Load next leaf to be processed
				if (todoPos > 0)
				{
					todoPos--;
					node = todo[todoPos].node;
					tMin = todo[todoPos].tMin;
					tMax = todo[todoPos].tMax;
				}
				else
				{
					break;
				}
			}
		}
		return hit;
	}
};


AcceleratorPtr makeDefaultAccelerator() {
	// TODO: uncomment or add the acceleration structure you have implemented
	return AcceleratorPtr(new KDTree());
	//return AcceleratorPtr(new BVHTree());
	//return AcceleratorPtr(new OctTree());
}

