
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// accelerators/kdtreeaccel.cpp*
#include "accelerators/kdtreeaccel.h"
#include "paramset.h"
#include "interaction.h"
#include "stats.h"
#include <algorithm>

namespace pbrt {

// KdTreeAccel Local Declarations
// P285 4.4.1 TREE REPRESENTATION
// Our initial implementation used a 16-byte node representation; when we reduced the size to 8 bytes
// we obtained nearly a 20 % speed increase
// KdAccelNode 为了优化内存访问，降低了可读性
struct KdAccelNode 
{
    // KdAccelNode Methods
    void InitLeaf(int *primNums, int np, std::vector<int> *primitiveIndices);
    void InitInterior(int axis, int ac, Float s)
    { 
        split = s;
        flags = axis;

        // 中间节点需要记录 KdTreeAccel::nodes[] 中两个子节点的位置, 其中一个紧邻它之后, 另一个用 aboveChild 来指示
        aboveChild |= (ac << 2);
    }

    Float SplitPos()  const { return split; } // 在这条坐标轴上划分的位置
    int nPrimitives() const { return nPrims >> 2; } // nPrims 的高三十位记录与叶子节点相交的图元数量
    int SplitAxis()   const { return flags & 3; } // flags 的最低两位记录了 split axis，0/1/2代表基于x/y/z 轴划分
    bool IsLeaf()     const { return (flags & 3) == 3; } // 当最低两位为 3 时， 代表这是一个叶子节点
    int AboveChild()  const { return aboveChild >> 2; }

    union 
    {
        Float split;                 // Interior

        int onePrimitive;            // Leaf, 如果只有一个图元与叶节点相交, 记录这个图元在 KdTreeAccel::primitives 中的索引, 当 onePrimitive 为0时表示没有与任何图元相交
        int primitiveIndicesOffset;  // Leaf, 否则在 KdTreeAccel::primitiveIndices 中连续记录这一组图元的下标, primitiveIndicesOffset 记录了 primitiveIndices 中的首个下标
    };

  private:
    // 这个结构的高三十位被 nPrims/aboveChild 使用, 最低两位被 flags 使用
    union 
    {
        int flags;       // Both
        int nPrims;      // Leaf
        int aboveChild;  // Interior
    };
};

// 记录坐标轴上包围盒的边界, 参考 Figure4.15
enum class EdgeType { Start, End };
struct BoundEdge 
{
    // BoundEdge Public Methods
    BoundEdge() {}
    BoundEdge(Float t, int primNum, bool starting) : t(t), primNum(primNum) 
    {
        type = starting ? EdgeType::Start : EdgeType::End;
    }

    Float t;
    int primNum;
    EdgeType type;
};

// KdTreeAccel Method Definitions
KdTreeAccel::KdTreeAccel(std::vector<std::shared_ptr<Primitive>> p,
                         int isectCost, int traversalCost, Float emptyBonus,
                         int maxPrims, int maxDepth)
    : isectCost(isectCost),
      traversalCost(traversalCost),
      maxPrims(maxPrims),
      emptyBonus(emptyBonus),
      primitives(std::move(p)) 
{
    // Build kd-tree for accelerator
    ProfilePhase _(Prof::AccelConstruction);

    nextFreeNode = nAllocedNodes = 0;
    // kdtree 建树时需要有一个最大的递归深度
    if (maxDepth <= 0)
        maxDepth = std::round(8 + 1.3f * Log2Int(int64_t(primitives.size())));

    // Compute bounds for kd-tree construction
    // optimize: 预先缓存一下图元的包围盒
    std::vector<Bounds3f> primBounds;
    primBounds.reserve(primitives.size());
    for (const std::shared_ptr<Primitive> &prim : primitives) 
    {
        Bounds3f b = prim->WorldBound();
        bounds = Union(bounds, b);
        primBounds.push_back(b);
    }

    // Allocate working memory for kd-tree construction
    // 记录x/y/z 三条轴上的包围盒的边界信息, 参考 Figure4.15
    // 一开始就按最大容量来分配, 这样只需要分配一次, 后续复用这段空间就行了
    std::unique_ptr<BoundEdge[]> edges[3];
    for (int i = 0; i < 3; ++i)
        edges[i].reset(new BoundEdge[2 * primitives.size()]);

    // 记录与两个子节点相交的图元, 也是一开始就按最大容量分配的内存
    // 递归划分左子树时, 不能覆盖掉右子树的数据(prims1), 所以这里 prims 分配了一个更大的空间
    std::unique_ptr<int[]> prims0(new int[primitives.size()]);
    std::unique_ptr<int[]> prims1(new int[(maxDepth + 1) * primitives.size()]);

    // Initialize _primNums_ for kd-tree construction
    // primNums 记录与**当前节点**相交的图元的索引, 初始时所有图元都和根节点相交
    std::unique_ptr<int[]> primNums(new int[primitives.size()]);
    for (size_t i = 0; i < primitives.size(); ++i) 
        primNums[i] = i;

    // Start recursive construction of kd-tree
    buildTree(
        0, 
        bounds, primBounds, 
        primNums.get(), primitives.size(),
        maxDepth, // 传递最大递归深度, 在递归调用 buildTree 中, 当 depth == 0 时, 停止递归
        edges, prims0.get(), prims1.get());
}

void KdAccelNode::InitLeaf(int *primNums, int np,
                           std::vector<int> *primitiveIndices) {
    flags = 3;
    nPrims |= (np << 2);

    // Store primitive ids for leaf node
    if (np == 0) // 空节点
        onePrimitive = 0;
    else if (np == 1) // 只与一个图元相交
        onePrimitive = primNums[0];
    else {
        primitiveIndicesOffset = primitiveIndices->size();
        for (int i = 0; i < np; ++i) primitiveIndices->push_back(primNums[i]);
    }
}

KdTreeAccel::~KdTreeAccel() { FreeAligned(nodes); }

void KdTreeAccel::buildTree(
    int nodeNum, 
    const Bounds3f &nodeBounds, const std::vector<Bounds3f> &allPrimBounds,
    int *primNums, int nPrimitives, 
    int depth,
    const std::unique_ptr<BoundEdge[]> edges[3],
    int *prims0, int *prims1, int badRefines) 
{
    CHECK_EQ(nodeNum, nextFreeNode);

    // Get next free node from _nodes_ array
    // kdtree 在建树时无法预先知道要分配多大的内存, 需要中间不断去扩张容量
    if (nextFreeNode == nAllocedNodes) 
    {
        int nNewAllocNodes = std::max(2 * nAllocedNodes, 512); // 每次按两倍扩容
        KdAccelNode *n = AllocAligned<KdAccelNode>(nNewAllocNodes);

        if (nAllocedNodes > 0) {
            memcpy(n, nodes, nAllocedNodes * sizeof(KdAccelNode));
            FreeAligned(nodes);
        }
        nodes = n;
        nAllocedNodes = nNewAllocNodes;
    }
    ++nextFreeNode;

    // Initialize leaf node if termination criteria met
    // 这个节点下的图元太少, 或者达到最大递归深度, 就不继续分割了
    if (nPrimitives <= maxPrims || depth == 0) 
    {
        // 与当前节点相交的图元, 图元数量, 要写入的 KdTreeAccel::primitiveIndices
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // Initialize interior node and continue recursion

    // Choose split axis position for interior node
    // 仍然使用 Section 4.3.2 中的 SAH 来评估不同划分方式的开销
    // 在x/y/z轴上分别选取一个划分位置并计算开销, 选择开销最低的那个来划分空间
    int bestAxis = -1, bestOffset = -1; // bestOffset: splitOffset in bestAxis

    Float bestCost = Infinity;
    Float oldCost = isectCost * Float(nPrimitives); // defalutCost, 不进行划分, 当作叶节点的开销

    Float totalSA = nodeBounds.SurfaceArea();
    Float invTotalSA = 1 / totalSA;

    Vector3f d = nodeBounds.pMax - nodeBounds.pMin;

    // Choose which axis to split along
    int axis = nodeBounds.MaximumExtent();
    int retries = 0; // retrySplit 的次数

retrySplit:

    // Initialize edges for _axis_
    // kdtree 总是在一条坐标轴上包围盒的边界处进行划分, 参考 Figure 4.15
    for (int i = 0; i < nPrimitives; ++i) 
    {
        int pn = primNums[i];
        const Bounds3f &bounds = allPrimBounds[pn];
        edges[axis][2 * i] = BoundEdge(bounds.pMin[axis], pn, true);
        edges[axis][2 * i + 1] = BoundEdge(bounds.pMax[axis], pn, false);
    }

    // Sort _edges_ for _axis_
    // 按从负到正的方向排序
    std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
              [](const BoundEdge &e0, const BoundEdge &e1) -> bool 
              {
                  if (e0.t == e1.t)
                      return (int)e0.type < (int)e1.type;
                  else
                      return e0.t < e1.t;
              });

    // Compute cost of all splits for _axis_ to find best
    int nBelow = 0, nAbove = nPrimitives;
    for (int i = 0; i < 2 * nPrimitives; ++i) 
    {
        // if(i % 2 == 1)
        if (edges[axis][i].type == EdgeType::End) 
            --nAbove;

        Float edgeT = edges[axis][i].t;
        if (edgeT > nodeBounds.pMin[axis] && edgeT < nodeBounds.pMax[axis]) 
        {
            // Compute cost for split at _i_th edge

            // Compute child surface areas for split at _edgeT_
            int otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;

            // Vector3f d = nodeBounds.pMax - nodeBounds.pMin;
            // 注意这里算的是表面积, 而不是体积
            Float belowSA = 2 * ( d[otherAxis0] * d[otherAxis1] + (edgeT - nodeBounds.pMin[axis]) * (d[otherAxis0] + d[otherAxis1]) );
            Float aboveSA = 2 * ( d[otherAxis0] * d[otherAxis1] + (nodeBounds.pMax[axis] - edgeT) * (d[otherAxis0] + d[otherAxis1]) );

            Float pBelow = belowSA * invTotalSA;
            Float pAbove = aboveSA * invTotalSA;

            // 如果有一个子节点是空的, 就可以跳过这个节点, 提升遍历性能
            // 所以 eb(bonus) 在这里可以看做是一个"奖励", 降低了 cost
            Float eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0; // emptyBonus 的默认值是 0.5f
            Float cost =
                traversalCost +
                isectCost * (1 - eb) * (pBelow * nBelow + pAbove * nAbove);

            // Update best split if this is lowest cost so far
            if (cost < bestCost) 
            {
                bestCost = cost;
                bestAxis = axis;
                bestOffset = i;
            }
        }

        if (edges[axis][i].type == EdgeType::Start) 
            ++nBelow;
    }
    // 这个 CHEKC 对应初始的 'int nBelow = 0, nAbove = nPrimitives;'
    CHECK(nBelow == nPrimitives && nAbove == 0);

    // Create leaf if no good splits were found
    // 为什么会出现 Figure4.16 这种情况, 不应该所有 primBounds 都与 nodeBounds 包含或相交吗???
    if (bestAxis == -1 && retries < 2) 
    {
        ++retries;
        axis = (axis + 1) % 3;
        goto retrySplit;
    }

    // 即便进行划分的开销(bestCost)比默认的开销(oldCost)还要大, 仍允许划分, 因为在后续递归生成子树的过程中可能会生成更好的结果
    // 但是要记录一下这种情况发生的次数
    if (bestCost > oldCost) ++badRefines;

    // 当节点数量太少时, 还是生成叶节点算了
    // 又或者向下递归的过程中, 如果进行划分老是比不划分的开销还大(badRefines == 3), 就不再划分, 而是直接生成叶节点
    if ((bestCost > 4 * oldCost && nPrimitives < 16) || bestAxis == -1 ||
        badRefines == 3) 
    {
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // Classify primitives with respect to split
    // 计算与两个子节点相交的图元索引, 这个索引对应的是 KdTreeAccel::primitives 和 allPrimBounds
    int n0 = 0, n1 = 0;

    for (int i = 0;              i < bestOffset;      ++i)
        if (edges[bestAxis][i].type == EdgeType::Start)
            prims0[n0++] = edges[bestAxis][i].primNum;

    for (int i = bestOffset + 1; i < 2 * nPrimitives; ++i)
        if (edges[bestAxis][i].type == EdgeType::End)
            prims1[n1++] = edges[bestAxis][i].primNum;

    // Recursively initialize children nodes
    // 直接对 nodeBounds 划分产生新的 subNodeBounds
    Float tSplit = edges[bestAxis][bestOffset].t;
    Bounds3f bounds0 = nodeBounds, bounds1 = nodeBounds;
    bounds0.pMax[bestAxis] = bounds1.pMin[bestAxis] = tSplit;

    // 划分左子树时, 虽然传入的 primNums 和 prims0 是同一个, 但在 'Classify primitives with respect to split' 时, primNums 已经无用了, 可以覆盖掉它的数据
    // 左子树划分完, 开始划分右子树时, 还需要用到 prims1 的数据, prims1 的数据不能被覆盖, 所以传入的是 'prims1 + nPrimitives'

    // 'nodeNum + 1' -> 这个子节点紧邻父节点, 不需要额外记录, 对另一个子节点, 就需要 'InitInterior' 记录一下了
    buildTree(nodeNum + 1, bounds0, allPrimBounds, prims0, n0, depth - 1, edges,
              prims0, prims1 + nPrimitives, badRefines); // 注意这里传入的是 'prims1 + nPrimitives'

    int aboveChild = nextFreeNode; // 深度遍历完前一个节点后, 产生新的 nextFreeNode value
    nodes[nodeNum].InitInterior(bestAxis, aboveChild, tSplit);

    buildTree(aboveChild,  bounds1, allPrimBounds, prims1, n1, depth - 1, edges,
              prims0, prims1 + nPrimitives, badRefines);
}



bool KdTreeAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const 
{
    ProfilePhase p(Prof::AccelIntersect);

    // 参考 Figure 4.17 的四个过程, 仔细关注 ray 的 [tMin, tMax] 的变化, 通过 [tMin, tMax] 可以更快的结束遍历

    // Compute initial parametric range of ray inside kd-tree extent
    Float tMin, tMax;
    if (!bounds.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }

    // Prepare to traverse kd-tree for ray
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    PBRT_CONSTEXPR int maxTodo = 64;
    KdToDo todo[maxTodo]; // 又是一个通过数组模拟的 stack, 可以对比一下 'BVHAccel::Intersect'
    int todoPos = 0;

    // Traverse kd-tree nodes in order for ray
    bool hit = false;
    const KdAccelNode *node = &nodes[0];
    
    // 理解这个 while 循环时, 不妨先当成遍历二叉树, 搞清楚递归, 回溯, 终止的流程, 然后再考虑 kdtree 的空间结构

    while (node != nullptr) 
    {
        // Bail out if we found a hit closer than the current node
        // 如果已经找到了一个比当前访问节点还要近的交点, 就不需要再查找了, 跳出循环
        if (ray.tMax < tMin) 
            break;

        if (!node->IsLeaf()) // 如果是中间节点, 就继续访问它最近的子节点, 向下递归, 直到找到了最近的叶节点为止
        {
            // Process kd-tree interior node

            // Compute parametric distance along ray to split plane
            int axis = node->SplitAxis();
            Float tPlane = (node->SplitPos() - ray.o[axis]) * invDir[axis];

            // Get node children pointers for ray
            const KdAccelNode *firstChild, *secondChild;

            // 参考 Figure 4.18, 先访问更近的那个子节点
            int belowFirst =
                (ray.o[axis] < node->SplitPos()) ||
                (ray.o[axis] == node->SplitPos() && ray.d[axis] <= 0); // 当光线原点在划分平面上, 则根据光线方向来选择

            if (belowFirst) 
            {
                firstChild = node + 1;
                secondChild = &nodes[node->AboveChild()];
            } 
            else 
            {
                firstChild = &nodes[node->AboveChild()];
                secondChild = node + 1;
            }

            // Advance to next child node, possibly enqueue other child
            // 参考 Figure4.19, 前两种情况下就不需要考虑访问另一个节点了, 否则把另一个节点加到 todoList 里
            if (tPlane > tMax || tPlane <= 0)
                node = firstChild;
            else if (tPlane < tMin)
                node = secondChild;
            else 
            {
                // Enqueue _secondChild_ in todo list
                todo[todoPos].node = secondChild;
                // 暂存 secondChild 的 [tMin, tMax], 然后更新当今节点的 tMax
                todo[todoPos].tMin = tPlane;
                todo[todoPos].tMax = tMax;
                ++todoPos;

                node = firstChild;
                tMax = tPlane;
            }
        } 
        else // 是叶节点的话, 遍历与这个节点相关的图元, 一旦找到了相交的图元, 就更新 ray.tMax, 看能不能跳出循环, 否则回溯到最近的中间节点去
        {
            // Check for intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) 
            {
                const std::shared_ptr<Primitive> &p =
                    primitives[node->onePrimitive];

                // Check one primitive inside leaf node
                if (p->Intersect(ray, isect)) 
                    hit = true;
            } 
            else 
            {
                for (int i = 0; i < nPrimitives; ++i) 
                {
                    int index =
                        primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &p = primitives[index];

                    // Check one primitive inside leaf node
                    if (p->Intersect(ray, isect)) 
                        hit = true;
                }
            }

            // Grab next node to process from todo list
            if (todoPos > 0) 
            {
                --todoPos;
                node = todo[todoPos].node;

                // 开始访问另一节点前, 先更新 ray 在该阶段的 [tMin, tMax]
                tMin = todo[todoPos].tMin;
                tMax = todo[todoPos].tMax;
            } 
            else
                break;
        }
    }

    return hit;
}

bool KdTreeAccel::IntersectP(const Ray &ray) const {
    ProfilePhase p(Prof::AccelIntersectP);
    // Compute initial parametric range of ray inside kd-tree extent
    Float tMin, tMax;
    if (!bounds.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }

    // Prepare to traverse kd-tree for ray
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    PBRT_CONSTEXPR int maxTodo = 64;
    KdToDo todo[maxTodo];
    int todoPos = 0;
    const KdAccelNode *node = &nodes[0];
    while (node != nullptr) {
        if (node->IsLeaf()) {
            // Check for shadow ray intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const std::shared_ptr<Primitive> &p =
                    primitives[node->onePrimitive];
                if (p->IntersectP(ray)) {
                    return true;
                }
            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    int primitiveIndex =
                        primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &prim =
                        primitives[primitiveIndex];
                    if (prim->IntersectP(ray)) {
                        return true;
                    }
                }
            }

            // Grab next node to process from todo list
            if (todoPos > 0) {
                --todoPos;
                node = todo[todoPos].node;
                tMin = todo[todoPos].tMin;
                tMax = todo[todoPos].tMax;
            } else
                break;
        } else {
            // Process kd-tree interior node

            // Compute parametric distance along ray to split plane
            int axis = node->SplitAxis();
            Float tPlane = (node->SplitPos() - ray.o[axis]) * invDir[axis];

            // Get node children pointers for ray
            const KdAccelNode *firstChild, *secondChild;
            int belowFirst =
                (ray.o[axis] < node->SplitPos()) ||
                (ray.o[axis] == node->SplitPos() && ray.d[axis] <= 0);
            if (belowFirst) {
                firstChild = node + 1;
                secondChild = &nodes[node->AboveChild()];
            } else {
                firstChild = &nodes[node->AboveChild()];
                secondChild = node + 1;
            }

            // Advance to next child node, possibly enqueue other child
            if (tPlane > tMax || tPlane <= 0)
                node = firstChild;
            else if (tPlane < tMin)
                node = secondChild;
            else {
                // Enqueue _secondChild_ in todo list
                todo[todoPos].node = secondChild;
                todo[todoPos].tMin = tPlane;
                todo[todoPos].tMax = tMax;
                ++todoPos;
                node = firstChild;
                tMax = tPlane;
            }
        }
    }
    return false;
}



std::shared_ptr<KdTreeAccel> CreateKdTreeAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) 
{
    int isectCost = ps.FindOneInt("intersectcost", 80);
    int travCost = ps.FindOneInt("traversalcost", 1);
    Float emptyBonus = ps.FindOneFloat("emptybonus", 0.5f);
    int maxPrims = ps.FindOneInt("maxprims", 1);
    int maxDepth = ps.FindOneInt("maxdepth", -1);
    return std::make_shared<KdTreeAccel>(std::move(prims), isectCost, travCost, emptyBonus,
                                         maxPrims, maxDepth);
}

}  // namespace pbrt
