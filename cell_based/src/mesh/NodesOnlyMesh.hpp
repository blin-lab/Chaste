/*

Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef NODESONLYMESH_HPP_
#define NODESONLYMESH_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>
#include "BoxCollection.hpp"
#include "MutableMesh.hpp"

/**
 * Mesh class for storing lists of nodes (no elements). This inherits from MutableMesh
 * because we want to be able to add and delete nodes.
 */
template<unsigned SPACE_DIM>
class NodesOnlyMesh: public MutableMesh<SPACE_DIM, SPACE_DIM>
{

protected:
    /**
     * A pointer to a box collection. Used to calculate neighbourhood information
     * for nodes in the mesh.
     */
    BoxCollection<SPACE_DIM>* mpBoxCollection;

private:

    /** The global number of nodes in the mesh. */
    unsigned mTotalNumNodes;

    /** Vector of pointer to halo nodes used by this process. */
    std::vector<Node<SPACE_DIM>* > mHaloNodes;

    /**
     * Defines connectivity in NodesOnlyMesh. Two nodes are connected
     * if their centres are less than mMaximumInteractionDistance apart.
     */
    double mMaximumInteractionDistance;

    /** A map from node global index to local index used by this process. */
    std::map<unsigned, unsigned> mNodesMapping;

    /** A map from halo node global index to local index used by this process. */
    std::map<unsigned, unsigned> mHaloNodesMapping;

    /** A counter of the number of fresh indices used on this process. */
    unsigned mIndexCounter;

    /** A vector of the global indices that have been freed on the process from deleting a node */
    //std::vector<unsigned> mDeletedGlobalNodeIndices;

    friend class TestNodesOnlyMesh;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archives the member variables of the object which have to be preserved
     * during its lifetime.
     *
     * Note that we must archive any member variables FIRST so that this
     * method can call a ReMesh (to convert from TrianglesMeshReader input
     * format into our native format).
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<MutableMesh<SPACE_DIM, SPACE_DIM> >(*this);
        archive & mMaximumInteractionDistance;
    }

    /**
     * Calculate the next unique global index available on this
     * process. Uses a hashing function to ensure that a unqiue
     * index is given to every node.
     *
     * For example for 3 process they will have access to the following
     * integers:
     *
     * Proc 0:  0   3   6   9   12  ...
     *
     * Proc 1:  1   4   7   10   13  ...
     *
     * Proc 2:  2   5   8   11   14  ...
     *
     * Deleted node indices can be locally re-used.
     *
     * Deleted node inidces of nodes that have *moved* process cannot be re-used.
     *
     * @return the next available index.
     */
    unsigned GetNextAvailableIndex();

public:

    /**
     * Default constructor to initialise BoxCollection to NULL.
     */
    NodesOnlyMesh();

    /**
     * Over-written destructor to delete pointer to BoxCollection
     */
    virtual ~NodesOnlyMesh();

    /**
     * Construct the mesh using only nodes. No mesh is created, but the nodes are stored.
     * The original vector of nodes is deep-copied: new node objects are made with are
     * independent of the pointers in the input so that they can be safely deleted.
     *
     * If this is the only way of constructing a mesh of this type, then we can be certain that
     * elements and boundary elements are always unused.
     *
     * @param rNodes a vector of pointers to nodes
     * @param maxInteractionDistance the distance that defines node neighbours in CalculateNodePairs
     * @param domainPadding the amount of padding space added to the edge of the spatial domain on construction. Should be larger than max movement distance of a node in one step.
     */
    void ConstructNodesWithoutMesh(const std::vector<Node<SPACE_DIM>*>& rNodes, double maxInteractionDistance, double domainPadding = 2.0);

    /**
     * A Helper method to enable you to construct a nodes-only mesh by stripping the nodes
     * TetrahedralMesh, this calls the ConstructNodesWithoutMesh method with the nodes
     *
     * If this is the only way of constructing a mesh of this type, then we can be certain that
     * elements and boundary elements are always unused.
     *
     * @param rGeneratingMesh any mesh with nodes, used to generate the NodesOnlyMesh
     * @param maxInteractionDistance the distance that defines node neighbours in CalculateNodePairs
     * @param domainPadding the amount of padding space added to the edge of the spatial domain on construction. Should be larger than max movement distance of a node in one step.
     */
    void ConstructNodesWithoutMesh(const AbstractMesh<SPACE_DIM,SPACE_DIM>& rGeneratingMesh, double maxInteractionDistance, double domainPadding = 2.0);

    /**
     * Overridden Clear() method for NodesOnlyMesh.
     */
    void Clear();

    /**
     * Overridden solve node mapping method
     *
     * @param index the global index of the node
     *
     * @return the local index of the node.
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Get the local number of nodes that are actually in use.
     * Does not include halo nodes.
     */
    unsigned GetNumNodes() const;

    /**
     * Get the global number of nodes.
     */
    unsigned GetGlobalNumNodes() const;

    /**
     * @return mMaxInteractionDistance
     */
    double GetMaximumInteractionDistance();

    /**
     * Get mpBoxCollection
     *
     * @return mpBoxCollection
     */
    BoxCollection<SPACE_DIM>* GetBoxCollection();

    /**
     * Clear the BoxCollection
     */
    void ClearBoxCollection();

    /**
     * Set up the box collection. Overridden in subclasses to implement periodicity.
     *
     * @param cutOffLength the cut off length for node neighbours
     * @param domainSize the size of the domain containing the nodes.
     */
    virtual void SetUpBoxCollection(double cutOffLength, c_vector<double, 2*SPACE_DIM> domainSize);

    /**
     * Calculate pairs of nodes using the BoxCollection
     *
     * @param rNodePairs reference to the set of node pairs to populate.
     * @param rNodeNeighbours reference to the list of neighbouring nodes for each node.
     */
    void CalculateNodePairs(std::set<std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*> >& rNodePairs, std::map<unsigned, std::set<unsigned> >& rNodeNeighbours);

    /**
     * Overridden ReMesh() method. Since only Nodes are stored, this method just cleans up mNodes by
     * removing nodes marked as deleted and reallocating mNodes to 'fill the gaps'.
     *
     * @param rMap a reference to a NodeMap which associated the indices of the old mesh
     * with the new mesh. It should have the same size as mNodes.
     */
    void ReMesh(NodeMap& rMap);

    /**
     * Overridden AddNode() method.
     *
     * @param pNewNode  pointer to the new node
     */
    unsigned AddNode(Node<SPACE_DIM>* pNewNode);

    /**
     * Overridden DeleteNode() method.
     *
     * @param index is the index of the node to be deleted
     */
    void DeleteNode(unsigned index);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodesOnlyMesh)

#endif /*NODESONLYMESH_HPP_*/
