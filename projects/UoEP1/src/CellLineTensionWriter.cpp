#include "CellLineTensionWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::CellLineTensionWriter() :
        AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("celllinetension.dat"), mAreaElasticityParameter(1.0), // These parameters are Case I in Farhadifar's paper
        mPerimeterContractilityParameter(0.04), mLineTensionParameter(0.12), mBoundaryLineTensionParameter(0.12) // this parameter as such does not exist in Farhadifar's model.
{
    this->mVtkCellDataName = "Cell Tension";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(
        CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> *pCellPopulation)
{

    int total_tension = 0.0;

    VertexBasedCellPopulation<SPACE_DIM> *p_cell_population =
            static_cast<VertexBasedCellPopulation<SPACE_DIM>*>(pCellPopulation);

    c_vector<double, SPACE_DIM> line_tension_contribution = zero_vector<double>(SPACE_DIM);

    unsigned cell_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    std::set<unsigned> node_index = pCellPopulation->GetNeighbouringNodeIndices(cell_index);

    for (std::set<unsigned>::iterator iter_node = ++node_index.begin(); iter_node != node_index.end(); ++iter_node)
    {

        Node<SPACE_DIM> *p_this_node = p_cell_population->GetNode(*iter_node);

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices =
                p_cell_population->GetNode(*iter_node)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin(); iter != containing_elem_indices.end();
                ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<SPACE_DIM, SPACE_DIM> *p_element = p_cell_population->GetElement(*iter);
            unsigned num_nodes_elem = p_element->GetNumNodes();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(*iter_node);

            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem + local_index - 1) % num_nodes_elem;
            Node<SPACE_DIM> *p_previous_node = p_element->GetNode(previous_node_local_index);

            unsigned next_node_local_index = (local_index + 1) % num_nodes_elem;
            Node<SPACE_DIM> *p_next_node = p_element->GetNode(next_node_local_index);

            // Compute the line tension parameter for each of these edges - be aware that this is half of the actual
            // value for internal edges since we are looping over each of the internal edges twice
            double previous_edge_line_tension_parameter = GetLineTensionParameter(p_previous_node, p_this_node,
                                                                                  *p_cell_population);
            double next_edge_line_tension_parameter = GetLineTensionParameter(p_this_node, p_next_node,
                                                                              *p_cell_population);

            // Compute the gradient of each these edges, computed at the present node
            c_vector<double, SPACE_DIM> previous_edge_gradient =
                    -p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element,
                                                                                      previous_node_local_index);
            c_vector<double, SPACE_DIM> next_edge_gradient =
                    p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);

            // Add the force contribution from cell-cell and cell-boundary line tension (note the minus sign)
            line_tension_contribution -= previous_edge_line_tension_parameter * previous_edge_gradient
                    + next_edge_line_tension_parameter * next_edge_gradient;

            total_tension += (line_tension_contribution[0] + line_tension_contribution[1]) / 2.0;
        }

    }

    return (total_tension);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::GetLineTensionParameter(Node<SPACE_DIM> *pNodeA, Node<SPACE_DIM> *pNodeB,
                                                           VertexBasedCellPopulation<SPACE_DIM> &rVertexCellPopulation)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(), elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(), elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    // Since each internal edge is visited twice in the loop above, we have to use half the line tension parameter
    // for each visit.
    double line_tension_parameter_in_calculation = GetLineTensionParameter() / 2.0;

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        line_tension_parameter_in_calculation = GetBoundaryLineTensionParameter();
    }

    return line_tension_parameter_in_calculation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::SetBoundaryLineTensionParameter(double BoundaryLineTensionParameter)
{
    mBoundaryLineTensionParameter = BoundaryLineTensionParameter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::SetLineTensionParameter(double lineTensionParameter)
{
    mLineTensionParameter = lineTensionParameter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::GetLineTensionParameter()
{
    return mLineTensionParameter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::GetBoundaryLineTensionParameter()
{
    return mBoundaryLineTensionParameter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(
        CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> *pCellPopulation)
{
    // Nothing to do
}

// Explicit instantiation
template class CellLineTensionWriter<2, 2> ;


#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
//EXPORT_TEMPLATE_CLASS_2_INTERNAL(CellLineTensionWriter)
