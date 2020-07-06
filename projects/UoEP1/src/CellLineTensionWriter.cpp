#include "CellLineTensionWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::CellLineTensionWriter() :
AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("celllinetension.dat"),
mAreaElasticityParameter(1.0), // These parameters are Case I in Farhadifar's paper
mPerimeterContractilityParameter(0.04),
mLineTensionParameter(0.12),
mBoundaryLineTensionParameter(0.12) // this parameter as such does not exist in Farhadifar's model.
{
    this->mVtkCellDataName = "Cell Tension";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(
        CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> *pCellPopulation) //
{

    double total_tension = 0.0; // Declare the variable total_tension which will later be used to do the sum of all the line tensions.
    unsigned edge_count = 0;

    //We have an Abstract Cell Population and want to specify that what we will be using in this case is a Vertex Based cell population.
    //This conversion can be done using static cast.
    VertexBasedCellPopulation<SPACE_DIM> *p_cell_population =
            static_cast<VertexBasedCellPopulation<SPACE_DIM>*>(pCellPopulation); //

    c_vector<double, SPACE_DIM> line_tension_contribution = zero_vector<double>(SPACE_DIM);

    // Find the cell by: accessing members of p_Cell_Population through a pointer with argument pCell.
    // And find the element attached to that cell.
    VertexElement <ELEMENT_DIM, SPACE_DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(pCell);
    // Get this element, its index and its number of nodes
    unsigned num_nodes_elem = p_element->GetNumNodes();

    // Iterate over vertices in the cell population
    for (unsigned local_index=0; local_index<num_nodes_elem; local_index++)
    {

        Node<SPACE_DIM>* p_this_node = p_element->GetNode(local_index);
        // Get the previous and next nodes in this element
        unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
        Node<SPACE_DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);

        unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
        Node<SPACE_DIM>* p_next_node = p_element->GetNode(next_node_local_index);

        // Compute the line tension parameter for each of these edges - be aware that this is half of the actual
        // value for internal edges since we are looping over each of the internal edges twice
        double previous_edge_line_tension_parameter = GetLineTensionParameter(p_previous_node, p_this_node, *p_cell_population);
        double next_edge_line_tension_parameter = GetLineTensionParameter(p_this_node, p_next_node, *p_cell_population);

        // Compute the gradient of each these edges, computed at the present node
        c_vector<double, SPACE_DIM> previous_edge_gradient =
                -p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
        c_vector<double, SPACE_DIM> next_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);

        // Add the force contribution from cell-cell and cell-boundary line tension (note the minus sign)
        line_tension_contribution -= previous_edge_line_tension_parameter*previous_edge_gradient +
                next_edge_line_tension_parameter*next_edge_gradient;

        total_tension += sqrt((line_tension_contribution[0] * line_tension_contribution[0]) + (line_tension_contribution[1] * line_tension_contribution[1]));
        edge_count++;

    }

    return (total_tension);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::GetLineTensionParameter(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, VertexBasedCellPopulation<SPACE_DIM>& rVertexCellPopulation)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    // Since each internal edge is visited twice in the loop above, we have to use half the line tension parameter
    // for each visit.
    double line_tension_parameter_in_calculation = GetLineTensionParameter()/2.0;

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
