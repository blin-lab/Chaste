#include "CellLineTensionWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::CellLineTensionWriter() :
AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("celllinetension.dat"),
mAreaElasticityParameter(1.0), // These parameters are Case I in Farhadifar's paper
mPerimeterContractilityParameter(0.04),
mLineTensionParameter(0.12),
mBoundaryLineTensionParameter(0.15) // this parameter as such does not exist in Farhadifar's model.
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
    VertexElement<ELEMENT_DIM, SPACE_DIM>* cell_element = p_cell_population->GetElementCorrespondingToCell(pCell);
    int num_nodes = cell_element->GetNumNodes();


    MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>& mesh = p_cell_population->rGetMesh();

    // Next we iterate over the nodes of that cell to calculate the line tension node by node.
    for (int current_index = 0; current_index < num_nodes; current_index++)
    {
        Node<SPACE_DIM> *current_node = cell_element->GetNode(current_index);
        for (int next_index = 0; next_index < num_nodes; next_index++)
        {
            if (current_index == next_index)
                continue;

            Node<SPACE_DIM> *next_node = cell_element->GetNode(next_index);

            // Find the indices of the elements owned by each node
            std::set<unsigned> elements_containing_nodeA = current_node->rGetContainingElementIndices();
            std::set<unsigned> elements_containing_nodeB = next_node->rGetContainingElementIndices();

            // Find common elements
            std::set<unsigned> shared_elements;
            std::set_intersection(elements_containing_nodeA.begin(), elements_containing_nodeA.end(),
                                  elements_containing_nodeB.begin(), elements_containing_nodeB.end(),
                                  std::inserter(shared_elements, shared_elements.begin()));

            // Check that the nodes have a common edge
            if (shared_elements.empty())
                continue;

            double line_tension_parameter_in_calculation = GetLineTensionParameter();

            // If the edge corresponds to a single element, then the cell is on the boundary
            if (shared_elements.size() == 1)
            {
                line_tension_parameter_in_calculation = GetBoundaryLineTensionParameter();
            }

            // Compute the gradient of each these edges, computed at the present node
            double next_edge_length = mesh.GetDistanceBetweenNodes(current_node->GetIndex(), next_node->GetIndex());
            assert(next_edge_length > DBL_EPSILON);

            c_vector<double, SPACE_DIM> edge_gradient = mesh.GetVectorFromAtoB(current_node->rGetLocation(), next_node->rGetLocation()) / next_edge_length;


            // Add the force contribution from cell-cell and cell-boundary line tension (note the minus sign)
            line_tension_contribution -= line_tension_parameter_in_calculation * edge_gradient;

            total_tension += (line_tension_contribution[0] + line_tension_contribution[1])/2;
            edge_count++;
        }
    }

    return (edge_count);
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
