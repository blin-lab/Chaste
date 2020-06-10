#include "CellLineTensionWriter.hpp"
#include "AbstractCellPopulation.hpp"

using namespace std;




template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::CellLineTensionWriter()
: AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("celllinetension.dat")
  {
    this->mVtkCellDataName = "Tension";
  }

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{

    int total_tension = 0.0;

    unsigned cell_index = pCellPopulation->GetLocationIndexUsingCell();
    std::set<unsigned> node_indices = pCellPopulation->GetNeighbouringNodeIndices(cell_index);

    for (std::set<unsigned>::iterator iter = node_indices.begin();
         iter != node_indices.end();
         ++iter)
    {
        pNodeA = pCellPopulation->GetNode(iter);
        pNodeB = pCellPopulation->GetNode(iter);
        // Find the indices of the elements owned by each node
        std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices(); //look at the farhadifar force and the local index and node index stuff in line 136
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
    }

    total_tension += line_tension_parameter_in_calculation;
    return pCell->GetCellData()->SetItem("Line Tension", total_tension);}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    // Write location index corresponding to cell
    *this->mpOutStream << pCellPopulation->GetLocationIndexUsingCell(pCell) << " ";

    // Write cell location
    c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << cell_location[i] << " ";
    }

    // Write cell age
    *this->mpOutStream << pCell->GetCellData()->SetItem("Line Tension", total_tension) << " ";
}

// Explicit instantiation
template class CellLineTensionWriter<1,1>;
template class CellLineTensionWriter<1,2>;
template class CellLineTensionWriter<2,2>;
template class CellLineTensionWriter<1,3>;
template class CellLineTensionWriter<2,3>;
template class CellLineTensionWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellLineTensionWriter)
