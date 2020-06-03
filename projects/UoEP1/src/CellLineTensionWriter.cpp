#include "CellLineTensionWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::CellLineTensionWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("celllinetension.dat")
{
    this->mVtkCellDataName = "Tension";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellLineTensionWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    return pCell->GetForceContribution();
}

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
    *this->mpOutStream << pCell->GetForceContribution << " ";
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
