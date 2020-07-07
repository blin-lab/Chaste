#include "CircularBoundaryCondition.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CircularBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::CircularBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* pCellPopulation,
                                                                      c_vector<double, SPACE_DIM> centre,
                                                                      double radius,
                                                                      double distance)
    : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>(pCellPopulation),
      mCentreOfCircle(centre),
      mRadiusOfCircle(radius),
      mMaximumDistance(distance)
{
    assert(mRadiusOfCircle > 0.0);
    assert(mMaximumDistance > 0.0);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const c_vector<double, SPACE_DIM>& CircularBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::rGetCentreOfCircle() const
{
    return mCentreOfCircle;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CircularBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::GetRadiusOfCircle() const
{
    return mRadiusOfCircle;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CircularBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations)
{
    if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation)==nullptr)
    {
        EXCEPTION("PlaneBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
    }

    assert(SPACE_DIM == ELEMENT_DIM); // LCOV_EXCL_LINE
    assert(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation));


    if (SPACE_DIM != 1)
    {
        if (dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation))
        {
            {
                // Iterate over all nodes...
                unsigned num_nodes = this->mpCellPopulation->GetNumNodes();
                for (unsigned node_index = 0; node_index < num_nodes; node_index++)
                {
                    Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(node_index);
                    c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();

                    {
                    // Find the radial distance between this cell and the surface of the sphere
                        double radius = norm_2(node_location - mCentreOfCircle);
                        assert(radius != 0.0); //Can't project the centre to anywhere sensible

                        std::cout << "Line 60 has passed.\n";

                        // If the cell is too far from the surface of the sphere...
                        if (fabs(radius - mRadiusOfCircle) > mMaximumDistance)
                        {
                        // ...move the cell back onto the surface of the sphere
                            c_vector<double, SPACE_DIM> location_on_circle =
                                    mCentreOfCircle + mRadiusOfCircle*(node_location - mCentreOfCircle)/radius;

                            p_node->rGetModifiableLocation() = location_on_circle;
                        }
                    }
                }
            }
        }
    }
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool CircularBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

    // Iterate over the cell population
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        // Find the radial distance between this cell and the surface of the sphere
        c_vector<double,SPACE_DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
        double radius = norm_2(cell_location - mCentreOfCircle);

        // If the cell is too far from the surface of the sphere...
        if (fabs(radius - mRadiusOfCircle) > mMaximumDistance)
        {
            // ...then the boundary condition is not satisfied
            condition_satisfied = false;
            break;
        }
    }

    return condition_satisfied;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CircularBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CentreOfCircle>";
    for (unsigned index=0; index != SPACE_DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mCentreOfCircle[index] << ",";
    }
    *rParamsFile << mCentreOfCircle[SPACE_DIM-1] << "</CentreOfCircle>\n";

    *rParamsFile << "\t\t\t<RadiusOfCircle>" << mRadiusOfCircle << "</RadiusOfCircle>\n";
    *rParamsFile << "\t\t\t<MaximumDistance>" << mMaximumDistance << "</MaximumDistance>\n";
    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class CircularBoundaryCondition<1,1>;
template class CircularBoundaryCondition<1,2>;
template class CircularBoundaryCondition<2,2>;
template class CircularBoundaryCondition<1,3>;
template class CircularBoundaryCondition<2,3>;
template class CircularBoundaryCondition<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CircularBoundaryCondition)
