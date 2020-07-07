#ifndef CIRCULARBOUNDARY_HPP_
#define CIRCULARBOUNDARY_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A spherical cell population boundary condition class, which restricts nodes to lie
 * on the surface of a Circular in the domain. Although the name of this class suggests
 * it is specific to 3D, it is actually also implemented in 2D, for which it is really
 * a circle geometry boundary condition.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class CircularBoundaryCondition : public AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>
{
private:

    /** The centre of the circle. */
    c_vector<double, SPACE_DIM> mCentreOfCircle;

    /** The radius of the circle. */
    double mRadiusOfCircle;

    /** The maximum distance from the surface of the circle that cells may be. */
    double mMaximumDistance;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mMaximumDistance;
        //archive & mUseJiggledNodesOnPlane;
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param centre the centre of the circle
     * @param radius the radius of the circle
     * @param distance the maximum distance from the surface of the circle that cells may be (defaults to 1e-5)
     */
    CircularBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* pCellPopulation,
                                    c_vector<double, SPACE_DIM> centre,
                                    double radius,
                                    double distance=1e-5);

    /**
     * @return #mCentreOfCircle.
     */
    const c_vector<double, SPACE_DIM>& rGetCentreOfCircle() const;

    /**
     * @return #mRadiusOfCircle.
     */
    double GetRadiusOfCircle() const;

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CircularBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CircularBoundaryCondition.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const CircularBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;

    // Archive c_vectors one component at a time
    c_vector<double, SPACE_DIM> point = t->rGetCentreOfCircle();
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        ar << point[i];
    }

    // Archive other member variables
    double radius = t->GetRadiusOfCircle();
    ar << radius;
}

/**
 * De-serialize constructor parameters and initialize a CircularBoundaryCondition.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, CircularBoundaryCondition<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* p_cell_population;
    ar >> p_cell_population;

    // Retrieve c_vectors one component at a time
    c_vector<double, SPACE_DIM> point;
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        ar >> point[i];
    }

    // Retrieve other member variables
    double radius;
    ar >> radius;

    // Invoke inplace constructor to initialise instance
    ::new(t)CircularBoundaryCondition<ELEMENT_DIM, SPACE_DIM>(p_cell_population, point, radius);
}
}
} // namespace ...

#endif /*CIRCULARBOUNDARYCONDITION_HPP_*/
