#ifndef TESTECLIPSEVERTEX5_HPP_
#define TESTECLIPSEVERTEX5_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractCellProperty.hpp"
#include "AbstractForce.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "OffLatticeSimulation.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PlaneBasedCellKiller.hpp"

#include "FakePetscSetup.hpp"

class MotileCellProperty : public AbstractCellProperty
{
private:
	unsigned mColour;
	friend class boost::serialization::access;
	void serialize (Archive & archive, const unisgned int version)
	{
		archive &boost::serialization::base_object<AbstractCellProperty>(*this);
		archive & mColour;
	}
public:
	MotileCellProperty(unsigned colour=5)
		: AbstractCellProperty(),
		  mColour(colour)
	{
	}

	~MotileCellProperty()
	{}

	unisnged GetColour () const
	{
		return mColour;
	}
};
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MotileCellProperty)
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MotileCellProperty)

class TesEclipseVertex4 : public AbstractCellBasedTestSuite
{
public:
    void TestPeriodicMonolayer()
    {
        HoneycombVertexMeshGenerator generator(16, 16);    // Parameters are: cells across, cells up
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

	MAKE_PTR(MotileCellProperty, p_motile);
	MAKE_PTR(CellLabel, p_lable);

	MAKE_PTR(WildTypeCellMutationState, p_state);
	MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
	std::vector<CellPtr> cells;
	for (unsigned i=0; i<mesh.GetNumElements() ; i++)
	{
		FixedG1GenerationalCellCylceModel*p_model = new FixedG1GenerationalCellCycleModel();

		CellPropertyCollection collection;
		if (RandomNumberGenerator::Instance()->ranf() < 0.2)
		{
			collection.AddProperty(p_motile);
			collection.AddProperty(p_label);
		}

		CellPtr p_cell (new Cell(p_state, p_model, NULL, false, collection));
		p_cell->SetCellProliferativeType(p_diff_type);

		double birth_time = - RandomNumberGenerator::Instance()->ranf() *
					(p_model->GetStemCellG1Duration()
					    + p_model->GetSG2MDuration());

		p_cell->SetBirthTime(birth_time);
		cells.push_back(p_cell);
	}


        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

	cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EclipseVertex5");
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(5.0);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

	MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

	c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        point(1) = 40.0;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, point, normal));
        simulator.AddCellKiller(p_killer);

        simulator.Solve();

    }
} ;

#endif /* TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_HPP_ */
