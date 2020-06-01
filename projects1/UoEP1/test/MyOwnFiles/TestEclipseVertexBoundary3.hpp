#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PlaneBasedCellKiller.hpp"

#include "FakePetscSetup.hpp"

class TestRunningVertexBasedSimulationsTutorial : public AbstractCellBasedTestSuite
{
public:
    void TestMonolayer()
    {
        HoneycombVertexMeshGenerator generator(6, 8);    // Parameters are: cells across, cells up
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EclipseVertexBoundary3");
        simulator.SetEndTime(5.0);

        simulator.SetSamplingTimestepMultiple(50);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

	c_vector<double,2> point = zero_vector<double>(2);
	c_vector<double,2> normal =zero_vector<double>(2);
	point (1) = 4.0;
	normal (1) = -1.0;
	MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
	simulator.AddCellPopulationBoundaryCondition(p_bc1);

	c_vector<double,2> point1 = zero_vector<double>(2);
	c_vector<double,2> normal1 = zero_vector<double>(2);
	normal1 (1) = -1.0;
	MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point1, normal1));
	simulator.AddCellPopulationBoundaryCondition(p_bc2);

	point (1) = 4.0;
	normal (1) = 1.0;
	point1 (1) = 1.0;
	normal1 (1) = -1.0;

	simulator.Solve();

    }
};
