#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"

#include "UniformG1GenerationalCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "FarhadifarForce.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include "PlaneBoundaryCondition.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "CellLabel.hpp"
#include "RandomNumberGenerator.hpp"
#include "CellLineTensionWriter.hpp"
#include "FakePetscSetup.hpp"

class TestEclipseVertexTensionInSquareBoundary : public AbstractCellBasedTestSuite
{
public:
    void TestMonolayer()
    {
        HoneycombVertexMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellLineTensionWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EclipseTensionWriterSquare");
        simulator.SetEndTime(200.0);

        simulator.SetSamplingTimestepMultiple(50);

        MAKE_PTR(FarhadifarForce<2>, p_force);
        p_force->SetBoundaryLineTensionParameter(0.12);
        simulator.AddForce(p_force);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        c_vector<double, 2> point = zero_vector<double>(2);
        c_vector<double, 2> normal = zero_vector<double>(2);

        point (0) = 0.0;
        normal (0) = -1.0;                             //creates a normal to the plane at x = 0 with vector direction -1.
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal)); //Vector -1 means that on the x plane it will go one unit towards the negative side.
        simulator.AddCellPopulationBoundaryCondition(p_bc1);

        point(0) = 0.0;
        normal(0) = 0.0;
        point (1) = 3.0;
        normal (1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal)); //Direction of the vector is one unit towards the positive side on the x plane.
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        point(0) = 3.0;
        normal(0) = 1.0;
        point(1) = 0.0;
        normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);

        normal (1) = -1.0;
        point (1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));

        simulator.Solve();

    }
};
