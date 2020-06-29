
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"

#include "UniformG1GenerationalCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "FarhadifarForce.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include "PlaneBasedCellKiller.hpp"
#include "CellLabel.hpp"
#include "CircularBoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"
#include "CellLabelWriter.hpp"
#include "CellLineTensionWriter.hpp"
#include "FakePetscSetup.hpp"
#include "Debug.hpp"


class TestMyOwnBoundary1 : public AbstractCellBasedTestSuite
{
public:
    void TestMonolayer()
    {
        HoneycombVertexMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellLineTensionWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EclipseCircularTension");
        simulator.SetEndTime(35.0);

        simulator.SetSamplingTimestepMultiple(50);

        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        c_vector<double,2> centre = zero_vector<double>(2);
        centre (0) = 2.0;
        double radius = 10.0;

        MAKE_PTR_ARGS(CircularBoundaryCondition<2>, p_bc, (&cell_population, centre, radius));
        simulator.AddCellPopulationBoundaryCondition(p_bc);


        simulator.Solve();

    }
};
