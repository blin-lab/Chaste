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
#include "NagaiHondaCellTensionWriter.hpp"
#include "FakePetscSetup.hpp"

class TestEclipseVertexBoundary4 : public AbstractCellBasedTestSuite
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
        cell_population.AddCellWriter<NagaiHondaCellTensionWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EclipseNagaiTensionWriter");
        simulator.SetEndTime(45.0);

        simulator.SetSamplingTimestepMultiple(50);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        simulator.Solve();

    }
};
