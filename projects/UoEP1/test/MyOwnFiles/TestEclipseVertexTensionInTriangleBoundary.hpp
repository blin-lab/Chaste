/*
 * TestEclipseVertexTensionInTriangleBoundary.hpp
 *
 *  Created on: 7 jul. 2020
 *      Author: mirob
 */
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

#include "CellLineTensionWriter.hpp"
#include "NagaiHondaCellTensionWriter.hpp"

/* This force law assumes that cells possess a "target area" property which determines the size of each
 * cell in the simulation. In order to assign target areas to cells and update them in each time step, we need
 * the next header file.*/

#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "VolumeTrackingModifier.hpp"

#include "PlaneBoundaryCondition.hpp"

#include "FakePetscSetup.hpp"

class TestEclipseVertexTensionInTriangleBoundary : public AbstractCellBasedTestSuite
{
public:
    void TestMonolayer()
    {
        HoneycombVertexMeshGenerator generator(1, 1);    // Parameters are: cells across, cells up
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<NagaiHondaCellTensionWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EclipseNagaiHondaTensionWriterUnevenTriangleTrying");
        simulator.SetEndTime(30.0);
        simulator.SetSamplingTimestepMultiple(100);



        /* We must now create one or more force laws, which determine the mechanics of the vertices
         * of each cell in a cell population. For this test, we use one force law, based on the
         * Nagai-Honda mechanics, and pass it to the {{{OffLatticeSimulation}}}.
         * For a list of possible forces see subclasses of {{{AbstractForce}}}.
         * These can be found in the inheritance diagram, here, [class:AbstractForce AbstractForce].
         * Note that some of these forces are not compatible with vertex-based simulations see the specific class documentation for details,
         * if you try to use an incompatible class then you will receive a warning.
         */
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        c_vector<double, 2> point = zero_vector<double>(2);
        c_vector<double, 2> normal = zero_vector<double>(2);

        //Make a boundary so that it creates a triangle with coordinates (6 ,6) , (-3, 0) and (3, 0). The area of the triangle is of 18.97u^2 so rounded to 19.
        // Impose a wall at y = 0 so that y < 0 is left out of boundary
        point(0) = -2.0;
        normal(0) = -1.0;
        point(1) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);

        point(0) = -2.0;
        normal(0) = -1.0;
        point(1) = 0.0;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        point(0) = 3.0;
        normal(0) = 1.0;
        point(1) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);

        point(0) = 3.0;
        normal(0) = 1.0;
        point(1) = 0.0;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc4);
        simulator.Solve();

        simulator.Solve();

    }
};
