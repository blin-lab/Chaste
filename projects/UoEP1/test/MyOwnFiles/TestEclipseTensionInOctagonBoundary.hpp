/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
#include "CellVolumesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"

#include "PlaneStickyBoundaryCondition.hpp"

#include "FakePetscSetup.hpp"

class TestEclipseVertexTensionInSquareBoundary : public AbstractCellBasedTestSuite
{
public:
    void TestMonolayer() //throw (Exception)
    {
        HoneycombVertexMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
        MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<NagaiHondaCellTensionWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EclipseNagaiHondaTensionWriterOctagonActiveBoundary");
        simulator.SetEndTime(31.0);

        simulator.SetSamplingTimestepMultiple(50);

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

        /* A {{{NagaiHondaForce}}} assumes that each cell has a target area. The target areas of cells are used to determine pressure
         * forces on each vertex and eventually determine the size of each cell in the simulation. In order to assign target areas to cells
         * and update them in each time step we add a {{{SimpleTargetAreaModifier}}} to the simulation, which inherits from
         *  {{{AbstractTargetAreaModifier}}}.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        c_vector<double, 2> point = zero_vector<double>(2);
        c_vector<double, 2> normal = zero_vector<double>(2);

        /* Create an Octagon Boundary by stablishing a square boundary and
         * overlaying a square (rotated 90º) on top so that it cuts the square
         * resulting in a octagon boundary.
         */

        /* Square Point A lies at (0 , 0).
         * Bottom Boundary is parallel to x axis and going negative. So that the cells populate the positive part.
         */

        point(0) = 0.0;
        point(1) = 0.0;
        normal(0) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);

        /* Square Point B lies at (-1.5 , 0).
         * Left Boundary is parallel to y axis at x = -1.5 always. And going negative.
         */

        point(0) = -1.5;
        point(1) = 0.0;
        normal(0) = -1.0;
        normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        /* Square Point C lies at (0 , 5).
         * Top Boundary is parallel y = 5.0 always. And it goes upwards so positive. Cell populate the down part of boundary.
         */

        point(0) = 0.0;
        point(1) = 5.0;
        normal(0) = 0.0;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);

        /* Square Point D lies at (3.6 , 0).
         * Right Boundary is parallel to y axis at x = -1.5 always. Cells populate to left of boundary.
         */

        point (0) = 3.6;
        point (1) = 0.0;
        normal (0) = 1.0;
        normal (1) = 0.0;
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        /* Diamond Point A lies at (-1.5 , 1.5).
         * Bottom left boundary of octagon
         */

        point(0) = -1.5;
        point(1) = 1.5;
        normal(0) = -1.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_bc5, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc5);

        /* Diamond Point B lies at (-1.5 , 3.6).
         * Top-Left boundary of octagon.
         */

        point(0) = -1.5;
        point(1) = 3.6;
        normal(0) = -1.0;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_bc6, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc6);

        /* Diamond Point C lies at (3.6 , 3.5).
         * Top-Right boundary of octagon.
         */

        point(0) = 3.6;
        point(1) = 3.5;
        normal(0) = 1.0;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_bc7, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc7);

        /* Diamond Point B lies at (3.6, 1.5).
         * Bottom-Right boundary of octagon.
         */

        point(0) = 3.6;
        point(1) = 1.5;
        normal(0) = 1.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneStickyBoundaryCondition<2>, p_bc8, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc8);

        simulator.Solve();

    }
};

