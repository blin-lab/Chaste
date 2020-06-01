#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "AbstractCellBasedSimulationModifier.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "FakePetscSetup.hpp"

class CellBoundaryFillModifier : public AbstractCellBasedSimulationModifier<2,2>
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize (Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<2,2> >(*this);
    }

public:

    CellBoundaryFillModifier()
        : AbstractCellBasedSimulationModifier<2,2>()
    {}

    ~CellBoundaryFillModifier()
    {}

    void UpdateAtEndOfTimeStep(AbstractCellPopulation<2,2>& rCellPopulation)
    {
        UpdateCellData(rCellPopulation);
    }

    void SetupSolve(AbstractCellPopulation<2,2>& rCellPopulation, std::string outputDirectory)
    {

        UpdateCellData(rCellPopulation) ;
    }

    void UpdateCellData(AbstractCellPopulation<2,2>& rCellPopulation)
    {
        rCellPopulation.Update();

        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin() ;
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            double cell_height = rCellPopulation.GetLocationOfCellCentre(*cell_iter)[1];

            cell_iter->GetCellData()->SetItem("height", cell_height);
        }
    }

    void OutputSimulationModifierParameters(out_stream& rParamsFile)
    {
        AbstractCellBasedSimulationModifier<2>::OutputSimulationModifierParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CellBoundaryFillModifier)
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CellBoundaryFillModifier)

class TestEclipseVertexBoundaryWithFillModifier : public AbstractCellBasedTestSuite
{
public:
    void TestMonolayer()
    {
        HoneycombVertexMeshGenerator generator(6, 8);    // Parameters are: cells across, cells up
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EclipseVertexBoundaryWithModifiers2");
        simulator.SetEndTime(12.0);

        simulator.SetSamplingTimestepMultiple(50);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(1.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        simulator.AddForce(p_force);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        point(0) = 10.0;
        normal(0) = 0.5;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);
        point(0) = -0.0;
        point(1) = -0.0;
        normal(0) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
        point(1) = 10.0;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        MAKE_PTR(CellBoundaryFillModifier, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        simulator.Solve();

    }
};
