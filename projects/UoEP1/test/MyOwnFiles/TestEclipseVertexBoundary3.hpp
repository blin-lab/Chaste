#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "AbstractCellPopulationBoundaryCondition.hpp"
#include "OffLatticeSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "SmartPointers.hpp"
#include "FakePetscSetup.hpp"

class MyBoundaryCondition : public AbstractCellPopulationBoundaryCondition<2>
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<2> >(*this);
    }

public:
    MyBoundaryCondition(AbstractCellPopulation<2>* pCellPopulation)
        : AbstractCellPopulationBoundaryCondition<2>(pCellPopulation)
    {
    }

    void ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
    {
        for (AbstractCellPopulation<2>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
            Node<2>* p_node = this->mpCellPopulation->GetNode(node_index);
            double y_coordinate = p_node->rGetLocation()[1];

            if (y_coordinate > 5.0)
            {
                p_node->rGetModifiableLocation()[1] = 5.0;
            }
            else if (y_coordinate < 0.0)
            {
                p_node->rGetModifiableLocation()[1] = 0.0;
            }
        }
    }

    bool VerifyBoundaryCondition()
    {
        bool condition_satisfied = true;

        for (AbstractCellPopulation<2>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            c_vector<double, 2> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
            double y_coordinate = cell_location(1);

            if ((y_coordinate < 0.0) || (y_coordinate > 5.0))
            {
                condition_satisfied = false;
                break;
            }
        }
        return condition_satisfied;
    }

    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
    {
        AbstractCellPopulationBoundaryCondition<2>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MyBoundaryCondition)
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MyBoundaryCondition)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const MyBoundaryCondition * t, const unsigned int file_version)
        {
            const AbstractCellPopulation<2>* const p_cell_population = t->GetCellPopulation();
            ar << p_cell_population;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, MyBoundaryCondition * t, const unsigned int file_version)
        {
            AbstractCellPopulation<2>* p_cell_population;
            ar >> p_cell_population;

            ::new(t)MyBoundaryCondition(p_cell_population);
        }
    }
}

class TestCreatingAndUsingANewCellPopulationBoundaryConditionTutorial : public AbstractCellBasedTestSuite
{
public:

    void TestMyBoundaryCondition()
    {
        HoneycombMeshGenerator generator(7, 7);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        MyBoundaryCondition bc(&cell_population);

        bool population_satisfies_bc = bc.VerifyBoundaryCondition();
        TS_ASSERT_EQUALS(population_satisfies_bc, false);

        std::map<Node<2>*, c_vector<double, 2> > old_node_locations;
        for (AbstractMesh<2,2>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
                node_iter != p_mesh->GetNodeIteratorEnd();
                ++node_iter)
        {
            old_node_locations[&(*node_iter)] = node_iter->rGetLocation();
        }

        bc.ImposeBoundaryCondition(old_node_locations);

        population_satisfies_bc = bc.VerifyBoundaryCondition();
        TS_ASSERT_EQUALS(population_satisfies_bc, true);

        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "my_bc.arch";
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            AbstractCellPopulationBoundaryCondition<2>* const p_bc = new MyBoundaryCondition(NULL);
            output_arch << p_bc;
            delete p_bc;
        }
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCellPopulationBoundaryCondition<2>* p_bc;
            input_arch >> p_bc;

            delete p_bc;
        }
    }

    void TestOffLatticeSimulationWithMyBoundaryCondition()
    {
        HoneycombMeshGenerator generator(7, 7, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        MAKE_PTR_ARGS(MyBoundaryCondition, p_bc, (&cell_population));

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithMyBoundaryCondition");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(1.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(3);
        simulator.AddForce(p_linear_force);

        simulator.AddCellPopulationBoundaryCondition(p_bc);

        simulator.Solve();
    }
}
