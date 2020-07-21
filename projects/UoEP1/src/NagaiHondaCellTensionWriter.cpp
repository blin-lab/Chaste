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

#include "NagaiHondaCellTensionWriter.hpp"

#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
NagaiHondaCellTensionWriter<ELEMENT_DIM, SPACE_DIM>::NagaiHondaCellTensionWriter() :
AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("celltension.dat"),
mNagaiHondaMembraneSurfaceEnergyParameter(10.0), // This is 0.1 in the Nagai & Honda paper.
mNagaiHondaCellCellAdhesionEnergyParameter(0.5), // This corresponds to a value of 1.0 for
// the sigma parameter in the Nagai & Honda
// paper. In the paper, the sigma value is
// set to 0.01.
mNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0) // This is 0.01 in the Nagai & Honda paper.
{
    this->mVtkCellDataName = "Cell Tension";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NagaiHondaCellTensionWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(
        CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> *pCellPopulation) //
{

    double total_tension = 0.0; // Declare the variable total_tension which will later be used to do the sum of all the line tensions.
    unsigned edge_count = 0;

    //We have an Abstract Cell Population and want to specify that what we will be using in this case is a Vertex Based cell population.
    //This conversion can be done using static cast.
    VertexBasedCellPopulation<SPACE_DIM> *p_cell_population =
            static_cast<VertexBasedCellPopulation<SPACE_DIM>*>(pCellPopulation); //

    unsigned num_elements = p_cell_population->GetNumElements();

    // Begin by computing the area and perimeter of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_areas(num_elements);
    std::vector<double> element_perimeters(num_elements);
    std::vector<double> target_areas(num_elements);
    //Iterate over the elements of the mesh in order to go one by one and calculating the area and perimeter of each.
    for (typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
            elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
            ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        element_areas[elem_index] = p_cell_population->rGetMesh().GetVolumeOfElement(elem_index);
        element_perimeters[elem_index] = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(elem_index);
        try
        {
            // If we haven't specified a growth modifier, there won't be any target areas in the CellData array and CellData
            // will throw an exception that it doesn't have "target area" entries.  We add this piece of code to give a more
            // understandable message. There is a slight chance that the exception is thrown although the error is not about the
            // target areas.
            target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");
        }
        catch (Exception&)
        {
            EXCEPTION("You need to add an AbstractTargetAreaModifier to the simulation in order to use NagaiHondaForce");
        }
    }

    //Declare the tension variable that will be used later.
    c_vector<double, SPACE_DIM> membrane_surface_tension = zero_vector<double>(SPACE_DIM);

    // Find the cell by: accessing members of p_cell_population through a pointer with argument pCell.
    // And find the element attached to that cell.
    VertexElement <ELEMENT_DIM, SPACE_DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(pCell);
    // Get this element, its index and its number of nodes
    unsigned elem_index = p_element->GetIndex();
    unsigned num_nodes_elem = p_element->GetNumNodes();

    // Iterate over vertices in the cell population
    //Is this the right way of doing it?
    for (unsigned local_index=0; local_index<num_nodes_elem; local_index++)
    {

        // Get the previous and next nodes in this element
        unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;

        // Compute the gradient of each these edges, computed at the present node
        c_vector<double, SPACE_DIM> previous_edge_gradient = -p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
        c_vector<double, SPACE_DIM> next_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);


        // Add the force contribution from this cell's membrane surface tension (note the minus sign)
        c_vector<double, SPACE_DIM> element_perimeter_gradient = previous_edge_gradient + next_edge_gradient;
        double cell_target_perimeter = 2*sqrt(M_PI*target_areas[elem_index]);
        membrane_surface_tension -= 2*GetNagaiHondaMembraneSurfaceEnergyParameter()*(element_perimeters[elem_index] - cell_target_perimeter)*element_perimeter_gradient;


        total_tension += sqrt((membrane_surface_tension[0] * membrane_surface_tension[0]) + (membrane_surface_tension[1] * membrane_surface_tension[1]));
        edge_count++;

    }

    return (total_tension);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NagaiHondaCellTensionWriter<ELEMENT_DIM, SPACE_DIM>::GetAdhesionParameter(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB, VertexBasedCellPopulation<SPACE_DIM>& rVertexCellPopulation)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    double adhesion_parameter = GetNagaiHondaCellCellAdhesionEnergyParameter();

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        adhesion_parameter = GetNagaiHondaCellBoundaryAdhesionEnergyParameter();
    }

    return adhesion_parameter;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NagaiHondaCellTensionWriter<ELEMENT_DIM, SPACE_DIM>::GetNagaiHondaMembraneSurfaceEnergyParameter()
{
    return mNagaiHondaMembraneSurfaceEnergyParameter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NagaiHondaCellTensionWriter<ELEMENT_DIM, SPACE_DIM>::GetNagaiHondaCellCellAdhesionEnergyParameter()
{
    return mNagaiHondaCellCellAdhesionEnergyParameter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double NagaiHondaCellTensionWriter<ELEMENT_DIM, SPACE_DIM>::GetNagaiHondaCellBoundaryAdhesionEnergyParameter()
{
    return mNagaiHondaCellBoundaryAdhesionEnergyParameter;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NagaiHondaCellTensionWriter<ELEMENT_DIM, SPACE_DIM>::SetNagaiHondaMembraneSurfaceEnergyParameter(double membraneSurfaceEnergyParameter)
{
    mNagaiHondaMembraneSurfaceEnergyParameter = membraneSurfaceEnergyParameter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NagaiHondaCellTensionWriter<ELEMENT_DIM, SPACE_DIM>::SetNagaiHondaCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mNagaiHondaCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NagaiHondaCellTensionWriter<ELEMENT_DIM, SPACE_DIM>::SetNagaiHondaCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mNagaiHondaCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NagaiHondaCellTensionWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(
        CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> *pCellPopulation)
{

    *this->mpOutStream << pCellPopulation->GetLocationIndexUsingCell(pCell) << " ";

    c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << cell_location[i] << " ";
    }


    double total_tension = 0.0; // Declare the variable total_tension which will later be used to do the sum of all the line tensions.
    unsigned edge_count = 0;

    //We have an Abstract Cell Population and want to specify that what we will be using in this case is a Vertex Based cell population.
    //This conversion can be done using static cast.
    VertexBasedCellPopulation<SPACE_DIM> *p_cell_population =
            static_cast<VertexBasedCellPopulation<SPACE_DIM>*>(pCellPopulation); //

    unsigned num_elements = p_cell_population->GetNumElements();

    // Begin by computing the area and perimeter of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_areas(num_elements);
    std::vector<double> element_perimeters(num_elements);
    std::vector<double> target_areas(num_elements);
    //Iterate over the elements of the mesh in order to go one by one and calculating the area and perimeter of each.
    for (typename VertexMesh<ELEMENT_DIM,SPACE_DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
            elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
            ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        element_areas[elem_index] = p_cell_population->rGetMesh().GetVolumeOfElement(elem_index);
        element_perimeters[elem_index] = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(elem_index);
        try
        {
            // If we haven't specified a growth modifier, there won't be any target areas in the CellData array and CellData
            // will throw an exception that it doesn't have "target area" entries.  We add this piece of code to give a more
            // understandable message. There is a slight chance that the exception is thrown although the error is not about the
            // target areas.
            target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");
        }
        catch (Exception&)
        {
            EXCEPTION("You need to add an AbstractTargetAreaModifier to the simulation in order to use NagaiHondaForce");
        }
    }

    //Declare the tension variable that will be used later.
    c_vector<double, SPACE_DIM> membrane_surface_tension = zero_vector<double>(SPACE_DIM);

    // Find the cell by: accessing members of p_cell_population through a pointer with argument pCell.
    // And find the element attached to that cell.
    VertexElement <ELEMENT_DIM, SPACE_DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(pCell);
    // Get this element, its index and its number of nodes
    unsigned elem_index = p_element->GetIndex();
    unsigned num_nodes_elem = p_element->GetNumNodes();

    // Iterate over vertices in the cell population
    //Is this the right way of doing it?
    for (unsigned local_index=0; local_index<num_nodes_elem; local_index++)
    {

        // Get the previous and next nodes in this element
        unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;

        // Compute the gradient of each these edges, computed at the present node
        c_vector<double, SPACE_DIM> previous_edge_gradient = -p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
        c_vector<double, SPACE_DIM> next_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);


        // Add the force contribution from this cell's membrane surface tension (note the minus sign)
        c_vector<double, SPACE_DIM> element_perimeter_gradient = previous_edge_gradient + next_edge_gradient;
        double cell_target_perimeter = 2*sqrt(M_PI*target_areas[elem_index]);
        membrane_surface_tension -= 2*GetNagaiHondaMembraneSurfaceEnergyParameter()*(element_perimeters[elem_index] - cell_target_perimeter)*element_perimeter_gradient;


        total_tension += sqrt((membrane_surface_tension[0] * membrane_surface_tension[0]) + (membrane_surface_tension[1] * membrane_surface_tension[1]));
        edge_count++;

    }

    *this->mpOutStream << total_tension << " ";
}

// Explicit instantiation
template class NagaiHondaCellTensionWriter<2, 2> ;


#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
//EXPORT_TEMPLATE_CLASS_2_INTERNAL(NagaiHondaCellTensionWriter)

