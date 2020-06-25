
#ifndef CELLLINETENSIONWRITER_HPP_
#define CELLLINETENSIONWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"
#include "FarhadifarForce.hpp"


/**
 * A class written using the visitor pattern for writing cell ages to file.
 *
 * The output file is called cellages.dat by default. If VTK is switched on,
 * then the writer also specifies the VTK output for each cell, which is stored
 * in the VTK cell data "Ages" by default.
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellLineTensionWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    // Needed for serialization.
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mAreaElasticityParameter;
        archive & mPerimeterContractilityParameter;
        archive & mLineTensionParameter;
        archive & mBoundaryLineTensionParameter;
    }

protected:

    /**
     * The strength of the area term in the model. Corresponds to K_alpha in Farhadifar's paper.
     */
    double mAreaElasticityParameter;

    /**
     * The strength of the perimeter term in the model. Corresponds to Gamma_alpha in Farhadifar's paper.
     */
    double mPerimeterContractilityParameter;

    /**
     * The strength of the line tension term in the model. Lambda_{i,j} in Farhadifar's paper.
     */
    double mLineTensionParameter;

    /**
     * The strength of the line tension at the boundary. This term does correspond to Lambda_{i,j} in Farhadifar's paper.
     */
    double mBoundaryLineTensionParameter;


public:

    /**
     * Default constructor.
     */
    CellLineTensionWriter();

    /**
     * Overridden GetCellDataForVtkOutput() method.
     *
     * Get a double associated with a cell. This method reduces duplication
     * of code between the methods VisitCell() and AddVtkData().
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     *
     * @return data associated with the cell
     */
    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Overridden VisitCell() method.
     *
     * Visit a cell and write its age.
     *
     * Outputs a line of space-separated values of the form:
     * ...[location index] [x-pos] [y-pos] [z-pos] [cell age] ...
     * with [y-pos] and [z-pos] included for 2 and 3 dimensional simulations, respectively.
     *
     * This is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);



    void SetBoundaryLineTensionParameter(double BoundaryLineTensionParameter);

    void SetLineTensionParameter(double lineTensionParameter);

    double GetLineTensionParameter();

    double GetBoundaryLineTensionParameter();

    double GetLineTensionParameter(Node<SPACE_DIM> *pNodeA, Node<SPACE_DIM> *pNodeB,
                                   VertexBasedCellPopulation<SPACE_DIM> &rVertexCellPopulation);


};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellLineTensionWriter)

#endif /* CELLLINETENSIONWRITER_HPP_ */

