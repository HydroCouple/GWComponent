/*!
*  \file    gwmodel.h
*  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
*  \version 1.0.0
*  \section Description
*  This file and its associated files and libraries are free software;
*  you can redistribute it and/or modify it under the terms of the
*  Lesser GNU Lesser General Public License as published by the Free Software Foundation;
*  either version 3 of the License, or (at your option) any later version.
*  fvhmcompopnent.h its associated files is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
*  \date 2018
*  \pre
*  \bug
*  \todo
*  \warning
*/

#ifndef GWMODEL_H
#define GWMODEL_H

#include "gwcomponent.h"
#include "odesolver.h"

#include <QObject>
#include <QFileInfo>
#include <unordered_map>
#include <vector>

#ifdef USE_NETCDF
#include "threadsafenetcdf/threadsafencvar.h"
#endif

class ThreadSafeNcFile;
struct Element;
struct ElementCell;
struct ElementJunction;
class IBoundaryCondition;
class TimeSeries;

struct SolverUserData
{
    GWModel *model = nullptr;
    int variableType = 0;
    int variableIndex = -1;
};

typedef void (*RetrieveCouplingData)(GWModel *model, double dateTime);
typedef void (*WriteVariableToNetCDF)(size_t currentTime,
                                      ThreadSafeNcVar &variable,
                                      const std::vector<Element*>& elements);

class GWCOMPONENT_EXPORT GWModel: public QObject
{
    Q_OBJECT

    friend struct Element;
    friend struct ElementCell;
    friend class EdgeBC;
    friend class ElementCellSourceBC;
    friend class ChannelBC;

  public:


    enum AdvectionDiscretizationMode
    {

      /*!
       * \brief Upwind discretization scheme.
       */
      Upwind,

      /*!
       * \brief Central discretization scheme.
       */
      Central,

      /*!
       * \brief Spalding, D.B., 1972. A Novel Finite Difference Formulation for Differential
       * Expressions Involving Both First and Second Derivatives.
       * International Journal for Numerical Methods in Engineering 4:551–559.
       */
      Hybrid,

      /*!
       * \brief TVD schemes
       */
      TVD,
    };

    /*!
     * \brief STMModel - Constructor for the Computational engine for the Stream Temperature Model.
     */
    GWModel(GWComponent *component);

    /*!
     * \brief ~STMModel - Destructor for the Computational engine for the Stream Temperature Model.
     */
    ~GWModel();

    /*!
     * \brief minTimeStep - Minimum timestep for the model in seconds.
     * \return Returns the minimum timestep for the current model instance.
     */
    double minTimeStep() const;

    /*!
     * \brief setMinTimeStep - Sets the minimum timestep for the current model instance
     * \param timeStep - Minimum timestep in seconds.
     */
    void setMinTimeStep(double timeStep);

    /*!
     * \brief maxTimeStep
     * \return
     */
    double maxTimeStep() const;

    /*!
     * \brief setMaxTimeStep
     * \param timeStep
     */
    void setMaxTimeStep(double timeStep);

    /*!
     * \brief useAdaptiveTimeStep
     * \return
     */
    bool useAdaptiveTimeStep() const;

    /*!
     * \brief setUseAdaptiveTimeStep
     * \param use
     */
    void setUseAdaptiveTimeStep(bool use);

    /*!
     * \brief timeStepRelaxationFactor
     * \return
     */
    double timeStepRelaxationFactor() const;

    /*!
     * \brief setTimeStepRelaxationFactor - Adaptive time-step relaxation factor
     * \param tStepRelaxFactor
     */
    void setTimeStepRelaxationFactor(double tStepRelaxFactor);

    /*!
     * \brief currentTimeStep - Current time step in seconds
     * \return
     */
    double currentTimeStep() const;

    /*!
     * \brief startDateTime - Start date and time specified as a modified julian date
     * \return
     */
    double startDateTime() const;

    /*!
     * \brief setStartDate - Sets the value of the start date and time
     * \param dateTime - Start date and time specified as modified julian date
     */
    void setStartDateTime(double dateTime);

    /*!
     * \brief endDateTime - End datetime specified as
     * \return
     */
    double endDateTime() const;

    /*!
     * \brief setEndDateTime Sets the value of the end datetime
     * \param dateTime
     */
    void setEndDateTime(double dateTime);

    /*!
     * \brief outputInterval
     * \return
     */
    double outputInterval() const;

    /*!
     * \brief setOutputInterval
     * \param interval
     */
    void setOutputInterval(double interval);

    /*!
     * \brief currentDateTime
     * \return
     */
    double currentDateTime() const;

    /*!
     * \brief hydHeadSolver
     * \return
     */
    ODESolver *odeSolver() const ;

    /*!
     * \brief waterDensity
     * \return
     */
    double waterDensity() const;

    /*!
     * \brief setWaterDensity
     * \param value
     */
    void setWaterDensity(double value);

    /*!
     * \brief defaultSedimentDensity
     * \return
     */
    double defaultSedimentDensity() const;

    /*!
     * \brief setDefaultSedimentDensity
     * \param value
     */
    void setDefaultSedimentDensity(double value);

    /*!
     * \brief specificHeatCapacityWater
     * \return
     */
    double specificHeatCapacityWater() const;

    /*!
     * \brief setSpecificHeatCapacityWater
     * \param value
     */
    void setSpecificHeatCapacityWater(double value);

    /*!
     * \brief defaultSedimentSpecificHeatCapacity
     * \return
     */
    double defaultSedimentSpecificHeatCapacity() const;

    /*!
     * \brief setDefaultSedimentSpecificHeatCapacity
     * \param value
     */
    void setDefaultSedimentSpecificHeatCapacity(double value);

    /*!
     * \brief defaultPorosity
     * \return
     */
    double defaultPorosity() const;

    /*!
     * \brief setDefaultPorosity
     * \param porosity
     */
    void setDefaultPorosity(double porosity);

    /*!
     * \brief hydraulicConductivityX
     * \return
     */
    double defaultHydraulicConductivityX() const;

    /*!
     * \brief setHydraulicConductivityX
     * \param hydConX
     */
    void setDefaultHydraulicConductivityX(double hydConX);

    /*!
     * \brief hydraulicConductivityY
     * \return
     */
    double defaultHydraulicConductivityY() const;

    /*!
     * \brief setHydraulicConductivityY
     * \param hydConY
     */
    void setDefaultHydraulicConductivityY(double hydConY);

    /*!
     * \brief numSolutes
     * \return
     */
    int numSolutes() const;

    /*!
     * \brief setNumSolutes
     * \param numSolutes
     */
    void  setNumSolutes(int numSolutes);

    /*!
     * \brief setSoluteNuame
     * \param soluteIndex
     * \param soluteName
     */
    void setSoluteName(int soluteIndex, const std::string &soluteName);

    /*!
     * \brief solute
     * \param soluteIndex
     * \return
     */
    std::string solute(int soluteIndex) const;

    /*!
     * \brief simulateWaterAge
     * \return
     */
    bool simulateWaterAge() const;

    /*!
     * \brief setSimulateWaterAge
     * \param simulate
     */
    void setSimulateWaterAge(bool simulate);

    /*!
     * \brief verbose
     * \return
     */
    bool verbose() const;

    /*!
     * \brief setVerbose
     * \param verbose
     */
    void setVerbose(bool verbose);

    /*!
     * \brief printFrequency
     * \return
     */
    int printFrequency() const;

    /*!
     * \brief setPrintFrequency
     * \param printFreq
     */
    void setPrintFrequency(int printFreq);

    /*!
     * \brief flushToDiskFrequency
     * \return
     */
    int flushToDiskFrequency() const;

    /*!
     * \brief setFlushToDiskFrequency
     * \param diskFlushFrequency
     */
    void setFlushToDiskFrequency(int diskFlushFrequency);

    /*!
     * \brief numElementJunctions
     * \return
     */
    int numElementJunctions() const;

    /*!
     * \brief addControlVolumeNode
     * \param id
     * \param x
     * \param y
     * \param z
     * \return
     */
    ElementJunction *addElementJunction(const std::string &id, double x = 0, double y = 0, double z = 0);

    /*!
     * \brief removeElementJunction
     * \param id
     */
    void deleteElementJunction(const std::string &id);

    /*!
     * \brief removeElementJunction
     * \param id
     */
    void deleteElementJunction(int index);

    /*!
     * \brief getElementJunction
     * \param id
     * \return
     */
    ElementJunction *getElementJunction(const std::string &id) ;

    /*!
     * \brief getElementJunction
     * \param index
     * \return
     */
    ElementJunction *getElementJunction(int index) ;

    /*!
     * \brief numCellsPerElement
     * \return
     */
    int numLeftElementCells() const;

    /*!
     * \brief setNumCellsPerElement
     * \param value
     */
    void setNumLeftElementCells(int value);

    /*!
     * \brief numRightElementCells
     * \return
     */
    int numRightElementCells() const;

    /*!
     * \brief setNumRightElementCells
     * \param value
     */
    void setNumRightElementCells(int value);

    /*!
     * \brief numElementCells
     * \return
     */
    int numElementCells() const;

    /*!
     * \brief numElements
     * \return
     */
    int numElements() const;

    /*!
     * \brief addElement
     * \param id
     * \param fromElement
     * \param toElement
     * \return
     */
    Element *addElement(const std::string &id, ElementJunction *upStream, ElementJunction *downStream);

    /*!
     * \brief removeElement
     * \param id
     */
    void deleteElement(const std::string &id);

    /*!
     * \brief removeElement
     * \param index
     */
    void deleteElement(int index);

    /*!
     * \brief getElement
     * \param id
     * \return
     */
    Element *getElement(const std::string &id);

    /*!
     * \brief getElement
     * \param index
     * \return
     */
    Element *getElement(int index);

    /*!
     * \brief inputFile
     * \return
     */
    QFileInfo inputFile() const;

    /*!
     * \brief setInputFile
     * \param inputFile
     */
    void setInputFile(const QFileInfo &inputFile);

    /*!
     * \brief outputNetCDFFile
     * \return
     */
    QFileInfo outputNetCDFFile() const;

    /*!
     * \brief setOutputNetCDFFile
     * \param outputNetCDFFile
     */
    void setOutputNetCDFFile(const QFileInfo &outputNetCDFFile);

    /*!
     * \brief retrieveCouplingDataFunction
     * \return
     */
    RetrieveCouplingData retrieveCouplingDataFunction() const;

    /*!
     * \brief setRetrieveCouplingDataFunction
     */
    void setRetrieveCouplingDataFunction(RetrieveCouplingData retrieveCouplingDataFunction);

    /*!
     * \brief initialize
     * \param errors
     * \return
     */
    bool initialize(std::list<std::string> &errors);

    /*!
     * \brief update
     */
    void update();

    /*!
     * \brief finalize
     * \param errors
     * \return
     */
    bool finalize(std::list<std::string> &errors);

    /*!
     * \brief printStatus
     */
    void printStatus();

    /*!
     * \brief saveSave
     * \param filePath
     */
    void saveSave(const QFileInfo &filePath);

    /*!
     * \brief computeDerivedHydraulics
     */
    void calculatePreComputedHydHeadVariables();

    /*!
     * \brief calculatePreComputedTempVariables
     */
    void calculatePreComputedTempVariables();

    /*!
     * \brief calculatePreComputedSoluteVariables
     */
    void calculatePreComputedSoluteVariables();

  private:

    /*!
     * \brief initializeInputFiles
     * \param errors
     * \return
     */
    bool initializeInputFiles(std::list<std::string> &errors);

    /*!
     * \brief initializeTimeVariables
     * \param errors
     * \return
     */
    bool initializeTimeVariables(std::list<std::string> &errors);

    /*!
     * \brief initializeElements
     * \param errors
     * \return
     */
    bool initializeElements(std::list<std::string> &errors);

    /*!
     * \brief initializeSolver
     * \param errors
     * \return
     */
    bool initializeSolvers(std::list<std::string> &errors);

    /*!
     * \brief intializeOutputFiles
     * \param errors
     * \return
     */
    bool initializeOutputFiles(std::list<std::string> &errors);

    /*!
     * \brief initializeNetCDFOutputFile
     * \param errors
     * \return
     */
    bool initializeNetCDFOutputFile(std::list<std::string> &errors);

    /*!
     * \brief initializeBoundaryConditions
     * \param errors
     * \return
     */
    bool initializeBoundaryConditions(std::list<std::string> &errors);

    /*!
     * \brief prepareForNextTimeStep
     */
    void prepareForNextTimeStep();

    /*!
     * \brief applyInitialConditions
     */
    void applyInitialConditions();

    /*!
     * \brief applyBoundaryConditions
     */
    void applyBoundaryConditions(double dateTime);

    /*!
     * \brief computeTimeStep
     * \return
     */
    double computeTimeStep();

    /*!
     * \brief solveHydHead
     * \param timeStep
     */
    void solve(double timeStep);

    /*!
     * \brief computeDTDt
     * \param model
     * \param variableIndex
     * \param t
     * \param y
     * \param dydt
     */
    static void computeDYDt(double t, double y[], double dydt[], void *userData);

    /*!
     * \brief readInputFileOptionTag
     * \param line
     */
    bool readInputFileOptionTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileOutputTag
     * \param line
     */
    bool readInputFileOutputTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileSolutesTag
     * \param line
     */
    bool readInputFileSolutesTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileElementJunctionsTag
     * \param line
     */
    bool readInputFileElementJunctionsTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileElementsTag
     * \param line
     */
    bool readInputFileElementsTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileElementCellWidths
     * \param line
     * \param errorMessage
     * \return
     */
    bool readInputFileElementCellWidths(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileElementCellBottomElevTag
     * \param line
     */
    bool readInputFileElementCellHydraulics(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileElementCellBCHydHeadTag
     * \param line
     * \param errorMessage
     * \return
     */
    bool readInputFileElementCellInitConditions(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileElementCellBCHydHeadDerivTag
     * \param line
     * \param errorMessage
     * \return
     */
    bool readInputFileElementCellBoundaryConditions(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileElementCellBCFluxTag
     * \param line
     * \param errorMessage
     * \return
     */
    bool readInputFileElementCellSoluteBoundaryConditions(int variable, const QStringList &options, QString &errorMessage);

    /*!
     * \brief readInputFileElementCellSources
     * \param line
     * \param errorMessage
     * \return
     */
    bool readInputFileElementCellSources(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileTimeSeries
     * \param options
     * \param errorMessage
     * \return
     */
    bool readInputFileTimeSeriesTag(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileChannelBoundaryConditions
     * \param line
     * \param errorMessage
     * \return
     */
    bool readInputFileChannelBoundaryConditions(const QString &line, QString &errorMessage);

    /*!
     * \brief readInputFileChannelZCellFactors
     * \param line
     * \param errorMessage
     * \return
     */
    bool readInputFileChannelZCellFactors(const QString &line, QString &errorMessage);

    /*!
     * \brief writeOutput
     */
    void writeOutput();

    /*!
     * \brief writeNetCDFOutput
     */
    void writeNetCDFOutput();

    /*!
     * \brief closeOutputFiles
     */
    void closeOutputFiles();

    /*!
     * \brief closeOutputFiles
     */
    void closeOutputNetCDFFile();

    /*!
     * \brief relativePathtoAbsolute
     * \param inputFile
     * \return
     */
    QFileInfo relativePathToAbsolute(const QFileInfo& fileInfo);

    /*!
     * \brief findProfile
     * \param from
     * \param to
     * \param m_profile
     * \return
     */
    bool findProfile(Element *from, Element *to, std::vector<Element*> &profile);

    /*!
     * \brief calculateDistanceFromUpstreamJunction
     * \param Element
     */
    void calculateDistanceFromUpstreamJunction(Element *element);

    /*!
     * \brief breadthFirstSearchSetIndex
     * \param cell
     * \param bfsIndex
     */
    void breadthFirstSearchSetIndex(ElementCell *cell, int &bfsIndex);

  private:

    std::vector<std::string> m_solutes; // Names of the solutes.

    //Time variables
    double m_timeStep = 1, //seconds
    m_prevTimeStep = 1,
    m_startDateTime = 0, //Modified Julian Day
    m_endDateTime = 1, //Modified Julian Day
    m_currentDateTime = 0, //Modified Julian Day
    m_prevDateTime = 0, //Modified Julian Day
    m_maxTimeStep = 1.0, //seconds
    m_minTimeStep = 0.1, //seconds
    m_outputInterval, //seconds
    m_nextOutputTime,//Julian Day
    m_timeStepRelaxationFactor,
    m_maxHead,
    m_minHead,
    m_maxTemp, //Tracks maximum temperature so far
    m_minTemp; //Tracks minimum temperature so far

    std::vector<double> m_maxSolute, //array for tracking maximum solute concentrations
    m_minSolute, //array for tracking minimum solute concentrations
    m_totalSoluteMassBalance, // Tracks total mass balance of solutes (kg)
    m_totalExternalSoluteMassBalance, //Tracks total mass balance from external sources (kg)
    m_solverCurrentValues,
    m_solverOutValues,
    m_currTemps,
    m_outTemps,
    m_solute_first_order_k,
    m_solute_kd,
    m_solute_molecular_diff;

    const static float dir[];

    int m_tempIndex;
    std::vector<int> m_soluteIndexes;

    std::unordered_map<std::string, QSharedPointer<TimeSeries>> m_timeSeries;

    int m_numInitFixedTimeSteps, //Number of initial fixed timeSteps of the minimum timestep to use when using the adaptive time step;
    m_numCurrentInitFixedTimeSteps, //Count number of initial minimum timesteps that have been used
    m_printFrequency, //Number of timesteps before printing
    m_currentPrintCount, // Number of timesteps
    m_flushToDiskFrequency, // Number of times to write output stored in memory to disk
    m_currentflushToDiskCount, //Number of timesteps that have been stored in memory so far since the last flush to disk
    m_addedSoluteCount,
    m_numLeftCellsPerElement,
    m_numRightCellsPerElement,
    m_totalCellsPerElement,
    m_numSolutes = 0,
    m_numBedZCells = 1;

    bool m_useAdaptiveTimeStep, //Use the adaptive time step option
    m_verbose, //Print simulation information to console
    m_solveHeatTransport,
    m_flushToDisk,
    m_simulateWaterAge = false;

    //Element junctions
    std::vector<ElementJunction*> m_elementJunctions;
    std::unordered_map<std::string, ElementJunction*> m_elementJunctionsById; //added for fast lookup using identifiers instead of indexes.

    //1D Computational elements
    std::vector<double> m_zcellFactors;
    std::vector<Element*> m_elements;
    std::unordered_map<std::string, Element*> m_elementsById; //added for fast lookup using identifiers instead of indexes.

    //Boundary conditions list
    std::vector<IBoundaryCondition*> m_boundaryConditions;

    //Solver Objects
    ODESolver *m_odeSolver;

    //Global water properties
    double m_waterDensity = 1000, //kg/m^3
    m_cp = 4184,// Water specific heat capacity default = 4187.0; // J/kg/C,
    m_sedDensity = 2650, //Sediment density //kg/m^3
    m_sedCp = 880, //Sediment specific heat capacity// J/kg/C,
    m_hydConX = 0.0, //Hydraulic conductivity in the X direction
    m_hydConY = 0.0001, //Hydraaulic conductivity in the Y direction
    m_hydConZ = 0.0001,
    m_specificYield = 0.2,
    m_specificStorage = 6.2e-3,
    m_porosity = 0.3,
    m_defaultCellWidth = 1.0,
    m_totalMassBalance = 0.0,
    m_totalHeatBalance = 0.0, //Track total heat accumulation from external heat sources (KJ)
    m_totalExternalHeatBalance = 0.0,
    m_dispersivityX = 0.05,
    m_dispersivityY = 0.05,
    m_dispersivityZ = 0.05,
    m_waterThermalConductivity = 0.606,
    m_sedThermalConductivity = 2.6;

    /*!
     * \brief m_advectionMode
     */
    AdvectionDiscretizationMode m_advectionMode;
    int m_TVDFluxLimiter;

    //File input and output
    QFileInfo m_inputFile, //Input filepath
    m_outputNetCDFFileInfo; //Output NetCDF filepath

#ifdef USE_NETCDF
    ThreadSafeNcFile *m_outputNetCDF = nullptr; //NetCDF output file object
    std::unordered_map<std::string, ThreadSafeNcVar> m_outNetCDFVariables;
    std::unordered_map<std::string, bool> m_outNetCDFVariablesOnOff;
    std::unordered_map<std::string, WriteVariableToNetCDF> m_outNetCDFVariablesIOFunctions;
    std::vector<std::string> m_optionalOutputVariables;
#endif

    static const std::unordered_map<std::string, int> m_inputFileFlags; //Input file flags
    static const std::unordered_map<std::string, int> m_hydraulicVariableFlags; //Hydraulic variable flags
    static const std::unordered_map<std::string, int> m_optionsFlags; //Input file flags
    static const std::unordered_map<std::string, int> m_solverTypeFlags; //Solver type flags
    static const std::unordered_map<std::string, int> m_linearSolverTypeFlags; //Solver type flags
    static const std::unordered_map<std::string, int> m_advectionFlags; //Advection type flags
    static const std::unordered_map<std::string, int> m_BCFlags; //Advection type flags
    static const std::unordered_map<std::string, int> m_edgeFlags; //Advection type flags
    static const std::unordered_map<std::string, int> m_srcFlags; //Advection type flags
    static const std::unordered_map<std::string, int> m_geomMultFlags; //Advection type flags
    static const std::unordered_map<std::string, int> m_chanVarFlags; //Advection type flags

    static const QRegExp m_dateTimeDelim;

    QRegExp m_delimiters; //Regex delimiter

    //Function to retrieve and apply coupling data
    RetrieveCouplingData m_retrieveCouplingDataFunction;

    //Parent component
    GWComponent *m_component;
};

#endif // GWMODEL_H
