#!/usr/bin/env python

#DOC-START imports
import sys, os, exfile
#DOC-END imports

#DOC-START load exfile
#Load the mesh information in the form of exregion format
exregion = exfile.Exregion("hetrogenouscylinder.exregion")
#DOC-END load exfile

#Additional information regading the mesh should be provided, the rudimentary exfile object does not provide these details
#The mesh is composed of quad with tri-quadratic basis

#Set the dimension of the mesh, here each of its elements are 3D
numberOfXi = 3

# Intialise OpenCMISS
from opencmiss.iron import iron

# Set problem parameters
#Use pressure to enforce incompressibililty constraint
usePressureBasis = True
#Atleast 3 Gauss quadrature points along an xi direction is required for finite elasticity problems
numberOfGaussXi = 3

#Setup field number handles
coordinateSystemUserNumber = 1
regionUserNumber = 1
basisUserNumber = 1
pressureBasisUserNumber = 2
meshUserNumber = 1
cellMLUserNumber = 1
decompositionUserNumber = 1
equationsSetUserNumber = 1
problemUserNumber = 1
#Mesh component numbers
quadraticMeshComponentNumber = 1
linearMeshComponentNumber = 2
#Fields
geometricFieldUserNumber = 1
fibreFieldUserNumber = 2
dependentFieldUserNumber = 3
cellMLModelsFieldUserNumber = 4
cellMLParametersFieldUserNumber = 5
cellMLIntermediateFieldUserNumber = 6
equationsSetFieldUserNumber = 7

# Set all diganostic levels on for testing
#iron.DiagnosticsSetOn(iron.DiagnosticTypes.ALL,[1,2,3,4,5],"Diagnostics",["DOMAIN_MAPPINGS_LOCAL_FROM_GLOBAL_CALCULATE"])

if(usePressureBasis):
    numberOfMeshComponents = 2
else:
    numberOfMeshComponents = 1
    
# Set the OpenCMISS random seed so that we can test this example by using the
# same parallel decomposition
numberOfRandomSeeds = iron.RandomSeedsSizeGet()
randomSeeds = [0]*numberOfRandomSeeds
randomSeeds[0] = 100
iron.RandomSeedsSet(randomSeeds)

# Get the number of computational nodes and this computational node number
computationEnvironment = iron.ComputationEnvironment()
numberOfComputationalNodes = computationEnvironment.NumberOfWorldNodesGet()
computationalNodeNumber = computationEnvironment.WorldNodeNumberGet()

if computationalNodeNumber == 0:
    if not os.path.exists("./results"):
        os.makedirs("./results")
    
# Create a 3D rectangular cartesian coordinate system
coordinateSystem = iron.CoordinateSystem()
coordinateSystem.CreateStart(coordinateSystemUserNumber)
coordinateSystem.DimensionSet(3)
coordinateSystem.CreateFinish()

# Create a region and assign the coordinate system to the region
region = iron.Region()
region.CreateStart(regionUserNumber,iron.WorldRegion)
region.LabelSet("Region")
region.coordinateSystem = coordinateSystem
region.CreateFinish()

# Define basis
basis = iron.Basis()
basis.CreateStart(basisUserNumber)
basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
basis.NumberOfXiSet(numberOfXi)
basis.InterpolationXiSet([iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*numberOfXi)
if(numberOfGaussXi>0):
    basis.QuadratureNumberOfGaussXiSet([numberOfGaussXi]*numberOfXi)
basis.CreateFinish()

if(usePressureBasis):
    # Define pressure basis
    pressureBasis = iron.Basis()
    pressureBasis.CreateStart(pressureBasisUserNumber)
    pressureBasis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
    pressureBasis.NumberOfXiSet(numberOfXi)
    pressureBasis.InterpolationXiSet([iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfXi)
    if(numberOfGaussXi>0):
        pressureBasis.QuadratureNumberOfGaussXiSet([numberOfGaussXi]*numberOfXi)
    pressureBasis.CreateFinish()

# Start the creation of input mesh in the region
mesh = iron.Mesh()
mesh.CreateStart(meshUserNumber, region, numberOfXi)
mesh.NumberOfComponentsSet(numberOfMeshComponents)
mesh.NumberOfElementsSet(exregion.num_elements)

# Define nodes for the mesh
nodes = iron.Nodes()
nodes.CreateStart(region, exregion.num_nodes)
nodes.CreateFinish()

#Specify the elementwise topology
quadraticElemnts = iron.MeshElements()
quadraticElemnts.CreateStart(mesh, quadraticMeshComponentNumber, basis)
for elem in exregion.elements:
    quadraticElemnts.NodesSet(elem.number, elem.nodes)
quadraticElemnts.CreateFinish()

if(usePressureBasis):
    #Specify the elementwise topology that suits the mesh element basis
    linearElements = iron.MeshElements()
    linearElements.CreateStart(mesh, linearMeshComponentNumber, pressureBasis)
    for elem in exregion.elements:
        linearNodes = list(elem.nodes[i - 1] for i in [1, 3, 7, 9, 19, 21, 25, 27])
        linearElements.NodesSet(elem.number, linearNodes)
    linearElements.CreateFinish()
mesh.CreateFinish()


# Create a decomposition for the mesh
decomposition = iron.Decomposition()
decomposition.CreateStart(decompositionUserNumber,mesh)
decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
decomposition.NumberOfDomainsSet(numberOfComputationalNodes)
decomposition.CreateFinish()

# Create a field for the geometry
geometricField = iron.Field()
geometricField.CreateStart(geometricFieldUserNumber,region)
geometricField.MeshDecompositionSet(decomposition)
geometricField.TypeSet(iron.FieldTypes.GEOMETRIC)
geometricField.VariableLabelSet(iron.FieldVariableTypes.U,"coordinates")
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
geometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
geometricField.CreateFinish()

#Set the geometric information from the exregion file
#Here we will record special nodes of interest, specifically boundary nodes where we wish to specify boundary values
left_boundary_nodes = []
right_boundary_nodes = []
# DOC-START define node coordinates
# Read the geometric field 
for nodeNumber in exregion.nodeids:
    nodeDomain = decomposition.NodeDomainGet(nodeNumber,1)
    if(nodeDomain == computationalNodeNumber):
        version = 1
        derivative = 1
        coord = []
        for component in range(1, numberOfXi + 1):
            component_name = ["1", "2", "3"][component - 1]
            value = exregion.node_value("coordinates", component_name, nodeNumber, derivative)
            geometricField.ParameterSetUpdateNodeDP(
                iron.FieldVariableTypes.U,
                iron.FieldParameterSetTypes.VALUES,
                version, derivative, nodeNumber, component, value)
            coord.append(value)
        #The cylinder has an axial length of 5 units and coord[2] contains the axial coordinate value        
        if coord[2] < 0.001:
            left_boundary_nodes.append(nodeNumber)
        elif coord[2] > 4.999:
            right_boundary_nodes.append(nodeNumber)
           

# Update the geometric field parameters 
geometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
geometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
# DOC-END define node coordinates

# Create a fibre field and attach it to the geometric field
fibreField = iron.Field()
fibreField.CreateStart(fibreFieldUserNumber,region)
fibreField.TypeSet(iron.FieldTypes.FIBRE)
fibreField.MeshDecompositionSet(decomposition)
fibreField.GeometricFieldSet(geometricField)
fibreField.VariableLabelSet(iron.FieldVariableTypes.U,"Fibre")
fibreField.CreateFinish()

# Create the dependent field
#Create the dependent field with 4 variables and the respective number of components
#1   U_Var_Type            4 components: 3 displacement (quad interpol) + 1 pressure (lin interpol))
#2   DELUDELN_Var_Type     4 components: 3 displacement (quad interpol) + 1 pressure (lin interpol))
#3   U1_Var_Type           6 components: 6 independent components of the strain tensor (quad interpol) [independent]
#4   U2_Var_Type           6 components: 6 independent components of the stress tensor (quad interpol) [dependent]

dependentField = iron.Field()
dependentField.CreateStart(dependentFieldUserNumber,region)
dependentField.VariableLabelSet(iron.FieldVariableTypes.U,"Dependent")
dependentField.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)  
dependentField.MeshDecompositionSet(decomposition)
dependentField.GeometricFieldSet(geometricField) 
dependentField.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
#Displacement, Gradient, Strain and Stress Tensors 
dependentField.NumberOfVariablesSet(4)
dependentField.VariableTypesSet([iron.FieldVariableTypes.U,iron.FieldVariableTypes.DELUDELN,iron.FieldVariableTypes.U1,iron.FieldVariableTypes.U2])
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U,4)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN,4)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U1,6)
dependentField.NumberOfComponentsSet(iron.FieldVariableTypes.U2,6)

#Assign the mesh from which the quantities are determined
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,quadraticMeshComponentNumber)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,quadraticMeshComponentNumber)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,quadraticMeshComponentNumber)
  
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,1,quadraticMeshComponentNumber)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,2,quadraticMeshComponentNumber)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,3,quadraticMeshComponentNumber)

if(usePressureBasis):
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,linearMeshComponentNumber)  
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,linearMeshComponentNumber)
else:
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,quadraticMeshComponentNumber)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,quadraticMeshComponentNumber)
    
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U1,1,quadraticMeshComponentNumber)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U1,2,quadraticMeshComponentNumber)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U1,3,quadraticMeshComponentNumber)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U1,4,quadraticMeshComponentNumber)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U1,5,quadraticMeshComponentNumber)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U1,6,quadraticMeshComponentNumber)

dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U2,1,quadraticMeshComponentNumber)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U2,2,quadraticMeshComponentNumber)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U2,3,quadraticMeshComponentNumber)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U2,4,quadraticMeshComponentNumber)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U2,5,quadraticMeshComponentNumber)
dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U2,6,quadraticMeshComponentNumber)

dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,1,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,2,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,3,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,4,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,5,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U1,6,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)

dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,1,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,2,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,3,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,4,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,5,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
dependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U2,6,iron.FieldInterpolationTypes.GAUSS_POINT_BASED)

if(usePressureBasis):
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,4,2)
    dependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,4,2)

dependentField.VariableLabelSet(iron.FieldVariableTypes.U1,"strain")    
dependentField.VariableLabelSet(iron.FieldVariableTypes.U2,"stress")

dependentField.CreateFinish()

# Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometricField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,
    dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3)
#Set hydrostatic pressure, the value depends on the constitutive law being used 
iron.Field.ComponentValuesInitialiseDP(dependentField,iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,4,-8.0)

# Create the equations_set
equationsSetField = iron.Field()
equationsSet = iron.EquationsSet()
equationsSetSpecification = [iron.EquationsSetClasses.ELASTICITY,
        iron.EquationsSetTypes.FINITE_ELASTICITY,
        iron.EquationsSetSubtypes.CONSTITUTIVE_LAW_IN_CELLML_EVALUATE]
equationsSet.CreateStart(equationsSetUserNumber, region, fibreField,
        equationsSetSpecification, equationsSetFieldUserNumber, equationsSetField)
equationsSet.CreateFinish()

equationsSet.DependentCreateStart(dependentFieldUserNumber, dependentField)
equationsSet.DependentCreateFinish()


#DOC-START create CellML environment
# Create the CellML environment
cellML = iron.CellML()
cellML.CreateStart(cellMLUserNumber, region)
# Import a Mooney-Rivlin material law from a file
mooneyRivlinModel = cellML.ModelImport("mooney_rivlin.xml")
#DOC-END create cellml environment

#DOC-START flag variables
# Now we have imported the model we are able to specify which variables from the model we want to set from openCMISS
cellML.VariableSetAsKnown(mooneyRivlinModel, "equations/E11")
cellML.VariableSetAsKnown(mooneyRivlinModel, "equations/E12")
cellML.VariableSetAsKnown(mooneyRivlinModel, "equations/E13")
cellML.VariableSetAsKnown(mooneyRivlinModel, "equations/E22")
cellML.VariableSetAsKnown(mooneyRivlinModel, "equations/E23")
cellML.VariableSetAsKnown(mooneyRivlinModel, "equations/E33")
# and variables to get from the CellML 
cellML.VariableSetAsWanted(mooneyRivlinModel, "equations/Tdev11")
cellML.VariableSetAsWanted(mooneyRivlinModel, "equations/Tdev12")
cellML.VariableSetAsWanted(mooneyRivlinModel, "equations/Tdev13")
cellML.VariableSetAsWanted(mooneyRivlinModel, "equations/Tdev22")
cellML.VariableSetAsWanted(mooneyRivlinModel, "equations/Tdev23")
cellML.VariableSetAsWanted(mooneyRivlinModel, "equations/Tdev33")
#DOC-END flag variables

#DOC-START create cellml finish
cellML.CreateFinish()
#DOC-END create cellml finish


# Start the creation of CellML <--> OpenCMISS field maps
cellML.FieldMapsCreateStart()
#Now we can set up the field variable component <--> CellML model variable mappings.
#DOC-START map strain components
#Map the strain components
cellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,1, iron.FieldParameterSetTypes.VALUES,mooneyRivlinModel,"equations/E11", iron.FieldParameterSetTypes.VALUES)
cellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,2, iron.FieldParameterSetTypes.VALUES,mooneyRivlinModel,"equations/E12", iron.FieldParameterSetTypes.VALUES)
cellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,3, iron.FieldParameterSetTypes.VALUES,mooneyRivlinModel,"equations/E13", iron.FieldParameterSetTypes.VALUES)
cellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,4, iron.FieldParameterSetTypes.VALUES,mooneyRivlinModel,"equations/E22", iron.FieldParameterSetTypes.VALUES)
cellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,5, iron.FieldParameterSetTypes.VALUES,mooneyRivlinModel,"equations/E23", iron.FieldParameterSetTypes.VALUES)
cellML.CreateFieldToCellMLMap(dependentField,iron.FieldVariableTypes.U1,6, iron.FieldParameterSetTypes.VALUES,mooneyRivlinModel,"equations/E33", iron.FieldParameterSetTypes.VALUES)
#DOC-END map strain components
    
#DOC-START map stress components
#Map the stress components
cellML.CreateCellMLToFieldMap(mooneyRivlinModel,"equations/Tdev11", iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,1,iron.FieldParameterSetTypes.VALUES)
cellML.CreateCellMLToFieldMap(mooneyRivlinModel,"equations/Tdev12", iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,2,iron.FieldParameterSetTypes.VALUES)
cellML.CreateCellMLToFieldMap(mooneyRivlinModel,"equations/Tdev13", iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,3,iron.FieldParameterSetTypes.VALUES)
cellML.CreateCellMLToFieldMap(mooneyRivlinModel,"equations/Tdev22", iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,4,iron.FieldParameterSetTypes.VALUES)
cellML.CreateCellMLToFieldMap(mooneyRivlinModel,"equations/Tdev23", iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,5,iron.FieldParameterSetTypes.VALUES)
cellML.CreateCellMLToFieldMap(mooneyRivlinModel,"equations/Tdev33", iron.FieldParameterSetTypes.VALUES,dependentField,iron.FieldVariableTypes.U2,6,iron.FieldParameterSetTypes.VALUES)
#DOC-END map stress components

#Finish the creation of cellML <--> OpenCMISS field maps
cellML.FieldMapsCreateFinish()

#DOC-START define CellML models field
#Create the CellML models field
cellMLModelsField = iron.Field()
cellML.ModelsFieldCreateStart(cellMLModelsFieldUserNumber,cellMLModelsField)
cellML.ModelsFieldCreateFinish()

# Stress, strain fields are evaluated at Gauss Points and any fields used in the computation should also have values at Gauss points
# Iterate through each element and set it for the elements gauss points

xidiv = 1.0/(numberOfGaussXi+1)
for elment in exregion.elements:
    elmentDomain = decomposition.ElementDomainGet(elment.number)
    if(elmentDomain == computationalNodeNumber):
        #Gauss point number counter    
        gaussPointIdx = 0
        #Assign model for each quadraturePoint:
        for xi in range(0,numberOfGaussXi):
            xi1 = (1.0+xi)*xidiv
            for xj in range(0,numberOfGaussXi):
                xi2 = (1.0+xj)*xidiv
                for xk in range(0,numberOfGaussXi):
                    xi3 = (1.0+xk)*xidiv
                    gaussPointIdx = gaussPointIdx + 1
                    cellMLModelsField.ParameterSetUpdateGaussPoint(iron.FieldVariableTypes.U,
                                                                   iron.FieldParameterSetTypes.VALUES,
                                                                   gaussPointIdx,elment.number,1,
                                                                   mooneyRivlinModel)
# Update the models field parameters 
cellMLModelsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
cellMLModelsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
#DOC-END define CellML models field

#DOC-START define CellML parameters and intermediate fields
#Create the CellML parameters field --- the strain field
cellMLParametersField = iron.Field()
cellML.ParametersFieldCreateStart(cellMLParametersFieldUserNumber,cellMLParametersField)
cellML.ParametersFieldCreateFinish()

#  Create the CellML intermediate field --- the stress field
cellMLIntermediateField = iron.Field()
cellML.IntermediateFieldCreateStart(cellMLIntermediateFieldUserNumber,cellMLIntermediateField)
cellML.IntermediateFieldCreateFinish()
#DOC-END define CellML parameters and intermediate fields

# Create equations
equations = iron.Equations()
equationsSet.EquationsCreateStart(equations)
equations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
equationsSet.EquationsCreateFinish()

#DOC-START define CellML finite elasticity problem
#Define the problem
problem = iron.Problem()
problemSpecification = [iron.ProblemClasses.ELASTICITY,
    iron.ProblemTypes.FINITE_ELASTICITY,
    iron.ProblemSubtypes.FINITE_ELASTICITY_WITH_CELLML]
problem.CreateStart(problemUserNumber, problemSpecification)
problem.CreateFinish()
#DOC-END define CellML finite elasticity problem

#Create the problem control loop
problem.ControlLoopCreateStart()
ControlLoop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],ControlLoop)
ControlLoop.TypeSet(iron.ProblemControlLoopTypes.SIMPLE)
problem.ControlLoopCreateFinish()

#Create the problem solvers
nonLinearSolver = iron.Solver()
linearSolver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,nonLinearSolver)
nonLinearSolver.OutputTypeSet(iron.SolverOutputTypes.MONITOR)
nonLinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD)
nonLinearSolver.NewtonLinearSolverGet(linearSolver)
#Use the DIRECT MUMPS solver
linearSolver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
#For large problems or problems with complex material behaviour, the direct solver may fail
#In such cases either preconditioners or the following solvers can be tried
#In case the matrix has zeros for some rows, SUPERLU is a good solver to try as it will report such errors
#linearSolver.LibraryTypeSet(iron.SolverLibraries.PASTIX)
#linearSolver.LibraryTypeSet(iron.SolverLibraries.SUPERLU)
problem.SolversCreateFinish()

#DOC-START define CellML solver
#Create the problem solver CellML equations
cellMLSolver = iron.Solver()
problem.CellMLEquationsCreateStart()
nonLinearSolver.NewtonCellMLSolverGet(cellMLSolver)
cellMLEquations = iron.CellMLEquations()
cellMLSolver.CellMLEquationsGet(cellMLEquations)
cellMLEquations.CellMLAdd(cellML)
problem.CellMLEquationsCreateFinish()
#DOC-END define CellML solver

#Create the problem solver equations
solver = iron.Solver()
solverEquations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,solver)
solver.SolverEquationsGet(solverEquations)
solverEquations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
equationsSetIndex = solverEquations.EquationsSetAdd(equationsSet)
problem.SolverEquationsCreateFinish()


# Prescribe boundary conditions (absolute nodal parameters)
boundaryConditions = iron.BoundaryConditions()
solverEquations.BoundaryConditionsCreateStart(boundaryConditions)

#Here we model axial stretch, by pulling along the axis at both the ends of the cylinder and holding them (Direchlet)
for nodeNumber in left_boundary_nodes:
    nodeDomain = decomposition.NodeDomainGet(nodeNumber, 1)
    if nodeDomain == computationalNodeNumber :
        boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nodeNumber, 1, iron.BoundaryConditionsTypes.FIXED, 0.0)
        boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nodeNumber, 2, iron.BoundaryConditionsTypes.FIXED, 0.0)
        boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nodeNumber, 3, iron.BoundaryConditionsTypes.FIXED, -1.0)

for nodeNumber in right_boundary_nodes:
    nodeDomain = decomposition.NodeDomainGet(nodeNumber, 1)
    if nodeDomain == computationalNodeNumber :
        boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nodeNumber, 1, iron.BoundaryConditionsTypes.FIXED, 0.0)
        boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nodeNumber, 2, iron.BoundaryConditionsTypes.FIXED, 0.0)
        boundaryConditions.AddNode(dependentField, iron.FieldVariableTypes.U, 1, 1, nodeNumber, 3, iron.BoundaryConditionsTypes.FIXED, 1.0)

solverEquations.BoundaryConditionsCreateFinish()


#The MUMPS solver may fail for large problems unable to allocate memory.
#The solver can be instructed to allocate large heap space apriori
#These settings are done here as the Solver is not assigned until problem Solver Equations are finalized
#linearSolver.MumpsSetIcntl(23,640000)
#linearSolver.MumpsSetIcntl(14,640000)
# Solve the problem
problem.Solve()

# Export the results, here we export them as standard exnode, exelem files
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("./results/AxialStretch","FORTRAN")
fields.ElementsExport("./results/AxialStretch","FORTRAN")
fields.Finalise()

