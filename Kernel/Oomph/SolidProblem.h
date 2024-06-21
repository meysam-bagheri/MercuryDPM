//Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
//For the list of developers, see <http://MercuryDPM.org/Team>.
//
//Redistribution and use in source and binary forms, with or without
//modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name MercuryDPM nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/**
 * For the material model, see
 * http://oomph-lib.maths.man.ac.uk/doc/solid/solid_theory/html/index.html
 **/

#ifndef SOLID_PROBLEM_H
#define SOLID_PROBLEM_H

// Generic oomph-lib headers
#include "generic.h"
#include "solid.h"
#include "constitutive.h"
#include "Math/Vector.h"
#include <array>
#include <cstdio>
#include "OomphHelpers.h"
#include "AnisotropicHookean.h"
#include "RefineableQDPVDElement.h"

//The mesh
#include "meshes/simple_cubic_mesh.h"
#include "meshes/tetgen_mesh.h"

#include "Oomph/ScaleCoupling/ScaleCoupledElement.h"

using namespace oomph;

#ifdef OOMPH_HAS_MPI
    #define OOMPH_MPI_PROCESSOR_NUM communicator_pt()->nproc()
#else
    #define OOMPH_MPI_PROCESSOR_NUM 1
#endif

#ifdef OOMPH_HAS_MPI
    #define OOMPH_MPI_PROCESSOR_ID communicator_pt()->my_rank()
#else
    #define OOMPH_MPI_PROCESSOR_ID 0
#endif

/**
 * Base class for (surface-coupled) problems with a solid body.
 *
 * Assumptions:
 *  - The element type is derived from QPVDElement, i.e. a solid.
 *  - The mesh type is derived from SolidMesh.
 *  - The constitutive law is derived from Hookean, with parameters elasticModulus, poissonRatio, density.
 *
 * The solid properties (elasticModulus, poissonRatio, density, constitutive_law_pt, body_force_fct, etc.) are stored in member variables, not in globals as typical in oomph-lib, and accessed by setters and getters.
 *
 * Additional functionality:
 * - setIsPinned allows one to define the pinned nodes using a function.
 */
template<class ELEMENT_TYPE>
class SolidProblem : public Problem
{
protected:
    //define element type (should be done by template)
    typedef ELEMENT_TYPE ELEMENT;

    // name of output files (user-defined)
    std::string name_;

    /// Elastic modulus (set via setSolid)
    double elasticModulus_ = 0;

    /// Poisson's ratio (set via setSolid)
    double poissonRatio_ = 0;

    /// Density
    double density_ = 0;

    /// Density
    double gravity_ = 0;

    /// Pointer to the body force function
    void (* body_force_fct)(const double& time, const Vector<double>& xi, Vector<double>& b) = nullptr;

    /// Pointer to constitutive law (should be set in constructor)
    ConstitutiveLaw* constitutive_law_pt  = nullptr;

    /// Pointer to solid mesh
    SolidMesh* Solid_mesh_pt = nullptr;

    /// Pointer to mesh of traction elements
    SolidMesh* Traction_mesh_pt = nullptr;

    /// Function to determine which nodal positions are pinned
    std::function<bool(SolidNode*, unsigned)> isPinned_;

public:

    /// Constructor: set default constitutive law and time stepper
    SolidProblem()
    {
        logger(INFO, "Set default constitutive law (AnisotropicHookean) and time stepper (Newmark<2>)");

        // Set Newton solver tolerance and maximum number of Newton iterations
        Max_newton_iterations = 20;
        Newton_solver_tolerance = 1e-10;
        Max_residuals = constants::inf;

        // Set constitutive law
        constitutive_law_pt = new AnisotropicHookean(&poissonRatio_, &elasticModulus_);

        // Allocate the timestepper
        add_time_stepper_pt(new Newmark<2>);

        //build_mesh();

        // setup boundary conditions
        //pin_and_assign_eqn_numbers();
    };

    /// set function for name_
    void setName(const std::string& name)
    {
        name_ = name;
        logger(INFO, "Name: %", name_);
    }

    /// get function for name_
    std::string getName() const
    {
        return name_;
    }

    /// set function for elasticModulus_
    void setElasticModulus(double elasticModulus)
    {
        elasticModulus_ = elasticModulus;
        logger(INFO, "Elastic Modulus: %", elasticModulus_);
    }

    /// get function for elasticModulus_
    double getElasticModulus() const
    {
        return elasticModulus_;
    }

    // add dissipation
    void addDissipation(double dissipation)
    {
        logger(INFO,"Adding dissipation %",dissipation);
        for (int i = 0; i < solid_mesh_pt()->nelement(); ++i)
        {
            dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(i))->dissipation_ = dissipation;
        }
    }

    void addAnisotropy(double Ex, double Ey, double Ez)
    {
        logger(INFO,"Adding anisotropy E = [% % %]",Ex, Ey, Ez);
        // make constitutive law anisotropic in x-direction
        elasticModulus_ = Ex;
        auto hooke_law_pt = dynamic_cast<AnisotropicHookean*>(constitutive_law_pt);
        hooke_law_pt->setAnisotropy({1.0, Ey/Ex, Ez/Ex});
    }

    /// set function for elasticModulus_
    void setOomphGravity(double gravity)
    {
        gravity_ = gravity;
        logger(INFO, "Elastic Modulus: %", gravity_);
    }

    /// get function for gravity_
    double getOomphGravity() const
    {
        return gravity_;
    }

    /// set function for poissonRatio_
    void setPoissonRatio(double poissonRatio)
    {
        poissonRatio_ = poissonRatio;
        logger(INFO, "Poisson Ratio: %", poissonRatio_);
    }

    /// get function for poissonRatio_
    double getPoissonRatio() const
    {
        return poissonRatio_;
    }

    /// set function for density_
    void setDensity(double density)
    {
        density_ = density;
        logger(INFO, "Density: %", density_);
    }

    /// get function for density_
    double getDensity() const
    {
        return density_;
    }

    /// set function for body_force_pt
    void setBodyForceAsGravity(double gravity = 9.8)
    {
        logger(INFO, "Setting oomph-gravity in z-direction");
        gravity_ = gravity;
        // define a static body force
        static double& Density = density_;
        static double& Gravity = gravity_;
        body_force_fct = [](const double& time, const Vector<double>& xi, Vector<double>& b) {
            b[0] = 0.0;
            b[1] = 0.0;
            b[2] = -Gravity * Density;
        };
    }

    /// set function for isPinned_
    void setIsPinned(std::function<bool(SolidNode*, unsigned)> isPinned)
    {
        isPinned_ = std::move(isPinned);
        logger(INFO, "Setting which positions on which nodes are pinned");
    }

    // set is_pinned such that a certain boundary is pinned
    void pinBoundary(unsigned b) {
        logger(INFO, "Pinning nodes on boundary %", b);
        isPinned_ = [b](SolidNode *n, unsigned d) {
            return n->is_on_boundary(b);
        };
    }

    // set is_pinned such that a certain boundary is pinned
    void pinBoundaries(std::vector<unsigned> b) {
        for (const auto a: b) {
            logger(INFO, "Pinning nodes on boundary %", a);
        }
        isPinned_ = [b](SolidNode *n, unsigned d) {
            for (const auto a : b) {
                if (n->is_on_boundary(a)) return true;
            }
            return false;
        };
    }

    /// Sets the dissipation coefficient for all elements.
    std::enable_if<std::is_base_of<RefineableQDPVDElement<3,2>, ELEMENT>::value, void> setDissipation(double dissipation)
    {
        for (int i = 0; i < solid_mesh_pt()->nelement(); ++i)
        {
            dynamic_cast<RefineableQDPVDElement<3,2>*>(solid_mesh_pt()->element_pt(i))->setDissipation(dissipation);
        }
    }

    /// set function for Newton_solver_tolerance
    void setNewtonSolverTolerance(double Newton_solver_tolerance)
    {
        this->Newton_solver_tolerance = Newton_solver_tolerance;
    }

    /// set function for Max_newton_iterations
    void setMaxNewtonIterations(unsigned Max_newton_iterations)
    {
        this->Max_newton_iterations = Max_newton_iterations;
    }

    /// set function for time step
    void setOomphTimeStep(double dt) {
        time_pt()->initialise_dt(dt);
    }

    /// get function for time step
    double getOomphTimeStep () const {
        return time_pt()->dt();
    }

    /// get function for current time
    double getOomphTime () const {
        return time_pt()->time();
    }

    /// Get function for the solid mesh pointer
    SolidMesh*& solid_mesh_pt() { return Solid_mesh_pt; }

    /// Get function for the traction mesh pointer
    SolidMesh*& traction_mesh_pt() { return Traction_mesh_pt; }

    /// Get function for the solid mesh pointer
    SolidMesh* const & solid_mesh_pt() const { return Solid_mesh_pt; }

    /// Get function for the traction mesh pointer
    SolidMesh* const & traction_mesh_pt() const { return Traction_mesh_pt; }

    /// Computes the domain size (min/max of the nodal positions in x/y/z)
    void getDomainSize(std::array<double, 3>& min, std::array<double, 3>& max) const
    {
        min[0] = min[1] = min[2] = constants::inf;
        max[0] = max[1] = max[2] = -constants::inf;
        logger.assert_always(solid_mesh_pt(),"mesh not found");
        for (unsigned i = 0; i < solid_mesh_pt()->nnode(); i++)
        {
            const auto n = solid_mesh_pt()->node_pt(i);
            for (int j = 0; j < 3; ++j)
            {
                min[j] = std::min(min[j], n->xi(j));
                max[j] = std::max(max[j], n->xi(j));
            }
        }
    }

    /**
     * - Asserts all parameters are set
     * - Assign pointers for all elements
     * - Builds global mesh
     * - Pins nodes
     * - Pins solid pressure
     * - Assigns equation numbers
     *
     * This function should be called before solve(), after creating the mesh, defining body force and constitutive law.
     */
    void prepareForSolve()
    {
        // check certain values are set
        logger.assert_always(!name_.empty(), "Set name via setName(..)");
        logger.assert_always(elasticModulus_>0, "Set elastic modulus via setElasticModulus(..)");
        logger.assert_always(density_>0, "Set density via setDensity(..)");
        logger.assert_always(solid_mesh_pt(), "Set solid mesh via e.g. setSolidCubicMesh(..)");

        // Assign constitutive_law_pt and body_force_fct_pt of each element
        logger(INFO, "Assign constitutive_law, body_force, density to all elements");
        for (unsigned i = 0; i < solid_mesh_pt()->nelement(); i++)
        {
            //Cast to a solid element
            ELEMENT* el_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(i));
            ///\todo UMBU Note, prepareForSolve did not set lambda correctly before
            // Set the constitutive law
            el_pt->constitutive_law_pt() = constitutive_law_pt;
            // Set body force
            el_pt->body_force_fct_pt() = body_force_fct;
            // Set density
            el_pt->lambda_sq_pt() = &density_;
        }

        // Construct the traction element mesh
        // Traction_mesh_pt = new SolidMesh;
        // create_traction_elements();

        //Build combined "global" mesh
        add_sub_mesh(solid_mesh_pt());
        if (traction_mesh_pt()) {
            add_sub_mesh(traction_mesh_pt());
            logger(INFO, "Built global mesh from solid mesh and traction mesh");
        } else {
            logger(INFO, "Built global mesh from solid mesh");
        }
        build_global_mesh();

        //Pin nodes
        if (isPinned_)
        {
            for ( unsigned n = 0; n < solid_mesh_pt()->nnode(); n++ )
            {
                SolidNode* n_pt = solid_mesh_pt()->node_pt( n );
                //Pin all nodes
                for ( unsigned i = 0; i < 3; i++ )
                {
                    if ( isPinned_( n_pt, i ) )
                    {
                        n_pt->pin_position( i );
                    }
                    else
                    {
                        n_pt->unpin_position( i );
                    }
                }
            }
        }
        countPinned();

        // Pin the redundant solid pressures (if any)
        PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(
            solid_mesh_pt()->element_pt());
        logger(INFO, "Pinned redundant nodal solid pressures");

        // Attach the boundary conditions to the mesh
        unsigned n_eq = assign_eqn_numbers();
        logger(INFO, "Assigned % equation numbers", n_eq);
    }

    /// returns statistics about pinned nodes to the console
    void countPinned()
    {
        // count pinned
        std::array<unsigned,3> countPinned {0,0,0};
        unsigned countAll = 0;
        for (unsigned n = 0; n < solid_mesh_pt()->nnode(); n++)
        {
            for (unsigned i = 0; i < 3; i++)
            {
                countPinned[i] += solid_mesh_pt()->node_pt(n)->position_is_pinned(i);
                countAll++;
            }
        }
        unsigned countPinnedAll = countPinned[0] + countPinned[1] +countPinned[2];
        logger(INFO, "Pinned % of % positions (% free): % in x, % in y, % in z", countPinnedAll, countAll, countAll - countPinnedAll, countPinned[0], countPinned[1], countPinned[2]);
    }

    /**
     * Solves a steady problem.
     */
    void solveSteady()
    {
        logger.assert_always(mesh_pt(), "Mesh pointer not set; did you call prepareForSolve?");
        logger(INFO, "Solve steady-state problem");
        actionsBeforeSolve();
        newton_solve();
        actionsAfterSolve();
        writeToVTK();
        saveSolidMesh();
    }

    virtual void actionsBeforeSolve() {}

    virtual void actionsAfterSolve() {}

    virtual void actionsBeforeOomphTimeStep() {}

    /**
     * Solves an unsteady problem.
     */
    void solveUnsteady(double timeMax, double dt, unsigned saveCount = 10)
    {
        logger.assert_always(mesh_pt(), "Mesh pointer not set; did you call prepareForSolve?");

        std::cout << "Solving oomph with dt=" << dt << " until timeMax=" << timeMax << std::endl;

        // Setup initial conditions. Default is a non-impulsive start
        this->time_pt()->initialise_dt(dt);
        assign_initial_values_impulsive(dt);
        actionsBeforeSolve();

        // This is the main loop over advancing time
        double& time = time_stepper_pt()->time();
        unsigned count = 0;
        const unsigned countMax = std::ceil(timeMax/dt);
        while (time < timeMax)
        {
            logger(INFO, "Time %s of %s (% of %)", time, timeMax, count, countMax);
            actionsBeforeOomphTimeStep();
            // solve the oomphProb for one time step (this also increments time)
            ///\todo UMBU changed from
            //adaptive_unsteady_newton_solve(dt,2e-7);
            unsteady_newton_solve(dt);
            // increase count
            ///\todo UMBU Changed counting back to the old way, otherwise we did not plot the same time step
            //count++;
            // write outputs of the oomphProb
            if (count++ % saveCount == 0 or time + dt > timeMax) {
                writeToVTK();
                // save in case code breaks
                // saveSolidMesh();
            }
            // abort if problem
            if (Max_res.back()==0) {
                saveSolidMesh();
                logger(ERROR, "Maximum residual is 0; exiting the code");
            }
        }
        saveSolidMesh();
        actionsAfterSolve();
    }

    /**
     * Solves an unsteady problem, returns successful if timeMaxMin has been reached
     */
    void solveUnsteadyForgiving(double timeMax, double timeMaxMin, double dt, unsigned saveCount = 10) {
        // solve
        try {
            solveUnsteady(timeMax, dt, saveCount);
        } catch(OomphLibError& error)  {
            //Store output if newton solver fails
            saveSolidMesh();
            if (time_stepper_pt()->time()-dt >= timeMaxMin) {
                // take it as successful if a fraction of the time evolution has finished
                logger(INFO,"Newton solver failed at t=% (tMax=%), but code will continue.",
                       time_stepper_pt()->time()-dt, timeMax);
                exit(0);
            } else {
                logger(ERROR,"Newton solver failed at t=% (tMax=%).",
                       time_stepper_pt()->time()-dt, timeMax);
            }
        }
    }

    /**
     * Removes old output files (vtu) with the same name as this problem.
     */
    void removeOldFiles() const
    {
        for (int i = 0; true; ++i)
        {
            std::string fileName = name_ + "FEM_" + std::to_string(i) + ".vtu";
            if (remove(fileName.c_str())) break;
            std::cout << "Deleted " << fileName << '\n';
        }
    }

    /**
     * Returns the x value at a certain xi value
     */
    void get_x(const Vector<double>& xi, Vector<double>& x) const
    {
        if (OOMPH_MPI_PROCESSOR_NUM>1) {
            logger(INFO, "get_x does not work with MPI");
            for (int i = 0; i < 3; ++i) {
                x[i] = xi[i];
            }
        } else {
            Vector<double> s(3);
            GeomObject *geom_obj_pt = nullptr;
            const unsigned long nelement = solid_mesh_pt()->nelement();
            for (unsigned long i = 0; i < nelement; i++) {
                auto el_pt = dynamic_cast<ELEMENT *>(solid_mesh_pt()->element_pt(i));
                el_pt->locate_zeta(xi, geom_obj_pt, s);
                if (geom_obj_pt) {
                    //logger(INFO,"Point % % % is in element % at % % %",
                    //       xi[0],xi[1],xi[2],i,s[0], s[1], s[2]);
                    el_pt->interpolated_x(s, x); //deformed coordinate
                    return;
                }
            }
            logger(ERROR, "x(xi) could not be found");
        }
    }

    /**
     * Returns the difference between the x and xi value in dimension d at a certain xi value
     */
    double getDeflection(Vector<double> xi, unsigned d) const
    {
        Vector<double> x(3);
        get_x(xi, x);
        return x[d] - xi[d];
    }

    // Boundary types (from cubic mesh)
    enum Boundary : unsigned
    {
        Z_MIN = 0, Y_MIN = 1, X_MAX = 2, Y_MAX = 3, X_MIN = 4, Z_MAX = 5
    };

    // Simple cubic mesh upgraded to become a solid mesh
    class SolidCubicMesh : public virtual RefineableSimpleCubicMesh<ELEMENT>, public virtual SolidMesh
    {
    public:
        /// \short Constructor:
        // nx, ny, nz: number of elements in the x, y, and z directions
        // xMax, yMax, zMax: dimensions of the cube (assume the center of the cube is at the origin)
        // timeStepper: defaults to Steady.
        SolidCubicMesh(const unsigned& nx, const unsigned& ny, const unsigned& nz,
                       const double& xMin, const double& xMax, const double& yMin,
                       const double& yMax, const double& zMin, const double& zMax,
                       TimeStepper* time_stepper_pt) :
                SimpleCubicMesh<ELEMENT>(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax, time_stepper_pt),
                RefineableSimpleCubicMesh<ELEMENT>(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax, time_stepper_pt),
                SolidMesh()
        {
            //Assign the initial lagrangian coordinates
            set_lagrangian_nodal_coordinates();
        }
    };

    // Create a solid cubic mesh
    void setSolidCubicMesh(const unsigned& nx, const unsigned& ny, const unsigned& nz,
                           const double& xMin, const double& xMax, const double& yMin,
                           const double& yMax, const double& zMin, const double& zMax)
    {
        // Create mesh
        solid_mesh_pt() = new SolidCubicMesh(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax, time_stepper_pt());

        logger(INFO, "Created %x%x% cubic mesh on domain [%,%]x[%,%]x[%,%]",
               nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax);
    }

    void saveSolidMesh()
    {
        std::string filename = name_;
#ifdef OOMPH_HAS_MPI
        // When the problem hasn't been distributed, the whole mesh is stored on every processor.
        // Prevent writing the same file multiple times, so only write for the first processor.
        if (!Problem_has_been_distributed && OOMPH_MPI_PROCESSOR_ID != 0)
            return;

        // When the problem has been distributed, add the processor id to the filename.
        if (Problem_has_been_distributed)
            filename += "_" + std::to_string(OOMPH_MPI_PROCESSOR_ID);
#endif
        filename += ".mesh";

        //if (useTetgen()) logger(ERROR, "Not implemented for Tetgen meshes");

        std::ofstream mesh(filename);
        auto solid_cubic_mesh_pt = dynamic_cast<SolidCubicMesh*>(solid_mesh_pt());
        if (solid_cubic_mesh_pt) {
            //logger.assert_always(solid_cubic_mesh_pt, "Mesh is not cubic");
            mesh << solid_cubic_mesh_pt->nx() << ' ';
            mesh << solid_cubic_mesh_pt->ny() << ' ';
            mesh << solid_cubic_mesh_pt->nz() << '\n';
        }
        for (int i = 0; i < solid_mesh_pt()->nnode(); ++i)
        {
            SolidNode* n = solid_mesh_pt()->node_pt(i);
            for (int j = 0; j < 3; ++j)
            {
                mesh << n->xi(j) << ' ' << n->x(j) << ' ' << n->position_is_pinned(j) << ' ';
            }
            mesh << '\n';
        }
        if (solid_cubic_mesh_pt) {
            logger(INFO, "Saved %x%x% mesh to %",
                   solid_cubic_mesh_pt->nx(), solid_cubic_mesh_pt->ny(), solid_cubic_mesh_pt->nz(), filename);
        } else {
            logger(INFO, "Saved mesh to % (% nodes)", filename, solid_mesh_pt()->nnode());
        }
    }

    void loadSolidMesh(std::string infileName, bool cubic=true)
    {
        logger(INFO, "Loading % (cubic=%)",infileName, cubic);
        //if (useTetgen()) logger(ERROR, "Not implemented for Tetgen meshes");

        std::ifstream mesh(infileName);
        logger.assert_always(mesh.good(),"Mesh file % could not be opened",infileName);

        if (cubic) {
            unsigned nx, ny, nz;
            mesh >> nx >> ny >> nz;
            logger(INFO, "Loaded %x%x% cubic mesh from %", nx, ny, nz, infileName);
            logger.assert_debug(nx > 1 and ny > 1 and nz > 1, "Mesh size invalid");
            Solid_mesh_pt = new SolidCubicMesh(nx, ny, nz, 0, 1, 0, 1, 0, 1, time_stepper_pt());
            // Assign physical properties to the elements before any refinement
            for (unsigned i = 0; i < solid_mesh_pt()->nelement(); i++)
            {
                //Cast to a solid element
                ELEMENT* el_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(i));
                // Set the constitutive law
                el_pt->constitutive_law_pt() = constitutive_law_pt;
            }
        }

        double xi, x;
        bool pin;
        logger(INFO, "Loading % nodes", solid_mesh_pt()->nnode());
        for (int i = 0; i < solid_mesh_pt()->nnode(); ++i)
        {
            SolidNode* n = solid_mesh_pt()->node_pt(i);
            for (int j = 0; j < 3; ++j)
            {
                mesh >> xi >> x >> pin;
                n->xi(j) = xi;
                n->x(j) = x;
                if (pin) n->pin_position(j);
            }
        }
    }

    // stores results in vtk file
    void writeToVTK()
    {
#ifdef OOMPH_HAS_MPI
        // When the problem hasn't been distributed, the whole mesh is stored on every processor.
        // Prevent writing the same file multiple times, so only write for the first processor.
        if (!Problem_has_been_distributed && OOMPH_MPI_PROCESSOR_ID != 0)
            return;
#endif

        //set local coordinates list
        std::vector<std::vector<double>> sList0;
        // order of nodes
        std::vector<unsigned> nList;
        // vtkFormat for given shape
        // https://kitware.github.io/vtk-examples/site/VTKFileFormats/
        unsigned vtkFormat;

        // define differently for tetrahedral and hexahedral (quad) elements
        if (std::is_base_of<SolidTElement<3, 2>, ELEMENT>::value)
        {
            vtkFormat = 10; // TETRA

            sList0 = {
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, 1},
                {0, 0, 0}
            };

            nList = {0, 1, 2, 3};
        }
        else if (std::is_base_of<SolidQElement<3, 2>, ELEMENT>::value)
        {
            vtkFormat = 12; // HEXAHEDRON

            sList0 = {
                {-1, -1, -1},
                {-1, -1, +1},
                {-1, +1, +1},
                {-1, +1, -1},
                {+1, -1, -1},
                {+1, -1, +1},
                {+1, +1, +1},
                {+1, +1, -1}
            };
            // order of nodes
            nList = {0, 4, 6, 2, 1, 5, 7, 3};
        } else {
            logger(ERROR,"Element type unknown");
        }

        // convert to Vector
        std::vector<Vector<double>> sList;
        for (auto s0 : sList0)
        {
            Vector<double> s(3);
            s[0] = s0[0];
            s[1] = s0[1];
            s[2] = s0[2];
            sList.push_back(s);
        }

        // number of cells
        unsigned nel = solid_mesh_pt()->nelement();

        // set up vtu Points
        struct Point
        {
            Vector<double> coordinate;
            struct Data
            {
                std::string name;
                std::vector<double> value;
            };
            std::vector<Data> data;
        };
        std::vector<Point> points;
        points.reserve(nel * sList.size());

        // set up vtu Cells
        struct Cell
        {
            std::vector<unsigned long> connectivity;
            unsigned long offset;
            unsigned type;
        };
        std::vector<Cell> cells;
        points.reserve(nel);

        //for all elements
        for (unsigned e = 0; e < nel; e++)
        {
            // Get pointer to element
            auto el_pt = dynamic_cast<ELEMENT*>(
                solid_mesh_pt()->element_pt(e));

            std::vector<unsigned long> connectivity;
            unsigned ix = 0; // node index
            for (const auto& s : sList)
            {
                // pointer to node info
                auto n = nList[ix];
                auto n_pt = dynamic_cast<SolidNode*>(el_pt->node_pt(n));
                // Get Eulerian coordinates
                Vector<double> x(3);
                el_pt->interpolated_x(s, x);
                // get velocity
                Vector<double> dxdt(3);
                el_pt->interpolated_dxdt(s, 1, dxdt);
                // get displacement
                Vector<double> xi(3);
                el_pt->interpolated_xi(s, xi);
                // getBodyForce
                Vector<double> body_force(3);
                auto bodyForceFct = el_pt->body_force_fct_pt();
                if (bodyForceFct) bodyForceFct(time(), xi, body_force);
                // get coupling residual (fails)
                //std::vector<double> couplingResidual
                //    {el_pt->get_nodal_coupling_residual(n, 0),
                //     el_pt->get_nodal_coupling_residual(n, 1),
                //     el_pt->get_nodal_coupling_residual(n, 2)};
                // get stress/pressure
                DenseMatrix<double> sigma(3,3);
                el_pt->get_stress(s, sigma);
                // first invariant
                double pressure = (sigma(0,0)+sigma(1,1)+sigma(2,2))/3;
                //https://en.wikipedia.org/wiki/Von_Mises_yield_criterion
                double J2 = 0;
                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        J2 += 0.5*sigma(i,j)*sigma(i,j);
                    }
                }
                //std::cout << " " << J2;
                // get strain
                DenseMatrix<double> strain(3,3);
                el_pt->get_strain(s, strain);
                // get velocity (fails)
                std::vector<double> dudt
                    {el_pt->dnodal_position_dt(n, 0),
                     el_pt->dnodal_position_dt(n, 1),
                     el_pt->dnodal_position_dt(n, 2)};
                // get boundary (works)
                std::set<unsigned>* boundaries_pt;
                n_pt->get_boundaries_pt(boundaries_pt);
                double b = boundaries_pt ? *boundaries_pt->begin()+1 : 0;
                std::vector<double> pin {(double) n_pt->position_is_pinned(0),
                                         (double) n_pt->position_is_pinned(1),
                                         (double) n_pt->position_is_pinned(2)};
                points.push_back(
                    {x, {
                        {"Velocity", {dxdt[0], dxdt[1], dxdt[2]}},
                        {"Displacement", {x[0] - xi[0], x[1] - xi[1], x[2] - xi[2]}},
                        {"BodyForce", {body_force[0], body_force[1], body_force[2]}},
                        {"Pin", pin},
                        //{"CouplingResidual", couplingResidual},
                        {"Velocity2", dudt},
                        {"Pressure", {pressure}},
                        {"J2", {J2}},
                        {"Stress", {sigma(0,0), sigma(0,1), sigma(0,2),
                                    sigma(1,0), sigma(1,1), sigma(1,2),
                                    sigma(2,0), sigma(2,1), sigma(2,2)}},
//                             {"Strain", {strain(0,0), strain(0,1), strain(0,2),
//                                           strain(1,0), strain(1,1), strain(1,2),
//                                           strain(2,0), strain(2,1), strain(2,2)}},
                        //\todo undo this comment, but it causes a problem in the output
                        {"Boundary", {b}}
                    }});

                // add to connectivity
                connectivity.push_back(points.size() - 1);
                ix++;
            }
            cells.push_back({connectivity, points.size(), vtkFormat});
        }

        //open vtk file
        static unsigned count = 0;
#ifdef OOMPH_HAS_MPI
        std::string vtkFileName = name_ + "FEM_" + std::to_string(OOMPH_MPI_PROCESSOR_ID) + "_" + std::to_string(count++) + ".vtu";
#else
        std::string vtkFileName = name_ + "FEM_" + std::to_string(count++) + ".vtu";
#endif

        std::ofstream vtk(vtkFileName);

        //write vtk file
        vtk << "<?xml version=\"1.0\"?>\n"
               "<!-- time 10.548-->\n"
               "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
               "<UnstructuredGrid>\n"
               "<Piece NumberOfPoints=\""
            << points.size()
            << "\" NumberOfCells=\""
            << cells.size()
            << "\">\n"
               "\n"
               "<Points>\n"
               "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (auto point : points)
        {
            vtk << point.coordinate[0] << " " << point.coordinate[1] << " " << point.coordinate[2] << "\n";
        }
        vtk << "</DataArray>\n"
               "</Points>\n\n";
        if (not points.empty())
        {
            vtk << "<PointData  Vectors=\"vector\">\n";
            for (int i = 0; i < points.front().data.size(); ++i)
            {
                auto data = points.front().data[i];
                vtk << R"(<DataArray type="Float32" Name=")"
                    << points.front().data[i].name
                    << R"(" NumberOfComponents=")"
                    << points.front().data[i].value.size()
                    << R"(" format="ascii">)" << "\n";
                for (const Point& point : points)
                {
                    for (const auto& value : point.data[i].value)
                    {
                        //ensure values are not too small (VTK issue)
                        vtk << ((!isfinite(value) or fabs(value)<1e-33 or fabs(value)>1e33)?0.0:value) << " ";
                    }
                }
                vtk << "\n"
                       "</DataArray>\n";
            }
            vtk << "</PointData>\n"
                   "\n";
        }
        vtk << "<Cells>\n"
               "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (const Cell& cell : cells)
        {
            for (auto point : cell.connectivity)
            {
                vtk << point << " ";
            }
        }
        vtk << "</DataArray>\n"
               "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        for (const Cell& cell : cells)
        {
            vtk << cell.offset << " ";
        }
        vtk << "</DataArray>\n"
               "<DataArray type=\"UInt8\"  Name=\"types\" format=\"ascii\">\n";
        for (const Cell& cell : cells)
        {
            vtk << cell.type << " ";
        }
        vtk << "\n"
               "</DataArray>\n"
               "</Cells>\n"
               "\n"
               "</Piece>\n"
               "</UnstructuredGrid>\n"
               "</VTKFile>\n";
        vtk.close();
        logger(INFO, "Written %", vtkFileName);
    }

    /// See PVDEquationsBase<DIM>::get_energy
    void getMassMomentumEnergy(double& mass, Vector<double>& com, Vector<double>& linearMomentum, Vector<double>& angularMomentum, double& elasticEnergy,
                               double& kineticEnergy) const {
        // Initialise mass
        mass = 0;
        // Initialise center of mass
        com.initialise(0);
        // Initialise momentum
        linearMomentum.initialise(0);
        angularMomentum.initialise(0);
        // Initialise energy
        elasticEnergy = 0;
        kineticEnergy = 0;

        // For each element
        const unsigned long nelement = this->solid_mesh_pt()->nelement();
        for (unsigned long e = 0; e < nelement; e++) {

            // Get the pointer to the element
            ELEMENT* e_pt = dynamic_cast<ELEMENT*>(Solid_mesh_pt->element_pt(e));

            const unsigned DIM = e_pt->dim();

            //Find out how many integration points there are
            unsigned n_intpt = e_pt->integral_pt()->nweight();

            //Set the Vector to hold local coordinates
            Vector<double> s(DIM);

            //Find out how many nodes there are
            const unsigned n_node = e_pt->nnode();

            //Find out how many positional dofs there are
            const unsigned n_position_type = e_pt->nnodal_position_type();

            //Set up memory for the shape functions
            Shape psi(n_node, n_position_type);
            DShape dpsidxi(n_node, n_position_type, DIM);

            // Timescale ratio (non-dim density)
            double lambda_sq = e_pt->lambda_sq();

            //Loop over the integration points
            for (unsigned ipt = 0; ipt < n_intpt; ipt++) {
                //Assign local coordinate s
                for (unsigned i = 0; i < DIM; i++) { s[i] = e_pt->integral_pt()->knot(ipt, i); }

                //Get the integral weight
                double w = e_pt->integral_pt()->weight(ipt);

                //Evaluate the shape function and its derivatives, and get Jacobian
                double J = e_pt->dshape_lagrangian_at_knot(ipt, psi, dpsidxi);

                //Get mass and the coupling weight at the integration point
                double coupling_w = 0;
                //for (unsigned l = 0; l < n_node; l++)
                //{
                //    double nodal_coupling_w = dynamic_cast<SolidNode*>(e_pt->node_pt(l))->get_coupling_weight();
                //    double psi_ = psi(l);
                //    coupling_w += psi_ * nodal_coupling_w;
                //}

                //Get the coordinates of the integration point
                Vector<double> x(DIM, 0.0);
                e_pt->interpolated_x(s,x);

                //Storage for Lagrangian coordinates and velocity (initialised to zero)
                Vector<double> interpolated_xi(DIM, 0.0);
                Vector<double> veloc(DIM, 0.0);

                //Calculate lagrangian coordinates
                for (unsigned l = 0; l < n_node; l++) {
                    //Loop over positional dofs
                    for (unsigned k = 0; k < n_position_type; k++) {
                        //Loop over displacement components (deformed position)
                        for (unsigned i = 0; i < DIM; i++) {
                            //Calculate the Lagrangian coordinates
                            interpolated_xi[i] += e_pt->lagrangian_position_gen(l, k, i) * psi(l, k);

                            //Calculate the velocity components (if unsteady solve)
                            if (e_pt->is_inertia_enabled()) {
                                veloc[i] += e_pt->dnodal_position_gen_dt(l, k, i) * psi(l, k);
                            }
                        }
                    }
                }

                //Get isotropic growth factor
                double gamma = 1.0;
                e_pt->get_isotropic_growth(ipt, s, interpolated_xi, gamma);

                //Premultiply the undeformed volume ratio (from the isotropic
                // growth), the integral weights, the coupling weights, and the Jacobian
                double W = gamma * w * (1.0 - coupling_w) * J;

                DenseMatrix<double> sigma(DIM, DIM);
                DenseMatrix<double> strain(DIM, DIM);

                //Now calculate the stress tensor from the constitutive law
                e_pt->get_stress(s, sigma);

                // Add pre-stress
                for (unsigned i = 0; i < DIM; i++) {
                    for (unsigned j = 0; j < DIM; j++) {
                        sigma(i, j) += e_pt->prestress(i, j, interpolated_xi);
                    }
                }

                //get the strain
                e_pt->get_strain(s, strain);

                // Initialise
                double localElasticEnergy = 0;
                double velocitySquared = 0;

                // Compute integrals
                for (unsigned i = 0; i < DIM; i++) {
                    for (unsigned j = 0; j < DIM; j++) {
                        localElasticEnergy += sigma(i, j) * strain(i, j);
                    }
                    velocitySquared += veloc[i] * veloc[i];
                }

                // Mass
                mass += lambda_sq * W;
                // Linear momentum and angular momentum
                Vector<double> cross_product(DIM, 0);
                VectorHelpers::cross(interpolated_xi, veloc, cross_product);
                for (unsigned i = 0; i < DIM; i++) {
                    com[i] += lambda_sq * W * x[i];
                    linearMomentum[i] += lambda_sq * veloc[i] * W;
                    angularMomentum[i] += lambda_sq * cross_product[i] * W;
                }
                // Potential energy
                elasticEnergy += 0.5 * localElasticEnergy * W;
                // Kinetic energy
                kineticEnergy += lambda_sq * 0.5 * velocitySquared * W;
            }
        }
        for (unsigned i = 0; i < com.size(); i++) {
            com[i] /= mass;
        }
    }

};
#endif //SOLID_PROBLEM_H
