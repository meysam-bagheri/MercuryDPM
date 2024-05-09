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

#ifndef SCOUPLEDELEMENT_H
#define SCOUPLEDELEMENT_H
#include "Oomph/Coupling/CoupledSolidNodes.h"

namespace oomph
{

//=================start_wrapper==================================
/// Wrapper class for solid elements to be coupled with discrete
/// solid particles on the surfaces
//================================================================
template<class ELEMENT>
class SCoupledElement : public ELEMENT
{

public:
    
    /// Constructor: Call constructor of underlying element
    SCoupledElement()
    {};
    
    /// Destructor (empty)
    ~SCoupledElement()
    {};

    /// Construct the local node n and return a pointer to it.
    /// Needs to be overridden because a SurfaceCoupledElement uses CoupledSolidNode's
    Node* construct_node(const unsigned& n)
    {
        // Construct a solid node and assign it to the local node pointer vector.
        // The dimension and number of values are taken from internal element data
        // The number of solid pressure dofs are also taken from internal data
        // The number of timesteps to be stored comes from the problem!
        this->node_pt(n) = new CoupledSolidNode(
            this->lagrangian_dimension(),
            this->nnodal_lagrangian_type(),
            this->nodal_dimension(),
            this->nnodal_position_type(),
            this->required_nvalue(n));
        // Now return a pointer to the node, so that the mesh can find it
        return this->node_pt(n);
    }

    /// Construct the local node n and return
    /// a pointer to it. Additionally, create storage for `history'
    /// values as required by timestepper
    /// Needs to be overridden because a SurfaceCoupledElement uses CoupledSolidNode's
    Node* construct_node(const unsigned& n, TimeStepper* const& time_stepper_pt)
    {
        // Construct a solid node and assign it to the local node pointer vector
        // The dimension and number of values are taken from internal element data
        // The number of solid pressure dofs are also taken from internal data
        this->node_pt(n) = new CoupledSolidNode(
            time_stepper_pt,
            this->lagrangian_dimension(),
            this->nnodal_lagrangian_type(),
            this->nodal_dimension(),
            this->nnodal_position_type(),
            this->required_nvalue(n));
        // Now return a pointer to the node, so that the mesh can find it
        return this->node_pt(n);
    }

    /// Construct the local node n and return a pointer to it.
    /// in the case when it is a boundary node; that is it MAY be
    /// located on a Mesh boundary
    /// Needs to be overridden because a SurfaceCoupledElement uses CoupledSolidNode's
    Node* construct_boundary_node(const unsigned& n)
    {
        // Construct a solid node and assign it to the local node pointer vector.
        // The dimension and number of values are taken from internal element data
        // The number of solid pressure dofs are also taken from internal data
        // The number of timesteps to be stored comes from the problem!
        this->node_pt(n) = new BoundaryNode<CoupledSolidNode>(
            this->lagrangian_dimension(),
            this->nnodal_lagrangian_type(),
            this->nodal_dimension(),
            this->nnodal_position_type(),
            this->required_nvalue(n));
        // Now return a pointer to the node, so that the mesh can find it
        return this->node_pt(n);
    }

    /// Construct the local node n and return
    /// a pointer to it, in the case when the node MAY be located
    /// on a boundary. Additionally, create storage for `history'
    /// values as required by timestepper
    /// Needs to be overridden because a SurfaceCoupledElement uses CoupledSolidNode's
    Node* construct_boundary_node(const unsigned& n,
                                   TimeStepper* const& time_stepper_pt)
    {
        // Construct a solid node and assign it to the local node pointer vector
        // The dimension and number of values are taken from internal element data
        // The number of solid pressure dofs are also taken from internal data
        this->node_pt(n) = new BoundaryNode<CoupledSolidNode>(
            time_stepper_pt,
            this->lagrangian_dimension(),
            this->nnodal_lagrangian_type(),
            this->nodal_dimension(),
            this->nnodal_position_type(),
            this->required_nvalue(n));
        // Now return a pointer to the node, so that the mesh can find it
        return this->node_pt(n);
    }


    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
        // Call the generic residuals function with flag set to 0 using a dummy matrix argument
        ELEMENT::fill_in_generic_contribution_to_residuals_pvd(
            residuals, GeneralisedElement::Dummy_matrix, 0);
        
        // Add point source contribution
        add_external_coupling_forces_to_residuals(residuals);
    }
    
    /// Add the element's contribution to its residual vector and Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
        //Call the generic routine with the flag set to 1
        ELEMENT::fill_in_generic_contribution_to_residuals_pvd(residuals, jacobian, 1);
        
        // Add point source contribution
        add_external_coupling_forces_to_residuals(residuals);
    }
    
    void get_momentum_and_energy(double& mass, Vector<double>& lin_mo, Vector<double>& ang_mo, double& pot_en,
                                 double& kin_en)
    {
        const unsigned DIM = this->dim();
        // Initialise mass
        mass = 0;
        // Initialise momentum
        lin_mo.initialise(0);
        ang_mo.initialise(0);
        // Initialise energy
        pot_en = 0;
        kin_en = 0;
        
        //Set the value of n_intpt
        unsigned n_intpt = this->integral_pt()->nweight();
        
        //Set the Vector to hold local coordinates
        Vector<double> s(DIM);
        
        //Find out how many nodes there are
        const unsigned n_node = this->nnode();
        
        //Find out how many positional dofs there are
        const unsigned n_position_type = this->nnodal_position_type();
        
        //Set up memory for the shape functions
        Shape psi(n_node, n_position_type);
        DShape dpsidxi(n_node, n_position_type, DIM);
        
        // Timescale ratio (non-dim density)
        double lambda_sq = this->lambda_sq();
        
        //Loop over the integration points
        for (unsigned ipt = 0; ipt < n_intpt; ipt++)
        {
            //Assign values of s
            for (unsigned i = 0; i < DIM; i++)
            { s[i] = this->integral_pt()->knot(ipt, i); }
            
            //Get the integral weight
            double w = this->integral_pt()->weight(ipt);
            
            //Call the derivatives of the shape functions and get Jacobian
            double J = this->dshape_lagrangian_at_knot(ipt, psi, dpsidxi);
            
            //Storage for Lagrangian coordinates and velocity (initialised to zero)
            Vector<double> interpolated_xi(DIM, 0.0);
            Vector<double> veloc(DIM, 0.0);
            
            //Calculate lagrangian coordinates
            for (unsigned l = 0; l < n_node; l++)
            {
                //Loop over positional dofs
                for (unsigned k = 0; k < n_position_type; k++)
                {
                    //Loop over displacement components (deformed position)
                    for (unsigned i = 0; i < DIM; i++)
                    {
                        //Calculate the Lagrangian coordinates
                        interpolated_xi[i] += this->lagrangian_position_gen(l, k, i) * psi(l, k);
                        
                        //Calculate the velocity components (if unsteady solve)
                        if (this->Unsteady)
                        {
                            veloc[i] += this->dnodal_position_gen_dt(l, k, i) * psi(l, k);
                        }
                    }
                }
            }
            
            //Get isotropic growth factor gamma
            double gamma = 1.0;
            this->get_isotropic_growth(ipt, s, interpolated_xi, gamma);
            
            //Premultiply the undeformed volume ratio (from the isotropic
            // growth), the integral weights, the coupling weights, and the Jacobian
            double W = gamma * w * J;
            
            DenseMatrix<double> sigma(DIM, DIM);
            DenseMatrix<double> strain(DIM, DIM);
            
            //Now calculate the stress tensor from the constitutive law
            this->get_stress(s, sigma);
            
            // Add pre-stress
            for (unsigned i = 0; i < DIM; i++)
            {
                for (unsigned j = 0; j < DIM; j++)
                {
                    sigma(i, j) += this->prestress(i, j, interpolated_xi);
                }
            }
            
            //get the strain
            this->get_strain(s, strain);
            
            // Initialise
            double local_pot_en = 0;
            double veloc_sq = 0;
            
            // Compute integrals
            for (unsigned i = 0; i < DIM; i++)
            {
                for (unsigned j = 0; j < DIM; j++)
                {
                    local_pot_en += sigma(i, j) * strain(i, j);
                }
                veloc_sq += veloc[i] * veloc[i];
            }
            
            // Mass
            mass += lambda_sq * W;
            // Linear momentum and angular momentum
            Vector<double> cross_product(DIM, 0);
            VectorHelpers::cross(interpolated_xi, veloc, cross_product);
            for (unsigned i = 0; i < DIM; i++)
            {
                lin_mo[i] += lambda_sq * veloc[i] * W;
                ang_mo[i] += lambda_sq * cross_product[i] * W;
            }
            // Potential energy
            pot_en += 0.5 * local_pot_en * W;
            // Kinetic energy
            kin_en += lambda_sq * 0.5 * veloc_sq * W;
        }
    }
    
    void set_nodal_coupling_residual(const bool& isCoupled, Vector <Vector<double>>& residual)
    {
        if (isCoupled)
        { nodal_coupling_residual = residual; }
        else
        { nodal_coupling_residual.clear(); }
    }
    
    double get_nodal_coupling_residual(const unsigned n, const unsigned d) {
        if (nodal_coupling_residual.size()>n && nodal_coupling_residual[n].size()>d)
            return nodal_coupling_residual[n][d];
        else
            return 0;
    }

private:
    
    /// Add the point source contribution to the residual vector
    void add_external_coupling_forces_to_residuals(Vector<double>& residuals)
    {
        // No further action if no coupling force
        if (nodal_coupling_residual.size() == 0) return;
        
        // Find out how many nodes there are
        const unsigned n_node = this->nnode();
        const unsigned DIM = this->dim();
        
        //find out how many positional dofs there are
        unsigned n_position_type = this->nnodal_position_type();
        
        // Set up memory for the shape/test functions
        //Shape psi(n_node);
        
        // Loop over all nodes belonging to the element
        int local_eqn = 0;
        // Loop over the nodes of the element
        for (unsigned l = 0; l < n_node; l++)
        {
            // Loop of types of dofs
            for (unsigned k = 0; k < n_position_type; k++)
            {
                // Loop over the force components
                for (unsigned i = 0; i < DIM; i++)
                {
                    // Get the local equation
                    local_eqn = this->position_local_eqn(l, k, i);
                    // IF it's not a boundary condition and the interpolated force is nonzero
                    if (local_eqn >= 0)
                    {
                        // add the nodal coupling residual to the global residual vector
                        residuals[local_eqn] += nodal_coupling_residual[l][i];
                    }
                }
            }
        }
        logger(VERBOSE, "Apply nodal_coupling_residual element %", this);
    }
    
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
        const unsigned DIM = this->dim();
        Vector<double> x(DIM);
        Vector<double> dxdt(DIM);
        Vector<double> s(DIM);
        
        // Tecplot header info
        outfile << this->tecplot_zone_string(n_plot);
        
        // Loop over plot points
        unsigned num_plot_points = this->nplot_points(n_plot);
        for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
        {
            // Get local coordinates of plot point
            this->get_s_plot(iplot, n_plot, s);
            
            // Get Eulerian coordinates
            this->interpolated_x(s, x);
            SolidFiniteElement* el_pt = dynamic_cast<SolidFiniteElement*>(this);
            el_pt->interpolated_dxdt(s, 1, dxdt);
            
            // Dummy integration point
            unsigned ipt = 0;
            
            // Output x,y,z
            for (unsigned i = 0; i < DIM; i++)
            {
                outfile << x[i] << " ";
            }
            
            // Output velocity dnodal_position_gen_dt
            for (unsigned i = 0; i < DIM; i++)
            {
                outfile << dxdt[i] << " ";
            }
            
            outfile << std::endl;
        }
        
        // Write tecplot footer (e.g. FE connectivity lists)
        this->write_tecplot_zone_footer(outfile, n_plot);
        outfile << std::endl;
    }
    
    /// Nodal coupling forces
    Vector <Vector<double>> nodal_coupling_residual;
    
};


//===========start_face_geometry==============================================
/// FaceGeometry of wrapped element is the same as the underlying element
//============================================================================
template<class ELEMENT>
class FaceGeometry<SCoupledElement<ELEMENT> > :
    public virtual FaceGeometry<ELEMENT>
{
public:
    
    /// \short Constructor [this was only required explicitly
    /// from gcc 4.5.2 onwards...]
    FaceGeometry() : FaceGeometry<ELEMENT>()
    {}
    
};

}

#endif //SCOUPLEDELEMENT_H
