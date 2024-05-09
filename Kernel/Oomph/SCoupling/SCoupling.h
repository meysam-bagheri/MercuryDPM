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

#ifndef SURFACE_COUPLING_H
#define SURFACE_COUPLING_H

#include "Oomph/Coupling/BaseCoupling.h"
#include "SCoupledElement.h"
#include "Logger.h"
#include "chrono"
using namespace oomph;

template<class M, class O>
class SCoupling : public BaseCoupling<M, O>
{
public:
    typedef typename O::ELEMENT ELEMENT;

    /**
     * information needed to link surface-coupled oomph elements to DPM
     * \todo should we move this to SurfaceCoupledElement.h?
     */
    struct SCoupledElement
    {
        // pointer to the oomph bulk element
        ELEMENT* bulk_elem_pt;
        // index of the face on the coupled boundary
        int face_index;
        // local coordinates of a traction element's nodes
        Vector<Vector<double*> > surface_node_position_pt;
        // pointer to solid nodes of the traction element
        Vector<CoupledSolidNode*> node_pt;
        // local coordinates of a traction element's center
        Vector<double> center_local;
    };

    // default constructor
    SCoupling() = default;
    
//    // sets a cg width before solving a surface coupled OomphMercuryProblem
//    void solveSurfaceCoupling(unsigned nStep, const double& width)
//    {
//        // set the coarse-grainning width w.r.t. length scale
//        setCGWidth(width); // dimensional
//        // solve surface coupled OomphMercuryProblem
//        solveSurfaceCoupling(nStep);
//    }
    
    // solve surface coupled OomphMercuryProblem
    void solveSurfaceCoupling()
    {
        // compute nStep
        unsigned nStep = O::getOomphTimeStep()/M::getTimeStep();
        if (nStep==0) {
            // if oomph has a smaller time step than Mercury
            logger(INFO, "Set nStep %, change mercuryTimeStep from % to %",
                   nStep, M::getTimeStep(), O::getOomphTimeStep());
            nStep = 1;
            M::setTimeStep(O::getOomphTimeStep());
        } else {
            logger(INFO, "Set nStep %, change oomphTimeStep from % to %",
                   nStep, O::getOomphTimeStep(), nStep * M::getTimeStep());
            O::setOomphTimeStep(nStep * M::getTimeStep());
        }

        // call solve routine
        solveSurfaceCoupling(nStep);
    }
    
    // solve surface coupled OomphMercuryProblem
    void solveSurfaceCoupling(unsigned nStep)
    {
        // check whether time steps are set
        logger.assert_always(O::getOomphTimeStep()>0,"Oomph time step not initialised");
        logger.assert_always(M::getTimeStep()>0,"Mercury time step not initialised");
        // check whether isCoupled is set
        logger.assert_always(!coupledBoundaries_.empty(), "isCoupled needs to be set, e.g. via setIsCoupled([](unsigned b) { return b == Boundary::Y_MIN; })");

        // first part of the Mercury solve routine, containing the actions before time-stepping
        M::initialiseSolve();

        // read Mercury time step
        double mercury_dt = M::getTimeStep();
        logger(INFO, "Mercury time step: %", mercury_dt);
        
        // set oomph time step
        double oomph_dt = nStep * mercury_dt;
        logger(INFO, "Oomph time step: %", oomph_dt);
        
        // set oomph initial conditions
        O::assign_initial_values_impulsive(oomph_dt);
        
        // read oomph-mesh
        logger(INFO, "Set up oomph mesh: % elements (% with traction)",
               O::solid_mesh_pt()?O::solid_mesh_pt()->nelement():0, O::traction_mesh_pt()?O::traction_mesh_pt()->nelement():0);
        
        // get list of bulk elements along the surface-coupled boundaries
        getSCoupledElements();
        
        // create DPM triangle walls from bulk finite elements
        createDPMWallsFromFiniteElems();
        
        // this is the main loop advancing over time
        unsigned nDone = 0; //< last written file number
        while (M::getTime() < M::getTimeMax())
        {
            this->actionsBeforeOomphTimeStep();
            // solve the coupled problem for one time step
            computeOneTimeStepForSCoupling(nStep);
            // write outputs of the oomphProb; this is slaved to the vtk output of Mercury, i.e. an oomph-lib output get written everytime a Mercury vtk file gets written
            if (M::getParticlesWriteVTK() && M::getVtkWriter()->getFileCounter() > nDone)
            {
                O::writeToVTK();
                nDone = M::getVtkWriter()->getFileCounter();
            }
        }

        // close output files of mercuryProb
        M::finaliseSolve();
    }

    /**
     * Solves an unsteady problem, returns successful if timeMaxMin has been reached
     */
    void solveSurfaceCouplingForgiving(unsigned nStep, double timeMaxMin=-constants::inf) {
        // solve
        try {
            solveSurfaceCoupling(nStep);
        } catch(OomphLibError& error)  {
            //Store output if newton solver fails
            O::saveSolidMesh();
            M::finaliseSolve();
            double time = O::time_stepper_pt()->time() - nStep * M::getTimeStep();;
            double timeMax = M::getTimeMax();;
            if (time >= timeMaxMin) {
                // take it as successful if a fraction of the time evolution has finished
                logger(INFO,"Newton solver failed at t=% (tMax=%), but code will continue.", time, timeMax);
                exit(0);
            } else {
                logger(ERROR,"Newton solver failed at t=% (tMax=%).", time, timeMax);
            }
        }
    }

    // solve OomphMercuryProblem, but with fixed solid
    void solveSurfaceCouplingFixedSolid()
    {
        // first part of the Mercury solve routine, containing the actions before time-stepping
        M::initialiseSolve();
        logger.assert_always(!coupledBoundaries_.empty(), "isCoupled needs to be set, e.g. via setIsCoupled([](unsigned b) { return b == Boundary::Y_MIN; })");
        
        // read Mercury time step
        double mercury_dt = M::getTimeStep();
        logger(INFO, "Mercury time step: %", mercury_dt);
        
        // set oomph time step
        logger(INFO, "Solid position fixed");
        
        // set oomph_dt
        this->time_pt()->initialise_dt(0);
        // By default do a non-impulsive start and provide initial conditions
        this->assign_initial_values_impulsive(0);

        // read oomph-mesh
        logger(INFO, "Set up oomph mesh: % elements (% with traction)",
               O::solid_mesh_pt()?O::solid_mesh_pt()->nelement():0, O::traction_mesh_pt()?O::traction_mesh_pt()->nelement():0);
        
        // get list of bulk elements along the surface-coupled boundaries
        getSCoupledElements();
        
        // create DPM triangle walls from bulk finite elements
        createDPMWallsFromFiniteElems();
        
        // this is the main loop advancing over time
        unsigned nDone = 0; //< last written file number
        while (M::getTime() < M::getTimeMax())
        {
            this->actionsBeforeOomphTimeStep();
            M::computeOneTimeStep();
            //if (getParticlesWriteVTK() && getVtkWriter()->getFileCounter() > nDone) {
            //    writeToVTK();
            //    nDone = getVtkWriter()->getFileCounter();
            //}
        }
        // close output files of mercuryProb
        M::finaliseSolve();
    }

    /**
     * Create a triangle wall from a set of vertices
     * \param vertex
     * \todo Why we need to set GroupID. I don't remember anymore.
     * \return
     */
    TriangleWall* createTriangleWall(std::array<Vec3D, 3> vertex)
    {
        TriangleWall wall;
        auto species = M::speciesHandler.getObject(0);
        wall.setSpecies(species);
        wall.setGroupId(100);
        wall.setVertices(vertex[0], vertex[1], vertex[2]);
        auto w = M::wallHandler.copyAndAddObject(wall);
        return w;
    }

    /**
     * Update the position of a triangle wall (used in updateDPMWallsFromFiniteElems)
     * \todo why is it using setPrescribedPosition
     * \todo we need to set velocity
     */
    void updateTriangleWall(TriangleWall*& wall, std::array<Vec3D, 3> vertex)
    {
        double time0 = M::getTime();
        double dTime = O::getOomphTimeStep();
        std::array<Vec3D,3> vertex0 = wall->getVertices();
        std::array<Vec3D,3> dVertex = {
            vertex[0] - vertex0[0],
            vertex[1] - vertex0[1],
            vertex[2] - vertex0[2]};
        wall->setPrescribedPosition( [time0, dTime, vertex0, dVertex, wall]( double time ) {
            double f = ( time - time0 ) / dTime;
            std::array<Vec3D, 3> vertex = {
                vertex0[0] + dVertex[0] * f,
                vertex0[1] + dVertex[1] * f,
                vertex0[2] + dVertex[2] * f };
            wall->setVertices( vertex[0], vertex[1], vertex[2] );
            //logger(INFO,"p %",vertex[0]);
            return wall->getPosition();
        } );
        Vec3D velocity = (dVertex[0]+dVertex[1]+dVertex[2])/3./dTime;
        wall->setPrescribedVelocity([velocity] (double time) {
            //logger(INFO,"v %",velocity);
            return velocity;
        });
    }

    /**
     * Solve surface-coupled problem for one oomph time step (n mercury time steps)
     */
    void computeOneTimeStepForSCoupling(const unsigned& nStepsMercury)
    {
        auto t0 = std::chrono::system_clock::now();
        updateDPMWallsFromFiniteElems();
        auto t1 = std::chrono::system_clock::now();
        BaseCoupling<M,O>::solveMercury(nStepsMercury);
        auto t2 = std::chrono::system_clock::now();
        if (solidFeelsParticles_) {
            updateTractionOnFiniteElems();
        }
        auto t3 = std::chrono::system_clock::now();
        BaseCoupling<M,O>::solveOomph();
        auto t4 = std::chrono::system_clock::now();
        if (logSurfaceCoupling) logger(INFO, "time % Elapsed time: FEM->DEM %, DEM %, DEM->FEM %, FEM %", M::getTime(),
               (t1-t0).count(), (t2-t1).count(), (t3-t2).count(), (t4-t3).count());
    }
    
    void createDPMWallsFromFiniteElems()
    {
        // loop over bulk elements at boundary
        for (auto sCoupledElement : sCoupledElements_)
        {
            // loop over nodes in element at boundary
            Vector<Vector<double*> > position_pt = sCoupledElement.surface_node_position_pt;
            // reordering vertices of oomph face element (default = {0,1,3,2})
            swap(position_pt[2], position_pt[3]);
            
            // get global coordinate at the center
            Vec3D center;
            Vector<double> x(3, 0.0);
            sCoupledElement.bulk_elem_pt->interpolated_x(sCoupledElement.center_local, x);
            center.setX(x[0]);
            center.setY(x[1]);
            center.setZ(x[2]);
            
            // create TriangleWalls from oomph face element
            const unsigned nTriangles = position_pt.size();
            unsigned n = 0;
            while (n < nTriangles)
            {
                // get vertices of TriangleWall (multiply vertex position with the length scale of the O)
                std::array<Vec3D, 3> vertex;
                // one vertex at the center
                vertex[0] = center;
                // two vertices from the O<element,TIMESTEPPER>
                vertex[1] = Vec3D(*position_pt[0][0],*position_pt[0][1],*position_pt[0][2]);
                vertex[2] = Vec3D(*position_pt[1][0],*position_pt[1][1],*position_pt[1][2]);
                
                // create triangle facet
                TriangleWall* w = createTriangleWall(vertex);
                triangleWalls_.push_back(w);
                
                // rotate forward by one element
                rotate(position_pt.begin(), position_pt.begin() + 1, position_pt.end());
                n++;
            }
        }
    }
    
    void updateDPMWallsFromFiniteElems()
    {
        // loop over bulk elements at boundary
        unsigned wallID = 0;
        for (auto sCoupledElement : sCoupledElements_)
        {
            // get members of a SCoupledElement
            Vector<Vector<double*> > position_pt = sCoupledElement.surface_node_position_pt;
            Vector<CoupledSolidNode*> node_pt = sCoupledElement.node_pt;
            
            // reordering vertices of oomph face element (default = {0,1,3,2})
            // \todo can we do this swap in getSCoupledElements, then we don't have to create a local copy position_pt all the time?
            swap(position_pt[2], position_pt[3]);
            
            // get global coordinate at the center
            Vector<double> x(3, 0.0);
            sCoupledElement.bulk_elem_pt->interpolated_x(sCoupledElement.center_local, x);
            Vec3D center {x[0], x[1], x[2]};
            
            // get number of TriangleWalls per oomph face element
            const unsigned nTriangles = position_pt.size();
            unsigned n = 0;
            while (n < nTriangles)
            {
                // get vertices of TriangleWall (multiply vertex position with the length scale of the O)
                std::array<Vec3D, 3> vertex;
                // one vertex at the center
                vertex[0] = center;
                // two vertices from the O<element,TIMESTEPPER>
                vertex[1] = Vec3D(*position_pt[0][0],
                                  *position_pt[0][1],
                                  *position_pt[0][2]);;
                vertex[2] = Vec3D(*position_pt[1][0],
                                  *position_pt[1][1],
                                  *position_pt[1][2]);
                
                // update vertices of triangle facet  (multiply vertex position with the length scale of the O<element,TIMESTEPPER>)
                updateTriangleWall(triangleWalls_[wallID], vertex);
                
                // rotate forward by one element
                rotate(position_pt.begin(), position_pt.begin() + 1, position_pt.end());
                rotate(node_pt.begin(), node_pt.begin() + 1, node_pt.end());
                n++;
                wallID++;
            }
        }
    }
    
    /**
     * sets nodal_coupling_residual in each scoupled element
     */
    void updateTractionOnFiniteElems()
    {
        // if construct mapping with FEM basis functions
        if (!BaseCoupling<M,O>::useCGMapping())
        {
            // tracks the id of the triangle walls (get incremented in computeSCouplingForcesFromTriangles)
            unsigned wallID = 0;
            
            // loop over scoupled elements
            for (auto sCoupledElement : sCoupledElements_)
            {
                // set up memory for nodal coupling force
                Vector<Vector<double> > nodalCouplingForces;
                // returns whether the element is coupled to particles
                bool elemIsCoupled = computeSCouplingForcesFromTriangles(sCoupledElement.bulk_elem_pt,
                                                                         sCoupledElement.surface_node_position_pt.size(),
                                                                         wallID,
                                                                         nodalCouplingForces);
                
                // assign nodal coupling force to the element to be used by element::fill_in_contribution_to_residuals(...)
                sCoupledElement.bulk_elem_pt->set_nodal_coupling_residual(elemIsCoupled, nodalCouplingForces);
            }
            logger(VERBOSE, "Update nodal_coupling_residual");
        }
        else
        {
//            // how many bulk elements in total
//            unsigned n_element = O::solid_mesh_pt()->nelement();
//
//            // loop over the bulk elements
//            for (unsigned e = 0; e < n_element; e++)
//            {
//                // get pointer to the bulk element
//                ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(O::solid_mesh_pt()->element_pt(e));
//
//                // set up memory for nodal coupling force
//                Vector<Vector<double> > nodalCouplingForces;
//                // whether the element is coupled to particles
//                bool elemIsCoupled = computeSCouplingForcesFromCG(elem_pt, nodalCouplingForces);
//
//                // assign nodal coupling force to the element to be used by element::fill_in_contribution_to_residuals(...)
//                elem_pt->set_nodal_coupling_residual(elemIsCoupled, nodalCouplingForces);
//            }
//            // reset the sum of the evaluated CG values
//            for (auto inter : M::interactionHandler)
//            {
//                if (inter->getI()->getName() == "TriangleWall")
//                {
////                    logger(INFO, "Sum of the CG values previously evaluated for this interaction %",
////                           inter->getPreviousTotalCGEval());
//                    inter->resetTotalCGEval();
//                }
//            }
        }
    }
    
    /**
     * Computes nodal scoupling forces from triangles.
     * \param[in] elem_pt pointer to the element for which forces are computed
     * \param[in] nTriangles number of triangles (equals number of nodes of the traction element)
     * \param[in,out] wallID
     * \param[out] nodalCouplingForces[l] = sum_i f * psi(l,cp) for each shape function l, interaction i, evaluated at contact point cp
     * \return
     */
    bool computeSCouplingForcesFromTriangles(ELEMENT* const elem_pt, const unsigned& nTriangles,
                                             unsigned& wallID, Vector<Vector<double> >& nodalCouplingForces)
    {
        // whether the element is coupled to particles
        bool elemIsCoupled = false;
        
        // get number of nodes in the bulk element
        const unsigned nnode = elem_pt->nnode();
        // get dimension of the problem
        const unsigned dim = elem_pt->dim();
        // initialize the coupling force vector
        nodalCouplingForces.resize(nnode, Vector<double>(dim, 0.0));
        
        // get number of TriangleWalls created from an oomph face element
        unsigned n = 0;
        // loop over TriangleWalls
        while (n++ < nTriangles)
        {
            // get pointer to the wall
            TriangleWall* w = triangleWalls_[wallID++];
            
            // skip the wall if not interactions with it
            if (w->getInteractions().size() == 0) continue;
            
            // Set up memory for the shape/test functions
            Shape psi(nnode);
            
            // loop over interactions with the wall
            for (auto inter : w->getInteractions())
            {
                if (!inter->getForce().isZero())
                {
                    // scale contact point and forces from DEM to FEM units
                    Vec3D xc = inter->getContactPoint();
                    Vec3D fc = inter->getForce();
                    
                    Vector<double> x(3, 0.0), f(3, 0.0);
                    x[0] = xc.getX();
                    x[1] = xc.getY();
                    x[2] = xc.getZ();
                    f[0] = fc.getX();
                    f[1] = fc.getY();
                    f[2] = fc.getZ();
                    
                    // get the local coordinate s if xc is located in the finite element
                    Vector<double> s(3, 0.0);
                    GeomObject* geom_obj_pt = 0;
                    elem_pt->locate_zeta(x, geom_obj_pt, s);
                    
                    // Get shape/test fcts
                    elem_pt->shape(s, psi);
                    // Loop over the test functions
                    for (unsigned l = 0; l < nnode; l++)
                    {
                        //Loop over the force components
                        for (unsigned i = 0; i < dim; i++)
                        {
                            // add contribution to the nodal coupling force
                            nodalCouplingForces[l][i] += f[i] * psi(l);
                        }
                    }
                    // set the flag to true
                    elemIsCoupled = true;
                }
            }
        }
        return elemIsCoupled;
    }
    
//    /**
//     * Computes nodal scoupling forces from triangles using coarse graining
//     * \param[in] elem_pt pointer to the element for which forces are computed
//     * \param[in,out] nodalCouplingForces coupling forces at nodal positions to be added to the residual
//     * \return
//     */
//    bool computeSCouplingForcesFromCG(ELEMENT*& elem_pt, Vector<Vector<double> >& nodalCouplingForces)
//    {
//        // whether the element is coupled to particles
//        bool elemIsCoupled = false;
//
//        // get number of nodes in the bulk element
//        const unsigned nnode = elem_pt->nnode();
//        // get dimension of the problem
//        const unsigned dim = elem_pt->dim();
//        // initialize the coupling force vector
//        nodalCouplingForces.resize(nnode, Vector<double>(dim, 0.0));
//
//        // get particles if there are in the bounding box
//        Vec3D min, max;
//        getElementBoundingBox(elem_pt, min, max);
//        Vector<BaseParticle*> pList;
//        BaseCoupling<M,O>::getParticlesInCell(min, max, pList);
//
//        // loop over shape functions at the contact points
//        for (const auto p : pList)
//        {
//            for (const auto inter : p->getInteractions())
//            {
//                if (inter->getI()->getName() == "TriangleWall")
//                {
//                    // get the number of integration points
//                    const unsigned n_intpt = elem_pt->integral_pt()->nweight();
//
//                    // loop over the integration points
//                    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
//                    {
//                        // set up memory for the local and global coordinates
//                        Vector<double> s(dim, 0.0), x(dim, 0.0);
//
//                        // get the value of the local coordinates at the integration point
//                        for (unsigned i = 0; i < dim; i++) s[i] = elem_pt->integral_pt()->knot(ipt, i);
//
//                        // get the value of the global coordinates at the integration point
//                        elem_pt->interpolated_x(s, x);
//
//                        // set CG coordinates
//                        CGCoordinates::XYZ coordinate;
//                        coordinate.setXYZ(Vec3D(x[0], x[1], x[2]));
//                        // evaluate the value of CG function around particle m at CGcoords \phi(\vec r_i-r_m)
//                        double phi = BaseCoupling<M,O>::getCGFunction().evaluateCGFunction(
//                            inter->getContactPoint(), coordinate);
//
//                        // add contributions to the coupling force
//                        if (!inter->getForce().isZero() && phi > 0.0)
//                        {
//                            // set the flag to true
//                            elemIsCoupled = true;
//
//                            // Set up memory for the shape/test functions
//                            Shape psi(nnode);
//
//                            // get the integral weight
//                            double w = elem_pt->integral_pt()->weight(ipt);
//                            // find the shape functions at the integration points r_i
//                            elem_pt->shape_at_knot(ipt, psi);
//
//                            // loop over the nodes
//                            for (unsigned l = 0; l < nnode; l++)
//                            {
//                                // CG mapping defined as \tilde{N_{l,m}}_ipt = w_ipt * \phi(\vec r_i-r_m) * N_l(r_i)
//                                double shape = w * phi * psi(l);
//                                inter->addCGEval(shape);
//                                shape /= inter->getPreviousTotalCGEval();
//                                Vec3D fc = inter->getForce() * shape;
//                                nodalCouplingForces[l][0] += fc.getX();
//                                nodalCouplingForces[l][1] += fc.getY();
//                                nodalCouplingForces[l][2] += fc.getZ();
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        return elemIsCoupled;
//    }
    
    /**
     * needed in computeSCouplingForcesFromCG to decide which particles affect a given element
     */
    void getElementBoundingBox(ELEMENT*& elem_pt, Vec3D& min, Vec3D& max)
    {
        // three arrays that contain the x, y and z coordinates of the bulk element
        Vector<double> listOfCoordX;
        Vector<double> listOfCoordY;
        Vector<double> listOfCoordZ;
        
        // get the x, y and z coordinates of the bulk element
        const unsigned nnode = elem_pt->nnode();
        for (unsigned n = 0; n < nnode; n++)
        {
            listOfCoordX.push_back(elem_pt->node_pt(n)->x(0));
            listOfCoordY.push_back(elem_pt->node_pt(n)->x(1));
            listOfCoordZ.push_back(elem_pt->node_pt(n)->x(2));
        }
        
        // get the bounding box of the bulk element
        min.X = *min_element(listOfCoordX.begin(), listOfCoordX.end());
        min.Y = *min_element(listOfCoordY.begin(), listOfCoordY.end());
        min.Z = *min_element(listOfCoordZ.begin(), listOfCoordZ.end());
        max.X = *max_element(listOfCoordX.begin(), listOfCoordX.end());
        max.Y = *max_element(listOfCoordY.begin(), listOfCoordY.end());
        max.Z = *max_element(listOfCoordZ.begin(), listOfCoordZ.end());
        
        // extend the bounding box if construct mapping with coarse graining
        logger.assert_always(M::particleHandler.getLargestParticle(), "No particles detected");
        //todo does it make sense with -2R here?
        min -= Vec3D(1.0, 1.0, 1.0) * ( BaseCoupling<M,O>::getCGWidth() - 2 * M::particleHandler.getLargestParticle()->getRadius() );
        max += Vec3D(1.0, 1.0, 1.0) * ( BaseCoupling<M,O>::getCGWidth() - 2 * M::particleHandler.getLargestParticle()->getRadius() );
    }
    
    /**
     * Get bulk elements along boundaries (SCoupling).
     * The results are stored in sCoupledElements_.
     */
    void getSCoupledElements()
    {
        // first clear bulk elements (for refineable problem)
        sCoupledElements_.clear();
        
        // loop over all boundaries
        for (unsigned b : coupledBoundaries_)
        {
            // we only need to couple the elements on the upper boundary
            logger(INFO,"Coupling boundary %", b);
            
            // number of bulk elements adjacent to boundary b
            unsigned n_element = this->solid_mesh_pt()->nboundary_element(b);
            
            // loop over the bulk elements adjacent to boundary b
            for (unsigned e = 0; e < n_element; e++)
            {
                // initialize the struct SCoupledElement
                SCoupledElement sCoupledElement;
                
                // get pointer to the bulk element that is adjacent to boundary b
                sCoupledElement.bulk_elem_pt = dynamic_cast<ELEMENT*>(this->solid_mesh_pt()->boundary_element_pt(b, e));
                
                // get the index of the face of element e along boundary b
                sCoupledElement.face_index = this->solid_mesh_pt()->face_index_at_boundary(b, e);
                
                // create temporary traction element
                SolidTractionElement<ELEMENT> trac_elem(sCoupledElement.bulk_elem_pt, sCoupledElement.face_index);
                
                // number of nodes on traction element
                unsigned n_node = trac_elem.nnode();
                
                // allocate space to store the local coordinates of the traction element's vertices
                sCoupledElement.surface_node_position_pt.reserve(n_node);
                
                // store the local coordinates of the traction element's center
                sCoupledElement.center_local.resize(3, 0.0);
                
                // assign addresses of nodal coordinates to the struct SCoupledElement
                for (unsigned n = 0; n < n_node; n++)
                {
                    // store pointer to solid nodes
                    sCoupledElement.node_pt.push_back(dynamic_cast<CoupledSolidNode*>(trac_elem.node_pt(n)));
                    
                    Vector<double> s(2, 0.0);
                    trac_elem.local_coordinate_of_node(n, s);
                    
                    Vector<double> s_bulk(3, 0.0);
                    trac_elem.get_local_coordinate_in_bulk(s, s_bulk);
                    
                    // pass the addresses of x coordinates at boundary to x_pt
                    Vector<double*> x_pt;
                    for (unsigned i = 0; i < 3; i++)
                    {
                        x_pt.push_back(&trac_elem.node_pt(n)->x(i));
                        sCoupledElement.center_local[i] += s_bulk[i];
                    }
                    /*
                                        // debug that the addresses of x coordinates are passed to surface_node_position_pt
                                        cout << "Address of x(1) of node on traction element\t " << &trac_elem.node_pt(n)->x(0) << endl;
                                        int nodeIDOfBulkElement = sCoupledElement.bulk_elem_pt->get_node_number(trac_elem.node_pt(n));
                                        cout << "Address of x(1) of the same node in bulk element "
                                             << &sCoupledElement.bulk_elem_pt->node_pt(nodeIDOfBulkElement)->x(0) << endl;
                                        cout << "Address of x(1) saved in Oomph interface object\t " << x_pt[0] << endl;
                    */
                    // save the addresses of Lagrange coordinates in each element at boundary
                    sCoupledElement.surface_node_position_pt.push_back(x_pt);
                }
                
                // local coordinate of the center of the face element
                for (unsigned i = 0; i < 3; i++) sCoupledElement.center_local[i] /= n_node;
                
                // save bulk elements and surface nodal position pointers
                sCoupledElements_.push_back(sCoupledElement);
            }
        }
    }

    // set is_pinned such that a certain boundary is pinned
    void coupleBoundary(unsigned b) {
        coupledBoundaries_ = {b};
    }

    // set is_pinned such that a certain boundary is pinned
    void coupleBoundaries(std::vector<unsigned> b) {
        coupledBoundaries_ = std::move(b);
    }

    void disableLogSurfaceCoupling() {
        logSurfaceCoupling = false;
    }

    void setSolidFeelsParticles(bool val) {
        solidFeelsParticles_ = val;
    }

    bool getSolidFeelsParticles() const {
        return solidFeelsParticles_;
    }

private:
    
    /// List of surface-coupled elements
    Vector<SCoupledElement> sCoupledElements_;
    
    /// List of triangle walls used for the surface coupling
    std::vector<TriangleWall*> triangleWalls_;
    
    /**
     * Function to determine whether boundary is coupled.
     * Needs to be set before solveSurfaceCoupling is called.
     */
    std::vector<unsigned> coupledBoundaries_;

    bool logSurfaceCoupling = true;

    /**
     * Set false for one-way coupling (solid does not feel particles), true for two-way coupling
     */
    bool solidFeelsParticles_ = true;
};

#endif //SURFACE_COUPLING_H
