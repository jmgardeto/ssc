/**
BSD-3-Clause
Copyright 2019 Alliance for Sustainable Energy, LLC
Redistribution and use in source and binary forms, with or without modification, are permitted provided
that the following conditions are met :
1.	Redistributions of source code must retain the above copyright notice, this list of conditions
and the following disclaimer.
2.	Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the distribution.
3.	Neither the name of the copyright holder nor the names of its contributors may be used to endorse
or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER, CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES
DEPARTMENT OF ENERGY, NOR ANY OF THEIR EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "csp_solver_cavity_receiver.h"
#include "csp_solver_core.h"
#include "csp_solver_util.h"
#include "sam_csp_util.h"

#include "Ambient.h"
#include "definitions.h"

C_cavity_receiver::C_cavity_receiver()
{
	
}

void C_cavity_receiver::genOctCavity(double height /*m*/, double width /*m*/)
{
    //function[PANEL1, PANEL2, PANEL3, PANEL4, FLOOR, COVER, TOPLIP, BOTLIP, APERTURE] = genOctCavity(height, width, lipTop, lipBot, varargin) % #codegen
    //    % creates geometry for a half - octagonal cavity reciever of specified width
    //    % andheight, with or without upper and lower lips

    size_t nPanels = 4;
    mv_rec_surfs.resize(nPanels + 5);

    mv_rec_surfs[PANEL1].type = 1;
    mv_rec_surfs[PANEL2].type = 1;
    mv_rec_surfs[PANEL3].type = 1;
    mv_rec_surfs[PANEL4].type = 1;
    mv_rec_surfs[FLOOR].type = 0;
    mv_rec_surfs[COVER].type = 0;
    mv_rec_surfs[TOPLIP].type = 1;
    mv_rec_surfs[BOTLIP].type = 1;
    mv_rec_surfs[APERTURE].type = 2;

    mv_rec_surfs[PANEL1].is_active_surf = true;
    mv_rec_surfs[PANEL2].is_active_surf = true;
    mv_rec_surfs[PANEL3].is_active_surf = true;
    mv_rec_surfs[PANEL4].is_active_surf = true;
    mv_rec_surfs[FLOOR].is_active_surf = false;
    mv_rec_surfs[COVER].is_active_surf = false;
    mv_rec_surfs[TOPLIP].is_active_surf = false;
    mv_rec_surfs[BOTLIP].is_active_surf = false;
    mv_rec_surfs[APERTURE].is_active_surf = false;

    // hardcode top and bottom lip for now
    double lipTop = 0.0;
    double lipBot = 0.0;

    double theta = CSP::pi / (double)nPanels;

    // matrix_t(nr, nc, val)
    // (x, y, z)
    mv_rec_surfs[FLOOR].vertices.resize_fill(nPanels + 1, 3, 0.0);
    mv_rec_surfs[COVER].vertices.resize_fill(nPanels + 1, 3, 0.0);

    for (size_t i = 0; i < nPanels+1; i++) {
        mv_rec_surfs[FLOOR].vertices(i, 0) = mv_rec_surfs[COVER].vertices(i, 0) = width * cos(i * theta)/2.0;
        mv_rec_surfs[FLOOR].vertices(i, 1) = mv_rec_surfs[COVER].vertices(i, 1) = width * sin(i * theta)/2.0;
        mv_rec_surfs[FLOOR].vertices(i, 2) = -height / 2.0;
        mv_rec_surfs[COVER].vertices(i, 2) = height / 2.0;
    }

    util::matrix_t<double> temp_total(4, 3, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> temp_p1(4, 3, std::numeric_limits<double>::quiet_NaN());

    for (size_t n = 0; n < nPanels; n++) {
        temp_p1 = mv_rec_surfs[FLOOR].vertices.row(n);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            temp_total(0, i) = temp_p1(0, i);
        }
        temp_p1 = mv_rec_surfs[COVER].vertices.row(n);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            temp_total(1, i) = temp_p1(0, i);
        }
        temp_p1 = mv_rec_surfs[COVER].vertices.row(n + 1);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            temp_total(2, i) = temp_p1(0, i);
        }
        temp_p1 = mv_rec_surfs[FLOOR].vertices.row(n + 1);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            temp_total(3, i) = temp_p1(0, i);
        }
        if (n == 0) {
            mv_rec_surfs[PANEL1].vertices = temp_total;
        }
        else if (n == 1) {
            mv_rec_surfs[PANEL2].vertices = temp_total;
        }
        else if (n == 2) {
            mv_rec_surfs[PANEL3].vertices = temp_total;
        }
        else if (n == 3) {
            mv_rec_surfs[PANEL4].vertices = temp_total;
        }
    }

    mv_rec_surfs[APERTURE].vertices.resize_fill(4, 3, std::numeric_limits<double>::quiet_NaN());
    mv_rec_surfs[TOPLIP].vertices.resize_fill(4, 3, std::numeric_limits<double>::quiet_NaN());
    mv_rec_surfs[BOTLIP].vertices.resize_fill(4, 3, std::numeric_limits<double>::quiet_NaN());

    if (lipTop <= 0.0) {
        temp_p1 = mv_rec_surfs[COVER].vertices.row(0);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            mv_rec_surfs[APERTURE].vertices(0, i) = temp_p1(0, i);
        }
        temp_p1 = mv_rec_surfs[COVER].vertices.row(nPanels);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            mv_rec_surfs[APERTURE].vertices(1, i) = temp_p1(0, i);
        }
    }

    if (lipBot <= 0.0) {
        temp_p1 = mv_rec_surfs[FLOOR].vertices.row(nPanels);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            mv_rec_surfs[APERTURE].vertices(2, i) = temp_p1(0, i);
        }
        temp_p1 = mv_rec_surfs[FLOOR].vertices.row(0);
        for (size_t i = 0; i < temp_p1.length(); i++) {
            mv_rec_surfs[APERTURE].vertices(3, i) = temp_p1(0, i);
        }
    }

    return;
}


void C_cavity_receiver::meshGeometry()
{
    /*% Mesh Surfaces
        [nodesGlobal, elemsPanel1, elemsPanel2, elemsPanel3, elemsPanel4, ...
        elemsFloor, elemsCover, elemsTopLip, elemsBotLip, elemsApert] = meshGeometry(elemSizes, meshTypes, ...
            PANEL1, PANEL2, PANEL3, PANEL4, FLOOR, COVER, TOPLIP, BOTLIP, APERTURE);*/


    for (size_t i = 0; i < mv_rec_surfs.size(); i++)
    {
        util::matrix_t<double> surf = mv_rec_surfs[i].vertices;
        double elemSize = mv_rec_surfs[i].e_size;
        size_t type = mv_rec_surfs[i].type;

        if (std::isnan(surf(0, 0))) {

            // nodes
            // varargout

            continue;
        }
        else {
            if (type == 0) {
                // Mesh with triangles
                double a = 1.23;
            }
            else if (type == 1) {
                // Mesh with quads
                meshMapped(surf, elemSize);
                //[nodes{ i }, elems] = meshMapped(SURF, elemSize);
            }
            else {
                // Mesh as a single element
                double c = 4.6;
            }
        }

        // Shift node IDs to account for previous element sets
    }


    /*for i = 1:nSurfs

        if isscalar(SURF) && isnan(SURF)
            % skip this surface
            nodes{ i } = NaN;
        varargout{ i } = NaN;
        else
            if type == 0 % mesh with triangles
                [nodes{ i }, elems] = meshPolygon(SURF, elemSize);

        elseif type == 1 % mesh with quads
            [nodes{ i }, elems] = meshMapped(SURF, elemSize);

            else% mesh as single element
                nodes{ i } = SURF;
        elems = 1:size(SURF, 1);

        end

            % shift node IDs to account for previous element sets
            varargout{ i } = elems + lastNodeID;
        lastNodeID = lastNodeID + size(nodes{ i }, 1);
        end
    end*/

    return;
}

void C_cavity_receiver::meshMapped(const util::matrix_t<double>& poly, double elemSize)
{
    double almostZero = 1.E-7;
    bool geom3D = true;

    // Confirm correct inputs
    size_t n_verts = poly.nrows();
    size_t n_dims = poly.ncols();

    if (n_dims == 3) {

        // Test for coplanar vertices
        if (n_verts != 4) {
            throw(C_csp_exception("meshMapped requires 4 vertices"));
        }
        else {
            util::matrix_t<double> less1_0(1, 3, std::numeric_limits<double>::quiet_NaN());
            util::matrix_t<double> less2_0(1, 3, std::numeric_limits<double>::quiet_NaN());
            for (size_t i = 0; i < 3; i++) {
                less1_0(0, i) = poly(1, i) - poly(0, i);
                less2_0(0, i) = poly(2, i) - poly(0, i);
            }
            util::matrix_t<double> n;
            crossproduct(less1_0, less2_0, n);
            util::matrix_t<double> n_hat;
            norm3Dvect(n, n_hat);
            double volume = 0.0;
            util::matrix_t<double> diff_local(1, 3, std::numeric_limits<double>::quiet_NaN());
            for (size_t i = 0; i < 4; i++) {
                // the triple product of any combination of vertices must be zero for the polygon to be planar
                for (size_t j = 0; j < 3; j++) {
                    diff_local(0, j) = poly(i, j) - poly(0, j);
                }
                volume = std::abs(dotprod3D(n, diff_local));

                if (volume > almostZero) {
                    throw(C_csp_exception("meshMapped polygon vertices not coplanar"));
                }
            }
        }
    }
    else {
        throw(C_csp_exception("meshMapped requires 3 dimensions for a vortex"));
    }

    return;
}

double C_cavity_receiver::dotprod3D(const util::matrix_t<double>& a, const util::matrix_t<double>& b)
{
    return a(0,0)*b(0,0) + a(0,1)*b(0,1) + a(0,2)*b(0,2);
}

void C_cavity_receiver::crossproduct(const util::matrix_t<double>& a_vert, const util::matrix_t<double>& b_vert, util::matrix_t<double>& cross)
{
    cross.resize_fill(1, 3, std::numeric_limits<double>::quiet_NaN());

    cross(0,0) = a_vert(0,1)*b_vert(0,2) - a_vert(0,2)*b_vert(0,1);
    cross(0,1) = a_vert(0,2)*b_vert(0,0) - a_vert(0,0)*b_vert(0,2);
    cross(0,2) = a_vert(0,0)*b_vert(0,1) - a_vert(0,1)*b_vert(0,0);
}

void C_cavity_receiver::norm3Dvect(const util::matrix_t<double>& vector_in, util::matrix_t<double>& norm_vect)
{
    norm_vect.resize_fill(1, 3, std::numeric_limits<double>::quiet_NaN());
    double magnitude = std::sqrt(std::pow(vector_in(0,0),2) + std::pow(vector_in(0,1),2) + std::pow(vector_in(0,2),2));
    for (size_t i = 0; i < 3; i++) {
        norm_vect(0, i) = vector_in(0, i) / magnitude;
    }
}

void C_cavity_receiver::init()
{
	ambient_air.SetFluid(ambient_air.Air);

    receiverHeight = 12; // Receiver opening height in meters
    receiverWidth = 14; // Reciever opening width in meters
    //topLipHeight = 1; // Height of top lip in meters
    //botLipHeight = 1; // Height of bottom lip in meters
    e_act_sol = 0.965; // Absorbtivity in short wave range for active surfaces
    e_pass_sol = 0.05; // Absorbtivity in short wave range for passive surfaces
    e_act_therm = 0.85; // Emissivity in long wave range for active surfaces
    e_pass_therm = 0.25; // Emissivity in long wave range for passive surfaces
    //T_HTFin = 290 + 273.15; // Inlet heat transfer fluid temperature
    //T_HTFout = 575 + 273.15; // Outlet heat transfer fluid temperature
    //T_inf = 20 + 273.15; // Temperature of surroundings
    //UA_elemental = 4000; // Specified conductance from HTF to each element
    //flux_elemental = 388858.025; // Specified incident solar flux on each element
    //h = 0; // Convective heat transfer coefficients per element

    genOctCavity(receiverHeight, receiverWidth);

    meshGeometry();

	return;
}

void C_cavity_receiver::call(const C_csp_weatherreader::S_outputs& weather,
	const C_csp_solver_htf_1state& htf_state_in,
	const C_pt_receiver::S_inputs& inputs,
	const C_csp_solver_sim_info& sim_info)
{
	return;
}

void C_cavity_receiver::off(const C_csp_weatherreader::S_outputs& weather,
	const C_csp_solver_htf_1state& htf_state_in,
	const C_csp_solver_sim_info& sim_info)
{

	return;
}

void C_cavity_receiver::converged()
{
	// Check HTF props?
	//!MJW 9.8.2010 :: Call the property range check subroutine with the inlet and outlet HTF temps to make sure they're in the valid range
	//call check_htf(Coolant,T_salt_hot)
	//call check_htf(Coolant,T_salt_cold)

	if (m_mode == C_csp_collector_receiver::STEADY_STATE)
	{
		throw(C_csp_exception("Receiver should only be run at STEADY STATE mode for estimating output. It must be run at a different mode before exiting a timestep",
			"MSPT receiver converged method"));
	}

}


void C_cavity_receiver::calc_pump_performance(double rho_f, double mdot, double ffact, double& PresDrop_calc, double& WdotPump_calc)
{

	

}

double C_cavity_receiver::get_pumping_parasitic_coef()
{
	return  12.3;
}

double C_cavity_receiver::area_proj()
{
	return 1.23;
}
