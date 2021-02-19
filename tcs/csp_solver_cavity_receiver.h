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

#ifndef __csp_solver_cavity_receiver_
#define __csp_solver_cavity_receiver_

#include "csp_solver_util.h"

#include "htf_props.h"
#include "csp_solver_core.h"
#include "csp_solver_pt_receiver.h"

class C_cavity_receiver : public C_pt_receiver
{
public:

    class C_rec_surface
    {
    public:
        util::matrix_t<double> vertices;    // (nr, nc) -> (vertices, xyx)
        double e_size;
        size_t type;
        bool is_active_surf;

        C_rec_surface()
        {
            e_size = 2.5;
        }
    };

    enum surf_order
    {
        PANEL1,
        PANEL2,
        PANEL3,
        PANEL4,
        FLOOR,
        COVER,
        TOPLIP,
        BOTLIP,
        APERTURE
    };

private:

    double receiverHeight; //[m] Receiver opening height in meters
    double receiverWidth; //[m] Reciever opening width in meters
    //double topLipHeight;  //[m] Height of top lip in meters
    //double botLipHeight;  //[m] Height of bottom lip in meters
    double e_act_sol;     //[-] Absorbtivity in short wave range for active surfaces
    double e_pass_sol;    //[-] Absorbtivity in short wave range for passive surfaces
    double e_act_therm;   //[-] Emissivity in long wave range for active surfaces
    double e_pass_therm;  //[-] Emissivity in long wave range for passive surfaces
    //double T_HTFin;       // Inlet heat transfer fluid temperature
    //double T_HTFout;      // Outlet heat transfer fluid temperature
    //double T_inf;         // Temperature of surroundings
    //double UA_elemental;  // Specified conductance from HTF to each element
    //double flux_elemental;// Specified incident solar flux on each element
    //double h;             // Convective heat transfer coefficients per element

    std::vector<C_rec_surface> mv_rec_surfs;

public:

	// Methods
	C_cavity_receiver();

	~C_cavity_receiver() {};

	virtual void init();

	virtual void call(const C_csp_weatherreader::S_outputs& weather,
		const C_csp_solver_htf_1state& htf_state_in,
		const C_pt_receiver::S_inputs& inputs,
		const C_csp_solver_sim_info& sim_info);

	virtual void off(const C_csp_weatherreader::S_outputs& weather,
		const C_csp_solver_htf_1state& htf_state_in,
		const C_csp_solver_sim_info& sim_info);

	virtual void converged();

	void calc_pump_performance(double rho_f, double mdot, double ffact, double& PresDrop_calc, double& WdotPump_calc);

	virtual double get_pumping_parasitic_coef();

	virtual double area_proj();

    void genOctCavity(double height /*m*/, double width /*m*/);

    void meshGeometry();

    void meshMapped(const util::matrix_t<double>& poly, double elemSize);

    void crossproduct(const util::matrix_t<double>&, const util::matrix_t<double>&, util::matrix_t<double>& cross);

    void norm3Dvect(const util::matrix_t<double>&, util::matrix_t<double>& norm_vect);

    double dotprod3D(const util::matrix_t<double>&, const util::matrix_t<double>&);

    void sumcolumns(const util::matrix_t<double>&, util::matrix_t<double>&);

    void ave_columns(const util::matrix_t<double>&, util::matrix_t<double>&);

    void to2D(const util::matrix_t<double>& poly, const util::matrix_t<double>& center,
        const util::matrix_t<double>& normal, const util::matrix_t<double>& xaxis);
};

#endif // __csp_solver_cavity_receiver_

