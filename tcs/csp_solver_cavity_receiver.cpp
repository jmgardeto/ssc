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
#include <math.h>

bool sort_pair_ascending(pair<double,double> i, pair<double, double> j)
{
    if (i.first > j.first) {
        return true;
    }
    else if (i.first == j.first && i.second > j.second) {
        return true;
    }
    else {
        return false;
    }
}

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
    mv_rec_surfs[FLOOR].type = 2;
    mv_rec_surfs[COVER].type = 2;
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

    vector<util::matrix_t<double>> v_nodes;

    int lastNodeID = 0;

    for (size_t i = 0; i < mv_rec_surfs.size(); i++)
    {
        util::matrix_t<double> surf = mv_rec_surfs[i].vertices;
        double elemSize = mv_rec_surfs[i].e_size;
        size_t type = mv_rec_surfs[i].type;

        util::matrix_t<double> nodes;
        util::matrix_t<double> elems;
        if (std::isnan(surf(0, 0))) {

            nodes.resize_fill(1,1,std::numeric_limits<double>::quiet_NaN());
            elems.resize_fill(1,1,std::numeric_limits<double>::quiet_NaN());
            v_nodes.push_back(nodes);
        }
        else {
            if (type == 0) {
                // Mesh with triangles
                throw(C_csp_exception("Triangle meshes not currently supported. Need to add qhull project"));
                meshPolygon(surf, elemSize);
            }
            else if (type == 1) {
                // Mesh with quads
                meshMapped(surf, elemSize, nodes, elems);
                v_nodes.push_back(nodes);
            }
            else {
                // Mesh as a single element
                v_nodes.push_back(surf);
                elems.resize(surf.nrows(), 1);
                for (size_t i = 0; i < surf.nrows(); i++) {
                    elems(i,0) = i;
                }
            }

            // Shift node IDs to account for previous element sets
            add_constant_to_each_element(lastNodeID, elems);
            m_v_elems.push_back(elems);
            lastNodeID = lastNodeID + v_nodes[i].nrows();
        }
    }

    m_nodesGlobal.resize_fill(lastNodeID, 3, 0.0);
    lastNodeID = 0;
    int nodeCount = 0;
    for (size_t k = 0; k < mv_rec_surfs.size(); k++) {
        if ( std::isfinite(v_nodes[k](0,0) )) { 
            nodeCount = v_nodes[k].nrows();
            for (size_t i = 0; i < nodeCount; i++) {
                for (size_t j = 0; j < 3; j++) {
                    m_nodesGlobal(lastNodeID + i,j) = v_nodes[k](i,j);
                }
            }
            lastNodeID += nodeCount;
        }
    }

    return;
}

void C_cavity_receiver::meshPolygon(const util::matrix_t<double>& poly, double elemSize)
{
    double almostZero = 1.E-7;

    // Confirm correct inputs
    size_t n_verts = poly.nrows();
    size_t n_dims = poly.ncols();

    util::matrix_t<double> nHat;
    util::matrix_t<double> less1_0(1, 3, std::numeric_limits<double>::quiet_NaN());
    if (n_dims == 3) {

        if (n_verts < 3) {
            throw(C_csp_exception("meshPolygon requires at least 3 vertices"));
        }

        // test for coplanar vertices
        util::matrix_t<double> poly0 = poly.row(0);
        util::matrix_t<double> poly1 = poly.row(1);
        util::matrix_t<double> poly2 = poly.row(2);

        util::matrix_t<double> less2_0(1, 3, std::numeric_limits<double>::quiet_NaN());

        diffrows(poly2, poly1, less2_0);
        diffrows(poly1, poly0, less1_0);

        util::matrix_t<double> n;
        crossproduct(less1_0, less2_0, n);

        norm3Dvect(n, nHat);

        if (n_verts > 3) { // check that points are planar

            double volume = 0.0;
            util::matrix_t<double> diff_local(1, 3, std::numeric_limits<double>::quiet_NaN());
            // Only need to check points after the first 3
            for (size_t i = 3; i < n_verts; i++) {
                // the triple product of any combination of vertices must be zero for the polygon to be planar
                // only need to check after first three points used to calculate the normal
                for (size_t j = 0; j < 3; j++) {
                    diff_local(0, j) = poly(i, j) - poly(0, j);
                }
                volume = std::abs(dotprod3D(n, diff_local));

                if (volume > almostZero) {
                    throw(C_csp_exception("meshPolygon polygon vertices not coplanar"));
                }
            }
        }
    }
    else {
        throw(C_csp_exception("meshMapped requires 3 dimensions for a vortex"));
    }

    util::matrix_t<double> center;
    ave_columns(poly, center);

    util::matrix_t<double> poly_2D;
    util::matrix_t<double> poly_rt;
    to2D(poly, center, nHat, less1_0, poly_2D, poly_rt);

    util::matrix_t<double> max_vect;
    util::matrix_t<double> min_vect;
    min_max_vects_from_columns(poly_2D, max_vect, min_vect);

    size_t n_cols = poly_2D.ncols();
    util::matrix_t<double> bbox(2, n_cols, std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < n_cols; i++) {
        bbox(0, i) = min_vect(0, i);
        bbox(1, i) = max_vect(0, i);
    }

    util::matrix_t<double> max_less_min;
    diffrows(max_vect, min_vect, max_less_min);
    double maxDim = max_row_value(max_less_min);

    if (maxDim / elemSize < 3 || maxDim / elemSize > 30) {
        throw(C_csp_exception("meshPolygon: Element size not within the required range"));
    }

    // evenly distribute mesh points on edges
    util::matrix_t<double> pfix_local = poly_2D;
    util::matrix_t<double> pfix = pfix_local;
    for (size_t j = 0; j < n_verts; j++) {

        pfix_local = pfix;

        util::matrix_t<double> pointA = poly_2D.row(j);
        util::matrix_t<double> pointB;
        if (j < n_verts - 1) {
            pointB = poly_2D.row(j+1);
        }
        else{
            pointB = poly_2D.row(0);
        }

        util::matrix_t<double> BlessA;
        diffrows(pointB, pointA, BlessA);
        double edgeLength = mag_vect(BlessA);

        // Determine number of elements on each edge
        int edgeDivs = max(1, (int)std::round(edgeLength / elemSize));

        // divide up edges into mesh points
        if (edgeDivs > 1) {
            util::matrix_t<double> AtoB;
            scale_vect(BlessA, 1./edgeLength, AtoB);
            double segment = edgeLength / (double)edgeDivs;

            util::matrix_t<double> newPoints(edgeDivs - 1, 2, 0.0);

            for (size_t i = 0; i < edgeDivs - 1; i++) {
                for (size_t k = 0; k < 2; k++) {
                    newPoints(i, k) = pointA(0, k) + (i + 1) * segment * AtoB(0, k);
                }
            }

            size_t n_row_pfix = pfix_local.nrows();
            size_t n_row_newPoints = newPoints.nrows();
            size_t n_row_pfix_new = n_row_pfix + n_row_newPoints;

            pfix.resize_preserve(n_row_pfix_new, 2, std::numeric_limits<double>::quiet_NaN());
            for (size_t i = 0; i < n_row_newPoints; i++) {
                for (size_t k = 0; k < 2; k++) {
                    pfix(n_row_pfix + i,k) = newPoints(i,k);
                }
            }
        }
        else {
            pfix = pfix_local;
        }
    }

    // Call the mesh engine
    triMesh2D(elemSize, bbox, pfix, poly_2D);

    //% call the mesh engine
    //    [nodes2D, triangles] = triMesh2D(fd, @huniform, elemSize, bbox, pfix);


}

void C_cavity_receiver::triMesh2D(double h0, const util::matrix_t<double>& bbox, const util::matrix_t<double>& pfix,
    const util::matrix_t<double>& poly_2D)
{
    // function settings
    double dptol = .001;
    double ttol = .1;
    double Fscale = 1.2;
    double deltat = .2;
    double geps = .001 * h0;
    double deps = sqrt(pow(2., -52)) * h0;


    // 1. Create initial distribution in bounding box(equilateral triangles)
    std::vector<double> x_mg;
    for (double ix = bbox(0, 0); ix < bbox(1, 0); ix = ix + h0) {
        x_mg.push_back(ix);
    }
    std::vector<double> y_mg;
    for (double iy = bbox(0, 1); iy < bbox(1, 1); iy = iy + h0 * sqrt(3) / 2.) {
        y_mg.push_back(iy);
    }

    size_t n_x_mg = x_mg.size();
    size_t n_y_mg = y_mg.size();

    util::matrix_t<double> x(n_y_mg, n_x_mg, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> y(n_y_mg, n_x_mg, std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < n_y_mg; i++) {
        for (size_t j = 0; j < n_x_mg; j++) {
            x(i,j) = x_mg[j];
            y(i,j) = y_mg[i];
        }
    }

    for (size_t i = 0; i < x.nrows(); i++) {
        if (i % 2 == 1) {
            for (size_t j = 0; j < x.ncols(); j++) {
                x(i, j) = x(i, j) + h0 / 2.0;       // Shift odd (even in 1-based indices) rows
            }
        }
    }

    util::matrix_t<double> p(n_y_mg*n_x_mg, 2, std::numeric_limits<double>::quiet_NaN());
    for(size_t j = 0; j < x.ncols(); j++){
        for (size_t i = 0; i < x.nrows(); i++) {
            p(j*x.nrows() + i, 0) = x(i,j);
            p(j*x.nrows() + i, 1) = y(i,j);
        }
    }

    // 2. Remove points outside the region, apply the rejection method
    util::matrix_t<double> d;
    pointToPoly(p, poly_2D, d);

    std::vector<size_t> i_p;
    for (size_t i = 0; i < d.nrows(); i++) {
        if (d(i, 0) + h0 / 2 < geps) {
            i_p.push_back(i);
        }
    }

    util::matrix_t<double> p_temp(i_p.size(), 2, std::numeric_limits<double>::quiet_NaN());
    for (size_t k = 0; k < i_p.size(); k++) {
        size_t i = i_p[k];
        p_temp(k,0) = p(i,0);
        p_temp(k,1) = p(i,1);
    }
    p = p_temp;

    std::vector<size_t> i_p_include;
    bool is_unique = true;
    for (size_t i = 0; i < p.nrows(); i++) {
        is_unique = true;
        for (size_t k = 0; k < pfix.nrows(); k++) {
            if (p(i, 0) == pfix(k, 0)) {
                if (p(i, 1) == pfix(k, 1)) {
                    is_unique = false;
                    break;
                }
            }
        }
        if (is_unique) {
            i_p_include.push_back(i);
        }
    }

    size_t n_row_pfix = pfix.nrows();
    size_t n_row_p = i_p_include.size();
    size_t n_row_p_new = n_row_pfix + n_row_p;
    p_temp = pfix;
    p_temp.resize_preserve(n_row_p_new, pfix.ncols(), std::numeric_limits<double>::quiet_NaN());

    for (size_t i = 0; i < i_p_include.size(); i++) {
        p_temp(n_row_pfix + i, 0) = p(i_p_include[i],0);
        p_temp(n_row_pfix + i, 1) = p(i_p_include[i],1);
    }

    p = p_temp;     // unordered

    vector<pair<double, double>> v_pair_p;
    for (size_t i = 0; i < p.nrows(); i++) {
        v_pair_p.push_back(make_pair(p(i,0),p(i,1)));
    }

    // using modified sort() function to sort
    sort(v_pair_p.rbegin(), v_pair_p.rend(), sort_pair_ascending);

    for (size_t i = 0; i < p.nrows(); i++) {
        p(i,0) = v_pair_p[i].first;
        p(i, 1) = v_pair_p[i].second;
    }

    util::matrix_t<double> pold(p.nrows(), p.ncols(), std::numeric_limits<double>::quiet_NaN());
    int count = 0;

    while (true) {

        // 3. Retriangulation by the Delaunay algorithm
        pold = p;       // Save current positions

    }

    /*while true
        % 3. Retriangulation by the Delaunay algorithm
        if max(sqrt(sum((p - pold). ^ 2, 2)) / h0) > ttol% Any large movement ?
            pold = p;% Save current positions
            t = delaunayn(p);% List of triangles
            pmid = (p(t(:, 1), :) + p(t(:, 2), :) + p(t(:, 3), :)) / 3; % Compute centroids
            t = t(feval(fd, pmid, varargin{ : }) < -geps, :);% Keep interior triangles
            % 4. Describe each bar by a unique pair of nodes
            bars = [t(:, [1, 2]); t(:, [1, 3]); t(:, [2, 3])];% Interior bars duplicated
            bars = unique(sort(bars, 2), 'rows');% Bars as node pairs
            % 5. Graphical output of the current mesh
            if plotFormation
                trimesh(t, p(:, 1), p(:, 2), zeros(N, 1))
                view(2); axis equal; axis off; drawnow
                end
                end*/


    return;
}

double C_cavity_receiver::pointToLine(const util::matrix_t<double>& p, const util::matrix_t<double>& a,
    const util::matrix_t<double>& b)
{
    // find the distance between point pand line segment a - b
    double x = p(0,0);
    double y = p(0,1);
    double x1 = a(0,0);
    double y1 = a(0,1);
    double x2 = b(0,0);
    double y2 = b(0,1);

    double A = x - x1;
    double B = y - y1;
    double C = x2 - x1;
    double D = y2 - y1;

    double dott = A*C + B*D;
    double len_sq = C*C + D*D;

    double param = -1.0;
    if (len_sq != 0) {
        param = dott / len_sq;
    }

    double xx = std::numeric_limits<double>::quiet_NaN();
    double yy = std::numeric_limits<double>::quiet_NaN();

    if (param < 0) {
        xx = x1;
        yy = y1;
    }
    else if (param > 1) {
        xx = x2;
        yy = y2;
    }
    else {
        xx = x1 + param*C;
        yy = y1 + param*D;
    }

    double dx = x - xx;
    double dy = y - yy;

    double val = sqrt(dx*dx + dy*dy);
    return val;
}

void C_cavity_receiver::pointToPoly(const util::matrix_t<double>& p, const util::matrix_t<double>& POLY,
    util::matrix_t<double>& d)
{
    int n = p.nrows();
    int m = p.ncols();
    int N = POLY.nrows();
    int M = POLY.ncols();

    if (m == 2 && M == 2) {

        d.resize_fill(n, 1, 0.0);

        for (size_t i = 0; i < n; i++) {

            util::matrix_t<double> D(N, 1, 0.0);

            for (size_t j = 0; j < N; j++) {
                util::matrix_t<double> a = POLY.row(j);

                util::matrix_t<double> b;
                if (j < N-1) {
                    b = POLY.row(j+1);
                }
                else {
                    b = POLY.row(0);
                }

                D(j,0) = std::abs(pointToLine(p.row(i), a, b));
            }

            d(i,0) = min_val_first_colum(D);
        }

        double abc = 1.23;
        util::matrix_t<double> p_x = p.col(0);
        util::matrix_t<double> p_y = p.col(1);

        util::matrix_t<double> poly_x;
        util::matrix_t<double> poly_y;

        transpose_matrix_t(POLY.col(0), poly_x);
        transpose_matrix_t(POLY.col(1), poly_y);

        util::matrix_t<int> in;
        inpolygon(p_x, p_y, poly_x, poly_y, in);

        for (size_t i = 0; i < n; i++) {
            d(i,0) *= (-1.0*in(i,0) + (double)(!in(i,0)));
        }
    }
    else {
        throw(C_csp_exception("pointToPoly: incorrect dimensions"));
    }
}

void C_cavity_receiver::inpolygon(const util::matrix_t<double>& p_x, const util::matrix_t<double>& p_y,
    const util::matrix_t<double>& poly_x, const util::matrix_t<double>& poly_y,
    util::matrix_t<int>& is_in_polygon)
{
    util::matrix_t<double> x = p_x;
    util::matrix_t<double> y = p_y;

    // Last point in polygon should equal first
    util::matrix_t<double> vx = poly_x;
    util::matrix_t<double> vy = poly_y;
    if (poly_x(poly_x.nrows() - 1, 0) != poly_x(0, 0) || poly_y(poly_y.nrows() - 1, 0)) {
        vx.resize_preserve(poly_x.nrows()+1, 1, std::numeric_limits<double>::quiet_NaN());
        vx(poly_x.nrows(),0) = poly_x(0,0);
        vy.resize_preserve(poly_y.nrows()+1, 1, std::numeric_limits<double>::quiet_NaN());
        vy(poly_y.nrows(),0) = poly_y(0,0);
    }

    size_t n_verts = vx.nrows();
    size_t n_points = x.ncols();

    util::matrix_t<double> xx(n_verts, n_points, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> yy(n_verts, n_points, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> v_vx(n_verts, n_points, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> v_vy(n_verts, n_points, std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < n_verts; i++) {
        for (size_t j = 0; j < n_points; j++) {
            xx(i,j) = x(0,j);
            yy(i,j) = y(0,j);
            v_vx(i,j) = vx(i,0);
            v_vy(i,j) = vy(i,0);
        }
    }

    x = xx;
    y = yy;
    vx = v_vx;
    vy = v_vy;

    util::matrix_t<double> avx(n_verts - 1, 1, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> avy(n_verts - 1, 1, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> ScaleFactor(n_verts - 1, 1, std::numeric_limits<double>::quiet_NaN());

    for (size_t i = 0; i < n_verts - 1; i++) {
        avx(i,0) = std::abs(0.5*(vx(i,0) + vx(i+1,0)));
        avy(i,0) = std::abs(0.5*(vy(i,0) + vy(i+1,0)));
        ScaleFactor(i,0) = max(max(avx(i,0),avy(i,0)), avx(i,0)*avy(i,0));
    }

    for (size_t i = 0; i < n_verts; i++) {
        for (size_t j = 0; j < n_points; j++) {
            vx(i,j) -= x(i,j);
            vy(i,j) -= y(i,j);
        }
    }

    util::matrix_t<int> quad(n_verts, n_points, std::numeric_limits<double>::quiet_NaN());
    bool posX, posY, negX, negY;
    for (size_t i = 0; i < n_verts; i++) {
        for (size_t j = 0; j < n_points; j++) {
            posX = vx(i, j) > 1.E-10;
            posY = vy(i, j) > 1.E-10;
            negX = !posX;
            negY = !posY;
            quad(i,j) = (size_t)(negX && posY) + 2*(size_t)(negX && negY) + 3*(size_t)(posX && negY);
        }
    }

    double scaledEps_base = sqrt(pow(2., -52)) * 3.0;
    double scaledEps = std::numeric_limits<double>::quiet_NaN();


    util::matrix_t<double> theCrossProd(n_verts-1, n_points, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> dotProd(n_verts - 1, n_points, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<int> signCrossProd(n_verts - 1, n_points, -99);
    util::matrix_t<int> diffQuad(n_verts - 1, n_points, -99);
    for (size_t i = 0; i < n_verts - 1; i++) {
        for (size_t j = 0; j < n_points; j++) {

            // Compute the sign() of the cross productand dot product of adjacent vertices.
            theCrossProd(i,j) = vx(i,j)*vy(i+1,j) - vx(i+1,j)*vy(i,j);

            scaledEps = ScaleFactor(i, 0) * scaledEps_base;
            if (abs(theCrossProd(i, j)) < scaledEps) {
                signCrossProd(i,j) = 0;
            }
            else if (theCrossProd(i, j) > 1.E-10) {
                signCrossProd(i,j) = 1;
            }
            else if (theCrossProd(i, j) < -1.E-10) {
                signCrossProd(i, j) = -1;
            }
            else {
                signCrossProd(i, j) = 0;
            }

            dotProd(i,j) = vx(i,j)*vx(i+1,j) + vy(i,j)*vy(i+1,j);

            diffQuad(i,j) = quad(i+1,j) - quad(i,j);

            if (abs(diffQuad(i, j)) == 3) {
                diffQuad(i,j) /= -3;
            }
            if (abs(diffQuad(i, j)) == 2) {
                diffQuad(i,j) = 2 * signCrossProd(i,j);
            }
        }        
    }

    util::matrix_t<int> in(1, n_points, -1);

    sum_int_columns(diffQuad, in);

    util::matrix_t<int> on(1, n_points, -1);

    size_t is_true = 0;
    for (size_t j = 0; j < n_points; j++) {

        if (in(0, j) != 0) {
            in(0, j) = 1;
        }
        else {
            in(0, j) = 0;
        }

        is_true = 0;
        for (size_t i = 0; i < n_verts - 1; i++) {
            if (signCrossProd(i, j) == 0 && dotProd(i, j) <= 0.0) {
                is_true = 1;
            }
        }
        on(0,j) = is_true;
        if (in(0,j) == 1 || on(0,j) == 1) {
            in(0,j) = 1;
        }
        else {
            in(0,j) = 0;
        }
    }

    transpose_int_matrix_t(in, is_in_polygon);

    return;
}


void C_cavity_receiver::transpose_matrix_t(const util::matrix_t<double>& a, util::matrix_t<double>& b)
{
    size_t n_row = a.nrows();
    size_t n_col = a.ncols();
    b.resize_fill(n_col, n_row, std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < n_row; i++) {
        for (size_t j = 0; j < n_col; j++) {
            b(j,i) = a(i,j);
        }
    }
}

void C_cavity_receiver::transpose_int_matrix_t(const util::matrix_t<int>& a, util::matrix_t<int>& b)
{
    size_t n_row = a.nrows();
    size_t n_col = a.ncols();
    b.resize_fill(n_col, n_row, std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < n_row; i++) {
        for (size_t j = 0; j < n_col; j++) {
            b(j, i) = a(i, j);
        }
    }
}


void C_cavity_receiver::meshMapped(const util::matrix_t<double>& poly, double elemSize,
    util::matrix_t<double>& nodes, util::matrix_t<double>& quads)
{
    double almostZero = 1.E-7;

    // Confirm correct inputs
    size_t n_verts = poly.nrows();
    size_t n_dims = poly.ncols();

    util::matrix_t<double> less1_0(1, 3, std::numeric_limits<double>::quiet_NaN());
    util::matrix_t<double> n_hat;
    if (n_dims == 3) {

        // Test for coplanar vertices
        if (n_verts != 4) {
            throw(C_csp_exception("meshMapped requires 4 vertices"));
        }
        else {
            util::matrix_t<double> less2_0(1, 3, std::numeric_limits<double>::quiet_NaN());
            for (size_t i = 0; i < 3; i++) {
                less1_0(0, i) = poly(1, i) - poly(0, i);
                less2_0(0, i) = poly(2, i) - poly(0, i);
            }
            util::matrix_t<double> n;
            crossproduct(less1_0, less2_0, n);
            norm3Dvect(n, n_hat);
            double volume = 0.0;
            util::matrix_t<double> diff_local(1, 3, std::numeric_limits<double>::quiet_NaN());
            for (size_t i = 3; i < n_verts; i++) {
                // the triple product of any combination of vertices must be zero for the polygon to be planar
                // only need to check after first three points used to calculate the normal
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

    util::matrix_t<double> center;
    ave_columns(poly, center);

    util::matrix_t<double> poly_2D;
    util::matrix_t<double> poly_rt;
    to2D(poly, center, n_hat, less1_0, poly_2D, poly_rt);

    util::matrix_t<double> nodes2D;
    map(poly_2D, elemSize, nodes2D, quads);

    to3D(nodes2D, center, n_hat, less1_0, nodes);

    return;
}

void C_cavity_receiver::to3D(const util::matrix_t<double>& poly_xy, const util::matrix_t<double>& origin,
    const util::matrix_t<double>& normal, const util::matrix_t<double>& xaxis,
    util::matrix_t<double>& poly3d)
{
    size_t n = poly_xy.nrows();     // number of points to process

    util::matrix_t<double> nHat;
    norm3Dvect(normal, nHat);

    util::matrix_t<double> xHat;
    norm3Dvect(xaxis, xHat);

    util::matrix_t<double> yHat;
    crossproduct(nHat, xHat, yHat);

    poly3d.resize_fill(n, 3, 0.0);
    for (size_t i = 0; i < n; i++){
        for (size_t j = 0; j < 3; j++) {
            poly3d(i,j) = origin(0,j) + xHat(0,j)*poly_xy(i,0) + yHat(0,j)*poly_xy(i,1);
        }
    }

    return;
}

void C_cavity_receiver::map(const util::matrix_t<double>& poly2D, double elemSize,
    util::matrix_t<double>& nodes, util::matrix_t<double>& quads)
{
    util::matrix_t<double> A = poly2D.row(0);
    util::matrix_t<double> B = poly2D.row(1);
    util::matrix_t<double> C = poly2D.row(2);
    util::matrix_t<double> D = poly2D.row(3);

    util::matrix_t<double> AtoB;
    diffrows(B, A, AtoB);
    double lengthAB = mag_vect(AtoB);

    util::matrix_t<double> BtoC;
    diffrows(C, B, BtoC);
    double lengthBC = mag_vect(BtoC);

    util::matrix_t<double> CtoD;
    diffrows(D, C, CtoD);
    double lengthCD = mag_vect(CtoD);

    util::matrix_t<double> AtoD;
    diffrows(D, A, AtoD);
    double lengthDA = mag_vect(AtoD);

    double maxDim = max({ lengthAB, lengthBC, lengthCD, lengthDA });

    // test for a reasonable element size
    if (maxDim / elemSize < 0.5 || maxDim / elemSize > 250) {
        throw(C_csp_exception("Element size not within the required range"));
    }

    // find reasonable # of elements per edge for mapping
    int elemsM = std::max(1, (int) std::round(0.5/elemSize*(lengthAB + lengthCD)));
    int elemsN = std::max(1, (int) std::round(0.5/elemSize*(lengthBC + lengthDA)));

    // initialize nodeand element arrays
    nodes.resize_fill((elemsM+1)*(elemsN+1), 2, 0.0);
    quads.resize_fill(elemsM*elemsN, 4, 0.0);

    int nodeID = -1;     // last used node index - set to -1 for C++ 0-based index
    int elemID = -1;     // last used element index - set to -1 for C++ 0-based index

    // Matlab code also starts at index 0
    for (size_t m = 0; m < elemsM + 1; m++) {
        // create bridge EF between AB and CD
        util::matrix_t<double> scaledAtoB;
        scale_vect(AtoB, m/(double)elemsM, scaledAtoB);
        util::matrix_t<double> E;
        add_vect_rows(A, scaledAtoB, E);

        util::matrix_t<double> scaledCtoD;
        scale_vect(CtoD, m/(double)elemsM, scaledCtoD);
        util::matrix_t<double> F;
        diffrows(D, scaledCtoD, F);

        util::matrix_t<double> EtoF;
        diffrows(F, E, EtoF);

        // Matlab code also starts at index 0
        for (size_t n = 0; n < elemsN + 1; n++) {
            // walk along bridge EF defining nodes
            nodeID++;
            util::matrix_t<double> scaledEtoF;
            scale_vect(EtoF, n/(double)elemsN, scaledEtoF);
            util::matrix_t<double> E_plus_scaledEtoF;
            add_vect_rows(E, scaledEtoF, E_plus_scaledEtoF);

            for (size_t i = 0; i < E_plus_scaledEtoF.ncols(); i++) {
                nodes(nodeID, i) = E_plus_scaledEtoF(0, i);
            }

            if (m > 0 && n > 0) { // then define new element
                elemID++;

                // Retrieve nodes around element in CCW direction
                // these MUST go around the element in order(CW or CCW) in
                // order to work correctly with viewFactor(...)
                quads(elemID, 0) = nodeID - elemsN - 2;
                quads(elemID, 1) = nodeID - 1;
                quads(elemID, 2) = nodeID;
                quads(elemID, 3) = nodeID - elemsN - 1;
            }

        }
    }

    return;
}

void C_cavity_receiver::to2D(const util::matrix_t<double>& poly, const util::matrix_t<double>& center,
    const util::matrix_t<double>& normal, const util::matrix_t<double>& xaxis,
    util::matrix_t<double>& poly_xy, util::matrix_t<double>& poly_rt)
{
    size_t n = poly.nrows();

    util::matrix_t<double> nHat;
    norm3Dvect(normal, nHat);

    util::matrix_t<double> xHat;
    norm3Dvect(xaxis, xHat);

    util::matrix_t<double> yHat;
    crossproduct(nHat, xHat, yHat);

    poly_xy.resize_fill(n, 2, 0.0);
    poly_rt.resize_fill(n, 2, 0.0);

    for (size_t i = 0; i < n; i++) {
        util::matrix_t<double> point = poly.row(i);
        util::matrix_t<double> arm(1, 3);
        for (size_t j = 0; j < 3; j++) {
            arm(0,j) = point(0,j) - center(0,j);
        }
        double radius = mag_vect(arm);
        double xComp = dotprod3D(arm, xHat);    // x coordinate in 2D CS
        double yComp = dotprod3D(arm, yHat);    // y coordinate in 2D CS
        double theta = atan2(yComp, xComp);     //
        if (theta < 0.0) {
            theta = theta + 2. * CSP::pi;
        }
        poly_xy(i, 0) = xComp;
        poly_xy(i, 1) = yComp;
        poly_rt(i, 0) = radius;
        poly_rt(i, 1) = theta;
    }

    return;
}

void C_cavity_receiver::add_constant_to_each_element(double val, util::matrix_t<double>& a)
{
    for (size_t i = 0; i < a.nrows(); i++) {
        for (size_t j = 0; j < a.ncols(); j++) {
            a(i,j) += val;
        }
    }
}

double C_cavity_receiver::dotprod3D(const util::matrix_t<double>& a, const util::matrix_t<double>& b)
{
    return a(0,0)*b(0,0) + a(0,1)*b(0,1) + a(0,2)*b(0,2);
}

void C_cavity_receiver::scale_vect(const util::matrix_t<double>& a, double scale, util::matrix_t<double>& out_vect)
{
    out_vect = a;
    for (size_t i = 0; i < a.ncols(); i++) {
        out_vect(0,i) = a(0,i)*scale;
    }
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
    double magnitude = mag_vect(vector_in);
    for (size_t i = 0; i < 3; i++) {
        norm_vect(0, i) = vector_in(0, i) / magnitude;
    }
}

double C_cavity_receiver::mag_vect(const util::matrix_t<double>& vector_in)
{
    double sum_of_sq = 0.0;
    for (size_t i = 0; i < vector_in.ncols(); i++) {
        sum_of_sq += std::pow(vector_in(0, i), 2);
    }

    return std::sqrt(sum_of_sq);

    //return std::sqrt(std::pow(vector_in(0, 0), 2) + std::pow(vector_in(0, 1), 2) + std::pow(vector_in(0, 2), 2));
}

void C_cavity_receiver::sumcolumns(const util::matrix_t<double>& a, util::matrix_t<double>& summed)
{
    size_t ncols = a.ncols();
    summed.resize_fill(1, ncols, 0.0);

    for (size_t i = 0; i < a.nrows(); i++) {
        for (size_t j = 0; j < ncols; j++) {
            summed(0, j) += a(i, j);
        }
    }
}

void C_cavity_receiver::sum_int_columns(const util::matrix_t<int>& a, util::matrix_t<int>& summed)
{
    size_t ncols = a.ncols();
    summed.resize_fill(1, ncols, 0.0);

    for (size_t i = 0; i < a.nrows(); i++) {
        for (size_t j = 0; j < ncols; j++) {
            summed(0, j) += a(i, j);
        }
    }
}

void C_cavity_receiver::diffrows(const util::matrix_t<double>& a, const util::matrix_t<double>& b, util::matrix_t<double>& a_less_b)
{
    a_less_b.resize_fill(1, a.ncols(), std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < a.ncols(); i++) {
        a_less_b(0, i) = a(0, i) - b(0, i);
    }
}

void C_cavity_receiver::add_vect_rows(const util::matrix_t<double>& a, const util::matrix_t<double>& b, util::matrix_t<double>& a_plus_b)
{
    a_plus_b.resize_fill(1, a.ncols(), std::numeric_limits<double>::quiet_NaN());
    for (size_t i = 0; i < a.ncols(); i++) {
        a_plus_b(0,i) = a(0,i) + b(0,i);
    }
}

void C_cavity_receiver::ave_columns(const util::matrix_t<double>& a, util::matrix_t<double>& averaged)
{
    double nrows = (double)a.nrows();
    sumcolumns(a, averaged);
    for (size_t i = 0; i < 3; i++) {
        averaged(0, i) /= nrows;
    }
}

double C_cavity_receiver::max_row_value(const util::matrix_t<double>& a)
{
    double maxval = a(0, 0);
    for (size_t i = 1; i < a.ncols(); i++) {
        maxval = max(maxval, a(0,i));
    }
    return maxval;
}

double C_cavity_receiver::min_val_first_colum(const util::matrix_t<double>& a)
{
    double minval = a(0, 0);
    for (size_t i = 1; i < a.nrows(); i++) {
        minval = min(minval, a(i, 0));
    }
    return minval;
}

void C_cavity_receiver::min_max_vects_from_columns(const util::matrix_t<double>& a, util::matrix_t<double>& max_vect, util::matrix_t<double>& min_vect)
{
    size_t ncols = a.ncols();
    max_vect = a.row(0);
    min_vect = a.row(0);
    for (size_t i = 1; i < a.nrows(); i++) {
        for (size_t j = 0; j < ncols; j++) {
            max_vect(0,j) = max(max_vect(0,j), a(i,j));
            min_vect(0,j) = min(min_vect(0,j), a(i,j));
        }
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
